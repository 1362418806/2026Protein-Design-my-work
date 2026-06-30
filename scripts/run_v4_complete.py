from __future__ import annotations

import argparse
import json
import math
import re
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from synbio_gfp_v3.candidates import generate_candidates
from synbio_gfp_v3.config import load_config, resolve_config
from synbio_gfp_v3.data import load_exclusion, load_training_table, parse_fasta_like
from synbio_gfp_v3.features import featurize
from synbio_gfp_v3.literature_priors import add_literature_priors
from synbio_gfp_v3.model_cache import ModelCache
from synbio_gfp_v3.models import predict_bundle, train_predictive_models
from synbio_gfp_v3.scoring import add_objective_scores, add_proxy_scores, rank01
from synbio_gfp_v3.site_policy import (
    BARREL_CORE_PRIOR,
    CHROMOPHORE,
    CHROMOPHORE_POCKET_8A,
    ION_NETWORK,
    build_site_policy,
)
from synbio_gfp_v3.structure_filter import compute_pdb_metrics
from synbio_gfp_v3.utils import ensure_jsonable, setup_logger, set_seed, validate_sequence

AA20 = set("ACDEFGHIKLMNPQRSTVWY")
MUT_RE = re.compile(r"([ACDEFGHIKLMNPQRSTVWY])(\d{1,3})([ACDEFGHIKLMNPQRSTVWY])")


def _resolve_file(data_dir: Path, filename: str) -> Path:
    p = data_dir / filename
    if p.exists():
        return p
    candidates = list(data_dir.glob(filename))
    if candidates:
        return candidates[0]
    raise FileNotFoundError(f"Cannot find {filename} under {data_dir}")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Run complete SynBio GFP V4 pipeline")
    p.add_argument("--config", required=True)
    p.add_argument("--data-dir", required=True)
    p.add_argument("--team-name", required=True)
    p.add_argument("--feature-mode", choices=["simple", "esm"], default="simple")
    p.add_argument("--esm-model", default="esm2_t30_150M_UR50D")
    p.add_argument("--max-train-samples", type=int, default=20000)
    p.add_argument("--n-candidates", type=int, default=20000)
    p.add_argument("--out-dir", default="outputs/run_v4_complete")
    p.add_argument("--cache-dir", default="cache")
    p.add_argument("--references", default=None)
    p.add_argument("--exclusion-list", default=None)
    p.add_argument("--gfp-data", default=None)
    p.add_argument("--reference-pdb", default=None)
    p.add_argument("--predicted-pdb-dir", default=None)
    p.add_argument("--structure-metrics-csv", default=None)
    p.add_argument("--proteinmpnn-score-csv", default=None)
    p.add_argument("--skip-structure", action="store_true")
    p.add_argument("--final-k", type=int, default=None)
    p.add_argument("--lambda-std", type=float, default=None)
    p.add_argument("--allow-proxy-final", action="store_true", help="Allow final Top-K without real AF/PDB metrics. Useful only for debugging.")
    return p.parse_args()


def resolve_v4_config(raw: dict[str, Any], args: argparse.Namespace) -> dict[str, Any]:
    cfg = resolve_config(raw, args)
    cfg.setdefault("v4", {})
    cfg.setdefault("structure", {})
    cfg.setdefault("proteinmpnn", {})
    cfg.setdefault("af_branch", {})
    if args.final_k is not None:
        cfg["v4"]["final_k"] = args.final_k
    if args.lambda_std is not None:
        cfg["v4"]["lambda_std"] = args.lambda_std
    if args.allow_proxy_final:
        cfg["v4"]["require_real_structure_for_final"] = False
    if args.structure_metrics_csv:
        cfg["structure"]["structure_metrics_csv"] = args.structure_metrics_csv
    if args.proteinmpnn_score_csv:
        cfg["proteinmpnn"]["score_csv"] = args.proteinmpnn_score_csv
    if args.skip_structure:
        cfg["structure"]["enabled"] = False
    return cfg


def mutation_list_from_sequences(wt: str, seq: str) -> list[str]:
    return [f"{a}{i}{b}" for i, (a, b) in enumerate(zip(wt, seq), start=1) if a != b]


def parse_mutation_text(text: Any) -> list[str]:
    if not isinstance(text, str):
        return []
    return [f"{a}{int(pos)}{b}" for a, pos, b in MUT_RE.findall(text.upper())]


def mutation_position(mut: str) -> int | None:
    m = MUT_RE.fullmatch(mut)
    return int(m.group(2)) if m else None


def mutation_alt(mut: str) -> str | None:
    m = MUT_RE.fullmatch(mut)
    return m.group(3) if m else None


def minmax01(s: pd.Series, higher_is_better: bool = True, default: float = 0.5) -> pd.Series:
    x = pd.to_numeric(s, errors="coerce")
    if x.notna().sum() == 0:
        return pd.Series(default, index=s.index, dtype=float)
    lo = float(x.min())
    hi = float(x.max())
    if math.isclose(lo, hi):
        return pd.Series(default, index=s.index, dtype=float)
    y = (x - lo) / (hi - lo)
    if not higher_is_better:
        y = 1.0 - y
    return y.fillna(default).clip(0.0, 1.0)


def safe_num(df: pd.DataFrame, col: str, default: float = 0.0) -> pd.Series:
    if col not in df.columns:
        return pd.Series(default, index=df.index, dtype=float)
    return pd.to_numeric(df[col], errors="coerce").fillna(default)


def normalize_v3_columns(df: pd.DataFrame, scaffold: str) -> pd.DataFrame:
    out = df.copy()
    if "source_branch" not in out.columns:
        out["source_branch"] = out.get("source", "unknown")
    if "candidate_id" not in out.columns:
        out["candidate_id"] = [f"cand_{i:06d}" for i in range(len(out))]
    if "Seq_ID" not in out.columns:
        out["Seq_ID"] = [f"candidate_{i}" for i in range(len(out))]
    if "mutations" not in out.columns:
        if "mutation_list" in out.columns:
            out["mutations"] = out["mutation_list"].astype(str)
        else:
            out["mutations"] = [";".join(mutation_list_from_sequences(scaffold, s)) for s in out["sequence"].astype(str)]
    if "mutation_count" not in out.columns:
        out["mutation_count"] = [len(parse_mutation_text(x)) for x in out["mutations"]]
    out["mutation_count"] = pd.to_numeric(out["mutation_count"], errors="coerce").fillna(0).astype(int)
    return out


def add_v4_lcb_and_risk(df: pd.DataFrame, cfg: dict[str, Any]) -> pd.DataFrame:
    out = df.copy()
    lam = float(cfg.get("v4", {}).get("lambda_std", 0.75))
    out["brightness_mean"] = safe_num(out, "brightness_pred", default=0.0)
    out["brightness_std"] = safe_num(out, "brightness_uncertainty", default=0.0)
    # LCB is the conservative brightness objective used by V4 to avoid uncertainty-chasing candidates.
    out["brightness_lcb"] = out["brightness_mean"] - lam * out["brightness_std"]
    out["brightness_lcb_norm"] = minmax01(out["brightness_lcb"], higher_is_better=True, default=0.5)
    out["low_risk_score"] = (1.0 - safe_num(out, "low_brightness_risk", default=1.0)).clip(0.0, 1.0)
    return out


def add_epistasis_and_ddg_proxies(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    mut_count = safe_num(out, "mutation_count")
    protected = safe_num(out, "protected_mutation_count")
    pocket = safe_num(out, "pocket_mutation_count")
    ion = safe_num(out, "ion_network_mutation_count")
    core = safe_num(out, "core_mutation_count") if "core_mutation_count" in out.columns else pd.Series(0.0, index=out.index)
    # The pairwise term is intentionally convex because GFP brightness often collapses nonlinearly under multi-site epistasis.
    pairwise = np.maximum(0.0, (mut_count * (mut_count - 1) / 2.0 - 6.0) / 20.0)
    out["pairwise_epistasis_penalty"] = pairwise
    out["ddg_effective"] = (
        0.18 * mut_count
        + 0.85 * protected
        + 0.70 * pocket
        + 0.95 * ion
        + 0.30 * core
        + 0.55 * pairwise
    )
    out["ddg_total_v4"] = out["ddg_effective"] + 0.10 * safe_num(out, "low_brightness_risk")
    out["p_fluorescent_v4"] = 1.0 / (1.0 + np.exp(1.15 * (out["ddg_total_v4"] - 1.45)))
    out["epistasis_penalty"] = (0.22 * pairwise + 0.18 * protected + 0.16 * pocket + 0.22 * ion).clip(0.0, 1.0)
    out["epistasis_score"] = (1.0 - out["epistasis_penalty"]).clip(0.0, 1.0)
    return out


def add_site_counts(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    chrom_counts = []
    core_counts = []
    for text in out["mutations"]:
        muts = parse_mutation_text(text)
        positions = [mutation_position(m) for m in muts]
        chrom_counts.append(sum(p in CHROMOPHORE for p in positions if p is not None))
        core_counts.append(sum(p in BARREL_CORE_PRIOR for p in positions if p is not None))
    out["chromophore_mutation_count"] = chrom_counts
    out["core_mutation_count"] = core_counts
    return out


def load_optional_score_table(path: str | None) -> pd.DataFrame | None:
    if not path:
        return None
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Optional score file not found: {p}")
    return pd.read_excel(p) if p.suffix.lower() in {".xlsx", ".xls"} else pd.read_csv(p)


def merge_proteinmpnn_scores(df: pd.DataFrame, cfg: dict[str, Any]) -> pd.DataFrame:
    out = df.copy()
    pm_cfg = cfg.get("proteinmpnn", {})
    score_df = load_optional_score_table(pm_cfg.get("score_csv"))
    out["proteinmpnn_status"] = "not_provided"
    out["proteinmpnn_score_raw"] = np.nan
    out["proteinmpnn_score_norm"] = 0.5
    if score_df is None or score_df.empty:
        return out
    score_col = pm_cfg.get("score_column", "proteinmpnn_score")
    id_col = next((c for c in ["candidate_id", "Seq_ID", "seq_id", "sequence"] if c in score_df.columns and c in out.columns), None)
    if id_col is None or score_col not in score_df.columns:
        raise ValueError("ProteinMPNN score CSV must share candidate_id/Seq_ID/sequence and contain the configured score column.")
    small = score_df[[id_col, score_col]].rename(columns={score_col: "proteinmpnn_score_raw"})
    out = out.drop(columns=["proteinmpnn_score_raw"], errors="ignore").merge(small, on=id_col, how="left")
    lower_is_better = bool(pm_cfg.get("lower_is_better", True))
    out["proteinmpnn_score_norm"] = minmax01(out["proteinmpnn_score_raw"], higher_is_better=not lower_is_better, default=0.5)
    out["proteinmpnn_status"] = np.where(out["proteinmpnn_score_raw"].notna(), "score_imported", "not_found")
    return out


def find_candidate_pdb(row: pd.Series, pred_dir: str | None, idx: int) -> Path | None:
    if not pred_dir:
        return None
    base = Path(pred_dir)
    names = [
        f"{row.get('candidate_id', '')}.pdb",
        f"{row.get('Seq_ID', '')}.pdb",
        f"candidate_{idx}.pdb",
        f"cand_{idx:06d}.pdb",
    ]
    for name in names:
        if name and (base / name).exists():
            return base / name
    return None


def structure_geometry_from_metrics(df: pd.DataFrame) -> pd.Series:
    barrel = safe_num(df, "barrel_rmsd", default=np.nan)
    pocket = safe_num(df, "pocket_rmsd", default=np.nan)
    ion = safe_num(df, "ion_network_rmsd", default=np.nan)
    barrel_score = np.exp(-np.nan_to_num(barrel, nan=4.0) / 1.20)
    pocket_score = np.exp(-np.nan_to_num(pocket, nan=3.0) / 0.75)
    ion_score = np.exp(-np.nan_to_num(ion, nan=3.0) / 0.65)
    return pd.Series(0.40 * barrel_score + 0.38 * pocket_score + 0.22 * ion_score, index=df.index).clip(0.0, 1.0)


def add_structure_branch(df: pd.DataFrame, cfg: dict[str, Any], out_dir: Path, logger: Any) -> pd.DataFrame:
    out = df.copy()
    st_cfg = cfg.get("structure", {})
    out["has_real_structure"] = False
    out["structure_status"] = out.get("structure_status", "not_checked")
    for col in ["barrel_rmsd", "pocket_rmsd", "ion_network_rmsd"]:
        if col not in out.columns:
            out[col] = np.nan
    metrics_csv = st_cfg.get("structure_metrics_csv")
    if metrics_csv:
        metrics = load_optional_score_table(metrics_csv)
        id_col = next((c for c in ["candidate_id", "Seq_ID", "sequence"] if c in metrics.columns and c in out.columns), None)
        if id_col is None:
            raise ValueError("Structure metrics CSV must share candidate_id, Seq_ID, or sequence with candidate table.")
        keep_cols = [id_col] + [c for c in ["barrel_rmsd", "pocket_rmsd", "ion_network_rmsd", "structure_status", "structure_pass"] if c in metrics.columns]
        out = out.drop(columns=[c for c in keep_cols if c != id_col and c in out.columns], errors="ignore").merge(metrics[keep_cols], on=id_col, how="left")
        out["has_real_structure"] = out[["barrel_rmsd", "pocket_rmsd", "ion_network_rmsd"]].notna().any(axis=1)
        # Pandas fillna does not accept ndarray values; build a Series so imported metrics work across pandas versions.
        default_status = pd.Series(np.where(out["has_real_structure"], "metrics_csv_checked", "not_checked"), index=out.index)
        out["structure_status"] = out["structure_status"].where(out["structure_status"].notna(), default_status)
        logger.info("V4 STRUCTURE METRICS IMPORTED | file=%s | real=%d", metrics_csv, int(out["has_real_structure"].sum()))
    elif bool(st_cfg.get("enabled", True)) and st_cfg.get("reference_pdb") and st_cfg.get("predicted_pdb_dir"):
        ref_pdb = st_cfg.get("reference_pdb")
        records = []
        for idx, row in out.iterrows():
            pdb = find_candidate_pdb(row, st_cfg.get("predicted_pdb_dir"), idx)
            if pdb is None:
                continue
            m = compute_pdb_metrics(ref_pdb, pdb)
            records.append({"_idx": idx, **m})
        for rec in records:
            idx = rec.pop("_idx")
            for k, v in rec.items():
                out.loc[idx, k] = v
            out.loc[idx, "has_real_structure"] = True
        logger.info("V4 PDB STRUCTURE CHECK DONE | real=%d", len(records))
    # Real structure is required for final eligibility unless explicitly disabled; proxy score remains diagnostic only.
    out["structure_geometry_score"] = structure_geometry_from_metrics(out)
    no_real = ~out["has_real_structure"].astype(bool)
    out.loc[no_real, "structure_geometry_score"] = 0.05
    out.loc[no_real & out["structure_status"].isna(), "structure_status"] = "not_checked"
    return out


def add_v4_gates(df: pd.DataFrame, cfg: dict[str, Any]) -> pd.DataFrame:
    out = df.copy()
    g = cfg.get("v4", {}).get("basic_gate", {})
    sg = cfg.get("v4", {}).get("structure_gate", {})
    sel = cfg.get("selection", {})
    min_len = int(sel.get("min_len", 220))
    max_len = int(sel.get("max_len", 250))
    valid_rows = []
    for seq in out["sequence"].astype(str):
        v = validate_sequence(seq, min_len, max_len)
        valid_rows.append(all(v.values()))
    out["contest_sequence_valid"] = valid_rows
    # Basic gate is deliberately independent of AF/ProteinMPNN availability so it can define the structure-priority pool.
    out["v4_basic_gate"] = (
        out["contest_sequence_valid"].astype(bool)
        & (safe_num(out, "protected_mutation_count") <= float(g.get("max_protected_mutations", 0)))
        & (safe_num(out, "pocket_mutation_count") <= float(g.get("max_pocket_mutations", 0)))
        & (safe_num(out, "ion_network_mutation_count") <= float(g.get("max_ion_network_mutations", 0)))
        & (safe_num(out, "chromophore_mutation_count") <= float(g.get("max_chromophore_mutations", 0)))
        & (safe_num(out, "core_mutation_count") <= float(g.get("max_core_mutations", 1)))
        & (safe_num(out, "low_brightness_risk") <= float(g.get("max_low_brightness_risk", 0.35)))
        & (safe_num(out, "mutation_count") <= float(g.get("max_mutation_count", 9)))
    )
    out["structure_pass"] = (
        out["has_real_structure"].astype(bool)
        & (safe_num(out, "barrel_rmsd", default=999) <= float(sg.get("max_barrel_rmsd", 1.50)))
        & (safe_num(out, "pocket_rmsd", default=999) <= float(sg.get("max_pocket_rmsd", 0.85)))
        & (safe_num(out, "ion_network_rmsd", default=999) <= float(sg.get("max_ion_network_rmsd", 0.75)))
        & (safe_num(out, "structure_geometry_score") >= float(sg.get("min_structure_geometry_score", 0.60)))
        & (safe_num(out, "proteinmpnn_score_norm", default=0.5) >= float(sg.get("min_proteinmpnn_score_norm", 0.20)))
    )
    require_real = bool(cfg.get("v4", {}).get("require_real_structure_for_final", True))
    if require_real:
        out["v4_structure_gate"] = out["structure_pass"].astype(bool)
    else:
        # This debug path preserves production columns but lets users smoke-test final selection without AF/PDB results.
        out["v4_structure_gate"] = out["v4_basic_gate"].astype(bool)
    out["v4_final_eligible"] = out["v4_basic_gate"].astype(bool) & out["v4_structure_gate"].astype(bool)
    return out


def add_v4_objective(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["thermal_score"] = minmax01(safe_num(out, "tgp_surface_prior"), True, 0.5)
    out["stability_score"] = minmax01(safe_num(out, "stability_proxy"), True, 0.5)
    out["foldability_score"] = minmax01(safe_num(out, "foldability_proxy"), True, 0.5)
    out["literature_score"] = minmax01(safe_num(out, "literature_prior"), True, 0.5)
    out["mutation_burden_score"] = (1.0 - np.clip((safe_num(out, "mutation_count") - 2.0) / 10.0, 0.0, 1.0)).clip(0.0, 1.0)
    out["structure_score_v4"] = (0.75 * safe_num(out, "structure_geometry_score", 0.05) + 0.25 * safe_num(out, "proteinmpnn_score_norm", 0.5)).clip(0.0, 1.0)
    # V4 objective combines conservative brightness, low-risk probability, stability priors, structural plausibility, and epistasis control.
    out["v4_objective"] = (
        0.26 * out["brightness_lcb_norm"]
        + 0.17 * out["low_risk_score"]
        + 0.13 * out["thermal_score"]
        + 0.12 * out["stability_score"]
        + 0.08 * out["foldability_score"]
        + 0.10 * out["literature_score"]
        + 0.08 * out["structure_score_v4"]
        + 0.04 * out["epistasis_score"]
        + 0.02 * out["mutation_burden_score"]
    )
    rank_cols = {
        "rank_brightness_lcb_norm": "brightness_lcb_norm",
        "rank_low_risk_score": "low_risk_score",
        "rank_thermal_score": "thermal_score",
        "rank_stability_score": "stability_score",
        "rank_foldability_score": "foldability_score",
        "rank_literature_score": "literature_score",
        "rank_structure_geometry_score": "structure_geometry_score",
        "rank_epistasis_score": "epistasis_score",
        "rank_mutation_burden_score": "mutation_burden_score",
    }
    for rcol, scol in rank_cols.items():
        out[rcol] = out[scol].rank(method="average", ascending=False)
    out["rank_agg_score"] = out[list(rank_cols)].mean(axis=1)
    return out


def select_portfolio(df: pd.DataFrame, cfg: dict[str, Any]) -> pd.DataFrame:
    k = int(cfg.get("v4", {}).get("final_k", 6))
    pool = df[df["v4_final_eligible"].astype(bool)].sort_values("v4_objective", ascending=False).copy()
    if pool.empty:
        return pool.head(0).copy()
    max_site_reuse = int(cfg.get("v4", {}).get("portfolio", {}).get("max_shared_mutations_per_site", 2))
    selected_idx: list[int] = []
    site_counts: Counter[int] = Counter()
    source_counts: Counter[str] = Counter()
    layers: dict[int, str] = {}
    # Greedy portfolio selection keeps high objective while preventing one hot mutation from dominating all final slots.
    for idx, row in pool.iterrows():
        muts = parse_mutation_text(row.get("mutations", ""))
        positions = [mutation_position(m) for m in muts if mutation_position(m) is not None]
        if any(site_counts[p] >= max_site_reuse for p in positions):
            continue
        selected_idx.append(idx)
        for p in positions:
            site_counts[p] += 1
        source_counts[str(row.get("source_branch", "unknown"))] += 1
        layers[idx] = infer_selection_layer(row)
        if len(selected_idx) >= k:
            break
    if len(selected_idx) < k:
        for idx, row in pool.iterrows():
            if idx not in selected_idx:
                selected_idx.append(idx)
                layers[idx] = infer_selection_layer(row) + "_fallback"
            if len(selected_idx) >= k:
                break
    out = pool.loc[selected_idx[:k]].copy()
    out["selection_layer"] = [layers.get(i, "v4_objective") for i in out.index]
    out["Seq_ID"] = [f"Seq_{i+1}" for i in range(len(out))]
    return out.reset_index(drop=True)


def infer_selection_layer(row: pd.Series) -> str:
    if bool(row.get("has_real_structure", False)) and bool(row.get("structure_pass", False)):
        return "structure_verified"
    if float(row.get("mutation_count", 0)) <= 4 and float(row.get("brightness_lcb_norm", 0)) >= 0.70:
        return "conservative_high_lcb"
    if float(row.get("tgp_surface_prior", 0)) >= 0.75:
        return "tgp_surface_stabilized"
    if float(row.get("stability_score", 0)) >= 0.65:
        return "balanced_stability"
    return "v4_objective"


def export_structure_tasks(priority: pd.DataFrame, cfg: dict[str, Any], out_dir: Path) -> dict[str, str]:
    outputs: dict[str, str] = {}
    af_n = int(cfg.get("af_branch", {}).get("export_top_n", cfg.get("v4", {}).get("structure_top_n", 200)))
    pm_n = int(cfg.get("proteinmpnn", {}).get("export_top_n", cfg.get("v4", {}).get("structure_top_n", 200)))
    af = priority.head(af_n)
    pm = priority.head(pm_n)
    af_fasta = out_dir / "af2_af3_structure_tasks.fasta"
    af_jsonl = out_dir / "af2_af3_structure_tasks.jsonl"
    pm_fasta = out_dir / "proteinmpnn_design_tasks.fasta"
    pm_jsonl = out_dir / "proteinmpnn_design_tasks.jsonl"
    with af_fasta.open("w", encoding="utf-8") as ffa, af_jsonl.open("w", encoding="utf-8") as fjs:
        for _, row in af.iterrows():
            sid = row.get("candidate_id") or row.get("Seq_ID")
            ffa.write(f">{sid}|mutations={row.get('mutations','')}|v4_objective={row.get('v4_objective',np.nan):.6f}\n{row['sequence']}\n")
            fjs.write(json.dumps({"id": sid, "sequence": row["sequence"], "mutations": row.get("mutations", ""), "backend": cfg.get("af_branch", {}).get("backend", "af2_or_af3")}, ensure_ascii=False) + "\n")
    with pm_fasta.open("w", encoding="utf-8") as ffa, pm_jsonl.open("w", encoding="utf-8") as fjs:
        for _, row in pm.iterrows():
            sid = row.get("candidate_id") or row.get("Seq_ID")
            ffa.write(f">{sid}|mutations={row.get('mutations','')}\n{row['sequence']}\n")
            # ProteinMPNN receives candidate backbones after AF/PDB generation; this task file binds sequence IDs to downstream backbone names.
            fjs.write(json.dumps({"id": sid, "sequence": row["sequence"], "mutations": row.get("mutations", ""), "expected_backbone_pdb": f"{sid}.pdb"}, ensure_ascii=False) + "\n")
    outputs.update({"af_fasta": str(af_fasta), "af_jsonl": str(af_jsonl), "proteinmpnn_fasta": str(pm_fasta), "proteinmpnn_jsonl": str(pm_jsonl)})
    return outputs


def write_v4_reports(
    selected: pd.DataFrame,
    ranked: pd.DataFrame,
    priority: pd.DataFrame,
    gate_top: pd.DataFrame,
    diagnostics: pd.DataFrame,
    metrics: dict[str, Any],
    cfg: dict[str, Any],
    out_dir: Path,
) -> None:
    with (out_dir / "training_metrics_v4.json").open("w", encoding="utf-8") as f:
        json.dump(ensure_jsonable(metrics), f, indent=2, ensure_ascii=False)
    final_empty_reason = ""
    if selected.empty:
        if int(diagnostics.loc[0, "has_real_structure_pass"]) == 0 and bool(cfg.get("v4", {}).get("require_real_structure_for_final", True)):
            final_empty_reason = "Final Top-K is intentionally empty because no real AF/PDB structure metrics were imported."
        elif int(diagnostics.loc[0, "v4_final_eligible_pass"]) == 0:
            final_empty_reason = "Final Top-K is empty because no candidate passed both V4 basic and structure gates."
        else:
            final_empty_reason = "Final Top-K is empty after portfolio constraints."
    lines = [
        "# SynBio GFP V4 Complete Pipeline Report",
        "",
        "## Executive summary",
        "",
        f"- total_candidates: {len(ranked)}",
        f"- v4_basic_gate_pass: {int(diagnostics.loc[0, 'v4_basic_gate_pass'])}",
        f"- v4_structure_gate_pass: {int(diagnostics.loc[0, 'v4_structure_gate_pass'])}",
        f"- v4_final_eligible_pass: {int(diagnostics.loc[0, 'v4_final_eligible_pass'])}",
        f"- final_selected: {len(selected)}",
        f"- final_empty_reason: {final_empty_reason or 'NA'}",
        "",
        "## Model metrics inherited from V3 training stage",
        "",
        f"- brightness_r2: {metrics.get('brightness_r2')}",
        f"- low_brightness_auc: {metrics.get('low_brightness_auc')}",
        f"- n_train: {metrics.get('n_train')}",
        f"- n_val: {metrics.get('n_val')}",
        "",
        "## V4 design logic",
        "",
        "- LCB = brightness_mean - lambda_std * brightness_std.",
        "- Basic gate blocks protected-site, chromophore-pocket, ion-network, high-risk, and excessive-burden candidates before AF/ProteinMPNN spending.",
        "- AF2/AF3 branch exports the basic-gated Top-N candidates for real structure prediction and can import PDB-derived RMSD metrics or a structure metrics CSV.",
        "- ProteinMPNN branch exports the same structure-priority pool and can import sequence/backbone compatibility scores.",
        "- Final Top-K is generated only from candidates that pass both basic and structure gates unless --allow-proxy-final is used for debugging.",
        "",
        "## Top structure-priority candidates",
        "",
    ]
    show_cols = ["candidate_id", "Seq_ID", "source_branch", "mutations", "v4_objective", "brightness_lcb", "low_brightness_risk", "structure_status", "v4_basic_gate", "v4_final_eligible"]
    for _, row in priority.head(12)[[c for c in show_cols if c in priority.columns]].iterrows():
        lines.append(f"- {row.to_dict()}")
    lines.extend(["", "## Final selected candidates", ""])
    if selected.empty:
        lines.append("No final candidates were emitted in this run.")
    else:
        for _, row in selected.iterrows():
            lines.extend([
                f"### {row.get('Seq_ID')}",
                f"- candidate_id: {row.get('candidate_id')}",
                f"- selection_layer: {row.get('selection_layer')}",
                f"- mutations: {row.get('mutations')}",
                f"- source_branch: {row.get('source_branch')}",
                f"- v4_objective: {row.get('v4_objective')}",
                f"- brightness_lcb: {row.get('brightness_lcb')}",
                f"- low_brightness_risk: {row.get('low_brightness_risk')}",
                f"- structure_status: {row.get('structure_status')}",
                f"- barrel_rmsd / pocket_rmsd / ion_network_rmsd: {row.get('barrel_rmsd')} / {row.get('pocket_rmsd')} / {row.get('ion_network_rmsd')}",
                "",
            ])
    (out_dir / "v4_pipeline_report.md").write_text("\n".join(lines), encoding="utf-8")


def write_submission(selected: pd.DataFrame, team_name: str, out_dir: Path) -> None:
    sub = pd.DataFrame(columns=["Team_Name", "Seq_ID", "Sequence"])
    if not selected.empty:
        sub = pd.DataFrame({"Team_Name": [team_name] * len(selected), "Seq_ID": selected["Seq_ID"].tolist(), "Sequence": selected["sequence"].tolist()})
    sub.to_csv(out_dir / "submission_v4.csv", index=False)


def diagnostics_frame(df: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame([{
        "total_candidates": int(len(df)),
        "v4_basic_gate_pass": int(df.get("v4_basic_gate", pd.Series(False, index=df.index)).sum()),
        "v4_structure_gate_pass": int(df.get("v4_structure_gate", pd.Series(False, index=df.index)).sum()),
        "v4_final_eligible_pass": int(df.get("v4_final_eligible", pd.Series(False, index=df.index)).sum()),
        "has_real_structure_pass": int(df.get("has_real_structure", pd.Series(False, index=df.index)).sum()),
        "structure_pass_true": int(df.get("structure_pass", pd.Series(False, index=df.index)).sum()),
        "protected_mutation_count_zero": int((safe_num(df, "protected_mutation_count") == 0).sum()),
        "pocket_mutation_count_zero": int((safe_num(df, "pocket_mutation_count") == 0).sum()),
        "ion_network_mutation_count_zero": int((safe_num(df, "ion_network_mutation_count") == 0).sum()),
        "low_brightness_risk_le_035": int((safe_num(df, "low_brightness_risk") <= 0.35).sum()),
        "mutation_count_le_9": int((safe_num(df, "mutation_count") <= 9).sum()),
    }])


def mutation_frequency(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, row in df.iterrows():
        seed = row.get("seed", "unknown")
        source = row.get("source_branch", row.get("source", "unknown"))
        for m in parse_mutation_text(row.get("mutations", "")):
            rows.append({"seed": seed, "source_branch": source, "mutation": m})
    if not rows:
        return pd.DataFrame(columns=["seed", "source_branch", "mutation", "count"])
    return pd.DataFrame(rows).value_counts(["seed", "source_branch", "mutation"]).reset_index(name="count").sort_values("count", ascending=False)


def run_pipeline_v4(cfg: dict[str, Any]) -> dict[str, str]:
    out_dir = Path(cfg["paths"]["out_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)
    logger = setup_logger(out_dir, name="synbio_gfp_v4")
    set_seed(int(cfg.get("seed", 42)))
    start_time = time.time()
    (out_dir / "config_resolved_v4.json").write_text(json.dumps(ensure_jsonable(cfg), indent=2, ensure_ascii=False), encoding="utf-8")
    logger.info("V4 COMPLETE PIPELINE START")

    data_dir = Path(cfg["paths"]["data_dir"])
    refs = parse_fasta_like(_resolve_file(data_dir, cfg["paths"].get("reference_sequences_file", "AAseqs of 5 GFP proteins.txt")))
    scaffold_key = next((k for k in refs if "sfgfp" in k.lower() or "sf" in k.lower()), None) or sorted(refs)[0]
    scaffold = refs[scaffold_key]
    exclusion = load_exclusion(_resolve_file(data_dir, cfg["paths"].get("exclusion_list_file", "Exclusion_List.csv")))
    policy = build_site_policy(scaffold, cfg.get("site_policy", {}).get("extra_protected_positions"))
    logger.info("REFERENCE LOADED | scaffold=%s | length=%d", scaffold_key, len(scaffold))

    train_df, meta = load_training_table(
        _resolve_file(data_dir, cfg["paths"].get("gfp_data_file", "GFP_data.xlsx")),
        scaffold=scaffold,
        max_rows=cfg.get("training", {}).get("max_train_samples"),
        seed=int(cfg.get("seed", 42)),
    )
    brightness_col = meta["brightness_col"]
    cache = ModelCache(cfg["paths"].get("cache_dir", "cache"), logger=logger)
    feature_cfg = dict(cfg.get("features", {}))
    X_train, train_feature_key = cache.get_or_compute_features(train_df["full_sequence"].tolist(), feature_cfg, lambda: featurize(train_df["full_sequence"].tolist(), feature_cfg, logger))
    model_key = cache.model_key(train_feature_key, cfg.get("training", {}), meta)
    bundle = cache.load_model_bundle(model_key)
    if bundle is None:
        bundle = train_predictive_models(X_train, train_df[brightness_col].astype(float).to_numpy(), cfg, meta, logger)
        cache.save_model_bundle(model_key, bundle)
        cache.save_metrics(model_key, bundle["metrics"])
    metrics = bundle.get("metrics") or cache.load_metrics(model_key) or {}

    cand = generate_candidates(scaffold, policy, cfg)
    explain = [policy.explain_mutations(scaffold, s) for s in cand["sequence"]]
    cand = pd.concat([cand.reset_index(drop=True), pd.DataFrame(explain)], axis=1)
    cand = add_literature_priors(cand, scaffold, policy)
    cand = add_proxy_scores(cand)
    X_cand, _ = cache.get_or_compute_features(cand["sequence"].tolist(), feature_cfg, lambda: featurize(cand["sequence"].tolist(), feature_cfg, logger))
    pred = predict_bundle(bundle, X_cand)
    for k, v in pred.items():
        cand[k] = v
    cand = add_objective_scores(cand)
    cand = normalize_v3_columns(cand, scaffold)
    cand = add_site_counts(cand)
    cand = add_v4_lcb_and_risk(cand, cfg)
    cand = add_epistasis_and_ddg_proxies(cand)
    cand = merge_proteinmpnn_scores(cand, cfg)
    cand = add_structure_branch(cand, cfg, out_dir, logger)
    cand = add_v4_gates(cand, cfg)
    cand = add_v4_objective(cand)

    cand["in_exclusion_list"] = cand["sequence"].isin(exclusion)
    cand.loc[cand["in_exclusion_list"], "v4_basic_gate"] = False
    cand.loc[cand["in_exclusion_list"], "v4_final_eligible"] = False
    cand = add_v4_objective(cand)
    ranked = cand.sort_values("v4_objective", ascending=False).reset_index(drop=True)
    ranked.to_csv(out_dir / "all_ranked_candidates_v4.csv", index=False)

    top_n = int(cfg.get("v4", {}).get("structure_top_n", cfg.get("structure", {}).get("top_n", 200)))
    priority = ranked[ranked["v4_basic_gate"].astype(bool)].sort_values("v4_objective", ascending=False).head(top_n).copy()
    priority.to_csv(out_dir / "structure_priority_top200_v4.csv", index=False)
    gate_top = ranked.sort_values("v4_objective", ascending=False).head(top_n).copy()
    gate_top.to_csv(out_dir / "ranked_top200_v4.csv", index=False)
    # Backward-compatible filename for earlier V4 analysis notebooks.
    gate_top.to_csv(out_dir / "structure_gate_top200_v4.csv", index=False)

    selected = select_portfolio(ranked, cfg)
    selected.to_csv(out_dir / "final_top6_v4.csv", index=False)
    write_submission(selected, cfg.get("team_name", "YourTeamName"), out_dir)

    diagnostics = diagnostics_frame(ranked)
    diagnostics.to_csv(out_dir / "v4_gate_diagnostics.csv", index=False)
    mutation_frequency(ranked).to_csv(out_dir / "mutation_frequency_by_seed_v6.csv", index=False)
    ranked[~ranked["v4_basic_gate"].astype(bool)].head(2000).to_csv(out_dir / "rejected_by_basic_gate.csv", index=False)
    ranked[ranked["v4_basic_gate"].astype(bool) & ~ranked["v4_structure_gate"].astype(bool)].head(2000).to_csv(out_dir / "rejected_by_structure_gate.csv", index=False)
    # Backward-compatible explorer diagnostic file; it is now populated only when explorer-like sources exist.
    explorer_rej = ranked[ranked["source_branch"].astype(str).str.contains("explor", case=False, na=False) & ~ranked["v4_basic_gate"].astype(bool)].copy()
    explorer_rej.head(2000).to_csv(out_dir / "rejected_explorer_candidates.csv", index=False)
    export_paths = export_structure_tasks(priority, cfg, out_dir)
    write_v4_reports(selected, ranked, priority, gate_top, diagnostics, metrics, cfg, out_dir)
    logger.info("V4 COMPLETE PIPELINE DONE | seconds=%.2f | final=%d", time.time() - start_time, len(selected))
    return {
        "ranked_candidates": str(out_dir / "all_ranked_candidates_v4.csv"),
        "structure_priority": str(out_dir / "structure_priority_top200_v4.csv"),
        "final_top6": str(out_dir / "final_top6_v4.csv"),
        "submission": str(out_dir / "submission_v4.csv"),
        "diagnostics": str(out_dir / "v4_gate_diagnostics.csv"),
        "report": str(out_dir / "v4_pipeline_report.md"),
        **export_paths,
    }


def main() -> None:
    args = parse_args()
    cfg = resolve_v4_config(load_config(args.config), args)
    outputs = run_pipeline_v4(cfg)
    print("\nV4 complete outputs:")
    for k, v in outputs.items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
