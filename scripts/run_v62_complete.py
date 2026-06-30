from __future__ import annotations

import argparse
import json
import hashlib
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

from synbio_gfp_v6.candidates import generate_candidates
from synbio_gfp_v6.config import load_config, resolve_config
from synbio_gfp_v6.data import load_exclusion, load_training_table, parse_fasta_like
from synbio_gfp_v6.features import featurize
from synbio_gfp_v6.literature_priors import add_literature_priors
from synbio_gfp_v6.model_cache import ModelCache
from synbio_gfp_v6.models import predict_bundle, train_predictive_models
from synbio_gfp_v6.scoring import add_objective_scores, add_proxy_scores, rank01
from synbio_gfp_v6.site_policy import (
    BARREL_CORE_PRIOR,
    CHROMOPHORE,
    CHROMOPHORE_POCKET_8A,
    ION_NETWORK,
    build_site_policy,
)
from synbio_gfp_v6.structure_filter import compute_pdb_metrics
from synbio_gfp_v6.feedback import add_v6_feedback_scores, build_feedback_policy_from_sources, write_feedback_reports
from synbio_gfp_v6.utils import ensure_jsonable, setup_logger, set_seed, validate_sequence

AA20 = set("ACDEFGHIKLMNPQRSTVWY")
MUT_RE = re.compile(r"([ACDEFGHIKLMNPQRSTVWY])(\d{1,3})([ACDEFGHIKLMNPQRSTVWY])")

def sequence_sha1(seq: str) -> str:
    return hashlib.sha1(str(seq).strip().upper().encode("utf-8")).hexdigest()


def _resolve_file(data_dir: Path, filename: str) -> Path:
    p = data_dir / filename
    if p.exists():
        return p
    candidates = list(data_dir.glob(filename))
    if candidates:
        return candidates[0]
    raise FileNotFoundError(f"Cannot find {filename} under {data_dir}")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Run complete SynBio GFP V6.2 Adaptive traceable pipeline")
    p.add_argument("--config", required=True)
    p.add_argument("--data-dir", required=True)
    p.add_argument("--team-name", required=True)
    p.add_argument("--feature-mode", choices=["simple", "esm"], default="simple")
    p.add_argument("--esm-model", default="esm2_t30_150M_UR50D")
    p.add_argument("--max-train-samples", type=int, default=20000)
    p.add_argument("--n-candidates", type=int, default=20000)
    p.add_argument("--out-dir", default="outputs/run_v6_complete")
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
    p.add_argument("--previous-ranked-csv", default=None, help="Optional V5/V5.1/V6 ranked candidate CSV used to learn closed-loop mutation feedback before generating candidates.")
    p.add_argument("--feedback-csv", default=None, help="Optional feedback table CSV with candidate-level structure/ProteinMPNN outcomes.")
    p.add_argument("--brightness-teacher-csv", default=None, help="Optional V5.1/V6 ranked CSV used as an exact-sequence brightness safety teacher.")
    return p.parse_args()


def resolve_v6_config(raw: dict[str, Any], args: argparse.Namespace) -> dict[str, Any]:
    cfg = resolve_config(raw, args)
    cfg.setdefault("v5", {})
    cfg.setdefault("structure", {})
    cfg.setdefault("proteinmpnn", {})
    cfg.setdefault("af_branch", {})
    cfg.setdefault("v6", {})
    cfg.setdefault("v6", {}).setdefault("feedback", {})
    if args.previous_ranked_csv:
        cfg["v6"]["feedback"]["previous_ranked_csv"] = args.previous_ranked_csv
    if args.feedback_csv:
        cfg["v6"]["feedback"]["feedback_csv"] = args.feedback_csv
    if getattr(args, "brightness_teacher_csv", None):
        cfg.setdefault("v62", {}).setdefault("brightness_teacher", {})["previous_ranked_csv"] = args.brightness_teacher_csv
    if args.final_k is not None:
        cfg["v5"]["final_k"] = args.final_k
    if args.lambda_std is not None:
        cfg["v5"]["lambda_std"] = args.lambda_std
    if args.allow_proxy_final:
        cfg["v5"]["require_real_structure_for_final"] = False
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


def mutation_position_set(text: Any) -> set[int]:
    return {p for p in (mutation_position(m) for m in parse_mutation_text(text)) if p is not None}


def sequence_hamming(a: str, b: str) -> int:
    if len(a) != len(b):
        return max(len(a), len(b))
    return sum(x != y for x, y in zip(a, b))


def shared_mutation_sites(a: Any, b: Any) -> int:
    return len(mutation_position_set(a) & mutation_position_set(b))


def infer_v6_bucket(row: pd.Series) -> str:
    source = str(row.get("source_branch", row.get("source", ""))).lower()
    mut_count = float(row.get("mutation_count", 0) or 0)
    pmpnn_norm = float(row.get("proteinmpnn_score_norm", 0.5) or 0.5)
    feedback_score = float(row.get("feedback_policy_score", 0.5) or 0.5)
    if source.startswith("measured_feedback"):
        return "measured_feedback_repair"
    if source.startswith("feedback") or feedback_score >= 0.68:
        return "measured_feedback_repair"
    if mut_count <= 4 and float(row.get("low_brightness_risk", 1) or 1) <= 0.10:
        return "low_brightness_safe_conservative"
    if mut_count <= 4 and float(row.get("brightness_lcb_norm", 0) or 0) >= 0.68 and float(row.get("low_brightness_risk", 1) or 1) <= 0.25:
        return "brightness_lcb_exploitation"
    if "loop" in source or float(row.get("loop_rigidity_score", 0.0) or 0.0) >= 0.66:
        return "loop_rigidity"
    if bool(row.get("has_real_structure", False)) and float(row.get("structure_geometry_score", 0) or 0) >= 0.45:
        return "structure_verified_safe"
    if "surface_patch" in source:
        return "surface_patch_thermal"
    if "tgp" in source or float(row.get("tgp_surface_prior", 0) or 0) >= 0.72:
        return "surface_patch_thermal"
    if pmpnn_norm >= 0.70:
        return "pmpnn_compatible"
    if mut_count >= 6:
        return "diversity_reserve"
    return "balanced_stability"



def add_v62_brightness_teacher(df: pd.DataFrame, cfg: dict[str, Any], logger: Any | None = None) -> pd.DataFrame:
    out = df.copy()
    teacher_cfg = cfg.get("v62", {}).get("brightness_teacher", {})
    teacher_path = teacher_cfg.get("previous_ranked_csv") or cfg.get("v6", {}).get("feedback", {}).get("previous_ranked_csv")
    out["brightness_teacher_coverage"] = False
    out["brightness_lcb_current"] = safe_num(out, "brightness_lcb", default=0.0)
    out["low_brightness_risk_current"] = safe_num(out, "low_brightness_risk", default=1.0)
    if not teacher_path or not Path(str(teacher_path)).exists():
        if logger:
            logger.info("V6.2 BRIGHTNESS TEACHER | no exact-sequence teacher CSV provided; using current model only")
        out["brightness_lcb_ensemble"] = out["brightness_lcb_current"]
        out["low_brightness_risk_ensemble"] = out["low_brightness_risk_current"]
        return out
    teacher = load_optional_score_table(str(teacher_path))
    if teacher is None or teacher.empty or "sequence" not in teacher.columns:
        out["brightness_lcb_ensemble"] = out["brightness_lcb_current"]
        out["low_brightness_risk_ensemble"] = out["low_brightness_risk_current"]
        return out
    teacher = teacher.copy()
    teacher["_seq_key"] = teacher["sequence"].astype(str).str.strip().str.upper()
    t_lcb_col = next((c for c in ["brightness_lcb", "brightness_lcb_ensemble", "brightness_pred", "brightness_mean"] if c in teacher.columns), None)
    t_risk_col = next((c for c in ["low_brightness_risk", "low_brightness_risk_ensemble"] if c in teacher.columns), None)
    keep = ["_seq_key"] + ([t_lcb_col] if t_lcb_col else []) + ([t_risk_col] if t_risk_col else [])
    teacher = teacher[keep].drop_duplicates("_seq_key")
    rename = {}
    if t_lcb_col:
        rename[t_lcb_col] = "teacher_brightness_lcb"
    if t_risk_col:
        rename[t_risk_col] = "teacher_low_brightness_risk"
    teacher = teacher.rename(columns=rename)
    out["_seq_key"] = out["sequence"].astype(str).str.strip().str.upper()
    out = out.merge(teacher, on="_seq_key", how="left").drop(columns=["_seq_key"], errors="ignore")
    out["brightness_teacher_coverage"] = out.get("teacher_brightness_lcb", pd.Series(np.nan, index=out.index)).notna() | out.get("teacher_low_brightness_risk", pd.Series(np.nan, index=out.index)).notna()
    w_teacher = float(teacher_cfg.get("teacher_weight", 0.50))
    w_current = float(teacher_cfg.get("current_weight", 0.50))
    denom = max(1e-8, w_teacher + w_current)
    teacher_lcb = pd.to_numeric(out.get("teacher_brightness_lcb", pd.Series(np.nan, index=out.index)), errors="coerce")
    current_lcb = out["brightness_lcb_current"]
    out["brightness_lcb_ensemble"] = np.where(teacher_lcb.notna(), (w_teacher * teacher_lcb.fillna(0.0) + w_current * current_lcb) / denom, current_lcb)
    teacher_risk = pd.to_numeric(out.get("teacher_low_brightness_risk", pd.Series(np.nan, index=out.index)), errors="coerce")
    current_risk = out["low_brightness_risk_current"]
    out["low_brightness_risk_ensemble"] = np.where(teacher_risk.notna(), (w_teacher * teacher_risk.fillna(1.0) + w_current * current_risk) / denom, current_risk)
    # The ensemble safety columns become production columns so all downstream gates consume the conservative teacher when available.
    out["brightness_lcb"] = pd.to_numeric(out["brightness_lcb_ensemble"], errors="coerce").fillna(current_lcb)
    out["brightness_lcb_norm"] = minmax01(out["brightness_lcb"], higher_is_better=True, default=0.5)
    out["low_brightness_risk"] = pd.to_numeric(out["low_brightness_risk_ensemble"], errors="coerce").fillna(current_risk).clip(0.0, 1.0)
    out["low_risk_score"] = (1.0 - out["low_brightness_risk"]).clip(0.0, 1.0)
    if logger:
        logger.info("V6.2 BRIGHTNESS TEACHER | file=%s | exact_sequence_coverage=%d/%d", teacher_path, int(out["brightness_teacher_coverage"].sum()), len(out))
    return out


def add_v62_physical_thermal_proxy(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    muts = out.get("mutations", pd.Series("", index=out.index)).astype(str)
    acid_gain = []
    hydrophobic_alt = []
    exposed_cys_or_aromatic = []
    loop_like = []
    for text in muts:
        parsed = parse_mutation_text(text)
        denom = max(1, len(parsed))
        acid_gain.append(sum(1 for m in parsed if mutation_alt(m) in {"D", "E"}) / denom)
        hydrophobic_alt.append(sum(1 for m in parsed if mutation_alt(m) in {"L", "I", "V", "F", "W", "M"}) / denom)
        exposed_cys_or_aromatic.append(sum(1 for m in parsed if mutation_alt(m) in {"C", "W", "F", "Y"}) / denom)
        # Mild polar/turn-compatible substitutions get a small rigidity prior without globally rewarding Pro/Gly.
        loop_like.append(sum(1 for m in parsed if mutation_alt(m) in {"S", "T", "N", "D", "E", "A"}) / denom)
    out["surface_charge_score"] = np.clip(0.50 * safe_num(out, "tgp_surface_prior", 0.0) + 0.50 * pd.Series(acid_gain, index=out.index), 0.0, 1.0)
    out["hydrophobic_patch_penalty"] = np.clip(0.65 * pd.Series(hydrophobic_alt, index=out.index) + 0.35 * pd.Series(exposed_cys_or_aromatic, index=out.index), 0.0, 1.0)
    out["aggregation_safety_score"] = (1.0 - out["hydrophobic_patch_penalty"]).clip(0.0, 1.0)
    out["salt_bridge_potential"] = np.clip(0.50 + 0.35 * out["surface_charge_score"] - 0.20 * out["hydrophobic_patch_penalty"], 0.0, 1.0)
    out["loop_rigidity_score"] = np.clip(0.35 + 0.45 * pd.Series(loop_like, index=out.index) - 0.10 * np.maximum(0.0, safe_num(out, "mutation_count") - 5), 0.0, 1.0)
    structure_component = (0.70 * safe_num(out, "structure_geometry_score", 0.05) + 0.30 * safe_num(out, "proteinmpnn_score_norm", 0.5)).clip(0.0, 1.0)
    out["thermal_proxy_v62"] = (
        0.20 * structure_component
        + 0.18 * out["surface_charge_score"]
        + 0.15 * out["salt_bridge_potential"]
        + 0.15 * out["aggregation_safety_score"]
        + 0.12 * out["loop_rigidity_score"]
        + 0.10 * safe_num(out, "proteinmpnn_score_norm", 0.5)
        + 0.10 * safe_num(out, "foldability_proxy", 0.5)
    ).clip(0.0, 1.0)
    return out


def add_v62_objective(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["v62_structure_score"] = (0.60 * safe_num(out, "structure_geometry_score", 0.05) + 0.40 * safe_num(out, "proteinmpnn_score_norm", 0.5)).clip(0.0, 1.0)
    out["v62_feedback_score"] = safe_num(out, "feedback_policy_score", 0.5).clip(0.0, 1.0)
    out["v62_blacklist_penalty"] = np.minimum(safe_num(out, "motif_blacklist_hits", 0.0), 4.0) / 4.0
    mut_count = safe_num(out, "mutation_count", 0.0)
    out["v62_mutation_burden_penalty"] = np.where(mut_count <= 5, 0.0, np.minimum((mut_count - 5.0) / 4.0, 1.0))
    out["v62_objective"] = (
        0.30 * out["brightness_lcb_norm"]
        + 0.20 * out["low_risk_score"]
        + 0.18 * safe_num(out, "thermal_proxy_v62", 0.5)
        + 0.12 * out["v62_structure_score"]
        + 0.08 * safe_num(out, "proteinmpnn_score_norm", 0.5)
        + 0.07 * out["v62_feedback_score"]
        + 0.05 * safe_num(out, "literature_score", 0.5)
        - 0.06 * safe_num(out, "epistasis_penalty", 0.0)
        - 0.05 * out["v62_mutation_burden_penalty"]
        - 0.04 * out["v62_blacklist_penalty"]
    ).clip(0.0, 1.0)
    # Compatibility aliases allow existing V6 notebooks and deploy scripts to keep working while ranking by V6.2 objective.
    out["v6_objective"] = out["v62_objective"]
    out["v5_objective"] = out["v62_objective"]
    return out


def add_v6_attribution(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    # V6 attribution keeps the single-ESM-style baseline visible so report gains cannot be mistaken for real algorithmic uplift.
    out["esm_baseline_score"] = safe_num(out, "objective", default=0.0)
    out["esm_baseline_rank"] = out["esm_baseline_score"].rank(method="first", ascending=False).astype(int)
    out["v6_rank"] = safe_num(out, "v5_objective", default=0.0).rank(method="first", ascending=False).astype(int)
    out["rank_delta_vs_esm"] = out["esm_baseline_rank"] - out["v6_rank"]
    out["is_in_esm_top6"] = out["esm_baseline_rank"] <= 6
    out["is_in_esm_top20"] = out["esm_baseline_rank"] <= 20
    out["structure_score_valid"] = out.get("has_real_structure", pd.Series(False, index=out.index)).astype(bool)
    out["proteinmpnn_score_valid"] = pd.to_numeric(out.get("proteinmpnn_score_raw", pd.Series(np.nan, index=out.index)), errors="coerce").notna()
    out["v6_bucket"] = out.apply(infer_v6_bucket, axis=1)
    return out


def greedy_diverse_select(pool: pd.DataFrame, k: int, cfg: dict[str, Any], *, min_hamming: int, max_shared_sites_pair: int, max_site_reuse: int) -> list[int]:
    selected: list[int] = []
    selected_sequences: list[str] = []
    selected_mutations: list[Any] = []
    site_counts: Counter[int] = Counter()
    for idx, row in pool.sort_values("v5_objective", ascending=False).iterrows():
        seq = str(row.get("sequence", ""))
        muts = row.get("mutations", "")
        positions = mutation_position_set(muts)
        # The hard diversity checks are intentionally local and auditable: sequence distance, shared mutation sites, and global site reuse.
        if selected_sequences and min(sequence_hamming(seq, s) for s in selected_sequences) < min_hamming:
            continue
        if selected_mutations and max(shared_mutation_sites(muts, m) for m in selected_mutations) > max_shared_sites_pair:
            continue
        if any(site_counts[p0] >= max_site_reuse for p0 in positions):
            continue
        selected.append(idx)
        selected_sequences.append(seq)
        selected_mutations.append(muts)
        for p0 in positions:
            site_counts[p0] += 1
        if len(selected) >= k:
            break
    return selected


def select_structure_priority_v6(ranked: pd.DataFrame, cfg: dict[str, Any]) -> pd.DataFrame:
    v6_cfg = cfg.get("v6", {})
    top_n = int(cfg.get("v5", {}).get("structure_top_n", cfg.get("structure", {}).get("top_n", 200)))
    pool = ranked[ranked["v5_basic_gate"].astype(bool)].copy()
    if pool.empty:
        return pool.head(0).copy()
    quotas = v6_cfg.get("structure_portfolio", {}).get("quotas", {
        "feedback_guided_exploitation": 60,
        "brightness_lcb_exploitation": 50,
        "surface_charge_exploration": 35,
        "structure_verified_safe": 30,
        "pmpnn_compatible": 15,
        "diversity_reserve": 10,
    })
    min_hamming = int(v6_cfg.get("structure_portfolio", {}).get("min_pairwise_hamming", 1))
    max_shared_pair = int(v6_cfg.get("structure_portfolio", {}).get("max_shared_mutations_per_pair", 5))
    max_site_reuse = int(v6_cfg.get("structure_portfolio", {}).get("max_shared_mutations_per_site", 25))
    chosen: list[int] = []
    for bucket, quota in quotas.items():
        sub = pool[pool["v6_bucket"].astype(str) == str(bucket)]
        chosen.extend(i for i in greedy_diverse_select(sub, int(quota), cfg, min_hamming=min_hamming, max_shared_sites_pair=max_shared_pair, max_site_reuse=max_site_reuse) if i not in chosen)
        if len(chosen) >= top_n:
            break
    if len(chosen) < top_n:
        fallback = pool.drop(index=chosen, errors="ignore")
        chosen.extend(i for i in greedy_diverse_select(fallback, top_n - len(chosen), cfg, min_hamming=0, max_shared_sites_pair=999, max_site_reuse=999) if i not in chosen)
    out = pool.loc[chosen[:top_n]].copy()
    out["structure_priority_rank_v6"] = range(1, len(out) + 1)
    return out.reset_index(drop=True)


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


def add_v5_lcb_and_risk(df: pd.DataFrame, cfg: dict[str, Any]) -> pd.DataFrame:
    out = df.copy()
    lam = float(cfg.get("v5", {}).get("lambda_std", 0.75))
    out["brightness_mean"] = safe_num(out, "brightness_pred", default=0.0)
    out["brightness_std"] = safe_num(out, "brightness_uncertainty", default=0.0)
    # LCB is the conservative brightness objective used by V5 to avoid uncertainty-chasing candidates.
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
    out["ddg_total_v5"] = out["ddg_effective"] + 0.10 * safe_num(out, "low_brightness_risk")
    out["p_fluorescent_v5"] = 1.0 / (1.0 + np.exp(1.15 * (out["ddg_total_v5"] - 1.45)))
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
        logger.info("V5 STRUCTURE METRICS IMPORTED | file=%s | real=%d", metrics_csv, int(out["has_real_structure"].sum()))
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
        logger.info("V5 PDB STRUCTURE CHECK DONE | real=%d", len(records))
    # Real structure is required for final eligibility unless explicitly disabled; proxy score remains diagnostic only.
    out["structure_geometry_score"] = structure_geometry_from_metrics(out)
    no_real = ~out["has_real_structure"].astype(bool)
    out.loc[no_real, "structure_geometry_score"] = 0.05
    out.loc[no_real & out["structure_status"].isna(), "structure_status"] = "not_checked"
    return out


def add_v5_gates(df: pd.DataFrame, cfg: dict[str, Any]) -> pd.DataFrame:
    out = df.copy()
    g = cfg.get("v5", {}).get("basic_gate", {})
    sg = cfg.get("v5", {}).get("structure_gate", {})
    sel = cfg.get("selection", {})
    min_len = int(sel.get("min_len", 220))
    max_len = int(sel.get("max_len", 250))
    valid_rows = []
    for seq in out["sequence"].astype(str):
        v = validate_sequence(seq, min_len, max_len)
        valid_rows.append(all(v.values()))
    out["contest_sequence_valid"] = valid_rows
    # Basic gate is deliberately independent of AF/ProteinMPNN availability so it can define the structure-priority pool.
    out["v5_basic_gate"] = (
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
    require_real = bool(cfg.get("v5", {}).get("require_real_structure_for_final", True))
    if require_real:
        out["v5_structure_gate"] = out["structure_pass"].astype(bool)
    else:
        # This debug path preserves production columns but lets users smoke-test final selection without AF/PDB results.
        out["v5_structure_gate"] = out["v5_basic_gate"].astype(bool)
    out["v5_final_eligible"] = out["v5_basic_gate"].astype(bool) & out["v5_structure_gate"].astype(bool)
    return out


def add_v5_objective(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["thermal_score"] = minmax01(safe_num(out, "tgp_surface_prior"), True, 0.5)
    out["stability_score"] = minmax01(safe_num(out, "stability_proxy"), True, 0.5)
    out["foldability_score"] = minmax01(safe_num(out, "foldability_proxy"), True, 0.5)
    out["literature_score"] = minmax01(safe_num(out, "literature_prior"), True, 0.5)
    out["mutation_burden_score"] = (1.0 - np.clip((safe_num(out, "mutation_count") - 2.0) / 10.0, 0.0, 1.0)).clip(0.0, 1.0)
    out["structure_score_v5"] = (0.75 * safe_num(out, "structure_geometry_score", 0.05) + 0.25 * safe_num(out, "proteinmpnn_score_norm", 0.5)).clip(0.0, 1.0)
    # V5 objective combines conservative brightness, low-risk probability, stability priors, structural plausibility, and epistasis control.
    out["v5_objective"] = (
        0.26 * out["brightness_lcb_norm"]
        + 0.17 * out["low_risk_score"]
        + 0.13 * out["thermal_score"]
        + 0.12 * out["stability_score"]
        + 0.08 * out["foldability_score"]
        + 0.10 * out["literature_score"]
        + 0.08 * out["structure_score_v5"]
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



def add_v6_objective(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    # V6 increases the causal impact of measured structural compatibility and feedback-derived mutation evidence.
    out["v6_structure_score"] = (0.60 * safe_num(out, "structure_geometry_score", 0.05) + 0.40 * safe_num(out, "proteinmpnn_score_norm", 0.5)).clip(0.0, 1.0)
    out["v6_feedback_score"] = safe_num(out, "feedback_policy_score", 0.5).clip(0.0, 1.0)
    out["v6_blacklist_penalty"] = np.minimum(safe_num(out, "motif_blacklist_hits", 0.0), 4.0) / 4.0
    out["v6_objective"] = (
        0.22 * out["brightness_lcb_norm"]
        + 0.14 * out["low_risk_score"]
        + 0.12 * out["thermal_score"]
        + 0.10 * out["stability_score"]
        + 0.06 * out["foldability_score"]
        + 0.07 * out["literature_score"]
        + 0.16 * out["v6_structure_score"]
        + 0.09 * out["v6_feedback_score"]
        + 0.03 * out["epistasis_score"]
        + 0.02 * out["mutation_burden_score"]
        - 0.04 * out["v6_blacklist_penalty"]
    ).clip(0.0, 1.0)
    # Keep the V5 column as a compatibility alias for older notebooks and deploy scripts.
    out["v5_objective"] = out["v6_objective"]
    return out

def select_portfolio(df: pd.DataFrame, cfg: dict[str, Any]) -> pd.DataFrame:
    k = int(cfg.get("v5", {}).get("final_k", 6))
    pool = df[df["v5_final_eligible"].astype(bool)].sort_values("v6_objective", ascending=False).copy()
    hard_cfg = cfg.get("v62", {}).get("final_constraints", {})
    max_final_risk = float(hard_cfg.get("max_final_low_brightness_risk", 0.20))
    normal_max_mut = int(hard_cfg.get("normal_max_mutation_count", 5))
    pool = pool[safe_num(pool, "low_brightness_risk", 1.0) <= max_final_risk].copy()
    if pool.empty:
        return pool.head(0).copy()
    final_cfg = cfg.get("v6", {}).get("final_portfolio", {})
    min_hamming = int(final_cfg.get("min_pairwise_hamming", 5))
    max_shared_pair = int(final_cfg.get("max_shared_mutations_per_pair", 3))
    max_site_reuse = int(final_cfg.get("max_shared_mutations_per_site", 2))
    max_per_source = int(final_cfg.get("max_per_source_branch", 2))
    bucket_quotas = final_cfg.get("mechanism_quotas", {
        "feedback_guided_exploitation": 1,
        "brightness_lcb_exploitation": 1,
        "surface_charge_exploration": 1,
        "structure_verified_safe": 1,
        "pmpnn_compatible": 1,
        "diversity_reserve": 1,
    })

    selected_idx: list[int] = []
    selected_sequences: list[str] = []
    selected_mutations: list[Any] = []
    site_counts: Counter[int] = Counter()
    source_counts: Counter[str] = Counter()

    def can_add(idx: int, row: pd.Series, *, relax: bool = False) -> bool:
        seq = str(row.get("sequence", ""))
        muts = row.get("mutations", "")
        positions = mutation_position_set(muts)
        source = str(row.get("source_branch", row.get("source", "unknown")))
        if idx in selected_idx:
            return False
        if not relax and float(row.get("mutation_count", 0) or 0) > normal_max_mut and sum(float(pool.loc[i].get("mutation_count", 0) or 0) > normal_max_mut for i in selected_idx) >= 1:
            return False
        if not relax and source_counts[source] >= max_per_source:
            return False
        if selected_sequences and min(sequence_hamming(seq, s) for s in selected_sequences) < (0 if relax else min_hamming):
            return False
        if selected_mutations and max(shared_mutation_sites(muts, m) for m in selected_mutations) > (999 if relax else max_shared_pair):
            return False
        if not relax and any(site_counts[p0] >= max_site_reuse for p0 in positions):
            return False
        return True

    def add_idx(idx: int, row: pd.Series) -> None:
        selected_idx.append(idx)
        selected_sequences.append(str(row.get("sequence", "")))
        muts = row.get("mutations", "")
        selected_mutations.append(muts)
        source_counts[str(row.get("source_branch", row.get("source", "unknown")))] += 1
        for p0 in mutation_position_set(muts):
            site_counts[p0] += 1

    # First pass enforces mechanism-level quotas so Top-6 is not merely a cluster of near-identical high-rank candidates.
    for bucket, quota in bucket_quotas.items():
        if len(selected_idx) >= k:
            break
        sub = pool[pool["v6_bucket"].astype(str) == str(bucket)].sort_values("v6_objective", ascending=False)
        taken = 0
        for idx, row in sub.iterrows():
            if taken >= int(quota) or len(selected_idx) >= k:
                break
            if can_add(idx, row):
                add_idx(idx, row)
                taken += 1

    # Second pass fills remaining slots with the best candidates under the same diversity constraints.
    if len(selected_idx) < k:
        for idx, row in pool.iterrows():
            if len(selected_idx) >= k:
                break
            if can_add(idx, row):
                add_idx(idx, row)

    # Final fallback is explicit and auditable; it only triggers when the eligible pool is too narrow.
    if len(selected_idx) < k:
        for idx, row in pool.iterrows():
            if len(selected_idx) >= k:
                break
            if can_add(idx, row, relax=True):
                add_idx(idx, row)

    out = pool.loc[selected_idx[:k]].copy()
    out["selection_layer"] = [infer_selection_layer(row) for _, row in out.iterrows()]
    if len(out) < k:
        out["selection_layer"] = out["selection_layer"].astype(str) + "_underfilled"
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
    return "v5_objective"


def export_structure_tasks(priority: pd.DataFrame, cfg: dict[str, Any], out_dir: Path) -> dict[str, str]:
    outputs: dict[str, str] = {}
    af_n = int(cfg.get("af_branch", {}).get("export_top_n", cfg.get("v5", {}).get("structure_top_n", 200)))
    pm_n = int(cfg.get("proteinmpnn", {}).get("export_top_n", cfg.get("v5", {}).get("structure_top_n", 200)))
    af = priority.head(af_n)
    pm = priority.head(pm_n)
    af_fasta = out_dir / "af2_af3_structure_tasks.fasta"
    af_jsonl = out_dir / "af2_af3_structure_tasks.jsonl"
    pm_fasta = out_dir / "proteinmpnn_design_tasks.fasta"
    pm_jsonl = out_dir / "proteinmpnn_design_tasks.jsonl"
    with af_fasta.open("w", encoding="utf-8") as ffa, af_jsonl.open("w", encoding="utf-8") as fjs:
        for _, row in af.iterrows():
            sid = row.get("candidate_id") or row.get("Seq_ID")
            ffa.write(f">{sid}|mutations={row.get('mutations','')}|v5_objective={row.get('v5_objective',np.nan):.6f}\n{row['sequence']}\n")
            fjs.write(json.dumps({"id": sid, "sequence": row["sequence"], "mutations": row.get("mutations", ""), "backend": cfg.get("af_branch", {}).get("backend", "af2_or_af3")}, ensure_ascii=False) + "\n")
    with pm_fasta.open("w", encoding="utf-8") as ffa, pm_jsonl.open("w", encoding="utf-8") as fjs:
        for _, row in pm.iterrows():
            sid = row.get("candidate_id") or row.get("Seq_ID")
            ffa.write(f">{sid}|mutations={row.get('mutations','')}\n{row['sequence']}\n")
            # ProteinMPNN receives candidate backbones after AF/PDB generation; this task file binds sequence IDs to downstream backbone names.
            fjs.write(json.dumps({"id": sid, "sequence": row["sequence"], "mutations": row.get("mutations", ""), "expected_backbone_pdb": f"{sid}.pdb"}, ensure_ascii=False) + "\n")
    outputs.update({"af_fasta": str(af_fasta), "af_jsonl": str(af_jsonl), "proteinmpnn_fasta": str(pm_fasta), "proteinmpnn_jsonl": str(pm_jsonl)})
    return outputs


def write_v5_reports(
    selected: pd.DataFrame,
    ranked: pd.DataFrame,
    priority: pd.DataFrame,
    gate_top: pd.DataFrame,
    diagnostics: pd.DataFrame,
    metrics: dict[str, Any],
    cfg: dict[str, Any],
    out_dir: Path,
) -> None:
    with (out_dir / "training_metrics_v5.json").open("w", encoding="utf-8") as f:
        json.dump(ensure_jsonable(metrics), f, indent=2, ensure_ascii=False)
    final_empty_reason = ""
    if selected.empty:
        if int(diagnostics.loc[0, "has_real_structure_pass"]) == 0 and bool(cfg.get("v5", {}).get("require_real_structure_for_final", True)):
            final_empty_reason = "Final Top-K is intentionally empty because no real AF/PDB structure metrics were imported."
        elif int(diagnostics.loc[0, "v5_final_eligible_pass"]) == 0:
            final_empty_reason = "Final Top-K is empty because no candidate passed both V5 basic and structure gates."
        else:
            final_empty_reason = "Final Top-K is empty after portfolio constraints."
    lines = [
        "# SynBio GFP V6.2 Adaptive Pipeline Report",
        "",
        "## Executive summary",
        "",
        f"- total_candidates: {len(ranked)}",
        f"- v5_basic_gate_pass: {int(diagnostics.loc[0, 'v5_basic_gate_pass'])}",
        f"- v5_structure_gate_pass: {int(diagnostics.loc[0, 'v5_structure_gate_pass'])}",
        f"- v5_final_eligible_pass: {int(diagnostics.loc[0, 'v5_final_eligible_pass'])}",
        f"- final_selected: {len(selected)}",
        f"- final_empty_reason: {final_empty_reason or 'NA'}",
        f"- esm_top6_overlap: {int(selected.get('is_in_esm_top6', pd.Series(dtype=bool)).sum()) if not selected.empty else 0}",
        f"- valid_structure_scores: {int(ranked.get('structure_score_valid', pd.Series(False, index=ranked.index)).sum())}",
        f"- valid_proteinmpnn_scores: {int(ranked.get('proteinmpnn_score_valid', pd.Series(False, index=ranked.index)).sum())}",
        "",
        "## Model metrics inherited from V3 training stage",
        "",
        f"- brightness_r2: {metrics.get('brightness_r2')}",
        f"- low_brightness_auc: {metrics.get('low_brightness_auc')}",
        f"- n_train: {metrics.get('n_train')}",
        f"- n_val: {metrics.get('n_val')}",
        "",
        "## V5 design logic",
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
    show_cols = ["candidate_id", "Seq_ID", "source_branch", "mutations", "v5_objective", "brightness_lcb", "low_brightness_risk", "structure_status", "v5_basic_gate", "v5_final_eligible"]
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
                f"- v6_bucket: {row.get('v6_bucket')}",
                f"- esm_baseline_rank: {row.get('esm_baseline_rank')}",
                f"- rank_delta_vs_esm: {row.get('rank_delta_vs_esm')}",
                f"- v5_objective: {row.get('v5_objective')}",
                f"- brightness_lcb: {row.get('brightness_lcb')}",
                f"- low_brightness_risk: {row.get('low_brightness_risk')}",
                f"- structure_status: {row.get('structure_status')}",
                f"- barrel_rmsd / pocket_rmsd / ion_network_rmsd: {row.get('barrel_rmsd')} / {row.get('pocket_rmsd')} / {row.get('ion_network_rmsd')}",
                "",
            ])
    (out_dir / "v6_pipeline_report.md").write_text("\n".join(lines), encoding="utf-8")


def write_submission(selected: pd.DataFrame, team_name: str, out_dir: Path) -> None:
    sub = pd.DataFrame(columns=["Team_Name", "Seq_ID", "Sequence"])
    if not selected.empty:
        sub = pd.DataFrame({"Team_Name": [team_name] * len(selected), "Seq_ID": selected["Seq_ID"].tolist(), "Sequence": selected["sequence"].tolist()})
    sub.to_csv(out_dir / "submission_v6.csv", index=False)
    sub.to_csv(out_dir / "submission_v62.csv", index=False)


def diagnostics_frame(df: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame([{
        "total_candidates": int(len(df)),
        "v5_basic_gate_pass": int(df.get("v5_basic_gate", pd.Series(False, index=df.index)).sum()),
        "v5_structure_gate_pass": int(df.get("v5_structure_gate", pd.Series(False, index=df.index)).sum()),
        "v5_final_eligible_pass": int(df.get("v5_final_eligible", pd.Series(False, index=df.index)).sum()),
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



def write_v62_survival_reports(ranked: pd.DataFrame, out_dir: Path) -> None:
    rows = []
    for branch, sub in ranked.groupby(ranked.get("source_branch", ranked.get("source", pd.Series("unknown", index=ranked.index))).astype(str)):
        rows.append({
            "source_branch": branch,
            "n": int(len(sub)),
            "basic_gate_survival_rate": float(sub.get("v5_basic_gate", pd.Series(False, index=sub.index)).mean()),
            "top200_candidate_rate_proxy": float((safe_num(sub, "structure_priority_rank_v6", 999999) <= 200).mean()) if "structure_priority_rank_v6" in sub.columns else 0.0,
            "structure_pass_rate": float(sub.get("structure_pass", pd.Series(False, index=sub.index)).mean()),
            "pmpnn_score_valid_rate": float(pd.to_numeric(sub.get("proteinmpnn_score_raw", pd.Series(np.nan, index=sub.index)), errors="coerce").notna().mean()),
            "final_eligible_rate": float(sub.get("v5_final_eligible", pd.Series(False, index=sub.index)).mean()),
            "low_brightness_safe_rate_010": float((safe_num(sub, "low_brightness_risk", 1.0) <= 0.10).mean()),
            "low_brightness_safe_rate_020": float((safe_num(sub, "low_brightness_risk", 1.0) <= 0.20).mean()),
            "mean_v62_objective": float(safe_num(sub, "v62_objective", 0.0).mean()),
            "mean_thermal_proxy_v62": float(safe_num(sub, "thermal_proxy_v62", 0.0).mean()),
        })
    survival = pd.DataFrame(rows).sort_values(["final_eligible_rate", "mean_v62_objective"], ascending=False) if rows else pd.DataFrame()
    survival.to_csv(out_dir / "v62_feedback_survival_by_branch.csv", index=False)
    lines = ["# V6.2 Feedback Survival Audit", ""]
    if survival.empty:
        lines.append("No source-branch survival records were available.")
    else:
        for _, row in survival.iterrows():
            lines.append(f"- {row['source_branch']}: n={row['n']}, basic={row['basic_gate_survival_rate']:.3f}, structure={row['structure_pass_rate']:.3f}, final={row['final_eligible_rate']:.3f}, low_risk<=0.10={row['low_brightness_safe_rate_010']:.3f}")
    (out_dir / "v62_feedback_survival_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_pipeline_v6(cfg: dict[str, Any]) -> dict[str, str]:
    out_dir = Path(cfg["paths"]["out_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)
    logger = setup_logger(out_dir, name="synbio_gfp_v6")
    set_seed(int(cfg.get("seed", 42)))
    start_time = time.time()
    (out_dir / "config_resolved_v5.json").write_text(json.dumps(ensure_jsonable(cfg), indent=2, ensure_ascii=False), encoding="utf-8")
    logger.info("V6.2 EVIDENCE-SAFE PIPELINE START")

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
    cand = add_v5_lcb_and_risk(cand, cfg)
    cand = add_v62_brightness_teacher(cand, cfg, logger)
    cand = add_epistasis_and_ddg_proxies(cand)
    cand = merge_proteinmpnn_scores(cand, cfg)
    cand = add_structure_branch(cand, cfg, out_dir, logger)
    cand = add_v5_gates(cand, cfg)
    cand = add_v5_objective(cand)
    feedback_policy = build_feedback_policy_from_sources(cand, cfg)
    cand = add_v6_feedback_scores(cand, feedback_policy)
    cand = add_v62_physical_thermal_proxy(cand)
    cand = add_v62_objective(cand)
    cand = add_v6_attribution(cand)

    cand["in_exclusion_list"] = cand["sequence"].isin(exclusion)
    cand.loc[cand["in_exclusion_list"], "v5_basic_gate"] = False
    cand.loc[cand["in_exclusion_list"], "v5_final_eligible"] = False
    cand = add_v5_objective(cand)
    cand = add_v6_feedback_scores(cand, feedback_policy)
    cand = add_v62_physical_thermal_proxy(cand)
    cand = add_v62_objective(cand)
    cand = add_v6_attribution(cand)
    write_feedback_reports(feedback_policy, out_dir)
    cand["sequence_sha1"] = cand["sequence"].astype(str).str.strip().str.upper().map(sequence_sha1)
    ranked = cand.sort_values("v6_objective", ascending=False).reset_index(drop=True)
    ranked.to_csv(out_dir / "all_ranked_candidates_v6.csv", index=False)
    ranked.to_csv(out_dir / "all_ranked_candidates_v62.csv", index=False)
    # Freeze the candidate universe so Stage 4 never silently regenerates a different pool.
    ranked.to_csv(out_dir / "candidate_universe_v62.csv", index=False)
    manifest_cols = [c for c in [
        "candidate_id", "Seq_ID", "sequence", "sequence_sha1", "sequence_length",
        "source_branch", "source", "generator_mode", "mutations", "mutation_count",
        "v62_objective", "v6_objective", "brightness_lcb", "brightness_lcb_norm",
        "low_brightness_risk", "thermal_proxy_v62", "v6_bucket"
    ] if c in ranked.columns]
    ranked[manifest_cols].to_csv(out_dir / "candidate_manifest_v62.csv", index=False)

    top_n = int(cfg.get("v5", {}).get("structure_top_n", cfg.get("structure", {}).get("top_n", 200)))
    priority = select_structure_priority_v6(ranked, cfg).head(top_n).copy()
    if "sequence_sha1" not in priority.columns:
        priority["sequence_sha1"] = priority["sequence"].astype(str).str.strip().str.upper().map(sequence_sha1)
    priority.to_csv(out_dir / "structure_priority_top200_v6.csv", index=False)
    priority.to_csv(out_dir / "structure_priority_top200_v62.csv", index=False)
    gate_top = ranked.sort_values("v5_objective", ascending=False).head(top_n).copy()
    gate_top.to_csv(out_dir / "ranked_top200_v6.csv", index=False)
    # Backward-compatible filename for earlier V5 analysis notebooks.
    gate_top.to_csv(out_dir / "structure_gate_top200_v6.csv", index=False)

    selected = select_portfolio(ranked, cfg)
    selected.to_csv(out_dir / "final_top6_v6.csv", index=False)
    selected.to_csv(out_dir / "final_top6_v62.csv", index=False)
    write_submission(selected, cfg.get("team_name", "YourTeamName"), out_dir)

    diagnostics = diagnostics_frame(ranked)
    diagnostics.to_csv(out_dir / "v6_gate_diagnostics.csv", index=False)
    mutation_frequency(ranked).to_csv(out_dir / "mutation_frequency_by_seed_v6.csv", index=False)
    write_v62_survival_reports(ranked, out_dir)
    ranked[~ranked["v5_basic_gate"].astype(bool)].head(2000).to_csv(out_dir / "rejected_by_basic_gate.csv", index=False)
    ranked[ranked["v5_basic_gate"].astype(bool) & ~ranked["v5_structure_gate"].astype(bool)].head(2000).to_csv(out_dir / "rejected_by_structure_gate.csv", index=False)
    # Backward-compatible explorer diagnostic file; it is now populated only when explorer-like sources exist.
    explorer_rej = ranked[ranked["source_branch"].astype(str).str.contains("explor", case=False, na=False) & ~ranked["v5_basic_gate"].astype(bool)].copy()
    explorer_rej.head(2000).to_csv(out_dir / "rejected_explorer_candidates.csv", index=False)
    export_paths = export_structure_tasks(priority, cfg, out_dir)
    write_v5_reports(selected, ranked, priority, gate_top, diagnostics, metrics, cfg, out_dir)
    logger.info("V6.2 EVIDENCE-SAFE PIPELINE DONE | seconds=%.2f | final=%d", time.time() - start_time, len(selected))
    return {
        "ranked_candidates": str(out_dir / "all_ranked_candidates_v6.csv"),
        "structure_priority": str(out_dir / "structure_priority_top200_v6.csv"),
        "final_top6": str(out_dir / "final_top6_v62.csv"),
        "submission": str(out_dir / "submission_v62.csv"),
        "diagnostics": str(out_dir / "v6_gate_diagnostics.csv"),
        "report": str(out_dir / "v6_pipeline_report.md"),
        **export_paths,
    }


def main() -> None:
    args = parse_args()
    cfg = resolve_v6_config(load_config(args.config), args)
    outputs = run_pipeline_v6(cfg)
    print("\nV6.2 evidence-safe outputs:")
    for k, v in outputs.items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
