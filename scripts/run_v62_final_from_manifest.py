#!/usr/bin/env python3
"""V6.2 final selector from a frozen Stage-1 candidate universe.

Stage 4 must not regenerate candidates. This script loads the exact Stage-1
candidate table, merges structure and ProteinMPNN evidence by sequence_sha1 first,
then emits final Top-6 only from evidence-backed candidates.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from run_v62_complete import (  # noqa: E402
    add_epistasis_and_ddg_proxies,
    add_v5_gates,
    add_v5_objective,
    add_v6_attribution,
    add_v62_objective,
    add_v62_physical_thermal_proxy,
    diagnostics_frame,
    ensure_jsonable,
    load_config,
    minmax01,
    mutation_frequency,
    safe_num,
    select_portfolio,
    setup_logger,
    structure_geometry_from_metrics,
    write_submission,
    write_v5_reports,
    write_v62_survival_reports,
)


def sequence_sha1(seq: str) -> str:
    return hashlib.sha1(str(seq).strip().upper().encode("utf-8")).hexdigest()


def load_table(path: str | Path) -> pd.DataFrame:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(p)
    return pd.read_excel(p) if p.suffix.lower() in {".xlsx", ".xls"} else pd.read_csv(p)


def ensure_keys(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if "sequence" not in out.columns:
        raise ValueError("Frozen candidate table must contain sequence.")
    out["sequence"] = out["sequence"].astype(str).str.strip().str.upper()
    if "sequence_sha1" not in out.columns:
        out["sequence_sha1"] = out["sequence"].map(sequence_sha1)
    if "candidate_id" not in out.columns:
        out["candidate_id"] = [f"cand_{i:06d}" for i in range(len(out))]
    if "Seq_ID" not in out.columns:
        out["Seq_ID"] = out["candidate_id"]
    return out


def merge_structure_metrics(cand: pd.DataFrame, metrics_csv: str | Path) -> pd.DataFrame:
    out = cand.copy()
    metrics = ensure_keys(load_table(metrics_csv))
    keep = [c for c in [
        "sequence_sha1", "candidate_id", "Seq_ID", "sequence", "barrel_rmsd", "pocket_rmsd",
        "ion_network_rmsd", "mean_plddt", "pdb_path", "has_real_structure", "structure_status", "structure_pass"
    ] if c in metrics.columns]
    metrics = metrics[keep].copy()
    # Prefer exact-sequence structure rows, then rows with valid RMSD/pass evidence.
    metrics["_pass_rank"] = metrics.get("structure_pass", pd.Series(False, index=metrics.index)).astype(str).str.lower().eq("true").astype(int)
    metrics["_rmsd_rank"] = metrics[[c for c in ["barrel_rmsd", "pocket_rmsd", "ion_network_rmsd"] if c in metrics.columns]].notna().sum(axis=1)
    metrics = metrics.sort_values(["sequence_sha1", "_pass_rank", "_rmsd_rank"], ascending=[True, False, False]).drop_duplicates("sequence_sha1")
    drop_cols = [c for c in ["barrel_rmsd", "pocket_rmsd", "ion_network_rmsd", "mean_plddt", "pdb_path", "has_real_structure", "structure_status", "structure_pass"] if c in out.columns]
    out = out.drop(columns=drop_cols, errors="ignore").merge(metrics.drop(columns=["candidate_id", "Seq_ID", "sequence", "_pass_rank", "_rmsd_rank"], errors="ignore"), on="sequence_sha1", how="left")
    out["structure_status"] = out.get("structure_status", pd.Series("not_checked", index=out.index)).fillna("not_checked")

    rmsd_cols = ["barrel_rmsd", "pocket_rmsd", "ion_network_rmsd"]
    for col in rmsd_cols:
        out[col] = pd.to_numeric(out.get(col, pd.Series(np.nan, index=out.index)), errors="coerce")

    status_real = out["structure_status"].astype(str).eq("rmsd_checked")
    rmsd_complete = out[rmsd_cols].notna().all(axis=1)
    imported_real = out.get("has_real_structure", pd.Series(False, index=out.index)).astype(str).str.lower().eq("true")
    # Final evidence requires reference-based RMSD rows; confidence-proxy metrics must not be upgraded to real structure evidence.
    out["has_real_structure"] = status_real & rmsd_complete & imported_real
    out["structure_pass"] = out.get("structure_pass", pd.Series(False, index=out.index)).astype(str).str.lower().eq("true") & out["has_real_structure"].astype(bool)
    out["structure_geometry_score"] = structure_geometry_from_metrics(out)
    out.loc[~out["has_real_structure"].astype(bool), "structure_geometry_score"] = 0.05
    return out


def merge_proteinmpnn(cand: pd.DataFrame, proteinmpnn_csv: str | Path | None) -> pd.DataFrame:
    out = cand.copy()
    out["proteinmpnn_status"] = out.get("proteinmpnn_status", "not_provided")
    out["proteinmpnn_score_raw"] = pd.to_numeric(out.get("proteinmpnn_score_raw", pd.Series(np.nan, index=out.index)), errors="coerce")
    out["proteinmpnn_score_norm"] = pd.to_numeric(out.get("proteinmpnn_score_norm", pd.Series(0.5, index=out.index)), errors="coerce").fillna(0.5)
    if not proteinmpnn_csv or not Path(proteinmpnn_csv).exists():
        return out
    pm = load_table(proteinmpnn_csv)
    if pm.empty:
        return out
    if "sequence" in pm.columns and "sequence_sha1" not in pm.columns:
        pm["sequence_sha1"] = pm["sequence"].astype(str).str.strip().str.upper().map(sequence_sha1)
    score_col = next((c for c in ["proteinmpnn_score", "proteinmpnn_score_raw", "score"] if c in pm.columns), None)
    if score_col is None:
        return out
    if "sequence_sha1" in pm.columns:
        key = "sequence_sha1"
    else:
        key = next((c for c in ["candidate_id", "Seq_ID", "seq_id"] if c in pm.columns and c in out.columns), None)
    if key is None:
        return out
    small = pm[[key, score_col]].drop_duplicates(key).rename(columns={score_col: "proteinmpnn_score_raw_imported"})
    out = out.merge(small, on=key, how="left")
    imported = pd.to_numeric(out["proteinmpnn_score_raw_imported"], errors="coerce")
    out["proteinmpnn_score_raw"] = imported.combine_first(out["proteinmpnn_score_raw"])
    out["proteinmpnn_score_norm"] = minmax01(out["proteinmpnn_score_raw"], higher_is_better=False, default=0.5)
    out["proteinmpnn_status"] = np.where(out["proteinmpnn_score_raw"].notna(), "score_imported", out["proteinmpnn_status"])
    out = out.drop(columns=["proteinmpnn_score_raw_imported"], errors="ignore")
    return out


def write_underfilled_recovery(ranked: pd.DataFrame, out_dir: Path, k: int) -> None:
    need = ranked.copy()
    need["needs_structure_recovery"] = ~need.get("structure_pass", pd.Series(False, index=need.index)).astype(bool)
    need = need[need["needs_structure_recovery"] & need.get("v5_basic_gate", pd.Series(False, index=need.index)).astype(bool)]
    need = need.sort_values("v6_objective", ascending=False).head(max(40, k * 20))
    cols = [c for c in ["candidate_id", "Seq_ID", "sequence", "sequence_sha1", "source_branch", "mutations", "v6_objective", "low_brightness_risk", "mutation_count", "structure_status", "structure_pass"] if c in need.columns]
    need[cols].to_csv(out_dir / "final_underfilled_needs_recovery.csv", index=False)


def main() -> None:
    p = argparse.ArgumentParser(description="Select V6.2 final Top-6 from frozen Stage-1 candidates.")
    p.add_argument("--config", default=str(ROOT / "configs" / "v62_adaptive.yaml"))
    p.add_argument("--frozen-candidates-csv", required=True)
    p.add_argument("--structure-metrics-csv", required=True)
    p.add_argument("--proteinmpnn-score-csv", default="")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--team-name", default="YourTeamName")
    args = p.parse_args()

    cfg = load_config(args.config)
    cfg.setdefault("v5", {})["require_real_structure_for_final"] = True
    cfg.setdefault("paths", {})["out_dir"] = args.out_dir
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    logger = setup_logger(out_dir, name="synbio_gfp_v62_final")

    cand = ensure_keys(load_table(args.frozen_candidates_csv))
    cand = merge_structure_metrics(cand, args.structure_metrics_csv)
    cand = merge_proteinmpnn(cand, args.proteinmpnn_score_csv)
    cand = add_epistasis_and_ddg_proxies(cand) if "epistasis_penalty" not in cand.columns else cand
    cand = add_v5_gates(cand, cfg)
    cand = add_v5_objective(cand)
    cand = add_v62_physical_thermal_proxy(cand)
    cand = add_v62_objective(cand)
    cand = add_v6_attribution(cand)
    ranked = cand.sort_values("v6_objective", ascending=False).reset_index(drop=True)
    ranked.to_csv(out_dir / "all_ranked_candidates_v62.csv", index=False)
    ranked.to_csv(out_dir / "all_ranked_candidates_v6.csv", index=False)

    selected = select_portfolio(ranked, cfg)
    selected.to_csv(out_dir / "final_top6_v62.csv", index=False)
    selected.to_csv(out_dir / "final_top6_v61.csv", index=False)
    selected.to_csv(out_dir / "final_top6_v6.csv", index=False)
    write_submission(selected, args.team_name, out_dir)

    diagnostics = diagnostics_frame(ranked)
    diagnostics.to_csv(out_dir / "v6_gate_diagnostics.csv", index=False)
    mutation_frequency(ranked).to_csv(out_dir / "mutation_frequency_by_seed_v6.csv", index=False)
    write_v62_survival_reports(ranked, out_dir)
    write_underfilled_recovery(ranked, out_dir, int(cfg.get("v5", {}).get("final_k", 6)))
    priority = ranked.head(int(cfg.get("v5", {}).get("structure_top_n", 200))).copy()
    write_v5_reports(selected, ranked, priority, priority, diagnostics, {"stage4_from_frozen": True}, cfg, out_dir)

    summary = {
        "frozen_candidates": args.frozen_candidates_csv,
        "structure_metrics": args.structure_metrics_csv,
        "proteinmpnn_scores": args.proteinmpnn_score_csv,
        "ranked_rows": int(len(ranked)),
        "final_rows": int(len(selected)),
        "final_all_real_structure": bool(len(selected) > 0 and selected.get("has_real_structure", pd.Series(False)).astype(bool).all()),
        "final_all_structure_pass": bool(len(selected) > 0 and selected.get("structure_pass", pd.Series(False)).astype(bool).all()),
    }
    (out_dir / "v62_final_summary.json").write_text(json.dumps(ensure_jsonable(summary), indent=2, ensure_ascii=False), encoding="utf-8")
    logger.info("V6.2 frozen final done | final=%d", len(selected))
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
