#!/usr/bin/env python3
"""Select a V7.2 final Top-6 from a frozen candidate universe.

The selection mathematics intentionally reuses the V7.1 implementation. This wrapper only
changes artifact names and removes legacy V6/V6.2 compatibility outputs so every matrix
experiment is self-contained and traceable.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

import run_v7_final_from_manifest as v7


def write_submission_v72(selected: pd.DataFrame, team_name: str, out_dir: Path) -> Path:
    """Write the sole canonical V7.2 submission artifact for one matrix experiment."""
    submission = pd.DataFrame(columns=["Team_Name", "Seq_ID", "Sequence"])
    if not selected.empty:
        submission = pd.DataFrame(
            {
                "Team_Name": [team_name] * len(selected),
                "Seq_ID": selected["Seq_ID"].tolist(),
                "Sequence": selected["sequence"].tolist(),
            }
        )
    output_path = out_dir / "submission_v72_for_upload.csv"
    submission.to_csv(output_path, index=False)
    return output_path


def write_v72_summary(selected: pd.DataFrame, ranked: pd.DataFrame, args: argparse.Namespace, out_dir: Path) -> None:
    """Persist a compact, experiment-local summary without legacy version aliases."""
    summary = {
        "frozen_candidates": args.frozen_candidates_csv,
        "structure_metrics": args.structure_metrics_csv,
        "proteinmpnn_scores": args.proteinmpnn_score_csv,
        "softcheck_scores": args.softcheck_csv,
        "ranked_rows": int(len(ranked)),
        "final_rows": int(len(selected)),
        "final_all_real_structure": bool(
            len(selected) > 0
            and selected.get("has_real_structure", pd.Series(False, index=selected.index)).astype(bool).all()
        ),
        "final_all_structure_pass": bool(
            len(selected) > 0
            and selected.get("structure_pass", pd.Series(False, index=selected.index)).astype(bool).all()
        ),
        "final_softcheck_available": int(
            selected.get("softcheck_available", pd.Series(False, index=selected.index)).astype(bool).sum()
        ) if not selected.empty else 0,
        "final_softcheck_pass": int(
            selected.get("oc2_pass", pd.Series(False, index=selected.index)).astype(bool).sum()
        ) if not selected.empty else 0,
        "ranked_softcheck_available": int(
            ranked.get("softcheck_available", pd.Series(False, index=ranked.index)).astype(bool).sum()
        ),
        "final_mean_biology_validity": float(v7.safe_num(selected, "v7_biology_validity_score", 0.0).mean())
        if not selected.empty else 0.0,
        "final_min_thermal_validity": float(v7.safe_num(selected, "v7_thermal_validity_score", 0.0).min())
        if not selected.empty else 0.0,
        "final_max_sequence_aggregation_risk": float(
            v7.safe_num(selected, "v7_sequence_aggregation_risk", 0.0).max()
        ) if not selected.empty else 0.0,
        "final_max_sequence_fragility": float(
            v7.safe_num(selected, "v7_sequence_fragility_score", 0.0).max()
        ) if not selected.empty else 0.0,
        "final_min_neighborhood_support": float(
            v7.safe_num(selected, "v7_neighborhood_support", 0.0).min()
        ) if not selected.empty else 0.0,
    }
    (out_dir / "v72_final_summary.json").write_text(
        json.dumps(v7.ensure_jsonable(summary), indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Select V7.2 final Top-6 from frozen Stage-1 candidates with optional soft validation."
    )
    parser.add_argument("--config", required=True)
    parser.add_argument("--frozen-candidates-csv", required=True)
    parser.add_argument("--structure-metrics-csv", required=True)
    parser.add_argument("--proteinmpnn-score-csv", default="")
    parser.add_argument("--softcheck-csv", default="")
    parser.add_argument("--data-dir", required=True)
    parser.add_argument("--exclusion-list-csv", default="")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--team-name", default="YourTeamName")
    args = parser.parse_args()

    cfg = v7.load_config(args.config)
    cfg.setdefault("v5", {})["require_real_structure_for_final"] = True
    cfg.setdefault("paths", {})["out_dir"] = args.out_dir
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    cand = v7.ensure_keys(v7.load_table(args.frozen_candidates_csv))
    exclusion = v7.load_exclusion_sequences(args.data_dir, args.exclusion_list_csv)
    cand["in_exclusion_list"] = cand["sequence"].isin(exclusion) if exclusion else False
    cand = v7.merge_structure_metrics(cand, args.structure_metrics_csv)
    cand = v7.merge_proteinmpnn(cand, args.proteinmpnn_score_csv)
    cand = v7.merge_softcheck(cand, args.softcheck_csv)
    # Recompute the complete epistasis proxy bundle unless all downstream objective fields are present.
    required_epistasis_fields = {"epistasis_penalty", "epistasis_score"}
    cand = v7.add_epistasis_and_ddg_proxies(cand) if not required_epistasis_fields.issubset(cand.columns) else cand
    cand = v7.add_v5_gates(cand, cfg)
    cand.loc[cand["in_exclusion_list"].astype(bool), ["v5_basic_gate", "v5_final_eligible"]] = False
    cand = v7.add_v5_objective(cand)
    cand = v7.add_v62_physical_thermal_proxy(cand)
    cand = v7.add_v62_objective(cand)
    cand = v7.add_v7_biology_soft_risk(cand, cfg)
    cand = v7.add_v7_confidence_and_objective(cand, cfg)
    cand = v7.add_v6_attribution(cand)

    ranked = cand.sort_values("v7_objective", ascending=False).reset_index(drop=True)
    ranked.to_csv(out_dir / "all_ranked_candidates_v72.csv", index=False)

    risk_cols = [
        c
        for c in [
            "candidate_id",
            "Seq_ID",
            "sequence_sha1",
            "source_branch",
            "mutations",
            "v7_objective",
            "v7_confidence_score",
            "v7_biology_validity_score",
            "v7_thermal_validity_score",
            "v7_solubility_score",
            "v7_sequence_aggregation_risk",
            "v7_sequence_fragility_score",
            "v7_neighborhood_support",
            "v7_soft_risk_penalty",
            "low_brightness_risk",
            "brightness_lcb",
            "brightness_lcb_norm",
            "thermal_proxy_v62",
            "mutation_count",
            "structure_pass",
            "softcheck_available",
            "oc2_pass",
        ]
        if c in ranked.columns
    ]
    ranked[risk_cols].head(500).to_csv(out_dir / "v72_soft_risk_top500.csv", index=False)

    selected = v7.select_portfolio(ranked, cfg)
    selected.to_csv(out_dir / "final_top6_v72.csv", index=False)
    selected[[c for c in risk_cols if c in selected.columns]].to_csv(
        out_dir / "v72_final_soft_risk_report.csv", index=False
    )
    write_submission_v72(selected, args.team_name, out_dir)

    diagnostics = v7.diagnostics_frame(ranked)
    diagnostics["softcheck_available_pass"] = int(
        ranked.get("softcheck_available", pd.Series(False, index=ranked.index)).astype(bool).sum()
    )
    diagnostics["softcheck_pass_true"] = int(
        ranked.get("oc2_pass", pd.Series(False, index=ranked.index)).astype(bool).sum()
    )
    diagnostics["in_exclusion_list"] = int(
        ranked.get("in_exclusion_list", pd.Series(False, index=ranked.index)).astype(bool).sum()
    )
    diagnostics.to_csv(out_dir / "v72_gate_diagnostics.csv", index=False)

    v7.write_underfilled_recovery(ranked, out_dir, int(cfg.get("v5", {}).get("final_k", 6)))
    underfilled = out_dir / "final_underfilled_needs_recovery.csv"
    if underfilled.exists():
        underfilled.rename(out_dir / "final_underfilled_needs_recovery_v72.csv")
    write_v72_summary(selected, ranked, args, out_dir)
    print(f"V7.2 frozen final done | final={len(selected)} | out_dir={out_dir}")


if __name__ == "__main__":
    main()
