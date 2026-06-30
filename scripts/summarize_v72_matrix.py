#!/usr/bin/env python3
"""Summarize the V7.2 matrix without modifying any experiment outputs.

All reported scores are proxies. In particular, thermal and brightness fields are not calibrated
wet-lab predictions and must not be interpreted as measured 72C retention or fluorescence.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def safe_num(frame: pd.DataFrame, column: str, default: float = math.nan) -> pd.Series:
    if column not in frame.columns:
        return pd.Series(default, index=frame.index, dtype=float)
    return pd.to_numeric(frame[column], errors="coerce").fillna(default)


def as_bool(values: pd.Series) -> pd.Series:
    return values.astype(str).str.strip().str.lower().isin({"1", "true", "yes"})


def mean_or_nan(series: pd.Series) -> float:
    value = pd.to_numeric(series, errors="coerce").mean()
    return float(value) if pd.notna(value) else math.nan


def min_or_nan(series: pd.Series) -> float:
    value = pd.to_numeric(series, errors="coerce").min()
    return float(value) if pd.notna(value) else math.nan


def max_or_nan(series: pd.Series) -> float:
    value = pd.to_numeric(series, errors="coerce").max()
    return float(value) if pd.notna(value) else math.nan


def clip01(value: float) -> float:
    return float(np.clip(value, 0.0, 1.0)) if math.isfinite(value) else 0.0


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def collect_run(run_dir: Path) -> tuple[dict[str, Any] | None, pd.DataFrame]:
    manifest_path = run_dir / "v72_experiment_manifest.json"
    final_dir = run_dir / "05_final"
    final_path = final_dir / "final_top6_v72.csv"
    audit_path = final_dir / "v72_submission_audit.json"
    metrics_path = run_dir / "02_structure" / "structure_metrics_v72.csv"
    if not manifest_path.is_file():
        return None, pd.DataFrame()
    manifest = read_json(manifest_path)
    if not final_path.is_file() or not audit_path.is_file() or not metrics_path.is_file():
        record = {
            "experiment_id": manifest.get("experiment_id", run_dir.name),
            "candidate_source": manifest.get("candidate_source", ""),
            "diagnostic_only": bool(manifest.get("diagnostic_only", False)),
            "completed": False,
            "audit_passed": False,
            "reason": "missing_final_or_audit_or_metrics",
        }
        return record, pd.DataFrame()

    final = pd.read_csv(final_path)
    metrics = pd.read_csv(metrics_path)
    audit = read_json(audit_path)
    if "sequence_sha1" not in final.columns and "sequence" in final.columns:
        final["sequence_sha1"] = final["sequence"].astype(str).str.strip().str.upper().map(
            lambda s: hashlib.sha1(s.encode("utf-8")).hexdigest()
        )
    metrics = metrics.drop_duplicates("sequence_sha1", keep="first")
    final = final.merge(
        metrics[
            [
                c
                for c in [
                    "sequence_sha1",
                    "pdb_exists_v72",
                    "pdb_sha256_v72",
                    "globalfit_status_v72",
                    "barrel_rmsd_globalfit_v72",
                    "pocket_rmsd_globalfit_v72",
                    "ion_network_rmsd_globalfit_v72",
                    "globalfit_geometry_score_v72",
                    "globalfit_local_gap_v72",
                ]
                if c in metrics.columns
            ]
        ],
        on="sequence_sha1",
        how="left",
    )

    final_rows = int(len(final))
    all_real = bool(as_bool(final.get("has_real_structure", pd.Series(False, index=final.index))).all()) if final_rows else False
    all_pass = bool(as_bool(final.get("structure_pass", pd.Series(False, index=final.index))).all()) if final_rows else False
    audit_passed = bool(audit.get("passed", False))
    hard_ok = bool(audit_passed and final_rows == 6 and all_real and all_pass)

    mean_lcb_norm = mean_or_nan(safe_num(final, "brightness_lcb_norm"))
    min_lcb = min_or_nan(safe_num(final, "brightness_lcb"))
    mean_lcb = mean_or_nan(safe_num(final, "brightness_lcb"))
    max_low_risk = max_or_nan(safe_num(final, "low_brightness_risk"))
    mean_low_risk = mean_or_nan(safe_num(final, "low_brightness_risk"))
    brightness_safety = (
        0.60 * clip01(mean_lcb_norm)
        + 0.40 * clip01(1.0 - max_low_risk if math.isfinite(max_low_risk) else 0.0)
    )

    mean_thermal_proxy = mean_or_nan(safe_num(final, "thermal_proxy_v62"))
    min_thermal_proxy = min_or_nan(safe_num(final, "thermal_proxy_v62"))
    mean_thermal_validity = mean_or_nan(safe_num(final, "v7_thermal_validity_score"))
    mean_solubility = mean_or_nan(safe_num(final, "v7_solubility_score"))
    max_aggregation = max_or_nan(safe_num(final, "v7_sequence_aggregation_risk"))
    max_fragility = max_or_nan(safe_num(final, "v7_sequence_fragility_score"))
    max_soft_risk = max_or_nan(safe_num(final, "v7_soft_risk_penalty"))
    thermal_proxy = (
        0.50 * clip01(mean_thermal_proxy)
        + 0.25 * clip01(mean_solubility)
        + 0.25 * clip01(1.0 - max_aggregation if math.isfinite(max_aggregation) else 0.0)
    )

    mean_plddt = mean_or_nan(safe_num(final, "mean_plddt"))
    max_local_barrel = max_or_nan(safe_num(final, "barrel_rmsd"))
    max_local_pocket = max_or_nan(safe_num(final, "pocket_rmsd"))
    max_local_ion = max_or_nan(safe_num(final, "ion_network_rmsd"))
    max_global_barrel = max_or_nan(safe_num(final, "barrel_rmsd_globalfit_v72"))
    max_global_pocket = max_or_nan(safe_num(final, "pocket_rmsd_globalfit_v72"))
    max_global_ion = max_or_nan(safe_num(final, "ion_network_rmsd_globalfit_v72"))
    mean_global_geometry = mean_or_nan(safe_num(final, "globalfit_geometry_score_v72"))
    mean_global_gap = mean_or_nan(safe_num(final, "globalfit_local_gap_v72"))
    mean_pmpnn = mean_or_nan(safe_num(final, "proteinmpnn_score_norm"))
    min_pmpnn = min_or_nan(safe_num(final, "proteinmpnn_score_norm"))
    mean_structure_geometry = mean_or_nan(safe_num(final, "structure_geometry_score"))
    plddt_score = clip01((mean_plddt - 70.0) / 25.0) if math.isfinite(mean_plddt) else 0.0
    structure_evidence = (
        0.35 * plddt_score
        + 0.30 * clip01(mean_global_geometry)
        + 0.20 * clip01(mean_pmpnn)
        + 0.15 * clip01(mean_structure_geometry)
    )

    biology_validity = mean_or_nan(safe_num(final, "v7_biology_validity_score"))
    pre_consensus_confidence = (
        0.30 * brightness_safety
        + 0.30 * thermal_proxy
        + 0.30 * structure_evidence
        + 0.10 * clip01(biology_validity)
    )
    pre_consensus_confidence *= 1.0 if hard_ok else 0.0

    record = {
        "experiment_id": manifest.get("experiment_id", run_dir.name),
        "candidate_source": manifest.get("candidate_source", ""),
        "diagnostic_only": bool(manifest.get("diagnostic_only", False)),
        "completed": True,
        "audit_passed": audit_passed,
        "audit_errors": len(audit.get("errors", [])),
        "audit_warnings": len(audit.get("warnings", [])),
        "hard_ok": hard_ok,
        "submission_rows": final_rows,
        "final_all_real_structure": all_real,
        "final_all_structure_pass": all_pass,
        "structure_msa_mode": manifest.get("parameters", {}).get("msa_mode", ""),
        "num_recycle": manifest.get("parameters", {}).get("num_recycle", ""),
        "num_models": manifest.get("parameters", {}).get("num_models", ""),
        "structure_top_n": manifest.get("parameters", {}).get("structure_top_n", ""),
        "mean_brightness_lcb": mean_lcb,
        "min_brightness_lcb": min_lcb,
        "mean_brightness_lcb_norm": mean_lcb_norm,
        "mean_low_brightness_risk": mean_low_risk,
        "max_low_brightness_risk": max_low_risk,
        "brightness_safety_proxy_v72": brightness_safety,
        "mean_thermal_proxy_v62": mean_thermal_proxy,
        "min_thermal_proxy_v62": min_thermal_proxy,
        "mean_thermal_validity_proxy_v7": mean_thermal_validity,
        "mean_solubility_proxy_v7": mean_solubility,
        "max_sequence_aggregation_risk": max_aggregation,
        "max_sequence_fragility": max_fragility,
        "max_soft_risk_penalty": max_soft_risk,
        "thermal_robustness_proxy_v72": thermal_proxy,
        "mean_plddt": mean_plddt,
        "max_barrel_rmsd_local": max_local_barrel,
        "max_pocket_rmsd_local": max_local_pocket,
        "max_ion_network_rmsd_local": max_local_ion,
        "max_barrel_rmsd_globalfit_v72": max_global_barrel,
        "max_pocket_rmsd_globalfit_v72": max_global_pocket,
        "max_ion_network_rmsd_globalfit_v72": max_global_ion,
        "mean_globalfit_geometry_score_v72": mean_global_geometry,
        "mean_globalfit_local_gap_v72": mean_global_gap,
        "mean_proteinmpnn_score_norm": mean_pmpnn,
        "min_proteinmpnn_score_norm": min_pmpnn,
        "mean_structure_geometry_score": mean_structure_geometry,
        "structure_evidence_proxy_v72": structure_evidence,
        "mean_biology_validity_proxy_v7": biology_validity,
        "pre_consensus_confidence_proxy_v72": pre_consensus_confidence,
        "submission_csv": str(final_dir / "submission_v72_for_upload.csv"),
        "submission_sha256": file_sha256(final_dir / "submission_v72_for_upload.csv"),
        "final_csv": str(final_path),
        "final_sha256": file_sha256(final_path),
        "metrics_csv": str(metrics_path),
        "metrics_sha256": file_sha256(metrics_path),
    }
    return record, final


def write_markdown(records: pd.DataFrame, path: Path, recommendation: str) -> None:
    lines = [
        "# V7.2 Matrix Summary",
        "",
        "All confidence, brightness, and thermal fields below are model proxies, not calibrated wet-lab probabilities.",
        "Thermal values do not claim to predict 72C/10min fluorescence retention directly.",
        "",
        f"Recommended matrix proxy winner: `{recommendation or 'none'}`",
        "",
        "| Experiment | Hard gate | Confidence proxy | Brightness safety | Thermal proxy | Structure evidence | Consensus |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    if records.empty:
        lines.append("| None | - | - | - | - | - | - |")
    else:
        order = records.sort_values(["hard_ok", "matrix_confidence_proxy_v72"], ascending=[False, False])
        for _, row in order.iterrows():
            lines.append(
                "| {experiment_id} | {hard_ok} | {matrix_confidence_proxy_v72:.4f} | "
                "{brightness_safety_proxy_v72:.4f} | {thermal_robustness_proxy_v72:.4f} | "
                "{structure_evidence_proxy_v72:.4f} | {submission_consensus_proxy_v72:.4f} |".format(
                    experiment_id=row.get("experiment_id", ""),
                    hard_ok=bool(row.get("hard_ok", False)),
                    matrix_confidence_proxy_v72=float(row.get("matrix_confidence_proxy_v72", 0.0) or 0.0),
                    brightness_safety_proxy_v72=float(row.get("brightness_safety_proxy_v72", 0.0) or 0.0),
                    thermal_robustness_proxy_v72=float(row.get("thermal_robustness_proxy_v72", 0.0) or 0.0),
                    structure_evidence_proxy_v72=float(row.get("structure_evidence_proxy_v72", 0.0) or 0.0),
                    submission_consensus_proxy_v72=float(row.get("submission_consensus_proxy_v72", 0.0) or 0.0),
                )
            )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Summarize an isolated V7.2 matrix.")
    parser.add_argument("--matrix-root", required=True)
    args = parser.parse_args()
    matrix_root = Path(args.matrix_root)
    run_root = matrix_root / "runs"
    if not run_root.is_dir():
        raise FileNotFoundError(run_root)

    records: list[dict[str, Any]] = []
    final_tables: dict[str, pd.DataFrame] = {}
    for run_dir in sorted(path for path in run_root.iterdir() if path.is_dir()):
        record, final = collect_run(run_dir)
        if record is not None:
            records.append(record)
            if not final.empty:
                final_tables[str(record["experiment_id"])] = final

    summary = pd.DataFrame(records)
    if summary.empty:
        summary.to_csv(matrix_root / "matrix_submission_metrics_v72.csv", index=False)
        write_markdown(summary, matrix_root / "matrix_summary_v72.md", "")
        return

    eligible_ids = summary.loc[
        summary["completed"] & summary["hard_ok"] & ~summary["diagnostic_only"], "experiment_id"
    ].astype(str).tolist()
    source_count = max(1, summary.loc[summary["experiment_id"].astype(str).isin(eligible_ids), "candidate_source"].nunique())

    occurrence_rows: list[dict[str, Any]] = []
    for experiment_id in eligible_ids:
        frame = final_tables.get(experiment_id, pd.DataFrame())
        for _, row in frame.iterrows():
            occurrence_rows.append(
                {
                    "experiment_id": experiment_id,
                    "candidate_source": summary.loc[
                        summary["experiment_id"].astype(str).eq(experiment_id), "candidate_source"
                    ].iloc[0],
                    "sequence_sha1": str(row.get("sequence_sha1", "")),
                    "candidate_id": str(row.get("candidate_id", "")),
                    "Seq_ID": str(row.get("Seq_ID", "")),
                    "sequence": str(row.get("sequence", "")),
                }
            )
    occurrence = pd.DataFrame(occurrence_rows)
    if not occurrence.empty:
        grouped = (
            occurrence.groupby("sequence_sha1", as_index=False)
            .agg(
                occurrences=("experiment_id", "nunique"),
                source_diversity=("candidate_source", "nunique"),
                candidate_id=("candidate_id", "first"),
                Seq_ID=("Seq_ID", "first"),
                sequence=("sequence", "first"),
            )
        )
        denominator = max(1, len(eligible_ids))
        grouped["exact_recurrence_fraction_v72"] = grouped["occurrences"] / denominator
        grouped["source_diversity_fraction_v72"] = grouped["source_diversity"] / source_count
        grouped["candidate_consensus_proxy_v72"] = (
            0.75 * grouped["exact_recurrence_fraction_v72"]
            + 0.25 * grouped["source_diversity_fraction_v72"]
        )
    else:
        grouped = pd.DataFrame(
            columns=[
                "sequence_sha1",
                "occurrences",
                "source_diversity",
                "candidate_id",
                "Seq_ID",
                "sequence",
                "exact_recurrence_fraction_v72",
                "source_diversity_fraction_v72",
                "candidate_consensus_proxy_v72",
            ]
        )
    grouped.to_csv(matrix_root / "matrix_candidate_consensus_v72.csv", index=False)

    consensus_map = grouped.set_index("sequence_sha1")["candidate_consensus_proxy_v72"].to_dict() if not grouped.empty else {}
    submission_consensus: list[float] = []
    for _, record in summary.iterrows():
        experiment_id = str(record["experiment_id"])
        frame = final_tables.get(experiment_id, pd.DataFrame())
        values = [float(consensus_map.get(str(sha), 0.0)) for sha in frame.get("sequence_sha1", pd.Series(dtype=str))]
        submission_consensus.append(float(np.mean(values)) if values else 0.0)
    summary["submission_consensus_proxy_v72"] = submission_consensus
    summary["matrix_confidence_proxy_v72"] = (
        0.85 * pd.to_numeric(summary["pre_consensus_confidence_proxy_v72"], errors="coerce").fillna(0.0)
        + 0.15 * pd.to_numeric(summary["submission_consensus_proxy_v72"], errors="coerce").fillna(0.0)
    ) * summary["hard_ok"].astype(float)

    eligible = summary[summary["completed"] & summary["hard_ok"] & ~summary["diagnostic_only"]].copy()
    recommended = ""
    if not eligible.empty:
        recommended = str(
            eligible.sort_values(
                ["matrix_confidence_proxy_v72", "brightness_safety_proxy_v72", "thermal_robustness_proxy_v72"],
                ascending=False,
            ).iloc[0]["experiment_id"]
        )

    summary = summary.sort_values(["hard_ok", "matrix_confidence_proxy_v72"], ascending=[False, False]).reset_index(drop=True)
    summary.to_csv(matrix_root / "matrix_submission_metrics_v72.csv", index=False)
    write_markdown(summary, matrix_root / "matrix_summary_v72.md", recommended)
    recommendation_payload = {
        "recommended_experiment_id": recommended,
        "selection_rule": (
            "Highest V7.2 matrix confidence proxy among non-diagnostic experiments that passed the "
            "audit and primary structure gate. This is a reporting recommendation, not a calibrated phenotype probability."
        ),
        "eligible_experiments": eligible["experiment_id"].astype(str).tolist(),
    }
    (matrix_root / "matrix_recommendation_v72.json").write_text(
        json.dumps(recommendation_payload, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    (matrix_root / "MATRIX_RECOMMENDATION_V72.txt").write_text(
        (recommended or "NO_ELIGIBLE_EXPERIMENT") + "\n", encoding="utf-8"
    )

    for _, record in summary.iterrows():
        experiment_id = str(record["experiment_id"])
        final_dir = matrix_root / "runs" / experiment_id / "05_final"
        per_submission = {
            key: (None if pd.isna(value) else value)
            for key, value in record.to_dict().items()
        }
        (final_dir / "submission_metrics_v72.json").write_text(
            json.dumps(per_submission, indent=2, ensure_ascii=False), encoding="utf-8"
        )

    print(f"Wrote V7.2 matrix summary: {matrix_root / 'matrix_submission_metrics_v72.csv'}")
    print(f"Recommended matrix proxy experiment: {recommended or 'none'}")


if __name__ == "__main__":
    main()
