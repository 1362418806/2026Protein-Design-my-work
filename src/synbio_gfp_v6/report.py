from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from .utils import ensure_jsonable


def build_design_reason(row: pd.Series) -> str:
    parts = [f"Layer={row.get('selection_layer', 'NA')}"]
    parts.append(f"mutations={row.get('mutation_count', 'NA')}")
    if row.get("literature_events"):
        parts.append(f"literature_events={row.get('literature_events')}")
    if row.get("source"):
        parts.append(f"source={row.get('source')}")
    if row.get("structure_status"):
        parts.append(f"structure={row.get('structure_status')}")
    return " | ".join(parts)


def write_reports(selected: pd.DataFrame, ranked: pd.DataFrame, metrics: dict[str, Any], config: dict[str, Any], out_dir: str | Path) -> None:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    selected_report = selected.copy()
    selected_report["design_reason"] = selected_report.apply(build_design_reason, axis=1)
    cols = [c for c in [
        "Seq_ID", "selection_layer", "source", "mutation_count", "mutation_list",
        "brightness_pred", "brightness_uncertainty", "brightness_lcb", "low_brightness_risk",
        "stability_proxy", "foldability_proxy", "literature_prior", "fitness_landscape_prior",
        "tgp_surface_prior", "pocket_prior", "superfolder_protection_score", "structure_score",
        "protected_site_mutations", "surface_mutations", "loop_mutations", "risk_explanation",
        "design_reason", "sequence",
    ] if c in selected_report.columns]
    selected_report[cols].to_csv(out / "selected_candidate_report.csv", index=False)
    with open(out / "training_metrics.json", "w", encoding="utf-8") as f:
        json.dump(ensure_jsonable(metrics), f, indent=2, ensure_ascii=False)
    report_lines = [
        "# SynBio GFP V3 Pipeline Test Report",
        "",
        "## Model metrics",
        "",
        f"- brightness_r2: {metrics.get('brightness_r2')}",
        f"- low_brightness_auc: {metrics.get('low_brightness_auc')}",
        f"- n_train: {metrics.get('n_train')}",
        f"- n_val: {metrics.get('n_val')}",
        "",
        "## Selected candidates",
        "",
    ]
    for _, row in selected_report.iterrows():
        report_lines.extend([
            f"### {row.get('Seq_ID')}",
            f"- selection_layer: {row.get('selection_layer')}",
            f"- mutation_count: {row.get('mutation_count')}",
            f"- mutation_list: {row.get('mutation_list')}",
            f"- brightness_pred / LCB: {row.get('brightness_pred')} / {row.get('brightness_lcb')}",
            f"- low_brightness_risk: {row.get('low_brightness_risk')}",
            f"- stability_proxy: {row.get('stability_proxy')}",
            f"- literature_prior: {row.get('literature_prior')}",
            f"- structure_status: {row.get('structure_status')}",
            f"- design_reason: {row.get('design_reason')}",
            f"- risk_explanation: {row.get('risk_explanation')}",
            "",
        ])
    (out / "v3_pipeline_report.md").write_text("\n".join(report_lines), encoding="utf-8")


def write_submission(selected: pd.DataFrame, team_name: str, out_dir: str | Path) -> None:
    sub = pd.DataFrame({
        "Team_Name": [team_name] * len(selected),
        "Seq_ID": selected["Seq_ID"].tolist(),
        "Sequence": selected["sequence"].tolist(),
    })
    sub.to_csv(Path(out_dir) / "submission.csv", index=False)
