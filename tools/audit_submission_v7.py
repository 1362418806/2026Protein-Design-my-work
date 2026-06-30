#!/usr/bin/env python3
"""Final submission legality and evidence audit for V7."""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

AA20_RE = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")


def load_table(path: str | Path) -> pd.DataFrame:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(p)
    return pd.read_excel(p) if p.suffix.lower() in {".xlsx", ".xls"} else pd.read_csv(p)


def load_exclusion(path: str | Path | None, data_dir: str | Path | None) -> set[str]:
    candidates: list[Path] = []
    if path:
        candidates.append(Path(path))
    if data_dir:
        candidates.append(Path(data_dir) / "Exclusion_List.csv")
    for p in candidates:
        if p.exists():
            df = load_table(p)
            col = "Sequence" if "Sequence" in df.columns else ("sequence" if "sequence" in df.columns else None)
            if col:
                return set(df[col].astype(str).str.strip().str.upper())
    return set()


def load_config(path: str | Path | None) -> dict[str, Any]:
    if not path:
        return {}
    p = Path(path)
    if not p.exists():
        return {}
    with p.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def mutation_count_from_row(row: pd.Series) -> int:
    if "mutation_count" in row and pd.notna(row["mutation_count"]):
        try:
            return int(float(row["mutation_count"]))
        except Exception:
            pass
    muts = str(row.get("mutations", ""))
    return len(re.findall(r"[ACDEFGHIKLMNPQRSTVWY]\d{1,3}[ACDEFGHIKLMNPQRSTVWY]", muts.upper()))


def main() -> None:
    p = argparse.ArgumentParser(description="Audit V7 final submission and selected-candidate evidence.")
    p.add_argument("--submission-csv", required=True)
    p.add_argument("--final-csv", required=True)
    p.add_argument("--data-dir", default="")
    p.add_argument("--exclusion-list-csv", default="")
    p.add_argument("--config", default="")
    p.add_argument("--out-json", required=True)
    p.add_argument("--out-md", default="")
    p.add_argument("--fail-on-error", action="store_true")
    args = p.parse_args()

    cfg = load_config(args.config)
    audit_cfg = cfg.get("v7", {}).get("final_audit", {})
    max_low = float(audit_cfg.get("max_final_low_brightness_risk", 0.20))
    preferred_count = int(audit_cfg.get("preferred_low_risk_count", 4))
    preferred_max = float(audit_cfg.get("preferred_low_risk_max", 0.10))
    max_normal_mut = int(audit_cfg.get("max_normal_mutation_count", 5))
    allow_one_six = bool(audit_cfg.get("allow_one_mutation_count_6", True))
    warn_max_agg = float(audit_cfg.get("warn_max_sequence_aggregation_risk", 0.55))
    warn_max_frag = float(audit_cfg.get("warn_max_sequence_fragility", 0.55))
    warn_min_bio = float(audit_cfg.get("warn_min_biology_validity_score", 0.45))
    hard_fail_bio = bool(audit_cfg.get("hard_fail_biology_soft_risk", False))

    sub = load_table(args.submission_csv)
    final = load_table(args.final_csv)
    exclusion = load_exclusion(args.exclusion_list_csv, args.data_dir)
    errors: list[str] = []
    warnings: list[str] = []

    required_sub_cols = {"Team_Name", "Seq_ID", "Sequence"}
    missing_cols = sorted(required_sub_cols - set(sub.columns))
    if missing_cols:
        errors.append(f"submission_missing_columns={missing_cols}")
    if len(sub) != 6:
        errors.append(f"submission_row_count={len(sub)} != 6")
    sequences = sub.get("Sequence", pd.Series([], dtype=str)).astype(str).str.strip().str.upper().tolist()
    if len(sequences) != len(set(sequences)):
        errors.append("submission_sequences_not_unique")
    for i, seq in enumerate(sequences, start=1):
        if not (220 <= len(seq) <= 250):
            errors.append(f"row_{i}_length={len(seq)} outside 220-250")
        if not seq.startswith("M"):
            errors.append(f"row_{i}_does_not_start_with_M")
        if not AA20_RE.fullmatch(seq):
            errors.append(f"row_{i}_has_invalid_amino_acids")
        if seq in exclusion:
            errors.append(f"row_{i}_sequence_in_exclusion_list")

    if len(final) != 6:
        errors.append(f"final_row_count={len(final)} != 6")
    bool_cols = ["has_real_structure", "structure_pass"]
    for col in bool_cols:
        if col not in final.columns:
            errors.append(f"final_missing_{col}")
        elif not final[col].astype(str).str.lower().eq("true").all():
            errors.append(f"final_{col}_not_all_true")
    if "structure_status" in final.columns and not final["structure_status"].astype(str).eq("rmsd_checked").all():
        errors.append("final_structure_status_not_all_rmsd_checked")
    if "in_exclusion_list" in final.columns and final["in_exclusion_list"].astype(str).str.lower().eq("true").any():
        errors.append("final_contains_exclusion_list_hit")
    for col in ["protected_mutation_count", "pocket_mutation_count", "ion_network_mutation_count"]:
        if col in final.columns and (pd.to_numeric(final[col], errors="coerce").fillna(999) != 0).any():
            errors.append(f"final_{col}_not_zero")
    if "low_brightness_risk" in final.columns:
        risk = pd.to_numeric(final["low_brightness_risk"], errors="coerce").fillna(1.0)
        if (risk > max_low).any():
            errors.append(f"final_low_brightness_risk_exceeds_{max_low}")
        if int((risk <= preferred_max).sum()) < preferred_count:
            warnings.append(f"only_{int((risk <= preferred_max).sum())}_final_candidates_have_low_brightness_risk<={preferred_max}")
    mut_counts = [mutation_count_from_row(row) for _, row in final.iterrows()]
    too_high = [m for m in mut_counts if m > max_normal_mut]
    if allow_one_six:
        if sum(m == max_normal_mut + 1 for m in mut_counts) > 1 or any(m > max_normal_mut + 1 for m in mut_counts):
            errors.append(f"final_mutation_count_policy_failed={mut_counts}")
    elif too_high:
        errors.append(f"final_mutation_count_policy_failed={mut_counts}")

    soft_available = int(final.get("softcheck_available", pd.Series(False, index=final.index)).astype(str).str.lower().eq("true").sum()) if not final.empty else 0
    soft_pass = int(final.get("oc2_pass", pd.Series(False, index=final.index)).astype(str).str.lower().eq("true").sum()) if not final.empty else 0
    if soft_available and soft_pass < soft_available:
        warnings.append(f"softcheck_available={soft_available}, softcheck_pass={soft_pass}")

    biology_warnings: list[str] = []
    if "v7_sequence_aggregation_risk" in final.columns:
        max_agg = float(pd.to_numeric(final["v7_sequence_aggregation_risk"], errors="coerce").fillna(1.0).max())
        if max_agg > warn_max_agg:
            biology_warnings.append(f"max_sequence_aggregation_risk={max_agg:.3f} > {warn_max_agg:.3f}")
    if "v7_sequence_fragility_score" in final.columns:
        max_frag = float(pd.to_numeric(final["v7_sequence_fragility_score"], errors="coerce").fillna(1.0).max())
        if max_frag > warn_max_frag:
            biology_warnings.append(f"max_sequence_fragility={max_frag:.3f} > {warn_max_frag:.3f}")
    if "v7_biology_validity_score" in final.columns:
        min_bio = float(pd.to_numeric(final["v7_biology_validity_score"], errors="coerce").fillna(0.0).min())
        if min_bio < warn_min_bio:
            biology_warnings.append(f"min_biology_validity_score={min_bio:.3f} < {warn_min_bio:.3f}")
    if biology_warnings:
        if hard_fail_bio:
            errors.extend(biology_warnings)
        else:
            warnings.extend(biology_warnings)

    summary = {
        "passed": not errors,
        "errors": errors,
        "warnings": warnings,
        "submission_rows": int(len(sub)),
        "final_rows": int(len(final)),
        "exclusion_loaded": bool(exclusion),
        "final_softcheck_available": soft_available,
        "final_softcheck_pass": soft_pass,
        "biology_soft_risk_hard_fail": hard_fail_bio,
        "mutation_counts": mut_counts,
    }
    out_json = Path(args.out_json)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    if args.out_md:
        lines = ["# V7 Final Audit", "", f"- passed: {summary['passed']}", f"- errors: {len(errors)}", f"- warnings: {len(warnings)}", f"- final_softcheck_available: {soft_available}", f"- final_softcheck_pass: {soft_pass}", ""]
        if errors:
            lines += ["## Errors", *[f"- {e}" for e in errors], ""]
        if warnings:
            lines += ["## Warnings", *[f"- {w}" for w in warnings], ""]
        Path(args.out_md).write_text("\n".join(lines), encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    if args.fail_on_error and errors:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
