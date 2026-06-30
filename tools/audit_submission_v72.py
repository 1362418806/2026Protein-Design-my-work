#!/usr/bin/env python3
"""Fail-closed V7.2 submission audit for one independent matrix experiment."""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path
from typing import Any

import pandas as pd

AA_RE = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")


def sequence_sha1(sequence: str) -> str:
    return hashlib.sha1(str(sequence).strip().upper().encode("utf-8")).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def as_bool(value: Any) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes"}


def read_csv(path: Path, label: str) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"{label} does not exist: {path}")
    return pd.read_csv(path)


def load_exclusion(path: Path) -> set[str]:
    frame = read_csv(path, "Exclusion list")
    col = "Sequence" if "Sequence" in frame.columns else "sequence" if "sequence" in frame.columns else None
    if col is None:
        raise ValueError("Exclusion list has no Sequence/sequence column.")
    values = set(frame[col].astype(str).str.strip().str.upper())
    if not values:
        raise ValueError("Exclusion list is empty; refusing fail-open audit.")
    return values


def write_markdown(report: dict[str, Any], path: Path) -> None:
    lines = ["# V7.2 Submission Audit", "", f"- passed: `{report['passed']}`", ""]
    for section in ("errors", "warnings", "checks"):
        lines.append(f"## {section.title()}")
        entries = report.get(section, [])
        if not entries:
            lines.append("- None")
        else:
            for entry in entries:
                lines.append(f"- {entry}")
        lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Audit a V7.2 final submission and its structure evidence.")
    parser.add_argument("--submission-csv", required=True)
    parser.add_argument("--final-csv", required=True)
    parser.add_argument("--structure-metrics-csv", required=True)
    parser.add_argument("--data-dir", default="")
    parser.add_argument("--exclusion-list-csv", default="")
    parser.add_argument("--out-json", required=True)
    parser.add_argument("--out-md", required=True)
    parser.add_argument("--fail-on-error", action="store_true")
    args = parser.parse_args()

    errors: list[str] = []
    warnings: list[str] = []
    checks: list[str] = []

    submission_path = Path(args.submission_csv)
    final_path = Path(args.final_csv)
    metrics_path = Path(args.structure_metrics_csv)
    exclusion_path = Path(args.exclusion_list_csv) if args.exclusion_list_csv else Path(args.data_dir) / "Exclusion_List.csv"

    report: dict[str, Any] = {
        "submission_csv": str(submission_path),
        "final_csv": str(final_path),
        "structure_metrics_csv": str(metrics_path),
        "exclusion_list_csv": str(exclusion_path),
        "passed": False,
        "errors": errors,
        "warnings": warnings,
        "checks": checks,
    }

    try:
        submission = read_csv(submission_path, "Submission CSV")
        final = read_csv(final_path, "Final CSV")
        metrics = read_csv(metrics_path, "Structure metrics CSV")
        exclusion = load_exclusion(exclusion_path)
        report["submission_sha256"] = sha256_file(submission_path)
        report["final_sha256"] = sha256_file(final_path)
        report["structure_metrics_sha256"] = sha256_file(metrics_path)
        report["exclusion_list_sha256"] = sha256_file(exclusion_path)
        checks.append(f"Loaded exclusion list with {len(exclusion)} entries.")
    except Exception as exc:
        errors.append(str(exc))
        submission = pd.DataFrame()
        final = pd.DataFrame()
        metrics = pd.DataFrame()
        exclusion = set()

    required_submission = {"Team_Name", "Seq_ID", "Sequence"}
    required_final = {"Seq_ID", "sequence"}
    required_metrics = {
        "sequence_sha1",
        "pdb_path",
        "has_real_structure",
        "structure_status",
        "structure_pass",
    }
    if not submission.empty and not required_submission.issubset(submission.columns):
        errors.append(f"Submission is missing required columns: {sorted(required_submission.difference(submission.columns))}")
    if not final.empty and not required_final.issubset(final.columns):
        errors.append(f"Final table is missing required columns: {sorted(required_final.difference(final.columns))}")
    if not metrics.empty and not required_metrics.issubset(metrics.columns):
        errors.append(f"Structure metrics are missing required columns: {sorted(required_metrics.difference(metrics.columns))}")

    if not submission.empty:
        if len(submission) != 6:
            errors.append(f"Submission has {len(submission)} rows; exactly 6 are required.")
        sequences = submission.get("Sequence", pd.Series(dtype=str)).astype(str).str.strip().str.upper()
        if sequences.duplicated().any():
            errors.append("Submission contains duplicate sequences.")
        for idx, seq in sequences.items():
            if not seq.startswith("M"):
                errors.append(f"Submission row {idx} does not start with M.")
            if not 220 <= len(seq) <= 250:
                errors.append(f"Submission row {idx} length={len(seq)} is outside [220, 250].")
            if not AA_RE.fullmatch(seq):
                errors.append(f"Submission row {idx} contains non-standard amino-acid characters.")
            if seq in exclusion:
                errors.append(f"Submission row {idx} exactly matches Exclusion_List.")
        checks.append("Validated sequence format, length, uniqueness, and exclusion membership.")

    if not submission.empty and not final.empty and required_submission.issubset(submission.columns) and required_final.issubset(final.columns):
        final_pairs = list(
            zip(
                final["Seq_ID"].astype(str).tolist(),
                final["sequence"].astype(str).str.strip().str.upper().tolist(),
            )
        )
        submission_pairs = list(
            zip(
                submission["Seq_ID"].astype(str).tolist(),
                submission["Sequence"].astype(str).str.strip().str.upper().tolist(),
            )
        )
        if final_pairs != submission_pairs:
            errors.append("Submission Seq_ID/Sequence rows do not exactly match final_top6_v72.csv in order.")
        else:
            checks.append("Submission and final Top-6 rows match exactly.")

    if not final.empty and not metrics.empty and required_final.issubset(final.columns) and required_metrics.issubset(metrics.columns):
        metrics = metrics.copy()
        metrics["sequence_sha1"] = metrics["sequence_sha1"].astype(str)
        metrics = metrics.sort_values("sequence_sha1").drop_duplicates("sequence_sha1", keep="first")
        metric_by_hash = metrics.set_index("sequence_sha1", drop=False)
        missing_structure = 0
        for _, row in final.iterrows():
            seq = str(row["sequence"]).strip().upper()
            sha = sequence_sha1(seq)
            if sha not in metric_by_hash.index:
                errors.append(f"Final sequence {row['Seq_ID']} has no metrics row for sequence_sha1={sha}.")
                missing_structure += 1
                continue
            m = metric_by_hash.loc[sha]
            if not as_bool(m.get("has_real_structure", False)):
                errors.append(f"Final sequence {row['Seq_ID']} is not marked has_real_structure=true.")
            if str(m.get("structure_status", "")) != "rmsd_checked":
                errors.append(f"Final sequence {row['Seq_ID']} has structure_status={m.get('structure_status')!r}, expected rmsd_checked.")
            if not as_bool(m.get("structure_pass", False)):
                errors.append(f"Final sequence {row['Seq_ID']} has structure_pass!=true.")
            pdb_path = Path(str(m.get("pdb_path", "") or ""))
            if not pdb_path.is_file():
                errors.append(f"Final sequence {row['Seq_ID']} PDB path is missing: {pdb_path}")
                continue
            if "pdb_sha256_v72" in metrics.columns:
                expected_hash = str(m.get("pdb_sha256_v72", "") or "")
                observed_hash = sha256_file(pdb_path)
                if not expected_hash:
                    errors.append(f"Final sequence {row['Seq_ID']} lacks pdb_sha256_v72.")
                elif expected_hash != observed_hash:
                    errors.append(f"Final sequence {row['Seq_ID']} PDB SHA256 does not match metrics.")
        if missing_structure == 0:
            checks.append("Verified final structure rows, PDB files, and V7.2 PDB hashes.")

    report["submission_rows"] = int(len(submission))
    report["final_rows"] = int(len(final))
    report["metrics_rows"] = int(len(metrics))
    report["errors"] = errors
    report["warnings"] = warnings
    report["checks"] = checks
    report["passed"] = not errors

    out_json = Path(args.out_json)
    out_md = Path(args.out_md)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_md.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    write_markdown(report, out_md)
    print(json.dumps(report, indent=2, ensure_ascii=False))
    if args.fail_on_error and errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
