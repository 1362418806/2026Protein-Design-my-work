#!/usr/bin/env python3
"""V6.1 adaptive structure resolver with automatic PDB registry and reuse-first workflow.

PDB files may be reused only when the amino-acid sequence is an exact match. Metrics are
always recomputed for the current run so old thresholds or old references cannot leak into
final decisions.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PDB_SUFFIXES = {".pdb", ".ent"}
CANDIDATE_CSV_PATTERNS = [
    "all_ranked_candidates*.csv",
    "structure_priority*.csv",
    "ranked_top*.csv",
    "final_top*.csv",
    "selected_candidates*.csv",
]


def sequence_sha1(seq: str) -> str:
    return hashlib.sha1(str(seq).strip().upper().encode("utf-8")).hexdigest()


def normalize_id(value: Any) -> str:
    x = str(value or "").strip()
    x = x.split("|", 1)[0]
    x = re.sub(r"_relaxed_rank_.*$", "", x)
    x = re.sub(r"_unrelaxed_rank_.*$", "", x)
    x = re.sub(r"_rank_.*$", "", x)
    x = re.sub(r"_model_.*$", "", x)
    x = re.sub(r"_pred_.*$", "", x)
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", x).strip("_")


def split_roots(text: str) -> list[Path]:
    out: list[Path] = []
    for item in re.split(r"[:;,]", text or ""):
        item = item.strip()
        if item:
            out.append(Path(item).expanduser())
    return out


def load_priority(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    if "sequence" not in df.columns:
        raise ValueError(f"Priority CSV must contain a sequence column: {path}")
    df = df.copy()
    if "candidate_id" not in df.columns:
        df["candidate_id"] = [f"cand_{i:06d}" for i in range(len(df))]
    if "Seq_ID" not in df.columns:
        df["Seq_ID"] = df["candidate_id"]
    df["sequence"] = df["sequence"].astype(str).str.strip().str.upper()
    df["sequence_sha1"] = df["sequence"].map(sequence_sha1)
    df["sequence_length"] = df["sequence"].str.len()
    return df


def collect_sequence_ids(search_roots: list[Path]) -> dict[str, dict[str, str]]:
    id_map: dict[str, dict[str, str]] = {}
    for root in search_roots:
        if not root.exists():
            continue
        for pat in CANDIDATE_CSV_PATTERNS:
            for csv_path in root.rglob(pat):
                try:
                    df = pd.read_csv(csv_path, usecols=lambda c: c in {"candidate_id", "Seq_ID", "seq_id", "id", "sequence"})
                except Exception:
                    continue
                if "sequence" not in df.columns:
                    continue
                for _, row in df.iterrows():
                    seq = str(row.get("sequence", "")).strip().upper()
                    if not re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]{50,400}", seq):
                        continue
                    payload = {"sequence": seq, "sequence_sha1": sequence_sha1(seq), "sequence_length": str(len(seq)), "source_csv": str(csv_path)}
                    for col in ["candidate_id", "Seq_ID", "seq_id", "id"]:
                        key = normalize_id(row.get(col, ""))
                        if key:
                            id_map[key] = payload
    return id_map


def discover_pdb_registry(search_roots: list[Path], id_map: dict[str, dict[str, str]]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    seen: set[tuple[str, str]] = set()
    for root in search_roots:
        if not root.exists():
            continue
        for pdb in sorted(root.rglob("*")):
            if not pdb.is_file() or pdb.suffix.lower() not in PDB_SUFFIXES:
                continue
            id_candidates = [normalize_id(pdb.stem)]
            id_candidates.extend(normalize_id(part.name) for part in pdb.parents if part != root.parent)
            matched_id = ""
            matched_payload: dict[str, str] | None = None
            for cand_id in id_candidates:
                if cand_id in id_map:
                    matched_id = cand_id
                    matched_payload = id_map[cand_id]
                    break
            if matched_payload is None:
                continue
            key = (matched_payload["sequence_sha1"], str(pdb.resolve()))
            if key in seen:
                continue
            seen.add(key)
            rows.append({
                "sequence_sha1": matched_payload["sequence_sha1"],
                "sequence_length": matched_payload["sequence_length"],
                "matched_id": matched_id,
                "pdb_path": str(pdb.resolve()),
                "source_csv": matched_payload.get("source_csv", ""),
                "search_root": str(root),
            })
    return pd.DataFrame(rows)


def link_or_copy(src: Path, dst: Path, copy_pdbs: bool = False) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists() or dst.is_symlink():
        dst.unlink()
    # Symlinks keep reused PDBs cheap; --copy-pdbs makes a portable run directory when needed.
    if copy_pdbs:
        shutil.copy2(src, dst)
    else:
        dst.symlink_to(src.resolve())


def write_fasta(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for _, row in df.iterrows():
            cid = normalize_id(row.get("candidate_id", row.get("Seq_ID", "candidate")))
            seq = str(row["sequence"]).strip().upper()
            handle.write(f">{cid}|sequence_sha1={row['sequence_sha1']}\n")
            for i in range(0, len(seq), 80):
                handle.write(seq[i:i + 80] + "\n")


def run_cmd(cmd: list[str], log_path: Path, dry_run: bool = False) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    print("$ " + " ".join(shlex.quote(str(x)) for x in cmd), flush=True)
    with log_path.open("a", encoding="utf-8", buffering=1) as log:
        log.write("$ " + " ".join(shlex.quote(str(x)) for x in cmd) + "\n")
        if dry_run:
            log.write("[dry_run] skipped\n")
            return
        # Stream backend logs line-by-line so SSH interruption does not hide progress.
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)
        assert proc.stdout is not None
        for line in proc.stdout:
            print(line, end="", flush=True)
            log.write(line)
        code = proc.wait()
        log.write(f"[exit_code] {code}\n")
        if code != 0:
            raise RuntimeError(f"Command failed with exit code {code}; see {log_path}")


def main() -> None:
    p = argparse.ArgumentParser(description="Resolve V6.1 structures by exact-sequence PDB reuse plus missing prediction.")
    p.add_argument("--priority-csv", required=True)
    p.add_argument("--reference-pdb", required=True)
    p.add_argument("--prediction-dir", required=True)
    p.add_argument("--metrics-csv", required=True)
    p.add_argument("--work-dir", required=True)
    p.add_argument("--pdb-search-roots", default=os.environ.get("PDB_SEARCH_ROOTS", ""))
    p.add_argument("--structure-policy", choices=["reuse-only", "reuse-first", "predict-missing", "metrics-only"], default=os.environ.get("STRUCTURE_POLICY", "reuse-first"))
    p.add_argument("--runner", choices=["auto", "colabfold", "external", "metrics-only"], default=os.environ.get("STRUCTURE_RUNNER", "auto"))
    p.add_argument("--colabfold-extra-args", default=os.environ.get("COLABFOLD_EXTRA_ARGS", "--num-models 1 --num-recycle 1"))
    p.add_argument("--external-command-template", default=os.environ.get("EXTERNAL_STRUCTURE_COMMAND", ""))
    p.add_argument("--copy-pdbs", action="store_true")
    p.add_argument("--dry-run", action="store_true")
    args = p.parse_args()

    priority_csv = Path(args.priority_csv)
    prediction_dir = Path(args.prediction_dir)
    metrics_csv = Path(args.metrics_csv)
    work_dir = Path(args.work_dir)
    prediction_dir.mkdir(parents=True, exist_ok=True)
    work_dir.mkdir(parents=True, exist_ok=True)

    priority = load_priority(priority_csv)
    roots = split_roots(args.pdb_search_roots)
    roots.extend([priority_csv.parents[1] if len(priority_csv.parents) > 1 else priority_csv.parent, prediction_dir])
    roots = [r for i, r in enumerate(roots) if r and r.exists() and r not in roots[:i]]

    id_map = collect_sequence_ids(roots)
    # Current priority candidates are authoritative for this run, so they override historical ID mappings.
    for _, row in priority.iterrows():
        payload = {"sequence": row["sequence"], "sequence_sha1": row["sequence_sha1"], "sequence_length": str(row["sequence_length"]), "source_csv": str(priority_csv)}
        for col in ["candidate_id", "Seq_ID"]:
            key = normalize_id(row.get(col, ""))
            if key:
                id_map[key] = payload

    registry = discover_pdb_registry(roots, id_map)
    registry_path = work_dir / "pdb_registry.csv"
    registry.to_csv(registry_path, index=False)

    registry_by_hash: dict[str, Path] = {}
    if not registry.empty:
        for _, row in registry.iterrows():
            registry_by_hash.setdefault(str(row["sequence_sha1"]), Path(row["pdb_path"]))

    reused_rows: list[dict[str, Any]] = []
    missing_rows: list[dict[str, Any]] = []
    for _, row in priority.iterrows():
        cid = normalize_id(row.get("candidate_id", row.get("Seq_ID", "candidate")))
        seq_hash = str(row["sequence_sha1"])
        src_pdb = registry_by_hash.get(seq_hash)
        normalized = prediction_dir / f"{cid}.pdb"
        if src_pdb and src_pdb.exists():
            link_or_copy(src_pdb, normalized, copy_pdbs=args.copy_pdbs)
            reused_rows.append({**row.to_dict(), "reused_pdb_path": str(normalized), "source_pdb_path": str(src_pdb)})
        else:
            missing_rows.append(row.to_dict())

    reused_manifest = work_dir / "reused_pdb_manifest.csv"
    missing_csv = work_dir / "missing_pdb_candidates.csv"
    missing_fasta = work_dir / "colabfold_missing_tasks.fasta"
    pd.DataFrame(reused_rows).to_csv(reused_manifest, index=False)
    missing_df = pd.DataFrame(missing_rows)
    missing_df.to_csv(missing_csv, index=False)
    if not missing_df.empty:
        write_fasta(missing_df, missing_fasta)
    else:
        missing_fasta.write_text("", encoding="utf-8")

    if args.structure_policy in {"reuse-first", "predict-missing"} and not missing_df.empty:
        runner = args.runner
        if runner == "metrics-only":
            print("Missing PDBs detected, but runner=metrics-only; no predictions will be generated.")
        else:
            backend_cmd = [
                sys.executable, str(ROOT / "tools" / "run_structure_backend_v6.py"),
                "--runner", runner,
                "--fasta-paths", str(missing_fasta),
                "--af-output-dir", str(prediction_dir),
                "--priority-csv", str(missing_csv),
                "--metrics-csv", str(work_dir / "missing_structure_metrics.csv"),
                "--work-dir", str(work_dir / "adapter_predict_missing"),
                "--reference-pdb", str(args.reference_pdb),
                "--colabfold-extra-args", args.colabfold_extra_args,
            ]
            if args.external_command_template:
                backend_cmd += ["--external-command-template", args.external_command_template]
            run_cmd(backend_cmd, work_dir / "resolve_structures_predict_missing.log", dry_run=args.dry_run)
    elif args.structure_policy == "reuse-only" and not missing_df.empty:
        print(f"reuse-only policy: {len(missing_df)} missing structures will remain missing.")
    elif args.structure_policy == "metrics-only":
        print("metrics-only policy: only currently available PDBs will be scored.")

    full_fasta = work_dir / "all_structure_tasks.fasta"
    write_fasta(priority, full_fasta)
    metric_cmd = [
        sys.executable, str(ROOT / "tools" / "run_structure_backend_v6.py"),
        "--runner", "metrics-only",
        "--fasta-paths", str(full_fasta),
        "--af-output-dir", str(prediction_dir),
        "--priority-csv", str(priority_csv),
        "--metrics-csv", str(metrics_csv),
        "--work-dir", str(work_dir / "adapter_recompute_metrics"),
        "--reference-pdb", str(args.reference_pdb),
    ]
    run_cmd(metric_cmd, work_dir / "resolve_structures_recompute_metrics.log", dry_run=args.dry_run)

    report = {
        "priority_candidates": int(len(priority)),
        "registry_rows": int(len(registry)),
        "reused_pdbs": int(len(reused_rows)),
        "missing_before_prediction": int(len(missing_rows)),
        "structure_policy": args.structure_policy,
        "metrics_csv": str(metrics_csv),
        "prediction_dir": str(prediction_dir),
    }
    (work_dir / "structure_reuse_report.md").write_text(
        "# V6.1 Structure Reuse Report\n\n" + "\n".join(f"- {k}: {v}" for k, v in report.items()) + "\n",
        encoding="utf-8",
    )
    (work_dir / "structure_reuse_report.json").write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
