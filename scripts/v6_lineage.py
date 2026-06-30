#!/usr/bin/env python3
"""Create file-level lineage tables for SynBio GFP V6 runs."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
from typing import Iterable

TRACK_SUFFIXES = {".csv", ".json", ".jsonl", ".fasta", ".fa", ".faa", ".pdb", ".log", ".md", ".yaml", ".yml", ".txt"}


def sha256_file(path: Path, block_size: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(block_size), b""):
            h.update(chunk)
    return h.hexdigest()


def infer_stage(path: Path, run_dir: Path) -> str:
    rel = path.relative_to(run_dir)
    first = rel.parts[0] if rel.parts else "root"
    mapping = {
        "00_preflight": "preflight",
        "01_stage1": "stage1_candidate_generation",
        "02_structure": "structure_prediction_and_metrics",
        "03_proteinmpnn": "proteinmpnn_score_only",
        "04_final": "final_reranking",
        "lineage": "lineage",
        "logs": "logs",
    }
    return mapping.get(first, first)


def iter_files(run_dir: Path) -> Iterable[Path]:
    for path in sorted(run_dir.rglob("*")):
        if not path.is_file():
            continue
        if path.suffix.lower() in TRACK_SUFFIXES or path.name.startswith("v5_") or path.name.startswith("v6_"):
            yield path


def main() -> None:
    p = argparse.ArgumentParser(description="Index V5 run artifacts with hashes and stage labels.")
    p.add_argument("--run-dir", required=True)
    p.add_argument("--out-dir", default=None)
    args = p.parse_args()
    run_dir = Path(args.run_dir).resolve()
    out_dir = Path(args.out_dir).resolve() if args.out_dir else run_dir / "lineage"
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for path in iter_files(run_dir):
        try:
            rel = path.relative_to(run_dir)
        except ValueError:
            rel = path
        rows.append({
            "stage": infer_stage(path, run_dir),
            "relative_path": str(rel),
            "absolute_path": str(path),
            "size_bytes": path.stat().st_size,
            "mtime_epoch": int(path.stat().st_mtime),
            "sha256": sha256_file(path),
        })

    csv_path = out_dir / "v6_artifact_index.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["stage", "relative_path", "absolute_path", "size_bytes", "mtime_epoch", "sha256"])
        writer.writeheader()
        writer.writerows(rows)

    summary = {
        "run_dir": str(run_dir),
        "artifact_count": len(rows),
        "stages": sorted({r["stage"] for r in rows}),
        "artifact_index_csv": str(csv_path),
    }
    (out_dir / "v6_lineage_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"Indexed {len(rows)} artifacts under {run_dir}")
    print(csv_path)


if __name__ == "__main__":
    main()
