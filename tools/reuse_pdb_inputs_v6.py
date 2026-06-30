#!/usr/bin/env python3
"""Normalize external ColabFold/AF PDB outputs and build V6 structure_metrics.csv.

This adapter is intentionally narrow: it never generates new structures. It only imports
existing PDBs through a stable candidate_id contract so V6 can be handed off without
binding the pipeline to one previous run directory.
"""
from __future__ import annotations

import argparse
import re
import shutil
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from synbio_gfp_v6.structure_filter import compute_pdb_metrics

PDB_SUFFIXES = {".pdb", ".ent"}


def normalize_id(value: str) -> str:
    x = str(value).strip().split("|", 1)[0]
    x = re.sub(r"_relaxed_rank_.*$", "", x)
    x = re.sub(r"_unrelaxed_rank_.*$", "", x)
    x = re.sub(r"_rank_.*$", "", x)
    x = re.sub(r"_model_.*$", "", x)
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", x).strip("_")


def load_priority(path: Path | None) -> list[str]:
    if not path or not path.exists():
        return []
    df = pd.read_csv(path)
    for col in ["candidate_id", "Seq_ID", "seq_id", "id"]:
        if col in df.columns:
            return [normalize_id(x) for x in df[col].astype(str).tolist()]
    return []


def discover_pdbs(pdb_dir: Path) -> dict[str, Path]:
    index: dict[str, Path] = {}
    for path in sorted(pdb_dir.rglob("*")):
        if not path.is_file() or path.suffix.lower() not in PDB_SUFFIXES:
            continue
        keys = [normalize_id(path.stem)]
        keys.extend(normalize_id(part.name) for part in path.parents if part != pdb_dir.parent)
        for key in keys:
            if key and key not in index:
                index[key] = path
    return index


def match_pdb(candidate_id: str, pdb_index: dict[str, Path]) -> Path | None:
    cid = normalize_id(candidate_id)
    if cid in pdb_index:
        return pdb_index[cid]
    for key, path in pdb_index.items():
        if key.startswith(cid) or cid.startswith(key):
            return path
    return None


def link_or_copy(src: Path, dst: Path, copy: bool = False) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists() or dst.is_symlink():
        dst.unlink()
    # Symlinks keep large PDB reuse cheap; --copy-pdbs is available for portable handoff bundles.
    if copy:
        shutil.copy2(src, dst)
    else:
        dst.symlink_to(src.resolve())


def main() -> None:
    p = argparse.ArgumentParser(description="Reuse external PDBs and emit V6 structure metrics.")
    p.add_argument("--input-pdb-dir", required=True)
    p.add_argument("--priority-csv", default=None)
    p.add_argument("--reference-pdb", required=True)
    p.add_argument("--normalized-pdb-dir", required=True)
    p.add_argument("--metrics-csv", required=True)
    p.add_argument("--copy-pdbs", action="store_true")
    args = p.parse_args()

    input_dir = Path(args.input_pdb_dir)
    ref_pdb = Path(args.reference_pdb)
    out_pdb_dir = Path(args.normalized_pdb_dir)
    metrics_csv = Path(args.metrics_csv)
    if not input_dir.exists():
        raise FileNotFoundError(f"Input PDB directory does not exist: {input_dir}")
    if not ref_pdb.exists():
        raise FileNotFoundError(f"Reference PDB does not exist: {ref_pdb}")

    pdb_index = discover_pdbs(input_dir)
    priority_ids = load_priority(Path(args.priority_csv) if args.priority_csv else None)
    if not priority_ids:
        priority_ids = sorted(pdb_index)

    rows: list[dict[str, Any]] = []
    for cid in priority_ids:
        pdb = match_pdb(cid, pdb_index)
        if pdb is None:
            rows.append({"candidate_id": cid, "structure_status": "missing_reuse_pdb", "structure_pass": False, "pdb_path": ""})
            continue
        normalized = out_pdb_dir / f"{cid}.pdb"
        link_or_copy(pdb, normalized, copy=args.copy_pdbs)
        try:
            metrics = compute_pdb_metrics(ref_pdb, normalized)
            rows.append({"candidate_id": cid, **metrics, "structure_pass": True, "pdb_path": str(normalized)})
        except Exception as exc:
            rows.append({"candidate_id": cid, "barrel_rmsd": np.nan, "pocket_rmsd": np.nan, "ion_network_rmsd": np.nan, "structure_status": f"metric_failed:{type(exc).__name__}", "structure_pass": False, "pdb_path": str(normalized)})

    metrics_csv.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows)
    df.to_csv(metrics_csv, index=False)
    print(f"Imported PDBs: {sum(df['structure_pass'].fillna(False).astype(bool))}/{len(df)}")
    print(f"Wrote structure metrics: {metrics_csv}")


if __name__ == "__main__":
    main()
