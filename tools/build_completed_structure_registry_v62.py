#!/usr/bin/env python3
"""Build a completed-structure registry from V5.1/V6/V6.1/V6.2 output roots.

The scanner is intentionally conservative: a PDB is considered reusable only when
it can be attached to an exact amino-acid sequence through a candidate/priority/final
CSV or an existing structure_metrics.csv row that carries both sequence and pdb_path.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path
from typing import Any

import pandas as pd

PDB_SUFFIXES = {".pdb", ".ent"}
CANDIDATE_CSV_PATTERNS = [
    "candidate_universe*.csv",
    "all_ranked_candidates*.csv",
    "structure_priority*.csv",
    "ranked_top*.csv",
    "final_top*.csv",
    "selected_candidates*.csv",
    "submission*.csv",
]
METRICS_PATTERNS = ["structure_metrics*.csv", "missing_*structure_metrics*.csv"]


def sequence_sha1(seq: str) -> str:
    return hashlib.sha1(str(seq).strip().upper().encode("utf-8")).hexdigest()


def normalize_id(value: Any) -> str:
    x = str(value or "").strip().split("|", 1)[0]
    x = re.sub(r"_relaxed_rank_.*$", "", x)
    x = re.sub(r"_unrelaxed_rank_.*$", "", x)
    x = re.sub(r"_rank_.*$", "", x)
    x = re.sub(r"_model_.*$", "", x)
    x = re.sub(r"_pred_.*$", "", x)
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", x).strip("_")


def split_roots(text: str) -> list[Path]:
    roots = []
    for item in re.split(r"[:;,]", text or ""):
        item = item.strip()
        if item:
            roots.append(Path(item).expanduser())
    return roots


def is_sequence(seq: str) -> bool:
    return bool(re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]{50,400}", str(seq).strip().upper()))


def collect_candidate_maps(roots: list[Path]) -> tuple[dict[str, dict[str, str]], dict[str, dict[str, str]], list[dict[str, Any]]]:
    id_map: dict[str, dict[str, str]] = {}
    hash_map: dict[str, dict[str, str]] = {}
    source_rows: list[dict[str, Any]] = []
    for root in roots:
        if not root.exists():
            continue
        for pat in CANDIDATE_CSV_PATTERNS:
            for csv_path in root.rglob(pat):
                try:
                    df = pd.read_csv(csv_path)
                except Exception:
                    continue
                seq_col = "sequence" if "sequence" in df.columns else ("Sequence" if "Sequence" in df.columns else None)
                if seq_col is None:
                    continue
                for _, row in df.iterrows():
                    seq = str(row.get(seq_col, "")).strip().upper()
                    if not is_sequence(seq):
                        continue
                    sha = str(row.get("sequence_sha1") or sequence_sha1(seq))
                    payload = {
                        "sequence": seq,
                        "sequence_sha1": sha,
                        "sequence_length": str(len(seq)),
                        "source_csv": str(csv_path),
                        "root": str(root),
                    }
                    hash_map.setdefault(sha, payload)
                    for col in ["candidate_id", "Seq_ID", "seq_id", "id"]:
                        key = normalize_id(row.get(col, ""))
                        if key:
                            id_map[key] = payload
                    source_rows.append({**payload, "csv_type": pat})
    return id_map, hash_map, source_rows


def collect_metrics_rows(roots: list[Path]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for root in roots:
        if not root.exists():
            continue
        for pat in METRICS_PATTERNS:
            for csv_path in root.rglob(pat):
                try:
                    df = pd.read_csv(csv_path)
                except Exception:
                    continue
                if "sequence" not in df.columns:
                    continue
                for _, row in df.iterrows():
                    seq = str(row.get("sequence", "")).strip().upper()
                    if not is_sequence(seq):
                        continue
                    pdb_path = str(row.get("pdb_path", "") or "")
                    status = str(row.get("structure_status", "") or "")
                    passed = str(row.get("structure_pass", "") or "").lower() == "true"
                    sha = str(row.get("sequence_sha1") or sequence_sha1(seq))
                    rows.append({
                        "sequence": seq,
                        "sequence_sha1": sha,
                        "sequence_length": len(seq),
                        "candidate_id": row.get("candidate_id", ""),
                        "Seq_ID": row.get("Seq_ID", ""),
                        "pdb_path": pdb_path,
                        "structure_status": status,
                        "structure_pass": passed,
                        "metrics_csv": str(csv_path),
                        "root": str(root),
                    })
    return rows


def discover_pdbs(roots: list[Path], id_map: dict[str, dict[str, str]], hash_map: dict[str, dict[str, str]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    seen: set[tuple[str, str]] = set()
    for root in roots:
        if not root.exists():
            continue
        for pdb in sorted(root.rglob("*")):
            if not pdb.is_file() or pdb.suffix.lower() not in PDB_SUFFIXES:
                continue
            matched_payload = None
            matched_id = ""
            candidates = [normalize_id(pdb.stem)]
            candidates.extend(normalize_id(part.name) for part in pdb.parents if part != root.parent)
            for key in candidates:
                if key in id_map:
                    matched_payload = id_map[key]
                    matched_id = key
                    break
            if matched_payload is None:
                continue
            sha = matched_payload["sequence_sha1"]
            key = (sha, str(pdb.resolve()))
            if key in seen:
                continue
            seen.add(key)
            rows.append({
                "sequence_sha1": sha,
                "sequence_length": matched_payload["sequence_length"],
                "sequence": matched_payload["sequence"],
                "matched_id": matched_id,
                "pdb_path": str(pdb.resolve()),
                "source_csv": matched_payload.get("source_csv", ""),
                "search_root": str(root),
                "registry_source": "pdb_scan_id_match",
            })
    return rows


def main() -> None:
    p = argparse.ArgumentParser(description="Find completed reusable structures from historical GFP runs.")
    p.add_argument("--search-roots", required=True, help="Colon/comma/semicolon separated output roots to scan.")
    p.add_argument("--out-csv", required=True)
    p.add_argument("--summary-json", default="")
    args = p.parse_args()

    roots = split_roots(args.search_roots)
    id_map, hash_map, source_rows = collect_candidate_maps(roots)
    pdb_rows = discover_pdbs(roots, id_map, hash_map)
    metrics_rows = collect_metrics_rows(roots)

    rows = pdb_rows[:]
    for m in metrics_rows:
        pdb_path = m.get("pdb_path") or ""
        if pdb_path and Path(pdb_path).exists():
            rows.append({
                "sequence_sha1": m["sequence_sha1"],
                "sequence_length": m["sequence_length"],
                "sequence": m["sequence"],
                "matched_id": normalize_id(m.get("candidate_id") or m.get("Seq_ID") or ""),
                "pdb_path": str(Path(pdb_path).resolve()),
                "source_csv": m.get("metrics_csv", ""),
                "search_root": m.get("root", ""),
                "registry_source": "structure_metrics_pdb_path",
                "structure_status": m.get("structure_status", ""),
                "structure_pass": m.get("structure_pass", False),
            })

    reg = pd.DataFrame(rows)
    if not reg.empty:
        reg = reg.drop_duplicates(["sequence_sha1", "pdb_path"]).sort_values(["sequence_sha1", "pdb_path"])
    Path(args.out_csv).parent.mkdir(parents=True, exist_ok=True)
    reg.to_csv(args.out_csv, index=False)

    summary = {
        "search_roots": [str(r) for r in roots],
        "candidate_csv_sequence_rows": len(source_rows),
        "id_map_entries": len(id_map),
        "unique_sequence_hashes": len(hash_map),
        "reusable_pdb_rows": int(len(reg)),
        "unique_reusable_sequences": int(reg["sequence_sha1"].nunique()) if not reg.empty else 0,
        "out_csv": str(args.out_csv),
    }
    if args.summary_json:
        Path(args.summary_json).parent.mkdir(parents=True, exist_ok=True)
        Path(args.summary_json).write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
