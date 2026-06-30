#!/usr/bin/env python3
"""Traceable ProteinMPNN score_only adapter for SynBio GFP V6.

V6 runs ProteinMPNN in directory mode because many ProteinMPNN helper versions expect
--input_path to be a PDB directory rather than a single PDB file. The output contract is
always proteinmpnn_scores.csv with candidate_id alignment for final reranking.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import shlex
import subprocess
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


def safe_name(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("_")


def normalize_id(value: str) -> str:
    x = str(value).strip().split("|", 1)[0]
    x = re.sub(r"_relaxed_rank_.*$", "", x)
    x = re.sub(r"_unrelaxed_rank_.*$", "", x)
    x = re.sub(r"_rank_.*$", "", x)
    x = re.sub(r"_model_.*$", "", x)
    return safe_name(x)


def run_cmd(cmd: list[str], log_path: Path) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    command_line = "\n$ " + " ".join(shlex.quote(str(x)) for x in cmd)
    print(command_line, flush=True)
    with log_path.open("a", encoding="utf-8", buffering=1) as log:
        log.write(command_line + "\n")
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1, universal_newlines=True)
        assert proc.stdout is not None
        for line in proc.stdout:
            log.write(line)
            log.flush()
            print(line, end="", flush=True)
        code = proc.wait()
        log.write(f"[exit_code] {code}\n")
        print(f"[exit_code] {code}", flush=True)
        if code != 0:
            raise subprocess.CalledProcessError(code, cmd)


def load_priority(priority_csv: Path | None) -> dict[str, dict[str, Any]]:
    if not priority_csv or not priority_csv.exists():
        return {}
    rows: dict[str, dict[str, Any]] = {}
    with priority_csv.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            for key in ("candidate_id", "Seq_ID", "id", "seq_id"):
                if row.get(key):
                    rows[normalize_id(row[key])] = row
                    break
    return rows


def extract_score(npz: Path) -> tuple[float, float, str]:
    z = np.load(npz, allow_pickle=True)
    if "global_score" in z:
        arr = z["global_score"]
        key = "global_score"
    elif "score" in z:
        arr = z["score"]
        key = "score"
    else:
        return math.nan, math.nan, "missing_score_key"
    return float(np.mean(arr)), float(np.std(arr)), f"score_only_imported_npz:{key}"


def main() -> None:
    p = argparse.ArgumentParser(description="Run ProteinMPNN score_only for V6 PDBs and emit a traceable CSV.")
    p.add_argument("--pdb-dir", required=True)
    p.add_argument("--out-csv", required=True)
    p.add_argument("--work-dir", required=True)
    p.add_argument("--proteinmpnn-dir", default=os.environ.get("PROTEINMPNN_DIR", "/hyperai/home/tools/ProteinMPNN"))
    p.add_argument("--priority-csv", default=None)
    p.add_argument("--chain-id", default="A")
    p.add_argument("--python-bin", default=sys.executable)
    p.add_argument("--num-seq-per-target", type=int, default=1)
    p.add_argument("--sampling-temp", type=float, default=0.1)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--skip-existing", action="store_true", default=False)
    args = p.parse_args()

    pdb_dir = Path(args.pdb_dir)
    out_csv = Path(args.out_csv)
    work_dir = Path(args.work_dir)
    proteinmpnn_dir = Path(args.proteinmpnn_dir)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    work_dir.mkdir(parents=True, exist_ok=True)
    if not pdb_dir.exists():
        raise FileNotFoundError(f"PDB directory does not exist: {pdb_dir}")
    if not proteinmpnn_dir.exists():
        raise FileNotFoundError(f"ProteinMPNN directory does not exist: {proteinmpnn_dir}")

    pdbs = sorted(pdb_dir.rglob("*.pdb"))
    if not pdbs:
        raise FileNotFoundError(f"No PDB files found under {pdb_dir}")
    priority = load_priority(Path(args.priority_csv) if args.priority_csv else None)
    parsed = work_dir / "parsed_pdbs.jsonl"
    assigned = work_dir / "assigned_pdbs.jsonl"
    log_path = work_dir / "proteinmpnn_score_only.log"

    if not args.skip_existing or not list(work_dir.rglob("*.npz")):
        # Parse the whole PDB directory once; this avoids the single-file parser incompatibility seen in V5 handoff runs.
        run_cmd([args.python_bin, str(proteinmpnn_dir / "helper_scripts" / "parse_multiple_chains.py"), "--input_path", str(pdb_dir), "--output_path", str(parsed)], log_path)
        run_cmd([args.python_bin, str(proteinmpnn_dir / "helper_scripts" / "assign_fixed_chains.py"), "--input_path", str(parsed), "--output_path", str(assigned), "--chain_list", args.chain_id if args.chain_id != "auto" else "A"], log_path)
        run_cmd([args.python_bin, str(proteinmpnn_dir / "protein_mpnn_run.py"), "--jsonl_path", str(parsed), "--chain_id_jsonl", str(assigned), "--out_folder", str(work_dir), "--num_seq_per_target", str(args.num_seq_per_target), "--sampling_temp", str(args.sampling_temp), "--score_only", "1", "--seed", str(args.seed)], log_path)

    npz_files = sorted(work_dir.rglob("*.npz"))
    rows: list[dict[str, Any]] = []
    seen: set[str] = set()
    for npz in npz_files:
        m = re.search(r"(cand_\d+)", str(npz))
        cid = normalize_id(m.group(1) if m else npz.stem)
        meta = priority.get(cid, {})
        candidate_id = meta.get("candidate_id", cid)
        score, score_std, status = extract_score(npz)
        rows.append({"candidate_id": candidate_id, "pdb_id": cid, "proteinmpnn_score": score, "proteinmpnn_score_std": score_std, "proteinmpnn_status": status, "proteinmpnn_npz": str(npz)})
        seen.add(str(candidate_id))
    # Keep missing rows explicit so final diagnostics can distinguish missing scores from valid low scores.
    for cid in priority:
        candidate_id = priority[cid].get("candidate_id", cid)
        if str(candidate_id) not in seen:
            rows.append({"candidate_id": candidate_id, "pdb_id": cid, "proteinmpnn_score": math.nan, "proteinmpnn_score_std": math.nan, "proteinmpnn_status": "missing_score_npz", "proteinmpnn_npz": ""})
    df = pd.DataFrame(rows).drop_duplicates("candidate_id", keep="first")
    df.to_csv(out_csv, index=False)
    (work_dir / "proteinmpnn_summary.json").write_text(json.dumps({"pdb_count": len(pdbs), "npz_count": len(npz_files), "rows": len(df), "out_csv": str(out_csv)}, indent=2), encoding="utf-8")
    print(f"Wrote {len(df)} rows to {out_csv}")
    print(df.head())


if __name__ == "__main__":
    main()
