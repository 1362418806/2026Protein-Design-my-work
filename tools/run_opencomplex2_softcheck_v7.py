#!/usr/bin/env python3
"""V7 orthogonal soft validator for GFP candidates.

This adapter is intentionally backend-agnostic. It can import an existing
OpenComplex2/AF-family metrics CSV, score an existing PDB directory, or invoke an
external command template that writes PDB files. The emitted scores are soft
validation evidence only; they must not replace the primary real-structure gate.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import shlex
import signal
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import pandas as pd

PDB_SUFFIXES = (".pdb", ".ent")
AA20_RE = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")
# 1-indexed sfGFP reference regions; these match the existing V6 structure backend contract.
BARREL_POSITIONS = list(range(1, 240))
POCKET_POSITIONS = [39, 63, 64, 65, 66, 67, 68, 69, 96, 145, 148, 203, 205, 222]
ION_NETWORK_POSITIONS = [17, 30, 32, 115, 122]


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


def load_table(path: str | Path) -> pd.DataFrame:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(p)
    return pd.read_excel(p) if p.suffix.lower() in {".xlsx", ".xls"} else pd.read_csv(p)


def ensure_candidate_keys(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if "sequence" not in out.columns:
        raise ValueError("Candidate CSV must contain a sequence column.")
    out["sequence"] = out["sequence"].astype(str).str.strip().str.upper()
    out = out[out["sequence"].map(lambda s: bool(AA20_RE.fullmatch(s)))].copy()
    if "sequence_sha1" not in out.columns:
        out["sequence_sha1"] = out["sequence"].map(sequence_sha1)
    if "candidate_id" not in out.columns:
        out["candidate_id"] = [f"cand_{i:06d}" for i in range(len(out))]
    if "Seq_ID" not in out.columns:
        out["Seq_ID"] = out["candidate_id"]
    out["candidate_id"] = out["candidate_id"].map(normalize_id)
    return out


def write_fasta(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for _, row in df.iterrows():
            cid = normalize_id(row.get("candidate_id", row.get("Seq_ID", "candidate")))
            seq = str(row["sequence"]).strip().upper()
            handle.write(f">{cid}|sequence_sha1={row['sequence_sha1']}\n")
            for i in range(0, len(seq), 80):
                handle.write(seq[i:i + 80] + "\n")


def parse_ca_atoms(path: Path, chain_id: str | None = None) -> dict[int, tuple[np.ndarray, float]]:
    atoms: dict[int, tuple[np.ndarray, float]] = {}
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith("ATOM") or line[12:16].strip() != "CA":
                continue
            chain = line[21].strip()
            if chain_id and chain != chain_id:
                continue
            try:
                resseq = int(line[22:26])
                xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])], dtype=np.float64)
                bfactor = float(line[60:66]) if len(line) >= 66 else math.nan
            except ValueError:
                continue
            atoms[resseq] = (xyz, bfactor)
    return atoms


def kabsch_rmsd(ref_xyz: np.ndarray, mob_xyz: np.ndarray) -> float:
    if ref_xyz.shape != mob_xyz.shape or ref_xyz.ndim != 2 or ref_xyz.shape[0] < 3:
        return math.nan
    ref0 = ref_xyz - ref_xyz.mean(axis=0)
    mob0 = mob_xyz - mob_xyz.mean(axis=0)
    cov = mob0.T @ ref0
    u, _, vt = np.linalg.svd(cov)
    d = np.sign(np.linalg.det(u @ vt))
    rot = u @ np.diag([1.0, 1.0, d]) @ vt
    aligned = mob0 @ rot
    return float(np.sqrt(np.mean(np.sum((aligned - ref0) ** 2, axis=1))))


def subset_rmsd(ref_atoms: dict[int, tuple[np.ndarray, float]], pred_atoms: dict[int, tuple[np.ndarray, float]], positions: Sequence[int]) -> float:
    common = [p for p in positions if p in ref_atoms and p in pred_atoms]
    if len(common) < 3:
        return math.nan
    ref_xyz = np.stack([ref_atoms[p][0] for p in common])
    pred_xyz = np.stack([pred_atoms[p][0] for p in common])
    return kabsch_rmsd(ref_xyz, pred_xyz)


def mean_plddt(pred_atoms: dict[int, tuple[np.ndarray, float]]) -> float:
    vals = [bf for _, bf in pred_atoms.values() if not math.isnan(bf)]
    return float(np.mean(vals)) if vals else math.nan


def discover_pdbs(pdb_dir: Path) -> dict[str, Path]:
    index: dict[str, Path] = {}
    if not pdb_dir.exists():
        return index
    for path in sorted(pdb_dir.rglob("*")):
        if not path.is_file() or path.suffix.lower() not in PDB_SUFFIXES:
            continue
        candidates = [path.stem]
        candidates.extend(parent.name for parent in path.parents if parent != pdb_dir.parent)
        for item in candidates:
            key = normalize_id(item)
            if key and key not in index:
                index[key] = path
    return index


def match_pdb(row: pd.Series, pdb_index: dict[str, Path]) -> Path | None:
    candidates = [row.get("candidate_id", ""), row.get("Seq_ID", ""), row.get("sequence_sha1", "")]
    for item in candidates:
        key = normalize_id(item)
        if key in pdb_index:
            return pdb_index[key]
    sha = str(row.get("sequence_sha1", ""))
    for key, path in pdb_index.items():
        if sha and sha in key:
            return path
    return None


def score_from_metrics(row: dict[str, Any], args: argparse.Namespace) -> tuple[float, bool, str]:
    barrel = float(row.get("oc2_barrel_rmsd", math.nan))
    pocket = float(row.get("oc2_pocket_rmsd", math.nan))
    ion = float(row.get("oc2_ion_network_rmsd", math.nan))
    plddt = float(row.get("oc2_mean_plddt", math.nan))
    if any(math.isnan(x) for x in [barrel, pocket, ion]):
        return 0.0, False, "missing_rmsd"
    barrel_score = 1.0 - min(1.0, barrel / max(1e-8, args.max_barrel_rmsd))
    pocket_score = 1.0 - min(1.0, pocket / max(1e-8, args.max_pocket_rmsd))
    ion_score = 1.0 - min(1.0, ion / max(1e-8, args.max_ion_network_rmsd))
    plddt_score = 0.50 if math.isnan(plddt) else min(1.0, max(0.0, plddt / 100.0))
    score = float(np.clip(0.40 * barrel_score + 0.25 * pocket_score + 0.20 * ion_score + 0.15 * plddt_score, 0.0, 1.0))
    passed = (
        barrel <= args.max_barrel_rmsd
        and pocket <= args.max_pocket_rmsd
        and ion <= args.max_ion_network_rmsd
        and (math.isnan(plddt) or plddt >= args.min_mean_plddt)
        and score >= args.min_softcheck_score
    )
    return score, bool(passed), "ok" if passed else "softcheck_gate_failed"


def run_external(command_template: str, fasta: Path, out_dir: Path, log_path: Path, timeout_sec: int) -> None:
    if not command_template:
        return
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd_text = command_template.format(fasta=shlex.quote(str(fasta)), out_dir=shlex.quote(str(out_dir)))
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a", encoding="utf-8", buffering=1) as log:
        log.write("$ " + cmd_text + "\n")
        # Run through bash so OpenComplex2 users can pass the native command line without adapter-specific parsing.
        proc = subprocess.Popen(["bash", "-lc", cmd_text], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, start_new_session=True)
        assert proc.stdout is not None
        start = time.time()
        for line in proc.stdout:
            print(line, end="", flush=True)
            log.write(line)
            if timeout_sec > 0 and time.time() - start > timeout_sec and proc.poll() is None:
                log.write(f"[terminate] softcheck wall-time exceeded {timeout_sec}s\n")
                try:
                    os.killpg(proc.pid, signal.SIGTERM)
                except ProcessLookupError:
                    pass
                time.sleep(5)
                if proc.poll() is None:
                    try:
                        os.killpg(proc.pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
                break
        code = proc.wait()
        log.write(f"[exit_code] {code}\n")
        if code != 0:
            raise RuntimeError(f"External softcheck command failed with exit code {code}. See {log_path}")


def import_softcheck_csv(path: Path) -> pd.DataFrame:
    df = load_table(path)
    if "sequence" in df.columns and "sequence_sha1" not in df.columns:
        df["sequence_sha1"] = df["sequence"].astype(str).str.strip().str.upper().map(sequence_sha1)
    if "sequence_sha1" not in df.columns:
        raise ValueError("Softcheck import CSV must contain sequence_sha1 or sequence.")
    # Normalize common column aliases while preserving backend-specific metadata.
    rename = {
        "barrel_rmsd": "oc2_barrel_rmsd",
        "pocket_rmsd": "oc2_pocket_rmsd",
        "ion_network_rmsd": "oc2_ion_network_rmsd",
        "mean_plddt": "oc2_mean_plddt",
        "pdb_path": "oc2_pdb_path",
        "structure_status": "oc2_status",
        "structure_pass": "oc2_pass",
        "softcheck_score": "oc2_softcheck_score",
    }
    df = df.rename(columns={k: v for k, v in rename.items() if k in df.columns and v not in df.columns})
    return df


def main() -> None:
    p = argparse.ArgumentParser(description="Run or import V7 OpenComplex2-style soft validation metrics.")
    p.add_argument("--candidates-csv", required=True)
    p.add_argument("--reference-pdb", required=True)
    p.add_argument("--out-csv", required=True)
    p.add_argument("--work-dir", required=True)
    p.add_argument("--runner", choices=["import", "metrics-only", "external"], default="metrics-only")
    p.add_argument("--import-csv", default="")
    p.add_argument("--pdb-dir", default="", help="Existing OpenComplex2 PDB directory for metrics-only mode.")
    p.add_argument("--external-command-template", default="", help="External command with {fasta} and {out_dir} placeholders.")
    p.add_argument("--top-n", type=int, default=80)
    p.add_argument("--timeout-sec", type=int, default=21600)
    p.add_argument("--max-barrel-rmsd", type=float, default=2.50)
    p.add_argument("--max-pocket-rmsd", type=float, default=1.20)
    p.add_argument("--max-ion-network-rmsd", type=float, default=1.20)
    p.add_argument("--min-mean-plddt", type=float, default=70.0)
    p.add_argument("--min-softcheck-score", type=float, default=0.45)
    args = p.parse_args()

    work_dir = Path(args.work_dir)
    out_csv = Path(args.out_csv)
    work_dir.mkdir(parents=True, exist_ok=True)
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    if args.runner == "import":
        if not args.import_csv:
            raise ValueError("--import-csv is required when --runner import.")
        imported = import_softcheck_csv(Path(args.import_csv))
        imported.to_csv(out_csv, index=False)
        print(json.dumps({"mode": "import", "rows": int(len(imported)), "out_csv": str(out_csv)}, indent=2))
        return

    candidates = ensure_candidate_keys(load_table(args.candidates_csv)).head(int(args.top_n)).copy()
    fasta = work_dir / "softcheck_tasks.fasta"
    write_fasta(candidates, fasta)
    oc2_out_dir = work_dir / "opencomplex2_outputs"
    if args.runner == "external":
        run_external(args.external_command_template, fasta, oc2_out_dir, work_dir / "logs" / "opencomplex2_softcheck.log", args.timeout_sec)
        pdb_dir = oc2_out_dir
    else:
        pdb_dir = Path(args.pdb_dir) if args.pdb_dir else oc2_out_dir

    ref_atoms = parse_ca_atoms(Path(args.reference_pdb))
    pdb_index = discover_pdbs(pdb_dir)
    rows: list[dict[str, Any]] = []
    for _, row in candidates.iterrows():
        pdb = match_pdb(row, pdb_index)
        payload: dict[str, Any] = {
            "candidate_id": row.get("candidate_id"),
            "Seq_ID": row.get("Seq_ID"),
            "sequence": row.get("sequence"),
            "sequence_sha1": row.get("sequence_sha1"),
            "oc2_status": "missing_pdb",
            "oc2_pdb_path": "",
            "oc2_barrel_rmsd": math.nan,
            "oc2_pocket_rmsd": math.nan,
            "oc2_ion_network_rmsd": math.nan,
            "oc2_mean_plddt": math.nan,
        }
        if pdb is not None:
            atoms = parse_ca_atoms(pdb)
            payload.update({
                "oc2_status": "rmsd_checked",
                "oc2_pdb_path": str(pdb),
                "oc2_barrel_rmsd": subset_rmsd(ref_atoms, atoms, BARREL_POSITIONS),
                "oc2_pocket_rmsd": subset_rmsd(ref_atoms, atoms, POCKET_POSITIONS),
                "oc2_ion_network_rmsd": subset_rmsd(ref_atoms, atoms, ION_NETWORK_POSITIONS),
                "oc2_mean_plddt": mean_plddt(atoms),
            })
        score, passed, reason = score_from_metrics(payload, args)
        payload["oc2_softcheck_score"] = score
        payload["oc2_pass"] = passed
        payload["oc2_reason"] = reason
        rows.append(payload)
    result = pd.DataFrame(rows)
    result.to_csv(out_csv, index=False)
    summary = {
        "runner": args.runner,
        "candidates": int(len(candidates)),
        "rmsd_checked": int(result["oc2_status"].astype(str).eq("rmsd_checked").sum()),
        "oc2_pass": int(result["oc2_pass"].astype(bool).sum()),
        "out_csv": str(out_csv),
    }
    (work_dir / "softcheck_summary_v7.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
