#!/usr/bin/env python3
"""Traceable structure backend adapter for SynBio GFP V6.

The adapter is intentionally backend-agnostic: any runner that produces PDB files can be
used, and the final contract is always the same structure_metrics.csv consumed by V5.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import queue
import re
import shlex
import shutil
import signal
import subprocess
import sys
import threading
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, Optional, Sequence, Tuple

import numpy as np

DEFAULT_BARREL_RANGES = "1-239"
DEFAULT_POCKET_POSITIONS = "65,66,67,68,69,70,95,96,145,146,147,148,203,204"
DEFAULT_ION_NETWORK_POSITIONS = "39,40,41,42,43,44,96,97,98,99,164,165,166,167"
PDB_SUFFIXES = (".pdb", ".ent")

def sequence_sha1(seq: str) -> str:
    return hashlib.sha1(str(seq).strip().upper().encode("utf-8")).hexdigest()


@dataclass(frozen=True)
class FastaRecord:
    seq_id: str
    sequence: str
    description: str


def safe_name(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("_")


def normalize_id(value: str) -> str:
    x = value.strip()
    x = x.split("|", 1)[0]
    x = re.sub(r"_relaxed_rank_.*$", "", x)
    x = re.sub(r"_unrelaxed_rank_.*$", "", x)
    x = re.sub(r"_rank_.*$", "", x)
    x = re.sub(r"_model_.*$", "", x)
    x = re.sub(r"_pred_.*$", "", x)
    return safe_name(x)


def parse_fasta(path: Path) -> list[FastaRecord]:
    records: list[FastaRecord] = []
    header: str | None = None
    seq_chunks: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append(make_record(header, seq_chunks))
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line)
    if header is not None:
        records.append(make_record(header, seq_chunks))
    if not records:
        raise ValueError(f"No FASTA records found: {path}")
    return records


def make_record(header: str, seq_chunks: Sequence[str]) -> FastaRecord:
    # Keep candidate_id stable by taking the token before metadata such as "|mutations=".
    seq_id = normalize_id(header.split()[0])
    seq = "".join(seq_chunks).replace(" ", "").upper()
    if not re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+", seq):
        raise ValueError(f"Invalid amino-acid FASTA sequence for {seq_id}")
    return FastaRecord(seq_id=seq_id, sequence=seq, description=header)


def write_fasta(records: Sequence[FastaRecord], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for rec in records:
            handle.write(f">{rec.seq_id}\n")
            for i in range(0, len(rec.sequence), 80):
                handle.write(rec.sequence[i:i + 80] + "\n")


def discover_pdbs(pdb_dir: Path) -> dict[str, Path]:
    out: dict[str, Path] = {}
    if not pdb_dir.exists():
        return out
    for path in sorted(pdb_dir.rglob("*")):
        if not path.is_file() or path.suffix.lower() not in PDB_SUFFIXES:
            continue
        candidates = [path.stem] + [p.name for p in path.parents if p != pdb_dir.parent]
        for item in candidates:
            key = normalize_id(item)
            if key and key not in out:
                out[key] = path
    return out


def match_pdb(seq_id: str, pdb_index: dict[str, Path]) -> Path | None:
    sid = normalize_id(seq_id)
    if sid in pdb_index:
        return pdb_index[sid]
    for key, path in pdb_index.items():
        if key.startswith(sid) or sid.startswith(key):
            return path
    return None


def terminate_process_group(proc: subprocess.Popen[Any], log: Any, reason: str) -> None:
    # Kill the whole backend process group so nested ColabFold/MMseqs workers do not survive a timeout.
    if proc.poll() is not None:
        return
    log.write(f"[terminate] {reason}\n")
    log.flush()
    try:
        os.killpg(proc.pid, signal.SIGTERM)
    except ProcessLookupError:
        return
    except Exception as exc:
        log.write(f"[terminate_warning] SIGTERM failed: {exc}\n")
    deadline = time.time() + 15.0
    while proc.poll() is None and time.time() < deadline:
        time.sleep(0.5)
    if proc.poll() is None:
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        except Exception as exc:
            log.write(f"[terminate_warning] SIGKILL failed: {exc}\n")


def run_cmd(cmd: Sequence[str], log_path: Path, dry_run: bool = False, timeout_sec: int = 0) -> int:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    command_line = "\n$ " + " ".join(shlex.quote(str(x)) for x in cmd)
    print(command_line, flush=True)
    with log_path.open("a", encoding="utf-8", buffering=1) as log:
        log.write(command_line + "\n")
        if timeout_sec and timeout_sec > 0:
            log.write(f"[timeout_sec] {timeout_sec}\n")
        if dry_run:
            log.write("[dry_run] skipped\n")
            print("[dry_run] skipped", flush=True)
            return 0

        # Start a new process group so timeout handling can clean up ColabFold child processes safely.
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            universal_newlines=True,
            start_new_session=True,
        )
        assert proc.stdout is not None
        output_queue: queue.Queue[str | None] = queue.Queue()

        def _reader() -> None:
            try:
                for child_line in proc.stdout:  # type: ignore[union-attr]
                    output_queue.put(child_line)
            finally:
                output_queue.put(None)

        reader_thread = threading.Thread(target=_reader, name="backend-output-reader", daemon=True)
        reader_thread.start()
        start = time.time()
        timed_out = False
        reader_done = False
        while not reader_done or proc.poll() is None or not output_queue.empty():
            try:
                line = output_queue.get(timeout=0.5)
            except queue.Empty:
                line = None
            if line is None:
                reader_done = True if proc.poll() is not None else reader_done
            else:
                log.write(line)
                log.flush()
                print(line, end="", flush=True)
            if timeout_sec and timeout_sec > 0 and proc.poll() is None and time.time() - start > timeout_sec:
                timed_out = True
                terminate_process_group(proc, log, f"wall-time exceeded {timeout_sec}s")
                break
        return_code = 124 if timed_out else proc.wait()
        reader_thread.join(timeout=2.0)
        exit_line = f"[exit_code] {return_code}"
        log.write(exit_line + "\n")
        print(exit_line, flush=True)
        return return_code


def split_positions(spec: str) -> list[int]:
    vals: list[int] = []
    for item in spec.split(","):
        item = item.strip()
        if not item:
            continue
        if "-" in item:
            left, right = item.split("-", 1)
            vals.extend(range(int(left), int(right) + 1))
        else:
            vals.append(int(item))
    return sorted(set(vals))


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


def load_priority(priority_csv: Path | None) -> dict[str, dict[str, Any]]:
    if not priority_csv or not priority_csv.exists():
        return {}
    rows: dict[str, dict[str, Any]] = {}
    with priority_csv.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            for key in ("candidate_id", "Seq_ID", "seq_id", "id"):
                if row.get(key):
                    rows[normalize_id(row[key])] = row
    return rows


def resolve_colabfold_batch(explicit_runner: str | None) -> str | None:
    """Resolve an explicit ColabFold command/path before falling back to PATH discovery."""
    requested = (explicit_runner or "").strip()
    if requested:
        candidate = Path(requested).expanduser()
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return str(candidate.resolve())
        command = shutil.which(requested)
        if command:
            return command
        raise FileNotFoundError(
            f"--colabfold-batch is not executable or discoverable: {requested}"
        )
    return shutil.which("colabfold_batch")


def dispatch_colabfold(records: Sequence[FastaRecord], args: argparse.Namespace, work_dir: Path) -> None:
    runner = resolve_colabfold_batch(args.colabfold_batch)
    if not runner:
        raise FileNotFoundError("colabfold_batch was not found in PATH or --colabfold-batch.")
    fasta = work_dir / "colabfold_missing.fasta"
    write_fasta(records, fasta)
    cmd = [runner]
    if args.colabfold_extra_args:
        cmd.extend(shlex.split(args.colabfold_extra_args))
    cmd.extend([str(fasta), str(args.af_output_dir)])
    code = run_cmd(cmd, work_dir / "logs" / "colabfold_batch.log", args.dry_run, timeout_sec=args.backend_timeout_sec)
    if code != 0:
        raise RuntimeError(f"colabfold_batch failed with exit code {code}")


def dispatch_external(records: Sequence[FastaRecord], args: argparse.Namespace, work_dir: Path) -> None:
    if not args.external_command_template:
        raise ValueError("--external-command-template is required for runner=external")
    fasta = work_dir / "external_missing.fasta"
    write_fasta(records, fasta)
    rendered = args.external_command_template.format(input_fasta=str(fasta), output_dir=str(args.af_output_dir), work_dir=str(work_dir))
    cmd = ["bash", "-lc", rendered]
    code = run_cmd(cmd, work_dir / "logs" / "external_backend.log", args.dry_run, timeout_sec=args.backend_timeout_sec)
    if code != 0:
        raise RuntimeError(f"external backend failed with exit code {code}")


def build_metrics(records: Sequence[FastaRecord], args: argparse.Namespace) -> None:
    output_dir = Path(args.af_output_dir)
    metrics_csv = Path(args.metrics_csv)
    metrics_csv.parent.mkdir(parents=True, exist_ok=True)
    pdb_index = discover_pdbs(output_dir)
    priority = load_priority(Path(args.priority_csv) if args.priority_csv else None)

    ref_atoms = None
    if args.reference_pdb:
        ref_atoms = parse_ca_atoms(Path(args.reference_pdb), args.reference_chain)
        if len(ref_atoms) < 3:
            raise ValueError(f"Reference PDB has too few CA atoms: {args.reference_pdb}")

    barrel_positions = split_positions(args.barrel_positions)
    pocket_positions = split_positions(args.pocket_positions)
    ion_positions = split_positions(args.ion_network_positions)

    fields = [
        "candidate_id", "Seq_ID", "sequence", "sequence_sha1", "mutations", "v5_objective",
        "barrel_rmsd", "pocket_rmsd", "ion_network_rmsd", "mean_plddt",
        "pdb_path", "has_real_structure", "structure_status", "structure_pass",
    ]
    with metrics_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for rec in records:
            row0 = priority.get(normalize_id(rec.seq_id), {})
            pdb_path = match_pdb(rec.seq_id, pdb_index)
            barrel = pocket = ion = plddt = math.nan
            status = "missing_pdb"
            if pdb_path:
                pred_atoms = parse_ca_atoms(pdb_path, args.predicted_chain)
                plddt = mean_plddt(pred_atoms)
                if ref_atoms is not None:
                    # RMSD metrics are computed on CA atoms for biologically relevant GFP residue subsets.
                    barrel = subset_rmsd(ref_atoms, pred_atoms, barrel_positions)
                    pocket = subset_rmsd(ref_atoms, pred_atoms, pocket_positions)
                    ion = subset_rmsd(ref_atoms, pred_atoms, ion_positions)
                    status = "rmsd_checked"
                elif args.allow_confidence_proxy:
                    # Confidence proxy is kept explicit so downstream lineage can reject it for final decisions.
                    barrel = max(0.0, (90.0 - plddt) / args.proxy_barrel_scale) if not math.isnan(plddt) else math.nan
                    pocket = max(0.0, (90.0 - plddt) / args.proxy_pocket_scale) if not math.isnan(plddt) else math.nan
                    ion = max(0.0, (90.0 - plddt) / args.proxy_ion_scale) if not math.isnan(plddt) else math.nan
                    status = "confidence_proxy_checked"
                else:
                    # PDBs are indexed for traceability, but final gate remains closed without reference-based metrics.
                    status = "pdb_without_reference"
            passed = all(not math.isnan(v) for v in [barrel, pocket, ion]) and barrel <= args.max_barrel_rmsd and pocket <= args.max_pocket_rmsd and ion <= args.max_ion_network_rmsd
            writer.writerow({
                "candidate_id": row0.get("candidate_id", rec.seq_id),
                "Seq_ID": row0.get("Seq_ID", rec.seq_id),
                "sequence": rec.sequence,
                "sequence_sha1": sequence_sha1(rec.sequence),
                "mutations": row0.get("mutations", ""),
                "v5_objective": row0.get("v5_objective", ""),
                "barrel_rmsd": "" if math.isnan(barrel) else f"{barrel:.6f}",
                "pocket_rmsd": "" if math.isnan(pocket) else f"{pocket:.6f}",
                "ion_network_rmsd": "" if math.isnan(ion) else f"{ion:.6f}",
                "mean_plddt": "" if math.isnan(plddt) else f"{plddt:.6f}",
                "pdb_path": str(pdb_path) if pdb_path else "",
                "has_real_structure": str(bool(pdb_path and status == "rmsd_checked")),
                "structure_status": status,
                "structure_pass": str(bool(passed)),
            })
    print(f"Wrote structure metrics: {metrics_csv}")


def main() -> None:
    p = argparse.ArgumentParser(description="Run a V5 structure backend and emit structure_metrics.csv.")
    p.add_argument("--fasta-paths", "--fasta_paths", dest="fasta_paths", required=True)
    p.add_argument("--af-output-dir", required=True)
    p.add_argument("--metrics-csv", required=True)
    p.add_argument("--priority-csv", default=None)
    p.add_argument("--work-dir", default="outputs/v5_structure_adapter")
    p.add_argument("--runner", choices=["auto", "colabfold", "external", "metrics-only"], default="auto")
    p.add_argument("--skip-existing", action="store_true", default=True)
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--colabfold-batch", default=None)
    p.add_argument("--colabfold-extra-args", default="--num-models 1 --num-recycle 1 --msa-mode mmseqs2_uniref_env")
    p.add_argument("--backend-timeout-sec", type=int, default=int(os.environ.get("V62_BACKEND_TIMEOUT_SEC", "7200")))
    p.add_argument("--external-command-template", default=None, help="Shell template with {input_fasta}, {output_dir}, and {work_dir} placeholders.")
    p.add_argument("--reference-pdb", default=None)
    p.add_argument("--reference-chain", default=None)
    p.add_argument("--predicted-chain", default=None)
    p.add_argument("--barrel-positions", default=DEFAULT_BARREL_RANGES)
    p.add_argument("--pocket-positions", default=DEFAULT_POCKET_POSITIONS)
    p.add_argument("--ion-network-positions", default=DEFAULT_ION_NETWORK_POSITIONS)
    p.add_argument("--max-barrel-rmsd", type=float, default=1.50)
    p.add_argument("--max-pocket-rmsd", type=float, default=0.85)
    p.add_argument("--max-ion-network-rmsd", type=float, default=0.75)
    p.add_argument("--allow-confidence-proxy", action="store_true")
    p.add_argument("--proxy-barrel-scale", type=float, default=30.0)
    p.add_argument("--proxy-pocket-scale", type=float, default=40.0)
    p.add_argument("--proxy-ion-scale", type=float, default=40.0)
    args = p.parse_args()

    records: list[FastaRecord] = []
    for item in args.fasta_paths.split(","):
        records.extend(parse_fasta(Path(item.strip())))
    work_dir = Path(args.work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    Path(args.af_output_dir).mkdir(parents=True, exist_ok=True)

    pdb_index = discover_pdbs(Path(args.af_output_dir))
    missing = [rec for rec in records if match_pdb(rec.seq_id, pdb_index) is None]
    (work_dir / "v5_structure_records.json").write_text(json.dumps({"total": len(records), "missing": len(missing), "runner": args.runner}, indent=2), encoding="utf-8")
    print(f"Structure records: total={len(records)} missing={len(missing)}")

    runner = args.runner
    if runner == "auto":
        # Auto mode must honor --colabfold-batch; silent metrics-only fallback invalidates V7 structure evidence.
        if missing:
            colabfold_runner = resolve_colabfold_batch(args.colabfold_batch)
            if colabfold_runner:
                runner = "colabfold"
            elif args.external_command_template:
                runner = "external"
            else:
                raise RuntimeError(
                    "runner=auto found missing structures but no executable ColabFold runner or external backend. "
                    "Set --colabfold-batch, add colabfold_batch to PATH, set --external-command-template, "
                    "or select --runner metrics-only explicitly."
                )
        else:
            runner = "metrics-only"
    if missing and runner == "colabfold":
        dispatch_colabfold(missing, args, work_dir)
    elif missing and runner == "external":
        dispatch_external(missing, args, work_dir)
    elif missing and runner == "metrics-only":
        print("metrics-only selected; no new predictions will be generated.")

    build_metrics(records, args)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
