#!/usr/bin/env python3
"""V6.2 recovery-safe structure resolver.

This resolver implements the artifact contract required by the V6.2 evidence-safe
pipeline: write recovery FASTA/state files before any long-running backend call,
reuse exact-sequence historical PDBs, predict missing structures in chunks, and
always recompute current-run RMSD metrics from normalized PDBs.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import queue
import re
import shlex
import shutil
import signal
import subprocess
import sys
import threading
import time
from pathlib import Path
from typing import Any

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PDB_SUFFIXES = {".pdb", ".ent"}
CANDIDATE_CSV_PATTERNS = [
    "candidate_universe*.csv",
    "all_ranked_candidates*.csv",
    "structure_priority*.csv",
    "ranked_top*.csv",
    "final_top*.csv",
    "selected_candidates*.csv",
]


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
    roots: list[Path] = []
    for item in re.split(r"[:;,]", text or ""):
        item = item.strip()
        if item:
            roots.append(Path(item).expanduser())
    return roots


def load_priority(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    if "sequence" not in df.columns:
        raise ValueError(f"Priority CSV must contain a sequence column: {path}")
    out = df.copy()
    if "candidate_id" not in out.columns:
        out["candidate_id"] = [f"cand_{i:06d}" for i in range(len(out))]
    if "Seq_ID" not in out.columns:
        out["Seq_ID"] = out["candidate_id"]
    out["candidate_id"] = out["candidate_id"].map(normalize_id)
    out["Seq_ID"] = out["Seq_ID"].astype(str)
    out["sequence"] = out["sequence"].astype(str).str.strip().str.upper()
    out["sequence_sha1"] = out["sequence"].map(sequence_sha1)
    out["sequence_length"] = out["sequence"].str.len()
    return out


def collect_sequence_ids(search_roots: list[Path]) -> dict[str, dict[str, str]]:
    id_map: dict[str, dict[str, str]] = {}
    for root in search_roots:
        if not root.exists():
            continue
        for pat in CANDIDATE_CSV_PATTERNS:
            for csv_path in root.rglob(pat):
                try:
                    df = pd.read_csv(csv_path, usecols=lambda c: c in {"candidate_id", "Seq_ID", "seq_id", "id", "sequence", "sequence_sha1"})
                except Exception:
                    continue
                if "sequence" not in df.columns:
                    continue
                for _, row in df.iterrows():
                    seq = str(row.get("sequence", "")).strip().upper()
                    if not re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]{50,400}", seq):
                        continue
                    payload = {
                        "sequence": seq,
                        "sequence_sha1": str(row.get("sequence_sha1") or sequence_sha1(seq)),
                        "sequence_length": str(len(seq)),
                        "source_csv": str(csv_path),
                    }
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
    src_resolved = src.resolve()
    if dst.exists() or dst.is_symlink():
        try:
            if dst.resolve() == src_resolved:
                # The normalized target already points to the desired PDB; keep it intact on reruns.
                return
        except RuntimeError:
            # Broken or recursive links are unsafe, so they are replaced below.
            pass
        dst.unlink()
    if copy_pdbs:
        shutil.copy2(src_resolved, dst)
    else:
        dst.symlink_to(src_resolved)


def write_fasta(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for _, row in df.iterrows():
            cid = normalize_id(row.get("candidate_id", row.get("Seq_ID", "candidate")))
            seq = str(row["sequence"]).strip().upper()
            handle.write(f">{cid}|sequence_sha1={row['sequence_sha1']}\n")
            for i in range(0, len(seq), 80):
                handle.write(seq[i:i + 80] + "\n")


def write_state(path: Path, **payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload.setdefault("updated_at", time.strftime("%Y-%m-%dT%H:%M:%S"))
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")


def _terminate_process_group(proc: subprocess.Popen[Any], log: Any, reason: str) -> None:
    # Terminate the entire stage process group to avoid orphaned ColabFold/MMseqs workers after timeouts.
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


def run_cmd(cmd: list[str], log_path: Path, dry_run: bool = False, timeout_sec: int = 0) -> int:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    print("$ " + " ".join(shlex.quote(str(x)) for x in cmd), flush=True)
    with log_path.open("a", encoding="utf-8", buffering=1) as log:
        log.write("$ " + " ".join(shlex.quote(str(x)) for x in cmd) + "\n")
        if timeout_sec and timeout_sec > 0:
            log.write(f"[timeout_sec] {timeout_sec}\n")
        if dry_run:
            log.write("[dry_run] skipped\n")
            return 0
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1, start_new_session=True)
        assert proc.stdout is not None
        output_queue: queue.Queue[str | None] = queue.Queue()

        def _reader() -> None:
            try:
                for child_line in proc.stdout:  # type: ignore[union-attr]
                    output_queue.put(child_line)
            finally:
                output_queue.put(None)

        reader_thread = threading.Thread(target=_reader, name="resolver-output-reader", daemon=True)
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
                print(line, end="", flush=True)
                log.write(line)
                log.flush()
            if timeout_sec and timeout_sec > 0 and proc.poll() is None and time.time() - start > timeout_sec:
                timed_out = True
                _terminate_process_group(proc, log, f"wall-time exceeded {timeout_sec}s")
                break
        code = 124 if timed_out else proc.wait()
        reader_thread.join(timeout=2.0)
        log.write(f"[exit_code] {code}\n")
        return int(code)


def chunk_dataframe(df: pd.DataFrame, chunk_size: int) -> list[pd.DataFrame]:
    if df.empty:
        return []
    size = max(1, int(chunk_size))
    return [df.iloc[i:i + size].copy() for i in range(0, len(df), size)]


def build_colabfold_extra_args(extra_args: str, single_sequence_default: bool, msa_mode: str) -> str:
    args = (extra_args or "").strip()
    if "--msa-mode" in args:
        return args
    if single_sequence_default:
        # Single-sequence mode is only retained as an explicit emergency/debug option, not a final-evidence default.
        return (args + " --msa-mode single_sequence").strip()
    if msa_mode:
        # GFP candidates require MSA/template-quality evidence; default to ColabFold's remote MSA mode.
        return (args + f" --msa-mode {msa_mode}").strip()
    return args


def recompute_metrics(priority_csv: Path, full_fasta: Path, prediction_dir: Path, metrics_csv: Path, work_dir: Path, reference_pdb: str, dry_run: bool) -> int:
    metric_cmd = [
        sys.executable, str(ROOT / "tools" / "run_structure_backend_v6.py"),
        "--runner", "metrics-only",
        "--fasta-paths", str(full_fasta),
        "--af-output-dir", str(prediction_dir),
        "--priority-csv", str(priority_csv),
        "--metrics-csv", str(metrics_csv),
        "--work-dir", str(work_dir / "adapter_recompute_metrics"),
        "--reference-pdb", str(reference_pdb),
    ]
    return run_cmd(metric_cmd, work_dir / "resolve_structures_recompute_metrics.log", dry_run=dry_run)


def main() -> None:
    p = argparse.ArgumentParser(description="Resolve V6.2 structures with exact-sequence reuse and recovery-safe chunks.")
    p.add_argument("--priority-csv", required=True)
    p.add_argument("--reference-pdb", required=True)
    p.add_argument("--prediction-dir", required=True)
    p.add_argument("--metrics-csv", required=True)
    p.add_argument("--work-dir", required=True)
    p.add_argument("--pdb-search-roots", default=os.environ.get("PDB_SEARCH_ROOTS", ""))
    p.add_argument("--structure-policy", choices=["reuse-only", "reuse-first", "predict-missing", "metrics-only"], default=os.environ.get("STRUCTURE_POLICY", "reuse-first"))
    p.add_argument("--runner", choices=["auto", "colabfold", "external", "metrics-only"], default=os.environ.get("STRUCTURE_RUNNER", "auto"))
    p.add_argument("--colabfold-extra-args", default=os.environ.get("COLABFOLD_EXTRA_ARGS", "--num-models 1 --num-recycle 1"))
    p.add_argument("--colabfold-batch", default=os.environ.get("COLABFOLD_BATCH", ""))
    p.add_argument("--colabfold-msa-mode", default=os.environ.get("V62_COLABFOLD_MSA_MODE", "mmseqs2_uniref_env"))
    p.add_argument("--external-command-template", default=os.environ.get("EXTERNAL_STRUCTURE_COMMAND", ""))
    p.add_argument("--chunk-size", type=int, default=int(os.environ.get("STRUCTURE_CHUNK_SIZE", "12")))
    p.add_argument("--chunk-timeout-sec", type=int, default=int(os.environ.get("V62_COLABFOLD_CHUNK_TIMEOUT_SEC", "7200")))
    p.add_argument("--chunk-retries", type=int, default=int(os.environ.get("V62_COLABFOLD_CHUNK_RETRIES", "1")))
    p.add_argument("--retry-backoff-sec", type=int, default=int(os.environ.get("V62_COLABFOLD_RETRY_BACKOFF_SEC", "60")))
    p.add_argument("--single-sequence-default", action="store_true", default=os.environ.get("V62_SINGLE_SEQUENCE_DEFAULT", "0") == "1")
    p.add_argument("--continue-on-chunk-failure", action="store_true", default=os.environ.get("V62_CONTINUE_ON_CHUNK_FAILURE", "1") == "1")
    p.add_argument("--copy-pdbs", action="store_true")
    p.add_argument("--dry-run", action="store_true")
    args = p.parse_args()

    priority_csv = Path(args.priority_csv)
    prediction_dir = Path(args.prediction_dir)
    metrics_csv = Path(args.metrics_csv)
    work_dir = Path(args.work_dir)
    prediction_dir.mkdir(parents=True, exist_ok=True)
    work_dir.mkdir(parents=True, exist_ok=True)
    state_path = work_dir / "resolver_state.json"

    priority = load_priority(priority_csv)
    full_fasta = work_dir / "all_structure_tasks.fasta"
    # Write the full task FASTA before any backend call so interrupted runs can resume in metrics-only mode.
    write_fasta(priority, full_fasta)
    write_state(state_path, stage="initialized", priority_candidates=len(priority), full_fasta=str(full_fasta))

    roots = split_roots(args.pdb_search_roots)
    roots.extend([priority_csv.parents[1] if len(priority_csv.parents) > 1 else priority_csv.parent, prediction_dir])
    roots = [r for i, r in enumerate(roots) if r and r.exists() and r not in roots[:i]]

    id_map = collect_sequence_ids(roots)
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
    normalized_rows: list[dict[str, Any]] = []
    for _, row in priority.iterrows():
        cid = normalize_id(row.get("candidate_id", row.get("Seq_ID", "candidate")))
        seq_hash = str(row["sequence_sha1"])
        src_pdb = registry_by_hash.get(seq_hash)
        normalized = prediction_dir / f"{cid}.pdb"
        if src_pdb and src_pdb.exists():
            link_or_copy(src_pdb, normalized, copy_pdbs=args.copy_pdbs)
            reused_rows.append({**row.to_dict(), "reused_pdb_path": str(normalized), "source_pdb_path": str(src_pdb)})
            normalized_rows.append({"candidate_id": cid, "sequence_sha1": seq_hash, "normalized_pdb_path": str(normalized), "source_pdb_path": str(src_pdb), "status": "reused_or_existing"})
        else:
            missing_rows.append(row.to_dict())

    reused_manifest = work_dir / "reused_pdb_manifest.csv"
    missing_csv = work_dir / "missing_pdb_candidates.csv"
    normalized_manifest = work_dir / "normalized_pdb_manifest.csv"
    pd.DataFrame(reused_rows).to_csv(reused_manifest, index=False)
    missing_df = pd.DataFrame(missing_rows)
    missing_df.to_csv(missing_csv, index=False)
    pd.DataFrame(normalized_rows).to_csv(normalized_manifest, index=False)
    write_state(state_path, stage="reuse_scanned", registry_rows=len(registry), reused=len(reused_rows), missing=len(missing_rows), full_fasta=str(full_fasta))

    extra_args = build_colabfold_extra_args(args.colabfold_extra_args, args.single_sequence_default, args.colabfold_msa_mode)
    chunk_failures: list[dict[str, Any]] = []
    chunk_successes = 0
    if args.structure_policy in {"reuse-first", "predict-missing"} and not missing_df.empty and args.runner != "metrics-only":
        chunks = chunk_dataframe(missing_df, args.chunk_size)
        chunks_dir = work_dir / "chunks"
        # Create the chunk artifact directory before writing CSV/FASTA recovery files.
        chunks_dir.mkdir(parents=True, exist_ok=True)
        for i, chunk in enumerate(chunks, start=1):
            chunk_fasta = chunks_dir / f"missing_chunk_{i:03d}.fasta"
            chunk_csv = chunks_dir / f"missing_chunk_{i:03d}.csv"
            chunk_metrics = chunks_dir / f"missing_chunk_{i:03d}_metrics.csv"
            chunk.to_csv(chunk_csv, index=False)
            write_fasta(chunk, chunk_fasta)
            backend_cmd = [
                sys.executable, str(ROOT / "tools" / "run_structure_backend_v6.py"),
                "--runner", args.runner,
                "--fasta-paths", str(chunk_fasta),
                "--af-output-dir", str(prediction_dir),
                "--priority-csv", str(chunk_csv),
                "--metrics-csv", str(chunk_metrics),
                "--work-dir", str(work_dir / "adapter_predict_chunks" / f"chunk_{i:03d}"),
                "--reference-pdb", str(args.reference_pdb),
                "--colabfold-extra-args", extra_args,
                "--backend-timeout-sec", str(args.chunk_timeout_sec),
            ]
            if args.colabfold_batch:
                backend_cmd += ["--colabfold-batch", args.colabfold_batch]
            if args.external_command_template:
                backend_cmd += ["--external-command-template", args.external_command_template]

            code = 1
            attempts = max(1, int(args.chunk_retries) + 1)
            for attempt in range(1, attempts + 1):
                attempt_log = work_dir / "logs" / f"predict_chunk_{i:03d}_attempt_{attempt:02d}.log"
                write_state(
                    state_path,
                    stage="predicting_chunk",
                    current_chunk=i,
                    current_attempt=attempt,
                    completed_chunks=chunk_successes,
                    failed_chunks=chunk_failures,
                    total_chunks=len(chunks),
                    colabfold_extra_args=extra_args,
                )
                # The resolver timeout is deliberately longer than the backend timeout so the backend can clean up first.
                resolver_timeout = args.chunk_timeout_sec + 90 if args.chunk_timeout_sec > 0 else 0
                code = run_cmd(backend_cmd, attempt_log, dry_run=args.dry_run, timeout_sec=resolver_timeout)
                if code == 0:
                    break
                if attempt < attempts and args.retry_backoff_sec > 0:
                    time.sleep(args.retry_backoff_sec)

            if code == 0:
                chunk_successes += 1
            else:
                failure = {"chunk": i, "rows": int(len(chunk)), "exit_code": code, "fasta": str(chunk_fasta), "attempts": attempts}
                chunk_failures.append(failure)
                if not args.continue_on_chunk_failure:
                    write_state(state_path, stage="chunk_failed", chunk_failures=chunk_failures, completed_chunks=chunk_successes)
                    raise RuntimeError(f"Structure chunk {i} failed with exit code {code}")
            write_state(state_path, stage="predicting_chunks", completed_chunks=chunk_successes, failed_chunks=chunk_failures, total_chunks=len(chunks), colabfold_extra_args=extra_args)
    elif args.structure_policy == "metrics-only" or args.runner == "metrics-only":
        print("metrics-only selected; no new predictions will be generated.")
    elif args.structure_policy == "reuse-only" and not missing_df.empty:
        print(f"reuse-only policy: {len(missing_df)} missing structures will remain missing.")

    recompute_code = recompute_metrics(priority_csv, full_fasta, prediction_dir, metrics_csv, work_dir, args.reference_pdb, args.dry_run)
    if recompute_code != 0:
        raise RuntimeError(f"Metrics recomputation failed with exit code {recompute_code}")

    needs_recovery = pd.DataFrame()
    if metrics_csv.exists():
        metrics = pd.read_csv(metrics_csv)
        status = metrics.get("structure_status", pd.Series("", index=metrics.index)).astype(str)
        passed = metrics.get("structure_pass", pd.Series(False, index=metrics.index)).astype(str).str.lower().eq("true")
        needs_recovery = metrics[status.eq("missing_pdb") | ~passed].copy()
        needs_recovery.to_csv(work_dir / "needs_recovery_candidates.csv", index=False)

    report = {
        "priority_candidates": int(len(priority)),
        "registry_rows": int(len(registry)),
        "reused_pdbs": int(len(reused_rows)),
        "missing_before_prediction": int(len(missing_rows)),
        "chunk_successes": int(chunk_successes),
        "chunk_failures": chunk_failures,
        "needs_recovery_rows": int(len(needs_recovery)),
        "structure_policy": args.structure_policy,
        "colabfold_extra_args": extra_args,
        "chunk_timeout_sec": int(args.chunk_timeout_sec),
        "chunk_retries": int(args.chunk_retries),
        "metrics_csv": str(metrics_csv),
        "prediction_dir": str(prediction_dir),
        "full_fasta": str(full_fasta),
    }
    (work_dir / "structure_reuse_report.md").write_text("# V6.2 Structure Reuse Report\n\n" + "\n".join(f"- {k}: {v}" for k, v in report.items()) + "\n", encoding="utf-8")
    (work_dir / "structure_reuse_report.json").write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    write_state(state_path, stage="done", **report)
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
