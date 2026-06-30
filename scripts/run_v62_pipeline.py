#!/usr/bin/env python3
"""Traceable SynBio GFP V6.1 Adaptive runner.

Each stage can be run independently, while this orchestrator provides a one-click path.
Structure handling is delegated to the V6.1 resolver so reusable PDBs are consumed before
new ColabFold jobs are launched.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import shlex
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]


def now_id() -> str:
    return time.strftime("v62_%Y%m%d_%H%M%S")


def shell_quote(cmd: list[str]) -> str:
    return " ".join(shlex.quote(str(x)) for x in cmd)


def run_stage(name: str, cmd: list[str], log_path: Path, status_rows: list[dict[str, Any]], dry_run: bool = False) -> None:
    start = time.time()
    log_path.parent.mkdir(parents=True, exist_ok=True)
    row = {"stage": name, "status": "running", "seconds": "", "log_path": str(log_path), "command": shell_quote(cmd)}
    status_rows.append(row)
    print(f"[stage] {name}\n$ {shell_quote(cmd)}", flush=True)
    with log_path.open("a", encoding="utf-8", buffering=1) as log:
        log.write(f"[stage] {name}\n$ {shell_quote(cmd)}\n")
        if dry_run:
            row.update({"status": "dry_run", "seconds": f"{time.time() - start:.3f}"})
            log.write("[dry_run] skipped\n")
            return
        # Stream child output to the parent log in real time for SSH-safe monitoring.
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)
        assert proc.stdout is not None
        for line in proc.stdout:
            print(line, end="", flush=True)
            log.write(line)
        code = proc.wait()
        row.update({"status": "pass" if code == 0 else "fail", "seconds": f"{time.time() - start:.3f}"})
        log.write(f"[exit_code] {code}\n")
        if code != 0:
            raise RuntimeError(f"Stage failed: {name}; see {log_path}")


def write_status(status_rows: list[dict[str, Any]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = ["stage", "status", "seconds", "log_path", "command"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(status_rows)


def main() -> None:
    p = argparse.ArgumentParser(description="Run SynBio GFP V6.1 Adaptive pipeline.")
    p.add_argument("--data-dir", default=os.environ.get("DATA_DIR", "/input0/2026Protein Design"))
    p.add_argument("--team-name", default=os.environ.get("TEAM_NAME", "YourTeamName"))
    p.add_argument("--run-id", default=os.environ.get("RUN_ID", now_id()))
    p.add_argument("--feature-mode", choices=["simple", "esm"], default=os.environ.get("FEATURE_MODE", "simple"))
    p.add_argument("--esm-model", default=os.environ.get("ESM_MODEL", "esm2_t30_150M_UR50D"))
    p.add_argument("--max-train-samples", type=int, default=int(os.environ.get("MAX_TRAIN_SAMPLES", "20000")))
    p.add_argument("--n-candidates", type=int, default=int(os.environ.get("N_CANDIDATES", "20000")))
    p.add_argument("--config", default=str(ROOT / "configs" / "v62_adaptive.yaml"))
    p.add_argument("--cache-dir", default=os.environ.get("V62_CACHE_DIR", str(ROOT / "cache" / "v62")))
    p.add_argument("--output-root", default=os.environ.get("V62_OUTPUT_ROOT", str(ROOT / "outputs")))
    p.add_argument("--gpu-profile", default=os.environ.get("GPU_PROFILE", "auto"))
    p.add_argument("--structure-policy", choices=["reuse-only", "reuse-first", "predict-missing", "metrics-only"], default=os.environ.get("STRUCTURE_POLICY", "reuse-first"))
    p.add_argument("--structure-runner", choices=["auto", "colabfold", "external", "metrics-only", "skip"], default=os.environ.get("STRUCTURE_RUNNER", "auto"))
    p.add_argument("--reference-pdb", default=os.environ.get("V62_REFERENCE_PDB", os.environ.get("REF_PDB", str(ROOT / "reference" / "reference_gfp.pdb"))))
    p.add_argument("--pdb-search-roots", default=os.environ.get("PDB_SEARCH_ROOTS", ":".join([
        "/hyperai/home/synbio_gfp_v51_complete/outputs",
        "/hyperai/home/synbio_gfp_v6_complete/outputs",
        "/hyperai/home/synbio_gfp_v5_complete/outputs",
        "/hyperai/home/synbio_gfp_v4_complete/outputs",
    ])))
    p.add_argument("--structure-metrics-csv", default=os.environ.get("V62_STRUCTURE_METRICS_CSV", ""))
    p.add_argument("--proteinmpnn-score-csv", default=os.environ.get("V62_PROTEINMPNN_SCORE_CSV", ""))
    p.add_argument("--colabfold-extra-args", default=os.environ.get("COLABFOLD_EXTRA_ARGS", "--num-models 1 --num-recycle 1"))
    p.add_argument("--colabfold-batch", default=os.environ.get("COLABFOLD_BATCH", ""))
    p.add_argument("--colabfold-msa-mode", default=os.environ.get("V62_COLABFOLD_MSA_MODE", "mmseqs2_uniref_env"))
    p.add_argument("--structure-chunk-size", type=int, default=int(os.environ.get("STRUCTURE_CHUNK_SIZE", "12")))
    p.add_argument("--chunk-timeout-sec", type=int, default=int(os.environ.get("V62_COLABFOLD_CHUNK_TIMEOUT_SEC", "7200")))
    p.add_argument("--chunk-retries", type=int, default=int(os.environ.get("V62_COLABFOLD_CHUNK_RETRIES", "1")))
    p.add_argument("--external-structure-command", default=os.environ.get("EXTERNAL_STRUCTURE_COMMAND", ""))
    p.add_argument("--run-proteinmpnn", action="store_true", default=os.environ.get("RUN_PROTEINMPNN", "0") == "1")
    p.add_argument("--proteinmpnn-dir", default=os.environ.get("PROTEINMPNN_DIR", "/hyperai/home/tools/ProteinMPNN"))
    p.add_argument("--previous-ranked-csv", default=os.environ.get("V62_PREVIOUS_RANKED_CSV", os.environ.get("V6_PREVIOUS_RANKED_CSV", "")))
    p.add_argument("--feedback-csv", default=os.environ.get("V62_FEEDBACK_CSV", os.environ.get("V6_FEEDBACK_CSV", "")))
    p.add_argument("--brightness-teacher-csv", default=os.environ.get("V62_BRIGHTNESS_TEACHER_CSV", ""))
    p.add_argument("--allow-proxy-final", action="store_true", default=os.environ.get("ALLOW_PROXY_FINAL", "0") == "1")
    p.add_argument("--dry-run", action="store_true")
    args = p.parse_args()

    run_dir = Path(args.output_root) / args.run_id
    dirs = {
        "preflight": run_dir / "00_preflight",
        "stage1": run_dir / "01_stage1",
        "structure": run_dir / "02_structure",
        "proteinmpnn": run_dir / "03_proteinmpnn",
        "final": run_dir / "04_final",
        "lineage": run_dir / "lineage",
        "logs": run_dir / "logs",
    }
    for d in dirs.values():
        d.mkdir(parents=True, exist_ok=True)

    context = vars(args).copy()
    context.update({"v62_dir": str(ROOT), "run_dir": str(run_dir), "python": sys.executable})
    (run_dir / "v62_run_context.json").write_text(json.dumps(context, indent=2, ensure_ascii=False), encoding="utf-8")
    status_rows: list[dict[str, Any]] = []
    status_csv = dirs["lineage"] / "v62_stage_status.csv"
    commands_sh = dirs["lineage"] / "v62_commands.sh"
    command_lines: list[str] = ["#!/usr/bin/env bash", "set -euo pipefail", f"cd {shlex.quote(str(ROOT))}"]

    try:
        preflight_cmd = [sys.executable, str(ROOT / "deploy" / "00_preflight_v6.py"), "--v5-dir", str(ROOT), "--data-dir", args.data_dir, "--out-dir", str(dirs["preflight"])]
        command_lines.append(shell_quote(preflight_cmd))
        run_stage("00_preflight", preflight_cmd, dirs["logs"] / "00_preflight.log", status_rows, args.dry_run)
        write_status(status_rows, status_csv)

        stage1_cmd = [
            sys.executable, str(ROOT / "scripts" / "run_v62_complete.py"),
            "--config", args.config,
            "--data-dir", args.data_dir,
            "--team-name", args.team_name,
            "--feature-mode", args.feature_mode,
            "--esm-model", args.esm_model,
            "--max-train-samples", str(args.max_train_samples),
            "--n-candidates", str(args.n_candidates),
            "--out-dir", str(dirs["stage1"]),
            "--cache-dir", args.cache_dir,
        ]
        if args.previous_ranked_csv:
            stage1_cmd += ["--previous-ranked-csv", args.previous_ranked_csv]
        if args.feedback_csv:
            stage1_cmd += ["--feedback-csv", args.feedback_csv]
        if args.brightness_teacher_csv:
            stage1_cmd += ["--brightness-teacher-csv", args.brightness_teacher_csv]
        command_lines.append(shell_quote(stage1_cmd))
        run_stage("01_candidate_generation", stage1_cmd, dirs["logs"] / "01_stage1.log", status_rows, args.dry_run)
        write_status(status_rows, status_csv)

        priority_csv = dirs["stage1"] / "structure_priority_top200_v6.csv"
        structure_metrics = dirs["structure"] / "structure_metrics.csv"
        pdb_dir = dirs["structure"] / "predictions"
        pdb_dir.mkdir(parents=True, exist_ok=True)

        if args.structure_metrics_csv:
            shutil.copy2(args.structure_metrics_csv, structure_metrics)
            run_stage("02_structure_metrics_import", ["bash", "-lc", f"wc -l {shlex.quote(str(structure_metrics))}"], dirs["logs"] / "02_structure.log", status_rows, args.dry_run)
        else:
            resolver_cmd = [
                sys.executable, str(ROOT / "tools" / "resolve_structures_v62.py"),
                "--priority-csv", str(priority_csv),
                "--reference-pdb", args.reference_pdb,
                "--prediction-dir", str(pdb_dir),
                "--metrics-csv", str(structure_metrics),
                "--work-dir", str(dirs["structure"]),
                "--pdb-search-roots", args.pdb_search_roots,
                "--structure-policy", args.structure_policy,
                "--runner", "metrics-only" if args.structure_runner == "skip" else args.structure_runner,
                "--colabfold-extra-args", args.colabfold_extra_args,
                "--colabfold-msa-mode", args.colabfold_msa_mode,
                "--chunk-size", str(args.structure_chunk_size),
                "--chunk-timeout-sec", str(args.chunk_timeout_sec),
                "--chunk-retries", str(args.chunk_retries),
            ]
            if args.colabfold_batch:
                resolver_cmd += ["--colabfold-batch", args.colabfold_batch]
            if args.external_structure_command:
                resolver_cmd += ["--external-command-template", args.external_structure_command]
            command_lines.append(shell_quote(resolver_cmd))
            run_stage("02_resolve_structures", resolver_cmd, dirs["logs"] / "02_structure.log", status_rows, args.dry_run)
        write_status(status_rows, status_csv)

        proteinmpnn_scores = dirs["proteinmpnn"] / "proteinmpnn_scores.csv"
        if args.proteinmpnn_score_csv:
            shutil.copy2(args.proteinmpnn_score_csv, proteinmpnn_scores)
            run_stage("03_proteinmpnn_scores_import", ["bash", "-lc", f"wc -l {shlex.quote(str(proteinmpnn_scores))}"], dirs["logs"] / "03_proteinmpnn.log", status_rows, args.dry_run)
        elif args.run_proteinmpnn:
            pm_cmd = [
                sys.executable, str(ROOT / "tools" / "run_proteinmpnn_v6.py"),
                "--pdb-dir", str(pdb_dir),
                "--out-csv", str(proteinmpnn_scores),
                "--work-dir", str(dirs["proteinmpnn"] / "score_work"),
                "--proteinmpnn-dir", args.proteinmpnn_dir,
                "--priority-csv", str(priority_csv),
                "--chain-id", "auto",
            ]
            command_lines.append(shell_quote(pm_cmd))
            run_stage("03_proteinmpnn_score_only", pm_cmd, dirs["logs"] / "03_proteinmpnn.log", status_rows, args.dry_run)
        write_status(status_rows, status_csv)

        frozen_candidates = dirs["stage1"] / "candidate_universe_v62.csv"
        if not frozen_candidates.exists():
            frozen_candidates = dirs["stage1"] / "all_ranked_candidates_v6.csv"
        final_cmd = [
            sys.executable, str(ROOT / "scripts" / "run_v62_final_from_manifest.py"),
            "--config", args.config,
            "--frozen-candidates-csv", str(frozen_candidates),
            "--structure-metrics-csv", str(structure_metrics),
            "--out-dir", str(dirs["final"]),
            "--team-name", args.team_name,
        ]
        if proteinmpnn_scores.exists() or args.dry_run:
            final_cmd += ["--proteinmpnn-score-csv", str(proteinmpnn_scores)]
        command_lines.append(shell_quote(final_cmd))
        run_stage("04_final_scoring_portfolio", final_cmd, dirs["logs"] / "04_final.log", status_rows, args.dry_run)
        write_status(status_rows, status_csv)

        lineage_cmd = [sys.executable, str(ROOT / "scripts" / "v6_lineage.py"), "--run-dir", str(run_dir), "--out-dir", str(dirs["lineage"])]
        command_lines.append(shell_quote(lineage_cmd))
        run_stage("05_lineage_index", lineage_cmd, dirs["logs"] / "05_lineage.log", status_rows, args.dry_run)
        write_status(status_rows, status_csv)
    finally:
        commands_sh.write_text("\n".join(command_lines) + "\n", encoding="utf-8")
        try:
            commands_sh.chmod(0o755)
        except Exception:
            pass
        write_status(status_rows, status_csv)

    print("V6.2 evidence-safe run completed")
    print(f"run_dir: {run_dir}")
    print(f"stage_status: {status_csv}")
    print(f"commands: {commands_sh}")
    print(f"submission: {dirs['final'] / 'submission_v62.csv'}")


if __name__ == "__main__":
    main()
