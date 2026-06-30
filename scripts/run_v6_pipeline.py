#!/usr/bin/env python3
"""End-to-end traceable SynBio GFP V6 runner.

This orchestrator does not hide the stage commands: it records every command, status,
log path, input path, and output path so each final sequence can be traced backward.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import shlex
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]


def now_id() -> str:
    return time.strftime("v6_%Y%m%d_%H%M%S")


def shell_quote(cmd: list[str]) -> str:
    return " ".join(shlex.quote(str(x)) for x in cmd)


def run_stage(name: str, cmd: list[str], log_path: Path, status_rows: list[dict[str, Any]], dry_run: bool = False) -> None:
    start = time.time()
    log_path.parent.mkdir(parents=True, exist_ok=True)
    row = {"stage": name, "status": "running", "seconds": "", "log_path": str(log_path), "command": shell_quote(cmd)}
    status_rows.append(row)
    stage_header = f"[stage] {name}"
    command_line = "$ " + shell_quote(cmd)
    print(stage_header, flush=True)
    print(command_line, flush=True)
    with log_path.open("a", encoding="utf-8", buffering=1) as log:
        log.write(stage_header + "\n")
        log.write(command_line + "\n")
        if dry_run:
            row.update({"status": "dry_run", "seconds": f"{time.time() - start:.3f}"})
            log.write("[dry_run] skipped\n")
            print("[dry_run] skipped", flush=True)
            return

        # Stream child output line-by-line so the fixed nohup log and per-stage log stay live.
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            universal_newlines=True,
        )
        assert proc.stdout is not None
        for line in proc.stdout:
            log.write(line)
            log.flush()
            print(line, end="", flush=True)
        return_code = proc.wait()
        row.update({"status": "pass" if return_code == 0 else "fail", "seconds": f"{time.time() - start:.3f}"})
        exit_line = f"[exit_code] {return_code}"
        log.write(exit_line + "\n")
        print(exit_line, flush=True)
        if return_code != 0:
            raise RuntimeError(f"Stage failed: {name}; see {log_path}")


def write_status(status_rows: list[dict[str, Any]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = ["stage", "status", "seconds", "log_path", "command"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(status_rows)


def main() -> None:
    p = argparse.ArgumentParser(description="Run traceable SynBio GFP V6 pipeline.")
    p.add_argument("--data-dir", default=os.environ.get("DATA_DIR", "/hyperai/input/input0/2026Protein Design"))
    p.add_argument("--team-name", default=os.environ.get("TEAM_NAME", "YourTeamName"))
    p.add_argument("--run-id", default=os.environ.get("RUN_ID", now_id()))
    p.add_argument("--feature-mode", choices=["simple", "esm"], default=os.environ.get("FEATURE_MODE", "simple"))
    p.add_argument("--esm-model", default=os.environ.get("ESM_MODEL", "esm2_t30_150M_UR50D"))
    p.add_argument("--max-train-samples", type=int, default=int(os.environ.get("MAX_TRAIN_SAMPLES", "20000")))
    p.add_argument("--n-candidates", type=int, default=int(os.environ.get("N_CANDIDATES", "20000")))
    p.add_argument("--config", default=str(ROOT / "configs" / "v6_complete.yaml"))
    p.add_argument("--cache-dir", default=os.environ.get("V6_CACHE_DIR", str(ROOT / "cache" / "v6")))
    p.add_argument("--output-root", default=os.environ.get("V6_OUTPUT_ROOT", str(ROOT / "outputs")))
    p.add_argument("--structure-runner", choices=["auto", "colabfold", "external", "metrics-only", "skip"], default=os.environ.get("STRUCTURE_RUNNER", "metrics-only"))
    p.add_argument("--reference-pdb", default=os.environ.get("V6_REFERENCE_PDB", os.environ.get("REF_PDB", "")))
    p.add_argument("--reuse-pdb-dir", default=os.environ.get("V6_REUSE_PDB_DIR", ""), help="Directory with existing ColabFold/AF PDBs to import instead of rerunning structure prediction.")
    p.add_argument("--structure-metrics-csv", default=os.environ.get("V6_STRUCTURE_METRICS_CSV", ""), help="Existing structure_metrics.csv to copy into this run.")
    p.add_argument("--proteinmpnn-score-csv", default=os.environ.get("V6_PROTEINMPNN_SCORE_CSV", ""), help="Existing proteinmpnn_scores.csv to copy into this run.")
    p.add_argument("--allow-confidence-proxy", action="store_true")
    p.add_argument("--colabfold-extra-args", default=os.environ.get("COLABFOLD_EXTRA_ARGS", "--num-recycle 3 --num-models 1"))
    p.add_argument("--external-structure-command", default=os.environ.get("EXTERNAL_STRUCTURE_COMMAND", ""))
    p.add_argument("--run-proteinmpnn", action="store_true", default=os.environ.get("RUN_PROTEINMPNN", "0") == "1")
    p.add_argument("--proteinmpnn-dir", default=os.environ.get("PROTEINMPNN_DIR", "/hyperai/home/tools/ProteinMPNN"))
    p.add_argument("--allow-proxy-final", action="store_true")
    p.add_argument("--previous-ranked-csv", default=os.environ.get("V6_PREVIOUS_RANKED_CSV", ""), help="V5/V5.1/V6 ranked CSV used for closed-loop candidate generation.")
    p.add_argument("--feedback-csv", default=os.environ.get("V6_FEEDBACK_CSV", ""), help="Optional explicit feedback CSV used for mutation policy learning.")
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
    context.update({"v6_dir": str(ROOT), "run_dir": str(run_dir), "python": sys.executable})
    (run_dir / "v6_run_context.json").write_text(json.dumps(context, indent=2, ensure_ascii=False), encoding="utf-8")

    status_rows: list[dict[str, Any]] = []
    status_csv = run_dir / "lineage" / "v6_stage_status.csv"
    commands_sh = run_dir / "lineage" / "v6_commands.sh"
    command_lines: list[str] = ["#!/usr/bin/env bash", "set -euo pipefail", f"cd {shlex.quote(str(ROOT))}"]

    try:
        preflight_cmd = [sys.executable, str(ROOT / "deploy" / "00_preflight_v6.py"), "--v5-dir", str(ROOT), "--data-dir", args.data_dir, "--out-dir", str(dirs["preflight"])]
        command_lines.append(shell_quote(preflight_cmd))
        run_stage("00_preflight", preflight_cmd, dirs["logs"] / "00_preflight.log", status_rows, args.dry_run)
        write_status(status_rows, status_csv)

        stage1_cmd = [
            sys.executable, str(ROOT / "scripts" / "run_v6_complete.py"),
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
        command_lines.append(shell_quote(stage1_cmd))
        run_stage("01_stage1_candidate_generation", stage1_cmd, dirs["logs"] / "01_stage1.log", status_rows, args.dry_run)
        write_status(status_rows, status_csv)

        af_fasta = dirs["stage1"] / "af2_af3_structure_tasks.fasta"
        priority_csv = dirs["stage1"] / "structure_priority_top200_v6.csv"
        structure_metrics = dirs["structure"] / "structure_metrics.csv"
        pdb_dir = dirs["structure"] / "predictions"
        pdb_dir.mkdir(parents=True, exist_ok=True)

        if args.structure_metrics_csv:
            # V6 can consume a precomputed structure_metrics.csv as a stable handoff interface.
            import shutil
            shutil.copy2(args.structure_metrics_csv, structure_metrics)
            run_stage("02_structure_metrics_import", ["bash", "-lc", f"wc -l {shlex.quote(str(structure_metrics))}"], dirs["logs"] / "02_structure.log", status_rows, args.dry_run)
            write_status(status_rows, status_csv)
        elif args.reuse_pdb_dir:
            if not args.reference_pdb:
                raise ValueError("--reference-pdb is required when --reuse-pdb-dir is used.")
            reuse_cmd = [
                sys.executable, str(ROOT / "tools" / "reuse_pdb_inputs_v6.py"),
                "--input-pdb-dir", args.reuse_pdb_dir,
                "--priority-csv", str(priority_csv),
                "--reference-pdb", args.reference_pdb,
                "--normalized-pdb-dir", str(pdb_dir),
                "--metrics-csv", str(structure_metrics),
            ]
            command_lines.append(shell_quote(reuse_cmd))
            run_stage("02_reuse_pdb_metrics", reuse_cmd, dirs["logs"] / "02_structure.log", status_rows, args.dry_run)
            write_status(status_rows, status_csv)
        elif args.structure_runner != "skip":
            structure_cmd = [
                sys.executable, str(ROOT / "tools" / "run_structure_backend_v6.py"),
                "--runner", args.structure_runner,
                "--fasta-paths", str(af_fasta),
                "--af-output-dir", str(pdb_dir),
                "--priority-csv", str(priority_csv),
                "--metrics-csv", str(structure_metrics),
                "--work-dir", str(dirs["structure"] / "adapter_work"),
                "--colabfold-extra-args", args.colabfold_extra_args,
            ]
            if args.reference_pdb:
                structure_cmd += ["--reference-pdb", args.reference_pdb]
            if args.allow_confidence_proxy:
                structure_cmd += ["--allow-confidence-proxy"]
            if args.external_structure_command:
                structure_cmd += ["--external-command-template", args.external_structure_command]
            command_lines.append(shell_quote(structure_cmd))
            run_stage("02_structure_prediction_and_metrics", structure_cmd, dirs["logs"] / "02_structure.log", status_rows, args.dry_run)
            write_status(status_rows, status_csv)

        proteinmpnn_scores = dirs["proteinmpnn"] / "proteinmpnn_scores.csv"
        if args.proteinmpnn_score_csv:
            # Existing ProteinMPNN scores are imported by file path so previous score_only runs remain reusable across V6 handoffs.
            import shutil
            shutil.copy2(args.proteinmpnn_score_csv, proteinmpnn_scores)
            run_stage("03_proteinmpnn_scores_import", ["bash", "-lc", f"wc -l {shlex.quote(str(proteinmpnn_scores))}"], dirs["logs"] / "03_proteinmpnn.log", status_rows, args.dry_run)
            write_status(status_rows, status_csv)
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

        final_cmd = [
            sys.executable, str(ROOT / "scripts" / "run_v6_complete.py"),
            "--config", args.config,
            "--data-dir", args.data_dir,
            "--team-name", args.team_name,
            "--feature-mode", args.feature_mode,
            "--esm-model", args.esm_model,
            "--max-train-samples", str(args.max_train_samples),
            "--n-candidates", str(args.n_candidates),
            "--out-dir", str(dirs["final"]),
            "--cache-dir", args.cache_dir,
        ]
        if args.previous_ranked_csv:
            final_cmd += ["--previous-ranked-csv", args.previous_ranked_csv]
        if args.feedback_csv:
            final_cmd += ["--feedback-csv", args.feedback_csv]
        if structure_metrics.exists() or args.dry_run:
            final_cmd += ["--structure-metrics-csv", str(structure_metrics)]
        if proteinmpnn_scores.exists() or args.dry_run:
            final_cmd += ["--proteinmpnn-score-csv", str(proteinmpnn_scores)]
        if args.allow_proxy_final:
            final_cmd += ["--allow-proxy-final"]
        command_lines.append(shell_quote(final_cmd))
        run_stage("04_final_reranking", final_cmd, dirs["logs"] / "04_final.log", status_rows, args.dry_run)
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

    print("V6 run completed")
    print(f"run_dir: {run_dir}")
    print(f"stage_status: {status_csv}")
    print(f"commands: {commands_sh}")
    print(f"submission: {dirs['final'] / 'submission_v6.csv'}")


if __name__ == "__main__":
    main()
