#!/usr/bin/env python3
"""Run an isolated V7.2 evidence matrix.

V7.2 deliberately keeps each experiment in its own directory. Candidate sources are generated
once, copied into each experiment as a frozen input, and then tested under controlled structure
settings. Matrix scores are reporting-only; no V7.2 score is fed back into candidate generation.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shlex
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def shell_quote(command: list[str]) -> str:
    return " ".join(shlex.quote(str(token)) for token in command)


def load_yaml(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        value = yaml.safe_load(handle) or {}
    if not isinstance(value, dict):
        raise ValueError(f"Matrix config must be a mapping: {path}")
    return value


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")


def run_command(command: list[str], log_path: Path, dry_run: bool = False) -> None:
    """Run one stage with a persistent log and immediate terminal streaming."""
    log_path.parent.mkdir(parents=True, exist_ok=True)
    line = shell_quote(command)
    print(f"$ {line}", flush=True)
    with log_path.open("a", encoding="utf-8", buffering=1) as log:
        log.write(f"$ {line}\n")
        if dry_run:
            log.write("[dry_run] skipped\n")
            return
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        assert process.stdout is not None
        for output_line in process.stdout:
            print(output_line, end="", flush=True)
            log.write(output_line)
        code = process.wait()
        log.write(f"[exit_code] {code}\n")
        if code != 0:
            raise RuntimeError(f"Command failed with exit code {code}: {line}")


def stage_complete(path: Path) -> bool:
    return path.is_file() and path.stat().st_size > 0


def infer_stage1_artifacts(stage1_dir: Path) -> tuple[Path, Path]:
    universe = stage1_dir / "candidate_universe_v62.csv"
    priority = stage1_dir / "structure_priority_top200_v6.csv"
    if not universe.is_file():
        raise FileNotFoundError(f"Missing candidate universe: {universe}")
    if not priority.is_file():
        raise FileNotFoundError(f"Missing structure priority: {priority}")
    return universe, priority


def resolve_experiment_ids(config: dict[str, Any], requested: str) -> list[dict[str, Any]]:
    experiments = config.get("experiments", [])
    if not isinstance(experiments, list):
        raise ValueError("experiments must be a list")
    if not requested:
        return [item for item in experiments if bool(item.get("enabled", True))]
    requested_ids = {value.strip() for value in requested.split(",") if value.strip()}
    available = {str(item.get("id", "")) for item in experiments}
    unknown = requested_ids.difference(available)
    if unknown:
        raise ValueError(f"Unknown experiment ids: {sorted(unknown)}")
    return [item for item in experiments if str(item.get("id", "")) in requested_ids]


def source_config(config: dict[str, Any], source_id: str) -> dict[str, Any]:
    sources = config.get("candidate_sources", {})
    if source_id not in sources:
        raise KeyError(f"Unknown candidate_source={source_id!r}")
    value = sources[source_id]
    if not isinstance(value, dict):
        raise ValueError(f"candidate_sources.{source_id} must be a mapping")
    return value


def run_candidate_source(
    *,
    matrix_root: Path,
    source_id: str,
    source: dict[str, Any],
    args: argparse.Namespace,
    dry_run: bool,
) -> tuple[Path, Path]:
    """Generate one candidate source once and freeze V7.2-named source artifacts."""
    source_dir = matrix_root / "candidate_sources" / source_id
    stage1_dir = source_dir / "01_stage1"
    canonical_universe = stage1_dir / "candidate_universe_v72.csv"
    canonical_priority = stage1_dir / "structure_priority_top200_v72.csv"
    metadata_path = source_dir / "candidate_source_manifest_v72.json"

    if args.resume and stage_complete(canonical_universe) and stage_complete(canonical_priority) and stage_complete(metadata_path):
        return canonical_universe, canonical_priority

    stage1_dir.mkdir(parents=True, exist_ok=True)
    command = [
        sys.executable,
        str(ROOT / "scripts" / "run_v62_complete.py"),
        "--config",
        str(args.pipeline_config),
        "--data-dir",
        str(args.data_dir),
        "--team-name",
        str(args.team_name),
        "--feature-mode",
        str(source.get("feature_mode", "simple")),
        "--esm-model",
        str(source.get("esm_model", "esm2_t30_150M_UR50D")),
        "--max-train-samples",
        str(source.get("max_train_samples", args.max_train_samples)),
        "--n-candidates",
        str(source.get("n_candidates", args.n_candidates)),
        "--out-dir",
        str(stage1_dir),
        "--cache-dir",
        str(args.cache_dir),
    ]
    run_command(command, source_dir / "logs" / "01_stage1.log", dry_run)
    legacy_universe, legacy_priority = infer_stage1_artifacts(stage1_dir)
    # The generator is inherited from V6.2, but matrix contracts expose only V7.2 artifact names.
    shutil.copy2(legacy_universe, canonical_universe)
    shutil.copy2(legacy_priority, canonical_priority)
    manifest = {
        "schema_version": "v72-candidate-source-1",
        "candidate_source": source_id,
        "created_at": utc_now(),
        "source_parameters": source,
        "pipeline_config": str(args.pipeline_config),
        "candidate_universe": str(canonical_universe),
        "candidate_universe_sha256": sha256_file(canonical_universe) if not dry_run else "",
        "structure_priority": str(canonical_priority),
        "structure_priority_sha256": sha256_file(canonical_priority) if not dry_run else "",
        "generator_compat_inputs": {
            "legacy_universe": str(legacy_universe),
            "legacy_priority": str(legacy_priority),
        },
    }
    write_json(metadata_path, manifest)
    return canonical_universe, canonical_priority

def prepare_experiment_inputs(
    *,
    run_dir: Path,
    source_id: str,
    universe_path: Path,
    priority_path: Path,
    top_n: int,
    dry_run: bool,
) -> tuple[Path, Path]:
    """Copy source artifacts so each experiment owns an immutable candidate snapshot."""
    input_dir = run_dir / "01_stage1"
    input_dir.mkdir(parents=True, exist_ok=True)
    local_universe = input_dir / "candidate_universe_v72.csv"
    local_priority = input_dir / f"structure_priority_top{top_n}_v72.csv"

    if not (dry_run and local_universe.exists()):
        shutil.copy2(universe_path, local_universe)
    priority = pd.read_csv(priority_path)
    priority.head(top_n).to_csv(local_priority, index=False)

    manifest = {
        "schema_version": "v72-experiment-input-1",
        "candidate_source": source_id,
        "source_universe": str(universe_path),
        "source_universe_sha256": sha256_file(universe_path) if not dry_run else "",
        "source_priority": str(priority_path),
        "source_priority_sha256": sha256_file(priority_path) if not dry_run else "",
        "local_universe": str(local_universe),
        "local_universe_sha256": sha256_file(local_universe) if not dry_run else "",
        "local_priority": str(local_priority),
        "local_priority_sha256": sha256_file(local_priority) if not dry_run else "",
        "structure_top_n": int(top_n),
        "created_at": utc_now(),
    }
    write_json(input_dir / "input_manifest_v72.json", manifest)
    return local_universe, local_priority


def run_experiment(
    *,
    matrix_root: Path,
    experiment: dict[str, Any],
    universe_path: Path,
    priority_path: Path,
    args: argparse.Namespace,
    dry_run: bool,
) -> None:
    """Run one serial, resume-safe matrix experiment."""
    experiment_id = str(experiment["id"])
    candidate_source = str(experiment["candidate_source"])
    run_dir = matrix_root / "runs" / experiment_id
    final_dir = run_dir / "05_final"
    final_submission = final_dir / "submission_v72_for_upload.csv"
    audit_json = final_dir / "v72_submission_audit.json"

    if args.resume and stage_complete(final_submission) and stage_complete(audit_json):
        try:
            audit = json.loads(audit_json.read_text(encoding="utf-8"))
            if bool(audit.get("passed", False)):
                print(f"[resume] experiment already passed: {experiment_id}", flush=True)
                return
        except Exception:
            pass

    top_n = int(experiment.get("structure_top_n", 200))
    local_universe, local_priority = prepare_experiment_inputs(
        run_dir=run_dir,
        source_id=candidate_source,
        universe_path=universe_path,
        priority_path=priority_path,
        top_n=top_n,
        dry_run=dry_run,
    )

    structure_dir = run_dir / "02_structure"
    prediction_dir = structure_dir / "predictions"
    raw_metrics = structure_dir / "structure_metrics_raw.csv"
    v72_metrics = structure_dir / "structure_metrics_v72.csv"
    structure_dir.mkdir(parents=True, exist_ok=True)

    recycle = int(experiment.get("num_recycle", 1))
    models = int(experiment.get("num_models", 1))
    resolver_command = [
        sys.executable,
        str(ROOT / "tools" / "resolve_structures_v62.py"),
        "--priority-csv",
        str(local_priority),
        "--reference-pdb",
        str(args.reference_pdb),
        "--prediction-dir",
        str(prediction_dir),
        "--metrics-csv",
        str(raw_metrics),
        "--work-dir",
        str(structure_dir / "resolver"),
        "--pdb-search-roots",
        "",
        "--structure-policy",
        "predict-missing",
        "--runner",
        str(experiment.get("structure_runner", args.structure_runner)),
        "--colabfold-extra-args",
        f"--num-models {models} --num-recycle {recycle}",
        "--colabfold-msa-mode",
        str(experiment["msa_mode"]),
        "--chunk-size",
        str(experiment.get("chunk_size", args.chunk_size)),
        "--chunk-timeout-sec",
        str(experiment.get("chunk_timeout_sec", args.chunk_timeout_sec)),
        "--chunk-retries",
        str(experiment.get("chunk_retries", args.chunk_retries)),
        "--retry-backoff-sec",
        str(experiment.get("retry_backoff_sec", args.retry_backoff_sec)),
    ]
    if args.colabfold_batch:
        resolver_command += ["--colabfold-batch", str(args.colabfold_batch)]
    run_command(resolver_command, run_dir / "logs" / "02_structure.log", dry_run)

    globalfit_command = [
        sys.executable,
        str(ROOT / "tools" / "augment_structure_metrics_v72.py"),
        "--metrics-csv",
        str(raw_metrics),
        "--reference-pdb",
        str(args.reference_pdb),
        "--out-csv",
        str(v72_metrics),
    ]
    run_command(globalfit_command, run_dir / "logs" / "02_globalfit.log", dry_run)

    pmpnn_dir = run_dir / "03_proteinmpnn"
    pmpnn_scores = pmpnn_dir / "proteinmpnn_scores.csv"
    if bool(experiment.get("run_proteinmpnn", args.run_proteinmpnn)):
        pmpnn_command = [
            sys.executable,
            str(ROOT / "tools" / "run_proteinmpnn_v6.py"),
            "--pdb-dir",
            str(prediction_dir),
            "--out-csv",
            str(pmpnn_scores),
            "--work-dir",
            str(pmpnn_dir / "score_work"),
            "--proteinmpnn-dir",
            str(args.proteinmpnn_dir),
            "--priority-csv",
            str(local_priority),
            "--chain-id",
            "auto",
        ]
        run_command(pmpnn_command, run_dir / "logs" / "03_proteinmpnn.log", dry_run)

    softcheck_path = ""
    if bool(experiment.get("run_softcheck", False)):
        soft_dir = run_dir / "04_softcheck"
        softcheck_csv = soft_dir / "softcheck_metrics_v72.csv"
        soft_command = [
            sys.executable,
            str(ROOT / "tools" / "run_opencomplex2_softcheck_v7.py"),
            "--candidates-csv",
            str(local_priority),
            "--reference-pdb",
            str(args.reference_pdb),
            "--out-csv",
            str(softcheck_csv),
            "--work-dir",
            str(soft_dir / "work"),
            "--runner",
            str(experiment.get("softcheck_runner", "metrics-only")),
            "--top-n",
            str(experiment.get("softcheck_top_n", min(80, top_n))),
            "--timeout-sec",
            str(experiment.get("softcheck_timeout_sec", args.softcheck_timeout_sec)),
        ]
        if args.softcheck_command:
            soft_command += ["--external-command-template", args.softcheck_command]
        run_command(soft_command, run_dir / "logs" / "04_softcheck.log", dry_run)
        softcheck_path = str(softcheck_csv)

    final_command = [
        sys.executable,
        str(ROOT / "scripts" / "run_v72_final_from_manifest.py"),
        "--config",
        str(args.pipeline_config),
        "--frozen-candidates-csv",
        str(local_universe),
        "--structure-metrics-csv",
        str(v72_metrics),
        "--out-dir",
        str(final_dir),
        "--team-name",
        str(args.team_name),
        "--data-dir",
        str(args.data_dir),
        "--exclusion-list-csv",
        str(args.data_dir / "Exclusion_List.csv"),
    ]
    if stage_complete(pmpnn_scores):
        final_command += ["--proteinmpnn-score-csv", str(pmpnn_scores)]
    if softcheck_path:
        final_command += ["--softcheck-csv", softcheck_path]
    run_command(final_command, run_dir / "logs" / "05_final.log", dry_run)

    audit_command = [
        sys.executable,
        str(ROOT / "tools" / "audit_submission_v72.py"),
        "--submission-csv",
        str(final_submission),
        "--final-csv",
        str(final_dir / "final_top6_v72.csv"),
        "--structure-metrics-csv",
        str(v72_metrics),
        "--data-dir",
        str(args.data_dir),
        "--out-json",
        str(final_dir / "v72_submission_audit.json"),
        "--out-md",
        str(final_dir / "v72_submission_audit.md"),
        "--fail-on-error",
    ]
    run_command(audit_command, run_dir / "logs" / "06_audit.log", dry_run)

    manifest = {
        "schema_version": "v72-experiment-1",
        "experiment_id": experiment_id,
        "candidate_source": candidate_source,
        "diagnostic_only": bool(experiment.get("diagnostic_only", False)),
        "parameters": experiment,
        "created_at": utc_now(),
        "input_universe": str(local_universe),
        "input_universe_sha256": sha256_file(local_universe) if not dry_run else "",
        "input_priority": str(local_priority),
        "input_priority_sha256": sha256_file(local_priority) if not dry_run else "",
        "raw_structure_metrics": str(raw_metrics),
        "structure_metrics_v72": str(v72_metrics),
        "proteinmpnn_scores": str(pmpnn_scores) if stage_complete(pmpnn_scores) else "",
        "softcheck_scores": softcheck_path,
        "final_submission": str(final_submission),
        "final_top6": str(final_dir / "final_top6_v72.csv"),
        "audit_json": str(audit_json),
        "code_sha256": {
            "matrix_runner": sha256_file(Path(__file__)),
            "final_selector": sha256_file(ROOT / "scripts" / "run_v72_final_from_manifest.py"),
            "globalfit": sha256_file(ROOT / "tools" / "augment_structure_metrics_v72.py"),
            "audit": sha256_file(ROOT / "tools" / "audit_submission_v72.py"),
        },
    }
    write_json(run_dir / "v72_experiment_manifest.json", manifest)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the isolated SynBio GFP V7.2 matrix.")
    parser.add_argument("--matrix-config", default=str(ROOT / "configs" / "v72_matrix.yaml"))
    parser.add_argument("--pipeline-config", default=str(ROOT / "configs" / "v72_adaptive.yaml"))
    parser.add_argument("--data-dir", default=os.environ.get("DATA_DIR", "/input0/2026Protein Design"))
    parser.add_argument("--team-name", default=os.environ.get("TEAM_NAME", "YourTeamName"))
    parser.add_argument("--matrix-root", default="")
    parser.add_argument("--experiments", default="", help="Comma-separated experiment ids. Overrides --profile when provided.")
    parser.add_argument("--profile", default=os.environ.get("V72_MATRIX_PROFILE", "deadline_core"), help="Named profile from the matrix config.")
    parser.add_argument("--resume", action="store_true", default=True)
    parser.add_argument("--no-resume", dest="resume", action="store_false")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--max-train-samples", type=int, default=int(os.environ.get("MAX_TRAIN_SAMPLES", "20000")))
    parser.add_argument("--n-candidates", type=int, default=int(os.environ.get("N_CANDIDATES", "20000")))
    parser.add_argument("--cache-dir", default=os.environ.get("V72_CACHE_DIR", str(ROOT / "cache" / "v72")))
    parser.add_argument("--reference-pdb", default=os.environ.get("REF_PDB", str(ROOT / "reference" / "reference_gfp.pdb")))
    parser.add_argument("--structure-runner", default=os.environ.get("STRUCTURE_RUNNER", "auto"))
    parser.add_argument("--colabfold-batch", default=os.environ.get("COLABFOLD_BATCH", ""))
    parser.add_argument("--chunk-size", type=int, default=int(os.environ.get("STRUCTURE_CHUNK_SIZE", "12")))
    parser.add_argument("--chunk-timeout-sec", type=int, default=int(os.environ.get("V72_COLABFOLD_CHUNK_TIMEOUT_SEC", "7200")))
    parser.add_argument("--chunk-retries", type=int, default=int(os.environ.get("V72_COLABFOLD_CHUNK_RETRIES", "1")))
    parser.add_argument("--retry-backoff-sec", type=int, default=int(os.environ.get("V72_COLABFOLD_RETRY_BACKOFF_SEC", "60")))
    parser.add_argument("--run-proteinmpnn", action="store_true", default=os.environ.get("RUN_PROTEINMPNN", "1") == "1")
    parser.add_argument("--proteinmpnn-dir", default=os.environ.get("PROTEINMPNN_DIR", "/hyperai/home/tools/ProteinMPNN"))
    parser.add_argument("--softcheck-timeout-sec", type=int, default=int(os.environ.get("V72_SOFTCHECK_TIMEOUT_SEC", "1800")))
    parser.add_argument("--softcheck-command", default=os.environ.get("V72_SOFTCHECK_COMMAND", ""))
    args = parser.parse_args()
    args.data_dir = Path(args.data_dir)
    args.pipeline_config = Path(args.pipeline_config)
    args.reference_pdb = Path(args.reference_pdb)

    matrix_config_path = Path(args.matrix_config)
    config = load_yaml(matrix_config_path)
    requested_experiments = args.experiments
    if not requested_experiments and args.profile:
        profiles = config.get("profiles", {})
        if args.profile not in profiles:
            raise ValueError(f"Unknown matrix profile: {args.profile}")
        requested_experiments = ",".join(str(value) for value in profiles[args.profile])
    selected = resolve_experiment_ids(config, requested_experiments)
    if not selected:
        raise ValueError("No enabled experiments selected.")

    if args.matrix_root:
        matrix_root = Path(args.matrix_root)
    else:
        matrix_root = ROOT / "experiments" / f"v72_matrix_{time.strftime('%Y%m%d_%H%M%S')}"
    matrix_root.mkdir(parents=True, exist_ok=True)

    plan = {
        "schema_version": "v72-matrix-plan-1",
        "created_at": utc_now(),
        "matrix_root": str(matrix_root),
        "matrix_config": str(matrix_config_path),
        "pipeline_config": str(args.pipeline_config),
        "data_dir": str(args.data_dir),
        "reference_pdb": str(args.reference_pdb),
        "profile": args.profile,
        "selected_experiments": selected,
        "matrix_config_sha256": sha256_file(matrix_config_path),
        "pipeline_config_sha256": sha256_file(args.pipeline_config),
    }
    write_json(matrix_root / "matrix_plan_resolved_v72.json", plan)
    if args.dry_run:
        print(json.dumps({"dry_run": True, "matrix_root": str(matrix_root), "selected_experiments": selected}, indent=2))
        return

    required_sources: list[str] = []
    for experiment in selected:
        source_id = str(experiment["candidate_source"])
        if source_id not in required_sources:
            required_sources.append(source_id)

    sources: dict[str, tuple[Path, Path]] = {}
    for source_id in required_sources:
        sources[source_id] = run_candidate_source(
            matrix_root=matrix_root,
            source_id=source_id,
            source=source_config(config, source_id),
            args=args,
            dry_run=args.dry_run,
        )

    for experiment in selected:
        source_id = str(experiment["candidate_source"])
        universe, priority = sources[source_id]
        run_experiment(
            matrix_root=matrix_root,
            experiment=experiment,
            universe_path=universe,
            priority_path=priority,
            args=args,
            dry_run=args.dry_run,
        )

    summary_command = [
        sys.executable,
        str(ROOT / "scripts" / "summarize_v72_matrix.py"),
        "--matrix-root",
        str(matrix_root),
    ]
    run_command(summary_command, matrix_root / "logs" / "99_matrix_summary.log", args.dry_run)
    print(f"V7.2 matrix completed: {matrix_root}")


if __name__ == "__main__":
    main()
