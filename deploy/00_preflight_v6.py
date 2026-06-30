#!/usr/bin/env python3
"""Preflight checks for SynBio GFP V6 deployments."""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

REQUIRED_DATA_FILES = [
    "AAseqs of 5 GFP proteins.txt",
    "Exclusion_List.csv",
    "GFP_data.xlsx",
]
REQUIRED_PY_MODULES = ["numpy", "pandas", "scipy", "sklearn", "joblib", "yaml", "openpyxl", "tqdm", "Bio", "esm"]


def sha256_file(path: Path, block_size: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(block_size), b""):
            h.update(chunk)
    return h.hexdigest()


def command_output(cmd: list[str]) -> str:
    try:
        out = subprocess.check_output(cmd, stderr=subprocess.STDOUT, text=True, timeout=15)
        return out.strip()
    except Exception as exc:
        return f"ERROR:{type(exc).__name__}:{exc}"


def add(rows: list[dict[str, Any]], name: str, status: str, detail: str, path: str = "") -> None:
    rows.append({"check": name, "status": status, "detail": detail, "path": path})


def main() -> None:
    p = argparse.ArgumentParser(description="Run V6 environment and data preflight checks.")
    p.add_argument("--v5-dir", default=os.environ.get("V6_DIR", "/hyperai/home/synbio_gfp_v6_complete"))
    p.add_argument("--data-dir", default=os.environ.get("DATA_DIR", "/hyperai/input/input0/2026Protein Design"))
    p.add_argument("--out-dir", default=None)
    args = p.parse_args()

    v5_dir = Path(args.v5_dir)
    data_dir = Path(args.data_dir)
    out_dir = Path(args.out_dir or (v5_dir / "outputs" / "preflight_v6"))
    out_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, Any]] = []
    add(rows, "v5_dir_exists", "pass" if v5_dir.exists() else "fail", str(v5_dir.exists()), str(v5_dir))
    add(rows, "data_dir_exists", "pass" if data_dir.exists() else "fail", str(data_dir.exists()), str(data_dir))

    # Record the exact interpreter so dependency issues cannot be hidden by PATH ordering.
    expected_prefix = os.environ.get("V6_ENV_PREFIX", "/hyperai/home/envs/synbio_gfp_v6")
    add(rows, "python_executable", "pass" if sys.executable.startswith(expected_prefix) else "warn", sys.executable, sys.executable)

    for rel in REQUIRED_DATA_FILES:
        path = data_dir / rel
        if path.exists():
            add(rows, f"data_file:{rel}", "pass", f"size={path.stat().st_size};sha256={sha256_file(path)[:16]}", str(path))
        else:
            add(rows, f"data_file:{rel}", "fail", "missing", str(path))

    for rel in ["scripts/run_v6_complete.py", "scripts/run_v6_pipeline.py", "configs/v6_complete.yaml", "requirements-v6.txt"]:
        path = v5_dir / rel
        add(rows, f"project_file:{rel}", "pass" if path.exists() else "fail", "present" if path.exists() else "missing", str(path))

    for mod in REQUIRED_PY_MODULES:
        add(rows, f"python_module:{mod}", "pass" if importlib.util.find_spec(mod) else "warn", "importable" if importlib.util.find_spec(mod) else "not found")

    for exe in ["python", "colabfold_batch", "git", "nvidia-smi"]:
        found = shutil.which(exe)
        add(rows, f"executable:{exe}", "pass" if found else "warn", found or "not found", found or "")

    add(rows, "python_version", "pass", sys.version.replace("\n", " "))
    add(rows, "nvidia_smi", "info", command_output(["nvidia-smi", "--query-gpu=name,driver_version,compute_cap,memory.total", "--format=csv,noheader"]))

    with (out_dir / "v6_preflight_checks.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["check", "status", "detail", "path"])
        writer.writeheader()
        writer.writerows(rows)

    payload = {"rows": rows, "v5_dir": str(v5_dir), "data_dir": str(data_dir)}
    (out_dir / "v6_preflight_checks.json").write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    failures = [r for r in rows if r["status"] == "fail"]
    print(f"Wrote preflight results to {out_dir}")
    if failures:
        print(f"Preflight failures: {len(failures)}", file=sys.stderr)
        sys.exit(2)


if __name__ == "__main__":
    main()
