#!/usr/bin/env python3
from __future__ import annotations

import argparse
import importlib
import json
import os
import platform
import sys
from pathlib import Path
from typing import Any

REQUIRED_IMPORTS = {
    "numpy": "numpy",
    "pandas": "pandas",
    "sklearn": "scikit-learn",
    "yaml": "pyyaml",
    "openpyxl": "openpyxl",
    "joblib": "joblib",
}
OPTIONAL_IMPORTS = {
    "torch": "torch",
    "esm": "fair-esm",
    "Bio": "biopython",
}
DATA_FILES = [
    "AAseqs of 5 GFP proteins.txt",
    "Exclusion_List.csv",
    "GFP_data.xlsx",
]


def check_import(module: str, package: str) -> dict[str, Any]:
    try:
        mod = importlib.import_module(module)
        version = getattr(mod, "__version__", "unknown")
        return {"package": package, "module": module, "ok": True, "version": version, "error": ""}
    except Exception as exc:  # pragma: no cover - preflight is intentionally diagnostic.
        return {"package": package, "module": module, "ok": False, "version": "", "error": repr(exc)}


def main() -> None:
    parser = argparse.ArgumentParser(description="Preflight check for SynBio GFP V4 deployment")
    parser.add_argument("--project-root", default=str(Path(__file__).resolve().parents[1]))
    parser.add_argument("--data-dir", default="")
    parser.add_argument("--require-esm", action="store_true")
    parser.add_argument("--json-out", default="")
    args = parser.parse_args()

    project_root = Path(args.project_root).expanduser().resolve()
    data_dir = Path(args.data_dir).expanduser().resolve() if args.data_dir else None
    report: dict[str, Any] = {
        "python": sys.version,
        "platform": platform.platform(),
        "project_root": str(project_root),
        "data_dir": str(data_dir) if data_dir else "",
        "imports_required": [check_import(k, v) for k, v in REQUIRED_IMPORTS.items()],
        "imports_optional": [check_import(k, v) for k, v in OPTIONAL_IMPORTS.items()],
        "data_files": [],
        "cuda": {"available": False, "device_count": 0, "devices": []},
        "status": "ok",
        "actions": [],
    }

    try:
        import torch
        report["cuda"] = {
            "available": bool(torch.cuda.is_available()),
            "device_count": int(torch.cuda.device_count()),
            "devices": [torch.cuda.get_device_name(i) for i in range(torch.cuda.device_count())],
        }
    except Exception as exc:
        report["cuda"]["error"] = repr(exc)

    if data_dir:
        for name in DATA_FILES:
            p = data_dir / name
            report["data_files"].append({"name": name, "path": str(p), "exists": p.exists()})

    missing_required = [x for x in report["imports_required"] if not x["ok"]]
    missing_esm = [x for x in report["imports_optional"] if x["module"] == "esm" and not x["ok"]]
    missing_data = [x for x in report["data_files"] if not x.get("exists", True)]

    if missing_required:
        report["status"] = "fail"
        report["actions"].append("Install core dependencies with: pip install -r requirements-v4.txt")
    if args.require_esm and missing_esm:
        report["status"] = "fail"
        report["actions"].append("Install FAIR-ESM with: pip install fair-esm")
    if args.require_esm and not report["cuda"].get("available"):
        report["actions"].append("CUDA is not visible. ESM mode will run slowly on CPU; use simple mode for smoke tests.")
    if missing_data:
        report["status"] = "fail"
        report["actions"].append("Check --data-dir and the official file names.")

    text = json.dumps(report, indent=2, ensure_ascii=False)
    print(text)
    if args.json_out:
        Path(args.json_out).write_text(text, encoding="utf-8")
    if report["status"] != "ok":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
