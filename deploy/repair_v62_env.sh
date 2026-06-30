#!/usr/bin/env bash
# Repair or create the Python environment required by V6.2.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
V62_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
ENV_PREFIX="${V62_ENV_PREFIX:-${V61_ENV_PREFIX:-/hyperai/home/envs/synbio_gfp_v51}}"
REQ="$V62_DIR/requirements-v6.txt"
echo "[repair_v62_env] V62_DIR=$V62_DIR"
echo "[repair_v62_env] ENV_PREFIX=$ENV_PREFIX"
echo "[repair_v62_env] requirements=$REQ"
if [ -x "$ENV_PREFIX/bin/python" ]; then
  PY="$ENV_PREFIX/bin/python"
else
  echo "[repair_v62_env] No python found at $ENV_PREFIX/bin/python; creating env prefix."
  mkdir -p "$(dirname "$ENV_PREFIX")"
  if command -v conda >/dev/null 2>&1; then
    conda create -y -p "$ENV_PREFIX" python=3.10 pip
  elif command -v mamba >/dev/null 2>&1; then
    mamba create -y -p "$ENV_PREFIX" python=3.10 pip
  else
    python -m venv "$ENV_PREFIX"
  fi
  PY="$ENV_PREFIX/bin/python"
fi
"$PY" -m pip install --upgrade pip setuptools wheel
"$PY" -m pip install -r "$REQ"
"$PY" - <<'PYCHECK'
import sys
import pandas, openpyxl, sklearn, yaml
print("python:", sys.executable)
print("pandas:", pandas.__version__)
print("openpyxl:", openpyxl.__version__)
print("sklearn:", sklearn.__version__)
print("yaml ok")
PYCHECK
echo "[repair_v62_env] done"
