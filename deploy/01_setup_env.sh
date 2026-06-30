#!/usr/bin/env bash
set -euo pipefail

# This script installs the V4 Python stack only. AF2/AF3/ColabFold are excluded
# by design because they are frequently cluster-specific and should be provided
# as an external structure-metrics producer.
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
ENV_NAME="${ENV_NAME:-synbio-gfp-v4}"
PYTHON_BIN="${PYTHON_BIN:-python}"

if command -v conda >/dev/null 2>&1; then
  if ! conda env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then
    conda create -y -n "${ENV_NAME}" python=3.10 pip
  fi
  # shellcheck disable=SC1091
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate "${ENV_NAME}"
  PYTHON_BIN="python"
fi

"${PYTHON_BIN}" -m pip install --upgrade pip setuptools wheel
"${PYTHON_BIN}" -m pip install -r "${PROJECT_ROOT}/requirements-v4.txt"
"${PYTHON_BIN}" "${PROJECT_ROOT}/deploy/00_preflight.py" --project-root "${PROJECT_ROOT}" --require-esm
