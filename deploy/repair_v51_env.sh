#!/usr/bin/env bash
# Repair/install only the V6 Python dependencies into the fixed V6 env.
# This script does not touch candidate generation, ranking, PDBs, or previous outputs.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/env_v6.sh"

if [[ ! -x "${V6_PYTHON}" ]]; then
  echo "[repair_v6][ERROR] Missing V6 Python: ${V6_PYTHON}" >&2
  echo "Run: bash ${SCRIPT_DIR}/install_v6_env.sh" >&2
  exit 2
fi

REQ="${V6_DIR}/requirements-v6.txt"
if [[ ! -f "${REQ}" ]]; then
  echo "[repair_v6][ERROR] Missing requirements file: ${REQ}" >&2
  exit 2
fi

# Install through the exact interpreter to avoid contaminating base conda or the V5 env.
"${V6_PYTHON}" -m pip install --upgrade pip setuptools wheel
"${V6_PYTHON}" -m pip install --cache-dir "${PIP_CACHE_DIR}" -r "${REQ}"

# Verify all Stage 1 imports before the pipeline spends time on candidate generation.
"${V6_PYTHON}" - <<'PY'
import sys
mods = ["numpy", "pandas", "scipy", "sklearn", "joblib", "yaml", "openpyxl", "tqdm", "Bio", "esm"]
print("python:", sys.executable)
for mod in mods:
    __import__(mod)
    print("import ok:", mod)
PY
