#!/usr/bin/env bash
set -euo pipefail

# Debug only: emits proxy final Top-6 without AF/PDB metrics. Do not use this as final submission.
PROJECT_ROOT="${PROJECT_ROOT:-/hyperai/home/synbio_gfp_v4_complete}"
DATA_DIR="${DATA_DIR:-/hyperai/input/input0/2026Protein Design}"
OUT_DIR="${OUT_DIR:-${PROJECT_ROOT}/outputs/run_v4_proxy_debug}"
CACHE_DIR="${CACHE_DIR:-${PROJECT_ROOT}/cache_v4}"
TEAM_NAME="${TEAM_NAME:-YourTeamName}"
FEATURE_MODE="${FEATURE_MODE:-simple}"
ESM_MODEL="${ESM_MODEL:-esm2_t30_150M_UR50D}"
MAX_TRAIN_SAMPLES="${MAX_TRAIN_SAMPLES:-20000}"
N_CANDIDATES="${N_CANDIDATES:-20000}"

python "${PROJECT_ROOT}/scripts/run_v4_complete.py" \
  --config "${PROJECT_ROOT}/configs/v4_complete.yaml" \
  --data-dir "${DATA_DIR}" \
  --team-name "${TEAM_NAME}" \
  --feature-mode "${FEATURE_MODE}" \
  --esm-model "${ESM_MODEL}" \
  --max-train-samples "${MAX_TRAIN_SAMPLES}" \
  --n-candidates "${N_CANDIDATES}" \
  --out-dir "${OUT_DIR}" \
  --cache-dir "${CACHE_DIR}" \
  --allow-proxy-final
