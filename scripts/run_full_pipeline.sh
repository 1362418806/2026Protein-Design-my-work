#!/usr/bin/env bash
set -euo pipefail

# Portable path variables. Edit these seven variables when moving between machines.
PROJECT_ROOT="${PROJECT_ROOT:-/hyperai/home/synbio_gfp_v3}"
DATA_DIR="${DATA_DIR:-/hyperai/input/input0/2026Protein Design}"
OUT_DIR="${OUT_DIR:-${PROJECT_ROOT}/outputs/run_v3_esm}"
CACHE_DIR="${CACHE_DIR:-${PROJECT_ROOT}/cache}"
TEAM_NAME="${TEAM_NAME:-YourTeamName}"
FEATURE_MODE="${FEATURE_MODE:-esm}"
ESM_MODEL="${ESM_MODEL:-esm2_t30_150M_UR50D}"
MAX_TRAIN_SAMPLES="${MAX_TRAIN_SAMPLES:-20000}"
N_CANDIDATES="${N_CANDIDATES:-20000}"
REFERENCE_FILE="${REFERENCE_FILE:-AAseqs of 5 GFP proteins.txt}"
EXCLUSION_FILE="${EXCLUSION_FILE:-Exclusion_List.csv}"
GFP_DATA_FILE="${GFP_DATA_FILE:-GFP_data.xlsx}"
REFERENCE_PDB="${REFERENCE_PDB:-}"
PREDICTED_PDB_DIR="${PREDICTED_PDB_DIR:-}"

cmd=(
  python "${PROJECT_ROOT}/scripts/run_v3_pipeline.py"
  --config "${PROJECT_ROOT}/configs/v3.yaml"
  --data-dir "${DATA_DIR}"
  --team-name "${TEAM_NAME}"
  --feature-mode "${FEATURE_MODE}"
  --esm-model "${ESM_MODEL}"
  --max-train-samples "${MAX_TRAIN_SAMPLES}"
  --n-candidates "${N_CANDIDATES}"
  --out-dir "${OUT_DIR}"
  --cache-dir "${CACHE_DIR}"
  --references "${REFERENCE_FILE}"
  --exclusion-list "${EXCLUSION_FILE}"
  --gfp-data "${GFP_DATA_FILE}"
)

if [[ -n "${REFERENCE_PDB}" && -n "${PREDICTED_PDB_DIR}" ]]; then
  cmd+=(--reference-pdb "${REFERENCE_PDB}" --predicted-pdb-dir "${PREDICTED_PDB_DIR}")
fi

printf 'Running command:\n%s\n' "${cmd[*]}"
"${cmd[@]}"
