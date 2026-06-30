#!/usr/bin/env bash
set -euo pipefail

# Stage 1: train/score/gate candidates and export AF2/AF3 + ProteinMPNN task files.
# No final biological submission is emitted unless real structure metrics are imported later.
PROJECT_ROOT="${PROJECT_ROOT:-/hyperai/home/synbio_gfp_v4_complete}"
DATA_DIR="${DATA_DIR:-/hyperai/input/input0/2026Protein Design}"
OUT_DIR="${OUT_DIR:-${PROJECT_ROOT}/outputs/run_v4_stage1}"
CACHE_DIR="${CACHE_DIR:-${PROJECT_ROOT}/cache_v4}"
TEAM_NAME="${TEAM_NAME:-YourTeamName}"
FEATURE_MODE="${FEATURE_MODE:-esm}"
ESM_MODEL="${ESM_MODEL:-esm2_t30_150M_UR50D}"
MAX_TRAIN_SAMPLES="${MAX_TRAIN_SAMPLES:-20000}"
N_CANDIDATES="${N_CANDIDATES:-20000}"

python "${PROJECT_ROOT}/deploy/00_preflight.py" \
  --project-root "${PROJECT_ROOT}" \
  --data-dir "${DATA_DIR}" \
  $( [[ "${FEATURE_MODE}" == "esm" ]] && printf '%s' '--require-esm' || true )

python "${PROJECT_ROOT}/scripts/run_v4_complete.py" \
  --config "${PROJECT_ROOT}/configs/v4_complete.yaml" \
  --data-dir "${DATA_DIR}" \
  --team-name "${TEAM_NAME}" \
  --feature-mode "${FEATURE_MODE}" \
  --esm-model "${ESM_MODEL}" \
  --max-train-samples "${MAX_TRAIN_SAMPLES}" \
  --n-candidates "${N_CANDIDATES}" \
  --out-dir "${OUT_DIR}" \
  --cache-dir "${CACHE_DIR}"

echo "Stage 1 completed. Structure tasks:"
echo "  ${OUT_DIR}/af2_af3_structure_tasks.fasta"
echo "  ${OUT_DIR}/af2_af3_structure_tasks.jsonl"
echo "  ${OUT_DIR}/proteinmpnn_design_tasks.fasta"
echo "  ${OUT_DIR}/proteinmpnn_design_tasks.jsonl"
echo "Import real metrics in Stage 2 before using submission_v4.csv for competition." 
