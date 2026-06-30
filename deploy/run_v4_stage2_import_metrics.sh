#!/usr/bin/env bash
set -euo pipefail

# Stage 2: rerun V4 with imported real structure metrics and optional ProteinMPNN scores.
PROJECT_ROOT="${PROJECT_ROOT:-/hyperai/home/synbio_gfp_v4_complete}"
DATA_DIR="${DATA_DIR:-/hyperai/input/input0/2026Protein Design}"
OUT_DIR="${OUT_DIR:-${PROJECT_ROOT}/outputs/run_v4_final}"
CACHE_DIR="${CACHE_DIR:-${PROJECT_ROOT}/cache_v4}"
TEAM_NAME="${TEAM_NAME:-YourTeamName}"
FEATURE_MODE="${FEATURE_MODE:-esm}"
ESM_MODEL="${ESM_MODEL:-esm2_t30_150M_UR50D}"
MAX_TRAIN_SAMPLES="${MAX_TRAIN_SAMPLES:-20000}"
N_CANDIDATES="${N_CANDIDATES:-20000}"
STRUCTURE_METRICS_CSV="${STRUCTURE_METRICS_CSV:?Set STRUCTURE_METRICS_CSV=/path/to/structure_metrics.csv}"
PROTEINMPNN_SCORE_CSV="${PROTEINMPNN_SCORE_CSV:-}"

cmd=(
  python "${PROJECT_ROOT}/scripts/run_v4_complete.py"
  --config "${PROJECT_ROOT}/configs/v4_complete.yaml"
  --data-dir "${DATA_DIR}"
  --team-name "${TEAM_NAME}"
  --feature-mode "${FEATURE_MODE}"
  --esm-model "${ESM_MODEL}"
  --max-train-samples "${MAX_TRAIN_SAMPLES}"
  --n-candidates "${N_CANDIDATES}"
  --out-dir "${OUT_DIR}"
  --cache-dir "${CACHE_DIR}"
  --structure-metrics-csv "${STRUCTURE_METRICS_CSV}"
)
if [[ -n "${PROTEINMPNN_SCORE_CSV}" ]]; then
  cmd+=(--proteinmpnn-score-csv "${PROTEINMPNN_SCORE_CSV}")
fi

printf 'Running Stage 2 command:\n%s\n' "${cmd[*]}"
"${cmd[@]}"
echo "Final outputs: ${OUT_DIR}/final_top6_v4.csv and ${OUT_DIR}/submission_v4.csv"
