#!/usr/bin/env bash
# Run only V6 candidate generation and structure-task export.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/env_v6.sh"
RUN_ID="${RUN_ID:-v6_stage1_$(date +%Y%m%d_%H%M%S)}"
OUT_DIR="${V6_OUTPUT_ROOT}/${RUN_ID}/01_stage1"
mkdir -p "${OUT_DIR}"
cmd=("${V6_PYTHON}" "${V6_DIR}/scripts/run_v6_complete.py"
  --config "${V6_DIR}/configs/v6_complete.yaml"
  --data-dir "${DATA_DIR}"
  --team-name "${TEAM_NAME:-YourTeamName}"
  --feature-mode "${FEATURE_MODE:-simple}"
  --esm-model "${ESM_MODEL:-esm2_t30_150M_UR50D}"
  --max-train-samples "${MAX_TRAIN_SAMPLES:-20000}"
  --n-candidates "${N_CANDIDATES:-20000}"
  --out-dir "${OUT_DIR}"
  --cache-dir "${V6_CACHE_DIR}")
[[ -n "${V6_PREVIOUS_RANKED_CSV:-}" ]] && cmd+=(--previous-ranked-csv "${V6_PREVIOUS_RANKED_CSV}")
[[ -n "${V6_FEEDBACK_CSV:-}" ]] && cmd+=(--feedback-csv "${V6_FEEDBACK_CSV}")
"${cmd[@]}"
