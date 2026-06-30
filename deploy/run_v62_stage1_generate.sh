#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/env_v62.sh"
if [ "${V62_ENV_READY:-0}" != "1" ]; then
  echo "[run_v62_stage1_generate.sh][ERROR] V6.2 Python environment is not ready. Run: bash deploy/repair_v62_env.sh" >&2
  exit 2
fi
PYTHON_BIN="${V62_PYTHON:-python}"
RUN_ID="${RUN_ID:-v62_manual}"
RUN_DIR="$V62_OUTPUT_ROOT/$RUN_ID"
mkdir -p "$RUN_DIR/01_stage1"
cmd=("$PYTHON_BIN" "$V62_DIR/scripts/run_v62_complete.py" --config "${V62_CONFIG:-$V62_DIR/configs/v62_adaptive.yaml}" --data-dir "$DATA_DIR" --team-name "${TEAM_NAME:-YourTeamName}" --feature-mode "${FEATURE_MODE:-simple}" --esm-model "${ESM_MODEL:-esm2_t30_150M_UR50D}" --max-train-samples "${MAX_TRAIN_SAMPLES:-20000}" --n-candidates "${N_CANDIDATES:-20000}" --out-dir "$RUN_DIR/01_stage1" --cache-dir "$V62_CACHE_DIR")
[ -n "${V62_PREVIOUS_RANKED_CSV:-${V61_PREVIOUS_RANKED_CSV:-}}" ] && cmd+=(--previous-ranked-csv "${V62_PREVIOUS_RANKED_CSV:-${V61_PREVIOUS_RANKED_CSV:-}}")
[ -n "${V62_FEEDBACK_CSV:-${V61_FEEDBACK_CSV:-}}" ] && cmd+=(--feedback-csv "${V62_FEEDBACK_CSV:-${V61_FEEDBACK_CSV:-}}")
[ -n "${V62_BRIGHTNESS_TEACHER_CSV:-${V61_BRIGHTNESS_TEACHER_CSV:-}}" ] && cmd+=(--brightness-teacher-csv "${V62_BRIGHTNESS_TEACHER_CSV:-${V61_BRIGHTNESS_TEACHER_CSV:-}}")
"${cmd[@]}"
