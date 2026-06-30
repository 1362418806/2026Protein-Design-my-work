#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/env_v7.sh"
if [ "${V7_ENV_READY:-0}" != "1" ]; then
  echo "[run_v7_stage1_generate.sh][ERROR] V7 Python environment is not ready. Run: bash deploy/env_v7.sh" >&2
  exit 2
fi
PYTHON_BIN="${V7_PYTHON:-python}"
RUN_ID="${RUN_ID:-v7_manual}"
RUN_DIR="$V7_OUTPUT_ROOT/$RUN_ID"
mkdir -p "$RUN_DIR/01_stage1"
cmd=("$PYTHON_BIN" "$V7_DIR/scripts/run_v62_complete.py" --config "${V7_CONFIG:-$V7_DIR/configs/v7_adaptive.yaml}" --data-dir "$DATA_DIR" --team-name "${TEAM_NAME:-YourTeamName}" --feature-mode "${FEATURE_MODE:-simple}" --esm-model "${ESM_MODEL:-esm2_t30_150M_UR50D}" --max-train-samples "${MAX_TRAIN_SAMPLES:-20000}" --n-candidates "${N_CANDIDATES:-20000}" --out-dir "$RUN_DIR/01_stage1" --cache-dir "$V7_CACHE_DIR")
[ -n "${V7_PREVIOUS_RANKED_CSV:-${V61_PREVIOUS_RANKED_CSV:-}}" ] && cmd+=(--previous-ranked-csv "${V7_PREVIOUS_RANKED_CSV:-${V61_PREVIOUS_RANKED_CSV:-}}")
[ -n "${V7_FEEDBACK_CSV:-${V61_FEEDBACK_CSV:-}}" ] && cmd+=(--feedback-csv "${V7_FEEDBACK_CSV:-${V61_FEEDBACK_CSV:-}}")
[ -n "${V7_BRIGHTNESS_TEACHER_CSV:-${V61_BRIGHTNESS_TEACHER_CSV:-}}" ] && cmd+=(--brightness-teacher-csv "${V7_BRIGHTNESS_TEACHER_CSV:-${V61_BRIGHTNESS_TEACHER_CSV:-}}")
"${cmd[@]}"
