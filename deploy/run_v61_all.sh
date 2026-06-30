#!/usr/bin/env bash
# Start the complete V6.1 Adaptive pipeline in the background with a fixed latest log.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/env_v61.sh"
LOG_DIR="$V61_DIR/logs"
mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/v61_all_latest.log"
PID_FILE="$LOG_DIR/v61_all_latest.pid"
RUN_ID_FILE="$LOG_DIR/v61_all_latest.run_id"
RUN_ID="${RUN_ID:-v61_$(date +%Y%m%d_%H%M%S)}"
echo "$RUN_ID" > "$RUN_ID_FILE"
: > "$LOG_FILE"

cmd=(python "$V61_DIR/scripts/run_v61_pipeline.py"
  --data-dir "$DATA_DIR"
  --team-name "${TEAM_NAME:-YourTeamName}"
  --run-id "$RUN_ID"
  --feature-mode "${FEATURE_MODE:-simple}"
  --esm-model "${ESM_MODEL:-esm2_t30_150M_UR50D}"
  --max-train-samples "${MAX_TRAIN_SAMPLES:-20000}"
  --n-candidates "${N_CANDIDATES:-20000}"
  --config "${V61_CONFIG:-$V61_DIR/configs/v61_adaptive.yaml}"
  --cache-dir "$V61_CACHE_DIR"
  --output-root "$V61_OUTPUT_ROOT"
  --gpu-profile "${GPU_PROFILE:-auto}"
  --structure-policy "${STRUCTURE_POLICY:-reuse-first}"
  --structure-runner "${STRUCTURE_RUNNER:-auto}"
  --reference-pdb "${REF_PDB:-$V61_DIR/reference/reference_gfp.pdb}"
  --pdb-search-roots "${PDB_SEARCH_ROOTS:-}"
  --colabfold-extra-args "${COLABFOLD_EXTRA_ARGS:---num-models 1 --num-recycle 1}"
)
if [ "${RUN_PROTEINMPNN:-1}" = "1" ]; then cmd+=(--run-proteinmpnn); fi
if [ -n "${V61_PREVIOUS_RANKED_CSV:-}" ]; then cmd+=(--previous-ranked-csv "$V61_PREVIOUS_RANKED_CSV"); fi
if [ -n "${V61_FEEDBACK_CSV:-}" ]; then cmd+=(--feedback-csv "$V61_FEEDBACK_CSV"); fi
if [ -n "${V61_BRIGHTNESS_TEACHER_CSV:-}" ]; then cmd+=(--brightness-teacher-csv "$V61_BRIGHTNESS_TEACHER_CSV"); fi

nohup "${cmd[@]}" > "$LOG_FILE" 2>&1 &
echo $! > "$PID_FILE"
echo "[run_v61_all] started pid=$(cat "$PID_FILE") run_id=$RUN_ID"
echo "[run_v61_all] log=$LOG_FILE"
echo "[run_v61_all] tail -f $LOG_FILE"
