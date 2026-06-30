#!/usr/bin/env bash
# Start the complete V6.2 Evidence-Safe pipeline in the background with a fixed latest log.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/env_v62.sh"
if [ "${V62_ENV_READY:-0}" != "1" ]; then
  echo "[run_v62_all.sh][ERROR] V6.2 Python environment is not ready. Run: bash deploy/repair_v62_env.sh" >&2
  exit 2
fi
PYTHON_BIN="${V62_PYTHON:-python}"
LOG_DIR="$V62_DIR/logs"
mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/v62_all_latest.log"
PID_FILE="$LOG_DIR/v62_all_latest.pid"
RUN_ID_FILE="$LOG_DIR/v62_all_latest.run_id"
RUN_ID="${RUN_ID:-v62_$(date +%Y%m%d_%H%M%S)}"
echo "$RUN_ID" > "$RUN_ID_FILE"
: > "$LOG_FILE"
cmd=("$PYTHON_BIN" "$V62_DIR/scripts/run_v62_pipeline.py"
  --data-dir "$DATA_DIR"
  --team-name "${TEAM_NAME:-YourTeamName}"
  --run-id "$RUN_ID"
  --feature-mode "${FEATURE_MODE:-simple}"
  --esm-model "${ESM_MODEL:-esm2_t30_150M_UR50D}"
  --max-train-samples "${MAX_TRAIN_SAMPLES:-20000}"
  --n-candidates "${N_CANDIDATES:-20000}"
  --config "${V62_CONFIG:-$V62_DIR/configs/v62_adaptive.yaml}"
  --cache-dir "$V62_CACHE_DIR"
  --output-root "$V62_OUTPUT_ROOT"
  --gpu-profile "${GPU_PROFILE:-auto}"
  --structure-policy "${STRUCTURE_POLICY:-reuse-first}"
  --structure-runner "${STRUCTURE_RUNNER:-auto}"
  --reference-pdb "${REF_PDB:-$V62_DIR/reference/reference_gfp.pdb}"
  --pdb-search-roots "${PDB_SEARCH_ROOTS:-}"
  --colabfold-extra-args "${COLABFOLD_EXTRA_ARGS:---num-models 1 --num-recycle 1}"
  --structure-chunk-size "${STRUCTURE_CHUNK_SIZE:-12}"
  --colabfold-msa-mode "${V62_COLABFOLD_MSA_MODE:-mmseqs2_uniref_env}"
  --chunk-timeout-sec "${V62_COLABFOLD_CHUNK_TIMEOUT_SEC:-7200}"
  --chunk-retries "${V62_COLABFOLD_CHUNK_RETRIES:-1}"
)
[ -n "${COLABFOLD_BATCH:-}" ] && cmd+=(--colabfold-batch "${COLABFOLD_BATCH}")
if [ "${RUN_PROTEINMPNN:-1}" = "1" ]; then cmd+=(--run-proteinmpnn); fi
[ -n "${V62_PREVIOUS_RANKED_CSV:-${V61_PREVIOUS_RANKED_CSV:-}}" ] && cmd+=(--previous-ranked-csv "${V62_PREVIOUS_RANKED_CSV:-${V61_PREVIOUS_RANKED_CSV:-}}")
[ -n "${V62_FEEDBACK_CSV:-${V61_FEEDBACK_CSV:-}}" ] && cmd+=(--feedback-csv "${V62_FEEDBACK_CSV:-${V61_FEEDBACK_CSV:-}}")
[ -n "${V62_BRIGHTNESS_TEACHER_CSV:-${V61_BRIGHTNESS_TEACHER_CSV:-}}" ] && cmd+=(--brightness-teacher-csv "${V62_BRIGHTNESS_TEACHER_CSV:-${V61_BRIGHTNESS_TEACHER_CSV:-}}")
nohup "${cmd[@]}" > "$LOG_FILE" 2>&1 &
echo $! > "$PID_FILE"
echo "[run_v62_all] started pid=$(cat "$PID_FILE") run_id=$RUN_ID"
echo "[run_v62_all] log=$LOG_FILE"
echo "[run_v62_all] tail -f $LOG_FILE"
