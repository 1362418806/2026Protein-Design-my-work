#!/usr/bin/env bash
# Start the complete V7 Orthogonal-Validation pipeline in the background with a fixed latest log.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/env_v7.sh"
if [ "${V7_ENV_READY:-0}" != "1" ]; then
  echo "[run_v7_all.sh][ERROR] V7 Python environment is not ready." >&2
  exit 2
fi
PYTHON_BIN="${V7_PYTHON:-python}"
LOG_DIR="$V7_DIR/logs"
mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/v7_all_latest.log"
PID_FILE="$LOG_DIR/v7_all_latest.pid"
RUN_ID_FILE="$LOG_DIR/v7_all_latest.run_id"
RUN_ID="${RUN_ID:-v7_$(date +%Y%m%d_%H%M%S)}"
echo "$RUN_ID" > "$RUN_ID_FILE"
: > "$LOG_FILE"
cmd=("$PYTHON_BIN" "$V7_DIR/scripts/run_v7_pipeline.py"
  --data-dir "$DATA_DIR"
  --team-name "${TEAM_NAME:-YourTeamName}"
  --run-id "$RUN_ID"
  --feature-mode "${FEATURE_MODE:-simple}"
  --esm-model "${ESM_MODEL:-esm2_t30_150M_UR50D}"
  --max-train-samples "${MAX_TRAIN_SAMPLES:-20000}"
  --n-candidates "${N_CANDIDATES:-20000}"
  --config "${V7_CONFIG:-$V7_DIR/configs/v7_adaptive.yaml}"
  --cache-dir "$V7_CACHE_DIR"
  --output-root "$V7_OUTPUT_ROOT"
  --gpu-profile "${GPU_PROFILE:-auto}"
  --structure-policy "${STRUCTURE_POLICY:-reuse-first}"
  --structure-runner "${STRUCTURE_RUNNER:-auto}"
  --reference-pdb "${REF_PDB:-$V7_DIR/reference/reference_gfp.pdb}"
  --pdb-search-roots "${PDB_SEARCH_ROOTS:-}"
  --colabfold-extra-args "${COLABFOLD_EXTRA_ARGS:---num-models 1 --num-recycle 1}"
  --structure-chunk-size "${STRUCTURE_CHUNK_SIZE:-12}"
  --colabfold-msa-mode "${V7_COLABFOLD_MSA_MODE:-mmseqs2_uniref_env}"
  --chunk-timeout-sec "${V7_COLABFOLD_CHUNK_TIMEOUT_SEC:-7200}"
  --chunk-retries "${V7_COLABFOLD_CHUNK_RETRIES:-1}"
)
[ -n "${COLABFOLD_BATCH:-}" ] && cmd+=(--colabfold-batch "${COLABFOLD_BATCH}")
[ -n "${EXTERNAL_STRUCTURE_COMMAND:-}" ] && cmd+=(--external-structure-command "${EXTERNAL_STRUCTURE_COMMAND}")
if [ "${RUN_PROTEINMPNN:-1}" = "1" ]; then cmd+=(--run-proteinmpnn); fi
if [ "${RUN_V7_SOFTCHECK:-0}" = "1" ]; then
  cmd+=(--run-softcheck --softcheck-runner "${V7_SOFTCHECK_RUNNER:-metrics-only}" --softcheck-top-n "${V7_SOFTCHECK_TOP_N:-80}" --softcheck-timeout-sec "${V7_SOFTCHECK_TIMEOUT_SEC:-21600}")
fi
[ -n "${V7_SOFTCHECK_CSV:-}" ] && cmd+=(--softcheck-csv "${V7_SOFTCHECK_CSV}")
[ -n "${V7_SOFTCHECK_PDB_DIR:-}" ] && cmd+=(--softcheck-pdb-dir "${V7_SOFTCHECK_PDB_DIR}")
[ -n "${OPENCOMPLEX2_COMMAND:-${V7_SOFTCHECK_COMMAND:-}}" ] && cmd+=(--softcheck-command "${OPENCOMPLEX2_COMMAND:-${V7_SOFTCHECK_COMMAND:-}}")
[ -n "${V7_PREVIOUS_RANKED_CSV:-${V62_PREVIOUS_RANKED_CSV:-${V61_PREVIOUS_RANKED_CSV:-}}}" ] && cmd+=(--previous-ranked-csv "${V7_PREVIOUS_RANKED_CSV:-${V62_PREVIOUS_RANKED_CSV:-${V61_PREVIOUS_RANKED_CSV:-}}}")
[ -n "${V7_FEEDBACK_CSV:-${V62_FEEDBACK_CSV:-${V61_FEEDBACK_CSV:-}}}" ] && cmd+=(--feedback-csv "${V7_FEEDBACK_CSV:-${V62_FEEDBACK_CSV:-${V61_FEEDBACK_CSV:-}}}")
[ -n "${V7_BRIGHTNESS_TEACHER_CSV:-${V62_BRIGHTNESS_TEACHER_CSV:-${V61_BRIGHTNESS_TEACHER_CSV:-}}}" ] && cmd+=(--brightness-teacher-csv "${V7_BRIGHTNESS_TEACHER_CSV:-${V62_BRIGHTNESS_TEACHER_CSV:-${V61_BRIGHTNESS_TEACHER_CSV:-}}}")
nohup "${cmd[@]}" > "$LOG_FILE" 2>&1 &
echo $! > "$PID_FILE"
echo "[run_v7_all] started pid=$(cat "$PID_FILE") run_id=$RUN_ID"
echo "[run_v7_all] log=$LOG_FILE"
echo "[run_v7_all] tail -f $LOG_FILE"
