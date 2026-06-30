#!/usr/bin/env bash
# Launch the complete traceable V6 pipeline in the background by default.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/env_v6.sh"

RUN_ID="${RUN_ID:-v6_$(date +%Y%m%d_%H%M%S)}"
TEAM_NAME="${TEAM_NAME:-YourTeamName}"
FEATURE_MODE="${FEATURE_MODE:-simple}"
STRUCTURE_RUNNER="${STRUCTURE_RUNNER:-metrics-only}"
RUN_PROTEINMPNN="${RUN_PROTEINMPNN:-0}"
V6_BACKGROUND="${V6_BACKGROUND:-1}"
V6_RESET_LATEST_LOG="${V6_RESET_LATEST_LOG:-1}"

if [[ "${V6_IN_BACKGROUND:-0}" != "1" && "${V6_BACKGROUND}" == "1" ]]; then
  mkdir -p "${V6_FIXED_LOG_DIR}"
  if [[ "${V6_RESET_LATEST_LOG}" == "1" ]]; then
    : > "${V6_LATEST_LOG}"
  fi

  cat > "${V6_LATEST_ENV}" <<ENV
RUN_ID=${RUN_ID}
TEAM_NAME=${TEAM_NAME}
FEATURE_MODE=${FEATURE_MODE}
STRUCTURE_RUNNER=${STRUCTURE_RUNNER}
RUN_PROTEINMPNN=${RUN_PROTEINMPNN}
DATA_DIR=${DATA_DIR}
V6_DIR=${V6_DIR}
V6_OUTPUT_ROOT=${V6_OUTPUT_ROOT}
V6_CACHE_DIR=${V6_CACHE_DIR}
REF_PDB=${REF_PDB:-}
V6_REUSE_PDB_DIR=${V6_REUSE_PDB_DIR:-}
V6_STRUCTURE_METRICS_CSV=${V6_STRUCTURE_METRICS_CSV:-}
V6_PROTEINMPNN_SCORE_CSV=${V6_PROTEINMPNN_SCORE_CSV:-}
V6_PREVIOUS_RANKED_CSV=${V6_PREVIOUS_RANKED_CSV:-}
V6_FEEDBACK_CSV=${V6_FEEDBACK_CSV:-}
ENV
  printf '%s\n' "${RUN_ID}" > "${V6_LATEST_RUN_ID_FILE}"

  # The child process does the real work; nohup keeps it alive after SSH disconnects.
  nohup env \
    V6_IN_BACKGROUND=1 \
    RUN_ID="${RUN_ID}" \
    TEAM_NAME="${TEAM_NAME}" \
    FEATURE_MODE="${FEATURE_MODE}" \
    STRUCTURE_RUNNER="${STRUCTURE_RUNNER}" \
    RUN_PROTEINMPNN="${RUN_PROTEINMPNN}" \
    REF_PDB="${REF_PDB:-}" \
    V6_REUSE_PDB_DIR="${V6_REUSE_PDB_DIR:-}" \
    V6_STRUCTURE_METRICS_CSV="${V6_STRUCTURE_METRICS_CSV:-}" \
    V6_PROTEINMPNN_SCORE_CSV="${V6_PROTEINMPNN_SCORE_CSV:-}" \
    V6_PREVIOUS_RANKED_CSV="${V6_PREVIOUS_RANKED_CSV:-}" \
    V6_FEEDBACK_CSV="${V6_FEEDBACK_CSV:-}" \
    ALLOW_CONFIDENCE_PROXY="${ALLOW_CONFIDENCE_PROXY:-}" \
    ALLOW_PROXY_FINAL="${ALLOW_PROXY_FINAL:-}" \
    MAX_TRAIN_SAMPLES="${MAX_TRAIN_SAMPLES:-20000}" \
    N_CANDIDATES="${N_CANDIDATES:-20000}" \
    ESM_MODEL="${ESM_MODEL:-esm2_t30_150M_UR50D}" \
    COLABFOLD_EXTRA_ARGS="${COLABFOLD_EXTRA_ARGS:---num-recycle 3 --num-models 1}" \
    EXTERNAL_STRUCTURE_COMMAND="${EXTERNAL_STRUCTURE_COMMAND:-}" \
    bash "${BASH_SOURCE[0]}" >> "${V6_LATEST_LOG}" 2>&1 &

  child_pid=$!
  printf '%s\n' "${child_pid}" > "${V6_LATEST_PID}"
  echo "[V6] launched in background"
  echo "[V6] pid: ${child_pid}"
  echo "[V6] run_id: ${RUN_ID}"
  echo "[V6] fixed log: ${V6_LATEST_LOG}"
  echo "[V6] follow log: tail -f ${V6_LATEST_LOG}"
  echo "[V6] status csv: ${V6_OUTPUT_ROOT}/${RUN_ID}/lineage/v6_stage_status.csv"
  exit 0
fi

cmd=("${V6_PYTHON}" -u "${V6_DIR}/scripts/run_v6_pipeline.py"
  --data-dir "${DATA_DIR}"
  --team-name "${TEAM_NAME}"
  --run-id "${RUN_ID}"
  --feature-mode "${FEATURE_MODE}"
  --esm-model "${ESM_MODEL:-esm2_t30_150M_UR50D}"
  --max-train-samples "${MAX_TRAIN_SAMPLES:-20000}"
  --n-candidates "${N_CANDIDATES:-20000}"
  --cache-dir "${V6_CACHE_DIR}"
  --output-root "${V6_OUTPUT_ROOT}"
  --structure-runner "${STRUCTURE_RUNNER}")

[[ -n "${REF_PDB:-}" ]] && cmd+=(--reference-pdb "${REF_PDB}")
[[ -n "${V6_REUSE_PDB_DIR:-}" ]] && cmd+=(--reuse-pdb-dir "${V6_REUSE_PDB_DIR}")
[[ -n "${V6_STRUCTURE_METRICS_CSV:-}" ]] && cmd+=(--structure-metrics-csv "${V6_STRUCTURE_METRICS_CSV}")
[[ -n "${V6_PROTEINMPNN_SCORE_CSV:-}" ]] && cmd+=(--proteinmpnn-score-csv "${V6_PROTEINMPNN_SCORE_CSV}")
[[ -n "${V6_PREVIOUS_RANKED_CSV:-}" ]] && cmd+=(--previous-ranked-csv "${V6_PREVIOUS_RANKED_CSV}")
[[ -n "${V6_FEEDBACK_CSV:-}" ]] && cmd+=(--feedback-csv "${V6_FEEDBACK_CSV}")
[[ -n "${ALLOW_CONFIDENCE_PROXY:-}" ]] && cmd+=(--allow-confidence-proxy)
[[ -n "${ALLOW_PROXY_FINAL:-}" ]] && cmd+=(--allow-proxy-final)
[[ "${RUN_PROTEINMPNN}" == "1" ]] && cmd+=(--run-proteinmpnn)
[[ -n "${EXTERNAL_STRUCTURE_COMMAND:-}" ]] && cmd+=(--external-structure-command "${EXTERNAL_STRUCTURE_COMMAND}")

echo "[V6] foreground worker started"
echo "[V6] run_id: ${RUN_ID}"
echo "[V6] fixed log: ${V6_LATEST_LOG}"
echo "[V6] command: ${cmd[*]}"
"${cmd[@]}"
