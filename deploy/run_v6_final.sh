#!/usr/bin/env bash
# Run V6 final reranking after structure/proteinmpnn metrics are available.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/env_v6.sh"
RUN_ID="${RUN_ID:?Set RUN_ID to the existing V6 run directory name.}"
STRUCTURE_METRICS="${STRUCTURE_METRICS:-${V6_OUTPUT_ROOT}/${RUN_ID}/02_structure/structure_metrics.csv}"
PMPNN_SCORES="${PMPNN_SCORES:-${V6_OUTPUT_ROOT}/${RUN_ID}/03_proteinmpnn/proteinmpnn_scores.csv}"
OUT_DIR="${V6_OUTPUT_ROOT}/${RUN_ID}/04_final"
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
[[ -f "${STRUCTURE_METRICS}" ]] && cmd+=(--structure-metrics-csv "${STRUCTURE_METRICS}")
[[ -f "${PMPNN_SCORES}" ]] && cmd+=(--proteinmpnn-score-csv "${PMPNN_SCORES}")
[[ -n "${ALLOW_PROXY_FINAL:-}" ]] && cmd+=(--allow-proxy-final)
"${cmd[@]}"
