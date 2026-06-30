#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/env_v61.sh"
: "${RUN_ID:?Set RUN_ID to an existing run directory.}"
RUN_DIR="$V61_OUTPUT_ROOT/$RUN_ID"
cmd=(python "$V61_DIR/scripts/run_v61_complete.py" --config "${V61_CONFIG:-$V61_DIR/configs/v61_adaptive.yaml}" --data-dir "$DATA_DIR" --team-name "${TEAM_NAME:-YourTeamName}" --feature-mode "${FEATURE_MODE:-simple}" --esm-model "${ESM_MODEL:-esm2_t30_150M_UR50D}" --max-train-samples "${MAX_TRAIN_SAMPLES:-20000}" --n-candidates "${N_CANDIDATES:-20000}" --out-dir "$RUN_DIR/04_final" --cache-dir "$V61_CACHE_DIR" --structure-metrics-csv "$RUN_DIR/02_structure/structure_metrics.csv")
[ -f "$RUN_DIR/03_proteinmpnn/proteinmpnn_scores.csv" ] && cmd+=(--proteinmpnn-score-csv "$RUN_DIR/03_proteinmpnn/proteinmpnn_scores.csv")
[ -n "${V61_PREVIOUS_RANKED_CSV:-}" ] && cmd+=(--previous-ranked-csv "$V61_PREVIOUS_RANKED_CSV")
[ -n "${V61_FEEDBACK_CSV:-}" ] && cmd+=(--feedback-csv "$V61_FEEDBACK_CSV")
[ -n "${V61_BRIGHTNESS_TEACHER_CSV:-}" ] && cmd+=(--brightness-teacher-csv "$V61_BRIGHTNESS_TEACHER_CSV")
"${cmd[@]}"
