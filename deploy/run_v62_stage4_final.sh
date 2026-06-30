#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/env_v62.sh"
if [ "${V62_ENV_READY:-0}" != "1" ]; then
  echo "[run_v62_stage4_final.sh][ERROR] V6.2 Python environment is not ready. Run: bash deploy/repair_v62_env.sh" >&2
  exit 2
fi
PYTHON_BIN="${V62_PYTHON:-python}"
: "${RUN_ID:?Set RUN_ID to an existing run directory.}"
RUN_DIR="$V62_OUTPUT_ROOT/$RUN_ID"
FROZEN="${FROZEN_CANDIDATES_CSV:-$RUN_DIR/01_stage1/candidate_universe_v62.csv}"
[ -f "$FROZEN" ] || FROZEN="$RUN_DIR/01_stage1/all_ranked_candidates_v6.csv"
cmd=("$PYTHON_BIN" "$V62_DIR/scripts/run_v62_final_from_manifest.py" --config "${V62_CONFIG:-$V62_DIR/configs/v62_adaptive.yaml}" --frozen-candidates-csv "$FROZEN" --structure-metrics-csv "$RUN_DIR/02_structure/structure_metrics.csv" --out-dir "$RUN_DIR/04_final" --team-name "${TEAM_NAME:-YourTeamName}")
[ -f "$RUN_DIR/03_proteinmpnn/proteinmpnn_scores.csv" ] && cmd+=(--proteinmpnn-score-csv "$RUN_DIR/03_proteinmpnn/proteinmpnn_scores.csv")
"${cmd[@]}"
