#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/env_v61.sh"
: "${RUN_ID:?Set RUN_ID to an existing run directory.}"
RUN_DIR="$V61_OUTPUT_ROOT/$RUN_ID"
python "$V61_DIR/scripts/v6_lineage.py" --run-dir "$RUN_DIR" --out-dir "$RUN_DIR/lineage"
echo "report=$RUN_DIR/04_final/v6_pipeline_report.md"
echo "submission=$RUN_DIR/04_final/submission_v61.csv"
