#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/env_v62.sh"
if [ "${V62_ENV_READY:-0}" != "1" ]; then
  echo "[run_v62_report.sh][ERROR] V6.2 Python environment is not ready. Run: bash deploy/repair_v62_env.sh" >&2
  exit 2
fi
PYTHON_BIN="${V62_PYTHON:-python}"
: "${RUN_ID:?Set RUN_ID to an existing run directory.}"
RUN_DIR="$V62_OUTPUT_ROOT/$RUN_ID"
"$PYTHON_BIN" "$V62_DIR/scripts/v6_lineage.py" --run-dir "$RUN_DIR" --out-dir "$RUN_DIR/lineage"
echo "report=$RUN_DIR/04_final/v6_pipeline_report.md"
echo "submission=$RUN_DIR/04_final/submission_v62.csv"
