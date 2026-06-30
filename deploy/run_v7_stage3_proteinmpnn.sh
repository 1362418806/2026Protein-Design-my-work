#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/env_v7.sh"
if [ "${V7_ENV_READY:-0}" != "1" ]; then
  echo "[run_v7_stage3_proteinmpnn.sh][ERROR] V6.2 Python environment is not ready. Run: bash deploy/repair_v7_env.sh" >&2
  exit 2
fi
PYTHON_BIN="${V7_PYTHON:-python}"
: "${RUN_ID:?Set RUN_ID to an existing run directory.}"
RUN_DIR="$V7_OUTPUT_ROOT/$RUN_ID"
"$PYTHON_BIN" "$V7_DIR/tools/run_proteinmpnn_v6.py" --pdb-dir "$RUN_DIR/02_structure/predictions" --out-csv "$RUN_DIR/03_proteinmpnn/proteinmpnn_scores.csv" --work-dir "$RUN_DIR/03_proteinmpnn/score_work" --proteinmpnn-dir "${PROTEINMPNN_DIR:-/hyperai/home/tools/ProteinMPNN}" --priority-csv "$RUN_DIR/01_stage1/structure_priority_top200_v6.csv" --chain-id auto
