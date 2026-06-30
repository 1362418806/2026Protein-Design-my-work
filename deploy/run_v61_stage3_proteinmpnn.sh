#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/env_v61.sh"
: "${RUN_ID:?Set RUN_ID to an existing run directory.}"
RUN_DIR="$V61_OUTPUT_ROOT/$RUN_ID"
python "$V61_DIR/tools/run_proteinmpnn_v6.py" --pdb-dir "$RUN_DIR/02_structure/predictions" --out-csv "$RUN_DIR/03_proteinmpnn/proteinmpnn_scores.csv" --work-dir "$RUN_DIR/03_proteinmpnn/score_work" --proteinmpnn-dir "${PROTEINMPNN_DIR:-/hyperai/home/tools/ProteinMPNN}" --priority-csv "$RUN_DIR/01_stage1/structure_priority_top200_v6.csv" --chain-id auto
