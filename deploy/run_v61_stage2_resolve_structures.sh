#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/env_v61.sh"
: "${RUN_ID:?Set RUN_ID to an existing run directory.}"
RUN_DIR="$V61_OUTPUT_ROOT/$RUN_ID"
PRIORITY_CSV="${PRIORITY_CSV:-$RUN_DIR/01_stage1/structure_priority_top200_v6.csv}"
python "$V61_DIR/tools/resolve_structures_v61.py" --priority-csv "$PRIORITY_CSV" --reference-pdb "${REF_PDB:-$V61_DIR/reference/reference_gfp.pdb}" --prediction-dir "$RUN_DIR/02_structure/predictions" --metrics-csv "$RUN_DIR/02_structure/structure_metrics.csv" --work-dir "$RUN_DIR/02_structure" --pdb-search-roots "${PDB_SEARCH_ROOTS:-}" --structure-policy "${STRUCTURE_POLICY:-reuse-first}" --runner "${STRUCTURE_RUNNER:-auto}" --colabfold-extra-args "${COLABFOLD_EXTRA_ARGS:---num-models 1 --num-recycle 1}"
