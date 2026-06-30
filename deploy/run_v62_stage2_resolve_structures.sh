#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/env_v62.sh"
if [ "${V62_ENV_READY:-0}" != "1" ]; then
  echo "[run_v62_stage2_resolve_structures.sh][ERROR] V6.2 Python environment is not ready. Run: bash deploy/repair_v62_env.sh" >&2
  exit 2
fi
PYTHON_BIN="${V62_PYTHON:-python}"
: "${RUN_ID:?Set RUN_ID to an existing run directory.}"
RUN_DIR="$V62_OUTPUT_ROOT/$RUN_ID"
PRIORITY_CSV="${PRIORITY_CSV:-$RUN_DIR/01_stage1/structure_priority_top200_v6.csv}"
cmd=("$PYTHON_BIN" "$V62_DIR/tools/resolve_structures_v62.py"
  --priority-csv "$PRIORITY_CSV"
  --reference-pdb "${REF_PDB:-$V62_DIR/reference/reference_gfp.pdb}"
  --prediction-dir "$RUN_DIR/02_structure/predictions"
  --metrics-csv "$RUN_DIR/02_structure/structure_metrics.csv"
  --work-dir "$RUN_DIR/02_structure"
  --pdb-search-roots "${PDB_SEARCH_ROOTS:-}"
  --structure-policy "${STRUCTURE_POLICY:-reuse-first}"
  --runner "${STRUCTURE_RUNNER:-auto}"
  --chunk-size "${STRUCTURE_CHUNK_SIZE:-12}"
  --chunk-timeout-sec "${V62_COLABFOLD_CHUNK_TIMEOUT_SEC:-7200}"
  --chunk-retries "${V62_COLABFOLD_CHUNK_RETRIES:-1}"
  --retry-backoff-sec "${V62_COLABFOLD_RETRY_BACKOFF_SEC:-60}"
  --colabfold-msa-mode "${V62_COLABFOLD_MSA_MODE:-mmseqs2_uniref_env}"
  --colabfold-extra-args "${COLABFOLD_EXTRA_ARGS:---num-models 1 --num-recycle 1}"
)
[ -n "${COLABFOLD_BATCH:-}" ] && cmd+=(--colabfold-batch "$COLABFOLD_BATCH")
"${cmd[@]}"
