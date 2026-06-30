#!/usr/bin/env bash
# Launch a serial, resume-safe V7.2 matrix. The caller may pass additional Python arguments.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/env_v72.sh"
if [ ! -x "$V72_PYTHON" ]; then
  echo "[run_v72_matrix][ERROR] V72_PYTHON is not executable: $V72_PYTHON" >&2
  exit 2
fi
if [ ! -x "$COLABFOLD_BATCH" ]; then
  echo "[run_v72_matrix][ERROR] COLABFOLD_BATCH is not executable: $COLABFOLD_BATCH" >&2
  exit 2
fi
exec "$V72_PYTHON" "$V72_DIR/scripts/run_v72_matrix.py" \
  --matrix-config "$V72_MATRIX_CONFIG" \
  --pipeline-config "$V72_CONFIG" \
  --data-dir "$DATA_DIR" \
  --cache-dir "$V72_CACHE_DIR" \
  --profile "$V72_MATRIX_PROFILE" \
  --reference-pdb "$REF_PDB" \
  --colabfold-batch "$COLABFOLD_BATCH" \
  --structure-runner "$STRUCTURE_RUNNER" \
  "$@"
