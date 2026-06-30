#!/usr/bin/env bash
# Recompute V7.2 summary files for an existing matrix root without rerunning experiments.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/env_v72.sh"
if [ "$#" -ne 1 ]; then
  echo "Usage: $0 /absolute/path/to/v72_matrix_YYYYMMDD_HHMMSS" >&2
  exit 2
fi
exec "$V72_PYTHON" "$V72_DIR/scripts/summarize_v72_matrix.py" --matrix-root "$1"
