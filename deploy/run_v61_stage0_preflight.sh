#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/env_v61.sh"
RUN_ID="${RUN_ID:-v61_manual}"
RUN_DIR="$V61_OUTPUT_ROOT/$RUN_ID"
mkdir -p "$RUN_DIR/00_preflight"
python "$V61_DIR/deploy/00_preflight_v6.py" --v5-dir "$V61_DIR" --data-dir "$DATA_DIR" --out-dir "$RUN_DIR/00_preflight"
