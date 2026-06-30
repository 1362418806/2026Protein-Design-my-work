#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/env_v61.sh" >/dev/null
tail -f "$V61_DIR/logs/v61_all_latest.log"
