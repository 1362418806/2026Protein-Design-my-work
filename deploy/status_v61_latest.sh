#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/env_v61.sh" >/dev/null
PID_FILE="$V61_DIR/logs/v61_all_latest.pid"
RUN_ID_FILE="$V61_DIR/logs/v61_all_latest.run_id"
if [ -f "$PID_FILE" ]; then
  pid="$(cat "$PID_FILE")"
  if kill -0 "$pid" 2>/dev/null; then echo "running pid=$pid"; else echo "not running pid=$pid"; fi
else
  echo "no pid file"
fi
if [ -f "$RUN_ID_FILE" ]; then
  run_id="$(cat "$RUN_ID_FILE")"
  echo "run_id=$run_id"
  status="$V61_OUTPUT_ROOT/$run_id/lineage/v61_stage_status.csv"
  [ -f "$status" ] && cat "$status"
fi
