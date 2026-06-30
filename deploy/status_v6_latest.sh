#!/usr/bin/env bash
# Print the latest V5 background PID and stage status table if available.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/env_v6.sh"
RUN_ID="${RUN_ID:-}"
if [[ -z "${RUN_ID}" && -f "${V6_LATEST_RUN_ID_FILE}" ]]; then
  RUN_ID="$(cat "${V6_LATEST_RUN_ID_FILE}")"
fi
if [[ -f "${V6_LATEST_PID}" ]]; then
  pid="$(cat "${V6_LATEST_PID}")"
  if ps -p "${pid}" >/dev/null 2>&1; then
    echo "latest pid: ${pid} (running)"
  else
    echo "latest pid: ${pid} (not running)"
  fi
else
  echo "latest pid: missing"
fi
if [[ -n "${RUN_ID}" ]]; then
  echo "latest run_id: ${RUN_ID}"
  status_csv="${V6_OUTPUT_ROOT}/${RUN_ID}/lineage/v6_stage_status.csv"
  echo "status csv: ${status_csv}"
  [[ -f "${status_csv}" ]] && cat "${status_csv}" || true
fi
echo "fixed log: ${V6_LATEST_LOG}"
