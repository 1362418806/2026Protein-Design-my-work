#!/usr/bin/env bash
# Follow the fixed-name V5 background log.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/env_v6.sh"
if [[ ! -f "${V6_LATEST_LOG}" ]]; then
  echo "V5 latest log does not exist yet: ${V6_LATEST_LOG}" >&2
  exit 1
fi
tail -f "${V6_LATEST_LOG}"
