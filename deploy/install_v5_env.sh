#!/usr/bin/env bash
# Compatibility wrapper for older commands. V6 requirements live in requirements-v6.txt.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
echo "[install_v6] deploy/install_v5_env.sh is deprecated for V6; forwarding to install_v6_env.sh" >&2
exec bash "${SCRIPT_DIR}/install_v6_env.sh" "$@"
