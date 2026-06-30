#!/usr/bin/env bash
# Minimal environment helper. It reuses the V5.1/V6 environment when present and installs only missing Python deps.
set -euo pipefail
ENV_PREFIX="${V61_ENV_PREFIX:-/hyperai/home/envs/synbio_gfp_v51}"
if [ -f "$ENV_PREFIX/activate_v5.sh" ]; then
  source "$ENV_PREFIX/activate_v5.sh"
elif [ -f "$ENV_PREFIX/bin/activate" ]; then
  source "$ENV_PREFIX/bin/activate"
else
  echo "[install_v61_env] Existing env not found at $ENV_PREFIX. Creating lightweight venv." >&2
  python -m venv "$ENV_PREFIX"
  source "$ENV_PREFIX/bin/activate"
fi
python -m pip install --upgrade pip
python -m pip install -r "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/requirements-v6.txt"
echo "[install_v61_env] done: $ENV_PREFIX"
