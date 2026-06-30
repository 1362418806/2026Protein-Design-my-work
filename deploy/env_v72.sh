#!/usr/bin/env bash
# V7.2 isolated matrix environment. Safe to source from an interactive SSH shell.
_v72_echo_err() { echo "$*" >&2; }
_v72_project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export V72_DIR="${V72_DIR:-/hyperai/home/synbio_gfp_v7_2_complete}"
if [ ! -d "$V72_DIR" ]; then export V72_DIR="$_v72_project_dir"; fi
export DATA_DIR="${DATA_DIR:-/input0/2026Protein Design}"
export V72_ENV_PREFIX="${V72_ENV_PREFIX:-/hyperai/home/envs/synbio_gfp_v51}"
export V72_PYTHON="${V72_PYTHON:-$V72_ENV_PREFIX/bin/python}"
if [ ! -x "$V72_PYTHON" ]; then
  export V72_PYTHON="$(command -v python 2>/dev/null || true)"
fi
export PATH="$V72_ENV_PREFIX/bin:${PATH:-}"
export VIRTUAL_ENV="$V72_ENV_PREFIX"
export PYTHONPATH="$V72_DIR/src:${PYTHONPATH:-}"

export V72_OUTPUT_ROOT="${V72_OUTPUT_ROOT:-$V72_DIR/outputs}"
export V72_CACHE_DIR="${V72_CACHE_DIR:-$V72_DIR/cache/v72}"
export V72_CONFIG="${V72_CONFIG:-$V72_DIR/configs/v72_adaptive.yaml}"
export V72_MATRIX_CONFIG="${V72_MATRIX_CONFIG:-$V72_DIR/configs/v72_matrix.yaml}"
export V72_MATRIX_PROFILE="${V72_MATRIX_PROFILE:-deadline_core}"
export REF_PDB="${REF_PDB:-$V72_DIR/reference/reference_gfp.pdb}"

export FEATURE_MODE="${FEATURE_MODE:-simple}"
export STRUCTURE_RUNNER="${STRUCTURE_RUNNER:-auto}"
export STRUCTURE_CHUNK_SIZE="${STRUCTURE_CHUNK_SIZE:-12}"
export V72_COLABFOLD_CHUNK_TIMEOUT_SEC="${V72_COLABFOLD_CHUNK_TIMEOUT_SEC:-7200}"
export V72_COLABFOLD_CHUNK_RETRIES="${V72_COLABFOLD_CHUNK_RETRIES:-1}"
export V72_COLABFOLD_RETRY_BACKOFF_SEC="${V72_COLABFOLD_RETRY_BACKOFF_SEC:-60}"
export COLABFOLD_BATCH="${COLABFOLD_BATCH:-/hyperai/home/tools/localcolabfold/.pixi/envs/default/bin/colabfold_batch}"
export RUN_PROTEINMPNN="${RUN_PROTEINMPNN:-1}"
export PROTEINMPNN_DIR="${PROTEINMPNN_DIR:-/hyperai/home/tools/ProteinMPNN}"
export V72_SOFTCHECK_TIMEOUT_SEC="${V72_SOFTCHECK_TIMEOUT_SEC:-1800}"

mkdir -p "$V72_OUTPUT_ROOT" "$V72_CACHE_DIR" "$V72_DIR/logs" "$V72_DIR/experiments"
echo "[env_v72] V72_DIR=$V72_DIR"
echo "[env_v72] DATA_DIR=$DATA_DIR"
echo "[env_v72] python=${V72_PYTHON:-MISSING}"
echo "[env_v72] profile=$V72_MATRIX_PROFILE"
echo "[env_v72] structure_runner=$STRUCTURE_RUNNER chunk_size=$STRUCTURE_CHUNK_SIZE"
echo "[env_v72] colabfold_batch=$COLABFOLD_BATCH"
if [ -z "${V72_PYTHON:-}" ] || [ ! -x "${V72_PYTHON:-/dev/null}" ]; then
  _v72_echo_err "[env_v72][ERROR] V72_PYTHON is not executable."
fi
return 0 2>/dev/null || exit 0
