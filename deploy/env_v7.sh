#!/usr/bin/env bash
# V7 orthogonal-validation environment.
# Source-safe by design: do not enable or leak `set -e`, `set -u`, or pipefail into interactive shells.
_v7_echo_err() { echo "$*" >&2; }
_v7_project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export V7_DIR="${V7_DIR:-/hyperai/home/synbio_gfp_v7_complete}"
if [ ! -d "$V7_DIR" ]; then export V7_DIR="$_v7_project_dir"; fi
export DATA_DIR="${DATA_DIR:-/input0/2026Protein Design}"
export V7_ENV_PREFIX="${V7_ENV_PREFIX:-${V62_ENV_PREFIX:-${V61_ENV_PREFIX:-/hyperai/home/envs/synbio_gfp_v51}}}"
_v7_candidates="${V7_ENV_PREFIX}"
[ -n "${V62_ENV_PREFIX:-}" ] && _v7_candidates="$_v7_candidates ${V62_ENV_PREFIX}"
[ -n "${V61_ENV_PREFIX:-}" ] && _v7_candidates="$_v7_candidates ${V61_ENV_PREFIX}"
_v7_candidates="$_v7_candidates /hyperai/home/envs/synbio_gfp_v51 /hyperai/home/envs/synbio_gfp_v61"
_v7_required_modules="pandas openpyxl sklearn yaml"
_v7_first_python=""
export V7_ENV_READY=0
export V7_PYTHON=""
_v7_python_has_deps() {
  "$1" - <<'PYCHECK' >/dev/null 2>&1
import importlib.util
required = ["pandas", "openpyxl", "sklearn", "yaml"]
missing = [name for name in required if importlib.util.find_spec(name) is None]
raise SystemExit(0 if not missing else 1)
PYCHECK
}
for _prefix in $_v7_candidates; do
  _py="$_prefix/bin/python"
  if [ -x "$_py" ]; then
    [ -z "$_v7_first_python" ] && _v7_first_python="$_py"
    if [ "${V7_SKIP_ENV_CHECK:-0}" = "1" ] || _v7_python_has_deps "$_py"; then
      export V7_PYTHON="$_py"
      export V7_ENV_PREFIX="$_prefix"
      export PATH="$_prefix/bin:${PATH:-}"
      export VIRTUAL_ENV="$_prefix"
      export V7_ENV_READY=1
      break
    fi
  fi
done
if [ "$V7_ENV_READY" != "1" ]; then
  if [ -n "$_v7_first_python" ]; then
    export V7_PYTHON="$_v7_first_python"
    export PATH="$(dirname "$_v7_first_python"):${PATH:-}"
  elif command -v python >/dev/null 2>&1; then
    export V7_PYTHON="$(command -v python)"
  fi
fi
hash -r 2>/dev/null || true
export V7_OUTPUT_ROOT="${V7_OUTPUT_ROOT:-$V7_DIR/outputs}"
export V7_CACHE_DIR="${V7_CACHE_DIR:-$V7_DIR/cache/v7}"
export V7_CONFIG="${V7_CONFIG:-$V7_DIR/configs/v7_adaptive.yaml}"
export REF_PDB="${REF_PDB:-$V7_DIR/reference/reference_gfp.pdb}"
export PYTHONPATH="$V7_DIR/src:${PYTHONPATH:-}"
export GPU_PROFILE="${GPU_PROFILE:-auto}"
export FEATURE_MODE="${FEATURE_MODE:-simple}"
export STRUCTURE_POLICY="${STRUCTURE_POLICY:-reuse-first}"
export STRUCTURE_RUNNER="${STRUCTURE_RUNNER:-auto}"
export STRUCTURE_CHUNK_SIZE="${STRUCTURE_CHUNK_SIZE:-12}"
export V7_COLABFOLD_MSA_MODE="${V7_COLABFOLD_MSA_MODE:-mmseqs2_uniref_env}"
export V7_COLABFOLD_CHUNK_TIMEOUT_SEC="${V7_COLABFOLD_CHUNK_TIMEOUT_SEC:-7200}"
export V7_COLABFOLD_CHUNK_RETRIES="${V7_COLABFOLD_CHUNK_RETRIES:-1}"
export V7_COLABFOLD_RETRY_BACKOFF_SEC="${V7_COLABFOLD_RETRY_BACKOFF_SEC:-60}"
export COLABFOLD_EXTRA_ARGS="${COLABFOLD_EXTRA_ARGS:---num-models 1 --num-recycle 1}"
export RUN_PROTEINMPNN="${RUN_PROTEINMPNN:-1}"
export RUN_V7_SOFTCHECK="${RUN_V7_SOFTCHECK:-0}"
export V7_SOFTCHECK_RUNNER="${V7_SOFTCHECK_RUNNER:-metrics-only}"
export V7_SOFTCHECK_TOP_N="${V7_SOFTCHECK_TOP_N:-80}"
export V7_SOFTCHECK_TIMEOUT_SEC="${V7_SOFTCHECK_TIMEOUT_SEC:-21600}"
# Backward-compatible aliases for V6.2/V6.1 deploy snippets.
export V62_DIR="$V7_DIR"
export V62_OUTPUT_ROOT="$V7_OUTPUT_ROOT"
export V62_CACHE_DIR="$V7_CACHE_DIR"
export V62_CONFIG="$V7_CONFIG"
export V62_PYTHON="$V7_PYTHON"
export V62_ENV_READY="$V7_ENV_READY"
export V61_DIR="$V7_DIR"
export V61_OUTPUT_ROOT="$V7_OUTPUT_ROOT"
export V61_CACHE_DIR="$V7_CACHE_DIR"
export V61_CONFIG="$V7_CONFIG"
mkdir -p "$V7_OUTPUT_ROOT" "$V7_CACHE_DIR" "$V7_DIR/logs"
echo "[env_v7] V7_DIR=$V7_DIR"
echo "[env_v7] DATA_DIR=$DATA_DIR"
echo "[env_v7] V7_ENV_PREFIX=$V7_ENV_PREFIX"
echo "[env_v7] python=${V7_PYTHON:-MISSING}"
echo "[env_v7] V7_ENV_READY=$V7_ENV_READY"
echo "[env_v7] GPU_PROFILE=$GPU_PROFILE FEATURE_MODE=$FEATURE_MODE STRUCTURE_POLICY=$STRUCTURE_POLICY"
echo "[env_v7] STRUCTURE_CHUNK_SIZE=$STRUCTURE_CHUNK_SIZE V7_COLABFOLD_MSA_MODE=$V7_COLABFOLD_MSA_MODE"
echo "[env_v7] RUN_V7_SOFTCHECK=$RUN_V7_SOFTCHECK V7_SOFTCHECK_RUNNER=$V7_SOFTCHECK_RUNNER"
if [ "$V7_ENV_READY" != "1" ]; then
  _v7_echo_err "[env_v7][ERROR] No usable V7 Python env found with required modules: $_v7_required_modules"
  _v7_echo_err "[env_v7][ERROR] Tried: $_v7_candidates"
fi
return 0 2>/dev/null || exit 0
