#!/usr/bin/env bash
# V6.2 evidence-safe environment.
# Source-safe by design: do not enable or leak `set -e`, `set -u`, or pipefail into interactive shells.
_v62_echo_err() { echo "$*" >&2; }
_v62_project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export V62_DIR="${V62_DIR:-/hyperai/home/synbio_gfp_v62_complete}"
if [ ! -d "$V62_DIR" ]; then export V62_DIR="$_v62_project_dir"; fi
export DATA_DIR="${DATA_DIR:-/input0/2026Protein Design}"
export V62_ENV_PREFIX="${V62_ENV_PREFIX:-${V61_ENV_PREFIX:-/hyperai/home/envs/synbio_gfp_v51}}"
_v62_candidates="${V62_ENV_PREFIX}"
[ -n "${V61_ENV_PREFIX:-}" ] && _v62_candidates="$_v62_candidates ${V61_ENV_PREFIX}"
_v62_candidates="$_v62_candidates /hyperai/home/envs/synbio_gfp_v51 /hyperai/home/envs/synbio_gfp_v61"
_v62_required_modules="pandas openpyxl sklearn yaml"
_v62_first_python=""
export V62_ENV_READY=0
export V62_PYTHON=""
_v62_python_has_deps() {
  "$1" - <<'PYCHECK' >/dev/null 2>&1
import importlib.util
required = ["pandas", "openpyxl", "sklearn", "yaml"]
missing = [name for name in required if importlib.util.find_spec(name) is None]
raise SystemExit(0 if not missing else 1)
PYCHECK
}
for _prefix in $_v62_candidates; do
  _py="$_prefix/bin/python"
  if [ -x "$_py" ]; then
    [ -z "$_v62_first_python" ] && _v62_first_python="$_py"
    if [ "${V62_SKIP_ENV_CHECK:-0}" = "1" ] || _v62_python_has_deps "$_py"; then
      export V62_PYTHON="$_py"
      export V62_ENV_PREFIX="$_prefix"
      export PATH="$_prefix/bin:${PATH:-}"
      export VIRTUAL_ENV="$_prefix"
      export V62_ENV_READY=1
      break
    fi
  fi
done
if [ "$V62_ENV_READY" != "1" ]; then
  if [ -n "$_v62_first_python" ]; then
    export V62_PYTHON="$_v62_first_python"
    export PATH="$(dirname "$_v62_first_python"):${PATH:-}"
  elif command -v python >/dev/null 2>&1; then
    export V62_PYTHON="$(command -v python)"
  fi
fi
hash -r 2>/dev/null || true
export V62_OUTPUT_ROOT="${V62_OUTPUT_ROOT:-$V62_DIR/outputs}"
export V62_CACHE_DIR="${V62_CACHE_DIR:-$V62_DIR/cache/v62}"
export V62_CONFIG="${V62_CONFIG:-$V62_DIR/configs/v62_adaptive.yaml}"
export REF_PDB="${REF_PDB:-$V62_DIR/reference/reference_gfp.pdb}"
export PYTHONPATH="$V62_DIR/src:${PYTHONPATH:-}"
export GPU_PROFILE="${GPU_PROFILE:-auto}"
export FEATURE_MODE="${FEATURE_MODE:-simple}"
export STRUCTURE_POLICY="${STRUCTURE_POLICY:-reuse-first}"
export STRUCTURE_RUNNER="${STRUCTURE_RUNNER:-auto}"
export STRUCTURE_CHUNK_SIZE="${STRUCTURE_CHUNK_SIZE:-12}"
export V62_COLABFOLD_MSA_MODE="${V62_COLABFOLD_MSA_MODE:-mmseqs2_uniref_env}"
export V62_SINGLE_SEQUENCE_DEFAULT="${V62_SINGLE_SEQUENCE_DEFAULT:-0}"
export V62_COLABFOLD_CHUNK_TIMEOUT_SEC="${V62_COLABFOLD_CHUNK_TIMEOUT_SEC:-7200}"
export V62_COLABFOLD_CHUNK_RETRIES="${V62_COLABFOLD_CHUNK_RETRIES:-1}"
export V62_COLABFOLD_RETRY_BACKOFF_SEC="${V62_COLABFOLD_RETRY_BACKOFF_SEC:-60}"
export COLABFOLD_EXTRA_ARGS="${COLABFOLD_EXTRA_ARGS:---num-models 1 --num-recycle 1}"
export RUN_PROTEINMPNN="${RUN_PROTEINMPNN:-1}"
export V61_DIR="$V62_DIR"
export V61_OUTPUT_ROOT="$V62_OUTPUT_ROOT"
export V61_CACHE_DIR="$V62_CACHE_DIR"
export V61_CONFIG="$V62_CONFIG"
mkdir -p "$V62_OUTPUT_ROOT" "$V62_CACHE_DIR" "$V62_DIR/logs"
echo "[env_v62] V62_DIR=$V62_DIR"
echo "[env_v62] DATA_DIR=$DATA_DIR"
echo "[env_v62] V62_ENV_PREFIX=$V62_ENV_PREFIX"
echo "[env_v62] python=${V62_PYTHON:-MISSING}"
echo "[env_v62] V62_ENV_READY=$V62_ENV_READY"
echo "[env_v62] GPU_PROFILE=$GPU_PROFILE FEATURE_MODE=$FEATURE_MODE STRUCTURE_POLICY=$STRUCTURE_POLICY"
echo "[env_v62] STRUCTURE_CHUNK_SIZE=$STRUCTURE_CHUNK_SIZE V62_COLABFOLD_MSA_MODE=$V62_COLABFOLD_MSA_MODE"
echo "[env_v62] COLABFOLD_EXTRA_ARGS=$COLABFOLD_EXTRA_ARGS"
if [ "$V62_ENV_READY" != "1" ]; then
  _v62_echo_err "[env_v62][ERROR] No usable V6.2 Python env found with required modules: $_v62_required_modules"
  _v62_echo_err "[env_v62][ERROR] Tried: $_v62_candidates"
  _v62_echo_err "[env_v62][ERROR] Run: bash deploy/repair_v62_env.sh"
fi
return 0 2>/dev/null || exit 0
