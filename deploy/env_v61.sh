#!/usr/bin/env bash
# SynBio GFP V6.1 adaptive runtime environment.
# Source this file before running one-click or stage-wise scripts.
set -euo pipefail

export V61_DIR="${V61_DIR:-/hyperai/home/synbio_gfp_v61_adaptive}"
if [ ! -d "$V61_DIR" ]; then
  # Allow sandbox/local testing from the checked-out package directory.
  export V61_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fi

# Dataset path priority: requested input0 path first, then historical HyperAI fallbacks.
if [ -z "${DATA_DIR:-}" ]; then
  if [ -d "/input0/2026Protein Design" ]; then
    export DATA_DIR="/input0/2026Protein Design"
  elif [ -d "/hyperai/input/input0/2026Protein Design" ]; then
    export DATA_DIR="/hyperai/input/input0/2026Protein Design"
  elif [ -d "/root/2026Protein Design" ]; then
    export DATA_DIR="/root/2026Protein Design"
  else
    export DATA_DIR="/input0/2026Protein Design"
  fi
fi

export V61_OUTPUT_ROOT="${V61_OUTPUT_ROOT:-$V61_DIR/outputs}"
export V61_CACHE_DIR="${V61_CACHE_DIR:-$V61_DIR/cache/v61}"
export REF_PDB="${REF_PDB:-$V61_DIR/reference/reference_gfp.pdb}"
export PYTHONPATH="$V61_DIR/src:${PYTHONPATH:-}"
export PATH="/hyperai/home/tools/localcolabfold/.pixi/envs/default/bin:${PATH}"

# Compute profile is adaptive and user-overridable. No algorithm path is hardcoded to 3090.
export GPU_PROFILE="${GPU_PROFILE:-auto}"
case "$GPU_PROFILE" in
  rtx3090_24g)
    export FEATURE_MODE="${FEATURE_MODE:-simple}"
    export ESM_BATCH_SIZE="${ESM_BATCH_SIZE:-2}"
    export STRUCTURE_TOP_N="${STRUCTURE_TOP_N:-160}"
    export COLABFOLD_EXTRA_ARGS="${COLABFOLD_EXTRA_ARGS:---num-models 1 --num-recycle 1}"
    ;;
  rtx5090)
    export FEATURE_MODE="${FEATURE_MODE:-esm}"
    export ESM_BATCH_SIZE="${ESM_BATCH_SIZE:-4}"
    export STRUCTURE_TOP_N="${STRUCTURE_TOP_N:-240}"
    export COLABFOLD_EXTRA_ARGS="${COLABFOLD_EXTRA_ARGS:---num-models 1 --num-recycle 3}"
    ;;
  a800)
    export FEATURE_MODE="${FEATURE_MODE:-esm}"
    export ESM_BATCH_SIZE="${ESM_BATCH_SIZE:-8}"
    export STRUCTURE_TOP_N="${STRUCTURE_TOP_N:-300}"
    export COLABFOLD_EXTRA_ARGS="${COLABFOLD_EXTRA_ARGS:---num-models 1 --num-recycle 3}"
    ;;
  auto|*)
    export FEATURE_MODE="${FEATURE_MODE:-simple}"
    export ESM_BATCH_SIZE="${ESM_BATCH_SIZE:-2}"
    export STRUCTURE_TOP_N="${STRUCTURE_TOP_N:-200}"
    export COLABFOLD_EXTRA_ARGS="${COLABFOLD_EXTRA_ARGS:---num-models 1 --num-recycle 1}"
    ;;
esac

export STRUCTURE_POLICY="${STRUCTURE_POLICY:-reuse-first}"
export RUN_PROTEINMPNN="${RUN_PROTEINMPNN:-1}"
export PDB_SEARCH_ROOTS="${PDB_SEARCH_ROOTS:-/hyperai/home/synbio_gfp_v51_complete/outputs:/hyperai/home/synbio_gfp_v6_complete/outputs:/hyperai/home/synbio_gfp_v5_complete/outputs:/hyperai/home/synbio_gfp_v4_complete/outputs}"
mkdir -p "$V61_OUTPUT_ROOT" "$V61_CACHE_DIR" "$V61_DIR/logs"

echo "[env_v61] V61_DIR=$V61_DIR"
echo "[env_v61] DATA_DIR=$DATA_DIR"
echo "[env_v61] GPU_PROFILE=$GPU_PROFILE FEATURE_MODE=$FEATURE_MODE STRUCTURE_POLICY=$STRUCTURE_POLICY"
