#!/usr/bin/env bash
# Shared environment variables for SynBio GFP V6.
# V6 is an independent code tree but reuses the existing V5.1 Python environment by default.

set -euo pipefail

export HYPERAI_HOME="${HYPERAI_HOME:-/hyperai/home}"
export V6_DIR="${V6_DIR:-${HYPERAI_HOME}/synbio_gfp_v6_complete}"
# Reuse the validated V5.1 environment unless a dedicated V6 env is explicitly provided.
export V6_ENV_PREFIX="${V6_ENV_PREFIX:-${HYPERAI_HOME}/envs/synbio_gfp_v51}"
export TOOLS_ROOT="${TOOLS_ROOT:-${HYPERAI_HOME}/tools}"
export THIRD_PARTY_ROOT="${THIRD_PARTY_ROOT:-${TOOLS_ROOT}/synbio_gfp_v6_third_party}"
export PIP_CACHE_DIR="${PIP_CACHE_DIR:-${HYPERAI_HOME}/pip_cache/synbio_gfp_v6}"

export DATA_DIR="${DATA_DIR:-/hyperai/input/input0/2026Protein Design}"
export V6_CACHE_DIR="${V6_CACHE_DIR:-${V6_DIR}/cache/v6}"
export V6_OUTPUT_ROOT="${V6_OUTPUT_ROOT:-${V6_DIR}/outputs}"
export V6_REFERENCE_DIR="${V6_REFERENCE_DIR:-${V6_DIR}/reference}"

export LOCALCOLABFOLD_DIR="${LOCALCOLABFOLD_DIR:-${TOOLS_ROOT}/localcolabfold}"
export LOCALCOLABFOLD_BIN="${LOCALCOLABFOLD_BIN:-${LOCALCOLABFOLD_DIR}/.pixi/envs/default/bin}"
export PROTEINMPNN_DIR="${PROTEINMPNN_DIR:-${TOOLS_ROOT}/ProteinMPNN}"
export OPENCOMPLEX_DIR="${OPENCOMPLEX_DIR:-${TOOLS_ROOT}/OpenComplex}"

export V6_FIXED_LOG_DIR="${V6_FIXED_LOG_DIR:-${V6_DIR}/logs}"
export V6_LATEST_LOG="${V6_LATEST_LOG:-${V6_FIXED_LOG_DIR}/v6_all_latest.log}"
export V6_LATEST_PID="${V6_LATEST_PID:-${V6_FIXED_LOG_DIR}/v6_all_latest.pid}"
export V6_LATEST_ENV="${V6_LATEST_ENV:-${V6_FIXED_LOG_DIR}/v6_all_latest.env}"
export V6_LATEST_RUN_ID_FILE="${V6_LATEST_RUN_ID_FILE:-${V6_FIXED_LOG_DIR}/v6_all_latest.run_id}"

# Optional closed-loop feedback inputs. Point these to V5/V5.1 final all_ranked_candidates CSVs.
export V6_PREVIOUS_RANKED_CSV="${V6_PREVIOUS_RANKED_CSV:-}"
export V6_FEEDBACK_CSV="${V6_FEEDBACK_CSV:-}"

export PYTHONPATH="${V6_DIR}/src:${PYTHONPATH:-}"
export PATH="${V6_ENV_PREFIX}/bin:${LOCALCOLABFOLD_BIN}:${PATH}"
export V6_PYTHON="${V6_PYTHON:-${V6_ENV_PREFIX}/bin/python}"
export V6_PIP="${V6_PIP:-${V6_ENV_PREFIX}/bin/pip}"
export PYTHONNOUSERSITE="${PYTHONNOUSERSITE:-1}"

mkdir -p "${V6_CACHE_DIR}" "${V6_OUTPUT_ROOT}" "${V6_REFERENCE_DIR}" "${THIRD_PARTY_ROOT}" "${PIP_CACHE_DIR}" "${V6_FIXED_LOG_DIR}"
