#!/usr/bin/env bash
# Install SynBio GFP V6 runtime into fixed /hyperai/home paths.
# The script avoids system-wide writes; it creates a relocatable conda/micromamba prefix.

set -euo pipefail

HYPERAI_HOME="${HYPERAI_HOME:-/hyperai/home}"
V6_DIR="${V6_DIR:-${HYPERAI_HOME}/synbio_gfp_v6_complete}"
V6_ENV_PREFIX="${V6_ENV_PREFIX:-${HYPERAI_HOME}/envs/synbio_gfp_v6}"
TOOLS_ROOT="${TOOLS_ROOT:-${HYPERAI_HOME}/tools}"
THIRD_PARTY_ROOT="${THIRD_PARTY_ROOT:-${TOOLS_ROOT}/synbio_gfp_v6_third_party}"
PIP_CACHE_DIR="${PIP_CACHE_DIR:-${HYPERAI_HOME}/pip_cache/synbio_gfp_v6}"
INSTALL_LOCALCOLABFOLD="${INSTALL_LOCALCOLABFOLD:-1}"
INSTALL_PROTEINMPNN="${INSTALL_PROTEINMPNN:-1}"
INSTALL_TORCH_MODE="${INSTALL_TORCH_MODE:-cu128}"  # supported: cu128, cpu, skip
INSTALL_JAX_CUDA="${INSTALL_JAX_CUDA:-1}"
PYTHON_VERSION="${PYTHON_VERSION:-3.10}"

log() { printf '[install_v6] %s\n' "$*"; }
fail() { printf '[install_v6][ERROR] %s\n' "$*" >&2; exit 1; }

ensure_writable_dir() {
  local d="$1"
  mkdir -p "$d" 2>/dev/null || true
  if [[ ! -w "$d" ]]; then
    # Try non-interactive sudo only if the server allows it; otherwise fail with a clear diagnosis.
    if command -v sudo >/dev/null 2>&1 && sudo -n true 2>/dev/null; then
      sudo mkdir -p "$d"
      sudo chown -R "$(id -u):$(id -g)" "$d"
      chmod -R u+rwX "$d"
    else
      fail "Directory is not writable and passwordless sudo is unavailable: $d"
    fi
  fi
  chmod -R u+rwX "$d" || true
}

ensure_writable_dir "${HYPERAI_HOME}/envs"
ensure_writable_dir "${TOOLS_ROOT}"
ensure_writable_dir "${THIRD_PARTY_ROOT}"
ensure_writable_dir "${PIP_CACHE_DIR}"

if [[ ! -d "$V6_DIR" ]]; then
  fail "V6_DIR does not exist: $V6_DIR. Unzip/copy synbio_gfp_v6_complete there first."
fi

if [[ -f "${V6_DIR}/deploy/env_v6.sh" ]]; then
  # shellcheck disable=SC1090
  source "${V6_DIR}/deploy/env_v6.sh"
fi

create_env_with_conda() {
  local conda_bin="$1"
  if [[ ! -x "${V6_ENV_PREFIX}/bin/python" ]]; then
    log "Creating conda env at ${V6_ENV_PREFIX}"
    "$conda_bin" create -y -p "$V6_ENV_PREFIX" "python=${PYTHON_VERSION}" pip
  else
    log "Using existing env at ${V6_ENV_PREFIX}"
  fi
}

create_env_with_micromamba() {
  local mm_root="${TOOLS_ROOT}/micromamba"
  local mm_bin="${mm_root}/bin/micromamba"
  if [[ ! -x "$mm_bin" ]]; then
    log "Installing micromamba under ${mm_root}"
    mkdir -p "${mm_root}/bin"
    if command -v curl >/dev/null 2>&1; then
      curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj -C "$mm_root" bin/micromamba
    elif command -v wget >/dev/null 2>&1; then
      wget -qO- https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj -C "$mm_root" bin/micromamba
    else
      fail "Neither curl nor wget is available to install micromamba."
    fi
  fi
  if [[ ! -x "${V6_ENV_PREFIX}/bin/python" ]]; then
    log "Creating micromamba env at ${V6_ENV_PREFIX}"
    "$mm_bin" create -y -p "$V6_ENV_PREFIX" -c conda-forge "python=${PYTHON_VERSION}" pip
  else
    log "Using existing env at ${V6_ENV_PREFIX}"
  fi
}

if command -v conda >/dev/null 2>&1; then
  create_env_with_conda "$(command -v conda)"
elif [[ -n "${CONDA_EXE:-}" && -x "${CONDA_EXE}" ]]; then
  create_env_with_conda "${CONDA_EXE}"
else
  create_env_with_micromamba
fi

PY="${V6_ENV_PREFIX}/bin/python"
PIP="${V6_ENV_PREFIX}/bin/pip"
[[ -x "$PY" ]] || fail "Python was not created at ${PY}"

log "Installing Python requirements"
"$PY" -m pip install --upgrade pip setuptools wheel
"$PIP" install --cache-dir "$PIP_CACHE_DIR" -r "${V6_DIR}/requirements-v6.txt"

if [[ "$INSTALL_TORCH_MODE" == "cu128" ]]; then
  log "Installing PyTorch CUDA 12.8 wheels"
  "$PIP" install --cache-dir "$PIP_CACHE_DIR" torch --index-url https://download.pytorch.org/whl/cu128 || {
    log "CUDA PyTorch install failed; falling back to CPU torch for pipeline smoke tests."
    "$PIP" install --cache-dir "$PIP_CACHE_DIR" torch --index-url https://download.pytorch.org/whl/cpu
  }
elif [[ "$INSTALL_TORCH_MODE" == "cpu" ]]; then
  log "Installing CPU PyTorch wheels"
  "$PIP" install --cache-dir "$PIP_CACHE_DIR" torch --index-url https://download.pytorch.org/whl/cpu
else
  log "Skipping PyTorch installation by request."
fi

if [[ "$INSTALL_JAX_CUDA" == "1" ]]; then
  log "Installing JAX CUDA12 wheels into the V6 env"
  "$PIP" install --cache-dir "$PIP_CACHE_DIR" --upgrade "jax[cuda12]" || log "JAX CUDA install failed; ColabFold env may still provide its own JAX."
fi

if [[ "$INSTALL_PROTEINMPNN" == "1" ]]; then
  if [[ ! -d "${TOOLS_ROOT}/ProteinMPNN/.git" ]]; then
    log "Cloning ProteinMPNN into ${TOOLS_ROOT}/ProteinMPNN"
    git clone https://github.com/dauparas/ProteinMPNN.git "${TOOLS_ROOT}/ProteinMPNN"
  else
    log "ProteinMPNN already exists at ${TOOLS_ROOT}/ProteinMPNN"
  fi
fi

if [[ "$INSTALL_LOCALCOLABFOLD" == "1" ]]; then
  if [[ ! -x "${TOOLS_ROOT}/localcolabfold/.pixi/envs/default/bin/colabfold_batch" && ! -x "${TOOLS_ROOT}/localcolabfold/colabfold-conda/bin/colabfold_batch" ]]; then
    log "Installing LocalColabFold into ${TOOLS_ROOT}/localcolabfold"
    rm -rf "${TOOLS_ROOT}/localcolabfold"
    git clone https://github.com/YoshitakaMo/localcolabfold.git "${TOOLS_ROOT}/localcolabfold"
    (cd "${TOOLS_ROOT}/localcolabfold" && bash install_colabbatch_linux.sh)
  else
    log "LocalColabFold already appears to be installed."
  fi
fi

log "Writing activation helper"
cat > "${V6_ENV_PREFIX}/activate_v6.sh" <<ACT
#!/usr/bin/env bash
source "${V6_DIR}/deploy/env_v6.sh"
export PATH="${V6_ENV_PREFIX}/bin:\$PATH"
ACT
chmod +x "${V6_ENV_PREFIX}/activate_v6.sh"

log "Running lightweight import verification"
"$PY" - <<'PY'
import sys
import numpy as np
import pandas as pd
import sklearn
import joblib
import yaml
import scipy
import openpyxl
import Bio
print('Python:', sys.version)
print('numpy:', np.__version__)
print('pandas:', pd.__version__)
print('sklearn:', sklearn.__version__)
print('joblib:', joblib.__version__)
print('scipy:', scipy.__version__)
print('openpyxl:', openpyxl.__version__)
print('Bio import OK')
print('yaml import OK')
try:
    import esm
    print('esm import OK')
except Exception as exc:
    raise RuntimeError(f'fair-esm import failed: {exc}')
PY

# Keep the old helper name as a compatibility shim for operators who still source activate_v5.sh.
cat > "${V6_ENV_PREFIX}/activate_v5.sh" <<ACT
#!/usr/bin/env bash
source "${V6_ENV_PREFIX}/activate_v6.sh"
ACT
chmod +x "${V6_ENV_PREFIX}/activate_v5.sh"

log "Done. Activate with: source ${V6_ENV_PREFIX}/activate_v6.sh"
