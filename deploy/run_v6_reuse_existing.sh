#!/usr/bin/env bash
# Run V6 with existing ColabFold PDBs and optional existing ProteinMPNN scores.
# This is the recommended handoff path after a V5 run has already produced PDBs.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/env_v6.sh"

: "${V6_REUSE_PDB_DIR:?Set V6_REUSE_PDB_DIR to the existing ColabFold/AF PDB directory.}"
: "${REF_PDB:?Set REF_PDB or V6_REFERENCE_PDB to the GFP reference PDB.}"

export STRUCTURE_RUNNER="${STRUCTURE_RUNNER:-skip}"
export RUN_PROTEINMPNN="${RUN_PROTEINMPNN:-0}"
export V6_BACKGROUND="${V6_BACKGROUND:-1}"

bash "${SCRIPT_DIR}/run_v6_all.sh"
