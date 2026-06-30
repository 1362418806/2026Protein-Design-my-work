#!/usr/bin/env bash
# Run only the V7 OpenComplex2-style soft validation stage for an existing RUN_ID.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/env_v7.sh"
PYTHON_BIN="${V7_PYTHON:-python}"
RUN_ID="${RUN_ID:-$(cat "$V7_DIR/logs/v7_all_latest.run_id" 2>/dev/null || true)}"
if [ -z "$RUN_ID" ]; then echo "RUN_ID is required" >&2; exit 2; fi
OUT="$V7_OUTPUT_ROOT/$RUN_ID"
mkdir -p "$OUT/04_softcheck"
"$PYTHON_BIN" "$V7_DIR/tools/run_opencomplex2_softcheck_v7.py" \
  --candidates-csv "$OUT/01_stage1/structure_priority_top200_v6.csv" \
  --reference-pdb "${REF_PDB:-$V7_DIR/reference/reference_gfp.pdb}" \
  --out-csv "$OUT/04_softcheck/softcheck_metrics_v7.csv" \
  --work-dir "$OUT/04_softcheck/softcheck_work" \
  --runner "${V7_SOFTCHECK_RUNNER:-metrics-only}" \
  --top-n "${V7_SOFTCHECK_TOP_N:-80}" \
  --timeout-sec "${V7_SOFTCHECK_TIMEOUT_SEC:-21600}" \
  ${V7_SOFTCHECK_PDB_DIR:+--pdb-dir "$V7_SOFTCHECK_PDB_DIR"} \
  ${OPENCOMPLEX2_COMMAND:+--external-command-template "$OPENCOMPLEX2_COMMAND"}
