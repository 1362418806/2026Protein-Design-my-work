#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/env_v7.sh"
PYTHON_BIN="${V7_PYTHON:-python}"
RUN_ID="${RUN_ID:-$(cat "$V7_DIR/logs/v7_all_latest.run_id" 2>/dev/null || true)}"
if [ -z "$RUN_ID" ]; then echo "RUN_ID is required" >&2; exit 2; fi
OUT="$V7_OUTPUT_ROOT/$RUN_ID"
mkdir -p "$OUT/05_final"
FROZEN="$OUT/01_stage1/candidate_universe_v62.csv"
[ -s "$FROZEN" ] || FROZEN="$OUT/01_stage1/all_ranked_candidates_v6.csv"
cmd=("$PYTHON_BIN" "$V7_DIR/scripts/run_v7_final_from_manifest.py"
  --config "${V7_CONFIG:-$V7_DIR/configs/v7_adaptive.yaml}"
  --frozen-candidates-csv "$FROZEN"
  --structure-metrics-csv "$OUT/02_structure/structure_metrics.csv"
  --out-dir "$OUT/05_final"
  --team-name "${TEAM_NAME:-YourTeamName}"
  --data-dir "$DATA_DIR")
[ -s "$OUT/03_proteinmpnn/proteinmpnn_scores.csv" ] && cmd+=(--proteinmpnn-score-csv "$OUT/03_proteinmpnn/proteinmpnn_scores.csv")
[ -s "$OUT/04_softcheck/softcheck_metrics_v7.csv" ] && cmd+=(--softcheck-csv "$OUT/04_softcheck/softcheck_metrics_v7.csv")
"${cmd[@]}"

# Final audit is part of Stage 5 so standalone reranking cannot bypass legality/evidence checks.
audit_cmd=("$PYTHON_BIN" "$V7_DIR/tools/audit_submission_v7.py"
  --submission-csv "$OUT/05_final/submission_v7_for_upload.csv"
  --final-csv "$OUT/05_final/final_top6_v7.csv"
  --data-dir "$DATA_DIR"
  --config "${V7_CONFIG:-$V7_DIR/configs/v7_adaptive.yaml}"
  --out-json "$OUT/05_final/v7_submission_audit.json"
  --out-md "$OUT/05_final/v7_submission_audit.md")
if [ "${V7_AUDIT_FAIL_ON_ERROR:-1}" = "1" ]; then
  audit_cmd+=(--fail-on-error)
fi
"${audit_cmd[@]}"
