# SynBio GFP V7.2 Independent Evidence Matrix

## Purpose

V7.2 is an isolated matrix runner derived from the V7.1 auto/artifact hotfix. It does **not**
write into the V7.1 project tree. Each experiment owns its candidate snapshot, PDB outputs,
ProteinMPNN scores, final Top-6, audit, and summary files.

The matrix extends test coverage while keeping the candidate/structure decision process explicit:

- controlled ColabFold MSA ablations on a fixed `simple_anchor` candidate universe;
- controlled recycle-count ablation on the same candidate universe;
- ESM-35M and ESM-150M candidate-source branches for source-diversity evidence;
- an optional `single_sequence` negative control, excluded from matrix recommendation;
- global-fit pocket/ion-network geometry audit using the beta-barrel transform;
- strict per-submission legality, PDB existence, PDB SHA256, and structure-evidence audit.

## Important interpretation boundary

`matrix_confidence_proxy_v72`, `brightness_safety_proxy_v72`, and
`thermal_robustness_proxy_v72` are **model proxies**, not calibrated wet-lab probabilities.
They do not claim to predict 72C/10min fluorescence retention directly.

V7.2 does not train on its own matrix outputs and does not feed any V7.2 result back into
candidate generation.

## Matrix profiles

- `quick_smoke`: one simple/uniref/recycle-1 experiment. Use this first after deployment.
- `deadline_core` (default): simple uniref/recycle-1, simple uniref_env/recycle-1,
  simple uniref/recycle-3, and an ESM-150M uniref_env/recycle-1 branch.
- `extended`: adds uniref_env/recycle-3, ESM-35M, and the single-sequence negative control.

The controlled simple branches share an identical frozen candidate source. Therefore MSA and
recycle comparisons do not confound candidate generation.

## Deployment

```bash
cd /hyperai/home
unzip -q /path/to/synbio_gfp_v7_2_deployable.zip
cd /hyperai/home/synbio_gfp_v7_2_complete
source deploy/env_v72.sh

"$V72_PYTHON" -m py_compile \
  scripts/run_v72_matrix.py \
  scripts/run_v72_final_from_manifest.py \
  scripts/summarize_v72_matrix.py \
  tools/augment_structure_metrics_v72.py \
  tools/audit_submission_v72.py

test -x "$COLABFOLD_BATCH"
test -f "$REF_PDB"
```

Run one smoke experiment first:

```bash
export V72_MATRIX_PROFILE=quick_smoke
MATRIX_ROOT="$V72_DIR/experiments/v72_matrix_$(date +%Y%m%d_%H%M%S)_smoke"
nohup bash deploy/run_v72_matrix.sh --matrix-root "$MATRIX_ROOT" \
  > "$V72_DIR/logs/v72_smoke.log" 2>&1 &
echo $! > "$V72_DIR/logs/v72_smoke.pid"
```

Confirm that the Stage2 log contains an actual `colabfold_batch` command and does not contain
`metrics-only selected; no new predictions will be generated`.

Then launch the deadline matrix with a fresh root:

```bash
export V72_MATRIX_PROFILE=deadline_core
MATRIX_ROOT="$V72_DIR/experiments/v72_matrix_$(date +%Y%m%d_%H%M%S)_deadline"
nohup bash deploy/run_v72_matrix.sh --matrix-root "$MATRIX_ROOT" \
  > "$V72_DIR/logs/v72_deadline.log" 2>&1 &
echo $! > "$V72_DIR/logs/v72_deadline.pid"
```

Resume safely after an interruption:

```bash
bash deploy/run_v72_matrix.sh --matrix-root "$MATRIX_ROOT" --profile deadline_core
```

Run remaining coverage after the deadline core finishes:

```bash
bash deploy/run_v72_matrix.sh --matrix-root "$MATRIX_ROOT" --profile extended
```

Completed experiments are skipped only when their V7.2 submission and audit JSON both exist and
the audit has `passed=true`.

## Primary outputs

```text
<MATRIX_ROOT>/
├── candidate_sources/<source>/01_stage1/
├── runs/<experiment>/
│   ├── 01_stage1/candidate_universe_v72.csv
│   ├── 02_structure/structure_metrics_v72.csv
│   ├── 03_proteinmpnn/proteinmpnn_scores.csv
│   └── 05_final/
│       ├── submission_v72_for_upload.csv
│       ├── final_top6_v72.csv
│       ├── v72_submission_audit.json
│       └── submission_metrics_v72.json
├── matrix_submission_metrics_v72.csv
├── matrix_candidate_consensus_v72.csv
├── matrix_summary_v72.md
├── matrix_recommendation_v72.json
└── MATRIX_RECOMMENDATION_V72.txt
```

The matrix recommendation is a reporting recommendation only. Any final competition submission
must still be manually chosen from a V7.2 experiment with `hard_ok=true` and `audit_passed=true`.
