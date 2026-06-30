# SynBio GFP V6 Closed-Loop Pipeline

V6 is an independent follow-up to V5/V5.1.  It keeps the validated V5.1 environment and structure/PDB interfaces, but fixes the main algorithmic issue: structure and ProteinMPNN feedback can now influence the next candidate generation round instead of only post-hoc reranking an already fixed Top-200 pool.

## What changed from V5/V5.1

1. `src/synbio_gfp_v6/feedback.py` converts previous ranked candidates, structure metrics, and ProteinMPNN scores into branch-level, mutation-level, and position-level feedback.
2. `src/synbio_gfp_v6/candidates.py` can generate feedback-conditioned candidates using successful motifs, repaired successful parents, surface-amplified variants, and guarded exploration.
3. `scripts/run_v6_complete.py` adds `v6_objective`, `feedback_policy_score`, `motif_blacklist_hits`, `v6_bucket`, and feedback reports.
4. Structure-priority selection now includes a `feedback_guided_exploitation` bucket.
5. Final Top-6 selection enforces mechanism quotas, source-branch caps, Hamming distance, shared-mutation limits, and site-reuse limits.

## Environment

V6 is an independent code tree, but it reuses the validated V5.1 environment by default:

```bash
source /hyperai/home/synbio_gfp_v6_complete/deploy/env_v6.sh
# default: V6_ENV_PREFIX=/hyperai/home/envs/synbio_gfp_v51
```

Use a dedicated V6 environment only if needed:

```bash
export V6_ENV_PREFIX=/hyperai/home/envs/synbio_gfp_v6
bash /hyperai/home/synbio_gfp_v6_complete/deploy/install_v6_env.sh
```

## Recommended closed-loop workflow

### 0. Copy package

```bash
cd /hyperai/home
unzip synbio_gfp_v6_complete.zip
source /hyperai/home/synbio_gfp_v6_complete/deploy/env_v6.sh
```

### 1. Use V5/V5.1 results as feedback seed

Point `V6_PREVIOUS_RANKED_CSV` to the best available V5/V5.1 final ranked table that already contains imported structure and ProteinMPNN columns if possible:

```bash
export V6_PREVIOUS_RANKED_CSV=/hyperai/home/synbio_gfp_v51_complete/outputs/<RUN_ID>/04_final/all_ranked_candidates_v51.csv
```

V5/V5.1 file names also work, for example:

```bash
export V6_PREVIOUS_RANKED_CSV=/hyperai/home/synbio_gfp_v5_complete/outputs/<RUN_ID>/04_final_reuse_pdb_fix_geometry_gate/all_ranked_candidates_v5.csv
```

### 2. Generate V6 feedback-conditioned Top-200 structure tasks

```bash
export RUN_ID=v6_round1_$(date +%Y%m%d_%H%M%S)
export TEAM_NAME=YourTeamName
export FEATURE_MODE=simple
export STRUCTURE_RUNNER=metrics-only
export RUN_PROTEINMPNN=0
export N_CANDIDATES=20000
bash /hyperai/home/synbio_gfp_v6_complete/deploy/run_v6_all.sh
```

Important outputs:

```text
outputs/<RUN_ID>/01_stage1/structure_priority_top200_v6.csv
outputs/<RUN_ID>/01_stage1/af2_af3_structure_tasks.fasta
outputs/<RUN_ID>/01_stage1/feedback/v6_feedback_policy.json
outputs/<RUN_ID>/01_stage1/feedback/mutation_feedback_v6.csv
outputs/<RUN_ID>/01_stage1/feedback/position_feedback_v6.csv
```

### 3. Predict structures and ProteinMPNN scores

Reuse existing V5 style commands, or set the runner directly:

```bash
export REF_PDB=/hyperai/home/synbio_gfp_v6_complete/reference/reference_gfp.pdb
export STRUCTURE_RUNNER=colabfold
export RUN_PROTEINMPNN=1
bash /hyperai/home/synbio_gfp_v6_complete/deploy/run_v6_all.sh
```

Alternatively, reuse PDBs/metrics:

```bash
export V6_STRUCTURE_METRICS_CSV=/path/to/structure_metrics.csv
export V6_PROTEINMPNN_SCORE_CSV=/path/to/proteinmpnn_scores.csv
bash /hyperai/home/synbio_gfp_v6_complete/deploy/run_v6_all.sh
```

### 4. Final selection

Final files are emitted under:

```text
outputs/<RUN_ID>/04_final/final_top6_v6.csv
outputs/<RUN_ID>/04_final/submission_v6.csv
outputs/<RUN_ID>/04_final/v6_pipeline_report.md
```

The final selector still requires real structure metrics by default.  `--allow-proxy-final` exists only for smoke tests and should not be used for serious submissions.

## Core diagnosis fixed by V6

V5/V5.1 did this:

```text
generate -> score -> structure/ProteinMPNN -> rerank existing candidates
```

V6 supports this:

```text
previous measured feedback -> update branch/mutation/site policy -> generate new candidates -> score -> structure priority -> structure/ProteinMPNN -> final portfolio
```

This is the intended causal path for objective validation improvement.
