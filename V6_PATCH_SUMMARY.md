# V6 Patch Summary

Implemented an independent V6 code tree from V5.1 with the following targeted changes.

## Core fixes

1. Closed-loop feedback generation
   - Added `src/synbio_gfp_v6/feedback.py`.
   - Converts previous V5/V5.1/V6 ranked candidate tables into mutation, position, and branch feedback.
   - Writes `feedback/v6_feedback_policy.json`, `mutation_feedback_v6.csv`, `position_feedback_v6.csv`, `branch_feedback_v6.csv`.

2. Feedback-conditioned generator
   - Replaced V5.1 static candidate generator with `src/synbio_gfp_v6/candidates.py`.
   - Supports `feedback_seed_repair`, `feedback_motif_recombine`, `feedback_surface_amplify`, and `feedback_guarded_explore` branches.
   - Keeps baseline conservative branches for safety.

3. Causal objective update
   - Added `v6_objective`, `v6_structure_score`, `v6_feedback_score`, and `v6_blacklist_penalty`.
   - Increased direct influence of structure / ProteinMPNN from V5.1’s weak post-hoc signal to a stronger V6 ranking component.

4. Top-200 structure-priority update
   - Added `feedback_guided_exploitation` bucket.
   - Retains brightness, surface-charge, structure-safe, ProteinMPNN, and diversity reserve buckets.

5. Mechanism-aware Top-6 selector
   - Added mechanism quotas, source-branch cap, Hamming distance, shared-mutation, and site-reuse constraints.
   - Avoids selecting six near-duplicate candidates from one local optimum.

6. Deployment/environment reuse
   - Added `deploy/env_v6.sh`.
   - V6 defaults to `V6_ENV_PREFIX=/hyperai/home/envs/synbio_gfp_v51`, matching the current server state.
   - Reference PDB is copied into `reference/reference_gfp.pdb`.

## Smoke tests run

- `python -m py_compile` on V6 scripts, tools, and package modules.
- `python scripts/run_v6_complete.py` on synthetic data with 200 candidates.
- `python scripts/run_v6_complete.py --previous-ranked-csv ...` verified feedback-conditioned candidate generation.
- `python scripts/run_v6_pipeline.py --dry-run` verified orchestrator command construction.
