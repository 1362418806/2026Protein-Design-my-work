# SynBio GFP V6.1 Adaptive Pipeline

V6.1 Adaptive keeps the V6 closed-loop generator and adds GPU-profile-aware execution,
automatic exact-sequence PDB reuse, V5.1-style brightness safety teaching, physical thermal
proxy scoring, feedback survival auditing, and stricter final portfolio selection.

## Quick start

```bash
cd /hyperai/home/synbio_gfp_v61_adaptive
source deploy/env_v61.sh

TEAM_NAME=YourTeamName \
DATA_DIR="/input0/2026Protein Design" \
GPU_PROFILE=auto \
STRUCTURE_POLICY=reuse-first \
RUN_PROTEINMPNN=1 \
bash deploy/run_v61_all.sh
```

Monitor:

```bash
tail -f /hyperai/home/synbio_gfp_v61_adaptive/logs/v61_all_latest.log
bash deploy/status_v61_latest.sh
```

## Stage-wise execution

```bash
RUN_ID=v61_manual bash deploy/run_v61_stage0_preflight.sh
RUN_ID=v61_manual bash deploy/run_v61_stage1_generate.sh
RUN_ID=v61_manual bash deploy/run_v61_stage2_resolve_structures.sh
RUN_ID=v61_manual bash deploy/run_v61_stage3_proteinmpnn.sh
RUN_ID=v61_manual bash deploy/run_v61_stage4_final.sh
RUN_ID=v61_manual bash deploy/run_v61_report.sh
```

## Structure policy

- `reuse-only`: only reuse exact-sequence PDBs; no prediction.
- `reuse-first`: default; reuse exact-sequence PDBs and predict missing structures.
- `predict-missing`: production mode when missing PDBs should be filled aggressively.
- `metrics-only`: smoke-test mode only; not recommended for final submission.

PDB reuse key is `sequence_sha1 + sequence_length`; candidate IDs are used only as hints.
`structure_metrics.csv` is recomputed every run.
