# SynBio GFP V6.2 Evidence-Safe Pipeline

V6.2 keeps the V6/V6.1 modeling stack but fixes the engineering contract that controls structure evidence, final selection, and feedback reuse.

## Main changes

1. Stage 1 freezes the candidate universe in `01_stage1/candidate_universe_v62.csv` and `candidate_manifest_v62.csv`.
2. Stage 2 writes `all_structure_tasks.fasta` and `resolver_state.json` before any ColabFold call, then reuses exact-sequence PDBs by `sequence_sha1`.
3. Stage 2 predicts missing structures in small MSA chunks and defaults to `--msa-mode mmseqs2_uniref_env`; `single_sequence` is no longer used for final evidence by default.
4. Stage 4 loads the frozen Stage-1 universe and never regenerates candidates.
5. Final Top-6 requires real structure metrics and `structure_pass=True`; underfilled final runs emit `final_underfilled_needs_recovery.csv` instead of proxy candidates.
6. Candidate generation adds `surface_patch`, `loop_rigidity`, and `measured_feedback_repair` sources.
7. Feedback policy defaults to measured-only evidence and no longer uses current-run proxy rows for mutation policy.

## Find reusable V5.1/V6 structures

```bash
cd /hyperai/home/synbio_gfp_v62_complete
source deploy/env_v62.sh

export PDB_SEARCH_ROOTS="/hyperai/home/synbio_gfp_v51_complete/outputs:/hyperai/home/synbio_gfp_v6_complete/outputs:/hyperai/home/synbio_gfp_v61_adaptive/outputs:/hyperai/home/synbio_gfp_v62_complete/outputs"

python tools/build_completed_structure_registry_v62.py \
  --search-roots "$PDB_SEARCH_ROOTS" \
  --out-csv "$V62_DIR/cache/v62/completed_structure_registry_v62.csv" \
  --summary-json "$V62_DIR/cache/v62/completed_structure_registry_v62_summary.json"
```

## Full run

```bash
cd /hyperai/home/synbio_gfp_v62_complete
source deploy/env_v62.sh

export STRUCTURE_POLICY="reuse-first"
export STRUCTURE_RUNNER="metrics-only"   # first audit reuse without new prediction
export PDB_SEARCH_ROOTS="/hyperai/home/synbio_gfp_v51_complete/outputs:/hyperai/home/synbio_gfp_v6_complete/outputs:/hyperai/home/synbio_gfp_v61_adaptive/outputs:/hyperai/home/synbio_gfp_v62_complete/outputs"

bash deploy/run_v62_all.sh
```

If final is underfilled, inspect:

```bash
RUN_ID=$(cat logs/v62_all_latest.run_id)
OUT="$V62_OUTPUT_ROOT/$RUN_ID"
cat "$OUT/02_structure/structure_reuse_report.md"
head "$OUT/04_final/final_underfilled_needs_recovery.csv"
```

Then only predict the highest-priority missing structures by switching:

```bash
export STRUCTURE_RUNNER="colabfold"
export STRUCTURE_CHUNK_SIZE=12
export V62_COLABFOLD_MSA_MODE="mmseqs2_uniref_env"
export V62_COLABFOLD_CHUNK_TIMEOUT_SEC=7200
export V62_COLABFOLD_CHUNK_RETRIES=1
export COLABFOLD_EXTRA_ARGS="--num-models 1 --num-recycle 1"
RUN_ID="$RUN_ID" bash deploy/run_v62_stage2_resolve_structures.sh
RUN_ID="$RUN_ID" bash deploy/run_v62_stage3_proteinmpnn.sh
RUN_ID="$RUN_ID" bash deploy/run_v62_stage4_final.sh
```
