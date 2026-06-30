# V6.2 ColabFold Structure Hotfix Summary

## Scope

This hotfix targets the V6.2 structure-generation path only. Candidate generation, final portfolio logic, and ProteinMPNN scoring are not algorithmically changed.

## Key changes

1. ColabFold now defaults to MSA evidence instead of `single_sequence`.
   - Default mode: `--msa-mode mmseqs2_uniref_env`.
   - Default extra args: `--num-models 1 --num-recycle 1`.
   - `single_sequence` remains available only when explicitly requested via `V62_SINGLE_SEQUENCE_DEFAULT=1` or explicit `COLABFOLD_EXTRA_ARGS`.

2. Missing structures are predicted in smaller chunks.
   - Default `STRUCTURE_CHUNK_SIZE=12`, matching the successful MSA pilot behavior.

3. ColabFold chunks are bounded and retryable.
   - `V62_COLABFOLD_CHUNK_TIMEOUT_SEC=7200` by default.
   - `V62_COLABFOLD_CHUNK_RETRIES=1` by default.
   - Timed-out backend process groups are terminated to avoid orphaned MMseqs/ColabFold workers.

4. Stage2 wrappers now expose ColabFold controls.
   - `COLABFOLD_BATCH` can explicitly point to the correct binary.
   - `V62_COLABFOLD_MSA_MODE` controls the default MSA mode.
   - `V62_COLABFOLD_RETRY_BACKOFF_SEC` controls retry backoff.

5. A Stage2 chunk artifact bug was fixed.
   - The chunk directory is created before writing chunk CSV/FASTA recovery artifacts.

## Recommended Stage2 recovery command

```bash
cd /hyperai/home/synbio_gfp_v62_complete
source deploy/env_v62.sh

export STRUCTURE_RUNNER="colabfold"
export STRUCTURE_POLICY="reuse-first"
export STRUCTURE_CHUNK_SIZE=12
export V62_COLABFOLD_MSA_MODE="mmseqs2_uniref_env"
export V62_COLABFOLD_CHUNK_TIMEOUT_SEC=7200
export V62_COLABFOLD_CHUNK_RETRIES=1
export COLABFOLD_EXTRA_ARGS="--num-models 1 --num-recycle 1"
# export COLABFOLD_BATCH="/path/to/colabfold_batch"  # Optional explicit binary.

RUN_ID="<existing_or_new_run_id>" bash deploy/run_v62_stage2_resolve_structures.sh
```

## Final evidence policy

Do not relax `structure_pass`. Do not use `single_sequence` outputs as final evidence unless explicitly running a diagnostic experiment outside the final submission path.

## v0.23 static-review adjustment

A final-selector evidence hardening fix was added after static review:

- `scripts/run_v62_final_from_manifest.py` now treats `has_real_structure=True` only when the merged structure row has `structure_status == "rmsd_checked"`, complete reference-based RMSD metrics, and imported real-structure evidence.
- This prevents confidence-proxy or malformed metrics rows from being upgraded to final real-structure evidence during sequence-SHA1 merging.
