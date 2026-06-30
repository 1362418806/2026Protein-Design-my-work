# V7 Patch Summary

## Added

1. `tools/run_opencomplex2_softcheck_v7.py`
   - Imports existing OC2 metrics, scores existing OC2 PDBs, or invokes an external command template.
   - Computes sfGFP barrel/pocket/ion-network RMSD against the reference PDB.
   - Emits `oc2_softcheck_score` and `oc2_pass` as soft validation evidence.

2. `scripts/run_v7_final_from_manifest.py`
   - Keeps frozen-manifest final selection.
   - Adds Exclusion List exact-match filtering at final-selection time.
   - Merges optional OC2/OpenComplex2 softcheck evidence by `sequence_sha1`.
   - Adds `v7_confidence_score`, `v7_orthogonal_structure_score`, and `v7_objective`.

3. `tools/audit_submission_v7.py`
   - Checks final submission legality, Exclusion List hits, primary structure evidence, protected/pocket/ion-network mutation counts, low-brightness risk policy, and mutation-burden policy.

4. `scripts/run_v7_pipeline.py` and deploy wrappers
   - Adds optional Stage 4 softcheck before final selection.
   - Adds mandatory final audit stage.

## Design principle

OpenComplex2 is a soft validator only. It improves confidence and exposes model disagreement, but it must never replace the primary real-structure gate.
