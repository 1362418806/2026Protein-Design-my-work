# SynBio GFP V7 Orthogonal-Validation Pipeline

V7 keeps the V6.2 evidence-safe contract and adds optional orthogonal soft validation, designed for OpenComplex2 (OC2) or any external AF-family/all-atom predictor that can emit PDB files.

## Core contract

- Stage 1 freezes the candidate universe (`candidate_universe_v62.csv`).
- Stage 2 resolves real structure evidence by `sequence_sha1` using reuse-first ColabFold/PDB metrics.
- Stage 3 imports or runs ProteinMPNN score-only compatibility checks.
- Stage 4 optionally runs OpenComplex2-style soft validation on a limited top-N pool.
- Stage 5 final selection uses only frozen candidates and requires primary real structure + `structure_pass=True`.
- Stage 6 runs final legality/evidence audit, including Exclusion List exact-match checking.

OC2 soft validation is never allowed to replace the primary ColabFold/PDB structure gate. It is an orthogonal confidence signal and a small ranking adjustment only.

## Main entrypoint

```bash
cd /hyperai/home/synbio_gfp_v7_complete
source deploy/env_v7.sh
export TEAM_NAME="YourTeamName"
export COLABFOLD_BATCH="/hyperai/home/tools/localcolabfold/.pixi/envs/default/bin/colabfold_batch"
bash deploy/run_v7_all.sh
```

## Enable OpenComplex2 soft validation

### Option A: run an external OC2 command template

The command must accept a FASTA and output directory. Use `{fasta}` and `{out_dir}` placeholders.

```bash
source deploy/env_v7.sh
export RUN_V7_SOFTCHECK=1
export V7_SOFTCHECK_RUNNER=external
export V7_SOFTCHECK_TOP_N=80
export OPENCOMPLEX2_COMMAND='python /path/to/opencomplex2_infer.py --input_fasta {fasta} --output_dir {out_dir}'
bash deploy/run_v7_all.sh
```

### Option B: score an existing OC2 PDB directory

```bash
source deploy/env_v7.sh
export RUN_V7_SOFTCHECK=1
export V7_SOFTCHECK_RUNNER=metrics-only
export V7_SOFTCHECK_PDB_DIR=/path/to/opencomplex2_pdbs
bash deploy/run_v7_all.sh
```

### Option C: import an existing softcheck CSV

```bash
source deploy/env_v7.sh
export V7_SOFTCHECK_CSV=/path/to/softcheck_metrics_v7.csv
bash deploy/run_v7_all.sh
```

The imported CSV should contain `sequence_sha1` or `sequence`, plus any of:

- `oc2_barrel_rmsd` / `barrel_rmsd`
- `oc2_pocket_rmsd` / `pocket_rmsd`
- `oc2_ion_network_rmsd` / `ion_network_rmsd`
- `oc2_mean_plddt` / `mean_plddt`
- `oc2_softcheck_score` / `softcheck_score`
- `oc2_pass` / `structure_pass`

## New V7 artifacts

- `04_softcheck/softcheck_metrics_v7.csv`
- `04_softcheck/softcheck_work/softcheck_summary_v7.json`
- `05_final/all_ranked_candidates_v7.csv`
- `05_final/final_top6_v7.csv`
- `05_final/submission_v7_for_upload.csv`
- `05_final/v7_final_summary.json`
- `05_final/v7_submission_audit.json`
- `05_final/v7_submission_audit.md`

Backward-compatible V6.2 file names are also written where practical.

## Interpretation

Use V7 when you do not have wet-lab feedback and want to reduce surrogate bias through model disagreement checks. A candidate that passes ColabFold/PDB evidence but fails OC2 should not be automatically discarded, but it should be treated as lower confidence. A candidate that passes both has stronger orthogonal structural support.

## V7.1 biology soft-risk layer

V7.1 adds a no-wet-lab-feedback risk layer before final portfolio selection. It does not replace the primary structure gate. The new layer estimates whether a candidate has hidden expression, solubility, aggregation, or thermal-retention risk that may not be captured by static RMSD alone.

New final columns include:

```text
v7_biology_validity_score
v7_thermal_validity_score
v7_solubility_score
v7_sequence_aggregation_risk
v7_sequence_fragility_score
v7_neighborhood_support
v7_soft_risk_penalty
v7_surrogate_confidence_score
v7_objective
```

The final selector still requires real structure evidence and `structure_pass=True`. OpenComplex2 or another external structure model is used only as a soft consensus signal. The final audit writes `v7_submission_audit.json` and `v7_submission_audit.md` and can be configured to hard-fail biology soft-risk warnings if desired.

Standalone final reranking now runs:

```bash
RUN_ID=<your_run_id> bash deploy/run_v7_stage5_final.sh
```

This produces both the final CSV files and the final audit report.
