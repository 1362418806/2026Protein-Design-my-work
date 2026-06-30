# V7.1 Deployable Algorithm Patch Summary

This patch keeps the V6.2/V7 evidence-safe contract intact and adds extra surrogate-validity checks for the no-wet-lab-feedback setting.

## What changed

1. Added V7 biology soft-risk features in `scripts/run_v7_final_from_manifest.py`:
   - approximate pH-7 net charge and charge drift from the sfGFP-like baseline
   - hydrophobic/aromatic/cysteine/proline-glycine composition checks
   - longest hydrophobic run and longest repeat run penalties
   - mutation-level fragile-alt, charge-flip, hydrophobic-alt and aromatic/cysteine-alt fractions
   - sequence aggregation risk, solubility score, sequence fragility score
   - mutation-neighborhood support and neighborhood fragility
   - thermal-validity and biology-validity scores

2. Integrated these features into V7 ranking:
   - `v7_confidence_score`
   - `v7_surrogate_confidence_score`
   - `v7_thermal_validity_score`
   - `v7_biology_validity_score`
   - `v7_soft_risk_penalty`
   - `v7_objective`

3. Kept hard gates unchanged:
   - real ColabFold/PDB structure evidence is still required
   - `structure_pass=True` is still required
   - OpenComplex2/OC2-style softcheck remains a soft validator only
   - confidence/proxy structures cannot become final hard evidence

4. Extended final audit:
   - final audit now reports biology-soft-risk warnings
   - optional hard fail can be enabled through `v7.final_audit.hard_fail_biology_soft_risk`

5. Made standalone Stage 5 safer:
   - `deploy/run_v7_stage5_final.sh` now runs both final selection and submission audit.

## New outputs

- `05_final/v7_soft_risk_top500.csv`
- `05_final/v7_final_soft_risk_report.csv`
- `05_final/v7_submission_audit.json`
- `05_final/v7_submission_audit.md`
- `05_final/v7_final_summary.json`

## Validation performed in sandbox

- Python compile check passed for new/modified V7 files.
- V7 dry-run command chain passed.
- Synthetic final-selector smoke test produced 6 final candidates with real-structure mock evidence.
- Synthetic final audit passed with only a low-risk portfolio warning.

## Recommended deployment mode

Use V7.1 as the deployable default. Keep `RUN_V7_SOFTCHECK=1` if OpenComplex2 or another external structure backend is available. If not available, V7.1 still uses biology soft-risk scoring and the original ColabFold/PDB hard gates.
