# SynBio GFP V4 Complete Pipeline

This package is a non-invasive V4 extension of the uploaded V3 project. The original V3 source tree is preserved; V4 is added as:

- `configs/v4_complete.yaml`
- `scripts/run_v4_complete.py`

## Main V4 additions

1. **LCB objective**: `brightness_lcb = brightness_mean - lambda_std * brightness_std`.
2. **Strict two-stage design**:
   - Stage A: V3-compatible generation, prediction, literature priors, and V4 basic gate.
   - Stage B: AF2/AF3/PDB/structure-metrics validation before final Top-6.
3. **ProteinMPNN branch**:
   - Exports `proteinmpnn_design_tasks.fasta` and `proteinmpnn_design_tasks.jsonl`.
   - Can import external ProteinMPNN scores through `--proteinmpnn-score-csv`.
4. **AF2/AF3 branch**:
   - Exports `af2_af3_structure_tasks.fasta` and `af2_af3_structure_tasks.jsonl`.
   - Can import real structure metrics through `--structure-metrics-csv` or parse PDBs from `--reference-pdb` + `--predicted-pdb-dir`.
5. **V4 gates and diagnostics**:
   - `v4_basic_gate`
   - `v4_structure_gate`
   - `v4_final_eligible`
   - `v4_gate_diagnostics.csv`
6. **Analysis outputs**:
   - `all_ranked_candidates_v4.csv`
   - `ranked_top200_v4.csv`
   - `structure_priority_top200_v4.csv`
   - `structure_gate_top200_v4.csv` (backward-compatible alias)
   - `final_top6_v4.csv`
   - `submission_v4.csv`
   - `mutation_frequency_by_seed_v6.csv`
   - `rejected_by_basic_gate.csv`
   - `rejected_by_structure_gate.csv`
   - `v4_pipeline_report.md`

## Smoke test

```bash
cd synbio_gfp_v4_complete
python tests/make_synthetic_data.py
python scripts/run_v4_complete.py \
  --config configs/v4_complete.yaml \
  --data-dir tests/synthetic_data \
  --team-name SmokeTeam \
  --feature-mode simple \
  --max-train-samples 250 \
  --n-candidates 1000 \
  --out-dir tests/smoke_v4 \
  --cache-dir tests/cache_v4
```

Expected behavior without real AF/PDB metrics: `final_top6_v4.csv` is intentionally empty, while `structure_priority_top200_v4.csv` and AF/ProteinMPNN task files are populated.

## Production run with structure metrics CSV

The structure metrics CSV should share one of `candidate_id`, `Seq_ID`, or `sequence`, and may contain:

- `barrel_rmsd`
- `pocket_rmsd`
- `ion_network_rmsd`
- `structure_status`
- `structure_pass`

```bash
python scripts/run_v4_complete.py \
  --config configs/v4_complete.yaml \
  --data-dir /path/to/data \
  --team-name YourTeamName \
  --feature-mode esm \
  --max-train-samples 20000 \
  --n-candidates 20000 \
  --out-dir outputs/run_v4_complete \
  --cache-dir cache \
  --structure-metrics-csv /path/to/structure_metrics.csv \
  --proteinmpnn-score-csv /path/to/proteinmpnn_scores.csv
```

## Debug-only proxy final

Use only for pipeline testing, not final biological selection:

```bash
python scripts/run_v4_complete.py ... --allow-proxy-final
```

## Server deployment patch

For HyperAI/Bohrium-style servers where AF2/AF3 is not compiled locally, use the deployment helpers under `deploy/`:

1. `deploy/00_preflight.py` checks Python packages, CUDA visibility, ESM availability, and official data files.
2. `deploy/run_v4_stage1_export.sh` runs V4 through candidate scoring and exports AF2/AF3 + ProteinMPNN task files.
3. `deploy/make_structure_metrics_template.py` creates a metrics CSV template for external AF/PDB results.
4. `deploy/run_v4_stage2_import_metrics.sh` imports real structure metrics and emits the final `submission_v4.csv`.
5. `deploy/run_v4_debug_proxy.sh` is only for deployment smoke tests and must not be used for final biological submission.

See `deploy/README_DEPLOYMENT.md` for the full deployment procedure.
