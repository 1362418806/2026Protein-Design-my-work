# Deployment patch changelog

## What changed

- Added `requirements-v4.txt` and `environment-v4.yml` for reproducible Python dependency setup.
- Added `deploy/00_preflight.py` to check packages, CUDA visibility, ESM availability, and official data files.
- Added Stage 1 / Stage 2 shell wrappers for AF2/AF3-missing servers:
  - `deploy/run_v4_stage1_export.sh`
  - `deploy/run_v4_stage2_import_metrics.sh`
  - `deploy/run_v4_debug_proxy.sh`
- Added `deploy/make_structure_metrics_template.py` for creating the real-structure metrics import template.
- Added `deploy/README_DEPLOYMENT.md` with HyperAI/Bohrium-style deployment instructions.
- Patched `scripts/run_v4_complete.py` structure-metrics import compatibility with recent pandas versions.

## Scope control

The patch does not modify V4 candidate generation, model training, scoring weights, basic gate thresholds, structure gate thresholds, or portfolio selection logic.
