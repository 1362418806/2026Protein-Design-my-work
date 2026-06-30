# SynBio GFP V6 Changelog

## V6 traceable release

- Added a clean `synbio_gfp_v6` Python namespace copied from the V4/V3 implementation for V6 isolation.
- Added `scripts/run_v6_complete.py`, using V6 output names and V6 config keys.
- Added `configs/v6_complete.yaml` with V6 gates, traceability keys, and ColabFold as the default structure backend label.
- Added `scripts/run_v6_pipeline.py` as a full-stage orchestrator:
  - `00_preflight`
  - `01_stage1`
  - `02_structure`
  - `03_proteinmpnn`
  - `04_final`
  - `05_lineage_index`
- Added `tools/run_structure_backend_v6.py` for backend-agnostic PDB回灌 and structure metrics generation.
- Added `tools/run_proteinmpnn_v6.py` for traceable ProteinMPNN score_only CSV generation.
- Added `scripts/v6_lineage.py` for artifact hashing and file-level lineage indexing.
- Added fixed-path environment setup under `/hyperai/home/envs/synbio_gfp_v6`.
- Added reusable third-party installation roots under `/hyperai/home/tools`.
- Added deployment helpers under `deploy/`:
  - `env_v6.sh`
  - `install_v5_env.sh`
  - `00_preflight_v6.py`
  - `run_v6_all.sh`
  - `run_v6_stage1.sh`
  - `run_v6_final.sh`
- Added `README_V5_TRACEABLE.md`.
