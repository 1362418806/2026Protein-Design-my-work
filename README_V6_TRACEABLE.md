# SynBio GFP V6 Traceable Pipeline

V6 is a production-oriented wrapper and code refresh built from the V4 workflow. Its goal is not to hide scientific decisions, but to make every stage reproducible, auditable, and easy to rerun.

## Fixed paths

```bash
V6_DIR=/hyperai/home/synbio_gfp_v61_complete
V6_ENV_PREFIX=/hyperai/home/envs/synbio_gfp_v6
TOOLS_ROOT=/hyperai/home/tools
DATA_DIR="/hyperai/input/input0/2026Protein Design"
V6_OUTPUT_ROOT=/hyperai/home/synbio_gfp_v61_complete/outputs
V6_CACHE_DIR=/hyperai/home/synbio_gfp_v61_complete/cache/v5
```

Source the shared environment file before running commands:

```bash
source /hyperai/home/synbio_gfp_v61_complete/deploy/env_v6.sh
```

## Install environment

```bash
cd /hyperai/home/synbio_gfp_v61_complete
bash deploy/install_v5_env.sh
source /hyperai/home/envs/synbio_gfp_v6/activate_v5.sh
```

The installer creates a fixed conda/micromamba prefix under `/hyperai/home/envs/synbio_gfp_v6`, uses `/hyperai/home/pip_cache/synbio_gfp_v6`, and places reusable tools under `/hyperai/home/tools`.

Optional knobs:

```bash
INSTALL_LOCALCOLABFOLD=1 INSTALL_PROTEINMPNN=1 INSTALL_TORCH_MODE=cu128 bash deploy/install_v5_env.sh
INSTALL_LOCALCOLABFOLD=0 INSTALL_PROTEINMPNN=0 INSTALL_TORCH_MODE=cpu bash deploy/install_v5_env.sh
```

## One-command traceable run

For a full lightweight run that does not launch structure prediction:

```bash
TEAM_NAME=YourTeamName \
FEATURE_MODE=simple \
STRUCTURE_RUNNER=metrics-only \
RUN_PROTEINMPNN=0 \
bash deploy/run_v6_all.sh
```

For ColabFold structure prediction after LocalColabFold is available:

```bash
TEAM_NAME=YourTeamName \
FEATURE_MODE=simple \
STRUCTURE_RUNNER=colabfold \
REF_PDB=/hyperai/home/synbio_gfp_v61_complete/reference/reference_gfp.pdb \
RUN_PROTEINMPNN=1 \
bash deploy/run_v6_all.sh
```

## Output layout

Each run is written as:

```text
outputs/<RUN_ID>/
  00_preflight/
  01_stage1/
  02_structure/
  03_proteinmpnn/
  04_final/
  lineage/
  logs/
```

Important files:

```text
01_stage1/all_ranked_candidates_v6.csv
01_stage1/structure_priority_top200_v6.csv
01_stage1/af2_af3_structure_tasks.fasta
02_structure/structure_metrics.csv
03_proteinmpnn/proteinmpnn_scores.csv
04_final/final_top6_v6.csv
04_final/submission_v6.csv
04_final/v6_gate_diagnostics.csv
lineage/v5_stage_status.csv
lineage/v5_commands.sh
lineage/v5_artifact_index.csv
```

## Structure backend contract

Any backend can be used if it emits PDB files. The V6 adapter converts PDBs into the common `structure_metrics.csv` schema:

```csv
candidate_id,Seq_ID,sequence,mutations,v5_objective,barrel_rmsd,pocket_rmsd,ion_network_rmsd,mean_plddt,pdb_path,structure_status,structure_pass
```

This allows ColabFold, OpenFold, OpenComplex, AlphaFold3, or manually generated PDB files to be回灌 into the same final reranking stage.

## Traceability rules

1. Every stage writes a log under `logs/`.
2. Every executed command is stored in `lineage/v5_commands.sh`.
3. Stage status is stored in `lineage/v5_stage_status.csv`.
4. All CSV/JSON/FASTA/PDB/log/report artifacts are hashed into `lineage/v5_artifact_index.csv`.
5. Final candidates can be traced backward through `candidate_id` into Stage 1, structure metrics, ProteinMPNN scores, and final diagnostics.
