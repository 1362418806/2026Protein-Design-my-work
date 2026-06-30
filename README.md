# SynBio GFP V6.1 Adaptive Pipeline

This package is a compute-adaptive upgrade of the V6 SynBio GFP closed-loop design pipeline.
It preserves V6 feedback-conditioned generation and adds:

- GPU profiles: `auto`, `rtx3090_24g`, `rtx5090`, `a800`.
- DATA_DIR priority for `/input0/2026Protein Design`.
- Exact-sequence PDB registry and reuse-first structure resolution.
- Recomputed structure metrics on every run.
- Optional V5.1 brightness-safety teacher via `--brightness-teacher-csv`.
- Physical thermal proxy scoring.
- Feedback survival audit by source branch.
- Stricter final Top-6 portfolio constraints.
- One-click and stage-wise execution scripts.

Read `README_V6.1_ADAPTIVE.md` and `SynBio_GFP_V6.1_Adaptive_维护文档.md` first.
