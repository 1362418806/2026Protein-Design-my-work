# V4 管线服务器部署说明

本补丁只补齐部署链路，不修改候选生成、模型训练、评分权重或最终选择算法主体。V4 的正确部署方式是 **两阶段运行**：先在当前服务器完成 ESM/候选/基础门控并导出结构任务，再把 AF2/AF3/ColabFold/ProteinMPNN 的真实结果以 CSV 或 PDB 目录形式导回 V4。

## 1. 环境安装

```bash
cd /hyperai/home/synbio_gfp_v4_complete
bash deploy/01_setup_env.sh
```

如果服务器无法安装 GPU 版 torch，可以先使用 `FEATURE_MODE=simple` 做 smoke test；正式 ESM 运行仍建议使用 CUDA。

## 2. 预检

```bash
python deploy/00_preflight.py \
  --project-root /hyperai/home/synbio_gfp_v4_complete \
  --data-dir "/hyperai/input/input0/2026Protein Design" \
  --require-esm
```

预检会检查核心 Python 包、FAIR-ESM、CUDA 可见性，以及官方数据包中的：

- `AAseqs of 5 GFP proteins.txt`
- `Exclusion_List.csv`
- `GFP_data.xlsx`

## 3. Stage 1：导出结构任务

```bash
PROJECT_ROOT=/hyperai/home/synbio_gfp_v4_complete \
DATA_DIR="/hyperai/input/input0/2026Protein Design" \
TEAM_NAME=YourTeamName \
FEATURE_MODE=esm \
N_CANDIDATES=20000 \
bash deploy/run_v4_stage1_export.sh
```

Stage 1 会生成：

- `outputs/run_v4_stage1/structure_priority_top200_v4.csv`
- `outputs/run_v4_stage1/af2_af3_structure_tasks.fasta`
- `outputs/run_v4_stage1/af2_af3_structure_tasks.jsonl`
- `outputs/run_v4_stage1/proteinmpnn_design_tasks.fasta`
- `outputs/run_v4_stage1/proteinmpnn_design_tasks.jsonl`

在没有真实结构指标时，`submission_v4.csv` 可能为空，这是预期行为，不是报错。

## 4. AF2/AF3 未编译时的推荐落地方式

不要在主 pipeline 内强行编译 AF2/AF3。推荐将 V4 视为结构任务调度器：

1. 用 Stage 1 导出的 `af2_af3_structure_tasks.fasta` 交给已有 AF2/AF3/ColabFold 环境。
2. 对得到的 PDB 计算或整理以下列：
   - `candidate_id` 或 `Seq_ID` 或 `sequence`
   - `barrel_rmsd`
   - `pocket_rmsd`
   - `ion_network_rmsd`
   - `structure_status`
   - `structure_pass`
3. 保存为 `structure_metrics.csv`。
4. Stage 2 导入该 CSV 生成最终 Top-6。

可先创建模板：

```bash
python deploy/make_structure_metrics_template.py \
  --priority-csv outputs/run_v4_stage1/structure_priority_top200_v4.csv \
  --out outputs/run_v4_stage1/structure_metrics_template.csv
```

## 5. Stage 2：导入真实结构指标并生成最终提交

```bash
PROJECT_ROOT=/hyperai/home/synbio_gfp_v4_complete \
DATA_DIR="/hyperai/input/input0/2026Protein Design" \
TEAM_NAME=YourTeamName \
FEATURE_MODE=esm \
STRUCTURE_METRICS_CSV=/path/to/structure_metrics.csv \
bash deploy/run_v4_stage2_import_metrics.sh
```

若已有 ProteinMPNN 分数，可追加：

```bash
PROTEINMPNN_SCORE_CSV=/path/to/proteinmpnn_scores.csv
```

最终输出：

- `outputs/run_v4_final/final_top6_v4.csv`
- `outputs/run_v4_final/submission_v4.csv`
- `outputs/run_v4_final/v4_gate_diagnostics.csv`
- `outputs/run_v4_final/v4_pipeline_report.md`

## 6. Debug-only 代理运行

仅用于确认管线部署，不用于最终提交：

```bash
PROJECT_ROOT=/hyperai/home/synbio_gfp_v4_complete \
DATA_DIR="/hyperai/input/input0/2026Protein Design" \
FEATURE_MODE=simple \
bash deploy/run_v4_debug_proxy.sh
```

该模式使用 `--allow-proxy-final`，会绕过真实结构门控。正式提交前必须改回 Stage 2 导入真实结构指标。
