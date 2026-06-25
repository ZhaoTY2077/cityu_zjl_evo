# RBP 跨灵长类进化保守性分析 — 设计规格

## 概述

在已有的人类 RBP 功能富集分析基础上，建立可扩展的序列进化分析管道，对比人类与灵长类动物的 RBP 序列，识别保守性和进化性特征。

**设计目标：**
- **Phase 1（本次实现）：** 人类 ↔ 黑猩猩 RBP 序列比对 + dN/dS 选择压力 + 保守性分析管道
- **Phase 2（未来）：** 扩展到更多灵长类（大猩猩、猩猩、猕猴、狨猴等）
- **Phase 3（未来）：** 交互式浏览工具

## 架构

### 技术栈

| 层 | 语言 | 职责 |
|----|------|------|
| **数据管道层** | Python (Biopython) | 序列预处理、格式转换、批量比对、dN/dS 计算、保守性打分 |
| **统计分析层** | R (tidyverse, ggplot2, clusterProfiler) | 统计汇总、出版级可视化、功能富集分析 |
| **外部工具** | MAFFT, PAML (codeml) | 序列比对、选择压力计算 |

### 数据流

```
原始 FASTA (human + chimp) + 同源映射表
        │
        ▼
┌─────────────────── Python ───────────────────┐
│  P1: 过滤 + 配对 (最长转录本)                   │
│  P2: MAFFT 成对比对 (cDNA + 蛋白)              │
│  P3: PAML codeml → dN/dS                      │
│  P4: 保守性打分 (p-distance, entropy)           │
└──────────────────────┬────────────────────────┘
         中间文件: CSV / FASTA / JSON
                       │
                       ▼
┌──────────────────── R ───────────────────────┐
│  R1: 数据导入与统计汇总                          │
│  R2: 保守性可视化 (ω分布/散点图/箱线图)          │
│  R3: 富集分析 (复用现有 enrich_chimp_100pct.R)  │
│  R4: 整合报告                                   │
└──────────────────────┬────────────────────────┘
                       ▼
               图表 / CSV / HTML 报告
```

### 目录结构

```
evo_project/
├── code/
│   ├── enrichment_analysis.R            # 现有 R 富集分析脚本
│   ├── enrichment_chimp_100pct.R        # 现有 R 黑猩猩富集分析脚本
│   ├── download_RBP_sequences.py        # 现有 Ensembl 序列下载脚本
│   ├── check_rbp_completeness.py
│   ├── download_missing_RBP.py
│   ├── download_remaining_RBP.py
│   └── evolutionary/                    # 进化分析管道代码
│       ├── config.yaml                  # 物种、路径、参数 — 唯一配置入口
│       ├── bin/                         # 外部工具脚本 (pal2nal 等)
│       ├── python/
│       │   ├── 01_filter_and_pair.py
│       │   ├── 02_run_alignment.py
│       │   ├── 03_calculate_dnds.py
│       │   └── 04_compute_conservation.py
│       └── r/
│           ├── 01_import_summary.R
│           ├── 02_visualize_conservation.R
│           ├── 03_enrichment_comparison.R
│           └── 04_render_report.R
├── data/
│   └── raw/                             # 原始下载数据 (FASTA, 同源映射表)
├── results/
│   ├── enrichment_analysis/             # 51 基因富集分析结果
│   └── enrichment_chimp_100pct/         # 391 基因 100% 保守 RBP 富集结果
├── docs/
│   └── superpowers/
│       ├── specs/                       # 设计规格文档
│       └── plans/                       # 实施计划文档
└── README.md
```

## Python 层模块设计

### P1: 序列预处理与配对

**输入：**
- `files/human_RBP_cdna.fasta`, `files/human_RBP_protein.fasta`
- `files/pan_troglodytes_RBP_cdna.fasta`, `files/pan_troglodytes_RBP_protein.fasta`
- `files/chimp_ortholog_mapping.txt`

**处理逻辑：**
1. 用 `BioPython SeqIO` 读取 FASTA
2. 对每个 ENSG 基因 ID，选择最长转录本 / 最长蛋白作为代表性序列
3. 过滤：翻译检查排除含提前终止密码子的序列；跳过非 `ortholog_one2one` 的映射
4. 配对输出：`data/paired_sequences/{gene_symbol}_hsa_ptr.fa`

**输出：**
- `data/paired_sequences/` — 每个 RBP 基因一对序列的 FASTA 文件
- `results/filtering_stats.csv` — 每个基因的处理状态

**未来扩展点：** `config.yaml` 的 `species` 列表决定输入文件命名，新增物种时映射规则从 `{species}_ortholog_mapping.txt` 读取。

### P2: 序列比对

**输入：** `data/paired_sequences/` 下每个基因的配对序列

**处理逻辑：**
1. 调用 `mafft` 进行成对蛋白比对（`--auto` 模式）
2. 调 `mafft` 进行成对 cDNA 比对（`--globalpair --maxiterate 1000`）
3. 对于不能翻译对齐的区域，用 `pal2nal` 工具将蛋白比对映射回密码子比对（为 PAML 准备）

**输出：**
- `data/alignments/protein/{gene}.aln`
- `data/alignments/cdna/{gene}.aln`
- `alignments_summary.csv` — 比对长度、gap 百分比

**未来扩展点：** 增加多物种时，比对模式从 pairwise 变为以人类为参考的多序列比对（Profile alignment）。

### P3: dN/dS 选择压力分析

**输入：** `data/alignments/cdna/` 下的密码子比对

**处理逻辑：**
1. 为每个基因生成 `codeml.ctl` 控制文件（成对模式 `runmode = -2`，`seqtype = 1`）
2. 调用 PAML `codeml`
3. 解析输出结果文件（`rst`, `mlc`），提取 dN, dS, ω
4. 计算显著性：似然比检验（零假设 ω=1 vs 自由模型）

**输出：**
- `results/dn_ds_results.csv` — 字段: `gene, dN, dS, omega, lnL_null, lnL_alt, p_value, significant`
- `results/dnds_summary.csv` — 全局统计

**未来扩展点：** 多物种时改用自由比模型（`model = 1`）+ 分支模型（`model = 2`），用系统发育树估算每个分支的 ω。

### P4: 保守性打分

**输入：** `data/alignments/protein/` 下的蛋白比对

**处理逻辑：**
1. **序列一致性：** 每个基因的 protein %identity、cDNA %identity（p-distance）
2. **位点熵值：** 基于比对计算每个位点的 Shannon entropy
3. **残基置换打分：** 基于 BLOSUM62 矩阵的平均置换得分（正值=保守替换，负值=非常规替换）
4. **滑动窗口保守性：** 窗口大小 10aa，滑动步长 1，输出每个窗口的平均保守性

**输出：**
- `results/conservation_scores.csv` — 字段: `gene, prot_identity, cdna_identity, min_entropy, mean_entropy, mean_blosum62_score, n_sites`
- `results/sitewise_conservation/{gene}.csv` — 每个位点的细粒度数据

**未来扩展点：** 多物种时用 Jensen-Shannon divergence 计算进化速率差异。

## R 层模块设计

### R1: 数据汇总

**输入：** Python 输出的所有 CSV 文件

**输出：** 控制台报告 + `results/analysis_summary.csv`

- 成功配对的基因数 / 比对成功数 / dN/dS 计算结果数
- 筛选流程中的样本量流失日志

### R2: 可视化

所有图表以 `theme_minimal()` 为基，颜色和字体对齐现有 `enrichment_analysis.R` 风格。

| 图表 | 类型 | 变量 |
|------|------|------|
| ω 分布 | 直方图 + 密度曲线 | dN/dS (omega)，标注 ω=1 阈值线 |
| 蛋白一致性 vs ω | 散点图 | x=prot_identity, y=omega, 颜色=显著性 |
| 滑动窗口保守性 | 热图/线图 | 按基因分面 |dN vs dS | 散点图 | x=dS, y=dN, 颜色=ω |
| 基因功能分类保守性 | 箱线图 | x=RBP_class, y=prot_identity |

### R3: 扩展的富集分析

复用现有 `enrichment_chimp_100pct.R` 的分析框架：

- **高保守基因集**（ω < 0.1, prot_identity > 99%）：做 GO/KEGG 富集
- **正选择基因集**（ω > 1, p < 0.05）：做 GO/KEGG 富集
- 与现有的 100% 保守基因富集结果对比

### R4: 整合报告（可选）

生成 R Markdown 格式的分析报告，包含所有图表和统计表。

## 配置入口

所有可变参数集中在 `config.yaml`，`species` 列表设计为可扩展：

```yaml
project:
  name: "RBP_evolutionary_analysis"
  root: "/saturn/zhaoty/evo_project/evolutionary_analysis"

species:
  reference: "human"
  targets:
    - code: "ptr"
      name: "Pan troglodytes"
      cdna: "../data/raw/pan_troglodytes_RBP_cdna.fasta"
      protein: "../data/raw/pan_troglodytes_RBP_protein.fasta"
      ortholog_map: "../data/raw/chimp_ortholog_mapping.txt"
  # 未来扩展:
  # - code: "ggo"
  #   name: "Gorilla gorilla"
  #   ...

alignment:
  mafft_opts: "--globalpair --maxiterate 1000"

dnds:
  paml_path: "codeml"

filtering:
  longest_transcript_only: true
  skip_non_one2one: true
```

更改 `targets` 列表即可增加新物种。

## 执行流程

```bash
# 顺序执行（Python → R）
cd /saturn/zhaoty/evo_project/evolutionary_analysis

# Python 管道
python python/01_filter_and_pair.py
python python/02_run_alignment.py
python python/03_calculate_dnds.py
python python/04_compute_conservation.py

# R 分析
Rscript r/01_import_summary.R
Rscript r/02_visualize_conservation.R
Rscript r/03_enrichment_comparison.R
# Rscript r/04_render_report.R  (可选)
```

每个步骤可以单独运行，调试时不需要重跑前面步骤。

## 未来扩展构架

| 阶段 | 内容 | 改动范围 |
|------|------|---------|
| **Phase 2** | 添加大猩猩、猕猴等 | 下载FASTA → 写同源映射 → 改 `config.yaml` → 重跑管道 |
| **Phase 2** | 多物种 dN/dS | P3 新增分支模型；R2 新增 `ggtree` 进化树标注 |
| **Phase 3** | 交互式展示 | 以 `results/` 下 CSV 为数据源，Shiny 或 Streamlit 实现 |

## 不纳入范围（YAGNI）

- 不做浏览器端的交互式分析（Phase 3 再做）
- 不做新的序列下载工具（现有 `files/` 中的脚本已覆盖）
- 不做并行计算框架（序列量级用单机足够）
- 不做实时 Web 服务
