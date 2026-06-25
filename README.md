# RBP 跨物种进化保守性分析

RNA结合蛋白（RBP）在人类与灵长类动物之间的序列保守性与进化特征分析。

## 项目结构

```
evo_project/
├── code/                              # 源代码
│   ├── enrichment_analysis.R          # 51 基因 GO/KEGG 富集分析
│   ├── enrichment_chimp_100pct.R      # 391 个 100% 保守 RBP 的富集分析
│   ├── download_RBP_sequences.py      # Ensembl API 序列下载工具
│   ├── check_rbp_completeness.py      # 序列完整性检查
│   ├── download_missing_RBP.py        # 缺失序列补充下载
│   ├── download_remaining_RBP.py      # 剩余序列下载
│   └── evolutionary/                  # 进化分析管道（开发中）
│       ├── config.yaml                # 管道配置
│       ├── bin/                       # 外部工具 (pal2nal 等)
│       ├── python/                    # Python 序列处理模块
│       └── r/                         # R 统计与可视化模块
├── data/
│   └── raw/                           # 原始数据
│       ├── human_RBP_cdna.fasta       # 人类 RBP cDNA 序列 (62,344 条)
│       ├── human_RBP_protein.fasta    # 人类 RBP 蛋白序列 (50,422 条)
│       ├── pan_troglodytes_RBP_cdna.fasta  # 黑猩猩 RBP cDNA 序列 (3,558 条)
│       ├── pan_troglodytes_RBP_protein.fasta # 黑猩猩 RBP 蛋白序列 (3,558 条)
│       ├── chimp_ortholog_mapping.txt # 人↔黑猩猩直系同源映射 (1,427 对)
│       ├── chimp_100pct_genes.txt     # 100% 保守的 391 个 RBP 基因
│       ├── RBP_geneid_list.txt        # 人类 RBP 基因 ENSG ID 列表
│       ├── human_no_chimp_ortholog_list.txt  # 无黑猩猩同源的人类 RBP
│       └── human_obsolete_RBP_ids.txt # 已废弃/合并的 RBP 基因 ID
├── results/
│   ├── enrichment_analysis/           # 51 基因富集分析结果
│   │   ├── GOBP_full.csv             # GO 生物过程富集结果
│   │   ├── GOMF_full.csv             # GO 分子功能富集结果
│   │   ├── GOCC_full.csv             # GO 细胞组分富集结果
│   │   ├── KEGG_full.csv             # KEGG 通路富集结果
│   │   └── *.pdf                     # 可视化图表
│   └── enrichment_chimp_100pct/      # 100% 保守 RBP 富集结果
│       ├── GOBP_full.csv
│       ├── GOMF_full.csv
│       ├── GOCC_full.csv
│       ├── KEGG_full.csv
│       └── *.pdf
├── docs/
│   └── superpowers/
│       ├── specs/                     # 设计规格文档
│       └── plans/                     # 实施计划文档
└── README.md
```

## 数据来源

- **人类 RBP 基因列表：** 从已知 RNA 结合蛋白数据库提取
- **序列数据：** 通过 Ensembl REST API（v110）下载（`code/download_RBP_sequences.py`）
- **人↔黑猩猩同源映射：** Ensembl Compara 管道获取

## 分析管道

### Phase 1: 功能富集分析（已完成）

| 脚本 | 分析内容 | 基因集 |
|------|---------|--------|
| `code/enrichment_analysis.R` | GO (BP/MF/CC) + KEGG 超几何富集 | 51 个人类 RBP 基因 |
| `code/enrichment_chimp_100pct.R` | GO (BP/MF/CC) + KEGG 超几何富集 | 391 个人类-黑猩猩 100% 保守 RBP 基因 |

**运行方法：**
```bash
cd /saturn/zhaoty/evo_project
Rscript code/enrichment_analysis.R
Rscript code/enrichment_chimp_100pct.R
```

### Phase 2: 序列进化分析（计划中）

人类↔黑猩猩 RBP 序列配对 → MAFFT 比对 → PAML dN/dS 选择压力 → 保守性打分 → 可视化

详见 `docs/superpowers/specs/2026-06-25-rbp-evolutionary-analysis-design.md`

### Phase 3: 多灵长类扩展（未来）

扩展到黑猩猩、大猩猩、猩猩、猕猴、狨猴等多物种比较。

## 结果摘要

- **51 基因分析：** 包含部分无黑猩猩同源物的基因，用于功能富集
- **391 基因分析：** 人与黑猩猩蛋白序列 100% 一致的 RBP 基因，揭示高度保守的 RNA 结合功能通路

## 环境依赖

**R (≥ 4.3):**
- org.Hs.eg.db, GO.db, KEGGREST, ggplot2, AnnotationDbi, reshape2, DOSE

**Python:**
- biopython, pyyaml, requests

**外部工具（进化分析用）:**
- MAFFT (序列比对)
- PAML codeml (dN/dS 计算)
- pal2nal (蛋白-密码子比对映射)

## 许可证

学术研究用途。
