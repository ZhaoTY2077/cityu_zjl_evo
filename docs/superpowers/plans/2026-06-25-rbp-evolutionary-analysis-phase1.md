# RBP 进化保守性分析管道 — Phase 1 实施计划

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 构建人类↔黑猩猩 RBP 序列比对、dN/dS 选择压力分析、保守性打分的可复现管道，并产出出版级可视化和富集分析报告。

**Architecture:** 混合管道 — Python (Biopython) 处理序列层面的批量计算，R (ggplot2/clusterProfiler) 负责统计分析和可视化。中间数据以 CSV/FASTA 文件传递，每个步骤可独立运行。

**Tech Stack:** Python (Biopython, PyYAML), R (ggplot2, clusterProfiler, tidyverse, rmarkdown), MAFFT, PAML (codeml), pal2nal

## Global Constraints

- 所有路径读取从 `config.yaml` 出发，不硬编码
- 基因过滤标准：选择最长转录本作为代表性序列；仅保留 `ortholog_one2one` 关系
- 中间文件保存在 `data/` 目录下，按类型分子目录
- 图形配色对齐现有 `enrichment_analysis.R` 风格（BP=#E64B35, MF=#4DBBD5, CC=#00A087, KEGG=#3C5488）
- 不引入新的大框架依赖（如 Dask、Spark）
- Python 脚本可单独运行，内部通过 `argparse` 接收 `--config` 参数

---

### Task 0: 环境准备和目录初始化

**Files:**
- Create: `evolutionary_analysis/config.yaml`
- Create: `evolutionary_analysis/data/` (及子目录)
- Create: `evolutionary_analysis/figures/`
- Create: `evolutionary_analysis/report/`
- None: 安装依赖

**Interfaces:**
- Produces: `config.yaml` — 一切后续任务读取的唯一配置入口

- [ ] **Step 1: 安装依赖**

```bash
# Python
pip install biopython pyyaml

# R 包
Rscript -e 'install.packages(c("tidyverse", "DT"), repos="https://cloud.r-project.org")'
Rscript -e 'if (!require("BiocManager", quietly=TRUE)) install.packages("BiocManager"); BiocManager::install("clusterProfiler")'

# 外部工具 (conda)
conda install -c bioconda mafft paml
# pal2nal 是一个 Perl 脚本，直接下载
cd /saturn/zhaoty/evo_project/evolutionary_analysis
wget https://raw.githubusercontent.com/malek-lab/pal2nal/master/pal2nal.pl -O bin/pal2nal.pl
chmod +x bin/pal2nal.pl
```

- [ ] **Step 2: 创建目录结构**

```bash
cd /saturn/zhaoty/evo_project
mkdir -p evolutionary_analysis/{python,r,data/{paired_sequences,alignments/{protein,cdna},results,sitewise_conservation},figures,report,bin}
```

- [ ] **Step 3: 创建 config.yaml**

```yaml
project:
  name: "RBP_evolutionary_analysis"
  root: "/saturn/zhaoty/evo_project/evolutionary_analysis"

species:
  reference: "human"
  ref_cdna: "../files/human_RBP_cdna.fasta"
  ref_protein: "../files/human_RBP_protein.fasta"
  targets:
    - code: "ptr"
      name: "Pan troglodytes"
      cdna: "../files/pan_troglodytes_RBP_cdna.fasta"
      protein: "../files/pan_troglodytes_RBP_protein.fasta"
      ortholog_map: "../files/chimp_ortholog_mapping.txt"
      ortholog_sep: "\t"          # 同源映射表分隔符

alignment:
  mafft_path: "mafft"
  mafft_opts_protein: "--auto"
  mafft_opts_cdna: "--globalpair --maxiterate 1000"
  pal2nal_path: "bin/pal2nal.pl"

dnds:
  paml_path: "codeml"
  temp_dir: "data/tmp_codeml"

filtering:
  longest_transcript_only: true
  skip_non_one2one: true
  min_protein_length: 50

paths:
  enst_mapping: "data/enst_to_ensg.csv"
```

- [ ] **Step 4: 验证工具可用性并提交**

```bash
cd /saturn/zhaoty/evo_project
mafft --version
codeml 2>&1 | head -3
python3 -c "import Bio; print('Biopython', Bio.__version__)"
Rscript -e 'library(ggplot2); library(clusterProfiler); cat("R packages ok\n")'
git add evolutionary_analysis/
git commit -m "feat: init evolutionary analysis directory structure and config"
```

---

### Task 1: 下载 Ensembl 转录本↔基因映射表

**Files:**
- Create: `evolutionary_analysis/python/download_enst_mapping.py`

**Interfaces:**
- Consumes: FASTA 文件中的 ENST/ENSPTRT 转录本 ID
- Produces: `data/enst_to_ensg.csv` (human) 和 `data/enstptr_to_ensptrg.csv` (chimp)

**问题说明：** 现有 FASTA 文件的序列头仅包含转录本 ID（`>ENST00000600027.5`），不包含基因 ID。需要通过 Ensembl API 获取转录本→基因的映射关系。

- [ ] **Step 1: 创建映射下载脚本**

```python
#!/usr/bin/env python3
"""
Download Ensembl transcript-to-gene mapping for transcript IDs found in FASTA files.
Uses REST API batch lookup (50 IDs/batch).
"""
import requests
import json
import time
import os
import csv
import re
import sys

BASE_URL = "https://rest.ensembl.org"
HEADERS = {"Content-Type": "application/json", "Accept": "application/json"}

def extract_transcript_ids(fasta_path):
    """Extract all unique transcript IDs (without version suffix) from a FASTA file."""
    ids = set()
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                tid = line[1:].strip()
                tid = re.sub(r'\.\d+$', '', tid)  # strip version
                ids.add(tid)
    return sorted(ids)

def batch_lookup(ids, batch_size=50):
    """Look up transcript IDs via Ensembl REST API in batches."""
    results = {}
    for i in range(0, len(ids), batch_size):
        batch = ids[i:i+batch_size]
        url = f"{BASE_URL}/lookup/id"
        resp = requests.post(url, headers=HEADERS, json={"ids": batch}, timeout=30)
        if resp.status_code == 200:
            data = resp.json()
            for tid, info in data.items():
                results[tid] = {
                    "gene_id": info.get("gene_id", ""),
                    "gene_name": info.get("display_name", ""),
                }
        else:
            print(f"  Batch {i//batch_size} failed: {resp.status_code}")
        time.sleep(0.5)  # rate limit
        if (i // batch_size + 1) % 10 == 0:
            print(f"  Progress: {i+len(batch)}/{len(ids)}")
    return results

def save_mapping(results, output_path):
    """Save mapping as CSV."""
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["transcript_id", "gene_id", "gene_symbol"])
        for tid, info in sorted(results.items()):
            writer.writerow([tid, info["gene_id"], info["gene_name"]])

if __name__ == "__main__":
    project_root = "/saturn/zhaoty/evo_project/evolutionary_analysis"

    print("=== Human ENST → ENSG ===")
    human_fasta = os.path.join(project_root, "../files/human_RBP_cdna.fasta")
    human_ids = extract_transcript_ids(human_fasta)
    print(f"  Found {len(human_ids)} unique ENST IDs")
    human_results = batch_lookup(human_ids)
    human_out = os.path.join(project_root, "data/enst_to_ensg.csv")
    save_mapping(human_results, human_out)
    print(f"  Saved {len(human_results)} mappings to {human_out}")

    print("\n=== Chimp ENSPTRT → ENSPTRG ===")
    chimp_fasta = os.path.join(project_root, "../files/pan_troglodytes_RBP_cdna.fasta")
    chimp_ids = extract_transcript_ids(chimp_fasta)
    print(f"  Found {len(chimp_ids)} unique ENSPTRT IDs")
    chimp_results = batch_lookup(chimp_ids)
    chimp_out = os.path.join(project_root, "data/enstptr_to_ensptrg.csv")
    save_mapping(chimp_results, chimp_out)
    print(f"  Saved {len(chimp_results)} mappings to {chimp_out}")
```

- [ ] **Step 2: 运行下载并检查结果**

```bash
cd /saturn/zhaoty/evo_project/evolutionary_analysis
python python/download_enst_mapping.py
# 预期: ~62000 human + ~3558 chimp mappings
wc -l data/enst_to_ensg.csv
wc -l data/enstptr_to_ensptrg.csv
head -5 data/enst_to_ensg.csv
```

- [ ] **Step 3: 提交**

```bash
cd /saturn/zhaoty/evo_project
git add evolutionary_analysis/python/download_enst_mapping.py \
       evolutionary_analysis/data/enst_to_ensg.csv \
       evolutionary_analysis/data/enstptr_to_ensptrg.csv
git commit -m "feat: add Ensembl transcript-to-gene mapping download"
```

---

### Task 2 (P1): 序列预处理与配对

**Files:**
- Create: `evolutionary_analysis/python/01_filter_and_pair.py`

**Interfaces:**
- Consumes:
  - `config.yaml` — 所有参数
  - `files/chimp_ortholog_mapping.txt` — 同源映射关系
  - `files/human_RBP_cdna.fasta`, `files/human_RBP_protein.fasta`
  - `files/pan_troglodytes_RBP_cdna.fasta`, `files/pan_troglodytes_RBP_protein.fasta`
  - `data/enst_to_ensg.csv`, `data/enstptr_to_ensptrg.csv`
- Produces:
  - `data/paired_sequences/{gene_symbol}_hsa_ptr.fa`
  - `data/results/filtering_stats.csv`

- [ ] **Step 1: 编写过滤配对脚本**

```python
#!/usr/bin/env python3
"""
P1: Filter and pair human-chimp RBP sequences.
- Reads ortholog mapping to identify one2one pairs
- Picks longest transcript/protein per gene (no premature stop codons)
- Writes paired FASTA files for downstream analysis
"""
import os
import sys
import csv
import re
import yaml
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def load_config():
    with open("config.yaml") as f:
        return yaml.safe_load(f)

def load_mapping(csv_path):
    """Load transcript→gene mapping CSV, return {transcript_id: {gene_id, gene_symbol}}."""
    mapping = {}
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            mapping[row["transcript_id"]] = {
                "gene_id": row["gene_id"],
                "gene_symbol": row["gene_symbol"],
            }
    return mapping

def load_ortholog_map(path, sep="\t"):
    """Load ortholog mapping, return {human_ensg: {chimp_ensg, chimp_symbol, orthology_type, perc_id}}."""
    orth = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter=sep)
        for row in reader:
            orth[row["human_gene_id"]] = {
                "chimp_gene_id": row["chimp_gene_id"],
                "chimp_symbol": row["chimp_gene_symbol"],
                "orthology_type": row["orthology_type"],
                "perc_id": float(row["perc_id"]),
            }
    return orth

def group_transcripts_by_gene(fasta_path, enst_mapping):
    """Parse FASTA and group sequences by gene ID using transcript mapping."""
    gene_groups = defaultdict(list)  # ensg → [(transcript_id, seq_record)]
    for record in SeqIO.parse(fasta_path, "fasta"):
        tid = record.id.split(".")[0]  # strip version
        if tid in enst_mapping:
            gene_id = enst_mapping[tid]["gene_id"]
            if gene_id:
                gene_groups[gene_id].append((tid, record))
        else:
            # Try without stripping version
            tid_full = record.id
            if tid_full in enst_mapping:
                gene_id = enst_mapping[tid_full]["gene_id"]
                if gene_id:
                    gene_groups[gene_id].append((tid_full, record))
    return gene_groups

def pick_longest(records):
    """Pick the longest sequence from a list of (id, SeqRecord) tuples."""
    if not records:
        return None
    return max(records, key=lambda x: len(x[1].seq))

def has_premature_stop(seq):
    """Check if a CDS sequence has a premature stop codon in any frame."""
    if len(seq) < 3:
        return True
    # Try first frame only (assuming correct CDS)
    prot = seq.translate(to_stop=True)
    cds_len = len(prot) * 3
    if cds_len < len(seq) - 3:  # more than 3 extra bases after stop
        return True
    return False

def main():
    cfg = load_config()
    root = os.path.abspath(cfg["project"]["root"])
    os.chdir(root)

    # Paths
    ortholog_file = os.path.abspath(cfg["species"]["targets"][0]["ortholog_map"])
    human_cdna = os.path.abspath(cfg["species"]["ref_cdna"])
    human_prot = os.path.abspath(cfg["species"]["ref_protein"])
    chimp_cdna = os.path.abspath(cfg["species"]["targets"][0]["cdna"])
    chimp_prot = os.path.abspath(cfg["species"]["targets"][0]["protein"])
    enst_map_file = os.path.abspath(cfg["paths"]["enst_mapping"])
    enstptr_map_file = "data/enstptr_to_ensptrg.csv"
    out_dir = "data/paired_sequences"
    result_dir = "data/results"

    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(result_dir, exist_ok=True)

    # Load mappings
    print("Loading transcript→gene mappings...")
    human_map = load_mapping(enst_map_file)
    chimp_map = load_mapping(enstptr_map_file)

    print("Loading ortholog mapping...")
    ortholog_map = load_ortholog_map(ortholog_file, sep=cfg["species"]["targets"][0].get("ortholog_sep", "\t"))
    skip_non_one2one = cfg["filtering"]["skip_non_one2one"]
    min_prot_len = cfg["filtering"].get("min_protein_length", 50)

    # Group transcripts by gene
    print("Grouping human transcripts by gene...")
    human_gene_groups = group_transcripts_by_gene(human_cdna, human_map)
    print(f"  {len(human_gene_groups)} human genes with transcripts")

    print("Grouping human protein sequences by gene...")
    human_prot_groups = group_transcripts_by_gene(human_prot, human_map)
    # For proteins, transcripts may have different IDs; use transcript mapping
    # Map ENSP → ENST → ENSG using known ENSP-to-ENST relationships
    # (Ensembl protein IDs are derived from transcript IDs)
    prot_from_tid = {}
    for record in SeqIO.parse(human_prot, "fasta"):
        # ENSP IDs start with ENSP for human; strip version
        pid = record.id.split(".")[0]
        # The transcript ID is stored in the sequence ID; for Ensembl protein FASTA,
        # the header may contain the transcript ID after "transcript:" or similar
        prot_from_tid[pid] = record

    print("Grouping chimp transcripts by gene...")
    chimp_gene_groups = group_transcripts_by_gene(chimp_cdna, chimp_map)
    print(f"  {len(chimp_gene_groups)} chimp genes with transcripts")

    # Process each human gene with a one2one ortholog
    stats = []
    for human_ensg, orth_info in ortholog_map.items():
        gene_symbol = human_map.get(next(iter(human_gene_groups.get(human_ensg, [])), [None, None])[0], {}).get("gene_symbol", "")
        chimp_ensg = orth_info["chimp_gene_id"]
        status = "skipped"

        # Filter: one2one only
        if skip_non_one2one and orth_info["orthology_type"] != "ortholog_one2one":
            continue

        # Get human transcripts for this gene
        human_transcripts = human_gene_groups.get(human_ensg, [])
        if not human_transcripts:
            stats.append([human_ensg, gene_symbol, chimp_ensg, orth_info["chimp_symbol"],
                          orth_info["perc_id"], "no_human_transcripts"])
            continue

        # Get chimp transcripts for this gene
        chimp_transcripts = chimp_gene_groups.get(chimp_ensg, [])
        if not chimp_transcripts:
            stats.append([human_ensg, gene_symbol, chimp_ensg, orth_info["chimp_symbol"],
                          orth_info["perc_id"], "no_chimp_transcripts"])
            continue

        # Pick longest sequences
        best_human = pick_longest(human_transcripts)
        best_chimp = pick_longest(chimp_transcripts)
        if not best_human or not best_chimp:
            continue

        # Filter: check premature stop codons in cDNA
        if has_premature_stop(best_human[1].seq):
            status = "human_premature_stop"
            # Still include but flag it
        if has_premature_stop(best_chimp[1].seq):
            status = "chimp_premature_stop"

        # Filter: minimum protein length
        human_prot_len = len(best_human[1].seq.translate(to_stop=True))
        chimp_prot_len = len(best_chimp[1].seq.translate(to_stop=True))
        if human_prot_len < min_prot_len or chimp_prot_len < min_prot_len:
            stats.append([human_ensg, gene_symbol, chimp_ensg, orth_info["chimp_symbol"],
                          orth_info["perc_id"], f"short_protein({human_prot_len},{chimp_prot_len})"])
            continue

        # Write paired FASTA (both sequences in one file)
        if not gene_symbol:
            gene_symbol = human_ensg.split(".")[0]
        safe_name = re.sub(r'[^a-zA-Z0-9_-]', '_', gene_symbol)
        out_file = os.path.join(out_dir, f"{safe_name}_hsa_ptr.fa")

        with open(out_file, "w") as f:
            f.write(f">{best_human[1].id} [gene={human_ensg}] [species=human]\n")
            seq_str = str(best_human[1].seq)
            for i in range(0, len(seq_str), 80):
                f.write(seq_str[i:i+80] + "\n")
            f.write(f">{best_chimp[1].id} [gene={chimp_ensg}] [species=chimp]\n")
            seq_str = str(best_chimp[1].seq)
            for i in range(0, len(seq_str), 80):
                f.write(seq_str[i:i+80] + "\n")

        stats.append([human_ensg, gene_symbol, chimp_ensg, orth_info["chimp_symbol"],
                      orth_info["perc_id"], status])
        if len(stats) % 200 == 0:
            print(f"  Processed {len(stats)} genes...")

    # Write filtering stats
    out_csv = os.path.join(result_dir, "filtering_stats.csv")
    with open(out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["human_ensg", "human_symbol", "chimp_ensg", "chimp_symbol",
                          "perc_id", "status"])
        writer.writerows(stats)

    paired_count = sum(1 for s in stats if s[5] in ("ok", "human_premature_stop", "chimp_premature_stop"))
    print(f"\nDone. Total one2one orthologs: {len(stats)}")
    print(f"  Successfully paired: {paired_count}")
    print(f"  Results saved to: {out_csv}")

if __name__ == "__main__":
    main()
```

- [ ] **Step 2: 运行配对脚本**

```bash
cd /saturn/zhaoty/evo_project/evolutionary_analysis
python python/01_filter_and_pair.py
# 预期: ~1300 one2one pairs, ~1100-1200 successfully paired after filtering
wc -l data/results/filtering_stats.csv
ls data/paired_sequences/ | wc -l
head -10 data/results/filtering_stats.csv
```

- [ ] **Step 3: 提交**

```bash
cd /saturn/zhaoty/evo_project
git add evolutionary_analysis/python/01_filter_and_pair.py \
       evolutionary_analysis/data/paired_sequences/ \
       evolutionary_analysis/data/results/filtering_stats.csv
git commit -m "feat(p1): sequence filtering and pairing for human-chimp RBP orthologs"
```

---

### Task 3 (P2): 序列比对

**Files:**
- Create: `evolutionary_analysis/python/02_run_alignment.py`

**Interfaces:**
- Consumes: `data/paired_sequences/*.fa`
- Produces: `data/alignments/protein/{gene}.aln`, `data/alignments/cdna/{gene}.aln`, `data/results/alignments_summary.csv`
- External: `mafft`, `pal2nal.pl`

- [ ] **Step 1: 编写比对脚本**

```python
#!/usr/bin/env python3
"""
P2: Run MAFFT alignments for all paired sequences.
- Protein alignment (MAFFT --auto)
- Back-translate protein alignment to codon alignment via pal2nal
"""
import os
import sys
import csv
import subprocess
import yaml
from concurrent.futures import ProcessPoolExecutor, as_completed

def load_config():
    with open("config.yaml") as f:
        return yaml.safe_load(f)

def get_paired_genes(paired_dir):
    """Get list of gene names from paired sequence files."""
    genes = []
    for fname in sorted(os.listdir(paired_dir)):
        if fname.endswith("_hsa_ptr.fa"):
            gene = fname.replace("_hsa_ptr.fa", "")
            genes.append(gene)
    return genes

def align_one_gene(gene, paired_dir, prot_align_dir, cdna_align_dir, mafft_path, pal2nal_path, mafft_opts_prot, mafft_opts_cdna):
    """Run MAFFT and pal2nal for one gene pair. Returns (gene, status, n_gaps, aln_length)."""
    paired_file = os.path.join(paired_dir, f"{gene}_hsa_ptr.fa")
    prot_out = os.path.join(prot_align_dir, f"{gene}.aln")
    cdna_out = os.path.join(cdna_align_dir, f"{gene}.aln")

    if not os.path.exists(paired_file):
        return (gene, "no_paired_file", 0, 0)

    try:
        # Step 1: Align protein sequences
        cmd_prot = [
            mafft_path, "--auto",
            "--out", prot_out,
            paired_file
        ]
        # Remove --out flag — MAFFT writes to stdout
        cmd_prot = [mafft_path] + mafft_opts_prot.split()
        with open(paired_file) as inf, open(prot_out, "w") as outf:
            result = subprocess.run(
                [mafft_path] + mafft_opts_prot.split(),
                stdin=inf, stdout=outf, stderr=subprocess.PIPE,
                text=True, timeout=120
            )
        if result.returncode != 0:
            return (gene, f"mafft_prot_error: {result.stderr[:100]}", 0, 0)

        # Step 2: Back-translate to codon alignment via pal2nal
        # pal2nal needs: protein alignment + unaligned cDNA sequences
        # It extracts CDS from the input sequences
        cmd_codon = [
            "perl", pal2nal_path, prot_out, paired_file,
            "-output", "fasta", "-nogap"
        ]
        with open(cdna_out, "w") as outf:
            result = subprocess.run(
                cmd_codon, stdout=outf, stderr=subprocess.PIPE,
                text=True, timeout=120
            )
        if result.returncode != 0:
            return (gene, f"pal2nal_error: {result.stderr[:100]}", 0, 0)

        # Compute alignment stats
        with open(prot_out) as f:
            lines = f.readlines()
        seqs = []
        current = ""
        for line in lines:
            if line.startswith(">"):
                if current:
                    seqs.append(current)
                current = ""
            else:
                current += line.strip()
        if current:
            seqs.append(current)

        if len(seqs) < 2 or not seqs[0] or not seqs[1]:
            return (gene, "empty_alignment", 0, 0)

        aln_len = max(len(s) for s in seqs)
        # Count gap positions in either sequence
        gaps = sum(1 for i in range(min(len(s) for s in seqs[:2]))
                   if seqs[0][i] == "-" or seqs[1][i] == "-")

        return (gene, "ok", gaps, aln_len)

    except subprocess.TimeoutExpired:
        return (gene, "timeout", 0, 0)
    except Exception as e:
        return (gene, f"error: {str(e)[:100]}", 0, 0)

def main():
    cfg = load_config()
    root = os.path.abspath(cfg["project"]["root"])
    os.chdir(root)

    paired_dir = "data/paired_sequences"
    prot_align_dir = "data/alignments/protein"
    cdna_align_dir = "data/alignments/cdna"
    result_dir = "data/results"
    os.makedirs(prot_align_dir, exist_ok=True)
    os.makedirs(cdna_align_dir, exist_ok=True)
    os.makedirs(result_dir, exist_ok=True)

    mafft_path = cfg["alignment"].get("mafft_path", "mafft")
    pal2nal_path = cfg["alignment"].get("pal2nal_path", "bin/pal2nal.pl")
    mafft_opts_prot = cfg["alignment"].get("mafft_opts_protein", "--auto")
    mafft_opts_cdna = cfg["alignment"].get("mafft_opts_cdna", "--globalpair --maxiterate 1000")

    genes = get_paired_genes(paired_dir)
    print(f"Found {len(genes)} gene pairs to align")

    results = []
    # Use ProcessPoolExecutor for parallel alignment (IO-bound, not CPU)
    # Limited to 4 workers to avoid overwhelming the system
    with ProcessPoolExecutor(max_workers=4) as executor:
        futures = {
            executor.submit(
                align_one_gene, g, paired_dir, prot_align_dir, cdna_align_dir,
                mafft_path, pal2nal_path, mafft_opts_prot, mafft_opts_cdna
            ): g for g in genes
        }
        for i, future in enumerate(as_completed(futures), 1):
            result = future.result()
            results.append(result)
            if i % 100 == 0:
                print(f"  Progress: {i}/{len(genes)}")

    # Write alignment summary
    out_csv = os.path.join(result_dir, "alignments_summary.csv")
    with open(out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["gene", "status", "gap_count", "alignment_length"])
        writer.writerows(results)

    ok_count = sum(1 for r in results if r[1] == "ok")
    print(f"\nDone. {ok_count}/{len(genes)} alignments successful")
    print(f"Results: {out_csv}")

if __name__ == "__main__":
    main()
```

- [ ] **Step 2: 运行比对**

```bash
cd /saturn/zhaoty/evo_project/evolutionary_analysis
python python/02_run_alignment.py
# 预期: ~90%+ success rate (closely related sequences)
wc -l data/results/alignments_summary.csv
ls data/alignments/protein/ | wc -l
ls data/alignments/cdna/ | wc -l
```

- [ ] **Step 3: 提交**

```bash
cd /saturn/zhaoty/evo_project
git add evolutionary_analysis/python/02_run_alignment.py \
       evolutionary_analysis/data/alignments/ \
       evolutionary_analysis/data/results/alignments_summary.csv
git commit -m "feat(p2): MAFFT protein alignment and pal2nal codon back-translation"
```

---

### Task 4 (P3): dN/dS 选择压力分析

**Files:**
- Create: `evolutionary_analysis/python/03_calculate_dnds.py`

**Interfaces:**
- Consumes: `data/alignments/cdna/*.aln` (密码子比对，FASTA 格式)
- Produces: `data/results/dn_ds_results.csv`
- External: `codeml` (PAML)

**注意：** PAML codeml 需要 PHYLIP 格式的输入。需要将 FASTA 格式的密码子比对转换为 PHYLIP。

- [ ] **Step 1: 编写 dN/dS 计算脚本**

```python
#!/usr/bin/env python3
"""
P3: Calculate dN/dS (omega) for each gene pair using PAML codeml.
- Converts codon alignment to PHYLIP format
- Runs codeml in pairwise mode (runmode = -2)
- Parses output to extract dN, dS, omega, likelihood
- Performs likelihood ratio test (H0: omega=1 vs H1: omega free)
"""
import os
import sys
import csv
import re
import shutil
import subprocess
import yaml
from concurrent.futures import ProcessPoolExecutor, as_completed
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment

def load_config():
    with open("config.yaml") as f:
        return yaml.safe_load(f)

def fasta_to_phylip(fasta_path, phylip_path):
    """Convert FASTA alignment to sequential PHYLIP format for PAML."""
    records = list(SeqIO.parse(fasta_path, "fasta"))
    if len(records) < 2:
        return False

    # Check alignment length is multiple of 3 (codon alignment)
    aln_len = len(records[0].seq)
    if aln_len % 3 != 0:
        return False

    # Truncate names to 10 chars (PHYLIP standard)
    names = []
    seqs = []
    for r in records:
        name = r.id[:10]
        # Ensure names are unique
        orig_name = name
        suffix = 1
        while name in names:
            name = f"{orig_name[:8]}{suffix}"
            suffix += 1
        names.append(name)
        seqs.append(str(r.seq))

    with open(phylip_path, "w") as f:
        f.write(f"  {len(records)} {aln_len}\n")
        for name, seq in zip(names, seqs):
            f.write(f"{name:<10} {seq}\n")
    return True

def write_codeml_ctl(ctl_path, phylip_path, out_path, run_mode="pairwise"):
    """Write codeml control file for pairwise dN/dS estimation."""
    # Create output directory if needed
    out_dir = os.path.dirname(out_path)
    os.makedirs(out_dir, exist_ok=True)

    with open(ctl_path, "w") as f:
        f.write(f"      seqfile = {phylip_path}\n")
        f.write(f"     outfile = {out_path}.mlc\n")
        f.write(f"       treefile = {phylip_path}.tree\n\n")
        f.write(f"      noisy = 3\n")
        f.write(f"      verbose = 1\n")
        f.write(f"      runmode = -2\n")  # pairwise comparison
        f.write(f"      seqtype = 1\n")   # codon-based
        f.write(f"    CodonFreq = 2\n")   # F3X4
        f.write(f"        clock = 0\n")
        f.write(f"       aaDist = 0\n")
        f.write(f"    model = 0\n")       # one omega for all branches
        f.write(f"      NSsites = 0\n")   # single omega
        f.write(f"        icode = 0\n")
        f.write(f"    fix_kappa = 0\n")
        f.write(f"      kappa = 2\n")
        f.write(f"    fix_omega = 0\n")
        f.write(f"      omega = 0.5\n")
        f.write(f"        ncatG = 3\n")

def write_null_ctl(ctl_path, phylip_path, out_path):
    """Write codeml control file with omega fixed to 1 (null model)."""
    out_dir = os.path.dirname(out_path)
    os.makedirs(out_dir, exist_ok=True)

    with open(ctl_path, "w") as f:
        f.write(f"      seqfile = {phylip_path}\n")
        f.write(f"     outfile = {out_path}.mlc_null\n")
        f.write(f"       treefile = {phylip_path}.tree\n\n")
        f.write(f"      noisy = 3\n")
        f.write(f"      verbose = 1\n")
        f.write(f"      runmode = -2\n")
        f.write(f"      seqtype = 1\n")
        f.write(f"    CodonFreq = 2\n")
        f.write(f"        clock = 0\n")
        f.write(f"       aaDist = 0\n")
        f.write(f"    model = 0\n")
        f.write(f"      NSsites = 0\n")
        f.write(f"        icode = 0\n")
        f.write(f"    fix_kappa = 0\n")
        f.write(f"      kappa = 2\n")
        f.write(f"    fix_omega = 1\n")   # fix omega=1 (null)
        f.write(f"      omega = 1\n")
        f.write(f"        ncatG = 3\n")

def write_pairwise_tree(tree_path, seq_names):
    """Write a simple NJ tree for pairwise comparison."""
    with open(tree_path, "w") as f:
        # Simple unrooted tree with two species and a branch
        # For pairwise, any tree works — codeml with runmode=-2 ignores it
        name1, name2 = seq_names[:2]
        f.write(f"({name1},{name2});\n")

def parse_mlc(mlc_path):
    """Parse codeml mlc output file to extract dN, dS, omega, lnL."""
    dN = dS = omega = lnL = None
    try:
        with open(mlc_path) as f:
            content = f.read()

        # Extract lnL
        lnl_match = re.search(r"lnL\(.+?\)\s*=\s*(-?\d+\.\d+)", content)
        if lnl_match:
            lnL = float(lnl_match.group(1))

        # Extract pairwise dN/dS — look for "pairwise" section
        # In runmode=-2 output, look for "dN/dS" or "omega" in the pairwise results
        pairwise_section = re.search(
            r"pairwise comparison.*?\n(.*?)(?=\n\s*\n|\Z)",
            content, re.DOTALL
        )
        if pairwise_section:
            pw = pairwise_section.group(1)
            # Extract dN, dS, omega values
            dn_match = re.search(r"dN\s*=\s*(\d+\.\d+)", pw)
            ds_match = re.search(r"dS\s*=\s*(\d+\.\d+)", pw)
            om_match = re.search(r"omega\s*=\s*(\d+\.\d+)|dN/dS\s*=\s*(\d+\.\d+)", pw)
            if dn_match: dN = float(dn_match.group(1))
            if ds_match: dS = float(ds_match.group(1))
            if om_match:
                omega = float(om_match.group(1) or om_match.group(2))

        # Alternative: parse from the "ml" section
        if dN is None:
            ml_section = re.search(r"ml\s+\(.+?\)(.*?)(?=\n\s*\n|\Z)", content, re.DOTALL)
            if ml_section:
                ml_text = ml_section.group(1)
                dn_match = re.search(r"dN\s*=\s*(\d+\.\d+)", ml_text)
                ds_match = re.search(r"dS\s*=\s*(\d+\.\d+)", ml_text)
                om_match = re.search(r"omega\s*=\s*(\d+\.\d+)", ml_text)
                if dn_match: dN = float(dn_match.group(1))
                if ds_match: dS = float(ds_match.group(1))
                if om_match: omega = float(om_match.group(1))

    except Exception as e:
        print(f"  Parse error for {mlc_path}: {e}")

    return dN, dS, omega, lnL

def run_codeml_for_gene(gene, cdna_align_dir, temp_dir, paml_path, result_dir):
    """Run codeml for one gene pair. Returns result tuple."""
    gene_result = {
        "gene": gene, "dN": "", "dS": "", "omega": "",
        "lnL_alt": "", "lnL_null": "", "p_value": "", "status": "ok"
    }

    cdna_aln = os.path.join(cdna_align_dir, f"{gene}.aln")
    if not os.path.exists(cdna_aln):
        gene_result["status"] = "no_alignment"
        return gene_result

    # Create temp directory for this gene
    gene_temp = os.path.join(temp_dir, gene)
    os.makedirs(gene_temp, exist_ok=True)

    phylip_path = os.path.join(gene_temp, f"{gene}.phy")
    ctl_alt = os.path.join(gene_temp, f"{gene}_alt.ctl")
    ctl_null = os.path.join(gene_temp, f"{gene}_null.ctl")
    tree_path = os.path.join(gene_temp, f"{gene}.tree")
    out_alt = os.path.join(gene_temp, gene)
    out_null = os.path.join(gene_temp, f"{gene}_null")

    try:
        # Convert to PHYLIP
        if not fasta_to_phylip(cdna_aln, phylip_path):
            gene_result["status"] = "phylip_conversion_failed"
            return gene_result

        # Read sequence names for tree file
        records = list(SeqIO.parse(cdna_aln, "fasta"))
        if len(records) < 2:
            gene_result["status"] = "too_few_sequences"
            return gene_result

        # Write tree file
        write_pairwise_tree(tree_path, [r.id[:10] for r in records])

        # Run alternative model (omega free)
        write_codeml_ctl(ctl_alt, phylip_path, out_alt, "pairwise")
        result_alt = subprocess.run(
            [paml_path, ctl_alt],
            cwd=gene_temp, capture_output=True, text=True, timeout=60
        )
        if result_alt.returncode != 0:
            gene_result["status"] = f"codeml_alt_failed"
            return gene_result

        # Run null model (omega=1)
        write_null_ctl(ctl_null, phylip_path, out_null)
        result_null = subprocess.run(
            [paml_path, ctl_null],
            cwd=gene_temp, capture_output=True, text=True, timeout=60
        )

        # Parse results
        dN, dS, omega, lnL_alt = parse_mlc(f"{out_alt}.mlc")
        _, _, _, lnL_null = parse_mlc(f"{out_null}.mlc_null")

        # Likelihood ratio test
        p_value = ""
        if lnL_alt is not None and lnL_null is not None:
            LRT = 2 * (lnL_alt - lnL_null)  # positive when alt is better
            # Chi-square test with 1 degree of freedom
            from scipy.stats import chi2
            p_value = 1 - chi2.cdf(max(0, LRT), 1)
            gene_result["p_value"] = f"{p_value:.6e}"

        gene_result["dN"] = f"{dN:.6f}" if dN is not None else ""
        gene_result["dS"] = f"{dS:.6f}" if dS is not None else ""
        gene_result["omega"] = f"{omega:.6f}" if omega is not None else ""
        gene_result["lnL_alt"] = f"{lnL_alt:.4f}" if lnL_alt is not None else ""
        gene_result["lnL_null"] = f"{lnL_null:.4f}" if lnL_null is not None else ""

    except subprocess.TimeoutExpired:
        gene_result["status"] = "timeout"
    except Exception as e:
        gene_result["status"] = f"error: {str(e)[:100]}"

    # Cleanup temp files (keep results)
    # shutil.rmtree(gene_temp, ignore_errors=True)

    return gene_result

def main():
    cfg = load_config()
    root = os.path.abspath(cfg["project"]["root"])
    os.chdir(root)

    cdna_align_dir = "data/alignments/cdna"
    temp_dir = cfg["dnds"].get("temp_dir", "data/tmp_codeml")
    paml_path = cfg["dnds"].get("paml_path", "codeml")
    result_dir = "data/results"
    os.makedirs(temp_dir, exist_ok=True)

    # Get all aligned genes from cdna directory
    genes = sorted(set(f.replace(".aln", "") for f in os.listdir(cdna_align_dir) if f.endswith(".aln")))
    print(f"Found {len(genes)} codon alignments for dN/dS analysis")

    results = []
    # Process sequentially to avoid PAML file conflicts
    for i, gene in enumerate(genes, 1):
        result = run_codeml_for_gene(gene, cdna_align_dir, temp_dir, paml_path, result_dir)
        results.append(result)
        if i % 100 == 0:
            ok = sum(1 for r in results if r["status"] == "ok")
            print(f"  Progress: {i}/{len(genes)} ({ok} ok)")

    # Write results
    out_csv = os.path.join(result_dir, "dn_ds_results.csv")
    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "gene", "dN", "dS", "omega", "lnL_alt", "lnL_null", "p_value", "status"
        ])
        writer.writeheader()
        writer.writerows(results)

    ok_count = sum(1 for r in results if r["status"] == "ok")
    print(f"\nDone. {ok_count}/{len(genes)} dN/dS calculations successful")
    print(f"Results: {out_csv}")

if __name__ == "__main__":
    main()
```

- [ ] **Step 2: 安装 scipy（用于似然比检验）并运行**

```bash
pip install scipy

cd /saturn/zhaoty/evo_project/evolutionary_analysis
python python/03_calculate_dnds.py
# 预期: ~1000+ successful dN/dS calculations
wc -l data/results/dn_ds_results.csv
head -10 data/results/dn_ds_results.csv
```

- [ ] **Step 3: 提交**

```bash
cd /saturn/zhaoty/evo_project
git add evolutionary_analysis/python/03_calculate_dnds.py \
       evolutionary_analysis/data/results/dn_ds_results.csv
git commit -m "feat(p3): pairwise dN/dS calculation with PAML codeml"
```

---

### Task 5 (P4): 保守性打分

**Files:**
- Create: `evolutionary_analysis/python/04_compute_conservation.py`

**Interfaces:**
- Consumes: `data/alignments/protein/*.aln`, `data/results/dn_ds_results.csv`
- Produces: `data/results/conservation_scores.csv`, `data/results/sitewise_conservation/{gene}.csv`

- [ ] **Step 1: 编写保守性打分脚本**

```python
#!/usr/bin/env python3
"""
P4: Compute conservation scores for each aligned gene pair.
- Percent identity (protein and inferred cDNA)
- Site-wise Shannon entropy
- Sliding-window conservation (10aa window)
- BLOSUM62-based substitution score
"""
import os
import sys
import csv
import math
import yaml
from collections import Counter
from Bio import SeqIO
from Bio.SubsMat.MatrixInfo import blosum62

def load_config():
    with open("config.yaml") as f:
        return yaml.safe_load(f)

def parse_alignment(aln_path):
    """Parse a FASTA alignment file. Returns sequence records."""
    return list(SeqIO.parse(aln_path, "fasta"))

def calc_percent_identity(seq1, seq2):
    """Calculate percent identity (excluding gap positions)."""
    aligned = 0
    identical = 0
    for a, b in zip(seq1, seq2):
        if a != "-" and b != "-":
            aligned += 1
            if a.upper() == b.upper():
                identical += 1
    if aligned == 0:
        return 0.0
    return (identical / aligned) * 100

def calc_shannon_entropy(seq1, seq2):
    """Calculate site-wise Shannon entropy for two aligned sequences."""
    # For two sequences: entropy = 0 if identical, -2*(0.5*log2(0.5)) = 1 if different
    sites = []
    for a, b in zip(seq1, seq2):
        if a == "-" or b == "-":
            continue
        if a.upper() == b.upper():
            sites.append(0.0)  # entropy 0 when identical
        else:
            # Two residues, each freq = 0.5
            sites.append(-2 * (0.5 * math.log2(0.5)))  # = 1.0
    if not sites:
        return 0.0, 0.0, 0  # mean, min, count
    return sum(sites) / len(sites), min(sites), len(sites)

def calc_blosum62_score(seq1, seq2):
    """Calculate average BLOSUM62 substitution score."""
    scores = []
    for a, b in zip(seq1, seq2):
        if a == "-" or b == "-":
            continue
        key = (a.upper(), b.upper())
        # Try both orderings
        score = blosum62.get(key, blosum62.get((b.upper(), a.upper())))
        if score is not None:
            scores.append(score)
    if not scores:
        return 0.0
    return sum(scores) / len(scores)

def sliding_window_identity(seq1, seq2, window=10, step=1):
    """Calculate percent identity in sliding windows."""
    positions = []
    identities = []
    for i in range(0, min(len(seq1), len(seq2)) - window + 1, step):
        w1 = seq1[i:i+window]
        w2 = seq2[i:i+window]
        pid = calc_percent_identity(w1, w2)
        positions.append(i)
        identities.append(pid)
    return positions, identities

def calc_cdna_identity(prot_aln_path, cdna_aln_path):
    """Calculate cDNA-level identity from codon alignment data."""
    try:
        cdna_records = list(SeqIO.parse(cdna_aln_path, "fasta"))
        if len(cdna_records) >= 2:
            return calc_percent_identity(str(cdna_records[0].seq), str(cdna_records[1].seq))
    except Exception:
        pass
    return 0.0

def main():
    cfg = load_config()
    root = os.path.abspath(cfg["project"]["root"])
    os.chdir(root)

    prot_align_dir = "data/alignments/protein"
    cdna_align_dir = "data/alignments/cdna"
    result_dir = "data/results"
    sitewise_dir = "data/results/sitewise_conservation"
    os.makedirs(sitewise_dir, exist_ok=True)

    # Get genes with successful protein alignments
    genes = sorted(set(f.replace(".aln", "") for f in os.listdir(prot_align_dir) if f.endswith(".aln")))
    print(f"Computing conservation scores for {len(genes)} genes...")

    results = []
    for i, gene in enumerate(genes, 1):
        prot_aln = os.path.join(prot_align_dir, f"{gene}.aln")
        cdna_aln = os.path.join(cdna_align_dir, f"{gene}.aln")

        if not os.path.exists(prot_aln):
            continue

        records = parse_alignment(prot_aln)
        if len(records) < 2:
            continue

        seq_human = str(records[0].seq)
        seq_chimp = str(records[1].seq)

        # Core metrics
        prot_id = calc_percent_identity(seq_human, seq_chimp)
        cdna_id = calc_cdna_identity(prot_aln, cdna_aln) or calc_percent_identity(
            str(list(SeqIO.parse(cdna_aln, "fasta"))[0].seq) if os.path.exists(cdna_aln) else "",
            str(list(SeqIO.parse(cdna_aln, "fasta"))[1].seq) if os.path.exists(cdna_aln) else ""
        ) if os.path.exists(cdna_aln) else 0.0

        mean_entropy, min_entropy, n_sites = calc_shannon_entropy(seq_human, seq_chimp)
        blosum_score = calc_blosum62_score(seq_human, seq_chimp)

        # Sliding window analysis
        win_pos, win_ids = sliding_window_identity(seq_human, seq_chimp, window=10, step=5)
        min_window_id = min(win_ids) if win_ids else 0
        mean_window_id = sum(win_ids) / len(win_ids) if win_ids else 0

        # Write site-wise data
        site_file = os.path.join(sitewise_dir, f"{gene}.csv")
        with open(site_file, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["position", "human_aa", "chimp_aa", "identical", "blosum62"])
            for pos, (a, b) in enumerate(zip(seq_human, seq_chimp)):
                if a == "-" or b == "-":
                    continue
                key = (a.upper(), b.upper())
                b62 = blosum62.get(key, blosum62.get((b.upper(), a.upper()), ""))
                writer.writerow([pos + 1, a, b, "1" if a.upper() == b.upper() else "0", b62])

        results.append({
            "gene": gene,
            "prot_identity": f"{prot_id:.2f}",
            "cdna_identity": f"{cdna_id:.2f}",
            "mean_entropy": f"{mean_entropy:.4f}",
            "min_entropy": f"{min_entropy:.4f}",
            "mean_blosum62": f"{blosum_score:.4f}",
            "n_sites": n_sites,
            "min_window_identity": f"{min_window_id:.2f}",
            "mean_window_identity": f"{mean_window_id:.2f}",
        })

        if (i + 1) % 200 == 0:
            print(f"  Progress: {i+1}/{len(genes)}")

    # Write conservation scores
    out_csv = os.path.join(result_dir, "conservation_scores.csv")
    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "gene", "prot_identity", "cdna_identity", "mean_entropy", "min_entropy",
            "mean_blosum62", "n_sites", "min_window_identity", "mean_window_identity"
        ])
        writer.writeheader()
        writer.writerows(results)

    print(f"\nDone. {len(results)} genes scored")
    print(f"Conservation scores: {out_csv}")
    print(f"Site-wise details: {sitewise_dir}")

if __name__ == "__main__":
    main()
```

- [ ] **Step 2: 运行保守性打分**

```bash
cd /saturn/zhaoty/evo_project/evolutionary_analysis
python python/04_compute_conservation.py
# 预期: ~1000+ genes scored
wc -l data/results/conservation_scores.csv
head -10 data/results/conservation_scores.csv
```

- [ ] **Step 3: 提交**

```bash
cd /saturn/zhaoty/evo_project
git add evolutionary_analysis/python/04_compute_conservation.py \
       evolutionary_analysis/data/results/conservation_scores.csv \
       evolutionary_analysis/data/results/sitewise_conservation/
git commit -m "feat(p4): conservation scoring - identity, entropy, BLOSUM62, sliding window"
```

---

### Task 6 (R1): 数据汇总统计

**Files:**
- Create: `evolutionary_analysis/r/01_import_summary.R`

**Interfaces:**
- Consumes: `data/results/filtering_stats.csv`, `data/results/alignments_summary.csv`, `data/results/dn_ds_results.csv`, `data/results/conservation_scores.csv`
- Produces: 控制台报告 + `data/results/analysis_summary.csv`

- [ ] **Step 1: 编写汇总脚本**

```r
#!/usr/bin/env Rscript
# R1: Import all pipeline results and generate summary statistics

library(tidyverse)

setwd("/saturn/zhaoty/evo_project/evolutionary_analysis")

# Read all result files
filtering <- read_csv("data/results/filtering_stats.csv", show_col_types = FALSE)
alignments <- read_csv("data/results/alignments_summary.csv", show_col_types = FALSE)
dnds <- read_csv("data/results/dn_ds_results.csv", show_col_types = FALSE)
conservation <- read_csv("data/results/conservation_scores.csv", show_col_types = FALSE)

# ---- Filtering summary ----
cat("\n========================================\n")
cat("        PIPELINE SUMMARY REPORT\n")
cat("========================================\n\n")

cat("--- Filtering ---\n")
cat("Total one2one orthologs in mapping: ", nrow(filtering), "\n")
status_counts <- filtering %>% count(status)
print(status_counts)
cat("\n")

# ---- Alignment summary ----
cat("--- Alignment ---\n")
aln_ok <- alignments %>% filter(status == "ok")
cat("Successful alignments: ", nrow(aln_ok), " / ", nrow(alignments), "\n")
if (nrow(aln_ok) > 0) {
  cat("Alignment length - mean:",
      round(mean(aln_ok$alignment_length, na.rm = TRUE), 0),
      "bp, range:",
      min(aln_ok$alignment_length, na.rm = TRUE), "-",
      max(aln_ok$alignment_length, na.rm = TRUE), "\n")
}
cat("\n")

# ---- dN/dS summary ----
cat("--- dN/dS Analysis ---\n")
dnds_ok <- dnds %>% filter(status == "ok")
cat("Successful dN/dS calculations: ", nrow(dnds_ok), " / ", nrow(dnds), "\n")
if (nrow(dnds_ok) > 0) {
  dnds_num <- dnds_ok %>%
    mutate(omega = as.numeric(omega),
           dN = as.numeric(dN),
           dS = as.numeric(dS))

  cat("Omega (dN/dS) distribution:\n")
  cat("  Mean:", round(mean(dnds_num$omega, na.rm = TRUE), 4), "\n")
  cat("  Median:", round(median(dnds_num$omega, na.rm = TRUE), 4), "\n")
  cat("  Range:",
      round(min(dnds_num$omega, na.rm = TRUE), 4), "-",
      round(max(dnds_num$omega, na.rm = TRUE), 4), "\n")

  # Count genes under positive selection
  sig_pos <- dnds_num %>%
    filter(omega > 1, as.numeric(p_value) < 0.05)
  cat("\n  Genes under positive selection (omega > 1, p < 0.05): ",
      nrow(sig_pos), "\n")

  # Highly conserved genes (omega near 0)
  hc <- dnds_num %>% filter(omega < 0.1)
  cat("  Highly conserved genes (omega < 0.1): ", nrow(hc), "\n")
}
cat("\n")

# ---- Conservation summary ----
cat("--- Conservation Scores ---\n")
if (nrow(conservation) > 0) {
  cons_num <- conservation %>%
    mutate(prot_id = as.numeric(prot_identity),
           cdna_id = as.numeric(cdna_identity))

  cat("Protein identity - mean:",
      round(mean(cons_num$prot_id, na.rm = TRUE), 2), "%\n")
  cat("Protein identity < 90%:",
      sum(cons_num$prot_id < 90, na.rm = TRUE), "genes\n")
  cat("Protein identity = 100%:",
      sum(cons_num$prot_id == 100, na.rm = TRUE), "genes\n")
}
cat("\n")

# ---- Write consolidated summary ----
summary_df <- bind_rows(
  filtering %>% mutate(source = "filtering"),
  alignments %>% mutate(source = "alignment"),
  dnds %>% mutate(source = "dnds"),
  conservation %>% mutate(source = "conservation")
)
write_csv(summary_df, "data/results/analysis_summary.csv")
cat("Consolidated summary saved to: data/results/analysis_summary.csv\n")
cat("========================================\n")
```

- [ ] **Step 2: 运行汇总**

```bash
cd /saturn/zhaoty/evo_project/evolutionary_analysis
Rscript r/01_import_summary.R
# 预期: 完整管道统计输出
```

- [ ] **Step 3: 提交**

```bash
cd /saturn/zhaoty/evo_project
git add evolutionary_analysis/r/01_import_summary.R \
       evolutionary_analysis/data/results/analysis_summary.csv
git commit -m "feat(r1): pipeline summary statistics and reporting"
```

---

### Task 7 (R2): 保守性可视化

**Files:**
- Create: `evolutionary_analysis/r/02_visualize_conservation.R`

**Interfaces:**
- Consumes: `data/results/dn_ds_results.csv`, `data/results/conservation_scores.csv`
- Produces: `figures/*.pdf` — 5 种图表

- [ ] **Step 1: 编写可视化脚本**

```r
#!/usr/bin/env Rscript
# R2: Conservation and selection pressure visualization
# Colors aligned with enrichment_analysis.R

library(tidyverse)

setwd("/saturn/zhaoty/evo_project/evolutionary_analysis")
fig_dir <- "figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Color palette (matching existing enrichment style)
col_bp <- "#E64B35"
col_mf <- "#4DBBD5"
col_cc <- "#00A087"
col_kegg <- "#3C5488"

# ---- Load data ----
dnds <- read_csv("data/results/dn_ds_results.csv", show_col_types = FALSE)
conservation <- read_csv("data/results/conservation_scores.csv", show_col_types = FALSE)

# Merge datasets
data <- dnds %>%
  filter(status == "ok") %>%
  mutate(
    omega = as.numeric(omega),
    dN = as.numeric(dN),
    dS = as.numeric(dS),
    p_val = as.numeric(p_value),
    significant = ifelse(!is.na(p_val) & p_val < 0.05 & !is.na(omega) & omega > 1, "Positive selection", "Neutral/Purifying")
  ) %>%
  left_join(
    conservation %>% mutate(
      prot_id = as.numeric(prot_identity),
      cdna_id = as.numeric(cdna_identity)
    ),
    by = "gene"
  )

# ---- 1. Omega distribution histogram ----
p1 <- ggplot(data, aes(x = omega)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 60, fill = col_bp, alpha = 0.7, color = "white", size = 0.3) +
  geom_density(color = "#333333", linewidth = 0.8) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#333333", linewidth = 0.7) +
  annotate("text", x = 1.2, y = Inf, label = "ω = 1 (neutral)",
           hjust = 0, vjust = 2, size = 3.5, color = "#333333") +
  scale_x_continuous(limits = c(0, min(max(data$omega, na.rm = TRUE) * 1.1, 3)),
                     oob = scales::squish) +
  labs(title = "Distribution of dN/dS (ω) across RBP orthologs",
       subtitle = paste0("Human-Chimpanzee, N = ", nrow(data), " genes"),
       x = "ω (dN/dS)", y = "Density") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
        panel.grid.minor = element_blank())
ggsave(file.path(fig_dir, "omega_distribution.pdf"), p1, width = 8, height = 5)

# ---- 2. Protein identity vs Omega scatter ----
p2 <- ggplot(data, aes(x = prot_id, y = omega)) +
  geom_point(aes(color = significant), alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Positive selection" = "#E41A1C",
                                 "Neutral/Purifying" = "#377EB8")) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", size = 0.5) +
  geom_smooth(method = "loess", color = "#333333", se = TRUE, alpha = 0.15, size = 0.7) +
  labs(title = "Protein Identity vs Selection Pressure",
       x = "Protein sequence identity (%)",
       y = "ω (dN/dS)",
       color = "") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom",
        panel.grid.minor = element_blank())
ggsave(file.path(fig_dir, "identity_vs_omega.pdf"), p2, width = 8, height = 6)

# ---- 3. dN vs dS scatter ----
p3 <- ggplot(data, aes(x = dS, y = dN)) +
  geom_point(aes(color = omega), alpha = 0.7, size = 2) +
  scale_color_gradient(low = "#377EB8", high = "#E41A1C", name = "ω") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50", size = 0.5) +
  annotate("text", x = Inf, y = Inf, label = "dN = dS (ω=1)",
           hjust = 1.1, vjust = 1.5, size = 3, color = "grey40") +
  labs(title = "Non-synonymous vs Synonymous Substitution Rates",
       x = "dS (synonymous)", y = "dN (non-synonymous)") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.minor = element_blank())
ggsave(file.path(fig_dir, "dN_vs_dS.pdf"), p3, width = 7, height = 6)

# ---- 4. Conservation score histogram ----
p4 <- ggplot(data, aes(x = prot_id)) +
  geom_histogram(bins = 50, fill = col_cc, alpha = 0.7, color = "white", size = 0.3) +
  geom_vline(xintercept = median(data$prot_id, na.rm = TRUE),
             linetype = "dashed", color = "#333333", size = 0.7) +
  annotate("text", x = median(data$prot_id, na.rm = TRUE) + 1,
           y = Inf, label = paste0("Median: ", round(median(data$prot_id, na.rm = TRUE), 1), "%"),
           hjust = 0, vjust = 2, size = 3.5) +
  labs(title = "Protein Sequence Conservation between Human and Chimpanzee",
       x = "Sequence identity (%)", y = "Number of genes") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.minor = element_blank())
ggsave(file.path(fig_dir, "conservation_histogram.pdf"), p4, width = 8, height = 5)

# ---- 5. omega vs dS with significance ----
data <- data %>% mutate(
  selection_class = case_when(
    !is.na(significant) & significant == "Positive selection" & !is.na(omega) & omega > 5 ~ "Strong positive",
    !is.na(significant) & significant == "Positive selection" ~ "Positive selection",
    !is.na(omega) & omega < 0.1 ~ "Strongly conserved",
    TRUE ~ "Neutral/Mild constraint"
  )
)

p5 <- ggplot(data, aes(x = dS, y = omega)) +
  geom_point(aes(color = selection_class), alpha = 0.7, size = 2) +
  scale_color_manual(values = c(
    "Strong positive" = "#FF0000",
    "Positive selection" = "#E41A1C",
    "Strongly conserved" = "#00A087",
    "Neutral/Mild constraint" = "#377EB8"
  )) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", size = 0.5) +
  labs(title = "Selection Pressure Classification",
       x = "dS (synonymous rate)", y = "ω (dN/dS)",
       color = "Selection class") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom",
        panel.grid.minor = element_blank()) +
  guides(color = guide_legend(nrow = 2))
ggsave(file.path(fig_dir, "selection_classification.pdf"), p5, width = 8, height = 6)

cat("\nAll conservation figures saved to:", fig_dir, "\n")
cat("  - omega_distribution.pdf\n")
cat("  - identity_vs_omega.pdf\n")
cat("  - dN_vs_dS.pdf\n")
cat("  - conservation_histogram.pdf\n")
cat("  - selection_classification.pdf\n")
```

- [ ] **Step 2: 运行可视化**

```bash
cd /saturn/zhaoty/evo_project/evolutionary_analysis
Rscript r/02_visualize_conservation.R
ls figures/*.pdf
```

- [ ] **Step 3: 提交**

```bash
cd /saturn/zhaoty/evo_project
git add evolutionary_analysis/r/02_visualize_conservation.R \
       evolutionary_analysis/figures/
git commit -m "feat(r2): conservation and selection pressure visualization"
```

---

### Task 8 (R3): 扩展的富集分析

**Files:**
- Create: `evolutionary_analysis/r/03_enrichment_comparison.R`
- Modify: 无需修改 `enrichment_chimp_100pct.R` — 通过 source() 复用其函数

**Interfaces:**
- Consumes: `data/results/dn_ds_results.csv`, `data/results/conservation_scores.csv`
- Produces: `figures/enrichment_*.pdf`, 控制台富集结果

- [ ] **Step 1: 编写富集对比脚本**

```r
#!/usr/bin/env Rscript
# R3: GO/KEGG enrichment for conserved and positively-selected RBP gene sets
# Reuses the hypergeometric test functions from enrichment_chimp_100pct.R

library(org.Hs.eg.db)
library(GO.db)
library(KEGGREST)
library(ggplot2)
library(tidyverse)

setwd("/saturn/zhaoty/evo_project/evolutionary_analysis")
fig_dir <- "figures"

# ---- Load analysis results ----
dnds <- read_csv("data/results/dn_ds_results.csv", show_col_types = FALSE)
conservation <- read_csv("data/results/conservation_scores.csv", show_col_types = FALSE)

data <- dnds %>%
  filter(status == "ok") %>%
  mutate(
    omega = as.numeric(omega),
    p_val = as.numeric(p_value)
  ) %>%
  left_join(conservation %>%
              mutate(prot_id = as.numeric(prot_identity)),
            by = "gene")

# ---- Define gene sets ----
# Get gene symbols for each set
all_genes_data <- data

# Highly conserved: omega < 0.1 AND prot_identity > 99%
conserved_set <- all_genes_data %>%
  filter(omega < 0.1, prot_id > 99) %>%
  pull(gene)
cat("Highly conserved genes (ω < 0.1, identity > 99%): ", length(conserved_set), "\n")

# Positively selected: omega > 1 and p < 0.05
positive_set <- all_genes_data %>%
  filter(omega > 1, !is.na(p_val), p_val < 0.05) %>%
  pull(gene)
cat("Positively selected genes (ω > 1, p < 0.05): ", length(positive_set), "\n")

# ---- Map symbols to Entrez IDs ----
map_symbols <- function(symbols) {
  if (length(symbols) == 0) return(character(0))
  map <- tryCatch(
    select(org.Hs.eg.db, keys = symbols, keytype = "SYMBOL",
           columns = "ENTREZID"),
    error = function(e) NULL
  )
  if (is.null(map)) return(character(0))
  unique(map$ENTREZID[!is.na(map$ENTREZID)])
}

conserved_entrez <- map_symbols(conserved_set)
positive_entrez <- map_symbols(positive_set)
cat("Conserved mapped to Entrez:", length(conserved_entrez), "\n")
cat("Positive mapped to Entrez:", length(positive_entrez), "\n")

# ---- Reuse hypergeometric test from existing enrichment_analysis.R ----
source("../enrichment_analysis.R")

# Override the genes variable with our sets
# Extract the enrichment function from the sourced script
do_go_enrichment <- function(gene_vec, universe_vec, ontology = "BP") {
  all_go_info <- select(GO.db, keys = ls(GOTERM), columns = c("GOID", "ONTOLOGY", "TERM"))
  go_ids_ont <- all_go_info$GOID[all_go_info$ONTOLOGY == ontology]
  message("  Ontology ", ontology, " — ", length(go_ids_ont), " terms to test")

  go2gene <- tryCatch(
    as.list(org.Hs.egGO2ALLEGS),
    error = function(e) NULL
  )
  if (is.null(go2gene)) return(NULL)

  gene_set <- unique(gene_vec)
  n_univ   <- length(unique(universe_vec))
  k        <- length(gene_set)

  .enrich_one <- function(goid) {
    term_genes <- go2gene[[goid]]
    if (is.null(term_genes)) return(NULL)
    term_genes <- intersect(term_genes, universe_vec)
    M <- length(term_genes)
    if (M < 2) return(NULL)

    x <- length(intersect(gene_set, term_genes))
    if (x < 1) return(NULL)

    pval <- phyper(x - 1, M, n_univ - M, k, lower.tail = FALSE)

    data.frame(
      GOID        = goid,
      Ontology    = ontology,
      Description = Term(GOTERM[[goid]]),
      GeneRatio   = paste0(x, "/", k),
      BgRatio     = paste0(M, "/", n_univ),
      pvalue      = pval,
      Count       = x,
      stringsAsFactors = FALSE
    )
  }

  res_list <- lapply(go_ids_ont, .enrich_one)
  res <- do.call(rbind, res_list[!sapply(res_list, is.null)])
  if (is.null(res) || nrow(res) == 0) return(NULL)

  res$p.adjust <- p.adjust(res$pvalue, method = "BH")
  res <- res[order(res$p.adjust), ]
  rownames(res) <- NULL
  res
}

universe <- keys(org.Hs.eg.db, keytype = "ENTREZID")

# Run enrichment for conserved set
if (length(conserved_entrez) >= 3) {
  cat("\n=== Enrichment for HIGHLY CONSERVED genes ===\n")
  cons_bp <- do_go_enrichment(conserved_entrez, universe, "BP")
  cons_mf <- do_go_enrichment(conserved_entrez, universe, "MF")
  cons_cc <- do_go_enrichment(conserved_entrez, universe, "CC")

  # Filter significant and plot
  sig_cons_bp <- cons_bp %>% filter(p.adjust < 0.05)
  cat("Conserved BP terms:", nrow(sig_cons_bp), "\n")

  if (!is.null(sig_cons_bp) && nrow(sig_cons_bp) > 0) {
    top_bp <- head(sig_cons_bp, 15)
    top_bp$Description <- factor(top_bp$Description, levels = rev(top_bp$Description))
    p <- ggplot(top_bp, aes(x = -log10(p.adjust), y = Description)) +
      geom_bar(stat = "identity", fill = "#00A087", alpha = 0.8) +
      labs(title = "GO BP: Highly Conserved RBPs (ω < 0.1, identity > 99%)",
           x = "-log10(p.adjust)", y = "") +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11))
    ggsave(file.path(fig_dir, "enrichment_conserved_BP.pdf"), p, width = 10, height = 6)
  }
}

# Run enrichment for positively selected set
if (length(positive_entrez) >= 3) {
  cat("\n=== Enrichment for POSITIVELY SELECTED genes ===\n")
  pos_bp <- do_go_enrichment(positive_entrez, universe, "BP")
  pos_mf <- do_go_enrichment(positive_entrez, universe, "MF")

  sig_pos_bp <- pos_bp %>% filter(p.adjust < 0.05)
  cat("Positively selected BP terms:", nrow(sig_pos_bp), "\n")

  if (!is.null(sig_pos_bp) && nrow(sig_pos_bp) > 0) {
    top_bp <- head(sig_pos_bp, 15)
    top_bp$Description <- factor(top_bp$Description, levels = rev(top_bp$Description))
    p <- ggplot(top_bp, aes(x = -log10(p.adjust), y = Description)) +
      geom_bar(stat = "identity", fill = "#E41A1C", alpha = 0.8) +
      labs(title = "GO BP: Positively Selected RBPs (ω > 1, p < 0.05)",
           x = "-log10(p.adjust)", y = "") +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11))
    ggsave(file.path(fig_dir, "enrichment_positive_BP.pdf"), p, width = 10, height = 6)
  }
}

cat("\nEnrichment comparison complete.\n")
```

- [ ] **Step 2: 运行富集分析**

```bash
cd /saturn/zhaoty/evo_project/evolutionary_analysis
Rscript r/03_enrichment_comparison.R
# 预期: 高保守基因和正选择基因的 GO 富集结果 + 图表
ls figures/enrichment_*.pdf
```

- [ ] **Step 3: 提交**

```bash
cd /saturn/zhaoty/evo_project
git add evolutionary_analysis/r/03_enrichment_comparison.R \
       evolutionary_analysis/figures/enrichment_*.pdf
git commit -m "feat(r3): GO enrichment for conserved and positively-selected RBP sets"
```

---

### Task 9 (R4): 整合报告（可选）

**Files:**
- Create: `evolutionary_analysis/r/04_render_report.R`
- Create: `evolutionary_analysis/report/analysis_report.Rmd`

**Interfaces:**
- Consumes: 所有 CSV + figures
- Produces: `report/analysis_report.html`

- [ ] **Step 1: 创建 R Markdown 报告**

```rmd
---
title: "RBP Evolutionary Conservation Analysis: Human vs Chimpanzee"
date: "`r Sys.Date()`"
output:
  html_document:
    toc: true
    toc_float: true
    theme: flatly
    code_folding: hide
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)
library(tidyverse)
setwd("/saturn/zhaoty/evo_project/evolutionary_analysis")
```

## Pipeline Summary

```{r summary}
dnds <- read_csv("data/results/dn_ds_results.csv", show_col_types = FALSE)
conservation <- read_csv("data/results/conservation_scores.csv", show_col_types = FALSE)

dnds_ok <- dnds %>% filter(status == "ok")
n_genes <- nrow(dnds_ok)
n_conserved <- dnds_ok %>% filter(as.numeric(omega) < 0.1) %>% nrow()
n_positive <- dnds_ok %>% filter(as.numeric(omega) > 1,
                                  as.numeric(p_value) < 0.05) %>% nrow()
```

- **Total ortholog pairs analyzed:** `r n_genes`
- **Highly conserved (ω < 0.1):** `r n_conserved`
- **Positive selection signals (ω > 1, p < 0.05):** `r n_positive`

## Selection Pressure Distribution

![](../figures/omega_distribution.pdf)

## Conservation vs Selection

![](../figures/identity_vs_omega.pdf)

## Substitution Rates

![](../figures/dN_vs_dS.pdf)

## Conservation Histogram

![](../figures/conservation_histogram.pdf)

## Selection Classification

![](../figures/selection_classification.pdf)

## Enrichment Analysis

### Highly Conserved RBPs
![](../figures/enrichment_conserved_BP.pdf)

### Positively Selected RBPs
![](../figures/enrichment_positive_BP.pdf)

## Top Genes of Interest

```{r top_genes}
data <- dnds_ok %>%
  mutate(omega = as.numeric(omega),
         dN = as.numeric(dN),
         dS = as.numeric(dS),
         p_val = as.numeric(p_value))

# Top 10 most conserved
cat("### Most Conserved (ω near 0)\n")
data %>% arrange(omega) %>% head(10) %>%
  select(gene, dN, dS, omega) %>% knitr::kable(digits = 4)

# Top 10 positive selection
cat("\n### Top Positive Selection Candidates\n")
data %>% filter(!is.na(p_val), p_val < 0.05) %>% arrange(desc(omega)) %>% head(10) %>%
  select(gene, dN, dS, omega, p_val) %>% knitr::kable(digits = 4)
```
```

- [ ] **Step 2: 渲染报告**

```bash
cd /saturn/zhaoty/evo_project/evolutionary_analysis
Rscript -e "rmarkdown::render('report/analysis_report.Rmd', output_file='analysis_report.html')"
```

- [ ] **Step 3: 提交**

```bash
cd /saturn/zhaoty/evo_project
git add evolutionary_analysis/r/04_render_report.R \
       evolutionary_analysis/report/
git commit -m "feat(r4): integrated R Markdown analysis report"
```

---

## 执行顺序

```
Task 0: 环境准备
    │
Task 1: 下载 Ensembl 转录本映射
    │
Task 2: P1 — 序列过滤配对 (01_filter_and_pair.py)
    │
    ▼
Task 3: P2 — MAFFT 比对 (02_run_alignment.py)
    │
    ├──────────────┐
    ▼              ▼
Task 4: P3      Task 5: P4
(dN/dS,       (保守性打分,
03_calc_dnds)  04_compute_conservation)
    │              │
    └──────┬───────┘
           ▼
Task 6: R1 — 数据汇总 (01_import_summary.R)
    │
    ▼
Task 7: R2 — 可视化 (02_visualize_conservation.R)
    │
    ▼
Task 8: R3 — 富集比较 (03_enrichment_comparison.R)
    │
    ▼
Task 9: R4 — 整合报告 (可选)
```

快速启动（全部依次执行）：
```bash
cd /saturn/zhaoty/evo_project/evolutionary_analysis
python python/01_filter_and_pair.py && \
python python/02_run_alignment.py && \
python python/03_calculate_dnds.py && \
python python/04_compute_conservation.py && \
Rscript r/01_import_summary.R && \
Rscript r/02_visualize_conservation.R && \
Rscript r/03_enrichment_comparison.R
```
