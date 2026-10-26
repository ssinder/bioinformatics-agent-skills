---
name: bulk-rnaseq
description: "This skill should be used when running bulk RNA-seq differential expression analysis: FastQC/trimming, Salmon quantification, counts matrix construction, DESeq2 statistics, GO/KEGG enrichment and publication figures. 用于 bulk 转录组（RNA-seq）差异表达分析全流程。"
license: MIT
compatibility: "Python 3.9+; numpy, scipy, pandas required. Optional: matplotlib (figures), Rscript + Bioconductor DESeq2 (r backend), fastqc/fastp/salmon (external tools)"
metadata:
  version: "1.0"
  last-reviewed: "2026-10-25"
  skill-author: sinder
---

## Overview

本技能把 bulk RNA-seq 的原始 FASTQ 一路处理到可发表的差异表达结论：先做质控与接头去除，再用 Salmon 做转录本定量并汇总成基因水平整数 counts，随后运行 DESeq2 风格负二项 Wald 检验得到差异基因列表，接着做 GO/KEGG 超几何富集，最后输出火山图、MA 图、PCA 图、热图与文字报告。

差异表达统计提供 **双后端**：`--backend python` 使用内置的纯 NumPy/SciPy 实现（size factor 归一化、Cox-Reid 离散度估计、趋势拟合、MAP 收缩、Wald 检验、BH 校正），`--backend r` 调用 Rscript 上的 Bioconductor DESeq2。两者输出同一套列名的结果表，可直接互换对比。

## When to Use This Skill

- 用户要做 bulk 转录组差异表达分析（关键词：RNA-seq、bulk RNA-seq、转录组、差异表达、DEG、DESeq2、Salmon、counts matrix）
- 用户已有 Salmon/STAR 的定量结果，需要构建 counts 矩阵并做差异表达
- 用户需要 GO / KEGG 富集分析（关键词：GO 富集、KEGG 富集、pathway enrichment、功能富集）
- 用户需要 RNA-seq 的标准发表级图表（火山图、MA 图、PCA、热图）
- 用户想在没有 R / Bioconductor 的环境里得到可靠的 DESeq2 风格结果

不适用：单细胞数据（用 `single-cell-rnaseq`）、仅做原始 reads 可视化、蛋白组与代谢组分析。

> **上游衔接**：如果还没有 FASTQ 文件，可先用 [`ena-download`](../ena-download/) 技能按
> SRA/ENA accession 下载原始测序数据（支持断点续传与 MD5 校验），再交给本技能做定量与差异分析。

## Quick Start

```bash
# 完整流程（需要 fastqc / fastp / salmon 在 PATH 中）
export PYTHONPATH="$PWD/packages"
python skills/bulk-rnaseq/scripts/run_pipeline.py \
    --input  /data/fastq \
    --output /results/rnaseq \
    --transcriptome /ref/transcriptome.fa \
    --tx2gene /ref/tx2gene.tsv \
    --metadata /data/metadata.csv \
    --go-annotation /ref/go_annotation.csv \
    --group-a treated --group-b control

# 跳过外部工具，直接从已有 Salmon 结果开始
python skills/bulk-rnaseq/scripts/run_pipeline.py \
    --input /data/fastq --output /results/rnaseq \
    --skip-qc --skip-trimming --skip-quant \
    --quant-dir /results/rnaseq/03_salmon \
    --tx2gene /ref/tx2gene.tsv \
    --metadata /data/metadata.csv \
    --group-a treated --group-b control --backend python

# 单独构建 counts 矩阵
python skills/bulk-rnaseq/scripts/build_counts_matrix.py \
    --quant-dir /results/rnaseq/03_salmon \
    --tx2gene /ref/tx2gene.tsv \
    --output /results/counts
```

## Workflow

1. **FastQC 质控** — 对每个 `.fastq.gz` 生成 HTML/ZIP 报告，可选 MultiQC 汇总。关注 per-base quality、adapter content、GC 分布。
2. **fastp 过滤** — 质量过滤（默认 Phred Q20）、PE 自动检测接头、最短长度 36 bp；同时输出 JSON/HTML 报告。
3. **Salmon 定量** — `salmon index` 构建转录组索引，`salmon quant` 带 `--gcBias --seqBias --validateMappings` 输出转录本级 TPM / NumReads / Length。
4. **counts 矩阵构建** — 用 `tx2gene.tsv` 把转录本汇总到基因；counts 取整（NumReads 四舍五入），TPM/FPKM 同步输出；生成元数据模板。
5. **DESeq2 差异表达** — size factor 归一化 → 基因离散度估计 → 趋势拟合 → MAP 离散度 → Wald 检验 → BH FDR。按 `padj` 与 `|log2FC|` 双阈值标注 Up / Down / NS。
6. **GO/KEGG 富集** — 单侧超几何检验 + BH FDR，按 BP / CC / MF 分别输出。
7. **可视化** — 火山图、MA 图、富集条形图、样本 PCA、Top DEG 热图（PDF + PNG，300 DPI）。
8. **报告** — 汇总实验设计、DEG 列表、富集结果、方法说明与产出清单。

## Parameters

| 参数 | 默认 | 说明 |
|---|---|---|
| `--input` | 必填 | 原始 FASTQ 目录 |
| `--output` | `<results>` | 输出根目录，下设 8 个步骤子目录 |
| `--transcriptome` | 无 | 转录组 FASTA，用于构建 Salmon 索引 |
| `--tx2gene` | 无 | 转录本-基因映射 TSV（列 `transcript_id` / `gene_id`） |
| `--metadata` | 无 | 样本元数据 CSV，需含 `sample` 与 `condition` 列 |
| `--go-annotation` | 无 | 本地 GO/KEGG 注释 CSV |
| `--backend` | `python` | 差异表达后端：`python`（内置）或 `r`（Rscript + DESeq2） |
| `--group-a` / `--group-b` | `treated` / `control` | 分子组与分母组标签 |
| `--design` | `~condition` | R 后端的设计公式 |
| `--alpha` | `0.05` | FDR 阈值 |
| `--lfc-threshold` | `1.0` | `|log2FC|` 阈值 |
| `--min-counts` | `10` | 低表达基因过滤阈值 |
| `--min-gene-set` / `--max-gene-set` | `5` / `500` | 富集分析的基因集大小范围 |
| `--threads` | `8` | 外部工具并行线程数 |
| `--skip-qc` / `--skip-trimming` / `--skip-quant` | 关 | 跳过对应步骤，复用已有中间结果 |

## Outputs

| 文件 | 说明 |
|---|---|
| `04_counts/counts.csv` | 基因水平整数 counts（行=基因，列=样本） |
| `04_counts/tpm.csv` / `fpkm.csv` | 基因水平 TPM / FPKM |
| `04_counts/metadata_template.csv` | 样本元数据模板（需人工确认 condition） |
| `05_deseq2/deseq2_results_full.csv` | 全基因结果：`baseMean`/`log2FoldChange`/`lfcSE`/`stat`/`pvalue`/`padj`/`significant`/`direction` |
| `05_deseq2/deseq2_results_significant.csv` | 显著差异基因子集 |
| `05_deseq2/size_factors.csv` | 每个样本的 size factor |
| `05_deseq2/dispersions.csv` | 每个基因的 MAP 离散度 |
| `05_deseq2/dispersion_trend.json` | 离散度趋势参数 `a0`/`a1`/`sigma2` |
| `06_enrichment/enrichment_{BP,CC,MF}.csv` | 各本体的富集结果 |
| `06_enrichment/enrichment_summary_all.csv` | 合并富集结果 |
| `07_figures/Figure{1..5}_*.{pdf,png}` | 火山图、MA 图、富集条形图、PCA、热图 |
| `08_report/analysis_report.txt` | 完整分析报告 |
| `run_manifest.json` | 运行元数据（输入路径、后端、对比方向） |

## Validation

- **counts 矩阵**：检查 `counts.csv` 各样本列和是否接近（差异过大说明测序深度差异大，需检查 size factor）；确认 gene ID 与注释文件 ID 空间一致。
- **size factor**：应集中在 1 附近（0.5–2）。出现极端值通常意味着某样本建库或测序异常。
- **火山图**：显著点应呈对称的双翼分布；MA 图应大致平坦。大量点集中在 `log2FC = 0` 附近是正常现象。
- **p 值分布**：期望绝大多数基因 `padj > 0.05`；若几乎全部显著，多半是分组标签写反或存在批次混杂。
- **离散度**：`dispersions.csv` 中落在下限（1e-8 附近）说明该基因数据量不足，其 p 值不可用。
- **富集**：检查 `Overlap_Genes` 与 `Term_Genes` 的比例与背景基因数；`min_gene_set` 太小会得到大量无意义的极小 term。
- **交叉验证**：同一数据用 `--backend python` 与 `--backend r` 各跑一次，比较 log2FC 排序与显著基因集合是否基本一致。

## References

- [references/methodology.md](references/methodology.md) — 各步骤方法学原理、统计推导与参数选择依据
- [references/troubleshooting.md](references/troubleshooting.md) — 常见报错、结果异常与排错手册
