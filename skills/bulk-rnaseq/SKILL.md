---
name: bulk-rnaseq
description: "Use when running bulk RNA-seq differential expression analysis: FastQC/trimming, Salmon quantification, counts matrix construction, DESeq2 statistics, GO/KEGG enrichment and publication figures. 用于 bulk 转录组（RNA-seq）差异表达分析全流程。"
license: MIT
compatibility: "Python 3.9+; numpy, scipy, pandas required. Optional: matplotlib (figures), Rscript + Bioconductor DESeq2 (r backend), fastqc/fastp/salmon (external tools)"
metadata:
  version: "1.0"
  last-reviewed: "2026-10-25"
  skill-author: sinder
  edam-topic: "topic_3170, topic_0203"
  edam-operation: "operation_3800, operation_3223, operation_3501"
  edam-input-data: "data_0928, data_3754"
  edam-output-data: "data_0928, data_3754"
  edam-input-format: "format_1930, format_3475"
  edam-output-format: "format_3752, format_3464"


---

## Overview

Takes bulk RNA-seq FASTQ files all the way to publishable differential expression results.

Differential expression statistics ship with **two backends**: `--backend python` uses the built-in pure NumPy/SciPy implementation, `--backend r` calls Bioconductor DESeq2 via Rscript. Both emit tables with identical column names, so they compare directly.

## When to Use This Skill

- Bulk transcriptome differential expression (RNA-seq, transcriptome, DEG, DESeq2, Salmon, counts matrix)
- Counts matrix from existing Salmon/STAR quantification, then DE testing
- GO or KEGG enrichment analysis (pathway / functional enrichment)
- Publication-grade RNA-seq figures (volcano, MA, PCA, heatmap)
- DESeq2-style results without R or Bioconductor

Not applicable to: single-cell data (use `single-cell-rnaseq`), raw-read visualisation, proteomics, metabolomics.

> **Upstream chaining**: without FASTQ yet, use [`ena-download`](../ena-download/) for SRA/ENA reads.

## Quick Start

```bash
export PYTHONPATH="$PWD/packages"   # needs fastqc / fastp / salmon on PATH
python skills/bulk-rnaseq/scripts/run_pipeline.py --input /data/fastq --output /results/rnaseq --transcriptome /ref/transcriptome.fa --tx2gene /ref/tx2gene.tsv --metadata /data/metadata.csv --go-annotation /ref/go_annotation.csv --group-a treated --group-b control
python skills/bulk-rnaseq/scripts/run_pipeline.py --input /data/fastq --output /results/rnaseq --skip-qc --skip-trimming --skip-quant --quant-dir /results/rnaseq/03_salmon --tx2gene /ref/tx2gene.tsv --metadata /data/metadata.csv --group-a treated --group-b control --backend python
python skills/bulk-rnaseq/scripts/build_counts_matrix.py --quant-dir /results/rnaseq/03_salmon --tx2gene /ref/tx2gene.tsv --output /results/counts
```


## Workflow

1. **FastQC** — report per `.fastq.gz`, optional MultiQC summary.
2. **fastp trimming** — quality filtering (Phred Q20 default), paired-end adapter auto-detection, 36 bp minimum length; JSON/HTML reports.
3. **Salmon quantification** — `salmon index` builds the index; `salmon quant` runs with `--gcBias --seqBias --validateMappings`, emitting transcript-level TPM / NumReads / Length.
4. **Counts matrix** — `tx2gene.tsv` maps transcripts to genes; counts rounded from NumReads, TPM/FPKM alongside, plus a metadata template.
5. **DESeq2-style DE** — size factor normalisation → gene-wise dispersion (Cox-Reid) → trend fitting → MAP shrinkage → Wald test → BH FDR; genes labelled Up / Down / NS by `padj`/`|log2FC|`.
6. **GO/KEGG enrichment** — one-sided hypergeometric test with BH FDR, reported separately for BP / CC / MF.
7. **Visualisation** — volcano, MA, enrichment bars, sample PCA, top-DEG heatmap (PDF + PNG, 300 DPI).
8. **Report** — design, DEG list, enrichment results, method notes, manifest.

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `--input` | required | Raw FASTQ directory |
| `--output` | `<results>` | Output root |
| `--transcriptome` | none | Transcriptome FASTA for the Salmon index |
| `--tx2gene` | none | Transcript-to-gene TSV (`transcript_id`/`gene_id`) |
| `--metadata` | none | Sample metadata CSV (`sample`,`condition`) |
| `--go-annotation` | none | Local GO/KEGG annotation CSV |
| `--backend` | `python` | `python` (built-in) or `r` (Rscript + DESeq2) |
| `--group-a` / `--group-b` | `treated` / `control` | Numerator / denominator groups |
| `--design` | `~condition` | Design formula (R backend) |
| `--alpha` | `0.05` | FDR threshold |
| `--lfc-threshold` | `1.0` | `|log2FC|` threshold |
| `--min-counts` | `10` | Low-expression gene filter |
| `--min-gene-set` / `--max-gene-set` | `5` / `500` | Gene-set size range for enrichment |
| `--threads` | `8` | Threads for external tools |
| `--skip-qc` / `--skip-trimming` / `--skip-quant` | off | Skip a step, reuse output |

## Outputs

- `04_counts/counts.csv` (gene-level integer counts), `04_counts/tpm.csv`, `fpkm.csv`, `04_counts/metadata_template.csv`
- `05_deseq2/deseq2_results_full.csv` (`baseMean`/`log2FoldChange`/`lfcSE`/`stat`/`pvalue`/`padj`/`significant`/`direction`), `05_deseq2/deseq2_results_significant.csv`, `05_deseq2/size_factors.csv`, `05_deseq2/dispersions.csv`, `05_deseq2/dispersion_trend.json`
- `06_enrichment/enrichment_{BP,CC,MF}.csv`, `06_enrichment/enrichment_summary_all.csv`
- `07_figures/Figure{1..5}_*.{pdf,png}`
- `08_report/analysis_report.txt`; root `run_manifest.json`

## Validation

- **Counts matrix** — column sums similar across samples; gene IDs share the annotation namespace.
- **Size factors** — cluster near 1 (0.5–2); extremes mean a library prep or sequencing problem.
- **Volcano / MA plots** — significant points form roughly symmetric wings, MA plot fairly flat; many points near `log2FC = 0` is normal.
- **p-value distribution** — most genes should have `padj > 0.05`; if nearly all are significant, group labels are swapped or batch effects confound it.
- **Dispersions** — genes at the lower bound (near 1e-8) have too few counts for a trustworthy p-value.
- **Enrichment** — check the `Overlap_Genes` / `Term_Genes` ratio and background gene count; too small a `min_gene_set` yields many tiny terms.
- **Cross-validation** — run the same dataset with both backends, then compare the log2FC ranking and significant gene sets; they should agree closely.

## References

- [references/methodology.md](references/methodology.md) — statistical rationale per step, derivations, parameter choices
- [references/troubleshooting.md](references/troubleshooting.md) — common errors, unexpected results, debugging
