---
name: single-cell-rnaseq
description: "Use when running single-cell RNA-seq analysis: reading 10x/AnnData matrices, QC filtering, normalization, HVG selection, PCA, clustering, UMAP/t-SNE and cell type annotation. 用于单细胞转录组（scRNA-seq）质控、降维、聚类与细胞注释全流程。"
license: MIT
compatibility: "Python 3.9+; numpy, scipy, pandas, scikit-learn required. Optional: matplotlib (figures), scanpy+anndata (scanpy backend), leidenalg+igraph (Leiden clustering), umap-learn (official UMAP)"
metadata:
  version: "1.0"
  last-reviewed: "2026-10-25"
  skill-author: sinder
  edam-topic: "topic_3170, topic_0203"
  edam-operation: "operation_0313, operation_3432, operation_3680"
  edam-input-data: "data_3112, data_0849"
  edam-output-data: "data_3112, data_0872"
  edam-input-format: "format_1475, format_3752"
  edam-output-format: "format_3752, format_3464"


---
## Overview

Turns a raw single-cell count matrix into cell populations and cell type labels: auto-detects the four input formats (10x HDF5 / MTX / dense CSV / h5ad), QC-filters on mitochondrial %, gene count and UMI count, then runs the normalisation → HVG → PCA → clustering → embedding → marker pipeline under Workflow. `--backend scanpy` calls official scanpy; `--backend numpy` uses a built-in NumPy/SciPy/scikit-learn fallback (own Louvain/UMAP spectral init + force optimisation) for equivalent results without scanpy.

## When to Use This Skill

- Single-cell transcriptome analysis (keywords: single-cell, scRNA-seq, cell atlas, annotation, UMAP, clustering)
- Raw 10x counts (`matrix.h5` / `filtered_feature_bc_matrix/`) needing standard QC + clustering
- Judging cell quality (keywords: doublet, mitochondrial %, QC, filtering threshold)
- Marker genes per population, to name cell types
- No scanpy, but reproducible reduction + clustering still needed
- Existing clustering results needing t-SNE / heatmaps / a report

Not applicable: bulk RNA-seq (`bulk-rnaseq`), spatial transcriptomics (CST/Slide-seq), single-cell ATAC, trajectory inference.

## Quick Start

```bash
export PYTHONPATH="$PWD/packages"
python skills/single-cell-rnaseq/scripts/run_pipeline.py --input /data/filtered_feature_bc_matrix --output /results/scrna --marker-file /ref/cell_marker_reference.csv
python skills/single-cell-rnaseq/scripts/run_pipeline.py --input /data/raw/filtered_feature_bc_matrix.h5 --output /results/scrna --backend numpy --n-neighbors 15 --resolution 1.0
python skills/single-cell-rnaseq/scripts/quality_control.py --input /data/filtered_feature_bc_matrix --output /results/scrna_qc --min-genes 200 --min-counts 500 --max-mt-pct 15
```

## Workflow

1. **Load** — sniff 10x HDF5 (`read_10x_h5`), MTX (`read_10x_mtx`), CSV/TSV (`read_dense_table`), h5ad (`read_h5ad`) → `SimpleAnnData` (anndata optional).
2. **QC metrics** — per-cell `total_counts` / `n_genes_by_counts` / `pct_counts_mt` / `pct_counts_hb` + logs.
3. **QC filter** — cells: `min_genes`, `min_counts`, `max_mt_pct`, opt. `max_hb_pct`; genes: `min_cells_per_gene`; report removal rate.
4. **Normalise** — library size → CP10K (`--target-sum`), then `log1p`.
5. **HVG** — Seurat-style binned z-score (20 bins, CV-corrected), top `--n-top-genes` genes.
6. **Scale + PCA** — per-gene z-score clipped to ±10; `--n-comps` via truncated SVD; variance ratio kept.
7. **Graph** — kNN (`NearestNeighbors`), Jaccard weights, UMAP-style prune threshold, symmetrised + normalised.
8. **Cluster** — Leiden first (`leidenalg` + `igraph`, RBConfiguration), else built-in Louvain; `--resolution` = granularity.
9. **Embeddings** — UMAP (`umap-learn` if available, else built-in spectral init + real force formula), t-SNE (sklearn `TSNE`).
10. **Markers + annotation** — Welch's t-test + BH correction; cell types by gene-overlap scoring vs marker reference (anonymous without).
11. **Figures + report** — UMAP cluster/cell-type plots, cluster-size bar, marker heatmap (PDF + PNG, 300 DPI), report, `run_manifest.json`.

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `--input` | required | 10x h5 / mtx dir / csv / h5ad |
| `--output` | `<results>/scrna` | Output root |
| `--marker-file` | none | Marker CSV (`class_label` / `subclass_label`) |
| `--backend` | `auto` | `auto` / `scanpy` / `numpy` (fallback) |
| `--min-genes` | `200` | Min genes per cell |
| `--min-counts` | `100` | Min UMIs per cell |
| `--max-mt-pct` | `20.0` | Mitochondrial % limit |
| `--max-hb-pct` | none | Haemoglobin % limit (off by default) |
| `--min-cells-per-gene` | `3` | Min cells per gene |
| `--target-sum` | `1e4` | Target total count (CP10K) |
| `--n-top-genes` | `2000` | Highly variable genes (PCA) |
| `--n-comps` | `50` | PCA components |
| `--n-neighbors` | `15` | kNN neighbours k |
| `--resolution` | `1.0` | Clustering resolution (larger = finer) |
| `--min-dist` | `0.3` | Min UMAP distance |
| `--n-markers` | `100` | Markers per cluster |
| `--seed` | `0` | Random seed |
| `--no-tsne` | off | Skip t-SNE (slow) |

## Outputs

7 numbered step subdirectories under `--output`:

- `01_qc/` — `obs_qc.csv` (per-cell counts/genes/mt%/hb%), `qc_summary.json` (counts before/after + thresholds)
- `02_normalized/` — `var_hvg.csv` (per-gene mean, HVG flag)
- `03_dimensionality/` — `embeddings.csv` (PCA/UMAP/t-SNE coords), `pca_variance.csv` (variance ratio)
- `04_clustered/` — cluster labels (to `obs`)
- `05_annotated/` — `cluster_annotations.csv` (cell type + overlap), `marker_genes.csv` (logFC/t value/p/padj), `cluster_summary.csv` (cell + marker counts)
- `06_figures/` — `Figure{1..4}_*.{pdf,png}`: UMAP cluster, UMAP cell-type, cluster-size, heatmap
- `07_report/` — `scrna_report.txt` (full report)

Plus `run_manifest.json` at the output root (input, backend, cell/gene counts, cluster count); `quality_control.py` writes `qc_metrics.csv` / `qc_filtered_cells.csv` / `qc_summary.json`.

## Validation

- **QC removal rate** — usually 5%–30%; above 50% ⇒ thresholds too strict (often `--min-genes` too high); below 1% ⇒ too loose, low-quality cells mixed in.
- **Mitochondrial %** — healthy cells usually < 10%; most cells > 30% ⇒ poorly dissociated or stressed sample, and raising `--max-mt-pct` then only admits noisy cells.
- **Cluster size** — clusters < 1% of the data are usually doublets or batch effects; raise `--resolution` to see whether they persist.
- **UMAP morphology** — clusters should form continuous, non-overlapping blobs; all cells in one point ⇒ `--n-top-genes` too small or PCA too few dims; extreme outlier cells ⇒ unfiltered low-quality cells.
- **Marker plausibility** — top markers should have clear biological meaning; all-ribosomal markers ⇒ no discriminative signal (resolution too low).
- **Embedding consistency** — same cluster in completely different places on UMAP and t-SNE ⇒ unstable neighbourhood graph; raise `--n-neighbors` and re-run.
- **Backend consistency** — `--backend scanpy` vs `--backend numpy` on the same data should give broadly consistent cluster counts and cluster-size distributions (within about ±1 cluster).

## References

- [Methodology](references/methodology.md) — rationale per step, formulas, parameter choices
- [Troubleshooting](references/troubleshooting.md) — common errors, unexpected results, debugging
