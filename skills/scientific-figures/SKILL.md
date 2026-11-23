---
name: scientific-figures
description: "Use when generating publication-ready scientific figures: ROC curves, confusion matrix heatmaps, docking binding-energy heatmaps, pLDDT confidence histograms, SHAP summary beeswarms, or class distribution charts. Renders each chart type with a unified journal-grade visual style. 用于生成顶刊级出版图表（ROC 曲线、混淆矩阵、热图、SHAP 摘要等）。"
license: MIT
compatibility: Python 3.9+；需要 matplotlib>=3.6（缺失时给出清晰安装提示而不抛 ImportError）
metadata:
  version: "1.0"
  last-reviewed: "2026-10-25"
  skill-author: sinder
  edam-topic: "topic_0622, topic_2269"
  edam-operation: "operation_2935"
  edam-input-data: "data_0928, data_0951"
  edam-output-data: "data_1710, data_0862"
  edam-input-format: "format_3752"
  edam-output-format: "format_1950, format_3464"

---

## Overview

Generates publication-grade figures through one multi-subcommand entry point. All
chart types share a single visual specification (Arial-first font stack, top/right
spines removed, grid disabled, colour-blind-safe palette, 300/600 dpi, tight bbox).

The unified `figure_style.py` module consolidates style logic that used to be
scattered across 15 separate scripts into reusable components, supports three
journal presets (`nature` / `cell` / `science`), and adds real statistical
capability: ROC significance testing, bootstrap confidence intervals and SHAP
feature-importance ranking.

## When to Use This Skill

- Standard charts for a paper or report — ROC curves, confusion matrices, heatmaps
- Several subfigures that must share exactly one style (a paper cannot mix Arial in
  one figure and Times in another)
- SHAP summary or feature-importance plots for machine-learning results
- A binding-energy heatmap for molecular docking results
- A pLDDT confidence distribution for structure-prediction results

## Quick Start

```bash
export PYTHONPATH="$PWD/packages"

# ROC curve (AUC significance testing + bootstrap CI)
python skills/scientific-figures/scripts/make_figures.py roc \
    --input data/roc.csv --output figures/fig8_roc \
    --bootstrap-ci --fill-best --stats-output figures/fig8_roc_stats.csv

# Confusion-matrix heatmap (multi-model)
python skills/scientific-figures/scripts/make_figures.py confusion \
    --input data/confusion.csv --output figures/fig9_cm

# Molecular-docking binding-energy heatmap
python skills/scientific-figures/scripts/make_figures.py heatmap \
    --input data/docking.csv --output figures/fig11_docking \
    --cmap coolwarm_r --colorbar-label "Binding Energy (kcal/mol)" \
    --x-label Receptor --y-label Ligand --percentile-clip

# pLDDT confidence histogram
python skills/scientific-figures/scripts/make_figures.py plddt \
    --input data/plddt.csv --output figures/fig3_plddt

# SHAP summary beeswarm + feature-importance bar chart
python skills/scientific-figures/scripts/make_figures.py shap \
    --input data/shap.csv --output figures/fig10_shap --top-n 15 \
    --importance-output figures/fig10_importance.csv

# Phylum-class two-level taxonomic distribution
python skills/scientific-figures/scripts/make_figures.py class-distribution \
    --input data/taxa.csv --output figures/fig1_taxa --top-n 20
```

## Workflow

1. **Load the style** — `apply_publication_style(style)` sets rcParams (font,
   spines, grid, dpi, palette)
2. **Read the data** — automatic BOM/delimiter sniffing, case-insensitive column
   name matching
3. **Dispatch by subcommand** — each is an independent plotting function sharing a
   unified save interface
4. **Statistical enrichment** (optional) — ROC: trapezoidal AUC, Z-test between
   AUCs (Hanley-McNeil standard error), bootstrap CI · Confusion matrix: globally
   unified colour scale, dual annotation of counts and percentages · SHAP: ranking
   by mean absolute SHAP value, standard deviation, Top-N selection
5. **Save** — `save_figure()` writes PNG (raster, for previewing and docs) and PDF
   (vector, for submission)

## Parameters

| Parameter | Description | Default |
| --- | --- | --- |
| `subcommand` | Chart type (positional): `roc` `confusion` `heatmap` `plddt` `shap` `class-distribution` `secondary-structure` `model-performance` `feature-importance` | required |
| `--input` | Input CSV | required |
| `--output` | Output path (no extension; one file per `--formats` entry) | required |
| `--style` | Journal preset: `default` / `nature` / `cell` / `science` | `nature` |
| `--formats` | Output formats | `png pdf` |
| `--dpi` | Resolution | `300` |
| `--top-n` | Top-N entries shown (SHAP / class distribution / feature importance) | `15` |
| `--bins` | Histogram bins (pLDDT) | `50` |
| `--cmap` | Heatmap colour map | `coolwarm_r` |
| `--colorbar-label` / `--x-label` / `--y-label` | Heatmap labels | `Value` / empty |
| `--percentile-clip` | Clip heatmap colour scale at 2%/98% percentiles | off |
| `--fill-best` | Fill area below the best model curve (ROC) | off |
| `--bootstrap-ci` | Bootstrap confidence interval of the AUC (ROC) | off |
| `--stats-output` | CSV of ROC significance test results | none |
| `--importance-output` | CSV of SHAP feature-importance values | none |
| `--row-labels` | Whether the first heatmap column is used as row labels | `auto` |

## Outputs

| File | Description |
| --- | --- |
| `<output>.png` | 300/600 dpi raster, for previews, PPT and Word embedding |
| `<output>.pdf` | Vector graphic, for submission when the journal requires vector art |
| `<output>_beeswarm.png/pdf` | SHAP beeswarm plot (`shap` subcommand) |
| `<output>_importance.png/pdf` | Feature-importance bar chart (`shap` subcommand) |
| `<output>_stats.csv` | Pairwise AUC Z-test results between ROC models (`--stats-output`) |
| `<output>_importance.csv` | `mean_abs_shap` ranking of all features (`--importance-output`) |

## Validation

- **ROC** — the AUC label on the figure must match the AUC in `--stats-output`; the
  random baseline must be annotated `AUC = 0.50`
- **ROC statistics** — `--n-samples` defaults to 500 and is the assumed sample size
  behind the Hanley-McNeil approximation. If your real sample size differs, pass it
  explicitly, otherwise the p-values are biased
- **Confusion matrix** — all subplots share one colour scale (global vmin/vmax), so
  cell colours **are comparable across models**. If shades look inconsistent, the
  data was split across multiple rows
- **Heatmap** — after `--percentile-clip`, extreme values land on the scale
  boundaries, which must be stated in the figure caption
- **pLDDT** — the three-bin colour scheme (>90 green / 70-90 blue / <70 yellow)
  matches the official AlphaFold bins
- **SHAP** — the beeswarm y-axis order must match the importance bar chart (both
  share the same `top_features` list)
- **General** — check the actual PNG pixel dimensions and whether any text is
  clipped; tight bbox is normally correct, but very long x-axis labels still need a
  visual check
- **Fonts** — when Arial is missing, matplotlib falls back to Helvetica/DejaVu Sans.
  Before final submission, confirm Arial is installed on the target machine,
  otherwise the layout shifts subtly

## References

- `references/methodology.md` — the style specification, the statistical method
  behind each chart type, and the reasoning behind palette choices
- `references/troubleshooting.md` — common errors, column-name mismatches, font and
  output problems
