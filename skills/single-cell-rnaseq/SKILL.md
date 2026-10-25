---
name: single-cell-rnaseq
description: "This skill should be used when running single-cell RNA-seq analysis: reading 10x/AnnData matrices, QC filtering, normalization, HVG selection, PCA, clustering, UMAP/t-SNE and cell type annotation. 用于单细胞转录组（scRNA-seq）质控、降维、聚类与细胞注释全流程。"
license: MIT
compatibility: "Python 3.9+; numpy, scipy, pandas, scikit-learn required. Optional: matplotlib (figures), scanpy+anndata (scanpy backend), leidenalg+igraph (Leiden clustering), umap-learn (official UMAP)"
metadata:
  version: "1.0"
  last-reviewed: "2026-10-25"
  skill-author: sinder
---

## Overview

本技能把单细胞 RNA-seq 数据从原始计数矩阵处理到可解释的细胞分群与细胞类型注释：自动识别 10x HDF5 / MTX / 稠密 CSV / h5ad 四类输入格式，做基于线粒体比例、基因数与 UMI 数的质控过滤，随后进行 CP10K + log1p 归一化、高变基因筛选、z-score 缩放、截断 SVD 主成分分析、kNN 邻接图构建与图聚类，最后输出 UMAP / t-SNE 二维嵌入、每簇 marker 基因表与可视化图表。

分析提供 **双后端**：`--backend scanpy` 调用官方 scanpy 实现，`--backend numpy` 使用内置的纯 NumPy/SciPy/scikit-learn 降级实现（自实现的 Louvain/UMAP 谱初始化与力优化），在没有 scanpy 的环境中仍能产出统计等价的结果。

## When to Use This Skill

- 用户要做单细胞转录组分析（关键词：single-cell、scRNA-seq、单细胞、细胞图谱、细胞注释、UMAP、聚类）
- 用户有 10x Genomics 原始计数（`matrix.h5` / `filtered_feature_bc_matrix/`）需要走标准 QC 与聚类流程
- 用户需要判断细胞质量（关键词：doublet、双细胞、线粒体比例、mitochondrial、QC、过滤阈值）
- 用户需要找每个细胞群的 marker 基因并据此命名细胞类型
- 用户环境里没有 scanpy，但仍需要可复现的降维与聚类结果
- 用户已有聚类结果，需要补充 t-SNE / marker 热图 / 分析报告

不适用：bulk RNA-seq（用 `bulk-rnaseq`）、空间转录组（CST/Slide-seq 需专门流程）、单细胞 ATAC、轨迹推断（trajectory inference）。

## Quick Start

```bash
export PYTHONPATH="$PWD/packages"

# 完整流程（自动检测输入格式）
python skills/single-cell-rnaseq/scripts/run_pipeline.py \
    --input  /data/filtered_feature_bc_matrix \
    --output /results/scrna \
    --marker-file /ref/cell_marker_reference.csv

# 10x HDF5 + 强制使用降级后端（无 scanpy 环境）
python skills/single-cell-rnaseq/scripts/run_pipeline.py \
    --input /data/raw/filtered_feature_bc_matrix.h5 \
    --output /results/scrna \
    --backend numpy --n-neighbors 15 --resolution 1.0

# 只做质控，不做降维聚类
python skills/single-cell-rnaseq/scripts/quality_control.py \
    --input /data/filtered_feature_bc_matrix \
    --output /results/scrna_qc \
    --min-genes 200 --min-counts 500 --max-mt-pct 15
```

## Workflow

1. **数据读取** — 自动嗅探输入类型：10x HDF5（`read_10x_h5`）、10x MTX 目录（`read_10x_mtx`）、稠密 CSV/TSV（`read_dense_table`）、h5ad（`read_h5ad`）；统一封装为轻量 `SimpleAnnData` 兼容容器，不强制依赖 anndata。
2. **质控指标** — 逐细胞计算 `total_counts` / `n_genes_by_counts` / `pct_counts_mt`（线粒体）/ `pct_counts_hb`（血红蛋白）及对数版本。
3. **质控过滤** — 按 `min_genes`、`min_counts`、`max_mt_pct`、可选 `max_hb_pct` 过滤细胞；按 `min_cells_per_gene` 过滤基因；输出过滤前后细胞数/基因数与去除率。
4. **归一化** — library size 归一化到 CP10K（`--target-sum`），随后 `log1p` 变换稳定方差。
5. **高变基因** — Seurat 风格分箱 z-score 法（20 个表达量分箱，对数期望值/标准差标准化，CV 校正），选出 `--n-top-genes` 个基因。
6. **缩放与 PCA** — 逐基因 z-score 缩放并截断到 ±10；用截断 SVD（稀疏输入）/ 稠密 SVD 得到 `--n-comps` 个主成分，同时记录解释方差比例。
7. **邻接图** — 在 PCA 空间用 sklearn `NearestNeighbors` 找 k 近邻，用 Jaccard 相似度算边权，经 UMAP 式 prune 阈值剪枝后对称化并行归一化。
8. **聚类** — 优先 Leiden（`leidenalg` + `igraph`，RBConfiguration），缺失时回退到内置 Louvain（模块度增益的局部移动 + 多层图聚合）。`--resolution` 控制粒度。
9. **降维嵌入** — UMAP（`umap-learn` 可用时用官方实现；否则用内置谱初始化 + UMAP 真实力公式优化）与 t-SNE（sklearn `TSNE`）。
10. **Marker 与注释** — Welch t 检验 + BH 校正对每簇排序 marker；用 marker 参考文件做基因重叠打分给出细胞类型标签，无参考文件时输出匿名标签。
11. **可视化与报告** — UMAP 分群图、UMAP 细胞类型图、簇大小条形图、marker 热图（PDF + PNG，300 DPI），并生成文字报告与 `run_manifest.json`。

## Parameters

| 参数 | 默认 | 说明 |
|---|---|---|
| `--input` | 必填 | 输入矩阵：10x h5 / 10x mtx 目录 / csv / h5ad |
| `--output` | `<results>/scrna` | 输出根目录，下设 7 个步骤子目录 |
| `--marker-file` | 无 | marker 参考 CSV（含 `class_label` / `subclass_label` / marker 列） |
| `--backend` | `auto` | `auto` / `scanpy` / `numpy`（降级实现） |
| `--min-genes` | `200` | 细胞至少检测到的基因数 |
| `--min-counts` | `100` | 细胞至少的 UMI 数 |
| `--max-mt-pct` | `20.0` | 线粒体基因比例上限（%） |
| `--max-hb-pct` | 无 | 血红蛋白基因比例上限（%），默认不过滤 |
| `--min-cells-per-gene` | `3` | 基因最少在多少细胞中表达 |
| `--target-sum` | `1e4` | 归一化目标总计数（CP10K） |
| `--n-top-genes` | `2000` | 参与 PCA 的高变基因数 |
| `--n-comps` | `50` | PCA 主成分数 |
| `--n-neighbors` | `15` | kNN 邻居数 k |
| `--resolution` | `1.0` | 聚类分辨率（越大簇越细） |
| `--min-dist` | `0.3` | UMAP 最小嵌入距离 |
| `--n-markers` | `100` | 每簇保留的 marker 基因数 |
| `--seed` | `0` | 随机种子 |
| `--no-tsne` | 关 | 跳过 t-SNE（较慢） |

## Outputs

| 文件 | 说明 |
|---|---|
| `01_qc/obs_qc.csv` | 每细胞质控指标（counts / genes / mt% / hb%） |
| `01_qc/qc_summary.json` | 过滤前后细胞数、基因数与所用阈值 |
| `02_normalized/var_hvg.csv` | 每基因的均值与高变标记 |
| `03_dimensionality/embeddings.csv` | PCA / UMAP / t-SNE 坐标合并表 |
| `03_dimensionality/pca_variance.csv` | 各主成分解释方差比例 |
| `04_clustered/` | 聚类标签（写入 `obs`） |
| `05_annotated/cluster_annotations.csv` | 每个 cluster 的细胞类型标签与 marker 重叠数 |
| `05_annotated/marker_genes.csv` | 每簇 marker 基因表（logFC / t 值 / p / padj） |
| `05_annotated/cluster_summary.csv` | 每簇细胞数与 marker 数 |
| `06_figures/Figure{1..4}_*.{pdf,png}` | UMAP 分群图、UMAP 细胞类型图、簇大小图、marker 热图 |
| `07_report/scrna_report.txt` | 完整分析报告 |
| `run_manifest.json` | 运行元数据（输入、后端、细胞/基因数、簇数） |
| `qc_metrics.csv` / `qc_filtered_cells.csv` / `qc_summary.json` | （`quality_control.py` 独立产出） |

## Validation

- **质控过滤率**：细胞去除率通常在 5%–30%。超过 50% 说明阈值过严（常见于 `--min-genes` 偏高），低于 1% 说明阈值过松，应检查是否混入了大量低质量细胞。
- **线粒体比例**：健康细胞通常 < 10%。若绝大多数细胞 > 30%，提示样本解离质量差或存在应激，此时上调 `--max-mt-pct` 只会引入噪声细胞。
- **簇大小**：出现占比 < 1% 的小簇通常是双细胞或批次效应，可通过调高 `--resolution` 观察其是否稳定存在。
- **UMAP 形态**：簇应呈连续、互不重叠的团块。若所有细胞挤成一点，通常是 `--n-top-genes` 过小或 PCA 维度不足；若出现极端离群细胞，检查是否有未过滤的低质量细胞。
- **marker 合理性**：每簇的 top marker 应有明确生物学含义。若某簇的 marker 全是核糖体蛋白基因，多半是簇内缺少区分信号（分辨率过低）。
- **降维一致性**：UMAP 与 t-SNE 上同一 cluster 的位置若完全不一致，说明邻接图不稳定，可增大 `--n-neighbors` 后重跑。
- **后端一致性**：同一数据用 `--backend scanpy` 与 `--backend numpy` 各跑一次，聚类数与簇大小分布应大体一致（±1 个簇以内）。

## References

- [references/methodology.md](references/methodology.md) — 各步骤方法学原理、公式推导与参数选择依据
- [references/troubleshooting.md](references/troubleshooting.md) — 常见报错、结果异常与排错手册
