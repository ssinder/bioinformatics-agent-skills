#!/usr/bin/env python3
"""单细胞 RNA-seq 端到端分析流程。

整合原仓库 ``scrna_python/scripts/`` 下的 6 个步骤：

1. **读取数据** —— 10x HDF5 / 10x MTX / 稠密 CSV / h5ad
2. **质量控制** —— 线粒体比例、基因数、UMI 数阈值过滤细胞与基因
3. **归一化** —— CP10K + log1p，高变基因选择，z-score 缩放
4. **降维** —— PCA（50 主成分）、kNN 邻接图、UMAP、t-SNE
5. **聚类** —— Leiden（无 leidenalg 时用内置 Louvain）
6. **细胞注释** —— 基于 marker 基因的重叠打分注释 + 每簇 marker 表

**双后端**：脚本检测 scanpy 是否可用，可用则走 scanpy 完整流程，
不可用则走内置的纯 NumPy/SciPy/sklearn 降级后端。两个后端都是真实算法实现。

Usage:
    python run_pipeline.py \\
        --input data/filtered_feature_bc_matrix.h5 \\
        --output results/scrna \\
        --min-genes 200 --min-counts 100 --max-mt-pct 20
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "packages"))

from biopipe_core import SkillError, get_config, get_logger  # noqa: E402
from biopipe_core.cli import cli_parser, run_cli  # noqa: E402
from biopipe_core.optional_deps import is_available  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))

from typing import Optional

from scalgorithms import SCANPY_AVAILABLE, run_pipeline  # noqa: E402
from scdata import read_matrix, summarize_sparsity  # noqa: E402

LOGGER = get_logger("biopipe.single-cell-rnaseq")

DIR_QC = "01_qc"
DIR_NORM = "02_normalized"
DIR_DIM = "03_dimensionality"
DIR_CLUSTER = "04_clustered"
DIR_ANNOT = "05_annotated"
DIR_FIG = "06_figures"
DIR_REPORT = "07_report"


def annotate_clusters(markers, adata, marker_file: Optional[Path],
                      top_n: int = 30) -> tuple[dict[str, dict], Optional[str]]:
    """基于 marker 基因重叠给 cluster 打细胞类型标签。

    优先使用用户提供的 marker 参考文件（CSV，需含 ``class_label``、
    ``subclass_label`` 与 marker 基因列）；未提供时退化为按 marker 基因
    数量排序的匿名标签。

    Args:
        markers: :func:`scalgorithms.rank_genes_groups` 的输出。
        adata: 数据对象（提供基因名）。
        marker_file: marker 参考文件路径。
        top_n: 每个 cluster 用于匹配的前 N 个 marker 基因。

    Returns:
        ``(annotations, source)``：cluster 到注释字典的映射，以及数据来源描述。
    """
    annotations: dict[str, dict] = {}
    if markers is None or len(markers) == 0:
        return annotations, None

    cluster_markers: dict[str, list[str]] = {}
    for cluster, grp in markers.groupby("cluster"):
        cluster_markers[str(cluster)] = grp["gene"].astype(str).tolist()[:top_n]

    if marker_file is not None and Path(marker_file).exists():
        import pandas as pd

        ref_df = pd.read_csv(marker_file)
        marker_col = None
        for cand in ("subclass.markers.combo", "markers", "marker_genes", "markers.combo"):
            if cand in ref_df.columns:
                marker_col = cand
                break
        if marker_col is None:
            LOGGER.warning("marker 参考文件缺少 marker 基因列，跳过注释: %s", marker_file)
            marker_file = None
        else:
            class_col = "class_label" if "class_label" in ref_df.columns else ref_df.columns[0]
            sub_col = "subclass_label" if "subclass_label" in ref_df.columns else ref_df.columns[1]
            ref_markers: dict[tuple[str, str], list[str]] = {}
            for _, row in ref_df.iterrows():
                raw = row.get(marker_col)
                if raw is None or (isinstance(raw, float) and pd.isna(raw)):
                    continue
                genes = [g.strip() for g in str(raw).split(",") if g.strip()]
                if genes:
                    ref_markers[(str(row[class_col]), str(row[sub_col]))] = genes

            for cluster, genes in cluster_markers.items():
                gene_set = set(genes)
                best = {"class": "Unknown", "subclass": "Unknown", "overlap": 0}
                for (cls, sub), ref_genes in ref_markers.items():
                    overlap = len(gene_set & set(ref_genes))
                    if overlap > best["overlap"]:
                        best = {"class": cls, "subclass": sub, "overlap": overlap}
                annotations[cluster] = best
            return annotations, f"marker 参考文件: {marker_file}"

    # 匿名注释：按 marker 数量与平均 |logFC| 打分
    if markers is not None and len(markers):
        agg = markers.groupby("cluster").agg(
            n_markers=("gene", "count"),
            mean_lfc=("logfoldchanges", lambda s: float(s.abs().mean())),
        )
        for cluster in cluster_markers:
            annotations[cluster] = {
                "class": "Unassigned",
                "subclass": f"Cluster_{cluster}",
                "overlap": 0,
                "n_markers": int(agg.loc[cluster, "n_markers"]) if cluster in agg.index else 0,
            }
    return annotations, "无 marker 参考文件，使用匿名标签"


def plot_results(adata, markers, annotations, fig_dir: Path, top_n_genes: int = 20) -> list[Path]:
    """生成单细胞分析的标准图表。

    产出：QC 小提琴图、UMAP cluster 图、UMAP 细胞类型图、
    cluster 大小条形图、marker 热图。

    Args:
        adata: 处理后的数据对象。
        markers: marker 基因表。
        annotations: cluster 注释字典。
        fig_dir: 输出目录。
        top_n_genes: 热图展示的基因数。

    Returns:
        产出文件路径列表。
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        LOGGER.warning("matplotlib 未安装，跳过可视化")
        return []

    import numpy as np
    import pandas as pd

    fig_dir.mkdir(parents=True, exist_ok=True)
    produced: list[Path] = []

    umap = adata.obsm.get("X_umap")
    if umap is None:
        LOGGER.warning("缺少 UMAP 坐标，跳过 UMAP 相关图")
        return produced

    clusters = adata.obs["cluster"].astype(str).to_numpy()
    uniq = sorted(set(clusters))

    # --- UMAP by cluster ---
    fig, ax = plt.subplots(figsize=(6, 5))
    cmap = plt.get_cmap("tab20")
    for i, c in enumerate(uniq):
        m = clusters == c
        ax.scatter(umap[m, 0], umap[m, 1], s=8, alpha=0.75,
                   color=cmap(i % 20), label=f"C{c}", edgecolors="none")
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_title("Clusters")
    ax.legend(frameon=False, fontsize=7, markerscale=2, loc="best")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for suffix in ("pdf", "png"):
        p = fig_dir / f"Figure1_UMAP_Clusters.{suffix}"
        fig.savefig(p, format=suffix)
        produced.append(p)
    plt.close(fig)
    LOGGER.info("UMAP cluster 图已生成")

    # --- UMAP by cell type ---
    if annotations:
        cell_types = {c: a.get("class", "Unknown") for c, a in annotations.items()}
        ct = np.array([cell_types.get(c, "Unknown") for c in clusters])
        fig, ax = plt.subplots(figsize=(6, 5))
        for i, t in enumerate(sorted(set(ct))):
            m = ct == t
            ax.scatter(umap[m, 0], umap[m, 1], s=8, alpha=0.75,
                       color=cmap(i % 20), label=str(t)[:24], edgecolors="none")
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")
        ax.set_title("Cell type annotation")
        ax.legend(frameon=False, fontsize=7, markerscale=2)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        for suffix in ("pdf", "png"):
            p = fig_dir / f"Figure2_UMAP_CellTypes.{suffix}"
            fig.savefig(p, format=suffix)
            produced.append(p)
        plt.close(fig)
        LOGGER.info("UMAP 细胞类型图已生成")

    # --- cluster 大小 ---
    counts = pd.Series(clusters).value_counts().sort_index()
    fig, ax = plt.subplots(figsize=(max(5, 0.4 * len(counts)), 4))
    ax.bar(counts.index.astype(str), counts.values, color="#3C8DBC")
    ax.set_xlabel("Cluster")
    ax.set_ylabel("Number of cells")
    ax.set_title("Cell distribution by cluster")
    for i, v in enumerate(counts.values):
        ax.text(i, v, str(v), ha="center", va="bottom", fontsize=7)
    for suffix in ("pdf", "png"):
        p = fig_dir / f"Figure3_Cluster_Sizes.{suffix}"
        fig.savefig(p, format=suffix)
        produced.append(p)
    plt.close(fig)
    LOGGER.info("cluster 分布图已生成")

    # --- marker 热图 ---
    if markers is not None and len(markers):
        import scipy.sparse as sp

        top_genes: list[str] = []
        for _cluster, grp in markers.groupby("cluster"):
            for gene in grp["gene"].astype(str).tolist()[:5]:
                if gene not in top_genes:
                    top_genes.append(gene)
        top_genes = top_genes[:top_n_genes]

        var_names = [str(v) for v in adata.var_names]
        gene_pos = {g: i for i, g in enumerate(var_names)}
        X = adata.X
        rows = []
        for g in top_genes:
            if g not in gene_pos:
                continue
            idx = gene_pos[g]
            vals = np.asarray(X[:, idx].todense()).ravel() if sp.issparse(X) else np.asarray(X[:, idx]).ravel()
            rows.append(vals)
        if rows:
            mat = np.vstack(rows)
            # 逐基因 z-score 后按 cluster 取均值做热图
            mat = mat - mat.mean(axis=1, keepdims=True)
            sd = mat.std(axis=1, keepdims=True)
            sd[sd == 0] = 1
            mat = mat / sd
            # 按 cluster 分组求均值
            order = sorted(set(clusters))
            cmat = np.vstack([
                np.array([mat[:, clusters == c].mean(axis=1) for c in order]).T
            ])
            fig, ax = plt.subplots(figsize=(1 + 0.55 * len(order), max(4, 0.22 * len(top_genes))))
            im = ax.imshow(cmat, aspect="auto", cmap="RdBu_r",
                           vmin=-2, vmax=2)
            ax.set_yticks(np.arange(len(top_genes)))
            ax.set_yticklabels(top_genes, fontsize=6)
            ax.set_xticks(np.arange(len(order)))
            ax.set_xticklabels([f"C{c}" for c in order], rotation=90, fontsize=7)
            ax.set_title("Marker gene expression (z-score by gene)")
            fig.colorbar(im, ax=ax, shrink=0.6, label="z-score")
            for suffix in ("pdf", "png"):
                p = fig_dir / f"Figure4_Marker_Heatmap.{suffix}"
                fig.savefig(p, format=suffix)
                produced.append(p)
            plt.close(fig)
            LOGGER.info("marker 热图已生成")

    return produced


def generate_report(output_root: Path, result, annotations, ann_source: str,
                    report_path: Path) -> Path:
    """生成单细胞分析报告。

    Args:
        output_root: 输出根目录。
        result: :class:`scalgorithms.PipelineResult`。
        annotations: cluster 注释字典。
        ann_source: 注释数据来源描述。
        report_path: 报告输出路径。

    Returns:
        报告路径。
    """
    lines: list[str] = []
    lines.append("=" * 78)
    lines.append("SINGLE-CELL RNA-SEQ ANALYSIS REPORT")
    lines.append("=" * 78)
    lines.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"Backend: {result.backend}")
    lines.append(f"Output root: {output_root}")

    qc = result.qc_stats
    lines.append("")
    lines.append("-" * 78)
    lines.append("1. DATA & QUALITY CONTROL")
    lines.append("-" * 78)
    lines.append(f"  Cells before QC: {qc['n_cells_before']}")
    lines.append(f"  Cells after QC:  {qc['n_cells_after']} "
                 f"(removed {qc['n_cells_removed']}, {qc['cell_removal_rate'] * 100:.2f}%)")
    lines.append(f"  Genes before:    {qc['n_genes_before']}")
    lines.append(f"  Genes after:     {qc['n_genes_after']} (removed {qc['n_genes_removed']})")
    lines.append("  Thresholds:")
    for key, value in qc["thresholds"].items():
        lines.append(f"    {key} = {value}")

    lines.append("")
    lines.append("-" * 78)
    lines.append("2. CLUSTERING")
    lines.append("-" * 78)
    lines.append(f"  Number of clusters: {len(set(result.adata.obs['cluster']))}")
    if result.cluster_summary is not None and len(result.cluster_summary):
        lines.append("")
        lines.append("  Cluster  Cells  Markers  Cell type")
        for _, row in result.cluster_summary.iterrows():
            cluster = str(row["cluster"])
            ann = annotations.get(cluster, {})
            lines.append(
                f"    {cluster:<7} {int(row['n_cells']):>5}  {int(row.get('n_markers', 0)):>7}  "
                f"{ann.get('class', '?')} / {ann.get('subclass', '?')}"
            )

    lines.append("")
    lines.append("-" * 78)
    lines.append("3. ANNOTATION SOURCE")
    lines.append("-" * 78)
    lines.append(f"  {ann_source}")

    lines.append("")
    lines.append("-" * 78)
    lines.append("4. TOP MARKER GENES PER CLUSTER")
    lines.append("-" * 78)
    if result.marker_table is not None and len(result.marker_table):
        for cluster, grp in result.marker_table.groupby("cluster"):
            top = grp.head(8)
            genes = ", ".join(
                f"{r['gene']}(LFC={r['logfoldchanges']:.2f}, padj={r['pvals_adj']:.2e})"
                for _, r in top.iterrows()
            )
            lines.append(f"  Cluster {cluster}: {genes}")

    lines.append("")
    lines.append("-" * 78)
    lines.append("5. METHODS")
    lines.append("-" * 78)
    lines.append(f"""
  Backend         : {result.backend}
  QC              : total_counts / n_genes_by_counts / pct_counts_mt / pct_counts_hb
  Normalization   : library-size normalization (CP10K) + log1p
  HVG             : Seurat-style binned z-score, {result.qc_stats.get('n_top_genes', 'n/a')}
  Scaling         : z-score, clipped at 10
  PCA             : truncated SVD on scaled HVG matrix
  Neighbors       : kNN graph with Jaccard edge weights (UMAP-style pruning)
  Clustering      : Leiden (igraph/leidenalg) or built-in Louvain fallback
  Embedding       : UMAP (umap-learn) or MDS + force-directed approximation; t-SNE (sklearn)
  Markers         : Welch t-test per cluster vs rest, BH FDR correction
""".rstrip())

    lines.append("")
    lines.append("-" * 78)
    lines.append("6. OUTPUT FILES")
    lines.append("-" * 78)
    for sub in (DIR_QC, DIR_NORM, DIR_DIM, DIR_CLUSTER, DIR_ANNOT, DIR_FIG):
        p = output_root / sub
        if p.exists():
            n = len([x for x in p.rglob("*") if x.is_file()])
            lines.append(f"  {sub}/  ({n} files)")
    lines.append("")
    lines.append("=" * 78)
    lines.append("END OF REPORT")
    lines.append("=" * 78)

    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(lines), encoding="utf-8")
    LOGGER.info("报告已生成: %s", report_path)
    return report_path


def main(args) -> int:
    """单细胞流程主函数。"""
    # 依赖守卫：numpy / pandas / scipy 缺失时给出明确的安装指引，
    # 而不是让 --help 之类的元信息查询也一并崩溃。
    for module in ("numpy", "pandas", "scipy.sparse"):
        if not is_available(module):
            raise SkillError(
                f"缺少依赖 {module}，单细胞流程无法运行。"
                f"\n请先安装: pip install {module}"
                f"\n或一次性安装: pip install 'biopipe-agent-skills[single-cell]'"
            )
    config = get_config(args.config)
    output_root = Path(args.output).expanduser() if args.output else config.output_dir / "scrna"
    output_root = output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    input_path = Path(args.input).expanduser()
    if not input_path.exists():
        raise SkillError(f"输入路径不存在: {input_path}")

    LOGGER.info("=" * 60)
    LOGGER.info("单细胞分析流程")
    LOGGER.info("=" * 60)
    LOGGER.info("输入: %s", input_path)
    LOGGER.info("输出根目录: %s", output_root)
    LOGGER.info(
        "后端: %s",
        "scanpy" if (args.backend == "auto" and SCANPY_AVAILABLE)
        else ("scanpy (强制)" if args.backend == "scanpy" else "numpy 降级后端"),
    )

    # --- 读取 ---
    LOGGER.info("[1/6] 读取数据")
    adata = read_matrix(input_path)
    LOGGER.info("  %d cells x %d genes", adata.n_obs, adata.n_vars)
    spars = summarize_sparsity(adata.X)
    LOGGER.info("  稀疏度 %.4f", spars["sparsity"])

    # --- 主流程 ---
    LOGGER.info("[2/6] 质控 -> 归一化 -> 降维 -> 聚类 -> 注释")
    result = run_pipeline(
        adata,
        backend=args.backend,
        min_genes=args.min_genes,
        min_counts=args.min_counts,
        max_mt_pct=args.max_mt_pct,
        max_hb_pct=args.max_hb_pct,
        min_cells_per_gene=args.min_cells_per_gene,
        target_sum=args.target_sum,
        n_top_genes=args.n_top_genes,
        n_comps=args.n_comps,
        n_neighbors=args.n_neighbors,
        resolution=args.resolution,
        min_dist=args.min_dist,
        n_marker_genes=args.n_markers,
        seed=args.seed,
        compute_tsne_embedding=not args.no_tsne,
        logger=lambda msg: LOGGER.info("  %s", msg),
    )

    qc = result.qc_stats
    LOGGER.info("质控后: %d cells x %d genes",
                qc["n_cells_after"], qc["n_genes_after"])

    # --- 分步保存 ---
    qc_dir = output_root / DIR_QC
    norm_dir = output_root / DIR_NORM
    dim_dir = output_root / DIR_DIM
    cluster_dir = output_root / DIR_CLUSTER
    ann_dir = output_root / DIR_ANNOT
    fig_dir = output_root / DIR_FIG
    for d in (qc_dir, norm_dir, dim_dir, cluster_dir, ann_dir, fig_dir):
        d.mkdir(parents=True, exist_ok=True)

    qc_dir.joinpath("qc_summary.json").write_text(
        json.dumps(qc, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    result.adata.obs.to_csv(qc_dir / "obs_qc.csv")
    LOGGER.info("[3/6] QC 结果已保存")

    result.adata.var.to_csv(norm_dir / "var_hvg.csv")
    LOGGER.info("[4/6] 归一化信息已保存")

    import numpy as np
    import pandas as pd

    emb_cols = [c for c in ("X_pca", "X_umap", "X_tsne") if c in result.adata.obsm]
    if emb_cols:
        emb = pd.DataFrame(
            np.column_stack([np.asarray(result.adata.obsm[c]) for c in emb_cols]),
            index=result.adata.obs_names,
            columns=[f"{c}_{i}" for c in emb_cols
                     for i in range(np.asarray(result.adata.obsm[c]).shape[1])],
        )
        emb.to_csv(dim_dir / "embeddings.csv")
    if "pca" in result.adata.uns:
        pd.DataFrame({
            "PC": [f"PC{i + 1}" for i in range(len(result.adata.uns["pca"]["variance_ratio"]))],
            "variance_ratio": result.adata.uns["pca"]["variance_ratio"],
        }).to_csv(dim_dir / "pca_variance.csv")
    LOGGER.info("[5/6] 降维结果已保存")

    # --- 注释 ---
    annotations, ann_source = annotate_clusters(
        result.marker_table, result.adata,
        Path(args.marker_file).expanduser() if args.marker_file else None,
    )
    if annotations:
        ann_df = pd.DataFrame(annotations).T
        ann_df.index.name = "cluster"
        ann_df.to_csv(ann_dir / "cluster_annotations.csv")
        result.adata.obs["cell_class"] = result.adata.obs["cluster"].astype(str).map(
            lambda c: annotations.get(c, {}).get("class", "Unknown")
        )
        result.adata.obs["cell_subclass"] = result.adata.obs["cluster"].astype(str).map(
            lambda c: annotations.get(c, {}).get("subclass", "Unknown")
        )
    LOGGER.info("[6/6] 细胞注释完成 (%s)", ann_source)

    if result.marker_table is not None and len(result.marker_table):
        result.marker_table.to_csv(ann_dir / "marker_genes.csv", index=False)
    if result.cluster_summary is not None and len(result.cluster_summary):
        result.cluster_summary.to_csv(ann_dir / "cluster_summary.csv", index=False)

    # --- 可视化 ---
    figures = plot_results(result.adata, result.marker_table, annotations,
                            fig_dir, top_n_genes=args.top_heatmap_genes)

    # --- 报告 ---
    report_path = output_root / DIR_REPORT / "scrna_report.txt"
    generate_report(output_root, result, annotations, ann_source, report_path)

    manifest = {
        "generated": datetime.now().isoformat(timespec="seconds"),
        "input": str(input_path),
        "output_root": str(output_root),
        "backend": result.backend,
        "scanpy_available": SCANPY_AVAILABLE,
        "n_cells": qc["n_cells_after"],
        "n_genes": qc["n_genes_after"],
        "n_clusters": int(len(set(result.adata.obs["cluster"]))),
        "figures": [str(p) for p in figures],
    }
    (output_root / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    LOGGER.info("=" * 60)
    LOGGER.info("流程完成。输出根目录: %s", output_root)
    LOGGER.info("=" * 60)
    return 0


if __name__ == "__main__":
    DESCRIPTION = (
        "单细胞 RNA-seq 端到端流程：读取 -> 质控 -> CP10K 归一化 -> HVG -> PCA -> "
        "kNN -> 聚类 -> UMAP/t-SNE -> 细胞注释"
    )

    def build_parser() -> argparse.ArgumentParser:
        """构造带业务参数的解析器（沿用 biopipe_core 的公共参数）。"""
        parser = cli_parser(prog="single-cell-rnaseq", description=DESCRIPTION,
                            skill="single-cell-rnaseq")

        g_in = parser.add_argument_group("输入")
        g_in.add_argument("--input", required=True,
                          help="输入矩阵（10x h5 / 10x mtx 目录 / csv / h5ad）")
        g_in.add_argument("--marker-file", help="marker 参考 CSV（class_label/subclass_label/markers）")

        g_out = parser.add_argument_group("输出")
        g_out.add_argument("--output", help="输出根目录（默认 <results>/scrna）")

        g_qc = parser.add_argument_group("质控")
        g_qc.add_argument("--min-genes", type=int, default=200, help="细胞至少检测到的基因数（默认 200）")
        g_qc.add_argument("--min-counts", type=int, default=100, help="细胞至少的 UMI 数（默认 100）")
        g_qc.add_argument("--max-mt-pct", type=float, default=20.0, help="线粒体比例上限%%（默认 20）")
        g_qc.add_argument("--max-hb-pct", type=float, default=None, help="血红蛋白比例上限%%（默认不过滤）")
        g_qc.add_argument("--min-cells-per-gene", type=int, default=3, help="基因最少表达细胞数（默认 3）")

        g_an = parser.add_argument_group("分析参数")
        g_an.add_argument("--backend", choices=["auto", "scanpy", "numpy"], default="auto",
                          help="分析后端：auto=自动检测，scanpy=强制 scanpy，numpy=降级实现")
        g_an.add_argument("--target-sum", type=float, default=1e4, help="归一化目标总计数（默认 1e4）")
        g_an.add_argument("--n-top-genes", type=int, default=2000, help="高变基因数（默认 2000）")
        g_an.add_argument("--n-comps", type=int, default=50, help="PCA 主成分数（默认 50）")
        g_an.add_argument("--n-neighbors", type=int, default=15, help="kNN 邻居数（默认 15）")
        g_an.add_argument("--resolution", type=float, default=1.0, help="聚类分辨率（默认 1.0）")
        g_an.add_argument("--min-dist", type=float, default=0.3, help="UMAP 最小距离（默认 0.3）")
        g_an.add_argument("--n-markers", type=int, default=100, help="每簇 marker 基因数（默认 100）")
        g_an.add_argument("--seed", type=int, default=0, help="随机种子（默认 0）")
        g_an.add_argument("--no-tsne", action="store_true", help="跳过 t-SNE（较慢）")
        g_an.add_argument("--top-heatmap-genes", type=int, default=20, help="热图展示基因数（默认 20）")

        return parser

    sys.exit(run_cli(main, parser=build_parser(), skill="single-cell-rnaseq"))
