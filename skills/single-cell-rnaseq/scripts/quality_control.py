#!/usr/bin/env python3
"""独立的单细胞质量控制（QC）步骤。

移植自原仓库 ``scrna_python/scripts/02_quality_control.py``，去掉了硬编码路径、
裸 print，并补齐 argparse 与日志。

QC 指标与原脚本一致（scanpy ``calculate_qc_metrics`` 的语义）：

- ``total_counts``：每个细胞的总 UMI 数
- ``n_genes_by_counts``：每个细胞检测到的基因数
- ``pct_counts_mt``：线粒体基因 UMI 占比（%）
- ``pct_counts_hb``：血红蛋白基因 UMI 占比（%）

过滤阈值默认与原脚本一致：``min_genes=200``、``min_counts=100``、
``max_mt_pct=20``、基因至少在 3 个细胞中表达。

Usage:
    python quality_control.py \\
        --input data/filtered_feature_bc_matrix.h5 \\
        --output results/01_qc
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "packages"))

from biopipe_core import SkillError, get_config, get_logger  # noqa: E402
from biopipe_core.cli import cli_parser, run_cli  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))

from scalgorithms import compute_qc_metrics, filter_cells_genes  # noqa: E402
from scdata import gene_is_hb, gene_is_mitochondrial, read_matrix, summarize_sparsity  # noqa: E402  # noqa: E402

LOGGER = get_logger("biopipe.single-cell-rnaseq")


def run_quality_control(
    input_path: Path,
    output_dir: Path,
    min_genes: int = 200,
    min_counts: int = 100,
    max_mt_pct: float = 20.0,
    max_hb_pct: float | None = None,
    min_cells_per_gene: int = 3,
    make_plots: bool = True,
) -> dict[str, Path]:
    """执行 QC 步骤。

    Args:
        input_path: 输入矩阵（10x h5 / mtx 目录 / csv / h5ad）。
        output_dir: 输出目录。
        min_genes: 细胞至少检测到的基因数。
        min_counts: 细胞至少的 UMI 数。
        max_mt_pct: 线粒体比例上限（%）。
        max_hb_pct: 血红蛋白比例上限（%）。
        min_cells_per_gene: 基因最少表达细胞数。
        make_plots: 是否生成 QC 小提琴图与散点图。

    Returns:
        产出文件路径字典。
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    LOGGER.info("读取数据: %s", input_path)
    adata = read_matrix(input_path)
    LOGGER.info("原始数据: %d cells x %d genes", adata.n_obs, adata.n_vars)

    spars = summarize_sparsity(adata.X)
    LOGGER.info("稀疏度: %.4f (非零元素 %d)", spars["sparsity"], spars["nnz"])

    adata.var_names_make_unique()
    mt_mask = gene_is_mitochondrial(adata.var_names)
    hb_mask = gene_is_hb(adata.var_names)
    LOGGER.info("线粒体基因 %d 个, 血红蛋白基因 %d 个", int(mt_mask.sum()), int(hb_mask.sum()))

    qc = compute_qc_metrics(adata, mt_mask, hb_mask)
    for key, values in qc.items():
        adata.obs[key] = values

    LOGGER.info("QC 指标统计:")
    LOGGER.info("  n_genes_by_counts: min=%.0f median=%.0f max=%.0f",
                qc["n_genes_by_counts"].min(),
                float(pd_median(qc["n_genes_by_counts"])),
                qc["n_genes_by_counts"].max())
    LOGGER.info("  total_counts:      min=%.0f median=%.0f max=%.0f",
                qc["total_counts"].min(),
                float(pd_median(qc["total_counts"])),
                qc["total_counts"].max())
    LOGGER.info("  pct_counts_mt:     min=%.2f%% median=%.2f%% max=%.2f%%",
                qc["pct_counts_mt"].min(),
                float(pd_median(qc["pct_counts_mt"])),
                qc["pct_counts_mt"].max())

    filtered, stats = filter_cells_genes(
        adata, qc,
        min_genes=min_genes, min_counts=min_counts, max_mt_pct=max_mt_pct,
        max_hb_pct=max_hb_pct, min_cells_per_gene=min_cells_per_gene,
    )

    LOGGER.info("过滤后: %d cells x %d genes (移除细胞 %d, 移除基因 %d)",
                stats["n_cells_after"], stats["n_genes_after"],
                stats["n_cells_removed"], stats["n_genes_removed"])
    LOGGER.info("细胞移除比例: %.2f%%", stats["cell_removal_rate"] * 100)

    outputs: dict[str, Path] = {}

    # QC 指标表（保留过滤前的全部细胞，便于画图与追溯）
    qc_table = adata.obs.copy()
    qc_table.to_csv(output_dir / "qc_metrics.csv")
    outputs["qc_metrics"] = output_dir / "qc_metrics.csv"

    filtered.obs.to_csv(output_dir / "qc_filtered_cells.csv")
    outputs["filtered_cells"] = output_dir / "qc_filtered_cells.csv"

    stats_path = output_dir / "qc_summary.json"
    stats_path.write_text(
        json.dumps(
            {
                **stats,
                "input": str(input_path),
                "sparsity": {k: v for k, v in spars.items() if k != "shape"},
                "n_mt_genes": int(mt_mask.sum()),
                "n_hb_genes": int(hb_mask.sum()),
                "qc_metric_stats": {
                    key: {
                        "min": float(values.min()),
                        "median": float(pd_median(values)),
                        "max": float(values.max()),
                    }
                    for key, values in qc.items()
                },
            },
            indent=2,
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )
    outputs["summary"] = stats_path

    if make_plots:
        produced = plot_qc(qc, stats, output_dir)
        outputs.update(produced)

    return outputs


def pd_median(values):
    """返回中位数（避免顶层导入 pandas 时的循环依赖）。"""
    import numpy as np

    return float(np.median(np.asarray(values)))


def plot_qc(qc: dict, stats: dict, output_dir: Path) -> dict[str, Path]:
    """生成 QC 可视化：小提琴图 + 散点图。

    Args:
        qc: QC 指标字典。
        stats: 过滤统计。
        output_dir: 输出目录。

    Returns:
        产出文件路径字典。
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        LOGGER.warning("matplotlib 未安装，跳过 QC 绘图")
        return {}

    import numpy as np

    outputs: dict[str, Path] = {}
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))

    # 小提琴图：过滤前的分布 + 阈值线
    panels = [
        ("n_genes_by_counts", "n_genes_by_counts", stats["thresholds"]["min_genes"], "#3C8DBC"),
        ("total_counts", "total_counts", stats["thresholds"]["min_counts"], "#4DBD74"),
        ("pct_counts_mt", "pct_counts_mt", stats["thresholds"]["max_mt_pct"], "#E64B35"),
    ]
    for ax, (key, label, threshold, color) in zip(axes, panels):
        data = np.asarray(qc[key], dtype=float)
        parts = ax.violinplot([data], showmeans=True, showextrema=False)
        for body in parts["bodies"]:
            body.set_facecolor(color)
            body.set_alpha(0.7)
        ax.axhline(threshold, color="black", linestyle="--", linewidth=0.8, label=f"threshold={threshold}")
        ax.set_xlabel(label)
        ax.set_ylabel("value")
        ax.legend(frameon=False, fontsize=8)
    fig.suptitle("QC distributions (before filtering)")
    fig.tight_layout()
    p = Path(output_dir) / "qc_violin.png"
    fig.savefig(p, dpi=150, bbox_inches="tight")
    plt.close(fig)
    outputs["violin"] = p
    LOGGER.info("QC 小提琴图已保存: %s", p)

    # 散点图：UMI vs 基因数，按线粒体比例着色
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    sc = ax.scatter(
        qc["total_counts"], qc["n_genes_by_counts"],
        c=qc["pct_counts_mt"], s=8, cmap="viridis", alpha=0.7,
    )
    fig.colorbar(sc, ax=ax, label="pct_counts_mt")
    ax.set_xlabel("total_counts")
    ax.set_ylabel("n_genes_by_counts")
    ax.set_title("UMI vs genes (colored by %mt)")
    p = Path(output_dir) / "qc_scatter.png"
    fig.savefig(p, dpi=150, bbox_inches="tight")
    plt.close(fig)
    outputs["scatter"] = p
    LOGGER.info("QC 散点图已保存: %s", p)

    return outputs


def main(args) -> int:
    """CLI 入口。"""
    input_path = Path(args.input).expanduser()
    if not input_path.exists():
        raise SkillError(f"输入路径不存在: {input_path}")

    output = Path(args.output).expanduser() if args.output else get_config().output_dir / "01_qc"

    outputs = run_quality_control(
        input_path=input_path,
        output_dir=output,
        min_genes=args.min_genes,
        min_counts=args.min_counts,
        max_mt_pct=args.max_mt_pct,
        max_hb_pct=args.max_hb_pct,
        min_cells_per_gene=args.min_cells_per_gene,
        make_plots=not args.no_plots,
    )
    for name, path in outputs.items():
        LOGGER.info("%s -> %s", name, path)
    return 0


if __name__ == "__main__":
    DESCRIPTION = "单细胞质量控制：计算 QC 指标并按阈值过滤细胞与基因"

    def build_parser() -> argparse.ArgumentParser:
        """构造带业务参数的解析器（沿用 biopipe_core 的公共参数）。"""
        parser = cli_parser(prog="single-cell-rnaseq", description=DESCRIPTION, skill="single-cell-rnaseq")
        parser.add_argument("--input", required=True, help="输入矩阵（10x h5 / mtx 目录 / csv / h5ad）")
        parser.add_argument("--output", help="输出目录（默认 <results>/01_qc）")
        parser.add_argument("--min-genes", type=int, default=200, help="细胞至少检测到的基因数（默认 200）")
        parser.add_argument("--min-counts", type=int, default=100, help="细胞至少的 UMI 数（默认 100）")
        parser.add_argument("--max-mt-pct", type=float, default=20.0, help="线粒体比例上限%%（默认 20）")
        parser.add_argument("--max-hb-pct", type=float, default=None, help="血红蛋白比例上限%%（默认不过滤）")
        parser.add_argument("--min-cells-per-gene", type=int, default=3, help="基因最少表达细胞数（默认 3）")
        parser.add_argument("--no-plots", action="store_true", help="跳过绘图")
        return parser

    sys.exit(run_cli(main, parser=build_parser(), skill="single-cell-rnaseq"))
