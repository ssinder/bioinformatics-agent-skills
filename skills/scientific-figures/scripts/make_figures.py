#!/usr/bin/env python3
"""出版级科学图表生成入口（多子命令）。

统一入口：``make_figures.py <子命令> --input ... --output ...``

支持的子命令
------------
============  ==========================================  ==============================
子命令        图类型                                       输入 CSV 结构
============  ==========================================  ==============================
roc           ROC 曲线（多模型 + AUC）                    Model, FPR, TPR, AUC
confusion     混淆矩阵热力图（多模型）                   Model, TN, FP, FN, TP
heatmap       分子对接结合能热图                          矩阵（首行/首列为标签）
plddt         结构预测置信度（pLDDT）直方图               含 pLDDT 数值的列
shap          SHAP 摘要蜂群图 + 特征重要性条形图          特征列 × 样本
class-distribution 门-纲双层级分类分布条形图              kingdom, phylum, class, count
secondary-structure  二级结构分布柱状图                    ss_name（每行一个残基）
model-performance   模型多指标对比柱状图                   Model, Accuracy, ...
feature-importance  特征重要性柱状图                        feature, mean_abs_shap, ...
============  ==========================================  ==============================

matplotlib 未安装时会给出清晰的安装提示，而不是 ImportError 堆栈。

用法示例::

    python make_figures.py roc --input roc.csv --output figures/roc
    python make_figures.py confusion --input cm.csv --output figures/cm
    python make_figures.py heatmap --input docking.csv --output figures/heatmap
    python make_figures.py plddt --input plddt.csv --output figures/plddt
    python make_figures.py shap --input shap.csv --output figures/shap
    python make_figures.py class-distribution --input taxa.csv --output figures/taxa
"""

from __future__ import annotations

import argparse
import csv
import math
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "packages"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from typing import Optional

from biopipe_core import (  # noqa: E402
    SkillError,
    cli_parser,
    get_config,
    get_logger,
    resolve_path,
    run_cli,
)
from figure_style import (  # noqa: E402
    MATPLOTLIB_AVAILABLE,
    PAPER_COLORS,
    SEMANTIC_COLORS,
    FigureStyleError,
    add_stats_box,
    apply_publication_style,
    get_color_palette,
    require_matplotlib,
    save_figure,
)

LOGGER = get_logger("biopipe.scientific-figures")

try:
    import numpy as np
except ImportError:  # pragma: no cover
    np = None  # type: ignore[assignment]


# ---------------------------------------------------------------- 数据读取


def read_table(path: Path) -> list[dict[str, str]]:
    """读取 CSV 为字典列表（自动处理 BOM 与分隔符嗅探）。"""
    from biopipe_core import read_csv_safe

    if not path.exists():
        raise SkillError(f"输入文件不存在: {path}")
    rows = read_csv_safe(path)
    if not rows:
        raise SkillError(f"输入文件没有数据: {path}")
    return rows


def to_float(value: Optional[str], default: float = 0.0) -> float:
    """安全转 float。"""
    if value is None:
        return default
    text = str(value).strip()
    if not text:
        return default
    try:
        return float(text)
    except ValueError:
        return default


def pick_column(rows: list[dict[str, str]], candidates: list[str]) -> Optional[str]:
    """从候选列名中挑出第一个实际存在的列。"""
    if not rows:
        return None
    available = {k.lower(): k for k in rows[0].keys()}
    for candidate in candidates:
        key = available.get(candidate.lower())
        if key:
            return key
    return None


# ---------------------------------------------------------------- 统计工具


def roc_auc(fpr: list[float], tpr: list[float]) -> float:
    """梯形法计算 AUC（从曲线点直接积分）。"""
    if len(fpr) < 2:
        return 0.5
    order = sorted(range(len(fpr)), key=lambda i: fpr[i])
    area = 0.0
    for a, b in zip(order[:-1], order[1:]):
        width = fpr[b] - fpr[a]
        height = (tpr[a] + tpr[b]) / 2.0
        area += width * height
    return max(0.0, min(1.0, area))


def auc_z_test(auc1: float, auc2: float, n_samples: int) -> tuple[float, float]:
    """两个 AUC 差异的 Z 检验（Hanley-McNeil 标准误近似）。

    Returns:
        ``(z 统计量, 双尾 p 值)``
    """
    if n_samples <= 0:
        return 0.0, 1.0

    def hanley_mcneil_se(auc: float, n: int) -> float:
        q1 = auc / (2.0 - auc)
        q2 = 2.0 * auc * auc / (1.0 + auc)
        numerator = auc * (1.0 - auc) + (n - 1) * (q1 - auc * auc + q2 - (auc * auc) ** 2)
        denominator = n - 1
        value = numerator / denominator
        return math.sqrt(value) if value > 0 else 0.0

    se1 = hanley_mcneil_se(auc1, n_samples)
    se2 = hanley_mcneil_se(auc2, n_samples)
    se_diff = math.sqrt(se1**2 + se2**2)
    if se_diff == 0:
        return 0.0, 1.0
    z = (auc1 - auc2) / se_diff
    p = 2.0 * (1.0 - _normal_cdf(abs(z)))
    return z, p


def _normal_cdf(x: float) -> float:
    """标准正态分布 CDF（用 erf 近似，避免 scipy 依赖）。"""
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def bootstrap_auc_ci(
    fpr: list[float],
    tpr: list[float],
    n_bootstrap: int = 200,
    confidence: float = 0.95,
) -> tuple[Optional[float], Optional[float]]:
    """用曲线点重采样做 AUC 的自助法置信区间。

    对 (FPR, TPR) 点对做有放回重采样后重算 AUC，取分位数。
    """
    if len(fpr) < 5:
        return None, None
    import random

    rng = random.Random(42)
    n = len(fpr)
    aucs: list[float] = []
    for _ in range(n_bootstrap):
        idx = [rng.randrange(n) for _ in range(n)]
        sample_fpr = [fpr[i] for i in idx]
        sample_tpr = [tpr[i] for i in idx]
        value = roc_auc(sample_fpr, sample_tpr)
        if 0.0 <= value <= 1.0:
            aucs.append(value)
    if len(aucs) < 10:
        return None, None
    aucs.sort()
    alpha = 1.0 - confidence
    lower = aucs[int(alpha / 2 * len(aucs))]
    upper = aucs[min(len(aucs) - 1, int((1 - alpha / 2) * len(aucs)))]
    return lower, upper


# ---------------------------------------------------------------- 各图类型


def plot_roc(rows: list[dict[str, str]], args: argparse.Namespace) -> list[Path]:
    """绘制多模型 ROC 曲线。"""
    require_matplotlib()
    import matplotlib.pyplot as plt

    model_col = pick_column(rows, ["Model", "model", "method", "classifier"])
    fpr_col = pick_column(rows, ["FPR", "fpr", "false_positive_rate"])
    tpr_col = pick_column(rows, ["TPR", "tpr", "true_positive_rate", "sensitivity"])
    auc_col = pick_column(rows, ["AUC", "auc"])
    if not model_col or not fpr_col or not tpr_col:
        raise SkillError("ROC 输入需要 Model / FPR / TPR 列")

    # 按模型分组
    grouped: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        grouped.setdefault(row.get(model_col, "unknown"), []).append(row)

    curves: dict[str, tuple[list[float], list[float], float]] = {}
    for model, model_rows in grouped.items():
        model_rows = sorted(model_rows, key=lambda r: to_float(r.get(fpr_col)))
        fpr = [to_float(r.get(fpr_col)) for r in model_rows]
        tpr = [to_float(r.get(tpr_col)) for r in model_rows]
        if auc_col:
            auc = to_float(model_rows[0].get(auc_col), roc_auc(fpr, tpr))
        else:
            auc = roc_auc(fpr, tpr)
        curves[model] = (fpr, tpr, auc)

    fig, ax = plt.subplots(figsize=(7, 6))
    colors = get_color_palette(len(curves), "okabe_ito")

    best_model = max(curves, key=lambda m: curves[m][2])
    for idx, (model, (fpr, tpr, auc)) in enumerate(sorted(curves.items())):
        label = f"{model} (AUC = {auc:.2f})"
        ax.plot(fpr, tpr, color=colors[idx % len(colors)], linewidth=2, label=label)
        if model == best_model and args.fill_best:
            ax.fill_between(fpr, tpr, alpha=0.10, color=colors[idx % len(colors)])

    ax.plot([0, 1], [0, 1], "k--", linewidth=1, label="Random (AUC = 0.50)")
    ax.set_xlabel("False Positive Rate (FPR)", fontsize=11, fontweight="bold")
    ax.set_ylabel("True Positive Rate (TPR)", fontsize=11, fontweight="bold")
    ax.legend(loc="lower right", fontsize=9, frameon=False)
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # 统计信息
    stats_lines = [f"Best model: {best_model}", f"AUC = {curves[best_model][2]:.3f}"]
    if args.bootstrap_ci:
        lower, upper = bootstrap_auc_ci(*curves[best_model][:2])
        if lower is not None and upper is not None:
            stats_lines.append(f"95% CI: [{lower:.3f}, {upper:.3f}]")
    add_stats_box(ax, "\n".join(stats_lines), loc="upper left", fontsize=8)

    fig.tight_layout()
    saved = save_figure(fig, args.output, formats=args.formats, dpi=args.dpi)
    plt.close(fig)

    # AUC 差异显著性检验表
    if args.stats_output:
        models = sorted(curves)
        stat_rows = []
        for i in range(len(models)):
            for j in range(i + 1, len(models)):
                auc1 = curves[models[i]][2]
                auc2 = curves[models[j]][2]
                z, p = auc_z_test(auc1, auc2, args.n_samples)
                sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
                stat_rows.append(
                    {
                        "Model 1": models[i],
                        "Model 2": models[j],
                        "AUC 1": round(auc1, 4),
                        "AUC 2": round(auc2, 4),
                        "Difference": round(abs(auc1 - auc2), 4),
                        "Z-statistic": round(z, 4),
                        "P-value": round(p, 6),
                        "Significance": sig,
                    }
                )
        _write_csv_rows(stat_rows, args.stats_output)
        LOGGER.info("已写出统计检验: %s", args.stats_output)

    return saved


def plot_confusion(rows: list[dict[str, str]], args: argparse.Namespace) -> list[Path]:
    """绘制混淆矩阵热力图（多模型子图）。"""
    require_matplotlib()
    import matplotlib.pyplot as plt

    model_col = pick_column(rows, ["Model", "model", "method", "classifier"])
    if not model_col:
        raise SkillError("混淆矩阵输入需要 Model 列")
    needed = {
        "TN": pick_column(rows, ["TN", "tn", "true_negative"]),
        "FP": pick_column(rows, ["FP", "fp", "false_positive"]),
        "FN": pick_column(rows, ["FN", "fn", "false_negative"]),
        "TP": pick_column(rows, ["TP", "tp", "true_positive"]),
    }
    missing = [k for k, v in needed.items() if not v]
    if missing:
        raise SkillError(f"混淆矩阵输入缺少列: {', '.join(missing)}")

    matrices: list[tuple[str, list[list[int]]]] = []
    all_values: list[int] = []
    for row in rows:
        model = row.get(model_col, "model")
        cm = [
            [int(to_float(row.get(needed["TN"]))), int(to_float(row.get(needed["FP"])))],
            [int(to_float(row.get(needed["FN"]))), int(to_float(row.get(needed["TP"])))],
        ]
        matrices.append((model, cm))
        all_values.extend(v for row_vals in cm for v in row_vals)

    n_models = len(matrices)
    cols = min(3, n_models)
    rows_n = (n_models + cols - 1) // cols
    fig = plt.figure(figsize=(3.3 * cols, 2.9 * rows_n))
    gs = fig.add_gridspec(rows_n, cols, wspace=0.65, hspace=0.65)

    vmin = min(all_values)
    vmax = max(all_values)
    if vmax == vmin:
        vmax = vmin + 1

    for idx, (model, cm) in enumerate(matrices):
        ax = fig.add_subplot(gs[idx // cols, idx % cols])
        total = sum(sum(r) for r in cm) or 1
        x_edges = [0, 1, 2]
        y_edges = [0, 1, 2]
        ax.pcolormesh(x_edges, y_edges, cm, cmap="Blues", vmin=vmin, vmax=vmax,
                      edgecolors="white", linewidth=0.5)
        for i in range(2):
            for j in range(2):
                val = cm[i][j]
                pct = val / total * 100
                color = "white" if val > (vmin + vmax) / 2 else "black"
                ax.text(j + 0.5, i + 0.5, f"{val}\n({pct:.1f}%)",
                        ha="center", va="center", color=color, fontsize=10,
                        fontweight="bold")
        ax.set_xticks([0.5, 1.5])
        ax.set_yticks([0.5, 1.5])
        ax.set_xticklabels(["Pred Neg", "Pred Pos"], fontsize=9)
        ax.set_yticklabels(["Actual Neg", "Actual Pos"], fontsize=9)
        ax.set_title(model, fontsize=10)
        ax.set_xlim(0, 2)
        ax.set_ylim(2, 0)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    saved = save_figure(fig, args.output, formats=args.formats, dpi=args.dpi)
    plt.close(fig)
    return saved


def plot_heatmap(rows: list[dict[str, str]], args: argparse.Namespace) -> list[Path]:
    """绘制矩阵热图（对接结合能 / 相关性矩阵）。"""
    require_matplotlib()
    import matplotlib.pyplot as plt

    if not rows:
        raise SkillError("热图输入没有数据")

    header = list(rows[0].keys())
    # 首列若是标签列（文本而非数字），视为行标签
    first_col = header[0]
    row_labels: list[str] = []
    matrix: list[list[float]] = []
    label_mode = args.row_labels

    for row in rows:
        values: list[float] = []
        for key in header[1:] if label_mode != "none" else header:
            values.append(to_float(row.get(key), float("nan")))
        if label_mode == "none":
            row_labels.append("")
        else:
            row_labels.append(str(row.get(first_col, "")))
        matrix.append(values)

    col_labels = header[1:] if label_mode != "none" else header

    # 用非 NaN 值统计色阶
    flat = [v for r in matrix for v in r if not math.isnan(v)]
    if not flat:
        raise SkillError("热图输入没有数值")
    if args.percentile_clip:
        flat_sorted = sorted(flat)
        lo_idx = int(len(flat_sorted) * 0.02)
        hi_idx = int(len(flat_sorted) * 0.98)
        vmin, vmax = flat_sorted[lo_idx], flat_sorted[hi_idx]
    else:
        vmin, vmax = min(flat), max(flat)
    if vmax == vmin:
        vmax = vmin + 1.0

    fig, ax = plt.subplots(figsize=(args.width, args.height))
    cmap = args.cmap
    image = ax.imshow(matrix, cmap=cmap, aspect="auto", vmin=vmin, vmax=vmax)

    cbar = plt.colorbar(image, ax=ax, shrink=0.8, pad=0.02)
    cbar.set_label(args.colorbar_label, fontsize=11)

    tick_step = max(1, len(col_labels) // 20)
    x_idx = list(range(0, len(col_labels), tick_step))
    ax.set_xticks(x_idx)
    ax.set_xticklabels([col_labels[i] for i in x_idx], rotation=45, ha="right", fontsize=8)

    y_step = max(1, len(row_labels) // 25)
    y_idx = list(range(0, len(row_labels), y_step))
    ax.set_yticks(y_idx)
    ax.set_yticklabels([row_labels[i] for i in y_idx], fontsize=8)

    ax.set_xlabel(args.x_label, fontsize=11, fontweight="bold")
    ax.set_ylabel(args.y_label, fontsize=11, fontweight="bold")

    stats_text = (
        f"Mean: {statistics.mean(flat):.2f}\n"
        f"Std: {statistics.pstdev(flat):.2f}\n"
        f"Range: {min(flat):.1f} to {max(flat):.1f}"
    )
    add_stats_box(ax, stats_text, loc="upper right", fontsize=8)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    saved = save_figure(fig, args.output, formats=args.formats, dpi=args.dpi)
    plt.close(fig)
    return saved


def plot_plddt(rows: list[dict[str, str]], args: argparse.Namespace) -> list[Path]:
    """绘制 pLDDT 结构置信度直方图。"""
    require_matplotlib()
    import matplotlib.pyplot as plt

    col = pick_column(rows, ["Global pLDDT", "plddt", "pLDDT", "confidence", "score"])
    if not col:
        # 退化处理：第一列全是数字就用第一列
        first = list(rows[0].keys())[0]
        try:
            to_float(rows[0].get(first))
            col = first
        except (TypeError, ValueError):
            col = None
    if not col:
        raise SkillError("pLDDT 输入需要数值列（如 Global pLDDT）")

    values = [to_float(r.get(col)) for r in rows]
    values = [v for v in values if v > 0]
    if not values:
        raise SkillError("pLDDT 输入没有有效数值")

    fig, ax = plt.subplots(figsize=(8, 5))
    n, bins, patches = ax.hist(
        values, bins=args.bins, color="#88CCEE", edgecolor="white", linewidth=0.5, alpha=0.85
    )
    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    for center, patch in zip(bin_centers, patches):
        if center >= 90:
            color = SEMANTIC_COLORS["high"]
        elif center >= 70:
            color = SEMANTIC_COLORS["medium"]
        else:
            color = SEMANTIC_COLORS["low"]
        patch.set_facecolor(color)

    mean_v = statistics.mean(values)
    median_v = statistics.median(values)
    ax.axvline(mean_v, color=SEMANTIC_COLORS["accent"], linestyle="--", linewidth=1.5,
               ymax=0.95, label=f"Mean: {mean_v:.1f}")
    ax.axvline(median_v, color="#88CCEE", linestyle=":", linewidth=1.5,
               ymax=0.95, label=f"Median: {median_v:.1f}")

    total = len(values)
    high_pct = sum(1 for v in values if v >= 90) / total * 100
    med_pct = sum(1 for v in values if 70 <= v < 90) / total * 100
    low_pct = sum(1 for v in values if v < 70) / total * 100

    ax.set_xlabel("pLDDT Confidence Score", fontsize=12, fontweight="bold")
    ax.set_ylabel("Frequency", fontsize=12, fontweight="bold")
    ax.legend(loc="upper right", fontsize=9, frameon=False, bbox_to_anchor=(1.0, 1.07))

    stats_text = (
        f"Categories:\n"
        f"High (>90): {high_pct:.1f}%\n"
        f"Medium (70-90): {med_pct:.1f}%\n"
        f"Low (<70): {low_pct:.1f}%\n"
        f"N: {total:,}"
    )
    add_stats_box(ax, stats_text, loc="upper left", fontsize=8)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    saved = save_figure(fig, args.output, formats=args.formats, dpi=args.dpi)
    plt.close(fig)
    return saved


def plot_shap(rows: list[dict[str, str]], args: argparse.Namespace) -> list[Path]:
    """绘制 SHAP 摘要蜂群图 + 特征重要性条形图。"""
    require_matplotlib()
    import matplotlib.pyplot as plt

    # 识别数值特征列（排除 ID 类列）
    id_col = pick_column(rows, ["id", "ID", "index", "sample", "name"])
    feature_cols: list[str] = []
    for key in rows[0].keys():
        if key == id_col:
            continue
        if all(_is_number(r.get(key)) for r in rows):
            feature_cols.append(key)
    if not feature_cols:
        raise SkillError("SHAP 输入需要至少一个数值特征列")

    # 计算平均绝对 SHAP 值
    importance: list[tuple[str, float, float]] = []
    for feature in feature_cols:
        values = [to_float(r.get(feature)) for r in rows]
        mean_abs = sum(abs(v) for v in values) / len(values)
        if len(values) > 1:
            std_abs = statistics.pstdev([abs(v) for v in values])
        else:
            std_abs = 0.0
        importance.append((feature, mean_abs, std_abs))
    importance.sort(key=lambda item: item[1], reverse=True)

    top_features = importance[: args.top_n]
    global_max = max((abs(to_float(r.get(f))) for r in rows for f, _, _ in top_features), default=1.0)
    global_max = global_max or 1.0

    import random

    # ---- 蜂群图 ----
    fig, ax = plt.subplots(figsize=(10, max(6, 0.45 * len(top_features) + 2)))
    cmap = plt.get_cmap(args.shap_cmap)
    for idx, (feature, _, _) in enumerate(top_features):
        values = [to_float(r.get(feature)) for r in rows]
        rng = random.Random(42 + idx)
        jitter = [rng.uniform(-0.35, 0.35) for _ in values]
        ax.scatter(
            values,
            [idx + offset for offset in jitter],
            c=values,
            cmap=cmap,
            s=18,
            alpha=0.6,
            edgecolors="none",
            vmin=-global_max,
            vmax=global_max,
        )
    ax.set_yticks(range(len(top_features)))
    ax.set_yticklabels([f for f, _, _ in top_features], fontsize=8)
    ax.axvline(x=0, color="gray", linestyle="--", linewidth=0.8, alpha=0.7)
    ax.set_xlabel("SHAP Value (impact on model output)", fontsize=11, fontweight="bold")
    ax.set_ylabel("Features", fontsize=11, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    mappable = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=-global_max, vmax=global_max))
    mappable.set_array([])
    cbar = plt.colorbar(mappable, ax=ax, shrink=0.6, aspect=20, pad=0.02)
    cbar.set_label("Feature Value Effect", fontsize=10)

    fig.tight_layout()
    saved = save_figure(fig, f"{args.output}_beeswarm", formats=args.formats, dpi=args.dpi)
    plt.close(fig)

    # ---- 特征重要性条形图 ----
    fig2, ax2 = plt.subplots(figsize=(8, max(6, 0.4 * len(top_features) + 2)))
    positions = range(len(top_features))
    colors = get_color_palette(len(top_features))
    bars = ax2.barh(
        list(positions),
        [v for _, v, _ in top_features],
        color=colors,
        edgecolor="white",
        linewidth=0.5,
        xerr=[s for _, _, s in top_features],
        error_kw={"ecolor": "#666666", "capsize": 2, "linewidth": 1},
        height=0.7,
    )
    ax2.set_yticks(list(positions))
    ax2.set_yticklabels([f for f, _, _ in top_features], fontsize=8)
    ax2.invert_yaxis()
    ax2.set_xlabel("Mean |SHAP Value|", fontsize=11, fontweight="bold")
    ax2.set_ylabel("Features", fontsize=11, fontweight="bold")
    for bar, (_feature, value, std) in zip(bars, top_features):
        ax2.text(
            value + std + 0.005,
            bar.get_y() + bar.get_height() / 2,
            f"{value:.3f}",
            va="center",
            fontsize=7,
            color="#333333",
        )
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    fig2.tight_layout()
    saved += save_figure(fig2, f"{args.output}_importance", formats=args.formats, dpi=args.dpi)
    plt.close(fig2)

    # ---- 特征重要性表 ----
    if args.importance_output:
        table_rows = [
            {
                "feature": feature,
                "mean_abs_shap": round(value, 6),
                "std_abs_shap": round(std, 6),
                "rank": idx + 1,
            }
            for idx, (feature, value, std) in enumerate(importance)
        ]
        _write_csv_rows(table_rows, args.importance_output)
        LOGGER.info("已写出特征重要性表: %s", args.importance_output)

    return saved


def plot_class_distribution(rows: list[dict[str, str]], args: argparse.Namespace) -> list[Path]:
    """绘制门-纲双层级分类分布条形图。"""
    require_matplotlib()
    import matplotlib.pyplot as plt

    phylum_col = pick_column(rows, ["phylum", "Phylum", "门"])
    class_col = pick_column(rows, ["class", "Class", "纲"])
    count_col = pick_column(rows, ["count", "Count", "数量", "n"])
    kingdom_col = pick_column(rows, ["kingdom", "Kingdom", "界"])
    if not class_col or not count_col:
        raise SkillError("分类分布输入需要 class 与 count 列")

    aggregated: dict[tuple[str, str, str], int] = {}
    for row in rows:
        key = (
            row.get(kingdom_col, "") if kingdom_col else "",
            row.get(phylum_col, "") if phylum_col else "",
            row.get(class_col, ""),
        )
        aggregated[key] = aggregated.get(key, 0) + int(to_float(row.get(count_col)))

    items = sorted(aggregated.items(), key=lambda kv: kv[1], reverse=True)
    total = sum(v for _, v in items)
    top_items = items[: args.top_n]

    fig, ax = plt.subplots(figsize=(11, max(5, 0.38 * len(top_items) + 2)))
    colors = get_color_palette(len(top_items))
    y_positions = list(range(len(top_items)))
    values = [v for _, v in top_items]
    labels = [k[2] for k, _ in top_items]

    bars = ax.barh(y_positions, values, color=colors, edgecolor="white",
                   linewidth=0.5, height=0.7)
    ax.invert_yaxis()

    max_val = max(values) if values else 1
    for idx, (_bar, (key, value)) in enumerate(zip(bars, top_items)):
        kingdom, phylum, class_name = key
        label = f"{value:,} [{kingdom}-{phylum}]" if phylum else f"{value:,}"
        ax.text(value + max_val * 0.005, idx, label, va="center", ha="left",
                fontsize=9, color="#333333")

    ax.set_yticks(y_positions)
    ax.set_yticklabels([name[:28] + ("..." if len(name) > 28 else "") for name in labels],
                       fontsize=10)
    ax.set_xlabel("Count", fontsize=11, fontweight="bold")
    ax.set_ylabel("Class", fontsize=11, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.set_xlim(0, max_val * 1.55)

    if phylum_col:
        import matplotlib.patches as mpatches

        seen: list[str] = []
        handles = []
        for key, _ in top_items:
            phylum = key[1]
            if phylum and phylum not in seen:
                seen.append(phylum)
                handles.append(mpatches.Patch(color=PAPER_COLORS[seen.index(phylum) % len(PAPER_COLORS)],
                                              label=phylum))
        if handles:
            ax.legend(handles=handles, loc="lower right", fontsize=8, title="Phylum",
                      title_fontsize=9, frameon=True, facecolor="white", edgecolor="#cccccc")

    add_stats_box(ax, f"Total: {total:,}\nShown: {len(top_items)}", loc="lower left", fontsize=8)

    fig.tight_layout()
    saved = save_figure(fig, args.output, formats=args.formats, dpi=args.dpi)
    plt.close(fig)
    return saved


def plot_secondary_structure(rows: list[dict[str, str]], args: argparse.Namespace) -> list[Path]:
    """绘制二级结构分布柱状图。"""
    require_matplotlib()
    import matplotlib.pyplot as plt

    col = pick_column(rows, ["ss_name", "secondary_structure", "structure", "二级结构"])
    if not col:
        raise SkillError("二级结构输入需要 ss_name 列")

    counts: dict[str, int] = {}
    for row in rows:
        name = str(row.get(col, "") or "unknown")
        counts[name] = counts.get(name, 0) + 1
    ordered = sorted(counts.items(), key=lambda kv: kv[1], reverse=True)

    fig, ax = plt.subplots(figsize=(7, 5))
    colors = get_color_palette(len(ordered))
    x_positions = list(range(len(ordered)))
    values = [v for _, v in ordered]
    bars = ax.bar(x_positions, values, color=colors, edgecolor="white",
                   linewidth=0.5, width=0.6)
    ax.set_xticks(x_positions)
    ax.set_xticklabels([name for name, _ in ordered], fontsize=9, rotation=15, ha="right")

    total = sum(values) or 1
    for idx, (_bar, value) in enumerate(zip(bars, values)):
        ax.text(idx, value + max(values) * 0.02, f"{value:,}\n({value / total * 100:.1f}%)",
                ha="center", va="bottom", fontsize=8)

    ax.set_xlabel("Secondary Structure", fontsize=11, fontweight="bold")
    ax.set_ylabel("Residues", fontsize=11, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.text(0.5, -0.02, f"Total residues: {total:,}", ha="center", fontsize=9, style="italic")
    fig.tight_layout()
    saved = save_figure(fig, args.output, formats=args.formats, dpi=args.dpi)
    plt.close(fig)
    return saved


def plot_model_performance(rows: list[dict[str, str]], args: argparse.Namespace) -> list[Path]:
    """绘制模型多指标对比柱状图。"""
    require_matplotlib()
    import matplotlib.pyplot as plt

    model_col = pick_column(rows, ["Model", "model", "method"])
    if not model_col:
        raise SkillError("模型性能输入需要 Model 列")

    metric_candidates = ["Accuracy", "Precision", "Recall", "F1-Score", "ROC-AUC"]
    metrics: list[tuple[str, str]] = []
    for candidate in metric_candidates:
        col = pick_column(rows, [candidate])
        if col:
            metrics.append((candidate, col))
    if not metrics:
        raise SkillError("模型性能输入需要至少一个指标列（Accuracy / Precision / ...）")

    models = [r.get(model_col, "") for r in rows]
    x = list(range(len(models)))
    width = 0.8 / len(metrics)
    colors = get_color_palette(len(metrics))

    fig, ax = plt.subplots(figsize=(max(8, len(models) * 1.6), 5))
    for idx, (label, col) in enumerate(metrics):
        values = [to_float(r.get(col)) for r in rows]
        offset = (idx - len(metrics) / 2 + 0.5) * width
        bars = ax.bar([p + offset for p in x], values, width, label=label,
                      color=colors[idx], edgecolor="white", linewidth=0.3)
        for bar, value in zip(bars, values):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.005,
                    f"{value:.2f}", ha="center", va="bottom", fontsize=5)

    ax.set_xlabel("Model", fontsize=11, fontweight="bold")
    ax.set_ylabel("Score", fontsize=11, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(models, fontsize=9, rotation=15, ha="right")
    all_values = [to_float(r.get(col)) for _, col in metrics for r in rows]
    if all_values:
        lower = min(all_values)
        ax.set_ylim(max(0.0, lower - 0.1), min(1.02, max(all_values) + 0.03))
    ax.legend(loc="upper right", fontsize=7, ncol=len(metrics), frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    saved = save_figure(fig, args.output, formats=args.formats, dpi=args.dpi)
    plt.close(fig)
    return saved


def plot_feature_importance(rows: list[dict[str, str]], args: argparse.Namespace) -> list[Path]:
    """绘制特征重要性柱状图。"""
    require_matplotlib()
    import matplotlib.pyplot as plt

    feature_col = pick_column(rows, ["feature", "Feature", "name", "特征"])
    value_col = pick_column(rows, ["mean_abs_shap", "importance", "value", "重要性"])
    if not feature_col or not value_col:
        raise SkillError("特征重要性输入需要 feature 与 mean_abs_shap 列")
    std_col = pick_column(rows, ["std_abs_shap", "std", "error"])

    items = []
    for row in rows:
        items.append(
            (
                str(row.get(feature_col, "")),
                to_float(row.get(value_col)),
                to_float(row.get(std_col)) if std_col else 0.0,
            )
        )
    items.sort(key=lambda item: item[1], reverse=True)
    items = items[: args.top_n]

    fig, ax = plt.subplots(figsize=(8, max(5, 0.4 * len(items) + 2)))
    colors = get_color_palette(len(items))
    x_positions = list(range(len(items)))
    values = [v for _, v, _ in items]
    errors = [e for _, _, e in items] if std_col else None

    bars = ax.bar(
        x_positions, values, color=colors, edgecolor="white", linewidth=0.5,
        width=0.7, yerr=errors,
        error_kw={"ecolor": "gray", "capsize": 3, "linewidth": 1} if errors else None,
    )
    ax.set_xticks(x_positions)
    ax.set_xticklabels([f for f, _, _ in items], rotation=45, ha="right", fontsize=7)

    for idx, (_bar, (_feature, value, std)) in enumerate(zip(bars, items)):
        ax.text(idx, value + std + max(values) * 0.02, f"{value:.3f}",
                ha="center", va="bottom", fontsize=7)

    ax.set_ylabel("Mean |SHAP Value|", fontsize=11, fontweight="bold")
    ax.set_xlabel("Feature", fontsize=11, fontweight="bold")
    ax.set_ylim(0, max(values) * 1.2 if values else 1)
    if errors:
        ax.errorbar([], [], yerr=[], fmt="none", ecolor="gray", capsize=3,
                    linewidth=1, label="Error bar: ±SD")
        ax.legend(loc="upper right", fontsize=8, framealpha=0.9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    saved = save_figure(fig, args.output, formats=args.formats, dpi=args.dpi)
    plt.close(fig)
    return saved


def _is_number(value: Optional[str]) -> bool:
    """判断字符串是否可解析为数字。"""
    if value is None:
        return False
    text = str(value).strip()
    if not text:
        return False
    try:
        float(text)
        return True
    except ValueError:
        return False


def _write_csv_rows(rows: list[dict], path: str) -> None:
    """写出 CSV（UTF-8 BOM）。"""
    from biopipe_core import write_csv

    if not rows:
        return
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    write_csv(rows, out, fieldnames=list(rows[0].keys()))


# ---------------------------------------------------------------- CLI


SUBCOMMANDS = {
    "roc": plot_roc,
    "confusion": plot_confusion,
    "heatmap": plot_heatmap,
    "plddt": plot_plddt,
    "shap": plot_shap,
    "class-distribution": plot_class_distribution,
    "secondary-structure": plot_secondary_structure,
    "model-performance": plot_model_performance,
    "feature-importance": plot_feature_importance,
}


def build_parser() -> argparse.ArgumentParser:
    """构建命令行解析器。"""
    parser = cli_parser(
        prog="scientific-figures",
        description="出版级科学图表生成（ROC / 混淆矩阵 / 热图 / pLDDT / SHAP / 分类分布等）",
    )
    parser.add_argument(
        "subcommand",
        choices=sorted(SUBCOMMANDS),
        help="图表类型",
    )
    parser.add_argument("--input", required=True, help="输入 CSV 路径")
    parser.add_argument("--output", required=True, help="输出路径（不含扩展名或含扩展名）")
    parser.add_argument("--style", default="nature",
                        choices=["default", "nature", "cell", "science"],
                        help="期刊样式预设")
    parser.add_argument("--formats", nargs="+", default=["png", "pdf"], help="输出格式")
    parser.add_argument("--dpi", type=int, default=300, help="输出分辨率")
    parser.add_argument("--title", default="", help="图表标题")
    # 通用调节
    parser.add_argument("--top-n", type=int, default=15, help="显示 Top N（SHAP / 分类分布）")
    parser.add_argument("--bins", type=int, default=50, help="直方图分箱数（pLDDT）")
    parser.add_argument("--cmap", default="coolwarm_r", help="热图色阶")
    parser.add_argument("--colorbar-label", default="Value", help="热图颜色条标签")
    parser.add_argument("--x-label", default="", help="热图 X 轴标签")
    parser.add_argument("--y-label", default="", help="热图 Y 轴标签")
    parser.add_argument("--width", type=float, default=12.0, help="热图宽度")
    parser.add_argument("--height", type=float, default=10.0, help="热图高度")
    parser.add_argument("--row-labels", default="auto", choices=["auto", "none"],
                        help="热图首列是否作为行标签")
    parser.add_argument("--percentile-clip", action="store_true",
                        help="热图色阶按 2%%/98%% 分位裁剪（抑制离群值）")
    parser.add_argument("--fill-best", action="store_true", help="ROC 图填充最佳模型曲线下方")
    parser.add_argument("--bootstrap-ci", action="store_true", help="ROC 图计算自助法置信区间")
    parser.add_argument("--n-samples", type=int, default=500, help="AUC Z 检验的样本量假设")
    parser.add_argument("--stats-output", help="ROC 统计检验结果 CSV 输出路径")
    parser.add_argument("--importance-output", help="SHAP 特征重要性表 CSV 输出路径")
    parser.add_argument("--shap-cmap", default="coolwarm", help="SHAP 图色阶")
    return parser


def main(args: argparse.Namespace) -> int:
    """主流程。"""
    if not MATPLOTLIB_AVAILABLE:
        raise FigureStyleError(
            "matplotlib 未安装，无法生成图表。请先安装：\n"
            "    pip install matplotlib\n"
            "    （或安装本技能的依赖：pip install -r skills/scientific-figures/requirements.txt）"
        )

    config = get_config(args.config)
    input_path = resolve_path(args.input)
    output_path = resolve_path(args.output)
    del output_path, config

    apply_publication_style(args.style)

    rows = read_table(input_path)
    LOGGER.info("读取 %s: %d 行, %d 列", input_path, len(rows), len(rows[0]))

    handler = SUBCOMMANDS[args.subcommand]
    saved = handler(rows, args)

    for path in saved:
        LOGGER.info("已保存图片: %s", path)

    LOGGER.info("图表生成完成（子命令: %s, 共 %d 个文件）", args.subcommand, len(saved))
    return 0


if __name__ == "__main__":
    sys.exit(
        run_cli(
            main,
            parser=build_parser(),
            skill="scientific-figures",
        )
    )
