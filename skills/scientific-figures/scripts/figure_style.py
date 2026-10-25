#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""出版级科学绘图样式模块。

移植自原仓库 ``data_visualization/scripts/common_style.py``，并做了以下增强：

- 统一的调色板（Okabe-Ito 色盲友好 / Paul Tol muted / 项目主色）
- 一致的 rcParams（Arial 优先字体、去除 top/right spines、关闭网格、
  300/600 dpi 保存、tight bbox）
- 期刊预设（``default`` / ``nature`` / ``cell`` / ``science``）
- 统一保存接口 :func:`save_figure`（PNG + PDF 双格式）
- 单图与多面板的画布创建接口
- matplotlib 缺失时的优雅降级（抛 :class:`FigureStyleError` 而非 ImportError）

用法::

    from figure_style import apply_publication_style, save_figure, get_color_palette

    apply_publication_style("nature")
    fig, ax = setup_single_figure(3.5, 2.8)
    ax.plot(x, y)
    save_figure(fig, "figures/fig1", formats=["png", "pdf"])
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Sequence

try:  # matplotlib 缺失时优雅降级
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    MATPLOTLIB_AVAILABLE = True
    _IMPORT_ERROR: str | None = None
except ImportError as exc:  # pragma: no cover - 依赖缺失路径
    MATPLOTLIB_AVAILABLE = False
    _IMPORT_ERROR = str(exc)
    mpl = None  # type: ignore[assignment]
    plt = None  # type: ignore[assignment]


class FigureStyleError(RuntimeError):
    """绘图相关错误（依赖缺失、数据不合法等）。"""


def require_matplotlib() -> None:
    """确认 matplotlib 可用，否则抛出带清晰提示的异常。"""
    if not MATPLOTLIB_AVAILABLE:
        raise FigureStyleError(
            "matplotlib 未安装，无法绘图。请先安装：\n"
            "    pip install matplotlib\n"
            f"原始错误: {_IMPORT_ERROR}"
        )


# ---------------------------------------------------------------- 调色板

# Okabe-Ito 色盲友好调色板（分类数据首选）
OKABE_ITO_COLORS = [
    "#E69F00",  # Orange
    "#56B4E9",  # Sky Blue
    "#009E73",  # Bluish Green
    "#F0E442",  # Yellow
    "#0072B2",  # Blue
    "#D55E00",  # Vermillion
    "#CC79A7",  # Reddish Purple
    "#000000",  # Black
]

# Paul Tol muted（最多 9 类）
TOL_MUTED = [
    "#332288", "#88CCEE", "#44AA99", "#117733",
    "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499",
]

# 项目主色
PRIMARY_COLORS = [
    "#88CCEE",  # Light Blue
    "#44AA99",  # Teal
    "#DDCC77",  # Gold
    "#CC6677",  # Rose
]

# 扩展色
EXTENDED_COLORS = [
    "#AA4499",  # Magenta
    "#882255",  # Purple
    "#999933",  # Olive
    "#117733",  # Dark Green
    "#88CCEE",
    "#44AA99",
    "#DDCC77",
    "#CC6677",
]

PAPER_COLORS = PRIMARY_COLORS + EXTENDED_COLORS

# 高对比（仅 3 类）
TOL_HIGH_CONTRAST = ["#004488", "#DDAA33", "#BB5566"]

# 二级结构配色
STRUCTURE_COLORS = {
    "Alpha helix": "#44AA99",
    "Beta sheet": "#88CCEE",
    "Loop/Disorder": "#999933",
}

# 语义色（用于分类 / 连续色阶）
SEMANTIC_COLORS = {
    "high": "#44AA99",
    "medium": "#88CCEE",
    "low": "#DDCC77",
    "negative": "#0072B2",
    "positive": "#D55E00",
    "accent": "#CC6677",
}

PALETTES = {
    "okabe_ito": OKABE_ITO_COLORS,
    "tol_muted": TOL_MUTED,
    "high_contrast": TOL_HIGH_CONTRAST,
    "paper": PAPER_COLORS,
    "primary": PRIMARY_COLORS,
}


def get_color_palette(n_colors: int, palette: str = "paper") -> list[str]:
    """按数量取调色板。

    Args:
        n_colors: 需要的颜色数量。
        palette: 调色板名（``paper`` / ``okabe_ito`` / ``tol_muted`` / ``high_contrast`` / ``primary``）。

    Returns:
        长度至少为 ``n_colors`` 的颜色列表。
    """
    selected = PALETTES.get(palette, PAPER_COLORS)
    if not selected:
        selected = PAPER_COLORS
    if n_colors <= len(selected):
        return selected[:n_colors]
    repeats = n_colors // len(selected) + 1
    return (selected * repeats)[:n_colors]


def sequential_cmap(name: str = "viridis"):
    """获取感知均匀的连续色阶。

    优先使用 matplotlib 内置；若不可用则回退到 viridis 的近似实现。
    """
    require_matplotlib()
    try:
        return plt.get_cmap(name)
    except (ValueError, AttributeError):
        return plt.get_cmap("viridis")


# ---------------------------------------------------------------- 样式


def get_base_style() -> dict[str, Any]:
    """返回出版级基础 rcParams。"""
    require_matplotlib()
    return {
        # 画布
        "figure.dpi": 100,
        "figure.facecolor": "white",
        "figure.edgecolor": "white",
        "figure.autolayout": False,
        "figure.constrained_layout.use": False,
        # 字体
        "font.size": 10,
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "mathtext.fontset": "dejavusans",
        # 坐标轴
        "axes.linewidth": 1.0,
        "axes.labelsize": 11,
        "axes.titlesize": 12,
        "axes.labelweight": "normal",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.spines.left": True,
        "axes.spines.bottom": True,
        "axes.edgecolor": "#333333",
        "axes.labelcolor": "#333333",
        "axes.axisbelow": True,
        "axes.prop_cycle": mpl.cycler(color=PAPER_COLORS),
        # 网格（默认关闭，保持简洁）
        "axes.grid": False,
        "grid.alpha": 0.0,
        "grid.linewidth": 0.0,
        "grid.linestyle": "none",
        # 刻度
        "xtick.major.size": 5,
        "xtick.minor.size": 2.5,
        "xtick.major.width": 0.8,
        "xtick.minor.width": 0.6,
        "xtick.labelsize": 9,
        "xtick.direction": "out",
        "xtick.color": "#333333",
        "ytick.major.size": 5,
        "ytick.minor.size": 2.5,
        "ytick.major.width": 0.8,
        "ytick.minor.width": 0.6,
        "ytick.labelsize": 9,
        "ytick.direction": "out",
        "ytick.color": "#333333",
        # 线条
        "lines.linewidth": 2.0,
        "lines.markersize": 6,
        "lines.markeredgewidth": 0.8,
        "lines.solid_capstyle": "round",
        # 图例
        "legend.fontsize": 9,
        "legend.frameon": False,
        "legend.framealpha": 0.0,
        "legend.loc": "best",
        "legend.edgecolor": "none",
        "legend.handlelength": 2.0,
        # 保存
        "savefig.dpi": 300,
        "savefig.format": "pdf",
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.1,
        "savefig.transparent": False,
        "savefig.facecolor": "white",
        "savefig.edgecolor": "white",
        # 图像
        "image.cmap": "viridis",
        "image.aspect": "auto",
        # 误差棒
        "errorbar.capsize": 4,
    }


# 期刊预设
JOURNAL_PRESETS: dict[str, dict[str, Any]] = {
    "nature": {
        "font.size": 7,
        "axes.labelsize": 8,
        "axes.titlesize": 9,
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,
        "legend.fontsize": 6,
        "savefig.dpi": 600,
        "lines.linewidth": 1.5,
        "lines.markersize": 4,
    },
    "cell": {
        "font.size": 8,
        "axes.labelsize": 9,
        "axes.titlesize": 10,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "savefig.dpi": 600,
        "lines.linewidth": 1.8,
    },
    "science": {
        "font.size": 9,
        "axes.labelsize": 10,
        "axes.titlesize": 11,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "savefig.dpi": 600,
    },
    "default": {},
}


def apply_publication_style(style_name: str = "default") -> None:
    """应用出版级样式。

    Args:
        style_name: ``default`` / ``nature`` / ``cell`` / ``science``。
    """
    require_matplotlib()
    style = get_base_style()
    preset = JOURNAL_PRESETS.get(style_name)
    if preset is None:
        raise FigureStyleError(
            f"未知样式 {style_name!r}，可选: {', '.join(sorted(JOURNAL_PRESETS))}"
        )
    style.update(preset)
    plt.rcParams.update(style)


def save_figure(
    fig: Any,
    filename: str | Path,
    formats: Sequence[str] = ("pdf", "png"),
    dpi: int = 300,
    **kwargs: Any,
) -> list[Path]:
    """保存出版级图片（默认 PNG + PDF 双格式）。

    Args:
        fig: matplotlib Figure。
        filename: 输出路径前缀（不含扩展名）或完整路径。
        formats: 输出格式列表。
        dpi: 分辨率。
        **kwargs: 透传给 ``fig.savefig``。

    Returns:
        成功写出的文件路径列表。
    """
    require_matplotlib()
    path = Path(filename)
    if path.suffix:
        base_name = path.stem
        formats = (path.suffix.lstrip("."),)
    else:
        base_name = path.name
    output_dir = path.parent if str(path.parent) else Path(".")
    output_dir.mkdir(parents=True, exist_ok=True)

    saved: list[Path] = []
    for fmt in formats:
        output_file = output_dir / f"{base_name}.{fmt}"
        save_kwargs: dict[str, Any] = {
            "dpi": dpi,
            "bbox_inches": "tight",
            "pad_inches": 0.1,
            "facecolor": "white",
            "edgecolor": "none",
            "transparent": False,
            "format": fmt,
        }
        save_kwargs.update(kwargs)
        try:
            fig.savefig(output_file, **save_kwargs)
            saved.append(output_file)
        except (OSError, ValueError) as exc:
            raise FigureStyleError(f"保存图片失败 {output_file}: {exc}") from exc
    return saved


def setup_single_figure(width: float = 3.5, height: float = 2.8) -> tuple[Any, Any]:
    """创建单幅标准出版尺寸画布。"""
    require_matplotlib()
    apply_publication_style("default")
    fig, ax = plt.subplots(figsize=(width, height))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    return fig, ax


def setup_multi_panel(
    nrows: int = 1,
    ncols: int = 1,
    width: float = 7.0,
    height: float | None = None,
    wspace: float = 0.3,
    hspace: float = 0.3,
) -> tuple[Any, Any]:
    """创建多面板画布。"""
    require_matplotlib()
    apply_publication_style("default")
    if height is None:
        height = width * 0.6

    fig, axes = plt.subplots(nrows, ncols, figsize=(width, height))
    if nrows == 1 and ncols == 1:
        axes = mpl_to_array(axes)
    else:
        axes = mpl_to_array(axes)

    for ax in axes.flat:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    fig.tight_layout()
    fig.subplots_adjust(wspace=wspace, hspace=hspace)
    return fig, axes


def mpl_to_array(axes: Any) -> Any:
    """把单个 Axes 包装成 2D 数组，统一遍历接口。"""
    import numpy as np

    return np.array([axes])


def add_stats_box(
    ax: Any,
    text: str,
    loc: str = "upper right",
    fontsize: float = 8,
) -> None:
    """在坐标轴上叠加统计信息文本框。"""
    positions = {
        "upper right": (0.98, 0.98, "right", "top"),
        "upper left": (0.02, 0.98, "left", "top"),
        "lower right": (0.98, 0.02, "right", "bottom"),
        "lower left": (0.02, 0.02, "left", "bottom"),
    }
    x, y, ha, va = positions.get(loc, positions["upper right"])
    ax.text(
        x,
        y,
        text,
        transform=ax.transAxes,
        fontsize=fontsize,
        verticalalignment=va,
        horizontalalignment=ha,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.9, edgecolor="none"),
    )


def label_bars(ax: Any, bars: Any, values: Sequence[float], fmt: str = "{:.0f}") -> None:
    """在柱状图末端标注数值。"""
    for bar, value in zip(bars, values):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height(),
            fmt.format(value),
            ha="center",
            va="bottom",
            fontsize=8,
            color="#333333",
        )
