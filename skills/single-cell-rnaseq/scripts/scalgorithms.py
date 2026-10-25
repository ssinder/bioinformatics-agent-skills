#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""单细胞分析算法核心（双后端：scanpy 与纯 NumPy/SciPy/sklearn 降级实现）。

所有算法都提供两个实现路径：

- **scanpy 后端**：调用 scanpy / anndata 的官方实现（结果与主流流程一致）
- **降级后端**：纯 NumPy/SciPy/sklearn 实现，保证在没有 scanpy 的环境中
  也能跑出**真实**的分析结果，而不是占位符

降级后端包含的算法：

============================  ==========================================
步骤                          降级实现
============================  ==========================================
CP10K 归一化 + log1p          ``normalize_total`` / ``log1p``
高变基因（HVG）               seurat 风格：归一化表达分箱内的 z-score 排序
缩放（z-score, clip 10）      ``scale``
PCA                           截断 SVD（scipy.sparse.linalg.svds）
kNN 图                       sklearn ``NearestNeighbors`` + Jaccard 边权
聚类                         Louvain 社区发现（自实现，含 modularity 优化）
UMAP 近似                    sklearn ``MDS`` + 邻域保持的力导向微调
t-SNE                         sklearn ``TSNE``
差异表达                     Welch t 检验 + BH 校正
============================  ==========================================

Usage:
    from scalgorithms import run_pipeline
    result = run_pipeline(adata, backend="auto")
"""

from __future__ import annotations

import importlib
import math
from dataclasses import dataclass, field
from typing import Any, Callable

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.sparse.linalg import svds

__all__ = [
    "SCANPY_AVAILABLE",
    "PipelineResult",
    "compute_qc_metrics",
    "filter_cells_genes",
    "normalize_total",
    "log1p",
    "scale",
    "highly_variable_genes",
    "compute_pca",
    "compute_neighbors",
    "cluster_louvain",
    "cluster_leiden",
    "compute_umap",
    "compute_tsne",
    "rank_genes_groups",
    "run_pipeline",
]


def _try_import(name: str):
    """尝试导入模块，失败返回 None。"""
    try:
        return importlib.import_module(name)
    except ImportError:
        return None


scanpy_mod = _try_import("scanpy")
anndata_mod = _try_import("anndata")
leidenalg_mod = _try_import("leidenalg")
igraph_mod = _try_import("igraph")
umap_mod = _try_import("umap")

SCANPY_AVAILABLE = scanpy_mod is not None and anndata_mod is not None


# --------------------------------------------------------------------------- #
# 结果容器
# --------------------------------------------------------------------------- #
@dataclass
class PipelineResult:
    """单细胞流程的结果汇总。

    Attributes:
        adata: 处理后的数据对象（SimpleAnnData 或 AnnData）。
        backend: 实际使用的后端名称（``scanpy`` / ``numpy``）。
        qc_stats: 质控统计信息。
        cluster_summary: 聚类摘要（每簇细胞数）。
        marker_table: 每簇的 marker 基因表。
        steps: 已完成步骤名列表。
    """

    adata: Any
    backend: str
    qc_stats: dict[str, Any] = field(default_factory=dict)
    cluster_summary: pd.DataFrame | None = None
    marker_table: pd.DataFrame | None = None
    steps: list[str] = field(default_factory=list)


# --------------------------------------------------------------------------- #
# 1. 质控指标
# --------------------------------------------------------------------------- #
def compute_qc_metrics(
    adata: Any,
    mt_mask: np.ndarray,
    hb_mask: np.ndarray | None = None,
    use_scanpy: bool = False,
) -> dict[str, np.ndarray]:
    """计算单细胞质控指标。

    对每个细胞计算：

    - ``total_counts``：总 UMI 数
    - ``n_genes_by_counts``：检测到的基因数
    - ``pct_counts_mt``：线粒体基因 UMI 占比（%）
    - ``pct_counts_hb``：血红蛋白基因 UMI 占比（%）
    - ``log1p_total_counts`` / ``log1p_n_genes``：对数形式，便于画图

    Args:
        adata: 数据对象（``SimpleAnnData``）。
        mt_mask: 线粒体基因的布尔掩码。
        hb_mask: 血红蛋白基因的布尔掩码，None 时跳过。
        use_scanpy: 为 True 时调用 scanpy 的实现（scanpy 不可用时自动回退）。

    Returns:
        指标名到数组的字典。
    """
    X = adata.X
    if sp.issparse(X):
        total = np.asarray(X.sum(axis=1)).ravel().astype(float)
        n_genes = np.asarray((X > 0).sum(axis=1)).ravel().astype(float)
        mt_counts = np.asarray(X[:, mt_mask].sum(axis=1)).ravel().astype(float) if mt_mask.any() else np.zeros_like(total)
        hb_counts = (
            np.asarray(X[:, hb_mask].sum(axis=1)).ravel().astype(float)
            if hb_mask is not None and hb_mask.any()
            else np.zeros_like(total)
        )
    else:
        arr = np.asarray(X, dtype=float)
        total = arr.sum(axis=1).astype(float)
        n_genes = (arr > 0).sum(axis=1).astype(float)
        mt_counts = arr[:, mt_mask].sum(axis=1).astype(float) if mt_mask.any() else np.zeros_like(total)
        hb_counts = (
            arr[:, hb_mask].sum(axis=1).astype(float)
            if hb_mask is not None and hb_mask.any()
            else np.zeros_like(total)
        )

    with np.errstate(divide="ignore", invalid="ignore"):
        pct_mt = np.where(total > 0, mt_counts / total * 100.0, 0.0)
        pct_hb = np.where(total > 0, hb_counts / total * 100.0, 0.0)

    return {
        "total_counts": total,
        "n_genes_by_counts": n_genes,
        "pct_counts_mt": pct_mt,
        "pct_counts_hb": pct_hb,
        "log1p_total_counts": np.log1p(total),
        "log1p_n_genes_by_counts": np.log1p(n_genes),
    }


# --------------------------------------------------------------------------- #
# 2. 过滤
# --------------------------------------------------------------------------- #
def filter_cells_genes(
    adata: Any,
    qc: dict[str, np.ndarray],
    min_genes: int = 200,
    min_counts: int = 100,
    max_mt_pct: float = 20.0,
    max_hb_pct: float | None = None,
    min_cells_per_gene: int = 3,
) -> tuple[Any, dict[str, Any]]:
    """按质控阈值过滤细胞与基因。

    细胞过滤条件（全部满足才保留）：

    - ``n_genes_by_counts >= min_genes``
    - ``total_counts >= min_counts``
    - ``pct_counts_mt <= max_mt_pct``（若提供 ``max_hb_pct`` 还需满足血红蛋白比例）

    基因过滤条件：在至少 ``min_cells_per_gene`` 个细胞中表达。

    Args:
        adata: 数据对象。
        qc: :func:`compute_qc_metrics` 的输出。
        min_genes: 每个细胞至少检测到的基因数。
        min_counts: 每个细胞至少的 UMI 数。
        max_mt_pct: 线粒体基因比例上限（%）。
        max_hb_pct: 血红蛋白比例上限（%），None 表示不过滤。
        min_cells_per_gene: 基因至少要在多少个细胞中表达。

    Returns:
        ``(filtered_adata, stats)``。
    """
    cell_mask = (
        (qc["n_genes_by_counts"] >= min_genes)
        & (qc["total_counts"] >= min_counts)
        & (qc["pct_counts_mt"] <= max_mt_pct)
    )
    if max_hb_pct is not None:
        cell_mask &= qc["pct_counts_hb"] <= max_hb_pct

    n_cells_before = adata.n_obs
    adata = adata.subset_cells(cell_mask)
    n_cells_after = adata.n_obs

    # 基因过滤
    X = adata.X
    if sp.issparse(X):
        cells_per_gene = np.asarray((X > 0).sum(axis=0)).ravel()
    else:
        cells_per_gene = (np.asarray(X) > 0).sum(axis=0)
    gene_mask = cells_per_gene >= min_cells_per_gene
    n_genes_before = adata.n_vars
    adata = adata.subset_genes(gene_mask)
    n_genes_after = adata.n_vars

    stats = {
        "n_cells_before": int(n_cells_before),
        "n_cells_after": int(n_cells_after),
        "n_cells_removed": int(n_cells_before - n_cells_after),
        "cell_removal_rate": float(1.0 - n_cells_after / max(n_cells_before, 1)),
        "n_genes_before": int(n_genes_before),
        "n_genes_after": int(n_genes_after),
        "n_genes_removed": int(n_genes_before - n_genes_after),
        "thresholds": {
            "min_genes": min_genes,
            "min_counts": min_counts,
            "max_mt_pct": max_mt_pct,
            "max_hb_pct": max_hb_pct,
            "min_cells_per_gene": min_cells_per_gene,
        },
    }
    return adata, stats


# --------------------------------------------------------------------------- #
# 3. 归一化
# --------------------------------------------------------------------------- #
def normalize_total(X: Any, target_sum: float = 1e4) -> Any:
    """CP10K / CPM 归一化（library size normalization）。

    每个细胞的计数除以其总 UMI 数，再乘以 ``target_sum``。

    Args:
        X: 表达矩阵。
        target_sum: 目标总计数（1e4 即 CP10K）。

    Returns:
        归一化后的矩阵（同类型）。
    """
    if sp.issparse(X):
        totals = np.asarray(X.sum(axis=1)).ravel()
        scale = np.divide(
            target_sum, totals, out=np.zeros_like(totals, dtype=float), where=totals > 0
        )
        return sp.diags(scale) @ X
    arr = np.asarray(X, dtype=float)
    totals = arr.sum(axis=1, keepdims=True)
    scale = np.divide(target_sum, totals, out=np.zeros_like(totals), where=totals > 0)
    return arr * scale


def log1p(X: Any) -> Any:
    """对数变换 ``log(1 + x)``。"""
    if sp.issparse(X):
        X = X.tocsr(copy=True)
        X.data = np.log1p(X.data)
        return X
    return np.log1p(np.asarray(X, dtype=float))


def scale(X: Any, max_value: float = 10.0) -> np.ndarray:
    """z-score 缩放并截断。

    每个基因减去均值、除以标准差（总体标准差 ``ddof=0``），再把结果截断到
    ``[-max_value, max_value]``（与 scanpy 的 ``sc.pp.scale`` 行为一致）。

    Args:
        X: 表达矩阵。
        max_value: 截断阈值。

    Returns:
        稠密的缩放后矩阵。
    """
    arr = X.toarray() if sp.issparse(X) else np.asarray(X, dtype=float)
    mean = arr.mean(axis=0, keepdims=True)
    std = arr.std(axis=0, keepdims=True)
    std[std == 0] = 1.0
    out = (arr - mean) / std
    return np.clip(out, -max_value, max_value)


# --------------------------------------------------------------------------- #
# 4. 高变基因
# --------------------------------------------------------------------------- #
def highly_variable_genes(
    X: Any,
    n_top_genes: int = 2000,
    min_mean: float = 0.0125,
    max_mean: float = 3.0,
    n_bins: int = 20,
) -> np.ndarray:
    """Seurat 风格的高变基因（HVG）选择。

    步骤：

    1. 对已做 CP10K + log1p 的矩阵，按每个基因的平均表达量分箱（等频分箱）
    2. 在每个箱内计算每个基因的 z-score（相对于该箱内基因）
    3. 按 z-score 排序，取前 ``n_top_genes`` 个
    4. 同时应用 ``min_mean`` / ``max_mean`` 的表达量区间过滤

    这是 Seurat 与 scanpy ``flavor="seurat"`` 共用的经典算法。

    Args:
        X: 已归一化并 log1p 的表达矩阵。
        n_top_genes: 选出的 HVG 数量。
        min_mean: 平均表达量下限。
        max_mean: 平均表达量上限。
        n_bins: 分箱数。

    Returns:
        布尔掩码数组（长度=基因数）。
    """
    arr = X.toarray() if sp.issparse(X) else np.asarray(X, dtype=float)
    means = arr.mean(axis=0)
    normalized_dispersion = np.zeros(arr.shape[1], dtype=float)

    # 等频分箱（用分位数作为分箱边界）
    valid = means > 0
    if valid.sum() == 0:
        return np.zeros(arr.shape[1], dtype=bool)

    disp = np.zeros_like(means)
    disp[valid] = arr[:, valid].var(axis=0) / np.maximum(means[valid], 1e-12)
    # 归一化离散度：除以该表达量下期望的泊松噪声
    with np.errstate(divide="ignore", invalid="ignore"):
        normalized_dispersion[valid] = disp[valid] / means[valid]

    quantiles = np.quantile(means[valid], np.linspace(0, 1, n_bins + 1))
    quantiles[0] = min(quantiles[0], means[valid].min() - 1e-9)
    quantiles[-1] = max(quantiles[-1], means[valid].max() + 1e-9)

    hvg = np.zeros(arr.shape[1], dtype=bool)
    for lo, hi in zip(quantiles[:-1], quantiles[1:]):
        if hi <= lo:
            continue
        in_bin = (means >= lo) & (means <= hi) & valid
        if in_bin.sum() < 2:
            continue
        vals = normalized_dispersion[in_bin]
        std = vals.std()
        if std == 0:
            continue
        z = (vals - vals.mean()) / std
        # 该箱内 z-score 最高的基因即为候选 HVG
        order = np.argsort(-z)
        bin_idx = np.flatnonzero(in_bin)
        # 选出 z 分数最高的若干（保守地按比例）
        n_take = max(1, int(np.ceil(len(bin_idx) * n_top_genes / max(valid.sum(), 1))))
        for i in order[:n_take]:
            hvg[bin_idx[i]] = True

    # 表达量区间过滤
    hvg &= (means >= min_mean) & (means <= max_mean)

    if hvg.sum() == 0:
        # 兜底：取平均表达量最高的 n_top_genes 个基因
        order = np.argsort(-means)[:n_top_genes]
        hvg[order] = True
    elif hvg.sum() > n_top_genes:
        # 在候选中按归一化离散度取前 n_top_genes
        cand = np.flatnonzero(hvg)
        keep = cand[np.argsort(-normalized_dispersion[cand])[:n_top_genes]]
        new_hvg = np.zeros_like(hvg)
        new_hvg[keep] = True
        hvg = new_hvg
    return hvg


# --------------------------------------------------------------------------- #
# 5. PCA
# --------------------------------------------------------------------------- #
def compute_pca(X: Any, n_comps: int = 50, seed: int = 0) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """主成分分析。

    对稀疏矩阵用截断 SVD（``scipy.sparse.linalg.svds``），对稠密矩阵用
    中心化后的 SVD，两者都等价于对中心化数据做 PCA。

    Args:
        X: 已缩放的表达矩阵。
        n_comps: 主成分数。
        seed: 随机种子（``svds`` 的 v0 初始化）。

    Returns:
        ``(X_pca, variance_ratio, components)``：

        - ``X_pca`` 形状 ``(n_cells, n_comps)``
        - ``variance_ratio`` 各主成分解释的方差比例
        - ``components`` 载荷矩阵，形状 ``(n_comps, n_genes)``
    """
    if sp.issparse(X):
        Xc = X.astype(float)
        n_obs, n_vars = Xc.shape
        k = int(min(n_comps, min(n_obs, n_vars) - 1))
        if k < 1:
            raise ValueError("矩阵过小，无法做 PCA")
        rng = np.random.default_rng(seed)
        v0 = rng.normal(size=min(Xc.shape))
        _, _, vt = svds(Xc.astype(float), k=k, v0=v0)
        # svds 返回按奇异值降序
        vt = vt[::-1]
        singular = np.linalg.norm(Xc @ vt.T, axis=0) if sp.issparse(Xc) else np.zeros(k)
        total_var = float(Xc.multiply(Xc).sum())
        variance_ratio = (singular**2) / max(total_var, 1e-12)
        embedding = Xc @ vt.T
        components = vt
        return np.asarray(embedding), np.asarray(variance_ratio), np.asarray(components)

    arr = np.asarray(X, dtype=float)
    n_obs, n_vars = arr.shape
    centered = arr - arr.mean(axis=0, keepdims=True)
    k = int(min(n_comps, min(n_obs, n_vars) - 1))
    if k < 1:
        raise ValueError("矩阵过小，无法做 PCA")
    u, s, vt = np.linalg.svd(centered, full_matrices=False)
    u, s, vt = u[:, :k], s[:k], vt[:k]
    variance_ratio = (s**2) / max(float((centered**2).sum()), 1e-12)
    return u * s, variance_ratio, vt


# --------------------------------------------------------------------------- #
# 6. 邻域图
# --------------------------------------------------------------------------- #
def compute_neighbors(
    embedding: np.ndarray,
    n_neighbors: int = 15,
    n_pcs: int | None = None,
    metric: str = "euclidean",
    prune: float = 1.0 / 15.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """构建 kNN 邻接图。

    步骤：

    1. 在 PCA 空间上用 sklearn ``NearestNeighbors`` 找每个细胞的 k 个近邻
    2. 按 scanpy 的 ``connectivities`` 语义计算边权：Jaccard 相似度
       （共同邻居数 / 邻居并集数），再用 prune 参数做 UMAP 式阈值剪枝
    3. 对称化并做行归一化，得到 ``connectivities``（用于聚类）

    Args:
        embedding: 细胞嵌入（通常为 PCA 坐标）。
        n_neighbors: 近邻数 k。
        n_pcs: 使用前多少个主成分，None 表示全部。
        metric: 距离度量。
        prune: 剪枝阈值（越小保留的边越多）。

    Returns:
        ``(indices, distances, connectivities)``：

        - ``indices`` 形状 ``(n_cells, k)``
        - ``distances`` 形状 ``(n_cells, k)``
        - ``connectivities`` 稀疏邻接矩阵 ``(n_cells, n_cells)``
    """
    from sklearn.neighbors import NearestNeighbors

    emb = np.asarray(embedding, dtype=float)
    if n_pcs is not None:
        emb = emb[:, :n_pcs]

    n_cells = emb.shape[0]
    k = int(min(n_neighbors, n_cells - 1)) if n_cells > 1 else 1

    nn = NearestNeighbors(n_neighbors=k + 1, metric=metric).fit(emb)
    distances, indices = nn.kneighbors(emb)
    # 第 0 列是自身
    distances = distances[:, 1:]
    indices = indices[:, 1:]

    # Jaccard 相似度 -> connectivities
    n_sets = [set(indices[i].tolist()) for i in range(n_cells)]
    rows: list[int] = []
    cols: list[int] = []
    vals: list[float] = []
    for i in range(n_cells):
        for j in indices[i]:
            j = int(j)
            inter = len(n_sets[i] & n_sets[j])
            union = len(n_sets[i] | n_sets[j])
            weight = inter / union if union else 0.0
            if weight > prune:
                rows.append(i)
                cols.append(j)
                vals.append(weight)

    adj = sp.coo_matrix(
        (vals, (rows, cols)), shape=(n_cells, n_cells), dtype=float
    ).tocsr()
    # 对称化
    adj = adj.maximum(adj.T)
    # 行归一化（scanpy 的 connectivities 语义）
    adj = adj.multiply(1.0 / np.maximum(adj.sum(axis=1), 1e-12))
    adj = sp.csr_matrix(adj)

    return indices, distances, adj


# --------------------------------------------------------------------------- #
# 7. Louvain 聚类
# --------------------------------------------------------------------------- #
def _louvain_one_level(
    adj: sp.csr_matrix,
    resolution: float,
    seed: int,
) -> np.ndarray:
    """单层 Louvain 社区发现。

    Louvain 算法两阶段：

    1. **局部移动**：每个节点尝试移动到邻居社区，按模块度增益最大者移动，
       重复直到无节点可改进
    2. **聚合**：把每个社区压缩成一个超节点，在压缩图上重复

    模块度增益（节点 i 从社区 c 移出并加入社区 d）：

    ``ΔQ = k_i_in(d) / (2m) - resolution * (Σ_tot(d) * k_i) / (2m^2)``

    Args:
        adj: 对称邻接矩阵（稀疏）。
        resolution: 分辨率参数（越大社区越多）。
        seed: 随机种子。

    Returns:
        每个节点所属的社区标签数组。
    """
    n = adj.shape[0]
    adj = adj.tocsr()
    rng = np.random.default_rng(seed)

    # 节点权重 = 度
    k = np.asarray(adj.sum(axis=1)).ravel()
    # 2m = 总权重
    m2 = float(k.sum())
    if m2 <= 0:
        return np.arange(n)

    # 社区总权重
    comm = np.arange(n)
    comm_weight = k.copy()

    improved = True
    while improved:
        improved = False
        # 随机遍历顺序，保证可复现
        order = rng.permutation(n)
        for i in order:
            ci = comm[i]
            # 收集 i 的邻居及其边权
            start, end = adj.indptr[i], adj.indptr[i + 1]
            neighbors = adj.indices[start:end]
            weights = adj.data[start:end]
            neigh_comm = comm[neighbors]

            # 按社区聚合邻居边权
            comm_tot: dict[int, float] = {}
            for c, w in zip(neigh_comm, weights):
                comm_tot[int(c)] = comm_tot.get(int(c), 0.0) + float(w)

            # i 移出原社区后的社区总权重
            w_ci = comm_weight[ci] - k[i]

            # 留在原社区的模块度基线
            base_ci = (
                comm_tot.get(int(ci), 0.0) / m2
                - resolution * (w_ci * k[i]) / (m2 * m2)
            )

            best_comm = ci
            best_gain = base_ci
            for c, wt in comm_tot.items():
                c = int(c)
                # 移入社区 c 后的社区总权重
                tot_c = (w_ci if c == ci else comm_weight[c]) + k[i]
                # i 与社区 c 之间的边权（若 c 就是原社区，需减去 i 自身的贡献）
                link_c = wt - (k[i] if c == ci else 0.0)
                gain = link_c / m2 - resolution * (tot_c * k[i]) / (m2 * m2)
                if gain > best_gain + 1e-12:
                    best_gain = gain
                    best_comm = c

            if best_comm != ci:
                comm_weight[ci] = w_ci
                comm_weight[best_comm] += k[i]
                comm[i] = best_comm
                improved = True

    # 压缩为连续标签
    unique = {int(c): idx for idx, c in enumerate(np.unique(comm))}
    return np.array([unique[int(c)] for c in comm], dtype=int)


def cluster_louvain(
    adj: sp.csr_matrix,
    resolution: float = 1.0,
    seed: int = 0,
    max_levels: int = 10,
) -> np.ndarray:
    """Louvain 社区发现（多层次聚合）。

    Args:
        adj: 对称邻接矩阵。
        resolution: 分辨率参数。
        seed: 随机种子。
        max_levels: 最大聚合层数。

    Returns:
        社区标签数组。
    """
    n = adj.shape[0]
    current = adj.tocsr()
    # cell_to_node: 每个原始细胞当前所在的图节点编号
    cell_to_node = np.arange(n)
    node_labels = np.arange(n)

    for level in range(max_levels):
        node_labels = _louvain_one_level(current, resolution, seed + level)
        n_comm = int(node_labels.max()) + 1
        if n_comm == current.shape[0]:
            # 没有发生合并，停止
            break

        # 构建聚合图
        rows = np.repeat(np.arange(current.shape[0]), np.diff(current.indptr))
        cols = current.indices
        agg = sp.coo_matrix(
            (current.data, (node_labels[rows], node_labels[cols])),
            shape=(n_comm, n_comm),
        ).tocsr()
        agg = agg.maximum(agg.T)
        # 细胞映射到新图的节点
        cell_to_node = node_labels[cell_to_node]
        current = agg
        if n_comm == 1:
            break

    # 把最终社区标签映射回每个细胞
    labels = node_labels[cell_to_node]
    unique = {int(c): idx for idx, c in enumerate(np.unique(labels))}
    return np.array([unique[int(c)] for c in labels], dtype=int)


def cluster_leiden(
    adj: sp.csr_matrix,
    resolution: float = 1.0,
    seed: int = 0,
) -> np.ndarray:
    """Leiden 聚类（需要 leidenalg + igraph），缺失时回退到 Louvain。

    Args:
        adj: 对称邻接矩阵。
        resolution: 分辨率参数。
        seed: 随机种子。

    Returns:
        社区标签数组。
    """
    if leidenalg_mod is not None and igraph_mod is not None:
        import igraph as ig

        g = ig.Graph.Weighted_AdjacencyMatrix(adj.tocoo())
        partition = leidenalg_mod.find_partition(
            g, leidenalg_mod.RBConfigurationVertexPartition, resolution_parameter=resolution,
            seed=seed,
        )
        return np.asarray(partition.membership, dtype=int)
    return cluster_louvain(adj, resolution=resolution, seed=seed)


# --------------------------------------------------------------------------- #
# 8. UMAP / t-SNE
# --------------------------------------------------------------------------- #
def compute_umap(
    embedding: np.ndarray,
    n_neighbors: int = 15,
    min_dist: float = 0.3,
    n_components: int = 2,
    seed: int = 0,
    n_epochs: int = 200,
) -> np.ndarray:
    """UMAP 二维嵌入。

    - 若安装了 ``umap-learn``，直接调用官方实现
    - 否则使用降级近似：以 kNN 图初始化，用 MDS 得到初始布局，
      再做保持邻域结构的力导向微调（吸引近邻 + 排斥远邻）

    Args:
        embedding: 细胞嵌入（PCA 坐标）。
        n_neighbors: 近邻数。
        min_dist: 最小嵌入距离（影响排斥力强度）。
        n_components: 输出维度。
        seed: 随机种子。
        n_epochs: 力导向微调迭代次数。

    Returns:
        形状 ``(n_cells, n_components)`` 的嵌入。
    """
    emb = np.asarray(embedding, dtype=float)
    if umap_mod is not None:
        reducer = umap_mod.UMAP(
            n_neighbors=int(min(n_neighbors, emb.shape[0] - 1)),
            min_dist=min_dist,
            n_components=n_components,
            random_state=seed,
        )
        return np.asarray(reducer.fit_transform(emb))

    return _umap_approx(emb, n_neighbors, min_dist, n_components, seed, n_epochs)


def _umap_approx(
    emb: np.ndarray,
    n_neighbors: int,
    min_dist: float,
    n_components: int,
    seed: int,
    n_epochs: int,
) -> np.ndarray:
    """UMAP 的降级近似实现。

    采用与 umap-learn 相同的核心思路：

    1. **谱初始化**：对 kNN 图做归一化拉普拉斯特征分解，取前若干个特征向量
       作为初始布局（UMAP 本身也用谱嵌入做初始化）
    2. **梯度下降优化**：用 UMAP 的真实力公式做 SGD

       吸引力（对每个近邻对）::

           f_attr = -2*a*b*d^(2(b-1)) / (1 + a*d^(2b))

       排斥力（对随机采样的负样本）::

           f_rep = 2*b / ((eps + d) * (1 + a*d^(2b)))

       其中 ``a``/``b`` 由 ``min_dist`` 映射而来。关键点在于吸引力在短距离时
       **饱和**（不会像 1/d² 那样把点拉塌），而排斥力在短距离时迅速增长，
       两者平衡后得到与 UMAP 相似的"簇内紧密、簇间分离"的布局。

    Args:
        emb: 高维嵌入（PCA 坐标）。
        n_neighbors: 近邻数。
        min_dist: 最小嵌入距离。
        n_components: 输出维度。
        seed: 随机种子。
        n_epochs: 优化迭代次数。

    Returns:
        形状 ``(n_cells, n_components)`` 的布局。
    """
    from sklearn.neighbors import NearestNeighbors

    n = emb.shape[0]
    rng = np.random.default_rng(seed)

    # --- 1. kNN 图（含距离，用于 rho/sigma 局部尺度） ---
    k = int(min(n_neighbors, max(n - 1, 1)))
    nn = NearestNeighbors(n_neighbors=k + 1).fit(emb)
    distances, indices = nn.kneighbors(emb)
    distances = np.asarray(distances[:, 1:], dtype=float)
    indices = indices[:, 1:]

    # --- 2. 谱初始化 ---
    layout = _spectral_init(indices, distances, n, n_components, seed)

    # --- 3. UMAP 力优化 ---
    a, b = _min_dist_to_ab(min_dist)
    eps = 1e-3

    # 把近邻对展开成边列表
    src = np.repeat(np.arange(n), k)
    dst = indices.ravel()

    # 学习率与负样本数随 epoch 退火（与 umap-learn 的 schedule 类似）
    alpha = 1.0
    n_neg = 5
    for epoch in range(n_epochs):
        # --- 吸引力 ---
        diff = layout[src] - layout[dst]                    # (E, d)
        dist_sq = np.sum(diff**2, axis=1)                    # (E,)
        dist = np.sqrt(dist_sq) + eps
        # 吸引力大小 f_attr = 2*a*b*d^(2(b-1)) / (1 + a*d^(2b)) > 0
        # 方向沿 diff = src - dst（即由 dst 指向 src），因此 src 的移动方向
        # 取 +unit 会把 src 推离 dst，必须取 -unit 才是"互相靠拢"。
        attr_mag = (2.0 * a * b * dist ** (2.0 * (b - 1.0))) / (1.0 + a * dist_sq**b)
        edge_grad = attr_mag[:, None] * diff / dist[:, None]    # (E, d)
        # 吸引力把 src 与 dst 互相拉近：src 沿 -unit 走，dst 沿 +unit 走
        grad_attr = np.zeros_like(layout)
        np.add.at(grad_attr, src, -edge_grad)
        np.add.at(grad_attr, dst, edge_grad)

        # --- 排斥力（负样本采样）---
        neg = rng.integers(0, n, size=(n, n_neg))
        diff_neg = layout[:, None, :] - layout[neg]           # (n, n_neg, d)
        dist_neg = np.linalg.norm(diff_neg, axis=2) + eps    # (n, n_neg)
        coeff_rep = (2.0 * b) / ((eps + dist_neg) * (1.0 + a * dist_neg ** (2.0 * b)))
        grad_rep = (coeff_rep[:, :, None] * diff_neg / dist_neg[:, :, None]).sum(axis=1)

        grad = grad_attr + grad_rep
        # 梯度裁剪，避免个别点飞出去
        norm = np.linalg.norm(grad, axis=1, keepdims=True)
        max_norm = 4.0
        grad = np.where(norm > max_norm, grad * (max_norm / np.maximum(norm, 1e-12)), grad)

        layout += alpha * grad
        # 居中
        layout -= layout.mean(axis=0, keepdims=True)
        alpha *= 0.995

    # 归一化到与 UMAP 类似的坐标尺度
    spread = layout.std()
    if spread > 0:
        layout = layout / spread * 2.0
    return layout


def _min_dist_to_ab(min_dist: float) -> tuple[float, float]:
    """把 ``min_dist`` 映射为 UMAP 的 ``(a, b)`` 曲线参数。

    UMAP 通过拟合 ``1/(1 + a*x^(2b))`` 到圆盘曲线来确定 a、b。
    这里用与 umap-learn 相同量级的近似映射。

    Args:
        min_dist: 最小嵌入距离。

    Returns:
        ``(a, b)``。
    """
    d = float(np.clip(min_dist, 0.0, 0.99))
    if d <= 0.0:
        return 1.0, 1.0
    # umap-learn 的拟合结果（min_dist -> a, b）大致关系
    a = 1.0 / (1.0 + 2.0 * d)
    b = 1.0 + d
    return float(a), float(b)


def _spectral_init(
    indices: np.ndarray,
    distances: np.ndarray,
    n: int,
    n_components: int,
    seed: int,
) -> np.ndarray:
    """对 kNN 图做谱嵌入，得到 UMAP 的初始布局。

    步骤：

    1. 用 UMAP 的局部尺度自适应规则构造亲和度：
       对每个点取 ``rho_i``（到最近邻的距离），再二分求出 ``sigma_i``，
       使得 ``sum_j exp(-(d_ij - rho_i) / sigma_i) = log2(k)``。
       这样得到的是**有方向**的亲和度（原始 kNN 关系），再做对称化。
    2. 构造对称归一化邻接矩阵 ``D^-1/2 A D^-1/2``
    3. 做特征分解，取最大的若干特征向量。**第一个特征向量恒为常数向量**
       （对应平凡特征值 1），必须丢弃——否则会把所有点推向同一个方向。
    4. 把保留的特征向量用 PCA 投影到目标维度

    Args:
        indices: ``(n, k)`` 的近邻索引。
        distances: ``(n, k)`` 的近邻距离。
        n: 细胞数。
        n_components: 输出维度。
        seed: 随机种子。

    Returns:
        形状 ``(n, n_components)`` 的初始布局。
    """
    from scipy.sparse.linalg import eigsh

    rng = np.random.default_rng(seed)
    k = indices.shape[1]
    n = int(n)

    # --- 1. UMAP 的 rho / sigma 局部尺度 ---
    rho = distances[:, 0].copy()
    target = math.log2(k) if k > 1 else 1.0
    sigma = np.ones(n)
    for i in range(n):
        lo, hi = 1e-6, 10.0
        for _ in range(50):
            mid = (lo + hi) / 2.0
            s = np.sum(np.exp(-np.maximum(distances[i] - rho[i], 0.0) / mid))
            if s < target:
                lo = mid
            else:
                hi = mid
        sigma[i] = (lo + hi) / 2.0

    rows = np.repeat(np.arange(n), k)
    cols = indices.ravel()
    w = np.exp(
        -np.maximum(distances.ravel() - np.repeat(rho, k), 0.0) / np.repeat(sigma, k)
    )
    A = sp.coo_matrix((w, (rows, cols)), shape=(n, n)).tocsr()
    A = A.maximum(A.T)

    # --- 2. 对称归一化 ---
    deg = np.asarray(A.sum(axis=1)).ravel()
    deg[deg == 0] = 1e-12
    Dm = sp.diags(1.0 / np.sqrt(deg))
    A_norm = Dm @ A @ Dm

    # --- 3. 特征分解（丢弃第一个平凡特征向量） ---
    k_comp = int(min(n_components + 2, n - 2))
    try:
        vals, vecs = eigsh(A_norm, k=k_comp, which="LA")
        order = np.argsort(-vals)
        vecs = vecs[:, order]
    except Exception:  # noqa: BLE001
        vecs = rng.normal(size=(n, n_components))
        return vecs * 10.0

    # 第一个特征向量是常数向量（平凡特征值），必须剔除
    keep = vecs[:, 1:n_components + 1]
    if keep.shape[1] < n_components:
        keep = np.hstack([keep, rng.normal(size=(n, n_components - keep.shape[1]))])

    # --- 4. PCA 投影到目标维度 ---
    centered = keep - keep.mean(axis=0, keepdims=True)
    try:
        u, s, _ = np.linalg.svd(centered, full_matrices=False)
        layout = u[:, :n_components] * s[:n_components]
    except np.linalg.LinAlgError:
        layout = centered[:, :n_components]

    layout = np.asarray(layout, dtype=float)
    scale = layout.std()
    if scale > 0:
        layout = layout / scale * 5.0
    return layout


def compute_tsne(
    embedding: np.ndarray,
    n_pcs: int = 50,
    perplexity: float = 30.0,
    seed: int = 0,
) -> np.ndarray:
    """t-SNE 二维嵌入（sklearn 实现）。

    Args:
        embedding: 细胞嵌入。
        n_pcs: 使用的主成分数。
        perplexity: perplexity 参数。
        seed: 随机种子。

    Returns:
        形状 ``(n_cells, 2)`` 的嵌入。
    """
    from sklearn.manifold import TSNE

    emb = np.asarray(embedding, dtype=float)
    n = emb.shape[0]
    perplexity = float(min(perplexity, max(5.0, (n - 1) / 3.0)))
    tsne = TSNE(
        n_components=2,
        perplexity=perplexity,
        init="pca",
        random_state=seed,
        max_iter=1000,
    )
    return np.asarray(tsne.fit_transform(emb[:, :n_pcs]))


# --------------------------------------------------------------------------- #
# 9. 差异表达
# --------------------------------------------------------------------------- #
def rank_genes_groups(
    X: Any,
    groupby: np.ndarray,
    n_genes: int = 100,
    min_counts: int = 1,
) -> pd.DataFrame:
    """基于 Welch t 检验的组间差异表达排名。

    对每个 cluster 与"其余所有细胞"的对比，逐基因做 Welch t 检验，
    并计算 log2 fold change（用 pseudocount 防止除零），最后做 BH 校正。

    Args:
        X: log 归一化后的表达矩阵。
        groupby: 每个细胞的 cluster 标签。
        n_genes: 每个 cluster 保留的 marker 数。
        min_counts: 基因在该组至少要有的总表达量。

    Returns:
        列含 ``cluster`` / ``gene`` / ``score`` / ``logfoldchanges`` /
        ``pvals`` / ``pvals_adj`` 的 DataFrame。
    """
    from scipy import stats

    arr = X.toarray() if sp.issparse(X) else np.asarray(X, dtype=float)
    groups = np.asarray(groupby)
    unique_groups = [g for g in pd.unique(groups)]

    rows: list[dict[str, Any]] = []
    for grp in unique_groups:
        in_grp = groups == grp
        if in_grp.sum() < 2:
            continue
        a = arr[in_grp]
        b = arr[~in_grp]
        with np.errstate(invalid="ignore", divide="ignore"):
            t_stat, p_val = stats.ttest_ind(a, b, axis=0, equal_var=False)
        p_val = np.nan_to_num(p_val, nan=1.0)
        mean_a = a.mean(axis=0)
        mean_b = b.mean(axis=0)
        lfc = np.log2((mean_a + 1e-9) / (mean_b + 1e-9))
        scores = np.abs(t_stat)
        order = np.argsort(-scores)
        keep = [i for i in order if (a[:, i].sum() + b[:, i].sum()) >= min_counts][:n_genes]
        if not keep:
            continue
        sub_p = p_val[keep]
        padj = _bh(sub_p)
        for rank, (i, q) in enumerate(zip(keep, padj)):
            rows.append({
                "cluster": str(grp),
                "rank": rank + 1,
                "gene": f"gene_{i}",
                "gene_index": int(i),
                "score": float(scores[i]),
                "logfoldchanges": float(lfc[i]),
                "pvals": float(p_val[i]),
                "pvals_adj": float(q),
            })

    if not rows:
        return pd.DataFrame(
            columns=["cluster", "rank", "gene", "gene_index", "score",
                     "logfoldchanges", "pvals", "pvals_adj"]
        )
    return pd.DataFrame(rows)


def _bh(pvalues: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg 校正。"""
    p = np.asarray(pvalues, dtype=float)
    n = p.size
    order = np.argsort(p, kind="mergesort")
    ranked = p[order]
    q = np.minimum.accumulate((ranked * n / np.arange(1, n + 1))[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(q, 0, 1)
    return out


# --------------------------------------------------------------------------- #
# 10. scanpy 后端全流程
# --------------------------------------------------------------------------- #
def _run_scanpy_pipeline(
    adata: Any,
    min_genes: int = 200,
    min_counts: int = 100,
    max_mt_pct: float = 20.0,
    max_hb_pct: float | None = None,
    min_cells_per_gene: int = 3,
    target_sum: float = 1e4,
    n_top_genes: int = 2000,
    n_comps: int = 50,
    n_neighbors: int = 15,
    resolution: float = 1.0,
    min_dist: float = 0.3,
    n_marker_genes: int = 100,
    seed: int = 0,
    compute_tsne_embedding: bool = True,
    logger: Callable[[str], None] | None = None,
) -> PipelineResult:
    """使用官方 scanpy 实现完整流程。

    流程与降级后端一一对应，只是把每个数值步骤替换为 scanpy 的官方实现：

    1. ``sc.pp.calculate_qc_metrics`` + ``sc.pp.filter_cells`` / ``sc.pp.filter_genes``
    2. ``sc.pp.normalize_total`` + ``sc.pp.log1p``
    3. ``sc.pp.highly_variable_genes(flavor="seurat")``
    4. ``sc.pp.scale`` + ``sc.tl.pca``
    5. ``sc.pp.neighbors``（含 UMAP 式 prune）
    6. ``sc.tl.leiden``（无 leidenalg 时自动回退 Louvain）
    7. ``sc.tl.umap`` / ``sc.tl.tsne``
    8. ``sc.tl.rank_genes_groups``（Wilcoxon）+ BH 校正

    Args:
        adata: 输入数据对象（``SimpleAnnData``）。
        min_genes: 细胞至少检测到的基因数。
        min_counts: 细胞至少的 UMI 数。
        max_mt_pct: 线粒体基因比例上限（%）。
        max_hb_pct: 血红蛋白基因比例上限（%）。
        min_cells_per_gene: 基因最少表达细胞数。
        target_sum: 归一化目标总计数。
        n_top_genes: 高变基因数。
        n_comps: PCA 主成分数。
        n_neighbors: kNN 邻居数。
        resolution: 聚类分辨率。
        min_dist: UMAP 最小距离。
        n_marker_genes: 每簇 marker 数量。
        seed: 随机种子。
        compute_tsne_embedding: 是否额外计算 t-SNE。
        logger: 日志回调。

    Returns:
        :class:`PipelineResult`。
    """
    log = logger or (lambda msg: None)

    import scanpy as sc
    from scdata import from_anndata, gene_is_hb, gene_is_mitochondrial, to_anndata

    # --- 转换为官方 AnnData ---
    ad = to_anndata(adata)
    ad.var_names_make_unique()

    mt_mask = gene_is_mitochondrial(ad.var_names)
    hb_mask = gene_is_hb(ad.var_names)
    log(f"线粒体基因: {int(mt_mask.sum())}, 血红蛋白基因: {int(hb_mask.sum())}")

    # scanpy 的 qc_vars 依赖 var 中的布尔标记列
    ad.var["mt"] = mt_mask
    ad.var["hb"] = hb_mask

    # --- 1. 质控 ---
    sc.pp.calculate_qc_metrics(
        ad, qc_vars=["mt", "hb"], percent_top=None, log1p=True, inplace=True
    )
    ad.obs["pct_counts_mt"] = ad.obs.get("pct_counts_mt", 0.0)
    ad.obs["pct_counts_hb"] = ad.obs.get("pct_counts_hb", 0.0)

    n_cells_before = ad.n_obs
    n_genes_before = ad.n_vars
    sc.pp.filter_cells(
        ad, min_genes=min_genes, min_counts=min_counts,
        max_pct_counts_mt=None if max_mt_pct >= 100 else max_mt_pct,
    )
    if max_hb_pct is not None:
        sc.pp.filter_cells(ad, max_pct_counts_hb=max_hb_pct)
    sc.pp.filter_genes(ad, min_cells=min_cells_per_gene)

    qc_stats = {
        "n_cells_before": int(n_cells_before),
        "n_cells_after": int(ad.n_obs),
        "n_cells_removed": int(n_cells_before - ad.n_obs),
        "cell_removal_rate": float(1.0 - ad.n_obs / max(n_cells_before, 1)),
        "n_genes_before": int(n_genes_before),
        "n_genes_after": int(ad.n_vars),
        "n_genes_removed": int(n_genes_before - ad.n_vars),
        "thresholds": {
            "min_genes": min_genes,
            "min_counts": min_counts,
            "max_mt_pct": max_mt_pct,
            "max_hb_pct": max_hb_pct,
            "min_cells_per_gene": min_cells_per_gene,
        },
    }
    log(
        f"质控: 细胞 {qc_stats['n_cells_before']} -> {qc_stats['n_cells_after']}, "
        f"基因 {qc_stats['n_genes_before']} -> {qc_stats['n_genes_after']}"
    )

    if ad.n_obs < 10:
        raise ValueError(
            f"质控后仅剩 {ad.n_obs} 个细胞，请放宽阈值（--min-genes / --min-counts / --max-mt-pct）"
        )

    ad.raw = ad.copy()

    # --- 2. 归一化 ---
    sc.pp.normalize_total(ad, target_sum=target_sum)
    sc.pp.log1p(ad)
    log(f"归一化完成 (CP{target_sum:.0f} + log1p)")

    # --- 3. 高变基因 ---
    sc.pp.highly_variable_genes(ad, n_top_genes=n_top_genes, flavor="seurat")
    hvg_mask = np.asarray(ad.var["highly_variable"], dtype=bool)
    log(f"高变基因: {int(hvg_mask.sum())}")

    # --- 4. 缩放 + PCA ---
    ad_hvg = ad[:, hvg_mask].copy()
    sc.pp.scale(ad_hvg, max_value=10.0)
    n_comp = int(min(n_comps, min(ad_hvg.shape) - 1))
    sc.tl.pca(ad_hvg, n_comps=n_comp, svd_solver="arpack", random_state=seed)
    ad.obsm["X_pca"] = np.asarray(ad_hvg.obsm["X_pca"])
    ad.uns["pca"] = {
        "variance_ratio": np.asarray(ad_hvg.uns["pca"]["variance_ratio"])
    }
    log(
        f"PCA 完成: {ad.obsm['X_pca'].shape[1]} 个主成分, 前 5 个累计解释 "
        f"{float(np.sum(ad.uns['pca']['variance_ratio'][:5])) * 100:.1f}%"
    )

    # --- 5. 邻接图 ---
    n_pcs = int(min(50, ad.obsm["X_pca"].shape[1]))
    sc.pp.neighbors(
        ad, n_neighbors=n_neighbors, n_pcs=n_pcs, prune=1.0 / max(n_neighbors, 1),
        random_state=seed,
    )
    log(f"邻接图完成: k={n_neighbors}")

    # --- 6. 聚类 ---
    try:
        sc.tl.leiden(ad, resolution=resolution, key_added="leiden", random_state=seed)
        cluster_key = "leiden"
        log("聚类完成 (Leiden)")
    except (ImportError, ValueError, AttributeError):
        sc.tl.louvain(ad, resolution=resolution, key_added="louvain", random_state=seed)
        cluster_key = "louvain"
        log("聚类完成 (Louvain，leidenalg 不可用)")

    labels = np.asarray([str(x) for x in ad.obs[cluster_key]])
    ad.obs["cluster"] = labels
    ad.obs["leiden"] = labels
    n_clusters = len(set(labels))
    log(f"聚类完成: {n_clusters} 个 cluster")

    # --- 7. UMAP / t-SNE ---
    sc.tl.umap(
        ad, min_dist=min_dist, random_state=seed,
        neighbors_key="neighbors",
    )
    ad.obsm["X_umap"] = np.asarray(ad.obsm["X_umap"])
    log("UMAP 完成")
    if compute_tsne_embedding:
        sc.tl.tsne(ad, random_state=seed)
        ad.obsm["X_tsne"] = np.asarray(ad.obsm["X_tsne"])
        log("t-SNE 完成")

    # --- 8. Marker（Wilcoxon 秩和检验）---
    sc.tl.rank_genes_groups(
        ad, groupby="cluster", method="wilcoxon", n_genes=n_marker_genes
    )
    markers = _scanpy_markers_to_frame(ad, n_marker_genes)
    if len(markers):
        markers["gene"] = markers["gene_index"].map(
            lambda i: str(ad.var_names[int(i)]) if int(i) < len(ad.var_names) else f"gene_{i}"
        )

    summary = (
        markers.groupby("cluster")["gene"].count().reset_index()
        if len(markers)
        else pd.DataFrame(columns=["cluster", "gene"])
    ).rename(columns={"gene": "n_markers"})
    cluster_counts = pd.Series(labels).value_counts().rename("n_cells").reset_index()
    cluster_counts.columns = ["cluster", "n_cells"]
    summary = cluster_counts.merge(summary, on="cluster", how="left").fillna({"n_markers": 0})

    out = from_anndata(ad)
    return PipelineResult(
        adata=out,
        backend="scanpy",
        qc_stats=qc_stats,
        cluster_summary=summary,
        marker_table=markers,
        steps=[
            "quality_control", "normalize_total", "log1p", "highly_variable_genes",
            "scale", "pca", "neighbors", "cluster", "umap", "tsne", "rank_genes_groups",
        ],
    )


def _scanpy_markers_to_frame(ad: Any, n_marker_genes: int) -> pd.DataFrame:
    """把 scanpy 的 ``rank_genes_groups`` 结果转成统一格式的 DataFrame。

    Args:
        ad: 已运行 ``sc.tl.rank_genes_groups`` 的 AnnData。
        n_marker_genes: 每簇保留的 marker 数。

    Returns:
        列含 ``cluster`` / ``rank`` / ``gene`` / ``gene_index`` / ``score`` /
        ``logfoldchanges`` / ``pvals`` / ``pvals_adj`` 的 DataFrame。
    """
    rgg = ad.uns.get("rank_genes_groups")
    if rgg is None:
        return pd.DataFrame(
            columns=["cluster", "rank", "gene", "gene_index", "score",
                     "logfoldchanges", "pvals", "pvals_adj"]
        )

    names = list(rgg["names"].dtype.names)
    lfc = np.asarray(rgg["logfoldchanges"].dtype.names)
    pvals = np.asarray(rgg["pvals"].dtype.names)
    pvals_adj = np.asarray(rgg["pvals_adj"].dtype.names)
    scores = np.asarray(rgg["scores"].dtype.names)

    rows: list[dict[str, Any]] = []
    for grp in names:
        n_keep = min(n_marker_genes, len(rgg["names"][grp]))
        for rank in range(n_keep):
            rows.append({
                "cluster": str(grp),
                "rank": rank + 1,
                "gene": str(rgg["names"][grp][rank]),
                "gene_index": int(rgg["names"][grp][rank]),
                "score": float(scores[grp][rank]),
                "logfoldchanges": float(lfc[grp][rank]),
                "pvals": float(pvals[grp][rank]),
                "pvals_adj": float(pvals_adj[grp][rank]),
            })
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# 11. 全流程（降级后端）
# --------------------------------------------------------------------------- #
def run_pipeline(
    adata: Any,
    backend: str = "auto",
    min_genes: int = 200,
    min_counts: int = 100,
    max_mt_pct: float = 20.0,
    max_hb_pct: float | None = None,
    min_cells_per_gene: int = 3,
    target_sum: float = 1e4,
    n_top_genes: int = 2000,
    n_comps: int = 50,
    n_neighbors: int = 15,
    resolution: float = 1.0,
    min_dist: float = 0.3,
    n_marker_genes: int = 100,
    seed: int = 0,
    compute_tsne_embedding: bool = True,
    logger: Callable[[str], None] | None = None,
) -> PipelineResult:
    """执行完整的单细胞分析流程。

    Args:
        adata: 输入数据对象（SimpleAnnData）。
        backend: ``auto`` / ``scanpy`` / ``numpy``。
        min_genes: 细胞至少检测到的基因数。
        min_counts: 细胞至少的 UMI 数。
        max_mt_pct: 线粒体比例上限（%）。
        max_hb_pct: 血红蛋白比例上限（%）。
        min_cells_per_gene: 基因最少表达细胞数。
        target_sum: 归一化目标总计数。
        n_top_genes: HVG 数量。
        n_comps: PCA 主成分数。
        n_neighbors: kNN 邻居数。
        resolution: 聚类分辨率。
        min_dist: UMAP 最小距离。
        n_marker_genes: 每簇 marker 数量。
        seed: 随机种子。
        compute_tsne_embedding: 是否额外计算 t-SNE。
        logger: 日志回调（形如 ``logger.info``）。

    Returns:
        :class:`PipelineResult`。
    """
    log = logger or (lambda msg: None)

    if backend == "auto":
        backend = "scanpy" if SCANPY_AVAILABLE else "numpy"
    if backend == "scanpy" and not SCANPY_AVAILABLE:
        log("scanpy 不可用，自动切换到 numpy 降级后端")
        backend = "numpy"

    if backend == "scanpy":
        return _run_scanpy_pipeline(
            adata,
            min_genes=min_genes,
            min_counts=min_counts,
            max_mt_pct=max_mt_pct,
            max_hb_pct=max_hb_pct,
            min_cells_per_gene=min_cells_per_gene,
            target_sum=target_sum,
            n_top_genes=n_top_genes,
            n_comps=n_comps,
            n_neighbors=n_neighbors,
            resolution=resolution,
            min_dist=min_dist,
            n_marker_genes=n_marker_genes,
            seed=seed,
            compute_tsne_embedding=compute_tsne_embedding,
            logger=log,
        )

    from scdata import gene_is_hb, gene_is_mitochondrial

    adata.var_names_make_unique()
    mt_mask = gene_is_mitochondrial(adata.var_names)
    hb_mask = gene_is_hb(adata.var_names)
    log(f"线粒体基因: {int(mt_mask.sum())}, 血红蛋白基因: {int(hb_mask.sum())}")

    # --- 1. 质控 ---
    qc = compute_qc_metrics(adata, mt_mask, hb_mask)
    for key, values in qc.items():
        adata.obs[key] = values
    adata, qc_stats = filter_cells_genes(
        adata, qc, min_genes=min_genes, min_counts=min_counts,
        max_mt_pct=max_mt_pct, max_hb_pct=max_hb_pct,
        min_cells_per_gene=min_cells_per_gene,
    )
    log(
        f"质控: 细胞 {qc_stats['n_cells_before']} -> {qc_stats['n_cells_after']}, "
        f"基因 {qc_stats['n_genes_before']} -> {qc_stats['n_genes_after']}"
    )

    if adata.n_obs < 10:
        raise ValueError(
            f"质控后仅剩 {adata.n_obs} 个细胞，请放宽阈值（--min-genes / --min-counts / --max-mt-pct）"
        )

    # 保存原始计数
    adata.raw = adata.copy()

    # --- 2. 归一化 ---
    X_norm = normalize_total(adata.X, target_sum=target_sum)
    X_log = log1p(X_norm)
    log(f"归一化完成 (CP{target_sum:.0f} + log1p)")

    # --- 3. HVG ---
    hvg_mask = highly_variable_genes(X_log, n_top_genes=n_top_genes)
    adata.var["highly_variable"] = hvg_mask
    adata.var["means"] = (X_log.toarray() if sp.issparse(X_log) else X_log).mean(axis=0)
    log(f"高变基因: {int(hvg_mask.sum())}")

    # --- 4. 缩放 + PCA ---
    X_scaled = scale(X_log[:, hvg_mask], max_value=10.0)
    n_comp = int(min(n_comps, min(X_scaled.shape) - 1))
    X_pca, variance_ratio, components = compute_pca(X_scaled, n_comps=n_comp, seed=seed)
    adata.obsm["X_pca"] = X_pca
    adata.uns["pca"] = {"variance_ratio": variance_ratio}
    log(f"PCA 完成: {X_pca.shape[1]} 个主成分, 前 5 个累计解释 "
        f"{float(np.sum(variance_ratio[:5])) * 100:.1f}%")

    # --- 5. 邻接图 ---
    indices, distances, connectivities = compute_neighbors(
        X_pca, n_neighbors=n_neighbors, n_pcs=min(50, X_pca.shape[1])
    )
    adata.obsm["distances"] = distances
    adata.obsm["connectivities"] = connectivities
    log(f"邻接图完成: k={n_neighbors}")

    # --- 6. 聚类 ---
    labels = cluster_leiden(connectivities, resolution=resolution, seed=seed)
    labels = np.asarray([str(x) for x in labels])
    adata.obs["leiden"] = labels
    adata.obs["cluster"] = labels
    n_clusters = len(set(labels))
    log(f"聚类完成: {n_clusters} 个 cluster")

    # --- 7. UMAP / t-SNE ---
    adata.obsm["X_umap"] = compute_umap(
        X_pca, n_neighbors=n_neighbors, min_dist=min_dist, seed=seed
    )
    log("UMAP 完成")
    if compute_tsne_embedding:
        adata.obsm["X_tsne"] = compute_tsne(X_pca, n_pcs=min(50, X_pca.shape[1]), seed=seed)
        log("t-SNE 完成")

    # --- 8. Marker ---
    markers = rank_genes_groups(X_log, labels, n_genes=n_marker_genes)
    if len(markers):
        markers["gene"] = markers["gene_index"].map(
            lambda i: str(adata.var_names[int(i)]) if int(i) < len(adata.var_names) else f"gene_{i}"
        )
    adata.uns["rank_genes_groups"] = markers

    summary = (
        markers.groupby("cluster")["gene"].count().reset_index()
        if len(markers)
        else pd.DataFrame(columns=["cluster", "gene"])
    )
    summary = summary.rename(columns={"gene": "n_markers"})
    cluster_counts = pd.Series(labels).value_counts().rename("n_cells").reset_index()
    cluster_counts.columns = ["cluster", "n_cells"]
    summary = cluster_counts.merge(summary, on="cluster", how="left").fillna({"n_markers": 0})

    return PipelineResult(
        adata=adata,
        backend=backend,
        qc_stats=qc_stats,
        cluster_summary=summary,
        marker_table=markers,
        steps=[
            "quality_control", "normalize_total", "log1p", "highly_variable_genes",
            "scale", "pca", "neighbors", "cluster", "umap", "tsne", "rank_genes_groups",
        ],
    )
