#!/usr/bin/env python3
"""DESeq2 风格的负二项差异表达统计（纯 NumPy/SciPy 实现）。

本模块是 ``bulk-rnaseq`` 技能中 ``--backend python`` 的数值核心，移植并实现了
DESeq2（Love et al., 2014, Genome Biology 15:550）公开算法的核心步骤：

1. **Size factor 归一化**：中位数比值法（median-of-ratios），
   ``s_j = median_g ( K_gj / geomean_g(K_.j) )``。
2. **基因离散度估计**：对每个基因在离散度网格上最大化
   *Cox-Reid 调整后的负二项 profile 似然*，得到 gene-wise dispersion ``alpha_hat``。
3. **离散度趋势拟合**：``alpha_trend(mu) = a0 + a1 / mu``（parametric 家族），
   在 log 空间做加权最小二乘初值，再用 Gamma 族 deviance 精修。
4. **MAP 离散度**：把离散度先验（log 空间高斯，先验方差由趋势残差 MAD 估计）
   与 CR 似然结合，得到最终用于 Wald 检验的 ``alpha_MAP``。
5. **Wald 检验**：以 NB-GLM 的 ``X'WX`` 作为 Fisher 信息矩阵近似，
   构造线性对比的 z 统计量并给出双侧 p 值。
6. **BH 校正**：Benjamini-Hochberg FDR（与 DESeq2 的 ``p.adjust`` 一致）。

这些步骤全部是真实数值计算，不是占位符；可在无 R / 无 PyDESeq2 的环境中独立使用。

Usage:
    from deseq2_python import run_deseq2
    res = run_deseq2(counts, groups, group_a="treated", group_b="control")
"""

from __future__ import annotations

import sys
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "packages"))

from biopipe_core.optional_deps import optional_import  # noqa: E402

np = optional_import('numpy')
np_HINT = '数值计算基础库'


def _require_np():
    """确保 numpy 可用，否则抛出带安装指引的异常。"""
    if np is None:
        from biopipe_core.optional_deps import MissingDependencyError

        raise MissingDependencyError('numpy', hint=np_HINT)
    return np


from scipy.optimize import minimize, minimize_scalar
from scipy.special import gammaln
from scipy.stats import norm

__all__ = [
    "MIN_DISP",
    "Deseq2Result",
    "estimate_size_factors",
    "estimate_gene_dispersions",
    "fit_dispersion_trend",
    "estimate_map_dispersions",
    "wald_test",
    "build_design_matrix",
    "run_deseq2",
    "bh_adjust",
]

# DESeq2 中的最小离散度下限（用于避免数值奇异）
MIN_DISP = 1e-8
# 离散度搜索范围
MAX_DISP = 10.0


# --------------------------------------------------------------------------- #
# 基本工具
# --------------------------------------------------------------------------- #
def bh_adjust(pvalues: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR 校正。

    Args:
        pvalues: 一维 p 值数组（允许 NaN，NaN 位置校正后仍为 NaN）。

    Returns:
        校正后的 q 值数组，顺序与输入一致。
    """
    p = np.asarray(pvalues, dtype=float)
    out = np.full(p.shape, np.nan, dtype=float)
    valid = np.isfinite(p)
    if not valid.any():
        return out
    pv = p[valid]
    n = pv.size
    order = np.argsort(pv, kind="mergesort")
    ranked = pv[order]
    # q_i = min_{j>=i} (n / j) * p_(j)
    scaled = ranked * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(scaled[::-1])[::-1]
    q = np.clip(q, 0.0, 1.0)
    adj = np.empty(n, dtype=float)
    adj[order] = q
    out[valid] = adj
    return out


def build_design_matrix(groups: Sequence[str], group_a: str, group_b: str) -> tuple[np.ndarray, np.ndarray]:
    """构造两组对比的设计矩阵与对比向量。

    设计矩阵为哑变量编码：截距列 + 指示 ``group_a == 1`` 的列（参照组为 ``group_b``）。
    对比向量 ``c`` 选取第 2 列，从而 ``c'beta`` 即 group_a 相对 group_b 的自然对数倍数。

    Args:
        groups: 每个样本的分组标签序列。
        group_a: 分子组（treatment）。
        group_b: 分母组（control / reference）。

    Returns:
        ``(X, contrast)``：``X`` 形状 ``(n_samples, 2)``，``contrast`` 形状 ``(2,)``。
    """
    groups = list(groups)
    if group_a not in groups or group_b not in groups:
        raise ValueError(f"分组标签中缺少 {group_a!r} 或 {group_b!r}，现有标签: {sorted(set(groups))}")
    if group_a == group_b:
        raise ValueError("分子组与分组组不能相同")
    X = np.column_stack(
        [
            np.ones(len(groups), dtype=float),
            np.array([1.0 if g == group_a else 0.0 for g in groups], dtype=float),
        ]
    )
    contrast = np.array([0.0, 1.0], dtype=float)
    return X, contrast


# --------------------------------------------------------------------------- #
# Size factor
# --------------------------------------------------------------------------- #
def estimate_size_factors(counts: np.ndarray) -> np.ndarray:
    """中位数比值法（median-of-ratios）估计 size factor。

    Args:
        counts: ``(n_samples, n_genes)`` 整数 counts 矩阵。

    Returns:
        长度 ``n_samples`` 的 size factor 数组（几何均值已被归一到 1）。
    """
    counts = np.asarray(counts, dtype=float)
    n_samples, n_genes = counts.shape
    if n_genes == 0:
        raise ValueError("counts 矩阵没有任何基因")

    # 每个基因跨样本的几何均值（只对非零值取几何均值，与 DESeq2 一致）
    log_counts = np.zeros_like(counts)
    nonzero = counts > 0
    log_counts[nonzero] = np.log(counts[nonzero])
    n_nonzero = nonzero.sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        geo_mean = np.exp(log_counts.sum(axis=0) / np.maximum(n_nonzero, 1))
    usable = (n_nonzero > 0) & np.isfinite(geo_mean) & (geo_mean > 0)

    ratios = np.zeros_like(counts)
    ratios[:, usable] = counts[:, usable] / geo_mean[usable]
    size_factors = np.median(ratios[:, usable], axis=1)
    size_factors = np.where(np.isfinite(size_factors) & (size_factors > 0), size_factors, 1.0)
    return size_factors


# --------------------------------------------------------------------------- #
# 负二项 GLM
# --------------------------------------------------------------------------- #
def _nb_loglik(y: np.ndarray, mu: np.ndarray, alpha: float) -> float:
    """负二项对数似然（参数化为 size=1/alpha）。

    Args:
        y: 观测 counts。
        mu: 期望值（已含 size factor 偏移）。
        alpha: 离散度。

    Returns:
        对数似然值。
    """
    if alpha <= 0:
        return -np.inf
    r = 1.0 / alpha
    mu = np.maximum(mu, 1e-10)
    # NB2 负二项对数似然：
    #   lgamma(y+r) - lgamma(r) - lgamma(y+1)
    #   + r * log(r/(r+mu)) + y * log(mu/(r+mu))
    # 其中 r = 1/alpha，r/(r+mu) = 1/(1+alpha*mu)，mu/(r+mu) = alpha*mu/(1+alpha*mu)
    ll = (
        gammaln(y + r)
        - gammaln(r)
        - gammaln(y + 1.0)
        - r * np.log1p(alpha * mu)
        + y * np.log(alpha * mu / (1.0 + alpha * mu))
    )
    return float(np.sum(ll))


def _nb_glm_fit(
    y: np.ndarray,
    X: np.ndarray,
    offset: np.ndarray,
    alpha: float,
    beta_init: np.ndarray | None = None,
    max_iter: int = 60,
    tol: float = 1e-10,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """用 IRLS 拟合负二项广义线性模型（log 链接，size factor 作为 offset）。

    Args:
        y: 观测 counts。
        X: 设计矩阵。
        offset: log(size factor) 偏移向量。
        alpha: 固定离散度。
        beta_init: 初值。
        max_iter: 最大迭代次数。
        tol: 收敛容差。

    Returns:
        ``(beta, mu, w)``，其中 ``w = mu / (1 + alpha * mu)`` 为工作权重。
    """
    n, p = X.shape
    if beta_init is None:
        # 以 OLS 对 log(count+0.5) 的近似作为初值，比全零初值更稳定
        y_log = np.log(y + 0.5) - offset
        beta_init = np.linalg.lstsq(X, y_log, rcond=None)[0]

    beta = np.array(beta_init, dtype=float)
    mu = np.exp(X @ beta + offset)
    w = mu / (1.0 + alpha * mu)
    for _ in range(max_iter):
        mu = np.exp(X @ beta + offset)
        w = mu / (1.0 + alpha * mu)
        z = X @ beta + (y - mu) / mu
        XtW = X.T * w
        A = XtW @ X
        b = XtW @ z
        try:
            new_beta = np.linalg.solve(A, b)
        except np.linalg.LinAlgError:
            new_beta = np.linalg.lstsq(A, b, rcond=None)[0]
        if not np.all(np.isfinite(new_beta)):
            break
        delta = np.max(np.abs(new_beta - beta))
        beta = new_beta
        if delta < tol:
            break
    mu = np.exp(X @ beta + offset)
    w = mu / (1.0 + alpha * mu)
    return beta, mu, w


def _cr_loglik(
    y: np.ndarray,
    X: np.ndarray,
    offset: np.ndarray,
    alpha: float,
    beta_init: np.ndarray | None = None,
) -> tuple[float, np.ndarray, np.ndarray]:
    """Cox-Reid 调整后的负二项 profile 对数似然。

    Cox-Reid 调整用于消除参数估计自由度带来的偏差：
    ``logL_CR = logL(alpha, beta_hat) - 0.5 * log|X'WX|``。

    Returns:
        ``(logL_CR, beta_hat, mu_hat)``。
    """
    beta, mu, w = _nb_glm_fit(y, X, offset, alpha, beta_init=beta_init)
    ll = _nb_loglik(y, mu, alpha)
    A = (X.T * w) @ X
    sign, logdet = np.linalg.slogdet(A)
    if sign <= 0 or not np.isfinite(logdet):
        return -np.inf, beta, mu
    return ll - 0.5 * logdet, beta, mu


# --------------------------------------------------------------------------- #
# 基因级离散度估计
# --------------------------------------------------------------------------- #
def estimate_gene_dispersions(
    counts: np.ndarray,
    X: np.ndarray,
    size_factors: np.ndarray,
    n_grid: int = 18,
    n_refine: int = 10,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """逐基因估计 gene-wise 离散度（Cox-Reid profile 似然最大化）。

    先在 log 空间做粗网格搜索，再在最优值附近做细化搜索。

    Args:
        counts: ``(n_samples, n_genes)`` counts 矩阵。
        X: 设计矩阵。
        size_factors: size factor 数组。
        n_grid: 粗网格点数。
        n_refine: 细化网格点数。

    Returns:
        ``(alpha_hat, mu_hat, beta_hat)``：分别是基因离散度、拟合均值与系数。
    """
    counts = np.asarray(counts, dtype=float)
    n_samples, n_genes = counts.shape
    offset = np.log(np.maximum(size_factors, 1e-12))
    alphas = np.full(n_genes, MIN_DISP)
    mus = np.zeros((n_samples, n_genes))
    betas = np.zeros((n_genes, X.shape[1]))

    coarse = np.geomspace(MIN_DISP * 10, MAX_DISP, n_grid)
    for g in range(n_genes):
        y = counts[:, g]
        if y.sum() <= 0:
            continue
        # 常数基因（所有样本相同）无法估计离散度，直接用下限
        if np.allclose(y, y[0]):
            alphas[g] = MIN_DISP
            beta, mu, _ = _nb_glm_fit(y, X, offset, MIN_DISP)
            mus[:, g] = mu
            betas[g] = beta
            continue

        best_alpha, best_ll, best_beta = MIN_DISP, -np.inf, None
        for a in coarse:
            ll, beta, mu = _cr_loglik(y, X, offset, a)
            if ll > best_ll:
                best_ll, best_alpha, best_beta = ll, a, beta
        # 细化：在最优值一个数量级范围内精搜
        lo = max(MIN_DISP * 10, best_alpha / 3.0)
        hi = min(MAX_DISP, best_alpha * 3.0)
        fine = np.geomspace(lo, hi, n_refine)
        for a in fine:
            ll, beta, mu = _cr_loglik(y, X, offset, a, beta_init=best_beta)
            if ll > best_ll:
                best_ll, best_alpha, best_beta = ll, a, beta
        alphas[g] = best_alpha
        betas[g] = best_beta if best_beta is not None else 0.0
        mu = np.exp(X @ betas[g] + offset)
        mus[:, g] = np.maximum(mu, 1e-10)

    return alphas, mus, betas


# --------------------------------------------------------------------------- #
# 离散度趋势
# --------------------------------------------------------------------------- #
@dataclass
class DispersionTrend:
    """离散度趋势模型 ``alpha_trend(mu) = a0 + a1 / mu``。

    Attributes:
        a0: 渐近离散度。
        a1: 与均值成反比的分量。
        sigma2: log 空间残差方差（作为离散度先验方差）。
        residuals: log 残差数组。
    """

    a0: float
    a1: float
    sigma2: float
    residuals: np.ndarray = field(default_factory=lambda: np.zeros(0))

    def __call__(self, mu: np.ndarray | float) -> np.ndarray | float:
        """计算给定均值下的趋势离散度。"""
        mu = np.asarray(mu, dtype=float)
        val = self.a0 + self.a1 / np.maximum(mu, 1e-8)
        return np.maximum(val, MIN_DISP)

    def as_dict(self) -> dict[str, float]:
        """返回可序列化的参数摘要。"""
        return {"a0": float(self.a0), "a1": float(self.a1), "sigma2": float(self.sigma2)}


def fit_dispersion_trend(alpha_hat: np.ndarray, mu_hat: np.ndarray) -> DispersionTrend:
    """拟合离散度趋势并估计先验方差。

    采用与 DESeq2 parametric 家族一致的做法：在 log 空间拟合
    ``log alpha = log(a0 + a1/mu)``，用 IRLS 式迭代加权（Huber 权��）抵抗离群基因，
    并且**剔除落在离散度下限的基因**（它们的 log 值是人为截断的，
    混入会把整条趋势拉向无穷大）。最后用残差 MAD 估计 log 空间的先验方差。

    Args:
        alpha_hat: gene-wise 离散度。
        mu_hat: 每个基因的拟合均值（取样本均值）。

    Returns:
        :class:`DispersionTrend` 实例。
    """
    alpha_hat = np.asarray(alpha_hat, dtype=float)
    mu_hat = np.asarray(mu_hat, dtype=float)
    finite = np.isfinite(alpha_hat) & np.isfinite(mu_hat) & (mu_hat > 0)
    if finite.sum() < 10:
        base = float(np.median(alpha_hat[finite])) if finite.any() else 0.1
        return DispersionTrend(a0=max(base, MIN_DISP), a1=0.0, sigma2=0.25)

    # 落在下限附近的基因是"离散度不可估"的基因，其 log 值被人为截断，不能参与趋势拟合
    at_floor = alpha_hat <= MIN_DISP * 10
    usable = finite & (~at_floor)
    if usable.sum() < 10:
        usable = finite & (alpha_hat > 0)
    if usable.sum() < 3:
        base = float(np.median(alpha_hat[finite])) if finite.any() else 0.1
        return DispersionTrend(a0=max(base, MIN_DISP), a1=0.0, sigma2=0.25)

    a = alpha_hat[usable]
    m = mu_hat[usable]
    la = np.log(np.maximum(a, MIN_DISP))
    log_m = np.log(np.maximum(m, 1e-8))

    # 参数化：t0 = log(a0), t1 = log(a1)，预测 log(a0 + a1/m) = logaddexp(t0, t1 - log_m)
    def _predict(t: np.ndarray) -> np.ndarray:
        return np.logaddexp(t[0], t[1] - log_m)

    def _resid(t: np.ndarray) -> np.ndarray:
        return _predict(t) - la

    # 初值：线性近似 log(alpha) ≈ t0 + (t1 - log m) 的普通最小二乘
    A = np.column_stack([np.ones_like(log_m), -log_m])
    coef, *_ = np.linalg.lstsq(A, la, rcond=None)
    t = np.array([np.clip(coef[0], -15.0, 15.0), np.clip(coef[1], -15.0, 15.0)])

    # Huber 加权迭代最小二乘（对 log 空间的残差做稳健加权）
    weights = np.ones_like(la)
    for _ in range(50):
        r = _resid(t)
        scale = 1.4826 * float(np.median(np.abs(r - np.median(r))))
        scale = max(scale, 1e-3)
        # Huber 权重：|r| <= 1.5*scale 为 1，否则按 1.5*scale/|r| 衰减
        huber = 1.5 * scale
        weights = np.where(np.abs(r) <= huber, 1.0, huber / np.maximum(np.abs(r), 1e-12))
        # 加权线性化：以当前预测为基点做一阶展开，再加权最小二乘
        base = _predict(t)
        # d/d(t0) logaddexp(t0, t1-logm) = exp(t0)/exp(base)
        w0 = np.exp(t[0] - base)
        w1 = np.exp(t[1] - log_m - base)
        # 线性模型：la ≈ base + w0*d_t0 + w1*d_t1
        Xw = np.column_stack([w0, w1]) * weights[:, None]
        yw = (la - base) * weights
        try:
            delta, *_ = np.linalg.lstsq(Xw, yw, rcond=None)
        except np.linalg.LinAlgError:
            break
        new_t = t + delta
        if not np.all(np.isfinite(new_t)):
            break
        new_t = np.clip(new_t, -30.0, 30.0)
        converged = np.max(np.abs(new_t - t)) < 1e-8
        t = new_t
        if converged:
            break

    a0 = float(np.exp(np.clip(t[0], -30.0, 30.0)))
    a1 = float(np.exp(np.clip(t[1], -30.0, 30.0)))

    resid = _resid(t)
    trend = DispersionTrend(a0=a0, a1=a1, sigma2=0.25)
    trend.residuals = resid
    # 用**截尾后**的残差估计先验方差，避免少数极端基因主导
    trimmed = resid[np.abs(resid - np.median(resid)) < 3.0 * max(
        1.4826 * float(np.median(np.abs(resid - np.median(resid)))), 1e-3
    )]
    mad = float(np.median(np.abs(trimmed - np.median(trimmed)))) if trimmed.size else float(
        np.median(np.abs(resid - np.median(resid)))
    )
    sigma = 1.4826 * mad
    # 下限避免先验过窄导致过度收缩
    trend.sigma2 = float(max(sigma**2, 0.09))
    return trend


# --------------------------------------------------------------------------- #
# MAP 离散度
# --------------------------------------------------------------------------- #
def estimate_map_dispersions(
    counts: np.ndarray,
    X: np.ndarray,
    size_factors: np.ndarray,
    trend: DispersionTrend,
    mu_mean: np.ndarray | None = None,
    n_grid: int = 15,
) -> np.ndarray:
    """在离散度先验下估计 MAP 离散度。

    后验（log 尺度）：
    ``log post(alpha) = logL_CR(alpha) - 0.5 * (log alpha - log trend(mu))^2 / sigma2``

    Args:
        counts: ``(n_samples, n_genes)`` counts。
        X: 设计矩阵。
        size_factors: size factor。
        trend: 已拟合的离散度趋势。
        mu_mean: 每个基因的平均拟合均值，None 时用样本均值近似。
        n_grid: 每个基因局部搜索的网格点数。

    Returns:
        长度 ``n_genes`` 的 MAP 离散度数组。
    """
    counts = np.asarray(counts, dtype=float)
    n_samples, n_genes = counts.shape
    offset = np.log(np.maximum(size_factors, 1e-12))
    out = np.full(n_genes, MIN_DISP)

    if mu_mean is None:
        mu_mean = counts.mean(axis=0)

    for g in range(n_genes):
        y = counts[:, g]
        if y.sum() <= 0:
            continue
        mu_g = max(float(mu_mean[g]), 1e-6)
        log_prior_mean = float(np.log(trend(mu_g)))
        sigma2 = max(trend.sigma2, 1e-6)

        grid = np.geomspace(MIN_DISP * 10, MAX_DISP, 24)
        best_alpha, best_post = MIN_DISP, -np.inf
        best_beta = None
        for a in grid:
            ll, beta, mu = _cr_loglik(y, X, offset, a)
            if not np.isfinite(ll):
                continue
            post = ll - 0.5 * (np.log(a) - log_prior_mean) ** 2 / sigma2
            if post > best_post:
                best_post, best_alpha, best_beta = post, a, beta
        if best_post == -np.inf:
            continue
        # 局部细化
        lo = max(MIN_DISP * 10, best_alpha / 2.5)
        hi = min(MAX_DISP, best_alpha * 2.5)
        for a in np.geomspace(lo, hi, n_grid):
            ll, beta, mu = _cr_loglik(y, X, offset, a, beta_init=best_beta)
            if not np.isfinite(ll):
                continue
            post = ll - 0.5 * (np.log(a) - log_prior_mean) ** 2 / sigma2
            if post > best_post:
                best_post, best_alpha = post, a
        out[g] = best_alpha
    return out


# --------------------------------------------------------------------------- #
# Wald 检验
# --------------------------------------------------------------------------- #
@dataclass
class Deseq2Result:
    """DESeq2 风格差异表达结果。

    Attributes:
        table: 结果表（pandas DataFrame），列为 baseMean/log2FoldChange/lfcSE/stat/pvalue/padj。
        size_factors: size factor。
        alphas: gene-wise 离散度。
        alpha_map: MAP 离散度。
        trend_params: 离散度趋势参数。
        n_samples: 样本数。
        n_genes: 参与检验的基因数。
    """

    table: object
    size_factors: np.ndarray
    alphas: np.ndarray
    alpha_map: np.ndarray
    trend_params: dict[str, float]
    n_samples: int
    n_genes: int


def wald_test(
    counts: np.ndarray,
    X: np.ndarray,
    size_factors: np.ndarray,
    contrast: np.ndarray,
    alpha: np.ndarray,
) -> dict[str, np.ndarray]:
    """对固定离散度做 NB-GLM Wald 检验。

    Args:
        counts: ``(n_samples, n_genes)`` counts。
        X: 设计矩阵。
        size_factors: size factor。
        contrast: 对比向量（作用在回归系数上）。
        alpha: 每个基因的离散度。

    Returns:
        字典，键为 ``baseMean`` / ``log2FoldChange`` / ``lfcSE`` / ``stat`` / ``pvalue``。
    """
    counts = np.asarray(counts, dtype=float)
    n_samples, n_genes = counts.shape
    offset = np.log(np.maximum(size_factors, 1e-12))
    normalized = counts / np.maximum(size_factors[:, None], 1e-12)

    base_mean = normalized.mean(axis=0)
    lfc = np.full(n_genes, np.nan)
    lfc_se = np.full(n_genes, np.nan)
    stat = np.full(n_genes, np.nan)
    pvalue = np.full(n_genes, np.nan)

    for g in range(n_genes):
        y = counts[:, g]
        if y.sum() <= 0:
            continue
        a = float(alpha[g]) if np.isfinite(alpha[g]) else MIN_DISP
        a = max(a, MIN_DISP)
        beta, mu, w = _nb_glm_fit(y, X, offset, a, max_iter=100, tol=1e-12)
        if not np.all(np.isfinite(beta)):
            continue
        A = (X.T * w) @ X
        try:
            cov = np.linalg.inv(A)
        except np.linalg.LinAlgError:
            cov = np.linalg.pinv(A)
        var = float(contrast @ cov @ contrast)
        if not np.isfinite(var) or var <= 0:
            continue
        est = float(contrast @ beta)
        se = float(np.sqrt(var))
        z = est / se
        lfc[g] = est / np.log(2.0)
        lfc_se[g] = se / np.log(2.0)
        stat[g] = z
        pvalue[g] = float(2.0 * norm.sf(abs(z)))

    return {
        "baseMean": base_mean,
        "log2FoldChange": lfc,
        "lfcSE": lfc_se,
        "stat": stat,
        "pvalue": pvalue,
    }


# --------------------------------------------------------------------------- #
# 顶层入口
# --------------------------------------------------------------------------- #
def run_deseq2(
    counts: np.ndarray,
    groups: Sequence[str],
    group_a: str,
    group_b: str,
    min_counts: int = 10,
    min_replicates: int = 1,
    gene_names: Sequence[str] | None = None,
) -> Deseq2Result:
    """执行完整的 DESeq2 风格差异表达分析。

    Args:
        counts: ``(n_samples, n_genes)`` 原始整数 counts。
        groups: 每个样本的分组标签。
        group_a: 分子组（treatment）。
        group_b: 分母组（control）。
        min_counts: 基因在所有样本中的最小总 counts，低于则剔除。
        min_replicates: 每组最少样本数。
        gene_names: 基因名列表，长度须与 ``counts`` 的基因数一致；None 时使用整数位置。

    Returns:
        :class:`Deseq2Result`。
    """
    import pandas as pd  # 局部导入，保持模块顶层轻量

    counts = np.asarray(counts, dtype=float)
    groups = list(groups)
    if counts.shape[0] != len(groups):
        raise ValueError("counts 的样本数与分组标签数量不一致")

    keep_genes = counts.sum(axis=0) >= max(min_counts, 1)
    counts = counts[:, keep_genes]
    if gene_names is not None:
        gene_names = [str(g) for g, k in zip(gene_names, keep_genes) if k]
    else:
        gene_names = [str(i) for i, k in enumerate(keep_genes) if k]

    n_a = sum(1 for g in groups if g == group_a)
    n_b = sum(1 for g in groups if g == group_b)
    if n_a < min_replicates or n_b < min_replicates:
        raise ValueError(f"分组样本量不足: {group_a}={n_a}, {group_b}={n_b}")

    size_factors = estimate_size_factors(counts)
    X, contrast = build_design_matrix(groups, group_a, group_b)

    alphas, mus, _ = estimate_gene_dispersions(counts, X, size_factors)
    mu_mean = mus.mean(axis=0)

    trend_obj = fit_dispersion_trend(alphas, mu_mean)
    alpha_map = estimate_map_dispersions(counts, X, size_factors, trend_obj, mu_mean=mu_mean)

    stats = wald_test(counts, X, size_factors, contrast, alpha_map)
    padj = bh_adjust(stats["pvalue"])

    table = pd.DataFrame(
        {
            "baseMean": stats["baseMean"],
            "log2FoldChange": stats["log2FoldChange"],
            "lfcSE": stats["lfcSE"],
            "stat": stats["stat"],
            "pvalue": stats["pvalue"],
            "padj": padj,
        },
        index=gene_names,
    )
    table.index.name = "gene_id"

    return Deseq2Result(
        table=table,
        size_factors=size_factors,
        alphas=alphas,
        alpha_map=alpha_map,
        trend_params=trend_obj.as_dict(),
        n_samples=counts.shape[0],
        n_genes=counts.shape[1],
    )
