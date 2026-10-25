"""分子动力学轨迹的核心分析算法（纯 numpy / scipy 实现）。

包含：

* :func:`kabsch_align` / :func:`rmsd` — Kabsch-Sander 最优刚体叠合与 RMSD
* :func:`rmsf` — 逐原子涨落（均方根涨落）
* :func:`radius_of_gyration` — 质量加权回旋半径
* :func:`pca_analysis` — 对坐标协方差做 SVD 主成分分析
* :func:`hbond_analysis` — 基于距离 + 角度的几何氢键判据
* :func:`dihedral_series` / :func:`ramachandran` — 二面角与 Ramachandran 分布
* :func:`rmsd_matrix` / :func:`hierarchical_cluster` — RMSD 矩阵与层次聚类
* :func:`interchain_contacts` — 链间相互作用
* :func:`sasa_estimate` — Shrake-Rupley 溶剂可及表面积

若 MDAnalysis 可用，入口脚本优先使用其内建实现；本模块作为不依赖
MDAnalysis 的等价回退。
"""

from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Sequence

import numpy as np

__all__ = [
    "kabsch_align",
    "rmsd",
    "rmsf",
    "radius_of_gyration",
    "pca_analysis",
    "hbond_analysis",
    "dihedral_series",
    "ramachandran",
    "rmsd_matrix",
    "hierarchical_cluster",
    "interchain_contacts",
    "sasa_estimate",
    "backbone_dihedral",
]


# ---------------------------------------------------------------------------
# Kabsch 叠合与 RMSD
# ---------------------------------------------------------------------------


def kabsch_align(
    mobile: np.ndarray,
    reference: np.ndarray,
    weights: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Kabsch-Sander 最优刚体叠合。

    求旋转矩阵 R 与平移 t，使得 ``mobile @ R + t`` 尽量接近 ``reference``。
    使用奇异值分解求解：设 ``C = (mobile_centered)^T @ W @ reference_centered``，
    ``SVD(C) = U S V^T``，则 ``R = V U^T``（行列式为负时翻转最后一个奇异向量）。

    Args:
        mobile: ``(n, 3)`` 待叠合坐标。
        reference: ``(n, 3)`` 参考坐标。
        weights: ``(n,)`` 权重（质量），None 表示等权。

    Returns:
        ``(rot, trans, rmsd)``：旋转矩阵 ``(3,3)``、平移向量 ``(3,)``、叠合后 RMSD。
    """
    mobile = np.asarray(mobile, dtype=np.float64)
    reference = np.asarray(reference, dtype=np.float64)
    if mobile.shape != reference.shape:
        raise ValueError("mobile 与 reference 形状不一致")
    if mobile.ndim != 2 or mobile.shape[1] != 3:
        raise ValueError("坐标数组形状必须是 (n, 3)")
    if mobile.shape[0] < 3:
        raise ValueError("叠合至少需要 3 个原子")

    if weights is None:
        w = np.ones(mobile.shape[0], dtype=np.float64)
    else:
        w = np.asarray(weights, dtype=np.float64)
        if w.shape[0] != mobile.shape[0]:
            raise ValueError("权重长度与原子数不一致")
    total = float(w.sum())
    if total <= 0:
        w = np.ones(mobile.shape[0], dtype=np.float64)
        total = float(w.sum())
    w = w / total

    m_center = (mobile * w[:, None]).sum(axis=0)
    r_center = (reference * w[:, None]).sum(axis=0)
    m_shift = mobile - m_center
    r_shift = reference - r_center

    # 协方差矩阵 H = Σ wᵢ p̃ᵢ q̃ᵢᵀ
    covariance = (m_shift * w[:, None]).T @ r_shift
    u, s, vt = np.linalg.svd(covariance)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    diag = np.diag([1.0, 1.0, d])
    rot = vt.T @ diag @ u.T

    aligned = m_shift @ rot.T + r_center
    diff = aligned - reference
    rmsd = float(np.sqrt(np.sum(w[:, None] * diff * diff)))

    # 平移向量：把叠合后的坐标还原为对原始 mobile 的变换
    trans = r_center - m_center @ rot.T
    return rot, trans, rmsd


def rmsd(
    coords_a: np.ndarray,
    coords_b: np.ndarray,
    align: bool = True,
    weights: np.ndarray | None = None,
) -> float:
    """计算两组坐标的 RMSD。

    Args:
        coords_a: ``(n, 3)``。
        coords_b: ``(n, 3)``。
        align: 是否先做最优刚体叠合。
        weights: 权重。

    Returns:
        RMSD（Å）。
    """
    if align:
        _, _, value = kabsch_align(coords_a, coords_b, weights)
        return value
    diff = np.asarray(coords_a, dtype=np.float64) - np.asarray(coords_b, dtype=np.float64)
    if weights is None:
        return float(np.sqrt(np.mean(np.sum(diff * diff, axis=1))))
    w = np.asarray(weights, dtype=np.float64)
    w = w / w.sum()
    return float(np.sqrt(np.sum(w * np.sum(diff * diff, axis=1))))


# ---------------------------------------------------------------------------
# RMSF / 回旋半径
# ---------------------------------------------------------------------------


def rmsf(coords: np.ndarray, align_reference: int | None = 0) -> np.ndarray:
    """逐原子均方根涨落。

    Args:
        coords: ``(n_frames, n_atoms, 3)``。
        align_reference: 若给定，先把所有帧叠合到该帧（Kabsch），再计算涨落。

    Returns:
        ``(n_atoms,)`` 每个原子的 RMSF（Å）。
    """
    coords = np.asarray(coords, dtype=np.float64)
    if coords.ndim != 3:
        raise ValueError("coords 形状必须是 (n_frames, n_atoms, 3)")
    if align_reference is not None and coords.shape[0] > 1:
        ref = coords[align_reference]
        aligned = np.empty_like(coords)
        aligned[0] = coords[0]
        for i in range(1, coords.shape[0]):
            rot, trans, _ = kabsch_align(coords[i], ref)
            aligned[i] = coords[i] @ rot.T + trans
        coords = aligned
    mean = coords.mean(axis=0)
    return np.sqrt(((coords - mean) ** 2).sum(axis=2).mean(axis=0))


def radius_of_gyration(coords: np.ndarray, weights: np.ndarray | None = None) -> np.ndarray:
    """逐帧质量加权回旋半径。

    Args:
        coords: ``(n_frames, n_atoms, 3)``。
        weights: ``(n_atoms,)`` 质量权重。

    Returns:
        ``(n_frames,)`` 每帧的 Rg（Å）。
    """
    coords = np.asarray(coords, dtype=np.float64)
    if weights is None:
        w = np.ones(coords.shape[1], dtype=np.float64)
    else:
        w = np.asarray(weights, dtype=np.float64)
    w = w / w.sum()
    center = (coords * w[None, :, None]).sum(axis=1)
    delta = coords - center[:, None, :]
    r2 = (w[None, :] * (delta ** 2).sum(axis=2)).sum(axis=1)
    return np.sqrt(r2)


# ---------------------------------------------------------------------------
# PCA
# ---------------------------------------------------------------------------


def pca_analysis(coords: np.ndarray, n_components: int = 3) -> dict:
    """对轨迹坐标做主成分分析（协方差矩阵的 SVD）。

    步骤：

    1. 选择所有帧作为样本，展平为 ``(n_frames, 3 * n_atoms)``；
    2. 减去均值；
    3. 对坐标做 Kabsch 叠合到第一帧（可选，由调用方决定是否预先叠合）；
    4. 协方差矩阵 ``C = X^T X / (n_frames - 1)``，做 SVD 得到主成分；
    5. 投影得到每个帧的主成分坐标与解释方差比例。

    Args:
        coords: ``(n_frames, n_atoms, 3)``。
        n_components: 返回的主成分个数。

    Returns:
        字典，含 ``eigenvalues``、``explained_variance_ratio``、``projection``
        ``(n_frames, n_components)``、``modes``（每个主成分的原子位移向量）。
    """
    coords = np.asarray(coords, dtype=np.float64)
    n_frames = coords.shape[0]
    if n_frames < 2:
        raise ValueError("PCA 至少需要 2 帧")

    flat = coords.reshape(n_frames, -1)
    centered = flat - flat.mean(axis=0, keepdims=True)
    covariance = centered.T @ centered / (n_frames - 1)
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]

    k = min(n_components, eigenvectors.shape[1])
    eigenvalues = eigenvalues[:k]
    eigenvectors = eigenvectors[:, :k]
    projection = centered @ eigenvectors

    total = float(np.sum(np.clip(eigenvalues, 0, None))) or 1.0
    ratio = np.clip(eigenvalues, 0, None) / total
    modes = eigenvectors.reshape(k, coords.shape[1], 3)

    return {
        "eigenvalues": eigenvalues,
        "explained_variance_ratio": ratio,
        "cumulative_variance": np.cumsum(ratio),
        "projection": projection,
        "modes": modes,
        "n_frames": n_frames,
    }


# ---------------------------------------------------------------------------
# 氢键（距离 + 角度几何判据）
# ---------------------------------------------------------------------------

# 常用原子半径（Å），用于范德华接触与 SASA
_VDW = {
    "H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80,
    "P": 1.80, "F": 1.47, "CL": 1.75, "BR": 1.85, "I": 1.98,
    "FE": 1.80, "ZN": 1.39, "MG": 1.73, "CA": 1.95, "NA": 2.27,
    "K": 2.75, "SE": 1.90,
}

# 常见氢键受体原子（按元素与杂化环境粗判）
_ACCEPTOR_ELEMENTS = {"N", "O", "S", "F"}
# 常见氢键供体（带 H 的 N/O）
_DONOR_ELEMENTS = {"N", "O"}


def _element_radius(element: str) -> float:
    return _VDW.get((element or "").strip().upper(), 1.70)


def hbond_analysis(
    coords: np.ndarray,
    elements: Sequence[str],
    donor_h_names: Sequence[str] = ("H", "HN", "HA", "HB2", "HB3", "HD1", "HD2", "HE", "HG1", "HG2", "HZ"),
    h_dist_max: float = 2.5,
    ha_dist_max: float = 3.5,
    angle_min: float = 120.0,
    exclude_pairs: Sequence[tuple[int, int]] = (),
    exclude_same_residue: bool = True,
) -> dict:
    """基于距离 + 角度几何判据的氢键分析。

    判据（几何近似，不做部分电荷计算）：

    * 供体 D 上存在氢原子 H（原子名以 ``H`` 开头）；
    * ``D···A`` 距离 <= ``ha_dist_max``（默认 3.5 Å）；
    * ``H···A`` 距离 <= ``h_dist_max``（默认 2.5 Å）；
    * 角度 ``D-H···A`` >= ``angle_min``（默认 120°）。

    Args:
        coords: ``(n_frames, n_atoms, 3)``。
        elements: 每个原子的元素符号，长度 = 原子数。
        donor_h_names: 视为氢原子的原子名（未列出的以 ``H`` 开头的也计入）。
        h_dist_max: H···A 最大距离。
        ha_dist_max: D···A 最大距离。
        angle_min: 最小 D-H···A 角度（度）。
        exclude_pairs: 需要排除的原子对（如 1-4 键）。
        exclude_same_residue: 是否排除同一残基内的接触。

    Returns:
        字典，含逐帧氢键数、存在性矩阵、平均存在率、检测到的氢键对列表。
    """
    coords = np.asarray(coords, dtype=np.float64)
    n_frames, n_atoms, _ = coords.shape
    elements = [str(e).strip().upper() for e in elements]

    h_names = {n.upper() for n in donor_h_names}
    is_h = np.array([e in ("H", "D") for e in elements], dtype=bool)
    acceptor_mask = np.array([e in _ACCEPTOR_ELEMENTS for e in elements], dtype=bool)
    donor_mask = np.array([e in _DONOR_ELEMENTS for e in elements], dtype=bool)

    excluded = set()
    for a, b in exclude_pairs:
        excluded.add((int(a), int(b)))
        excluded.add((int(b), int(a)))

    # 枚举 (供体重原子 D, 氢 H, 受体 A) 三元组：先用第一帧做拓扑筛选
    candidates = _expand_hydrogens(
        coords[0], donor_mask, is_h, acceptor_mask, elements, exclude_same_residue, excluded
    )

    presence = np.zeros((n_frames, len(candidates)), dtype=bool)
    for f in range(n_frames):
        pos = coords[f]
        for ci, (d, h, a) in enumerate(candidates):
            d_pos, h_pos, a_pos = pos[d], pos[h], pos[a]
            dist_ha = float(np.linalg.norm(h_pos - a_pos))
            if dist_ha > h_dist_max:
                continue
            dist_da = float(np.linalg.norm(d_pos - a_pos))
            if dist_da > ha_dist_max:
                continue
            v1 = d_pos - h_pos
            v2 = a_pos - h_pos
            denom = np.linalg.norm(v1) * np.linalg.norm(v2)
            if denom < 1e-8:
                continue
            cosang = float(np.dot(v1, v2) / denom)
            ang = math.degrees(math.acos(max(-1.0, min(1.0, cosang))))
            if ang < angle_min:
                continue
            presence[f, ci] = True

    per_frame = presence.sum(axis=1)
    occupancy = presence.mean(axis=0) if len(candidates) else np.zeros(0)
    bonds = []
    for ci, (d, h, a) in enumerate(candidates):
        occ = float(occupancy[ci]) if len(candidates) else 0.0
        if occ <= 0.0:
            continue
        bonds.append(
            {
                "donor_index": d,
                "hydrogen_index": h,
                "acceptor_index": a,
                "occupancy": round(occ, 4),
                "frames": int(presence[:, ci].sum()),
            }
        )
    bonds.sort(key=lambda b: (-b["occupancy"], b["donor_index"], b["acceptor_index"]))

    return {
        "per_frame_count": per_frame,
        "mean_count": float(per_frame.mean()) if n_frames else 0.0,
        "occupancy": bonds,
        "n_candidate_pairs": len(candidates),
        "criteria": {
            "h_dist_max": h_dist_max,
            "ha_dist_max": ha_dist_max,
            "angle_min": angle_min,
        },
    }


def _expand_hydrogens(
    ref: np.ndarray,
    donor_mask: np.ndarray,
    is_h: np.ndarray,
    acceptor_mask: np.ndarray,
    elements: list[str],
    exclude_same_residue: bool,
    excluded: set,
) -> list[tuple[int, int, int]]:
    """枚举 (供体重原子 D, 氢 H, 受体 A) 三元组。"""
    triples: list[tuple[int, int, int]] = []
    donor_indices = np.flatnonzero(donor_mask)
    h_indices = np.flatnonzero(is_h)
    acceptor_indices = np.flatnonzero(acceptor_mask)
    for d in donor_indices:
        attached = []
        for h in h_indices:
            if abs(int(h) - int(d)) > 2:
                continue
            if np.linalg.norm(ref[h] - ref[d]) > 1.35:
                continue
            attached.append(int(h))
        if not attached:
            continue
        for a in acceptor_indices:
            if int(a) == int(d):
                continue
            if exclude_same_residue and _same_residue(int(d), int(a), elements):
                continue
            if (int(d), int(a)) in excluded:
                continue
            for h in attached:
                triples.append((int(d), h, int(a)))
    return triples


_RESIDUE_KEYS: dict[int, tuple[str, int]] = {}


def _same_residue(i: int, j: int, elements: list[str]) -> bool:
    """粗略判断两个原子是否属于同一残基（依赖外部注册的残基映射）。"""
    key_i = _RESIDUE_KEYS.get(i)
    key_j = _RESIDUE_KEYS.get(j)
    if key_i is None or key_j is None:
        return False
    return key_i == key_j


def register_residue_map(residue_map: dict[int, tuple[str, int]]) -> None:
    """注册原子 -> 残基标识的映射，供氢键判据排除残基内接触。"""
    _RESIDUE_KEYS.clear()
    _RESIDUE_KEYS.update(residue_map)


# ---------------------------------------------------------------------------
# 二面角 / Ramachandran
# ---------------------------------------------------------------------------


def _dihedral(p0: np.ndarray, p1: np.ndarray, p2: np.ndarray, p3: np.ndarray) -> float:
    """计算 IUPAC 约定的二面角（度），范围 (-180, 180]。"""
    b0 = p0 - p1
    b1 = p2 - p1
    b2 = p3 - p2
    norm = np.linalg.norm(b1)
    if norm < 1e-8:
        return float("nan")
    axis = b1 / norm
    v = b0 - np.dot(b0, axis) * axis
    w = b2 - np.dot(b2, axis) * axis
    x = float(np.dot(v, w))
    cr = np.cross(axis, v)
    y = float(np.dot(cr, w))
    return math.degrees(math.atan2(y, x))


def backbone_dihedral(coords: np.ndarray, n_idx: int, ca_idx: int, c_idx: int, c_next: int) -> np.ndarray:
    """逐帧计算 phi 角 ``(C(i-1), N(i), CA(i), C(i))``。

    Args:
        coords: ``(n_frames, n_atoms, 3)``。
        n_idx: N 原子索引。
        ca_idx: CA 原子索引。
        c_idx: C 原子索引。
        c_next: 下一个残基的 C 原子索引（phi 用）。

    Returns:
        ``(n_frames,)`` 角度数组（度）。
    """
    out = np.full(coords.shape[0], np.nan)
    for f in range(coords.shape[0]):
        p = coords[f]
        out[f] = _dihedral(p[c_next], p[n_idx], p[ca_idx], p[c_idx])
    return out


def dihedral_series(
    coords: np.ndarray,
    atom_indices: Sequence[int],
    label: str = "chi",
) -> dict:
    """计算绕指定四原子序列的二面角序列。

    Args:
        coords: ``(n_frames, n_atoms, 3)``。
        atom_indices: 四个原子索引（依次为 i-1, i, i+1, i+2）。
        label: 角度标签。

    Returns:
        字典，含 ``angles`` 数组与统计量。
    """
    if len(atom_indices) != 4:
        raise ValueError("二面角需要恰好 4 个原子索引")
    a, b, c, d = atom_indices
    angles = np.full(coords.shape[0], np.nan)
    for f in range(coords.shape[0]):
        p = coords[f]
        angles[f] = _dihedral(p[a], p[b], p[c], p[d])
    valid = angles[~np.isnan(angles)]
    hist, edges = np.histogram(valid, bins=36, range=(-180.0, 180.0)) if valid.size else (
        np.zeros(36),
        np.linspace(-180.0, 180.0, 37),
    )
    return {
        "label": label,
        "angles": angles,
        "mean": float(valid.mean()) if valid.size else float("nan"),
        "std": float(valid.std()) if valid.size else float("nan"),
        "histogram": hist.tolist(),
        "bin_edges": edges.tolist(),
        "circular_mean": _circular_mean(valid),
    }


def _circular_mean(angles: np.ndarray) -> float:
    """角度的圆周均值（度）。"""
    if angles.size == 0:
        return float("nan")
    rad = np.radians(angles)
    return float(math.degrees(math.atan2(np.sin(rad).mean(), np.cos(rad).mean())))


# 常见二级结构区域的粗略判定（φ, ψ）区间
_RAMA_HELIX = [(-160.0, -20.0, -90.0, 0.0), (-120.0, 50.0, 60.0, 180.0)]
_RAMA_SHEET = [(-180.0, -50.0, 60.0, 180.0), (-150.0, -20.0, 90.0, 180.0)]


def _in_region(phi: float, psi: float, regions) -> bool:
    for p_lo, p_hi, s_lo, s_hi in regions:
        if p_lo <= phi <= p_hi and s_lo <= psi <= s_hi:
            return True
    return False


def ramachandran(
    phis: np.ndarray,
    psis: np.ndarray,
    bins: int = 18,
) -> dict:
    """Ramachandran 图分析。

    按粗略的 (φ, ψ) 区域划分为 alpha 区、beta 区、alphaL 区与 coil，
    并给出二维分布直方图（以 30° 为网格）。

    Args:
        phis: ``(n_frames,)`` phi 角。
        psis: ``(n_frames,)`` psi 角。
        bins: 每个轴的直方图 bin 数（默认 18，即 20° 一格）。

    Returns:
        字典，含区域占比与二维分布。
    """
    phis = np.asarray(phis, dtype=np.float64)
    psis = np.asarray(psis, dtype=np.float64)
    valid = ~np.isnan(phis) & ~np.isnan(psis)
    phi_v, psi_v = phis[valid], psis[valid]

    counts = {"alpha": 0, "beta": 0, "alphaL": 0, "coil": 0}
    for phi, psi in zip(phi_v, psi_v):
        if _in_region(phi, psi, _RAMA_HELIX):
            counts["alpha"] += 1
        elif _in_region(phi, psi, _RAMA_SHEET):
            counts["beta"] += 1
        elif -90.0 <= phi <= 30.0 and -60.0 <= psi <= 100.0:
            counts["alphaL"] += 1
        else:
            counts["coil"] += 1

    total = max(len(phi_v), 1)
    hist, xedges, yedges = np.histogram2d(
        phi_v, psi_v, bins=bins, range=[[-180.0, 180.0], [-180.0, 180.0]]
    )
    return {
        "n_valid": int(total),
        "fractions": {k: v / total for k, v in counts.items()},
        "counts": counts,
        "histogram": hist.tolist(),
        "phi_edges": xedges.tolist(),
        "psi_edges": yedges.tolist(),
        "outliers": int(
            np.sum(
                ~(
                    _in_region_arr(phi_v, psi_v, _RAMA_HELIX)
                    | _in_region_arr(phi_v, psi_v, _RAMA_SHEET)
                    | ((phi_v >= -90.0) & (phi_v <= 30.0) & (psi_v >= -60.0) & (psi_v <= 100.0))
                )
            )
        ),
    }


def _in_region_arr(phi: np.ndarray, psi: np.ndarray, regions) -> np.ndarray:
    mask = np.zeros(phi.shape, dtype=bool)
    for p_lo, p_hi, s_lo, s_hi in regions:
        mask |= (phi >= p_lo) & (phi <= p_hi) & (psi >= s_lo) & (psi <= s_hi)
    return mask


# ---------------------------------------------------------------------------
# 聚类
# ---------------------------------------------------------------------------


def rmsd_matrix(
    coords: np.ndarray,
    max_frames: int = 200,
    selection: np.ndarray | None = None,
) -> tuple[np.ndarray, list[int]]:
    """计算帧间 RMSD 矩阵。

    对每一对帧做 Kabsch 叠合后计算 RMSD。帧数较多时按均匀间隔抽样以控制
    O(N^2) 的计算量。

    Args:
        coords: ``(n_frames, n_atoms, 3)``。
        max_frames: 参与计算的最大帧数（超过则均匀抽样）。
        selection: 参与叠合的原子索引。

    Returns:
        ``(matrix, frame_indices)``：对称 RMSD 矩阵与对应的原始帧索引。
    """
    coords = np.asarray(coords, dtype=np.float64)
    n_frames = coords.shape[0]
    if n_frames > max_frames:
        idx = np.linspace(0, n_frames - 1, max_frames).astype(int)
        idx = np.unique(idx)
    else:
        idx = np.arange(n_frames)

    sub = coords[idx]
    if selection is not None:
        sub = sub[:, selection, :]

    m = len(idx)
    matrix = np.zeros((m, m), dtype=np.float64)
    for i in range(m):
        for j in range(i + 1, m):
            value = rmsd(sub[i], sub[j], align=True)
            matrix[i, j] = value
            matrix[j, i] = value
    return matrix, idx.tolist()


def hierarchical_cluster(
    matrix: np.ndarray,
    n_clusters: int | None = None,
    threshold: float | None = None,
    method: str = "average",
) -> dict:
    """对 RMSD 矩阵做层次聚类。

    优先使用 ``scipy.cluster.hierarchy``（linkage + fcluster）；
    若 scipy 不可用则退化为简单的阈值连通分量划分。

    Args:
        matrix: ``(n, n)`` 对称距离矩阵。
        n_clusters: 期望的簇数。
        threshold: 距离阈值（与 ``n_clusters`` 二选一）。
        method: linkage 方法（``average`` / ``complete`` / ``ward`` / ``single``）。

    Returns:
        字典，含每个帧的簇标签、簇数量与各簇成员。
    """
    matrix = np.asarray(matrix, dtype=np.float64)
    n = matrix.shape[0]
    if n == 0:
        return {"labels": [], "n_clusters": 0, "clusters": []}

    labels = None
    n_clusters_final = n_clusters or 0
    try:
        from scipy.cluster.hierarchy import fcluster, linkage  # type: ignore
        from scipy.spatial.distance import squareform  # type: ignore

        condensed = squareform(matrix, checks=False)
        tree = linkage(condensed, method=method)
        if threshold is not None:
            labels = fcluster(tree, t=threshold, criterion="distance")
        else:
            k = max(1, n_clusters or 2)
            labels = fcluster(tree, t=k, criterion="maxclust")
        n_clusters_final = int(len(np.unique(labels)))
    except Exception:  # noqa: BLE001 - scipy 不可用时退化实现
        threshold_value = threshold if threshold is not None else float(np.median(matrix))
        labels = np.zeros(n, dtype=int)
        current = 0
        for i in range(n):
            if labels[i]:
                continue
            current += 1
            labels[i] = current
            stack = [i]
            while stack:
                cur = stack.pop()
                for j in range(n):
                    if not labels[j] and matrix[cur, j] <= threshold_value:
                        labels[j] = current
                        stack.append(j)
        n_clusters_final = current

    clusters: dict[int, list[int]] = {}
    for idx, lab in enumerate(labels):
        clusters.setdefault(int(lab), []).append(idx)
    cluster_list = [
        {"id": int(lab), "size": len(members), "members": members}
        for lab, members in sorted(clusters.items())
    ]
    return {
        "labels": [int(x) for x in labels],
        "n_clusters": n_clusters_final,
        "clusters": cluster_list,
        "method": method,
        "threshold": threshold,
    }


# ---------------------------------------------------------------------------
# 链间相互作用
# ---------------------------------------------------------------------------


def interchain_contacts(
    coords: np.ndarray,
    chain_ids: Sequence[str],
    contact_dist: float = 4.5,
    min_occupancy: float = 0.5,
) -> dict:
    """统计链间原子接触。

    对每一对不同链，计算所有原子对的最小距离；距离小于阈值的记为一次接触，
    并按帧统计存在率。

    Args:
        coords: ``(n_frames, n_atoms, 3)``。
        chain_ids: 每个原子的链 ID。
        contact_dist: 接触距离阈值（Å）。
        min_occupancy: 报告时要求的最小存在率。

    Returns:
        字典，含每对链的接触数量与存在率。
    """
    coords = np.asarray(coords, dtype=np.float64)
    chain_arr = np.asarray([str(c) for c in chain_ids])
    if chain_arr.shape[0] != coords.shape[1]:
        raise ValueError(
            f"chain_ids 长度 ({chain_arr.shape[0]}) 与原子数 ({coords.shape[1]}) 不一致"
        )
    chains = sorted(set(chain_arr.tolist()))

    n_frames = coords.shape[0]
    occupancy_counter: dict[tuple[str, str], dict[tuple[int, int], int]] = {}

    for f in range(n_frames):
        pos = coords[f]
        for i in range(len(chains)):
            for j in range(i + 1, len(chains)):
                ci, cj = chains[i], chains[j]
                idx_i = np.flatnonzero(chain_arr == ci)
                idx_j = np.flatnonzero(chain_arr == cj)
                if idx_i.size == 0 or idx_j.size == 0:
                    continue
                d = np.linalg.norm(pos[idx_i][:, None, :] - pos[idx_j][None, :, :], axis=2)
                hit = np.argwhere(d < contact_dist)
                if hit.size == 0:
                    continue
                key = (ci, cj)
                occ = occupancy_counter.setdefault(key, {})
                for a_i, a_j in hit:
                    pair = (int(idx_i[a_i]), int(idx_j[a_j]))
                    occ[pair] = occ.get(pair, 0) + 1

    summary = []
    for (ci, cj), occ in sorted(occupancy_counter.items()):
        persistent = [(p, c / n_frames) for p, c in occ.items() if c / n_frames >= min_occupancy]
        if not persistent:
            continue
        persistent.sort(key=lambda x: -x[1])
        summary.append(
            {
                "chain_a": ci,
                "chain_b": cj,
                "n_contacts": len(persistent),
                "max_occupancy": round(persistent[0][1], 4),
                "contacts": [
                    {
                        "atom_a": p[0],
                        "atom_b": p[1],
                        "occupancy": round(occ, 4),
                    }
                    for p, occ in persistent[:50]
                ],
            }
        )
    return {
        "chains": chains,
        "contact_distance": contact_dist,
        "min_occupancy": min_occupancy,
        "pairs": summary,
        "n_frames": n_frames,
    }


# ---------------------------------------------------------------------------
# SASA（Shrake-Rupley）
# ---------------------------------------------------------------------------


def sasa_estimate(
    coords: np.ndarray,
    elements: Sequence[str],
    n_points: int = 100,
    probe_radius: float = 1.4,
) -> dict:
    """Shrake-Rupley 数值溶剂可及表面积。

    在每个原子的范德华球面上按黄金螺旋均匀撒点；对每个点统计它是否被其它
    原子的"扩展球"（半径 + 探针半径）覆盖；未被覆盖的点对应的球面元
    即为可及面积。

    Args:
        coords: ``(n_frames, n_atoms, 3)``。
        elements: 元素符号序列。
        n_points: 每个原子的球面采样点数（黄金螺旋）。
        probe_radius: 探针半径（Å，水为 1.4）。

    Returns:
        字典，含总 SASA、逐原子 SASA 与极性/非极性分解。
    """
    coords = np.asarray(coords, dtype=np.float64)
    n_frames, n_atoms, _ = coords.shape
    radii = np.array([_element_radius(e) for e in elements], dtype=np.float64)
    expanded = radii + probe_radius

    # 黄金螺旋球面采样点
    idx = np.arange(n_points, dtype=np.float64) + 0.5
    phi = np.arccos(1.0 - 2.0 * idx / n_points)
    theta = np.pi * (1.0 + 5.0 ** 0.5) * idx
    sphere = np.stack(
        [np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)], axis=1
    )

    per_frame_total = np.zeros(n_frames, dtype=np.float64)
    per_frame_polar = np.zeros(n_frames, dtype=np.float64)
    per_atom_last: np.ndarray | None = None

    for f in range(n_frames):
        pos = coords[f]
        accessible = np.zeros(n_atoms, dtype=np.float64)
        for a in range(n_atoms):
            center = pos[a]
            points = center + radii[a] * sphere
            covered = np.zeros(n_points, dtype=bool)
            delta = points[:, None, :] - pos[None, :, :]
            dist2 = (delta ** 2).sum(axis=2)
            # 自身不算遮挡
            dist2[:, a] = np.inf
            within = dist2 < (expanded[None, :] ** 2)
            covered = within.any(axis=1)
            n_acc = int((~covered).sum())
            area = 4.0 * math.pi * radii[a] ** 2 * (n_acc / n_points)
            accessible[a] = area
        per_frame_total[f] = accessible.sum()
        polar_mask = np.array(
            [str(e).strip().upper() in ("N", "O") for e in elements], dtype=bool
        )
        per_frame_polar[f] = accessible[polar_mask].sum()
        per_atom_last = accessible

    return {
        "per_frame_total": per_frame_total,
        "per_frame_polar": per_frame_polar,
        "per_frame_nonpolar": per_frame_total - per_frame_polar,
        "mean_total": float(per_frame_total.mean()) if n_frames else 0.0,
        "std_total": float(per_frame_total.std()) if n_frames else 0.0,
        "per_atom_last_frame": per_atom_last.tolist() if per_atom_last is not None else [],
        "n_points": n_points,
        "probe_radius": probe_radius,
    }


# ---------------------------------------------------------------------------
# 自检入口
# ---------------------------------------------------------------------------


def _main(argv: list[str] | None = None) -> int:
    """命令行自检：对合成轨迹跑通全部核心算法。"""
    import argparse

    sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "packages"))
    from biopipe_core import get_logger

    logger = get_logger("biopipe.molecular-dynamics.md_analysis")
    parser = argparse.ArgumentParser(
        prog="md_analysis", description="分子动力学核心算法自检"
    )
    parser.add_argument("--frames", type=int, default=20, help="合成轨迹帧数")
    parser.add_argument("--atoms", type=int, default=12, help="合成轨迹原子数")
    args = parser.parse_args(argv)

    rng = np.random.default_rng(42)
    coords = np.cumsum(rng.normal(scale=0.1, size=(args.frames, args.atoms, 3)), axis=1)
    elements = ["C"] * args.atoms

    # Kabsch 自洽性检验：施加已知旋转+平移后应能精确还原
    theta = 0.7
    rot_true = np.array(
        [
            [np.cos(theta), -np.sin(theta), 0.0],
            [np.sin(theta), np.cos(theta), 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    base = rng.normal(size=(60, 3))
    shifted = base @ rot_true.T + np.array([1.0, 2.0, 3.0])
    rot, trans, value = kabsch_align(base, shifted)
    logger.info(
        "Kabsch 自洽检验: RMSD=%.3e, 旋转误差=%.3e, det(R)=%+.1f",
        value,
        float(np.abs(rot - rot_true).max()),
        float(np.linalg.det(rot)),
    )

    rot, trans, value = kabsch_align(coords[-1], coords[0])
    logger.info("Kabsch 叠合 RMSD = %.4f Å（合成轨迹，应接近 0）", value)
    logger.info("RMSF 均值 = %.4f Å", rmsf(coords).mean())
    logger.info("Rg 均值 = %.4f Å", radius_of_gyration(coords).mean())
    pca = pca_analysis(coords)
    logger.info(
        "PCA 前 3 个主成分解释方差 = %s",
        np.round(pca["explained_variance_ratio"], 4).tolist(),
    )
    matrix, _ = rmsd_matrix(coords, max_frames=10)
    cluster = hierarchical_cluster(matrix, n_clusters=2)
    logger.info("层次聚类得到 %d 个簇", cluster["n_clusters"])
    logger.info("SASA 均值 = %.2f Å²", sasa_estimate(coords, elements, n_points=30)["mean_total"])
    return 0


if __name__ == "__main__":
    sys.exit(_main())
