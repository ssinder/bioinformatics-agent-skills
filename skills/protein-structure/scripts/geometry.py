#!/usr/bin/env python3
"""结构几何量计算（纯 Python，仅依赖标准库）。

提供：

- 回旋半径（radius of gyration, Rg）及其分量分解
- 陀螺张量（gyration tensor）主轴分解与形状描述符
- 瞬时偶极矩（由原子部分电荷估算）
- 端到端距离、包围盒、惯性椭球形状分类

所有计算都以质量加权的质心为原点，符合高分子物理的常规定义。

Usage:
    from geometry import radius_of_gyration, gyration_tensor, shape_descriptors
"""

from __future__ import annotations

import math
import sys
from collections.abc import Iterable, Sequence
from pathlib import Path
from typing import Optional

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "packages"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from biopipe_core import get_logger  # noqa: E402
from pdb_parser import parse_pdb_file  # noqa: E402

logger = get_logger("biopipe.protein-structure")

__all__ = [
    "radius_of_gyration",
    "gyration_tensor",
    "shape_descriptors",
    "instantaneous_dipole",
    "end_to_end_distance",
    "bounding_box",
]

# 元素 -> 原子质量（u）
ATOMIC_MASS = {
    "H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999, "S": 32.06,
    "P": 30.974, "F": 18.998, "CL": 35.45, "BR": 79.904, "I": 126.904,
    "MG": 24.305, "ZN": 65.38, "CA": 40.078, "NA": 22.990, "K": 39.098,
    "FE": 55.845, "MN": 54.938, "CU": 63.546, "SE": 78.971, "PBE": 207.2,
}

# 元素 -> 近似部分电荷（e），用于瞬时偶极估算（简化模型）
ELEMENT_PARTIAL_CHARGE = {
    "H": 0.31, "C": 0.05, "N": -0.30, "O": -0.42, "S": -0.20,
    "P": 0.35, "F": -0.32, "CL": -0.35, "BR": -0.30, "I": -0.25,
    "MG": 1.60, "ZN": 1.30, "CA": 1.20, "NA": 0.90, "K": 0.85,
    "FE": 1.80, "MN": 1.50, "CU": 1.20, "SE": -0.10,
}

Coord = Sequence[float]


def _mass_of(atom) -> float:
    """取原子的质量（u），未知元素按碳处理。"""
    return ATOMIC_MASS.get(atom.element.upper(), 12.011)


def _mass_center(atoms: Sequence) -> tuple[Coord, float]:
    """返回质量中心与总质量。"""
    total = sum(_mass_of(a) for a in atoms)
    if total <= 0:
        return (0.0, 0.0, 0.0), 0.0
    cx = sum(a.x * _mass_of(a) for a in atoms) / total
    cy = sum(a.y * _mass_of(a) for a in atoms) / total
    cz = sum(a.z * _mass_of(a) for a in atoms) / total
    return (cx, cy, cz), total


def radius_of_gyration(atoms: Sequence) -> dict:
    """计算回旋半径及分量。

    公式（质量加权）::

        Rg = sqrt( sum_i m_i * |r_i - R_cm|^2 / sum_i m_i )

    同时给出各坐标轴上的分量 Rg_x / Rg_y / Rg_z（以质心为原点），
    以及以 Cα 子集计算的结果（用于与文献中的 Rg 对齐）。

    Args:
        atoms: 原子序列（需具备 ``x``/``y``/``z``/``element`` 属性）。

    Returns:
        包含 ``rg``（Å）、``rg_components``、``n_atoms``、``total_mass``、
        ``rg_ca`` 的字典；无原子时所有数值为 0。
    """
    if not atoms:
        return {
            "rg": 0.0,
            "rg_components": {"x": 0.0, "y": 0.0, "z": 0.0},
            "n_atoms": 0,
            "total_mass": 0.0,
            "rg_ca": None,
            "center_of_mass": (0.0, 0.0, 0.0),
        }

    com, total_mass = _mass_center(atoms)
    sx = sy = sz = 0.0
    for atom in atoms:
        m = _mass_of(atom)
        dx, dy, dz = atom.x - com[0], atom.y - com[1], atom.z - com[2]
        sx += m * dx * dx
        sy += m * dy * dy
        sz += m * dz * dz

    rg = math.sqrt((sx + sy + sz) / total_mass) if total_mass else 0.0
    comp = {
        "x": math.sqrt(sx / total_mass) if total_mass else 0.0,
        "y": math.sqrt(sy / total_mass) if total_mass else 0.0,
        "z": math.sqrt(sz / total_mass) if total_mass else 0.0,
    }

    ca_atoms = [a for a in atoms if a.name.upper() == "CA"]
    rg_ca = None
    if len(ca_atoms) >= 2:
        ca_com, ca_mass = _mass_center(ca_atoms)
        acc = 0.0
        for atom in ca_atoms:
            m = _mass_of(atom)
            acc += m * (
                (atom.x - ca_com[0]) ** 2
                + (atom.y - ca_com[1]) ** 2
                + (atom.z - ca_com[2]) ** 2
            )
        rg_ca = math.sqrt(acc / ca_mass) if ca_mass else None

    return {
        "rg": rg,
        "rg_components": comp,
        "n_atoms": len(atoms),
        "total_mass": total_mass,
        "rg_ca": rg_ca,
        "center_of_mass": com,
    }


def gyration_tensor(atoms: Sequence) -> dict:
    """计算质量加权的陀螺张量并做主轴分解。

    陀螺张量 S = sum_i m_i (r_i - R_cm)(r_i - R_cm)^T，
    三个特征值对应主轴上的回旋半径平方（除以总质量）。
    特征分解采用 Jacobi 旋转法（纯 Python，无 numpy 依赖）。

    Returns:
        ``{"eigenvalues": [a, b, c], "axes": [[u1],[u2],[u3]], "rg": ...}``
        特征值已降序排列。
    """
    if not atoms:
        return {"eigenvalues": [0.0, 0.0, 0.0], "axes": [[1, 0, 0], [0, 1, 0], [0, 0, 1]], "rg": 0.0}

    com, total_mass = _mass_center(atoms)
    sxx = sxy = sxz = syy = syz = szz = 0.0
    for atom in atoms:
        m = _mass_of(atom)
        dx, dy, dz = atom.x - com[0], atom.y - com[1], atom.z - com[2]
        sxx += m * dx * dx
        sxy += m * dx * dy
        sxz += m * dx * dz
        syy += m * dy * dy
        syz += m * dy * dz
        szz += m * dz * dz

    matrix = [
        [sxx / total_mass, sxy / total_mass, sxz / total_mass],
        [sxy / total_mass, syy / total_mass, syz / total_mass],
        [sxz / total_mass, syz / total_mass, szz / total_mass],
    ]
    eigenvalues, eigenvectors = _jacobi_eigen(matrix)
    total = sum(eigenvalues)
    return {
        "eigenvalues": eigenvalues,
        "axes": eigenvectors,
        "rg": math.sqrt(total) if total > 0 else 0.0,
    }


def _jacobi_eigen(matrix: list[list[float]], max_sweeps: int = 100, tol: float = 1e-12):
    """对称 3x3 矩阵的 Jacobi 特征分解。

    Args:
        matrix: 3x3 对称矩阵（按行给出）。
        max_sweeps: 最大旋转次数。
        tol: 收敛阈值。

    Returns:
        ``(eigenvalues, eigenvectors)``，特征值降序，``eigenvectors[k]``
        是第 k 个特征值对应的单位主轴向量。
    """
    a = [row[:] for row in matrix]
    v = [[1.0 if i == j else 0.0 for j in range(3)] for i in range(3)]

    for _ in range(max_sweeps):
        off = a[0][1] ** 2 + a[0][2] ** 2 + a[1][2] ** 2
        if off < tol:
            break
        for p, q in ((0, 1), (0, 2), (1, 2)):
            if abs(a[p][q]) < 1e-18:
                continue
            theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q])
            t = (1.0 if theta >= 0 else -1.0) / (abs(theta) + math.sqrt(theta * theta + 1.0))
            c = 1.0 / math.sqrt(t * t + 1.0)
            s = t * c

            for k in range(3):
                akp, akq = a[k][p], a[k][q]
                a[k][p] = c * akp - s * akq
                a[k][q] = s * akp + c * akq
            for k in range(3):
                apk, aqk = a[p][k], a[q][k]
                a[p][k] = c * apk - s * aqk
                a[q][k] = s * apk + c * aqk
            for k in range(3):
                vkp, vkq = v[k][p], v[k][q]
                v[k][p] = c * vkp - s * vkq
                v[k][q] = s * vkp + c * vkq

    eigenvalues = [a[i][i] for i in range(3)]
    order = sorted(range(3), key=lambda i: eigenvalues[i], reverse=True)
    eigenvalues = [eigenvalues[i] for i in order]
    eigenvectors = [[v[i][j] for i in range(3)] for j in order]
    return eigenvalues, eigenvectors


def shape_descriptors(atoms: Sequence) -> dict:
    """由陀螺张量特征值导出形状描述符。

    特征值按 **降序** 记为 ``λ1 ≥ λ2 ≥ λ3``（分别对应最长、中等、最短主轴）：

    - ``asphericity`` b = λ1 - λ2（棒状程度；球形为 0）
    - ``acylindricity`` c = λ2 - λ3
    - ``relative_shape_anisotropy`` kappa = 0.5 * (b² + c²) / Rg⁴（0 = 球形，1 = 棒形）
    - ``shape_class``：rod / coil / sphere 的经验判据
    """
    gt = gyration_tensor(atoms)
    lam = gt["eigenvalues"]
    rg2 = sum(lam)
    l1, l2, l3 = lam
    b = l1 - l2
    c = l2 - l3
    kappa = 0.5 * (b * b + c * c) / (rg2 * rg2) if rg2 > 0 else 0.0

    if rg2 <= 0:
        shape = "undefined"
    elif b / rg2 < 0.01:
        shape = "sphere"
    elif c / rg2 < 0.05:
        shape = "rod"
    else:
        shape = "coil"

    return {
        "eigenvalues": lam,
        "rg": gt["rg"],
        "asphericity": b,
        "acylindricity": c,
        "relative_shape_anisotropy": kappa,
        "shape_class": shape,
    }


def instantaneous_dipole(atoms: Sequence) -> dict:
    """估算瞬时偶极矩（Debye）。

    使用元素级简化部分电荷（见 ``ELEMENT_PARTIAL_CHARGE``），
    对所有原子求和：``mu = sum_i q_i * r_i``（以质心为原点则电荷和为 0）。

    该值是数量级的粗略估计，仅用于不同体系/构象之间的相对比较，
    不应与 QM 计算的偶极矩直接对照。
    """
    if not atoms:
        return {"dipole": (0.0, 0.0, 0.0), "magnitude": 0.0, "n_charged": 0}

    com, _ = _mass_center(atoms)
    mx = my = mz = 0.0
    charged = 0
    for atom in atoms:
        q = ELEMENT_PARTIAL_CHARGE.get(atom.element.upper())
        if q is None:
            continue
        charged += 1
        mx += q * (atom.x - com[0])
        my += q * (atom.y - com[1])
        mz += q * (atom.z - com[2])

    magnitude = math.sqrt(mx * mx + my * my + mz * mz)
    return {
        "dipole": (mx, my, mz),
        "magnitude": magnitude,
        "n_charged": charged,
    }


def end_to_end_distance(atoms: Sequence) -> float:
    """首尾原子间距离（Å）。"""
    if len(atoms) < 2:
        return 0.0
    return math.dist(atoms[0].coord, atoms[-1].coord)


def bounding_box(atoms: Iterable) -> dict:
    """包围盒（Å）。"""
    xs: list[float] = []
    ys: list[float] = []
    zs: list[float] = []
    for atom in atoms:
        xs.append(atom.x)
        ys.append(atom.y)
        zs.append(atom.z)
    if not xs:
        return {"dx": 0.0, "dy": 0.0, "dz": 0.0, "volume": 0.0}
    dx = max(xs) - min(xs)
    dy = max(ys) - min(ys)
    dz = max(zs) - min(zs)
    return {"dx": dx, "dy": dy, "dz": dz, "volume": dx * dy * dz}


def _main(argv: Optional[list[str]] = None) -> int:
    """命令行入口：解析 PDB 并打印几何量。"""
    import argparse

    parser = argparse.ArgumentParser(
        prog="geometry", description="结构几何量计算（命令行自检工具）"
    )
    parser.add_argument("input", help="PDB 文件路径")
    args = parser.parse_args(argv)

    structure = parse_pdb_file(args.input)
    atoms = structure.atoms()
    logger.info("Rg: %.3f Å", radius_of_gyration(atoms)["rg"])
    shape = shape_descriptors(atoms)
    logger.info("形状: %s (κ=%.4f)", shape["shape_class"], shape["relative_shape_anisotropy"])
    logger.info("偶极矩: %.3f D", instantaneous_dipole(atoms)["magnitude"])
    return 0


if __name__ == "__main__":
    sys.exit(_main())
