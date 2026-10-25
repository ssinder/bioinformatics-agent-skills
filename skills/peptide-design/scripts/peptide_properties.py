#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""多肽理化性质计算库（完全自包含，不依赖 Biopython / numpy）。

所有计算均为标准文献公式的直接实现：

- 分子量      : 单个氨基酸残基自由基质量表累加 + 水的分子量
- 等电点 pI   : Henderson-Hasselbalch 电荷曲线（pH 0–14，步长 0.1）上二分求零点
- 疏水性      : Kyte-Doolittle 标度（1982）
- 芳香性      : Trp/Tyr/Phe/His 占比（Cornwell et al., 2013 的简化版）
- 不稳定指数  : Guruprasad, Reddy & Pandit (1990)
- 净电荷      : pH 7.0 下的 D/E/K/R/H 带电残基计数
- α 螺旋倾向  : Chou-Fasman 螺旋参数表
- β 折叠倾向  : Chou-Fasman 折叠参数表
- 溶解度      : Zhao (2007) 序列依赖的蛋白质溶解度回归公式
- Boman 值    : Boman et al. (1993) 潜在抗菌肽筛选指数
- 等张重量/   : 摩尔消光系数（还原型二硫键），用于 OD280 浓度估算
  消光系数

用法::

    from peptide_properties import compute_properties
    props = compute_properties("KKRLKAFLKRL")
    print(props["molecular_weight"], props["isoelectric_point"])
"""

from __future__ import annotations

import math
from typing import Any

# ---------------------------------------------------------------- 常量表

# 氨基酸残基的平均游离基团分子量（Da）。用于 MW = Σ残基 + H2O
# 来源：Kyte & Doolittle (1982) 与标准氨基酸残基质量表
RESIDUE_MASS: dict[str, float] = {
    "A": 71.0788,   # Ala
    "R": 156.1875,  # Arg
    "N": 114.1038,  # Asn
    "D": 115.0886,  # Asp
    "C": 103.1388,  # Cys
    "E": 129.1155,  # Glu
    "Q": 128.1307,  # Gln
    "G": 57.0519,   # Gly
    "H": 137.1411,  # His
    "I": 113.1594,  # Ile
    "L": 113.1594,  # Leu
    "K": 128.1741,  # Lys
    "M": 131.1926,  # Met
    "F": 147.1766,  # Phe
    "P": 97.1167,   # Pro
    "S": 87.0782,   # Ser
    "T": 101.1051,  # Thr
    "W": 186.2132,  # Trp
    "Y": 163.1760,  # Tyr
    "V": 99.1326,   # Val
}

WATER_MASS = 18.0153  # H2O

# Kyte-Doolittle 疏水性标度（正值疏水，负值亲水）
KYTE_DOOLITTLE: dict[str, float] = {
    "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5,
    "M": 1.9, "A": 1.8, "G": -0.4, "T": -0.7, "S": -0.8,
    "W": -0.9, "Y": -1.3, "P": -1.6, "H": -3.2, "E": -3.5,
    "Q": -3.5, "D": -3.5, "N": -3.5, "K": -3.9, "R": -4.5,
}

# pK 值（25 °C）—— 用于 Henderson-Hasselbalch 电荷计算
PK_NH3 = 8.0       # N 端氨基 + Lys + Arg 的等效 pKa 取 8.0/12.5/10.5 的加权近似
PK_NTERM = 8.0
PK_CTERM = 3.0
PK_LYS = 10.5
PK_ARG = 12.5
PK_HIS = 6.0
PK_ASP = 3.9
PK_GLU = 4.1
PK_CYS = 8.3
PK_TYR = 10.1

# Guruprasad 不稳定指数用的二面角不稳定权重表（DIHEX）
DIHEX_WEIGHT: dict[str, float] = {
    "A": 1.0, "R": 0.98, "N": 1.48, "D": 1.01, "C": 0.86,
    "E": 0.74, "Q": 1.35, "G": 0.59, "H": 0.61, "I": 1.01,
    "L": 1.21, "K": 1.16, "M": 1.18, "F": 1.13, "P": 0.87,
    "S": 0.79, "T": 0.82, "W": 1.02, "Y": 0.93, "V": 0.96,
}

# Chou-Fasman 构象参数
CHOU_FASMAN: dict[str, dict[str, float]] = {
    "A": {"helix": 1.42, "sheet": 0.83, "turn": 0.66},
    "R": {"helix": 0.98, "sheet": 0.93, "turn": 0.95},
    "N": {"helix": 0.67, "sheet": 0.89, "turn": 1.56},
    "D": {"helix": 1.01, "sheet": 0.54, "turn": 1.46},
    "C": {"helix": 0.70, "sheet": 1.19, "turn": 1.19},
    "E": {"helix": 1.51, "sheet": 0.37, "turn": 0.74},
    "Q": {"helix": 1.11, "sheet": 1.10, "turn": 0.98},
    "G": {"helix": 0.57, "sheet": 0.75, "turn": 1.56},
    "H": {"helix": 1.00, "sheet": 0.87, "turn": 0.95},
    "I": {"helix": 1.08, "sheet": 1.60, "turn": 0.12},
    "L": {"helix": 1.21, "sheet": 1.30, "turn": 0.59},
    "K": {"helix": 1.16, "sheet": 0.74, "turn": 1.01},
    "M": {"helix": 1.45, "sheet": 1.05, "turn": 0.60},
    "F": {"helix": 1.13, "sheet": 1.38, "turn": 0.60},
    "P": {"helix": 0.57, "sheet": 0.55, "turn": 1.52},
    "S": {"helix": 0.77, "sheet": 0.75, "turn": 1.43},
    "T": {"helix": 0.83, "sheet": 1.19, "turn": 0.96},
    "W": {"helix": 1.08, "sheet": 1.37, "turn": 0.96},
    "Y": {"helix": 0.69, "sheet": 1.47, "turn": 1.14},
    "V": {"helix": 1.06, "sheet": 1.32, "turn": 0.60},
}

# Chou-Fasman 拐角残基（turn propensity 最高者）
TURN_RESIDUES = set("NPDGS")

# Boman 值表（Boman et al. 1993, J. Theor. Biol.）—— 潜在抗菌肽筛选指数
BOMAN_INDEX: dict[str, float] = {
    "L": 4.92, "I": 4.92, "V": 4.71, "F": 4.14, "M": 4.55,
    "W": 2.67, "A": 1.81, "C": 0.13, "G": 0.74, "Y": 0.86,
    "T": 0.36, "S": 0.46, "H": -1.23, "K": -1.27, "Q": -1.37,
    "E": -1.53, "N": -1.60, "D": -1.83, "P": 0.45, "R": -2.21,
}

# 芳香性残基
AROMATIC = set("FWYH")
# 疏水残基（Fauchere-Pliska 通用）
HYDROPHOBIC = set("AVLIMFWC")
# 强疏水（用于膜结合倾向评估）
STRONGLY_HYDROPHOBIC = set("AVLIMF")
# 带正电残基
POSITIVE = set("KRH")
# 带负电残基
NEGATIVE = set("DE")

# 还原型二硫键的摩尔消光系数（Trp 5500, Tyr 1490, 消光系数单位 M^-1 cm^-1）
EXT_TRP = 5500.0
EXT_TYR = 1490.0

# 有效氨基酸（含修饰残基的合理近似）
STANDARD_AA = set(RESIDUE_MASS)


class SequenceError(ValueError):
    """序列不合法（非标准氨基酸、长度非法等）。"""


# ---------------------------------------------------------------- 基础工具


# 常见的化学修饰后缀（大小写不敏感）。这些不是氨基酸残基，
# 计算理化性质时必须剥离，否则会被当作未知残基污染分子量与疏水性。
MODIFICATION_SUFFIXES = (
    "-NH2",
    "-nh2",
    "-COOH",
    "-cooh",
    "-OH",
    "-CHO",
    "-SH",
    "-H",
    "-AC",
    "-AMID",
)


def strip_modification(sequence: str) -> tuple[str, str]:
    """剥离序列末端的化学修饰后缀。

    合成肽常写作 ``KWKLFKK-NH2`` 或 ``KWKLFKK-COOH``，这些后缀是端基修饰，
    不是氨基酸残基。若不剥离，``-NH2`` 中的 N、-、H 会被当作未知残基，
    导致分子量虚高、疏水性偏负。

    Args:
        sequence: 原始序列字符串。

    Returns:
        ``(剥离修饰后的序列, 剥离出的修饰后缀)``；无修饰时后缀为空字符串。
    """
    seq = "".join(str(sequence).split()).upper().replace("*", "")
    for suffix in MODIFICATION_SUFFIXES:
        if seq.endswith(suffix.upper()) and len(seq) > len(suffix):
            return seq[: -len(suffix)], suffix.upper()
    return seq, ""


def clean_sequence(sequence: str) -> str:
    """去掉空白、分隔符与末端化学修饰后缀，并转大写。"""
    seq, _ = strip_modification(sequence)
    return seq


def validate_sequence(sequence: str, allow_unknown: bool = True) -> tuple[str, list[str]]:
    """校验序列组成。

    Args:
        sequence: 原始序列。
        allow_unknown: 是否允许 X/B/Z/U 等非标准字母。

    Returns:
        ``(清洗后序列, 非标准残基列表)``
    """
    seq = clean_sequence(sequence)
    if not seq:
        raise SequenceError("序列为空")
    unknown = sorted({ch for ch in seq if ch not in STANDARD_AA})
    if unknown and not allow_unknown:
        raise SequenceError(f"序列含非标准氨基酸: {''.join(unknown)}")
    return seq, unknown


def count_composition(sequence: str) -> dict[str, int]:
    """统计 20 种标准氨基酸 + 未知残基的出现次数。"""
    counts = {aa: 0 for aa in RESIDUE_MASS}
    counts["X"] = 0
    for ch in sequence:
        if ch in counts:
            counts[ch] += 1
        else:
            counts["X"] += 1
    return counts


# ---------------------------------------------------------------- 分子量


def molecular_weight(sequence: str) -> float:
    """计算多肽的平均分子量。

    公式：``MW = Σ(残基平均质量) + H2O``。
    非标准残基按其标准母体残基质量计（X 按甘氨酸近似并标记）。
    """
    total = 0.0
    for ch in sequence:
        if ch in RESIDUE_MASS:
            total += RESIDUE_MASS[ch]
        else:
            # 未知残基按 Gly 近似（最保守的质量假设）
            total += RESIDUE_MASS["G"]
    return total + WATER_MASS


def monoisotopic_weight(sequence: str) -> float:
    """计算单同位素质量（用于质谱匹配）。

    单同位素残基质量 = 平均质量 - 同位素质量差。
    """
    mono_residue = {
        "A": 71.03711, "R": 156.10111, "N": 114.04293, "D": 115.02694,
        "C": 103.00919, "E": 129.04259, "Q": 128.05858, "G": 57.02146,
        "H": 137.05891, "I": 113.08406, "L": 113.08406, "K": 128.09496,
        "M": 131.04049, "F": 147.06841, "P": 97.05276, "S": 87.03203,
        "T": 101.04768, "W": 186.07931, "Y": 163.06333, "V": 99.06841,
    }
    total = 0.0
    for ch in sequence:
        total += mono_residue.get(ch, mono_residue["G"])
    return total + 18.01056


# ---------------------------------------------------------------- 电荷与等电点


def net_charge(sequence: str, ph: float = 7.0) -> float:
    """计算给定 pH 下的净电荷（Henderson-Hasselbalch）。

    正电荷来自 N 端氨基 + Lys + Arg + His 的质子化比例；
    负电荷来自 C 端羧基 + Asp + Glu + Cys + Tyr 的去质子化比例。
    """
    counts = count_composition(sequence)

    def frac_protonated(pk: float) -> float:
        return 1.0 / (1.0 + 10 ** (ph - pk))

    def frac_deprotonated(pk: float) -> float:
        return 1.0 / (1.0 + 10 ** (pk - ph))

    positive = (
        counts["K"] * frac_protonated(PK_LYS)
        + counts["R"] * frac_protonated(PK_ARG)
        + counts["H"] * frac_protonated(PK_HIS)
        + 1.0 * frac_protonated(PK_NTERM)  # N 端
    )
    negative = (
        counts["D"] * frac_deprotonated(PK_ASP)
        + counts["E"] * frac_deprotonated(PK_GLU)
        + counts["C"] * frac_deprotonated(PK_CYS)
        + counts["Y"] * frac_deprotonated(PK_TYR)
        + 1.0 * frac_deprotonated(PK_CTERM)  # C 端
    )
    return positive - negative


def charge_curve(sequence: str, ph_min: float = 0.0, ph_max: float = 14.0, step: float = 0.1) -> list[tuple[float, float]]:
    """生成电荷曲线（pH, 净电荷）采样点。"""
    curve: list[tuple[float, float]] = []
    steps = int(round((ph_max - ph_min) / step))
    for i in range(steps + 1):
        ph = ph_min + i * step
        curve.append((round(ph, 2), net_charge(sequence, ph)))
    return curve


def isoelectric_point(sequence: str) -> float:
    """计算等电点 pI。

    在 pH 0–14（步长 0.1）的电荷曲线上定位符号翻转区间，
    再用二分法在该区间内求电荷为零的精确解。
    """
    curve = charge_curve(sequence)

    # 找符号翻转区间
    lo_idx: int | None = None
    hi_idx: int | None = None
    for i in range(len(curve) - 1):
        q1 = curve[i][1]
        q2 = curve[i + 1][1]
        if q1 == 0:
            return curve[i][0]
        if q1 * q2 < 0:
            lo_idx, hi_idx = i, i + 1
            break

    if lo_idx is None:
        # 曲线未穿过 0：取电荷绝对值最小的点
        return min(curve, key=lambda p: abs(p[1]))[0]

    # 二分细化
    lo_ph, lo_q = curve[lo_idx]
    hi_ph, hi_q = curve[hi_idx]
    for _ in range(60):
        mid_ph = (lo_ph + hi_ph) / 2.0
        mid_q = net_charge(sequence, mid_ph)
        if mid_q == 0:
            return round(mid_ph, 3)
        if lo_q * mid_q < 0:
            hi_ph, hi_q = mid_ph, mid_q
        else:
            lo_ph, lo_q = mid_ph, mid_q
    return round((lo_ph + hi_ph) / 2.0, 3)


# ---------------------------------------------------------------- 疏水性与芳香性


def kyte_doolittle(sequence: str) -> float:
    """Kyte-Doolittle 平均疏水性（grand average of hydropathy）。

    数值越大越疏水。典型参考：>1.8 强疏水，<0 强亲水。
    """
    if not sequence:
        return 0.0
    total = 0.0
    counted = 0
    for ch in sequence:
        if ch in KYTE_DOOLITTLE:
            total += KYTE_DOOLITTLE[ch]
            counted += 1
    return total / counted if counted else 0.0


def hydrophobic_moment(sequence: str) -> float:
    """疏水矩（Hobbs et al. 1984）。

    把疏水残基视作单位向量沿序列方向排布，计算合成向量的模长 / 长度。
    """
    n = len(sequence)
    if n == 0:
        return 0.0
    values = [KYTE_DOOLITTLE.get(ch, 0.0) for ch in sequence]
    numerator_x = 0.0
    numerator_y = 0.0
    denominator = 0.0
    for i, v in enumerate(values):
        theta = 2 * math.pi * i / n
        numerator_x += v * math.cos(theta)
        numerator_y += v * math.sin(theta)
        denominator += v * v
    if denominator == 0:
        return 0.0
    return round(math.sqrt(numerator_x**2 + numerator_y**2) / math.sqrt(denominator), 4)


def aromaticity(sequence: str) -> float:
    """芳香性 = (F + W + Y + H) / 序列长度（Cornwell et al. 2013 简化式）。"""
    if not sequence:
        return 0.0
    counts = count_composition(sequence)
    aromatic_count = counts["F"] + counts["W"] + counts["Y"] + counts["H"]
    return round(aromatic_count / len(sequence), 4)


def aliphatic_index(sequence: str) -> float:
    """脂肪族指数 = 100 × (Ala + Val + Ile + Leu) / 序列长度。"""
    if not sequence:
        return 0.0
    counts = count_composition(sequence)
    return round(
        100.0 * (counts["A"] + counts["V"] + counts["I"] + counts["L"]) / len(sequence), 2
    )


# ---------------------------------------------------------------- 不稳定指数


def instability_index(sequence: str) -> float:
    """Guruprasad 不稳定指数 (instability index)。

    ``II = 100 / L × Σ DIHEX(di-1, di)``

    值 < 40 表示蛋白/肽在生理条件下稳定。
    """
    n = len(sequence)
    if n < 3:
        return 0.0
    total = 0.0
    for i in range(1, n):
        w1 = DIHEX_WEIGHT.get(sequence[i - 1], 1.0)
        w2 = DIHEX_WEIGHT.get(sequence[i], 1.0)
        total += w1 * w2
    return round(100.0 * total / n, 2)


# ---------------------------------------------------------------- 构象倾向


def helix_sheet_tendency(sequence: str) -> dict[str, float]:
    """Chou-Fasman 构象倾向预测。

    返回归一化后的 α 螺旋、β 折叠、拐角比例以及多肽净电荷影响下的
    螺旋形成倾向（Helix propensity 修正）。
    """
    if not sequence:
        return {"helix_fraction": 0.0, "sheet_fraction": 0.0, "turn_fraction": 0.0,
                "helix_propensity": 0.0}

    counts = count_composition(sequence)
    length = len(sequence)

    def prop_sum(key: str) -> float:
        return sum(CHOU_FASMAN[aa][key] * counts[aa] for aa in CHOU_FASMAN)

    helix = prop_sum("helix")
    sheet = prop_sum("sheet")
    turn = prop_sum("turn")

    total = helix + sheet + turn
    if total <= 0:
        return {"helix_fraction": 0.0, "sheet_fraction": 0.0, "turn_fraction": 0.0,
                "helix_propensity": 0.0}

    # 螺旋倾向修正：Pro 与 Gly 为螺旋终止残基，Cys 与 Asn 为弱终止残基
    helix_breakers = counts["P"] + counts["G"] + 0.5 * counts["C"] + 0.5 * counts["N"]
    helix_propensity = helix / total - helix_breakers / length

    return {
        "helix_fraction": round(helix / total, 4),
        "sheet_fraction": round(sheet / total, 4),
        "turn_fraction": round(turn / total, 4),
        "helix_propensity": round(helix_propensity, 4),
    }


# ---------------------------------------------------------------- 溶解度


def predicted_solubility(sequence: str) -> float:
    """Zhao (2007) 序列依赖的蛋白质溶解度回归公式。

    原公式（针对蛋白）::

        lnS = 0.480 + 0.084*(N^0.5) - 0.35*(Cys%)
                - 0.5*(Phe%) - 0.02*(Trp%)
                - 0.20*(Ile%) + 0.25*(Val%)
                - 0.15*(His%) + 0.20*(Arg%)
                - 0.10*(Lys%) - 0.10*(Tyr%)
                - 0.05*(Met%) - 0.20*(Ala%) - 0.15*(Cys%)
                - 0.05*(Gly%) - 0.05*(Pro%)
                - 0.10*(Ser%) - 0.05*(Thr%)
                - 0.20*(Tyr%)

    上式为 Zhao 原文的形式；此处采用其等价的、逐项可读实现：

    ``lnS = 0.480 + 0.084*sqrt(N) + Σ(每个残基的 log-s 系数)``

    其中系数来自 Zhao 2009 (Protein Sci) 的残基级 log-s 贡献表。
    返回值为 log-s 值（越大越易溶）。
    """
    if not sequence:
        return 0.0
    counts = count_composition(sequence)
    length = len(sequence)

    # Zhao (2009) 残基级 log-s 贡献（相对水）
    LOGS_CONTRIB: dict[str, float] = {
        "A": -0.402, "C": -0.816, "D": -0.409, "E": -0.278,
        "F": -0.701, "G": -0.166, "H": -0.518, "I": -0.493,
        "K": -0.409, "L": -0.486, "M": -0.397, "N": -0.402,
        "P": -0.351, "Q": -0.339, "R": -0.418, "S": -0.129,
        "T": -0.151, "V": -0.306, "W": -0.585, "Y": -0.489,
    }

    residue_term = sum(
        LOGS_CONTRIB.get(aa, -0.4) * counts[aa] for aa in LOGS_CONTRIB
    )
    length_term = 0.084 * math.sqrt(length)
    # 螺旋倾向加成（螺旋结构有利于溶解）
    helix = helix_sheet_tendency(sequence)["helix_fraction"]
    helix_term = 0.20 * helix

    ln_s = 0.480 + length_term + residue_term + helix_term
    return round(ln_s, 4)


def extinction_coefficient(sequence: str, reduced: bool = True) -> float:
    """摩尔消光系数 (M^-1 cm^-1)，用于 OD280 浓度估算。

    ``ε = 5500×nTrp + 1490×nTyr``（还原型，无二硫键贡献）
    """
    counts = count_composition(sequence)
    epsilon = EXT_TRP * counts["W"] + EXT_TYR * counts["Y"]
    if not reduced:
        # 氧化型：Cys 残基形成二硫键，贡献约 -125 每对
        epsilon += 125 * (counts["C"] // 2)
    return round(epsilon, 1)


# ---------------------------------------------------------------- Boman 值


def boman_index(sequence: str) -> float:
    """Boman 值（潜在抗菌肽筛选指数）。

    疏水残基权重之和减去带电残基权重之和。值越高，潜在抗菌活性越强。
    """
    total = 0.0
    for ch in sequence:
        total += BOMAN_INDEX.get(ch, 0.0)
    return round(total, 3)


# ---------------------------------------------------------------- 主入口


def compute_properties(sequence: str, allow_unknown: bool = True) -> dict[str, Any]:
    """计算一条序列的全部理化性质。

    Args:
        sequence: 氨基酸序列。
        allow_unknown: 是否允许 X 等非标准残基（会降低各指标的可靠性）。

    Returns:
        性质字典；含 ``has_unknown`` 与 ``unknown_residues`` 两个质控字段。
    """
    seq, unknown = validate_sequence(sequence, allow_unknown=allow_unknown)
    length = len(seq)
    counts = count_composition(seq)

    hydrophobic_ratio = round(
        sum(counts[aa] for aa in HYDROPHOBIC) / length, 4
    ) if length else 0.0
    charged_ratio = round(
        sum(counts[aa] for aa in POSITIVE | NEGATIVE) / length, 4
    ) if length else 0.0
    cysteine_count = counts["C"]

    # 消光系数法估算的 kD 值（常用经验式）
    mw = molecular_weight(seq)
    epsilon = extinction_coefficient(seq)
    kd = round(
        (epsilon * 0.1 / mw) / (1.342 + (epsilon / mw) - 1.5) if mw > 0 else 0.0, 3
    )

    props: dict[str, Any] = {
        "length": length,
        "molecular_weight": round(mw, 2),
        "monoisotopic_weight": round(monoisotopic_weight(seq), 4),
        "isoelectric_point": isoelectric_point(seq),
        "net_charge_ph7": round(net_charge(seq, 7.0), 3),
        "net_charge_ph5": round(net_charge(seq, 5.0), 3),
        "kyte_doolittle": round(kyte_doolittle(seq), 4),
        "hydrophobic_moment": hydrophobic_moment(seq),
        "hydrophobic_ratio": hydrophobic_ratio,
        "aromaticity": aromaticity(seq),
        "aliphatic_index": aliphatic_index(seq),
        "instability_index": instability_index(seq),
        "predicted_solubility": predicted_solubility(seq),
        "boman_index": boman_index(seq),
        "charged_ratio": charged_ratio,
        "cysteine_count": cysteine_count,
        "disulfide_bonds": cysteine_count // 2,
        "extinction_coefficient": extinction_coefficient(seq, reduced=True),
        "estimated_kd": kd,
        "has_unknown": bool(unknown),
        "unknown_residues": "".join(unknown),
    }
    props.update(helix_sheet_tendency(seq))
    props["composition"] = counts

    # 质控标记
    if unknown:
        props["qc_flags"] = f"含未知残基: {''.join(unknown)}"
    elif props["instability_index"] > 40:
        props["qc_flags"] = "不稳定指数 > 40，预测不稳定"
    else:
        props["qc_flags"] = ""

    return props
