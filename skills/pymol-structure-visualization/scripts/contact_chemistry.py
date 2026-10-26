#!/usr/bin/env python3
"""残基化学分类表 —— 链间接触类型判定的化学依据。

这个模块是整个技能的知识核心：**接触类型不是靠几何距离猜的，而是靠
参与接触的两侧残基的化学性质组合判定的。**

判定链条
--------
1. **原子层面**：氢键要求一侧是供体原子（N/O 上带 H）、另一侧是受体原子
   （N/O 上的孤对电子）。只靠距离无法区分氢键与任意极性接触。
2. **残基层面**：盐桥要求正电残基（LYS/ARG/HIS）与负电残基（ASP/GLU）；
   π-π 堆积要求两侧都是芳香族；阳离子-π 要求正电与芳香相遇。
3. **距离层面**：原子间距小于该类型的作用距离上限才算接触。

每种类型有各自的距离阈值 —— 氢键 3.5 Å、盐桥 4.0 Å、疏水 4.5 Å、
π-π 5.0 Å、阳离子-π 6.0 Å。

来源
----
从老仓库 ``pymol_show/scripts/interaction count.py`` 移植。该脚本把这些
分类表硬编码在函数调用里，且全部依赖 PyMOL 的 ``cmd`` 对象（必须先
``pymol.finish_launching()``），**无法 headless 运行**。

本模块把分类表抽成纯数据 + 纯函数，**零依赖**，可以在任何环境跑。
PyMOL 降级为可选的渲染层。

修掉的问题
----------
1. **分类表不一致** —— 老仓库里 ``donors`` 列表含 ASP/GLU（``hbond count.py``）
   而 ``interaction count.py`` 的 ``donors`` 不含；``acceptors`` 在两处也不同。
   同一套相互作用在不同脚本里被判成不同类型。现在统一到一份表。
2. **ASP/GLU 被列为供体是错的** —— 羧基氧没有可离散的氢，不构成供体。
   ASP/GLU 是纯受体。保留在 ``DONOR_AMBIGUOUS`` 里仅作对照说明。
3. **缺少 HIS 的电荷态区分** —— HIS 在 pH 7 约 10% 质子化，既可作供体
   也可作受体，但**不应**计入正电残基（老仓库的 ``positive`` 列表含 HIS，
   会虚增盐桥数）。这里把 HIS 归入 ``POSITIVE_PARTIAL``，默认不参与
   盐桥判定，除非显式开启。
4. **CYS 未计入任何列表** —— 疏水列表漏了 CYS。

本模块只依赖 Python 标准库。
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

__all__ = [
    "INTERACTION_TYPES",
    "DONORS",
    "ACCEPTORS",
    "POSITIVE",
    "POSITIVE_PARTIAL",
    "NEGATIVE",
    "HYDROPHOBIC",
    "AROMATIC",
    "IONIZABLE",
    "STANDARD_RESIDUES",
    "LIGAND_CONTACT",
    "ContactRule",
    "classify_pair",
    "classify_pair_multi",
    "is_donor_residue",
    "is_acceptor_residue",
    "is_positive_residue",
    "is_negative_residue",
    "is_hydrophobic_residue",
    "is_aromatic_residue",
    "is_standard_residue",
    "residue_kind",
    "distance_cutoff",
    "describe_rules",
]


# ---------------------------------------------------------------- 残基化学分类

#: **氢键供体**残基 —— 侧链能提供 N-H / O-H 键。
#:
#: 注意 **ASP / GLU 不在此列**：羧基（-COOH）上的氧虽然在质子化状态下
#: 带氢，但它是**弱供体**（pKa ≈ 4），在生理 pH 下几乎完全去质子化，
#: 作为氢键供体的贡献可以忽略。把它们当供体会虚增氢键数。
DONORS = frozenset({
    "ARG", "ASN", "GLN", "HIS", "LYS", "SER", "THR", "TRP", "TYR",
})

#: **氢键受体**残基 —— 侧链有带孤对电子的 N / O。
ACCEPTORS = frozenset({
    "ASN", "GLN", "HIS", "SER", "THR", "TYR",
    "ASP", "GLU",
})

#: **带正电**残基（生理 pH 下几乎完全质子化）—— 盐桥的供体侧。
POSITIVE = frozenset({"LYS", "ARG"})

#: **部分带正电**残基 —— HIS 的咪唑基 pKa ≈ 6.0，生理 pH 下约 10% 质子化。
#:
#: 老仓库把 HIS 直接算进正电残基，会**虚增盐桥数量**。这里单列出来，
#: 默认**不**参与盐桥判定（:func:`classify_pair` 的 ``include_histidine``
#: 参数控制）。
POSITIVE_PARTIAL = frozenset({"HIS"})

#: **带负电**残基 —— 盐桥的受体侧。
NEGATIVE = frozenset({"ASP", "GLU"})

#: **疏水**残基 —— 疏水接触的两侧。
HYDROPHOBIC = frozenset({
    "ALA", "CYS", "ILE", "LEU", "MET", "PHE", "PRO", "TRP", "VAL",
})

#: **芳香族**残基 —— π-π 堆积与阳离子-π 的参与者。
AROMATIC = frozenset({"PHE", "TRP", "TYR", "HIS"})

#: 可电离残基（带电 + 部分带电），用于"这个残基的电荷状态"判定。
IONIZABLE = POSITIVE | NEGATIVE | POSITIVE_PARTIAL

#: 标准氨基酸三字母码全集（用于校验输入）
STANDARD_RESIDUES = frozenset({
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
})

#: 非标准残基（配体、水、离子等 HETATM 残基）的接触类型标记。
#:
#: 配体接触**不参与**残基级化学分类 —— 我们不知道配体是芳香还是极性、
#: 是否带电，套用氨基酸的分类表会得出无意义的结论（比如把苯环配体
#: 当成"疏水残基"，把带电配体当成"负电残基"）。
#:
#: 但配体接触**也不该被静默丢弃** —— 配体结合位点往往就是研究热点。
#: 所以这类接触会被单独标记为 ``Ligand_contact``，距离照常报告。
LIGAND_CONTACT = "Ligand_contact"


# ---------------------------------------------------------------- 接触规则

@dataclass(frozen=True)
class ContactRule:
    """一种相互作用类型的判定规则。

    Attributes:
        name: 类型标识（CSV 里的 ``Interaction_Type`` 值）。
        label: 中文显示名。
        max_distance: 原子间距上限（Å）。
        description: 判定条件的人类可读描述。
    """

    name: str
    label: str
    max_distance: float
    description: str


#: 全部接触类型。顺序即报告与 CSV 里的排序顺序。
INTERACTION_TYPES: tuple[ContactRule, ...] = (
    ContactRule(
        name="Hydrogen_bond",
        label="氢键",
        max_distance=3.5,
        description="一侧为氢键供体残基、另一侧为受体残基，原子间距 ≤ 3.5 Å",
    ),
    ContactRule(
        name="Salt_bridge",
        label="盐桥",
        max_distance=4.0,
        description="正电残基（LYS/ARG）与负电残基（ASP/GLU）之间，原子间距 ≤ 4.0 Å",
    ),
    ContactRule(
        name="Hydrophobic",
        label="疏水接触",
        max_distance=4.5,
        description="两侧均为疏水残基，原子间距 ≤ 4.5 Å",
    ),
    ContactRule(
        name="Pi_Pi",
        label="π-π 堆积",
        max_distance=5.0,
        description="两侧均为芳香族残基，原子间距 ≤ 5.0 Å",
    ),
    ContactRule(
        name="Cation_Pi",
        label="阳离子-π",
        max_distance=6.0,
        description="正电残基与芳香族残基之间，原子间距 ≤ 6.0 Å",
    ),
)


# ---------------------------------------------------------------- 判定函数


def _norm(res_name: str) -> str:
    """残基名归一化：大写 + 去掉首尾空白。"""
    return (res_name or "").strip().upper()


def is_donor_residue(res_name: str, include_ambiguous: bool = False) -> bool:
    """该残基能否作为氢键供体。

    参数
    ----
    res_name
        残基三字母码。
    include_ambiguous
        是否把 ASP/GLU 也算作（弱）供体。默认 **False** ——
        羧基在生理 pH 下去质子化，不是有效供体。老仓库把 ASP/GLU
        列进供体表，会虚增氢键数量。
    """
    name = _norm(res_name)
    if name in DONORS:
        return True
    if include_ambiguous and name in ("ASP", "GLU"):
        return True
    return False


def is_acceptor_residue(res_name: str) -> bool:
    """该残基能否作为氢键受体。"""
    return _norm(res_name) in ACCEPTORS


def is_positive_residue(res_name: str, include_histidine: bool = False) -> bool:
    """该残基是否带正电。

    参数
    ----
    res_name
        残基三字母码。
    include_histidine
        是否把 HIS 算作正电。默认 **False** —— HIS 咪唑基 pKa ≈ 6.0，
        生理 pH 下仅约 10% 质子化。老仓库直接把它算进正电列表，
        会虚增盐桥数量。
    """
    name = _norm(res_name)
    if name in POSITIVE:
        return True
    if include_histidine and name in POSITIVE_PARTIAL:
        return True
    return False


def is_negative_residue(res_name: str) -> bool:
    """该残基是否带负电。"""
    return _norm(res_name) in NEGATIVE


def is_hydrophobic_residue(res_name: str) -> bool:
    """该残基是否疏水。"""
    return _norm(res_name) in HYDROPHOBIC


def is_aromatic_residue(res_name: str) -> bool:
    """该残基是否芳香族。"""
    return _norm(res_name) in AROMATIC


def is_standard_residue(res_name: str) -> bool:
    """该残基是否标准氨基酸（20 种）。

    非标准残基（配体、水、离子、修饰残基）无法套用氨基酸的化学分类表。
    """
    return _norm(res_name) in STANDARD_RESIDUES


def residue_kind(res_name: str) -> str:
    """残基的类别。

    返回
    ----
    ``"standard"``
        标准氨基酸 —— 参与完整的化学分类判定。
    ``"ligand"``
        非标准残基（配体等）—— 只标记接触，不做化学分类。
    ``"water"``
        水分子。
    """
    name = _norm(res_name)
    if name in STANDARD_RESIDUES:
        return "standard"
    if name in ("HOH", "WAT", "DOD", "H2O", "TIP"):
        return "water"
    return "ligand"


def classify_pair(
    res_a: str,
    res_b: str,
    distance: float,
    include_histidine: bool = False,
    include_ambiguous_donors: bool = False,
) -> Optional[str]:
    """判定两个残基之间属于哪种相互作用。

    参数
    ----
    res_a
        第一个残基的三字母码。
    res_b
        第二个残基的三字母码。
    distance
        参与接触的两个原子间的最小距离（Å）。
    include_histidine
        是否把 HIS 当作正电残基（影响盐桥与阳离子-π）。
    include_ambiguous_donors
        是否把 ASP/GLU 当作弱供体（影响氢键）。

    返回
    ----
    相互作用类型名（``INTERACTION_TYPES`` 里的 ``name``）；
    不满足任何类型时返回 ``None``。

    说明
    ----
    判定**按优先级顺序**进行，一个残基对只归入**第一个**匹配的类型。
    优先级反映的是"特异性"：氢键比疏水接触更具体，
    盐桥比一般的极性接触更具体。

    但有两个例外，它们是**叠加**而非互斥的：

    - **π-π 堆积**：两侧都是芳香族时，即使同时也是疏水的，
      也应报 π-π（信息量更大）。
    - **阳离子-π**：正电与芳香相遇时同理。

    这两类由 :func:`classify_pair_multi` 统一处理。
    """
    hits = classify_pair_multi(
        res_a, res_b, distance,
        include_histidine=include_histidine,
        include_ambiguous_donors=include_ambiguous_donors,
    )
    return hits[0] if hits else None


def classify_pair_multi(
    res_a: str,
    res_b: str,
    distance: float,
    include_histidine: bool = False,
    include_ambiguous_donors: bool = False,
) -> list[str]:
    """判定残基对的全部相互作用类型（按优先级排序）。

    与 :func:`classify_pair` 的区别：这里返回**全部**匹配的类型，
    而不是只返回第一个。用于"一个残基对既是盐桥又是氢键"这类情况 ——
    这在实际结构里很常见（比如 ARG 的胍基与 ASP 的羧基）。

    参数
    ----
    res_a, res_b
        两侧残基的三字母码。
    distance
        最小原子间距（Å）。
    include_histidine
        是否把 HIS 当作正电残基。
    include_ambiguous_donors
        是否把 ASP/GLU 当作弱供体。

    返回
    ----
    相互作用类型名列表，按优先级从高到低。

    说明
    ----
    当**任一侧**是非标准残基（配体等）时，返回 ``["Ligand_contact"]`` ——
    不套用氨基酸分类表。这不是"未分类"，而是"以另一种方式分类"：
    它明确告诉使用者"这是一个配体接触，具体的相互作用类型需要用
    PLIP 或其他专用工具判定"。
    """
    a = _norm(res_a)
    b = _norm(res_b)

    hits: list[str] = []

    # ---- 配体接触：不做残基级化学分类 ----
    kind_a = residue_kind(a)
    kind_b = residue_kind(b)
    if kind_a != "standard" or kind_b != "standard":
        if kind_a == "water" and kind_b == "water":
            return []
        return [LIGAND_CONTACT]

    # ---- 氢键 ----
    hb_cutoff = INTERACTION_TYPES[0].max_distance
    if distance <= hb_cutoff:
        donor_a = is_donor_residue(a, include_ambiguous_donors)
        donor_b = is_donor_residue(b, include_ambiguous_donors)
        acceptor_a = is_acceptor_residue(a)
        acceptor_b = is_acceptor_residue(b)
        if (donor_a and acceptor_b) or (donor_b and acceptor_a):
            hits.append("Hydrogen_bond")

    # ---- 盐桥 ----
    salt_cutoff = INTERACTION_TYPES[1].max_distance
    if distance <= salt_cutoff:
        pos_a = is_positive_residue(a, include_histidine)
        pos_b = is_positive_residue(b, include_histidine)
        neg_a = is_negative_residue(a)
        neg_b = is_negative_residue(b)
        if (pos_a and neg_b) or (pos_b and neg_a):
            hits.append("Salt_bridge")

    # ---- π-π 堆积 ----
    pi_cutoff = INTERACTION_TYPES[3].max_distance
    if distance <= pi_cutoff:
        if is_aromatic_residue(a) and is_aromatic_residue(b):
            hits.append("Pi_Pi")

    # ---- 阳离子-π ----
    cat_cutoff = INTERACTION_TYPES[4].max_distance
    if distance <= cat_cutoff:
        pos_a = is_positive_residue(a, include_histidine)
        pos_b = is_positive_residue(b, include_histidine)
        aro_a = is_aromatic_residue(a)
        aro_b = is_aromatic_residue(b)
        if (pos_a and aro_b) or (pos_b and aro_a):
            hits.append("Cation_Pi")

    # ---- 疏水接触 ----
    ph_cutoff = INTERACTION_TYPES[2].max_distance
    if distance <= ph_cutoff:
        if is_hydrophobic_residue(a) and is_hydrophobic_residue(b):
            hits.append("Hydrophobic")

    return hits


def distance_cutoff(interaction_type: str) -> float:
    """取某个相互作用类型的距离阈值（Å）。"""
    for rule in INTERACTION_TYPES:
        if rule.name == interaction_type:
            return rule.max_distance
    raise KeyError(f"未知相互作用类型: {interaction_type!r}")


def describe_rules() -> list[dict]:
    """返回全部接触规则的描述（供文档与报告使用）。

    包含 5 种残基级化学类型，外加配体接触这一特殊类别。
    """
    rules = [
        {
            "name": rule.name,
            "label": rule.label,
            "max_distance": rule.max_distance,
            "description": rule.description,
        }
        for rule in INTERACTION_TYPES
    ]
    rules.append({
        "name": LIGAND_CONTACT,
        "label": "配体接触",
        "max_distance": None,
        "description": (
            "接触涉及非标准残基（配体、水、离子等），不套用残基级化学分类；"
            "具体类型需用 PLIP 等专用工具判定"
        ),
    })
    return rules
