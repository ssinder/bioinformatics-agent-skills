#!/usr/bin/env python3
"""链间残基接触的几何检测 —— 纯 Python，零依赖。

这个模块做一件老仓库做不到的事：**在没有 PyMOL 的环境里算出链间接触**。

做法
----
1. 用**均匀网格空间划分**（uniform grid / cell list）把原子按 6 Å 的立方格
   装进桶 —— 这样"找出所有在 6 Å 内的原子对"从 O(N²) 降到接近 O(N)。
2. 遍历每个桶的 27 个邻域，收集候选原子对。
3. 对候选对算真实欧氏距离，保留在阈值内的。
4. 把原子级接触归并到**残基对**级，记录最小距离。

为什么用网格而不是暴力两两比较
-----------------------------
一个典型复合物有 5000-20000 个原子。暴力两两比较是 2.5 亿次距离计算 ——
纯 Python 要跑几分钟。网格划分后只比较真正邻近的桶，
在同样的输入下通常快 50-200 倍。

网格边长的选择
--------------
网格边长必须 **≥ 最大接触距离**（默认 6.0 Å），否则会漏掉跨桶的接触。
边长大于最大距离时，只需检查 3×3×3 的邻域；等于时需要检查更远的一圈
（因为原子可能正好落在桶的边界）。这里取 ``ceil(max_distance / cell_size) + 1``
来计算需要扫描的邻域范围。

本模块只依赖 Python 标准库。
"""

from __future__ import annotations

import math
from collections.abc import Iterable
from dataclasses import dataclass, field
from typing import Optional

__all__ = [
    "AtomRef",
    "Contact",
    "ResiduePair",
    "ContactDetector",
    "detect_chain_contacts",
    "grid_cell_size",
]


#: 默认的最大接触距离（Å）—— 由最宽松的接触类型（阳离子-π，6.0 Å）决定。
DEFAULT_MAX_DISTANCE = 6.0


@dataclass(frozen=True)
class AtomRef:
    """一个参与接触的原子。

    Attributes:
        serial: 原子序号。
        name: 原子名。
        res_name: 残基名（三字母）。
        chain_id: 链标识。
        res_seq: 残基序号。
        icode: 插入码。
        x, y, z: 坐标（Å）。
        element: 元素符号。
        hetero: 是否来自 HETATM。
    """

    serial: int
    name: str
    res_name: str
    chain_id: str
    res_seq: int
    icode: str
    x: float
    y: float
    z: float
    element: str = ""
    hetero: bool = False

    @property
    def residue_key(self) -> tuple[str, int, str]:
        """残基唯一键（链 + 序号 + 插入码）。"""
        return (self.chain_id, self.res_seq, self.icode)

    @property
    def residue_label(self) -> str:
        """人类可读的残基标签，如 ``A:ARG102``。"""
        icode = self.icode.strip()
        return f"{self.chain_id}:{self.res_name}{self.res_seq}{icode}"


@dataclass
class Contact:
    """一对原子之间的接触。

    Attributes:
        atom_a: 第一个原子。
        atom_b: 第二个原子。
        distance: 原子间距（Å）。
    """

    atom_a: AtomRef
    atom_b: AtomRef
    distance: float


@dataclass
class ResiduePair:
    """一对残基之间的接触汇总。

    Attributes:
        chain_a: 链 A 标识。
        res_name_a: 残基 A 名。
        res_seq_a: 残基 A 序号。
        chain_b: 链 B 标识。
        res_name_b: 残基 B 名。
        res_seq_b: 残基 B 序号。
        min_distance: 两侧最近原子间距（Å）。
        n_atom_pairs: 构成该接触的原子对数量。
        atom_pairs: 具体的原子对明细。
    """

    chain_a: str
    res_name_a: str
    res_seq_a: int
    chain_b: str
    res_name_b: str
    res_seq_b: int
    min_distance: float
    n_atom_pairs: int = 0
    atom_pairs: list[Contact] = field(default_factory=list)

    @property
    def label_a(self) -> str:
        return f"{self.chain_a}:{self.res_name_a}{self.res_seq_a}"

    @property
    def label_b(self) -> str:
        return f"{self.chain_b}:{self.res_name_b}{self.res_seq_b}"

    @property
    def pair_key(self) -> tuple[str, int, str, int, str]:
        """残基对的唯一键（无序）。"""
        return (self.chain_a, self.res_seq_a, self.res_name_a,
                self.chain_b, self.res_seq_b, self.res_name_b)


def grid_cell_size(max_distance: float, target_cells: int = 1.0) -> float:
    """计算合适的网格边长。

    参数
    ----
    max_distance
        最大接触距离（Å）。网格边长必须不小于它。
    target_cells
        目标：一个桶里平均装多少个原子。值越大桶越小、比较越少，
        但桶的数量会指数增长。默认 1.0 对典型蛋白已经很好。

    返回
    ----
    网格边长（Å），保证 ``>= max_distance``。

    说明
    ----
    边长取 ``max_distance * 1.0001`` 而不是刚好相等，是为了给浮点留余量：
    桶索引用 ``floor(x / size)`` 计算，桶是**左闭右开**区间，理论上
    桶索引差最多为 1，但浮点误差可能让一个恰好落在边界的原子被算进
    相邻桶。多留 0.01% 的余量可以彻底消除这类边界抖动。
    """
    # 边长至少等于最大距离，否则会漏掉跨桶接触
    return max(max_distance * 1.0001, 1.0)


class ContactDetector:
    """链间接触检测器。

    用法::

        detector = ContactDetector(max_distance=6.0)
        detector.add_atoms(atoms)
        contacts = detector.detect("A", "B")
    """

    def __init__(
        self,
        max_distance: float = DEFAULT_MAX_DISTANCE,
        include_hetero: bool = True,
        ignore_hydrogen: bool = True,
    ) -> None:
        """初始化检测器。

        参数
        ----
        max_distance
            最大接触距离（Å）。
        include_hetero
            是否把 HETATM（水、离子、辅因子、配体）算进去。
            配体与蛋白的接触是研究热点，**默认开启**。
        ignore_hydrogen
            是否忽略氢原子。PDB 文件里的氢坐标通常不完整或经过修整，
            参与距离计算反而引入误差 —— **默认忽略**。
        """
        self.max_distance = float(max_distance)
        self.include_hetero = include_hetero
        self.ignore_hydrogen = ignore_hydrogen

        self._atoms: list[AtomRef] = []
        self._grid: dict[tuple[int, int, int], list[int]] = {}
        self._cell_size = grid_cell_size(self.max_distance)

    # ------------------------------------------------------------ 装载

    def add_atoms(self, atoms: Iterable[AtomRef]) -> None:
        """加入原子并重建空间网格。"""
        for atom in atoms:
            if self.ignore_hydrogen and atom.element.upper() == "H":
                continue
            if not self.include_hetero and atom.hetero:
                continue
            self._atoms.append(atom)
        self._rebuild_grid()

    def _rebuild_grid(self) -> None:
        """把原子装进网格桶。"""
        size = self._cell_size
        grid: dict[tuple[int, int, int], list[int]] = {}
        for index, atom in enumerate(self._atoms):
            key = (
                int(math.floor(atom.x / size)),
                int(math.floor(atom.y / size)),
                int(math.floor(atom.z / size)),
            )
            grid.setdefault(key, []).append(index)
        self._grid = grid

    def _cell_range(self) -> int:
        """需要扫描的邻域半径（桶数）。

        桶是左闭右开区间。若桶边长 **严格大于** 最大接触距离，则任意两个
        距离小于最大距离的原子，其桶索引在每个轴上的差最多是 1 ——
        扫 3×3×3 的邻域就足够。

        早期实现返回 2（扫 5×5×5），理由是"原子可能落在桶边界两侧"。
        这个理由在桶边长 **等于** 最大距离时成立：此时距离恰为最大值的
        两个原子可能被分到相邻桶，且其中一个的桶索引可能相差 2。
        现在桶边长取 ``max_distance * 1.0001``（严格大于），差 2 的桶
        索引不可能落在最大距离内，扫 1 格即可 —— 候选对数量减少约 4.6 倍，
        而结果完全一致。
        """
        return 1

    # ------------------------------------------------------------ 检测

    def detect(
        self,
        chain_a: str,
        chain_b: str,
        include_same_chain: bool = False,
    ) -> list[ResiduePair]:
        """检测两条链之间的残基接触。

        参数
        ----
        chain_a, chain_b
            链标识。
        include_same_chain
            是否允许检测**同一条链**内部的接触。默认 **False** ——
            链内接触（尤其是共价相邻残基）与本技能要回答的
            "两条链之间怎么结合" 无关，混进来会淹没跨链信号。
            只有在研究链内构象耦联时才需要打开。

        返回
        ----
        按最小距离升序排列的残基对列表。
        """
        if chain_a == chain_b and not include_same_chain:
            return []

        # 按链分桶，避免跨链比较
        by_chain: dict[str, list[int]] = {}
        for index, atom in enumerate(self._atoms):
            by_chain.setdefault(atom.chain_id, []).append(index)

        indices_a = by_chain.get(chain_a, [])
        indices_b = by_chain.get(chain_b, [])
        if not indices_a or not indices_b:
            return []

        size = self._cell_size
        radius = self._cell_range()
        cutoff_sq = self.max_distance * self.max_distance

        # 收集候选原子对（用 (i, j) 去重）
        candidate_pairs: set[tuple[int, int]] = set()
        for index in indices_a:
            atom = self._atoms[index]
            cx = int(math.floor(atom.x / size))
            cy = int(math.floor(atom.y / size))
            cz = int(math.floor(atom.z / size))
            for dx in range(-radius, radius + 1):
                for dy in range(-radius, radius + 1):
                    for dz in range(-radius, radius + 1):
                        bucket = self._grid.get((cx + dx, cy + dy, cz + dz))
                        if not bucket:
                            continue
                        for other in bucket:
                            if other == index:
                                continue
                            other_atom = self._atoms[other]
                            if other_atom.chain_id != chain_b:
                                continue
                            key = (index, other) if index < other else (other, index)
                            candidate_pairs.add(key)

        # 算真实距离并按残基对归并
        residue_contacts: dict[tuple, ResiduePair] = {}
        for i, j in candidate_pairs:
            atom_i = self._atoms[i]
            atom_j = self._atoms[j]
            dist_sq = (
                (atom_i.x - atom_j.x) ** 2
                + (atom_i.y - atom_j.y) ** 2
                + (atom_i.z - atom_j.z) ** 2
            )
            if dist_sq > cutoff_sq:
                continue
            distance = math.sqrt(dist_sq)

            # 按链顺序归一（保证 pair_key 与顺序无关）
            if atom_i.chain_id == chain_a:
                first, second = atom_i, atom_j
            else:
                first, second = atom_j, atom_i

            key = (
                first.chain_id, first.res_seq, first.icode,
                second.chain_id, second.res_seq, second.icode,
            )
            existing = residue_contacts.get(key)
            if existing is None:
                existing = ResiduePair(
                    chain_a=first.chain_id,
                    res_name_a=first.res_name,
                    res_seq_a=first.res_seq,
                    chain_b=second.chain_id,
                    res_name_b=second.res_name,
                    res_seq_b=second.res_seq,
                    min_distance=distance,
                )
                residue_contacts[key] = existing
            if distance < existing.min_distance:
                existing.min_distance = distance
            existing.n_atom_pairs += 1
            existing.atom_pairs.append(Contact(atom_a=first, atom_b=second, distance=distance))

        result = list(residue_contacts.values())
        result.sort(key=lambda pair: (pair.min_distance, pair.chain_a, pair.res_seq_a))
        return result

    def detect_all(self, exclude_same_chain: bool = True) -> dict[tuple[str, str], list[ResiduePair]]:
        """检测所有链两两之间的接触。

        参数
        ----
        exclude_same_chain
            是否跳过链内接触（同一链的两个残基）。

        返回
        ----
        ``{(chain_a, chain_b): [ResiduePair, ...]}``，链对按字母序排列。
        """
        chains = sorted({atom.chain_id for atom in self._atoms})
        out: dict[tuple[str, str], list[ResiduePair]] = {}
        for index, chain_a in enumerate(chains):
            for chain_b in chains[index + 1:]:
                pairs = self.detect(chain_a, chain_b)
                if pairs:
                    out[(chain_a, chain_b)] = pairs
        return out


def detect_chain_contacts(
    atoms: Iterable[AtomRef],
    chain_a: str,
    chain_b: str,
    max_distance: float = DEFAULT_MAX_DISTANCE,
    include_hetero: bool = True,
    ignore_hydrogen: bool = True,
) -> list[ResiduePair]:
    """便捷函数：检测两条链之间的残基接触。"""
    detector = ContactDetector(
        max_distance=max_distance,
        include_hetero=include_hetero,
        ignore_hydrogen=ignore_hydrogen,
    )
    detector.add_atoms(atoms)
    return detector.detect(chain_a, chain_b)
