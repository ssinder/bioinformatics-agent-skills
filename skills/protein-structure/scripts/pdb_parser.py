#!/usr/bin/env python3
"""纯 Python PDB 解析器（不依赖 BioPython）。

按 PDB v3.3 原子记录的固定列定义逐行解析 ``ATOM`` / ``HETATM`` 记录，
构建 ``Structure -> Model -> Chain -> Residue -> Atom`` 的层级结构，
并提供原子类型统计、几何量（回旋半径 / 陀螺张量 / 瞬时偶极）等计算。

固定列（0-based 切片）::

    record      [0:6]
    serial      [6:11]
    name        [12:16]
    altLoc      [16]
    resName     [17:20]
    chainID     [21]
    resSeq      [22:26]
    iCode       [26]
    x,y,z       [30:38], [38:46], [46:54]
    occupancy   [54:60]
    tempFactor  [60:66]
    element     [76:78]

Usage:
    from pdb_parser import parse_pdb_file
    structure = parse_pdb_file("1crn.pdb")
    print(structure.summary())
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Iterator
from dataclasses import dataclass, field
from pathlib import Path

__all__ = [
    "Atom",
    "Residue",
    "Chain",
    "Model",
    "Structure",
    "parse_pdb_file",
    "parse_pdb_text",
    "dihedral_angle",
    "angle_between",
]

# 标准氨基酸三字母 -> 单字母
AA3_TO_AA1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "SEC": "U", "PYL": "O",
}

# 常见残基的骨架原子名（用于二级结构计算）
BACKBONE_ATOMS = ("N", "CA", "C", "O", "OXT")

# 各元素的有效共价半径（Å），用于 SASA/接触判据
VDW_RADII = {
    "H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80,
    "P": 1.80, "F": 1.47, "CL": 1.75, "BR": 1.85, "I": 1.98,
    "MG": 1.73, "ZN": 1.39, "CA": 1.95, "NA": 2.27, "K": 2.75,
    "FE": 1.80, "MN": 1.61, "CU": 1.40, "SE": 1.90,
}


def _safe_float(text: str, default: float = 0.0) -> float:
    """把定宽字段转成 float，失败时返回默认值。"""
    try:
        return float(text.strip())
    except (TypeError, ValueError):
        return default


def _safe_int(text: str, default: int = 0) -> int:
    """把定宽字段转成 int，失败时返回默认值。"""
    try:
        return int(text.strip())
    except (TypeError, ValueError):
        return default


@dataclass
class Atom:
    """单个原子。

    Attributes:
        serial: 原子序号。
        name: 原子名（已去除空白）。
        alt_loc: 交替位置标记。
        res_name: 残基名。
        chain_id: 链标识。
        res_seq: 残基序号。
        icode: 插入码。
        coord: ``(x, y, z)`` 坐标。
        occupancy: 占据率。
        bfactor: B 因子。
        element: 元素符号。
        hetatm: 是否来自 HETATM 记录。
    """

    serial: int
    name: str
    alt_loc: str
    res_name: str
    chain_id: str
    res_seq: int
    icode: str
    coord: tuple[float, float, float]
    occupancy: float
    bfactor: float
    element: str
    hetatm: bool = False

    @property
    def x(self) -> float:
        return self.coord[0]

    @property
    def y(self) -> float:
        return self.coord[1]

    @property
    def z(self) -> float:
        return self.coord[2]

    @property
    def is_hydrogen(self) -> bool:
        return self.element.upper() == "H"

    @property
    def is_backbone(self) -> bool:
        return self.name in BACKBONE_ATOMS

    @property
    def radius(self) -> float:
        """原子的范德华半径（未知元素取碳的值）。"""
        return VDW_RADII.get(self.element.upper(), 1.70)

    def distance_to(self, other: Atom) -> float:
        """到另一个原子的欧氏距离（Å）。"""
        return math.dist(self.coord, other.coord)


@dataclass
class Residue:
    """一个残基（按 ``(chain_id, res_seq, icode)`` 唯一）。"""

    res_name: str
    res_seq: int
    chain_id: str
    icode: str = ""
    atoms: list[Atom] = field(default_factory=list)
    hetatm: bool = False

    def get(self, name: str) -> Atom | None:
        """按原子名取原子（大小写不敏感），取不到返回 None。"""
        upper = name.upper()
        for atom in self.atoms:
            if atom.name.upper() == upper:
                return atom
        return None

    @property
    def atoms_dict(self) -> dict[str, Atom]:
        """原子名 -> 原子 的映射。"""
        return {a.name: a for a in self.atoms}

    @property
    def is_standard(self) -> bool:
        """是否为标准氨基酸（含 MSE 变体）。"""
        return self.res_name in AA3_TO_AA1

    @property
    def one_letter(self) -> str:
        """单字母代码，非标准残基用 ``X``。"""
        return AA3_TO_AA1.get(self.res_name, "X")

    @property
    def label(self) -> str:
        """人类可读标签，如 ``A:TYR42``。"""
        return f"{self.chain_id}:{self.res_name}{self.res_seq}{self.icode.strip()}"

    def backbone(self) -> dict[str, tuple[float, float, float]]:
        """返回骨架原子坐标字典（仅包含实际存在的原子）。"""
        out: dict[str, tuple[float, float, float]] = {}
        for name in ("N", "CA", "C", "O"):
            atom = self.get(name)
            if atom is not None:
                out[name] = atom.coord
        return out

    def centroid(self) -> tuple[float, float, float]:
        """残基所有原子的几何中心。"""
        if not self.atoms:
            return (0.0, 0.0, 0.0)
        n = len(self.atoms)
        return (
            sum(a.x for a in self.atoms) / n,
            sum(a.y for a in self.atoms) / n,
            sum(a.z for a in self.atoms) / n,
        )


@dataclass
class Chain:
    """一条链（同一 ``chain_id`` 的残基集合）。"""

    chain_id: str
    residues: list[Residue] = field(default_factory=list)

    def __iter__(self) -> Iterator[Residue]:
        return iter(self.residues)

    @property
    def n_residues(self) -> int:
        return len(self.residues)

    @property
    def n_atoms(self) -> int:
        return sum(len(r.atoms) for r in self.residues)

    def sequence(self, only_standard: bool = True) -> str:
        """输出单字母序列。"""
        return "".join(r.one_letter for r in self.residues if (r.is_standard or not only_standard))


@dataclass
class Model:
    """一个模型（PDB 可含多个 MODEL 记录，通常为 NMR 多模型）。"""

    model_id: int
    chains: list[Chain] = field(default_factory=list)

    def __iter__(self) -> Iterator[Chain]:
        return iter(self.chains)

    @property
    def n_chains(self) -> int:
        return len(self.chains)

    @property
    def n_residues(self) -> int:
        return sum(c.n_residues for c in self.chains)

    @property
    def n_atoms(self) -> int:
        return sum(c.n_atoms for c in self.chains)

    def chain(self, chain_id: str) -> Chain | None:
        """按 ID 取链。"""
        for c in self.chains:
            if c.chain_id == chain_id:
                return c
        return None


@dataclass
class Structure:
    """整个结构文件。"""

    structure_id: str
    models: list[Model] = field(default_factory=list)
    header: dict[str, str] = field(default_factory=dict)

    @property
    def model(self) -> Model:
        """第一个模型（没有模型时返回空模型）。"""
        return self.models[0] if self.models else Model(model_id=1)

    @property
    def n_models(self) -> int:
        return len(self.models)

    @property
    def n_chains(self) -> int:
        return self.model.n_chains

    @property
    def n_residues(self) -> int:
        return self.model.n_residues

    @property
    def n_atoms(self) -> int:
        return self.model.n_atoms

    def atoms(self) -> list[Atom]:
        """第一个模型的所有原子。"""
        return [a for chain in self.model for res in chain for a in res.atoms]

    def residues(self) -> list[Residue]:
        """第一个模型的所有残基。"""
        return [r for chain in self.model for r in chain]

    def summary(self) -> dict:
        """返回结构概要字典。"""
        chains = self.model.chains
        return {
            "structure_id": self.structure_id,
            "n_models": self.n_models,
            "n_chains": self.n_chains,
            "n_residues": self.n_residues,
            "n_atoms": self.n_atoms,
            "chain_ids": [c.chain_id for c in chains],
            "seq_length": max((c.n_residues for c in chains), default=0),
            "has_ca": any(r.get("CA") is not None for r in self.residues()),
            "element_counts": element_statistics(self.atoms()),
            "hetatm_residues": sorted(
                {r.res_name for r in self.residues() if r.hetatm}
            ),
        }


def element_statistics(atoms: Iterable[Atom]) -> dict[str, int]:
    """统计各元素的原子数量（按元素符号排序）。"""
    counts: dict[str, int] = {}
    for atom in atoms:
        key = atom.element.upper() or atom.name[0].upper()
        counts[key] = counts.get(key, 0) + 1
    return dict(sorted(counts.items()))


def angle_between(v1: tuple[float, float, float], v2: tuple[float, float, float]) -> float:
    """两个向量的夹角（度）。"""
    n1 = math.sqrt(sum(c * c for c in v1))
    n2 = math.sqrt(sum(c * c for c in v2))
    if n1 == 0.0 or n2 == 0.0:
        return 0.0
    cos = sum(a * b for a, b in zip(v1, v2)) / (n1 * n2)
    cos = max(-1.0, min(1.0, cos))
    return math.degrees(math.acos(cos))


def dihedral_angle(
    p0: tuple[float, float, float],
    p1: tuple[float, float, float],
    p2: tuple[float, float, float],
    p3: tuple[float, float, float],
) -> float:
    """计算四点二面角（度，范围 -180 ~ +180）。

    采用 IUPAC 约定：扭转角是 ``p0-p1-p2-p3`` 绕 ``p1→p2`` 轴的旋转角，
    从 ``p1`` 看向 ``p2`` 时为负（例如 α-螺旋的 φ ≈ -57°、ψ ≈ -47°）。

    实现基于投影法：把 ``p0-p1`` 与 ``p3-p2`` 都投影到垂直于主轴
    ``p1→p2`` 的平面上，投影向量的夹角即为扭转角；符号由
    ``(主轴 × 投影1) · 投影2`` 的正负决定。
    """
    b0 = (p0[0] - p1[0], p0[1] - p1[1], p0[2] - p1[2])
    b1 = (p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2])
    b2 = (p3[0] - p2[0], p3[1] - p2[1], p3[2] - p2[2])

    b1_len = math.sqrt(sum(c * c for c in b1))
    if b1_len < 1e-9:
        return 0.0
    axis = (b1[0] / b1_len, b1[1] / b1_len, b1[2] / b1_len)

    # v: b0 在垂直于 axis 平面内的投影（指向 p0 一侧）
    d0 = sum(b0[i] * axis[i] for i in range(3))
    v = (b0[0] - d0 * axis[0], b0[1] - d0 * axis[1], b0[2] - d0 * axis[2])
    # w: b2 在垂直于 axis 平面内的投影
    d2 = sum(b2[i] * axis[i] for i in range(3))
    w = (b2[0] - d2 * axis[0], b2[1] - d2 * axis[1], b2[2] - d2 * axis[2])

    x = sum(v[i] * w[i] for i in range(3))
    cr = (
        axis[1] * v[2] - axis[2] * v[1],
        axis[2] * v[0] - axis[0] * v[2],
        axis[0] * v[1] - axis[1] * v[0],
    )
    y = sum(cr[i] * w[i] for i in range(3))
    return math.degrees(math.atan2(y, x))


def _parse_atom_line(line: str, serial_hint: int) -> Atom | None:
    """解析一行 ATOM/HETATM 记录。"""
    if len(line) < 54:
        return None
    record = line[0:6].strip()
    if record not in ("ATOM", "HETATM"):
        return None

    name = line[12:16].strip()
    alt_loc = line[16:17].strip() if len(line) > 16 else ""
    res_name = line[17:20].strip()
    chain_id = line[21:22].strip() if len(line) > 21 else ""
    res_seq = _safe_int(line[22:26]) if len(line) >= 26 else 0
    icode = line[26:27] if len(line) >= 27 else ""
    coord = (
        _safe_float(line[30:38]),
        _safe_float(line[38:46]),
        _safe_float(line[46:54]),
    )
    occupancy = _safe_float(line[54:60], 1.0) if len(line) >= 60 else 1.0
    bfactor = _safe_float(line[60:66], 0.0) if len(line) >= 66 else 0.0
    element = ""
    if len(line) >= 78:
        element = line[76:78].strip().upper()
    if not element:
        # PDB 中 element 字段常留空，此时从原子名推断（首字母为数字时取第二位）
        raw = line[12:16]
        stripped = raw.strip()
        element = stripped[0].upper() if stripped else "X"
        if stripped and stripped[0].isdigit() and len(stripped) > 1:
            element = stripped[1].upper()

    return Atom(
        serial=_safe_int(line[6:11], serial_hint) or serial_hint,
        name=name,
        alt_loc=alt_loc,
        res_name=res_name,
        chain_id=chain_id,
        res_seq=res_seq,
        icode=icode,
        coord=coord,
        occupancy=occupancy,
        bfactor=bfactor,
        element=element,
        hetatm=(record == "HETATM"),
    )


def parse_pdb_text(text: str, structure_id: str = "structure") -> Structure:
    """从 PDB 文本解析出 :class:`Structure`。

    处理要点：

    - 只保留 ``ATOM`` / ``HETATM`` 记录；
    - 交替位置（altLoc）只取 ``' '``、``'A'`` 或第一个出现的构象；
    - 同一 ``(chain, resSeq, iCode)`` 的原子合并到同一残基；
    - ``MODEL`` / ``ENDMDL`` 记录切分多模型。
    """
    structure = Structure(structure_id=structure_id)
    current_model: Model | None = None
    model_index = 0
    residue_index: dict[tuple[str, int, str], Residue] = {}
    chain_index: dict[str, Chain] = {}
    serial_hint = 0

    for raw in text.splitlines():
        line = raw.rstrip("\n\r")
        if not line:
            continue
        record = line[0:6].strip()

        if record == "MODEL":
            model_index += 1
            current_model = Model(model_id=model_index)
            structure.models.append(current_model)
            residue_index = {}
            chain_index = {}
            continue
        if record == "ENDMDL":
            current_model = None
            continue
        if record in ("HEADER", "TITLE", "COMPND", "SOURCE", "REMARK", "EXPDTA"):
            key = record
            value = line[10:].strip() if len(line) > 10 else line.strip()
            if key in ("HEADER", "TITLE", "EXPDTA"):
                structure.header.setdefault(key, value)
            continue

        if record not in ("ATOM", "HETATM"):
            continue

        serial_hint += 1
        atom = _parse_atom_line(line, serial_hint)
        if atom is None:
            continue
        # 交替位置：只保留主构象
        if atom.alt_loc not in ("", "A", "1"):
            continue

        if current_model is None:
            model_index = 1
            current_model = Model(model_id=1)
            structure.models.append(current_model)

        chain = chain_index.get(atom.chain_id)
        if chain is None:
            chain = Chain(chain_id=atom.chain_id)
            chain_index[atom.chain_id] = chain
            current_model.chains.append(chain)

        key = (atom.chain_id, atom.res_seq, atom.icode)
        residue = residue_index.get(key)
        if residue is None:
            residue = Residue(
                res_name=atom.res_name,
                res_seq=atom.res_seq,
                chain_id=atom.chain_id,
                icode=atom.icode,
                hetatm=atom.hetatm,
            )
            residue_index[key] = residue
            chain.residues.append(residue)

        # 同一原子名重复出现（多模型或多重占据）时保留第一个
        if residue.get(atom.name) is None:
            residue.atoms.append(atom)

    return structure


def parse_pdb_file(path: str | Path) -> Structure:
    """解析单个 PDB 文件。"""
    p = Path(path)
    text = p.read_text(encoding="utf-8", errors="replace")
    return parse_pdb_text(text, structure_id=p.stem)


def _main(argv: list[str] | None = None) -> int:
    """命令行入口：解析 PDB 并打印结构概要。"""
    import argparse

    parser = argparse.ArgumentParser(
        prog="pdb_parser", description="纯 Python PDB 解析器（命令行自检工具）"
    )
    parser.add_argument("input", help="PDB 文件路径")
    args = parser.parse_args(argv)

    structure = parse_pdb_file(args.input)
    summary = structure.summary()
    for key, value in summary.items():
        logger.info("%s: %s", key, value)
    return 0


if __name__ == "__main__":
    import sys

    sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "packages"))
    from biopipe_core import get_logger

    logger = get_logger("biopipe.protein-structure")
    sys.exit(_main())
