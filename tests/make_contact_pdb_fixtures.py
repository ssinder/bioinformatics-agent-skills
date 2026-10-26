#!/usr/bin/env python3
"""生成链间接触分析的 PDB 测试夹具。

夹具是**手工构造**的两链复合物。残基采用标准氨基酸理想坐标，再通过
**平移**把 B 侧残基的接触原子精确送到 A 侧接触原子的指定距离处 ——
这样每种相互作用类型都能被精确复现。

关键设计
========

**每个残基对占据独立的空间区域**（沿 x 轴按 30 Å 间隔排开），
避免不同配对的原子意外靠近。老实现把 6 个配对全部按 x 轴 16 Å 间距排列，
导致不同残基的主链原子互相穿插，出现了 1.18 Å 这种物理上不可能的距离
（原子会重叠）。本实现用 **30 Å 间距**（大于两倍最大接触距离 6 Å）
彻底隔离各配对。

期望的相互作用
==============
============  ==================  ==================  ============
残基对          接触原子            距离（Å）           期望类型
============  ==================  ==================  ============
ARG10 - ASP60   NH1 ... OD1        2.9                氢键 + 盐桥
LYS20 - GLU70   NZ  ... OE1        3.0                氢键 + 盐桥
LYS25 - PHE90   NZ  ... CG         5.5                阳离子-π
PHE30 - TYR80   CG  ... CG         4.5                π-π 堆积
ILE40 - LEU100  CG1 ... CG         4.2                疏水接触
GLY50 - GLY110  CA  ... CA         3.8                无类型
============  ==================  ==================  ============

另外还包含：

- 一个 ``HETATM`` 配体（``LIG`` 苯环），用于验证配体接触被单独标记
- 一个氢原子，用于验证默认忽略氢
- 一个含两个 ``MODEL`` 的文件，用于验证只取第一个 MODEL

生成器会**自检**：

1. 期望的相互作用类型都出现
2. GLY-GLY 不产生任何类型
3. **所有原子两两间距 ≥ 2.0 Å** —— 防止原子重叠

夹具坏了要立刻知道，而不是让下游测试莫名其妙地失败。

用法::

    python tests/make_contact_pdb_fixtures.py
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
FIXTURES_DIR = REPO_ROOT / "tests" / "fixtures_contacts"

#: 不同残基对之间的 x 轴间距（Å）。必须远大于最大接触距离（6 Å），
#: 否则相邻配对的侧链会意外靠近，污染接触检测结果。
PAIR_SPACING = 30.0

#: 允许的最小原子间距（Å）。低于此值认为原子重叠，夹具无效。
#:
#: 注意这个阈值只适用于**不同残基之间**的原子。残基内部的共价键
#: （C-N 约 1.33 Å、C-O 约 1.23 Å、C-C 约 1.52 Å）本来就短得多，
#: 用同一阈值检查会把正常结构误判为"重叠"。
MIN_ATOM_SEPARATION = 2.0


# ---------------------------------------------------------------- 理想残基坐标
# 标准氨基酸的理想坐标（主链 N 位于原点，z 轴朝前）

ARG_LOCAL = {
    "N": (0.000, 0.000, 0.000), "CA": (1.458, 0.000, 0.000),
    "C": (2.009, 1.420, 0.000), "O": (1.251, 2.390, 0.000),
    "CB": (2.009, -0.518, 1.426), "CG": (3.478, -0.341, 1.685),
    "CD": (4.089, 0.828, 2.301), "NE": (5.443, 1.034, 2.712),
    "CZ": (5.891, 2.360, 3.327), "NH1": (5.062, 3.331, 3.662),
    "NH2": (6.986, 2.685, 3.598),
}

LYS_LOCAL = {
    "N": (0.000, 0.000, 0.000), "CA": (1.458, 0.000, 0.000),
    "C": (2.009, 1.420, 0.000), "O": (1.251, 2.390, 0.000),
    "CB": (2.009, -0.518, 1.426), "CG": (3.478, -0.341, 1.685),
    "CD": (3.877, 1.093, 2.160), "CE": (4.346, 0.824, 3.221),
    "NZ": (3.627, 1.733, 4.061),
}

ASP_LOCAL = {
    "N": (0.000, 0.000, 0.000), "CA": (1.458, 0.000, 0.000),
    "C": (2.009, 1.420, 0.000), "O": (1.251, 2.390, 0.000),
    "CB": (2.009, -0.518, 1.426), "CG": (3.477, -0.341, 1.685),
    "OD1": (3.541, -1.269, 2.523), "OD2": (4.522, 0.116, 1.891),
}

GLU_LOCAL = {
    "N": (0.000, 0.000, 0.000), "CA": (1.458, 0.000, 0.000),
    "C": (2.009, 1.420, 0.000), "O": (1.251, 2.390, 0.000),
    "CB": (2.009, -0.518, 1.426), "CG": (3.478, -0.341, 1.685),
    "CD": (3.877, 1.093, 2.160), "OE1": (3.541, 2.379, 2.834),
    "OE2": (5.029, 0.881, 2.774),
}

PHE_LOCAL = {
    "N": (0.000, 0.000, 0.000), "CA": (1.458, 0.000, 0.000),
    "C": (2.009, 1.420, 0.000), "O": (1.251, 2.390, 0.000),
    "CB": (2.009, -0.518, 1.426), "CG": (3.478, -0.341, 1.685),
    "CD1": (3.633, -1.362, 2.372), "CD2": (4.516, 0.341, 1.860),
    "CE1": (4.716, -1.850, 3.416), "CE2": (5.599, 0.871, 2.904),
    "CZ": (5.682, -0.162, 4.460),
}

TYR_LOCAL = {
    **PHE_LOCAL,
    "OH": (6.653, -0.516, 4.922),
}

ILE_LOCAL = {
    "N": (0.000, 0.000, 0.000), "CA": (1.458, 0.000, 0.000),
    "C": (2.009, 1.420, 0.000), "O": (1.251, 2.390, 0.000),
    "CB": (2.009, -0.518, 1.426), "CG1": (3.478, -0.341, 1.685),
    "CG2": (1.853, -1.373, 2.341), "CD1": (4.248, 0.451, 2.321),
}

LEU_LOCAL = {
    "N": (0.000, 0.000, 0.000), "CA": (1.458, 0.000, 0.000),
    "C": (2.009, 1.420, 0.000), "O": (1.251, 2.390, 0.000),
    "CB": (2.009, -0.518, 1.426), "CG": (3.478, -0.341, 1.685),
    "CD1": (3.510, -1.293, 2.774), "CD2": (4.516, 0.276, 2.072),
}

GLY_LOCAL = {
    "N": (0.000, 0.000, 0.000), "CA": (1.458, 0.000, 0.000),
    "C": (2.009, 1.420, 0.000), "O": (1.251, 2.390, 0.000),
}

#: 苯环配体（LIG）
LIG_LOCAL = {
    "C1": (0.000, 0.000, 0.000), "C2": (1.390, 0.000, 0.000),
    "C3": (2.085, 1.200, 0.000), "C4": (1.390, 2.400, 0.000),
    "C5": (0.000, 2.400, 0.000), "C6": (-0.695, 1.200, 0.000),
    "O1": (-0.695, -0.780, 0.000),
}

RESIDUE_TEMPLATES = {
    "ARG": ARG_LOCAL, "LYS": LYS_LOCAL, "ASP": ASP_LOCAL, "GLU": GLU_LOCAL,
    "PHE": PHE_LOCAL, "TYR": TYR_LOCAL, "ILE": ILE_LOCAL, "LEU": LEU_LOCAL,
    "GLY": GLY_LOCAL, "LIG": LIG_LOCAL,
}


# ---------------------------------------------------------------- 构造


def _translate(local: dict, offset: tuple[float, float, float]) -> dict:
    """给残基的局部坐标表加上平移量。"""
    return {
        name: (x + offset[0], y + offset[1], z + offset[2])
        for name, (x, y, z) in local.items()
    }


def _dist(a: tuple, b: tuple) -> float:
    """两点距离。"""
    return math.sqrt(sum((p - q) ** 2 for p, q in zip(a, b)))


def _atom_line(
    serial: int,
    name: str,
    res_name: str,
    chain: str,
    res_seq: int,
    xyz: tuple,
    element: str,
    record: str = "ATOM  ",
) -> str:
    """按 PDB 固定列宽格式生成一行原子记录。"""
    x, y, z = xyz
    return (
        f"{record:<6}{serial:>5} {name:<4}{'':1}{res_name:>3} {chain:>1}"
        f"{res_seq:>4}{'':1}   {x:>8.3f}{y:>8.3f}{z:>8.3f}"
        f"{1.00:>6.2f}{20.00:>6.2f}          {element:>2}"
    )


def _element_of(atom_name: str) -> str:
    """从原子名推断元素符号。"""
    name = atom_name.strip()
    if name[:1].isdigit():
        name = name[1:]
    return name[:1].upper()


def _candidate_directions() -> list[tuple[float, float, float]]:
    """候选接触方向的单位向量列表。

    在球面上均匀取点。用黄金角螺旋避免采样点聚集。
    """
    directions: list[tuple[float, float, float]] = []
    golden_angle = math.pi * (3.0 - math.sqrt(5.0))
    n_points = 200
    for i in range(n_points):
        z = 1.0 - (2.0 * i) / (n_points - 1)
        radius = math.sqrt(max(0.0, 1.0 - z * z))
        theta = golden_angle * i
        directions.append((radius * math.cos(theta), radius * math.sin(theta), z))
    return directions


def _min_cross_residue_distance(
    coords_a: dict,
    coords_b: dict,
) -> float:
    """两组坐标之间的最小原子间距。"""
    best = float("inf")
    for _, pa in coords_a.items():
        for _, pb in coords_b.items():
            d = _dist(pa, pb)
            if d < best:
                best = d
    return best


def _place_pair(
    local_a: dict,
    local_b: dict,
    atom_a: str,
    atom_b: str,
    distance: float,
    origin_a: tuple[float, float, float],
) -> tuple[dict, dict]:
    """摆放一对残基，使接触原子间距恰为 ``distance``。

    通过在球面上搜索方向，挑出**残基间最小原子间距最大**的那个方向。
    这样既保证接触原子距离精确，又避免两个残基的其他原子撞在一起
    （老实现固定用 +y 方向，导致 ARG 的 NH2 尖端撞上 GLU 的 OE2 尖端，
    出现 1.18 Å 的非物理距离）。

    参数
    ----
    local_a, local_b
        两侧残基的局部坐标表。
    atom_a, atom_b
        参与接触的原子名。
    distance
        期望的接触距离（Å）。
    origin_a
        A 侧残基的平移量。

    返回
    ----
    ``(coords_a, coords_b)`` —— 两组绝对坐标。
    """
    coords_a = _translate(local_a, origin_a)
    anchor = coords_a[atom_a]

    best_direction = (0.0, 1.0, 0.0)
    best_separation = -1.0

    for direction in _candidate_directions():
        target = tuple(anchor[i] + direction[i] * distance for i in range(3))
        offset_b = tuple(target[i] - local_b[atom_b][i] for i in range(3))
        coords_b = _translate(local_b, offset_b)
        separation = _min_cross_residue_distance(coords_a, coords_b)
        if separation > best_separation:
            best_separation = separation
            best_direction = direction

    target = tuple(anchor[i] + best_direction[i] * distance for i in range(3))
    offset_b = tuple(target[i] - local_b[atom_b][i] for i in range(3))
    coords_b = _translate(local_b, offset_b)
    return coords_a, coords_b


def _build() -> str:
    """构造两链复合物的 PDB 文本。"""
    lines: list[str] = [
        "HEADER    SYNTHETIC TEST COMPLEX",
        "TITLE     TWO CHAIN COMPLEX FOR INTERACTION CONTACT TESTING",
        "REMARK    FIXTURE GENERATED BY tests/make_contact_pdb_fixtures.py",
        "REMARK    EACH RESIDUE PAIR IS PLACED AT A PRECISE CONTACT DISTANCE",
        "REMARK    PAIRS ARE ISOLATED IN SPACE BY 30 A TO AVOID CROSS-TALK",
        "",
    ]

    # 每个配对：(A 残基, A 接触原子, B 残基, B 接触原子, 目标距离)
    # 方向由 `_place_pair` 自动搜索确定 —— 挑出让残基间最小原子间距最大的方向
    pair_specs = [
        ("ARG", 10, "NH1", "ASP", 60, "OD1", 2.9),
        ("LYS", 20, "NZ", "GLU", 70, "OE1", 3.0),
        ("LYS", 25, "NZ", "PHE", 90, "CG", 5.5),
        ("PHE", 30, "CG", "TYR", 80, "CG", 4.5),
        ("ILE", 40, "CG1", "LEU", 100, "CG", 4.2),
        ("GLY", 50, "CA", "GLY", 110, "CA", 3.8),
    ]

    residues: list[tuple[str, str, int, dict]] = []  # (chain, res_name, res_seq, coords)

    for index, (res_a, seq_a, atom_a, res_b, seq_b, atom_b, distance) in enumerate(pair_specs):
        # 每个配对占据独立的 x 区域
        origin_a = (index * PAIR_SPACING, 0.0, 0.0)

        coords_a, coords_b = _place_pair(
            local_a=RESIDUE_TEMPLATES[res_a],
            local_b=RESIDUE_TEMPLATES[res_b],
            atom_a=atom_a,
            atom_b=atom_b,
            distance=distance,
            origin_a=origin_a,
        )

        residues.append(("A", res_a, seq_a, coords_a))
        residues.append(("B", res_b, seq_b, coords_b))

    # 配体：苯环，靠近 TYR80。放在独立区域避免干扰其他配对。
    tyr80 = next(coords for chain, name, seq, coords in residues
                 if chain == "B" and name == "TYR" and seq == 80)
    lig_target = tuple(
        tyr80["CZ"][i] + (0.0, 1.0, 0.0)[i] * 4.0 for i in range(3)
    )
    lig_offset = tuple(lig_target[i] - LIG_LOCAL["C1"][i] for i in range(3))
    residues.append(("L", "LIG", 1, _translate(LIG_LOCAL, lig_offset)))

    # 输出原子
    serial = 1
    for chain, res_name, res_seq, coords in residues:
        record = "HETATM" if res_name == "LIG" else "ATOM  "
        for atom_name, xyz in coords.items():
            lines.append(_atom_line(
                serial, atom_name, res_name, chain, res_seq, xyz,
                _element_of(atom_name), record,
            ))
            serial += 1
        lines.append("TER")
    lines.append("")

    # 氢原子（默认应被忽略）—— 放在 ARG10 的 NH1 旁边
    arg10 = next(coords for chain, name, seq, coords in residues
                  if chain == "A" and name == "ARG" and seq == 10)
    lines.append(_atom_line(
        serial, "H", "ARG", "A", 10,
        (arg10["NH1"][0], arg10["NH1"][1] + 1.0, arg10["NH1"][2]),
        "H",
    ))
    serial += 1
    lines.append("END")
    return "\n".join(lines) + "\n"


def _build_multimodel() -> str:
    """构造含两个 MODEL 的 PDB（用于验证只取第一个 MODEL）。"""
    lines: list[str] = [
        "HEADER    MULTI MODEL TEST",
        "MODEL        1",
    ]
    serial = 1
    for chain, res_name, res_seq, origin in [
        ("A", "ALA", 1, (0.0, 0.0, 0.0)),
        ("B", "GLY", 1, (3.0, 0.0, 0.0)),
    ]:
        local = RESIDUE_TEMPLATES.get(res_name)
        if local is None:
            local = {
                "N": (0.000, 0.000, 0.000), "CA": (1.458, 0.000, 0.000),
                "C": (2.009, 1.420, 0.000), "O": (1.251, 2.390, 0.000),
            }
        for atom_name, xyz in _translate(local, origin).items():
            lines.append(_atom_line(
                serial, atom_name, res_name, chain, res_seq, xyz,
                _element_of(atom_name),
            ))
            serial += 1
    lines.append("ENDMDL")
    lines.append("MODEL        2")
    for chain, res_name, res_seq, origin in [
        ("A", "ALA", 1, (50.0, 0.0, 0.0)),
        ("B", "GLY", 1, (53.0, 0.0, 0.0)),
    ]:
        local = {
            "N": (0.000, 0.000, 0.000), "CA": (1.458, 0.000, 0.000),
            "C": (2.009, 1.420, 0.000), "O": (1.251, 2.390, 0.000),
        }
        for atom_name, xyz in _translate(local, origin).items():
            lines.append(_atom_line(
                serial, atom_name, res_name, chain, res_seq, xyz,
                _element_of(atom_name),
            ))
            serial += 1
    lines.append("ENDMDL")
    lines.append("END")
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------- 自检


def _parse_atoms(pdb_text: str) -> list:
    """从 PDB 文本解析出原子列表（用技能自己的解析器，保证一致）。"""
    sys.path.insert(0, str(REPO_ROOT / "skills" / "pymol-structure-visualization" / "scripts"))
    from contact_geometry import AtomRef

    atoms = []
    serial = 0
    for line in pdb_text.splitlines():
        if line[:6] not in ("ATOM  ", "HETATM"):
            continue
        serial += 1
        atoms.append(AtomRef(
            serial=serial,
            name=line[12:16].strip(),
            res_name=line[17:20].strip(),
            chain_id=line[21:22].strip() or "_",
            res_seq=int(line[22:26].strip() or 0),
            icode=line[26:27].strip(),
            x=float(line[30:38]), y=float(line[38:46]), z=float(line[46:54]),
            element=line[76:78].strip(),
            hetero=(line[:6] == "HETATM"),
        ))
    return atoms


def _self_check(pdb_text: str) -> None:
    """验证夹具确实产生了期望的相互作用类型。

    夹具坏了要立刻知道 —— 而不是让下游测试莫名其妙地失败。
    """
    sys.path.insert(0, str(REPO_ROOT / "skills" / "pymol-structure-visualization" / "scripts"))
    from contact_chemistry import classify_pair_multi
    from contact_geometry import ContactDetector

    atoms = _parse_atoms(pdb_text)

    detector = ContactDetector(max_distance=6.0)
    detector.add_atoms(atoms)

    expected = {
        ("A", "ARG", 10, "B", "ASP", 60): {"Hydrogen_bond", "Salt_bridge"},
        ("A", "LYS", 20, "B", "GLU", 70): {"Hydrogen_bond", "Salt_bridge"},
        ("A", "LYS", 25, "B", "PHE", 90): {"Cation_Pi"},
        ("A", "PHE", 30, "B", "TYR", 80): {"Pi_Pi"},
        ("A", "ILE", 40, "B", "LEU", 100): {"Hydrophobic"},
    }

    problems: list[str] = []
    found: dict[tuple, set[str]] = {}

    chains = sorted({a.chain_id for a in atoms})
    for i, chain_a in enumerate(chains):
        for chain_b in chains[i + 1:]:
            for pair in detector.detect(chain_a, chain_b):
                key = (pair.chain_a, pair.res_name_a, pair.res_seq_a,
                       pair.chain_b, pair.res_name_b, pair.res_seq_b)
                types = set(classify_pair_multi(
                    pair.res_name_a, pair.res_name_b, pair.min_distance))
                found[key] = types

    for key, wanted in expected.items():
        got = found.get(key)
        if got is None:
            problems.append(f"残基对 {key} 没有产生任何接触")
            continue
        missing = wanted - got
        if missing:
            problems.append(
                f"残基对 {key} 缺少类型 {sorted(missing)}（实际 {sorted(got)}）"
            )

    # GLY-GLY 不应产生任何类型
    gly_key = ("A", "GLY", 50, "B", "GLY", 110)
    if found.get(gly_key):
        problems.append(f"GLY-GLY 不该有类型，实际 {sorted(found[gly_key])}")

    # 配体接触应被单独标记
    lig_key = ("B", "TYR", 80, "L", "LIG", 1)
    if found.get(lig_key) != {"Ligand_contact"}:
        problems.append(
            f"配体接触应标记为 Ligand_contact，实际 {sorted(found.get(lig_key, []))}"
        )

    # **原子不得重叠** —— 不同残基之间任何两原子间距 < 2.0 Å 都说明几何坏了
    #
    # 只检查**跨残基**的原子对：残基内部的共价键（C-N 约 1.33 Å）本来就
    # 短得多，用同一阈值检查会把正常结构误判为"重叠"。
    # 预期接触的那一对残基本身距离就 < 2.0 Å，是设计意图，跳过。
    expected_residue_pairs = {
        frozenset({(chain_a, seq_a), (chain_b, seq_b)})
        for chain_a, _, seq_a, chain_b, _, seq_b in expected
    }
    expected_residue_pairs.add(frozenset({("B", 80), ("L", 1)}))

    overlap = []
    for i in range(len(atoms)):
        a = atoms[i]
        for j in range(i + 1, len(atoms)):
            b = atoms[j]
            if a.chain_id == b.chain_id and a.res_seq == b.res_seq:
                continue  # 同一残基内部，跳过
            if frozenset({(a.chain_id, a.res_seq), (b.chain_id, b.res_seq)}) in expected_residue_pairs:
                continue  # 预期接触对，跳过
            d = _dist((a.x, a.y, a.z), (b.x, b.y, b.z))
            if d < MIN_ATOM_SEPARATION:
                overlap.append(
                    f"{a.chain_id}:{a.res_name}{a.res_seq}:{a.name} - "
                    f"{b.chain_id}:{b.res_name}{b.res_seq}:{b.name} = {d:.3f} Å"
                )
    if overlap:
        problems.append(
            f"有 {len(overlap)} 对原子间距 < {MIN_ATOM_SEPARATION} Å（原子重叠）:\n    "
            + "\n    ".join(overlap[:10])
        )

    # **接触距离必须落在合理区间** —— 报告距离不应异常小
    for key in expected:
        chain_a, res_a, seq_a, chain_b, res_b, seq_b = key
        for pair in detector.detect(chain_a, chain_b):
            if pair.res_seq_a != seq_a or pair.res_seq_b != seq_b:
                continue
            if pair.min_distance < 2.0:
                problems.append(
                    f"残基对 {key} 的最小距离 {pair.min_distance:.3f} Å 不合理"
                )

    if problems:
        raise SystemExit("夹具自检失败:\n  " + "\n  ".join(problems))

    print("  夹具自检通过")


def main() -> int:
    """写出全部夹具。"""
    FIXTURES_DIR.mkdir(parents=True, exist_ok=True)

    complex_text = _build()
    _self_check(complex_text)

    files = {
        "two_chain_complex.pdb": complex_text,
        "multimodel.pdb": _build_multimodel(),
    }
    for name, content in files.items():
        path = FIXTURES_DIR / name
        path.write_text(content, encoding="utf-8")
        print(f"已写出 {path.relative_to(REPO_ROOT)}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
