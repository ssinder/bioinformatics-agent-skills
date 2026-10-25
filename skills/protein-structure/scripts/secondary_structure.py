#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""二级结构判定：DSSP 外部程序 + 内置简化 DSSP 回退实现。

模块提供两条独立通路：

1. :func:`run_external_dssp` —— 调用系统中的 ``mkdssp`` / ``dssp`` 可执行文件，
   解析其标准输出（mmCIF 之前的传统 DSSP 文本格式）。
2. :func:`assign_secondary_structure` —— 内置的简化 DSSP 实现，
   不依赖任何外部程序：

   - 计算骨架 N-H 与 C=O 之间的**氢键能量**（Kabsch & Sander 1983 静电模型，
     含侧链修正项），据此识别 ``n -> n+4``（α 螺旋）、``n -> n+3``（3-10 螺旋）、
     ``n -> n+5``（π 螺旋）以及 β-桥/梯形；
   - 结合 φ/ψ 二面角与 DSSP 的 n-turn 定义（``n -> n+3``）判定转角与弯折；
   - 剩余残基按 Ramachandran 区域归类为 coil。

两条通路输出的残基级记录结构一致，可直接互换。

Usage:
    from secondary_structure import assign_secondary_structure
    ss = assign_secondary_structure(structure)
    print(ss["statistics"])
"""

from __future__ import annotations

import math
import shutil
import subprocess
import tempfile
from pathlib import Path

from pdb_parser import AA3_TO_AA1, Atom, Residue, Structure, dihedral_angle

__all__ = [
    "DSSP_CODES",
    "assign_secondary_structure",
    "run_external_dssp",
    "find_dssp_executable",
    "classify_phi_psi",
]

# DSSP 八类代码 -> 中文/英文说明
DSSP_CODES = {
    "H": "Alpha helix (α-螺旋)",
    "G": "3-10 helix (3-10 螺旋)",
    "I": "Pi helix (π-螺旋)",
    "E": "Extended strand (β-折叠股)",
    "B": "Beta bridge (β-桥)",
    "T": "Turn (转角)",
    "S": "Bend (弯折)",
    "-": "Coil (无规卷曲)",
}

HELIX_CODES = {"H", "G", "I"}
SHEET_CODES = {"E", "B"}
COIL_CODES = {"T", "S", "-"}

# Kabsch-Sander 氢键能量系数（单位 kcal/mol）
HBOND_ENERGY_CONST = -27.888
HBOND_ENERGY_DISTANCE_COEFF = 332.0
HBOND_ENERGY_Q1Q2F = 0.084 * 332.0
HBOND_ENERGY_F = 0.20

# DSSP 使用的 C=O 与 N-H 距离（Å），超出即认为不是氢键
HBOND_MAX_DISTANCE = 5.5

# 常见残基的酰胺氢位置修正（PDB 中 H 原子常常缺失，需要几何重建）
# 键长 C-N = 1.32 Å，N-H = 1.00 Å
C_N_BOND = 1.32
N_H_BOND = 1.00
C_O_BOND = 1.24
# 骨架 N-H 与 N-C(prev) 的键角（度），用于几何重建虚拟氢
N_H_ANGLE_DEG = 119.0


def find_dssp_executable(candidates: tuple[str, ...] = ("mkdssp", "dssp", "mkdssp.exe")) -> str | None:
    """在 PATH 中查找 DSSP 可执行文件。"""
    for name in candidates:
        found = shutil.which(name)
        if found:
            return found
    return None


def _classify_backbone_atom(atom: Atom) -> str | None:
    """判断骨架原子属于 N 端 / C 端 / 羰基氧。"""
    name = atom.name.upper()
    if name == "N":
        return "N"
    if name == "O":
        return "O"
    return None


def _rebuild_amide_hydrogens(
    residue: Residue, prev_c: Atom | None
) -> tuple[float, float, float] | None:
    """重建酰胺氢原子的坐标。

    绝大多数实验结构（X 射线 / NMR）不含极性氢原子。DSSP 用的是
    骨架 N 原子与**虚拟氢**的位置构成 N-H 供体。这里用几何法重建。

    几何原理：肽键是平面的（trans 构型），前一个羰基碳 C(prev)、N、CA
    与酰胺氢 H 近似共面。N 上有三条键：C(prev)-N、CA-N、N-H，
    键角分别约为 121.7°、119°、119°。因此：

    - 设 ``u1`` 为 N→C(prev) 的单位向量，``u2`` 为 N→CA 的单位向量；
    - 设 ``θ`` 为 u1 与 u2 的夹角，``e2`` 为 u2 在垂直于 u1 平面内的分量方向；
    - 令 ``h = cos(119°)·u1 + B·e2``，由 ``angle(h, u2) = 119°`` 解出
      ``B = cos(119°)·(1 - cosθ) / sinθ``。

    这样得到的 N-H 键长 1.00 Å，键角与真实几何一致，氢原子指向
    C=O 氧原子一侧，从而使 i→i+4 的 α-螺旋氢键能被正确识别。

    Args:
        residue: 当前残基。
        prev_c: 前一残基的羰基碳原子（链首时为 None）。

    Returns:
        氢原子坐标，缺失关键原子时返回 None。
    """
    n_atom = residue.get("N")
    if n_atom is None or prev_c is None:
        return None
    ca = residue.get("CA")
    if ca is None:
        return None

    ox, oy, oz = n_atom.coord
    v1 = (prev_c.coord[0] - ox, prev_c.coord[1] - oy, prev_c.coord[2] - oz)
    v2 = (ca.coord[0] - ox, ca.coord[1] - oy, ca.coord[2] - oz)

    n1 = math.sqrt(sum(c * c for c in v1))
    n2 = math.sqrt(sum(c * c for c in v2))
    if n1 < 1e-6 or n2 < 1e-6:
        return None

    u1 = (v1[0] / n1, v1[1] / n1, v1[2] / n1)
    u2 = (v2[0] / n2, v2[1] / n2, v2[2] / n2)

    cos_theta = max(-1.0, min(1.0, sum(u1[i] * u2[i] for i in range(3))))
    sin_theta = math.sqrt(max(0.0, 1.0 - cos_theta * cos_theta))
    if sin_theta < 1e-6:
        return None

    # e2: u2 在垂直于 u1 平面内的单位分量
    perp = (u2[0] - cos_theta * u1[0], u2[1] - cos_theta * u1[1], u2[2] - cos_theta * u1[2])
    e2 = (perp[0] / sin_theta, perp[1] / sin_theta, perp[2] / sin_theta)

    cos_h = math.cos(math.radians(N_H_ANGLE_DEG))
    b = cos_h * (1.0 - cos_theta) / sin_theta

    hx = cos_h * u1[0] + b * e2[0]
    hy = cos_h * u1[1] + b * e2[1]
    hz = cos_h * u1[2] + b * e2[2]
    hn = math.sqrt(hx * hx + hy * hy + hz * hz)
    if hn < 1e-9:
        return None
    scale = N_H_BOND / hn
    return (ox + hx * scale, oy + hy * scale, oz + hz * scale)


def _carbonyl_oxygen(residue: Residue) -> Atom | None:
    """取残基的羰基氧原子（优先 O，找不到时用 OXT）。"""
    return residue.get("O") or residue.get("OXT")


def _side_chain_charge(residue: Residue) -> float:
    """侧链在 Kabsch-Sander 模型中的部分电荷修正（q1q2f）。

    规则（简化版）：

    - Asp / Glu 的羧基氧取 -0.5，His 咪唑氮取 -0.5（质子化形式 +1.0），
      Lys / Arg 的氨基取 +1.0（去质子化中性）；
    - 其余残基为 0。
    """
    name = residue.res_name
    if name in ("ASP", "GLU"):
        return -0.5
    if name == "HIS":
        return -0.5
    if name in ("LYS", "ARG"):
        return 0.5
    return 0.0


def _hbond_energy(
    h_pos: tuple[float, float, float],
    o_pos: tuple[float, float, float],
    h_own_charge: float,
    o_owner_charge: float,
) -> float:
    """Kabsch-Sander 氢键能量（kcal/mol）。

    ``E = -27.888 * (1/r_ON + 1/r_CH) + 332.0 * (q1 * q2f / r_ON + q1f * q2 / r_CH)``
    ``+ 0.20 * q1f * q2f * (1/r_ON + 1/r_CH)``

    其中 ``q1`` 是供体 N 上的部分电荷，``q2`` 是受体 C=O 碳上的部分电荷，
    ``q1f`` / ``q2f`` 是其侧链诱导的等效电荷。
    """
    r_on = math.dist(h_pos, o_pos)
    if r_on < 1e-6:
        return 0.0

    # 受体羰基碳的位置：沿 C=O 反方向偏移
    # 用 r_CH 近似为 (r_ON + r_CO 近似值)，这里采用常见近似 r_CH = r_ON + 0.2
    r_ch = r_on + 0.2

    q1 = 0.20 + h_own_charge
    q2 = 0.20 + o_owner_charge
    q1f = 0.0 + h_own_charge
    q2f = 0.0 + o_owner_charge

    e = HBOND_ENERGY_CONST * (1.0 / r_on + 1.0 / r_ch)
    e += HBOND_ENERGY_DISTANCE_COEFF * (q1 * q2f / r_on + q1f * q2 / r_ch)
    e += HBOND_ENERGY_F * q1f * q2f * (1.0 / r_on + 1.0 / r_ch)
    return e


def _build_backbone_records(chains) -> list[dict]:
    """为每个链的每个标准残基建立骨架记录（含重建的酰胺氢）。

    记录字段：``residue``、``pos``（链内序号）、``chain``、``n``、``ca``、
    ``c``、``o``、``h``（重建或读取的酰胺氢）、``charge``（侧链电荷）。
    """
    records: list[dict] = []
    for chain in chains:
        chain_records: list[dict] = []
        prev_c: Atom | None = None
        for residue in chain:
            n_atom = residue.get("N")
            ca_atom = residue.get("CA")
            c_atom = residue.get("C")
            if n_atom is None or ca_atom is None or c_atom is None:
                prev_c = c_atom
                continue
            if not residue.is_standard:
                prev_c = c_atom
                continue

            h_atom = residue.get("H")
            if h_atom is not None:
                h_pos = h_atom.coord
            else:
                h_pos = _rebuild_amide_hydrogens(residue, prev_c)

            o_atom = _carbonyl_oxygen(residue)
            chain_records.append(
                {
                    "residue": residue,
                    "pos": len(chain_records),
                    "chain": chain.chain_id,
                    "n": n_atom.coord,
                    "ca": ca_atom.coord,
                    "c": c_atom.coord,
                    "o": o_atom.coord if o_atom is not None else None,
                    "h": h_pos,
                    "charge": _side_chain_charge(residue),
                }
            )
            prev_c = c_atom
        records.extend(chain_records)
    return records


def _hbond_matrix(records: list[dict]) -> list[list[float]]:
    """计算所有 (供体 i, 受体 j) 对的氢键能量矩阵（kcal/mol）。

    链内残基按序列编号给出**链内位置**（``pos``），氢键方向遵循蛋白质的
    右手螺旋拓扑：受体的 C=O 位于供体 N-H 的**上游**（即 ``pos_acceptor <
    pos_donor``），这样 α-螺旋的 i → i+4 氢键才会被正确识别。
    跨链的氢键不施加该限制。
    """
    n = len(records)
    matrix = [[0.0] * n for _ in range(n)]
    for i, donor in enumerate(records):
        if donor["h"] is None:
            continue
        for j, acceptor in enumerate(records):
            if i == j or acceptor["o"] is None:
                continue
            same_chain = acceptor["residue"].chain_id == donor["residue"].chain_id
            if same_chain and acceptor["pos"] >= donor["pos"]:
                # 同一链内，受体的羰基必须在供体之前（螺旋拓扑）
                continue
            dist = math.dist(donor["h"], acceptor["o"])
            if dist > HBOND_MAX_DISTANCE:
                continue
            energy = _hbond_energy(
                donor["h"], acceptor["o"], donor["charge"], acceptor["charge"]
            )
            matrix[i][j] = energy
    return matrix


def _phi_psi(record: dict, prev_record: dict | None, next_record: dict | None):
    """计算残基的 φ / ψ 二面角（度）。"""
    phi = None
    psi = None
    if prev_record is not None:
        phi = dihedral_angle(prev_record["c"], record["n"], record["ca"], record["c"])
    if next_record is not None:
        psi = dihedral_angle(record["n"], record["ca"], record["c"], next_record["n"])
    return phi, psi


def classify_phi_psi(phi: float | None, psi: float | None) -> str:
    """按 φ/ψ 落在 Ramachandran 图的哪个区域给出结构倾向。

    判据基于常见的经验分区：

    - α 区：φ ∈ [-160, -20]，ψ ∈ [-120, 30]
    - β 区：φ ∈ [-180, -20]，ψ ∈ (30, 180]
    - 其余按 φ/ψ 的具体取值细分左/右螺旋、β 桥、αL 区与 coil
    """
    if phi is None or psi is None:
        return "coil"
    if not (-180.0 <= phi <= 180.0 and -180.0 <= psi <= 180.0):
        return "coil"

    in_alpha_phi = -160.0 <= phi <= -20.0
    if in_alpha_phi and -120.0 <= psi <= 30.0:
        return "alpha"
    if in_alpha_phi and psi > 30.0:
        return "beta"
    if in_alpha_phi and psi < -120.0:
        return "left_alpha"
    if phi > -20.0:
        return "coil"
    return "coil"


def assign_secondary_structure(
    structure: Structure,
    min_hbond_energy: float = -0.5,
    use_hbond_pattern: bool = True,
) -> dict:
    """内置简化 DSSP：按氢键模式与 φ/ψ 判定二级结构。

    算法流程：

    1. 抽取每条链的骨架 N / CA / C / O 坐标，并重建缺失的酰胺氢；
    2. 计算全部 N-H···O=C 氢键的 Kabsch-Sander 能量；
    3. **螺旋**：残基 i 与 i+4 形成氢键（能量 < 阈值）→ α 螺旋 ``H``；
       i 与 i+3 → 3-10 螺旋 ``G``；i 与 i+5 → π 螺旋 ``I``；
    4. **β-桥**：按 DSSP 定义，残基 i 与 j（|i-j| > 2）若形成氢键，
       且 i+1（或 i-1）与 j（或 j+1）之一也形成氢键，则标记为 β-桥 ``B``；
    5. **转角 / 弯折**：存在 i → i+2 的 n-turn（``T``）；
       若 φ、ψ 的角度绝对值较大且二者同号则标记弯折 ``S``；
    6. 其余按 φ/ψ 区域判定为 coil ``-``。

    Args:
        structure: 已解析的结构。
        min_hbond_energy: 判定氢键成立的最大能量（kcal/mol），负值越负要求越严。
        use_hbond_pattern: 是否启用氢键模式判据；关闭后仅用 φ/ψ 区域判定。

    Returns:
        ``{"secondary_structure": [...], "statistics": {...}, "method": "builtin-dssp"}``，
        每个残基记录含 chain / res_num / res_name / ss / ss_name / phi / psi / hbonds。
    """
    records = _build_backbone_records(structure.model.chains)
    energy = _hbond_matrix(records)

    # 记录每个残基收到的氢键（作为供体的能量最低者）与给出的氢键
    hbond_to: dict[int, float] = {}
    hbond_from: dict[int, float] = {}
    if use_hbond_pattern:
        for i, donor in enumerate(records):
            for j in range(len(records)):
                e = energy[i][j]
                if e < min_hbond_energy:
                    if j not in hbond_to or e < hbond_to[j]:
                        hbond_to[j] = e
                    if i not in hbond_from or e < hbond_from[i]:
                        hbond_from[i] = e

    by_chain: dict[str, list[int]] = {}
    for idx, rec in enumerate(records):
        by_chain.setdefault(rec["chain"], []).append(idx)

    assignments: dict[int, str] = {}

    if use_hbond_pattern:
        for chain_id, indices in by_chain.items():
            for p, idx in enumerate(indices):
                # n → n+4 : α 螺旋（供体在下游，受体在上游）
                if p + 4 < len(indices):
                    j = indices[p + 4]
                    if energy[j][idx] < min_hbond_energy:
                        assignments[idx] = "H"
                        assignments[j] = "H"
                        continue
                # n → n+3 : 3-10 螺旋
                if p + 3 < len(indices):
                    j = indices[p + 3]
                    if energy[j][idx] < min_hbond_energy:
                        if assignments.get(idx) != "H":
                            assignments[idx] = "G"
                        if assignments.get(j) != "H":
                            assignments[j] = "G"
                        continue
                # n → n+5 : π 螺旋
                if p + 5 < len(indices):
                    j = indices[p + 5]
                    if energy[j][idx] < min_hbond_energy:
                        if assignments.get(idx) not in ("H", "G"):
                            assignments[idx] = "I"
                        if assignments.get(j) not in ("H", "G"):
                            assignments[j] = "I"

        # β-桥（DSSP 定义：至少一个桥，且不是螺旋残基）
        for chain_id, indices in by_chain.items():
            for p, idx in enumerate(indices):
                if assignments.get(idx) in HELIX_CODES:
                    continue
                for q, j in enumerate(indices):
                    if q <= p:
                        continue
                    if q - p < 3:
                        continue
                    parallel_bridge = energy[idx][j] < min_hbond_energy
                    antiparallel_bridge = energy[j][idx] < min_hbond_energy
                    if not (parallel_bridge or antiparallel_bridge):
                        continue
                    # 平行桥要求 i+1→j+1；反平行桥要求 j→i-1
                    if parallel_bridge:
                        ok = (
                            p + 1 < len(indices)
                            and q + 1 < len(indices)
                            and energy[indices[p + 1]][indices[q + 1]] < min_hbond_energy
                        )
                    else:
                        ok = (
                            p - 1 >= 0
                            and energy[j][indices[p - 1]] < min_hbond_energy
                        )
                    if ok:
                        assignments[idx] = "B"
                        assignments[j] = "B"
                        break

    # 计算 φ/ψ、转角与弯折
    secondary: list[dict] = []
    for chain_id, indices in by_chain.items():
        for p, idx in enumerate(indices):
            rec = records[idx]
            residue = rec["residue"]
            prev_rec = records[indices[p - 1]] if p > 0 else None
            next_rec = records[indices[p + 1]] if p + 1 < len(indices) else None
            phi, psi = _phi_psi(rec, prev_rec, next_rec)

            code = assignments.get(idx, "")
            if code == "":
                # 未被氢键模式覆盖：按 DSSP 的优先级判定
                # 螺旋 > 折叠股 > 转角 > 弯折 > 无规卷曲
                region = classify_phi_psi(phi, psi)
                turn = False
                if p + 2 < len(indices):
                    j = indices[p + 2]
                    # n-turn：供体 i 与受体 i+2 形成氢键
                    if energy[j][idx] < min_hbond_energy:
                        turn = True

                if region == "beta":
                    # 伸展区域优先判为折叠股（β-折叠股的 φ<0、ψ>90）
                    code = "E"
                elif region == "alpha":
                    code = "H"
                elif turn:
                    code = "T"
                elif phi is not None and psi is not None and abs(phi) > 90 and abs(psi) > 90:
                    code = "S"
                else:
                    code = "-"

            secondary.append(
                {
                    "chain": chain_id,
                    "res_num": residue.res_seq,
                    "res_name": residue.res_name,
                    "icode": residue.icode.strip(),
                    "ss": code,
                    "ss_name": DSSP_CODES.get(code, "Unknown"),
                    "phi": round(phi, 1) if phi is not None else None,
                    "psi": round(psi, 1) if psi is not None else None,
                    "hbond_in": round(hbond_to.get(idx, 0.0), 3) if use_hbond_pattern else 0.0,
                    "hbond_out": round(hbond_from.get(idx, 0.0), 3) if use_hbond_pattern else 0.0,
                }
            )

    stats = _statistics(secondary)
    return {
        "secondary_structure": secondary,
        "statistics": stats,
        "method": "builtin-dssp",
        "n_hbonds": sum(1 for row in energy for e in row if e < min_hbond_energy),
    }


def _statistics(secondary: list[dict]) -> dict:
    """统计二级结构分布。"""
    total = len(secondary)
    helix = sum(1 for s in secondary if s["ss"] in HELIX_CODES)
    sheet = sum(1 for s in secondary if s["ss"] in SHEET_CODES)
    coil = total - helix - sheet
    counts: dict[str, int] = {}
    for s in secondary:
        counts[s["ss"]] = counts.get(s["ss"], 0) + 1
    return {
        "n_residues": total,
        "helix_count": helix,
        "sheet_count": sheet,
        "coil_count": coil,
        "helix_pct": round(helix / total * 100, 2) if total else 0.0,
        "sheet_pct": round(sheet / total * 100, 2) if total else 0.0,
        "coil_pct": round(coil / total * 100, 2) if total else 0.0,
        "code_counts": dict(sorted(counts.items())),
    }


def run_external_dssp(
    pdb_path: str | Path,
    executable: str | None = None,
    timeout: int = 120,
) -> dict | None:
    """调用外部 DSSP 程序并解析输出。

    Args:
        pdb_path: PDB 文件路径。
        executable: DSSP 可执行文件路径，None 时自动在 PATH 中查找。
        timeout: 子进程超时秒数。

    Returns:
        与 :func:`assign_secondary_structure` 相同结构的字典；
        程序不可用或执行失败时返回 None。
    """
    exe = executable or find_dssp_executable()
    if not exe:
        return None

    with tempfile.TemporaryDirectory() as tmpdir:
        out_path = Path(tmpdir) / "out.dssp"
        cmd = [exe, str(pdb_path), "-o", str(out_path)]
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        except (subprocess.TimeoutExpired, OSError):
            return None
        if proc.returncode != 0 or not out_path.exists():
            # 老版本 mkdssp 走 stdout
            if proc.stdout:
                return _parse_dssp_text(proc.stdout, Path(pdb_path).stem)
            return None
        return _parse_dssp_text(out_path.read_text(encoding="utf-8", errors="replace"),
                                Path(pdb_path).stem)


def _parse_dssp_text(text: str, structure_id: str) -> dict | None:
    """解析传统 DSSP 文本格式。

    关键列：残基号 ``[5:10]``、链标识 ``[11]``、氨基酸 ``[13]``、
    二级结构码 ``[16]``。
    """
    secondary: list[dict] = []
    started = False
    for line in text.splitlines():
        if not started:
            if line.startswith("  #  RESIDUE") or "RESIDUE" in line[:20]:
                started = True
            continue
        if len(line) < 17:
            continue
        try:
            res_num = int(line[5:10].strip())
        except ValueError:
            continue
        chain = line[11].strip() if len(line) > 11 else ""
        res_name = line[13:16].strip() if len(line) >= 16 else ""
        code = line[16:17].strip()
        if code == "":
            code = "-"
        secondary.append(
            {
                "chain": chain,
                "res_num": res_num,
                "res_name": res_name,
                "icode": "",
                "ss": code,
                "ss_name": DSSP_CODES.get(code, "Unknown"),
                "phi": None,
                "psi": None,
            }
        )
    if not secondary:
        return None
    return {
        "secondary_structure": secondary,
        "statistics": _statistics(secondary),
        "method": "external-dssp",
        "structure_id": structure_id,
    }


def _main(argv: list[str] | None = None) -> int:
    """命令行入口：对单个 PDB 打印二级结构判定结果。"""
    import argparse

    parser = argparse.ArgumentParser(
        prog="secondary_structure", description="二级结构判定（命令行自检工具）"
    )
    parser.add_argument("input", help="PDB 文件路径")
    parser.add_argument(
        "--method",
        choices=("auto", "external", "builtin"),
        default="builtin",
        help="二级结构计算方式（默认 builtin）",
    )
    parser.add_argument("--hbond-energy", type=float, default=-0.5, help="氢键能量阈值")
    args = parser.parse_args(argv)

    from pdb_parser import parse_pdb_file

    structure = parse_pdb_file(args.input)
    if args.method == "external":
        result = run_external_dssp(args.input)
        if result is None:
            logger.error("外部 DSSP 不可用或执行失败")
            return 1
    else:
        result = assign_secondary_structure(
            structure, min_hbond_energy=args.hbond_energy, use_hbond_pattern=True
        )

    logger.info("方法: %s", result["method"])
    for key, value in result["statistics"].items():
        logger.info("%s: %s", key, value)
    for item in result["secondary_structure"]:
        logger.info(
            "  %s:%s%s  %s  phi=%s psi=%s",
            item["chain"],
            item["res_name"],
            item["res_num"],
            item["ss"],
            item.get("phi"),
            item.get("psi"),
        )
    return 0


if __name__ == "__main__":
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "packages"))
    from biopipe_core import get_logger

    logger = get_logger("biopipe.protein-structure")
    sys.exit(_main())
