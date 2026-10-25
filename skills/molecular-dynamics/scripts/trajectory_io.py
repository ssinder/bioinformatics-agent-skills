"""轨迹读取：内置纯 numpy 的 DCD / PDB 轨迹解析器，MDAnalysis 可用时优先使用。

本模块提供统一的 ``load_trajectory`` 接口，返回一个 :class:`Trajectory` 对象：

* ``coords``   -- ``(n_frames, n_atoms, 3)`` float64 坐标数组
* ``topology`` -- 拓扑对象，提供 ``n_atoms`` / ``positions`` / ``residues`` 等
* ``times``    -- ``(n_frames,)`` 时间序列（ps）

DCD 为定长记录格式：文件开头是 84 字节 header 记录，随后是可选的
"每帧含 3 个 4 字节整数记录头 + 原子坐标" 的定长记录块。详见
``references/methodology.md``。
"""

from __future__ import annotations

import struct
import sys
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

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



__all__ = [
    "Trajectory",
    "SimpleTopology",
    "load_trajectory",
    "read_dcd",
    "read_pdb_trajectory",
    "read_xyz_trajectory",
    "dcd_header_info",
]


# ---------------------------------------------------------------------------
# 轻量拓扑容器
# ---------------------------------------------------------------------------


@dataclass
class Residue:
    """残基信息。"""

    resid: int
    resname: str
    chain: str
    index: int = 0
    atoms: list[int] = field(default_factory=list)


@dataclass
class Atom:
    """原子信息。"""

    index: int
    name: str
    resname: str
    resid: int
    chain: str
    element: str = ""


class SimpleTopology:
    """不依赖 MDAnalysis 的最小拓扑描述。"""

    def __init__(self, atoms: Sequence[Atom], residues: Sequence[Residue]):
        self.atoms = list(atoms)
        self.residues = list(residues)
        self._residue_index = {r.index: r for r in self.residues}

    @property
    def n_atoms(self) -> int:
        return len(self.atoms)

    @property
    def n_residues(self) -> int:
        return len(self.residues)

    @property
    def n_chains(self) -> int:
        return len({r.chain for r in self.residues})

    def chain_ids(self) -> list[str]:
        seen: list[str] = []
        for r in self.residues:
            if r.chain not in seen:
                seen.append(r.chain)
        return seen

    def residue_atoms(self, res_index: int) -> list[int]:
        res = self._residue_index.get(res_index)
        return list(res.atoms) if res else []

    def atom_mask_by_name(self, names: Iterable[str]) -> np.ndarray:
        """返回布尔掩码，标记原子名在给定集合中的原子。"""
        wanted = {n.strip().upper() for n in names}
        return np.array([a.name.strip().upper() in wanted for a in self.atoms], dtype=bool)


@dataclass
class Trajectory:
    """轨迹数据容器。"""

    coords: np.ndarray  # (n_frames, n_atoms, 3)
    topology: SimpleTopology
    times: np.ndarray  # (n_frames,)
    source: str = ""
    reader: str = "builtin"

    @property
    def n_frames(self) -> int:
        return int(self.coords.shape[0])

    @property
    def n_atoms(self) -> int:
        return int(self.coords.shape[1])

    def __repr__(self) -> str:  # pragma: no cover - 调试用
        return (
            f"Trajectory(frames={self.n_frames}, atoms={self.n_atoms}, "
            f"chains={self.topology.n_chains}, reader={self.reader})"
        )


# ---------------------------------------------------------------------------
# DCD 解析（纯 numpy / 标准库）
# ---------------------------------------------------------------------------


def _read_record_header(fh) -> int:
    """读取一个 Fortran 风格记录的 4 字节长度标记。"""
    raw = fh.read(4)
    if len(raw) < 4:
        raise EOFError("DCD 记录头意外结束")
    (length,) = struct.unpack("<i", raw)
    return length


def _read_record_body(fh, length: int) -> bytes:
    """读取已知长度的记录体，并校验尾部长度标记。"""
    raw = fh.read(length)
    if len(raw) < length:
        raise EOFError("DCD 记录数据意外结束")
    tail = fh.read(4)
    if len(tail) < 4:
        raise EOFError("DCD 记录尾标记意外结束")
    (tail_length,) = struct.unpack("<i", tail)
    if tail_length != length:
        raise ValueError(
            f"DCD 记录长度不匹配: 头 {length} != 尾 {tail_length}（文件可能已损坏）"
        )
    return raw


def _read_record(fh, dtype: str, count: int) -> np.ndarray:
    """读取一个 Fortran 风格记录，返回 ``count`` 个 ``dtype`` 元素。"""
    length = _read_record_header(fh)
    expected = count * struct.calcsize(dtype)
    if length != expected:
        # 部分写盘工具会写错长度标记，这里按期望长度兜底读取。
        raw = fh.read(expected)
        if len(raw) < expected:
            raise EOFError("DCD 记录数据意外结束")
        fh.read(4)  # 跳过尾部长度标记
        return np.frombuffer(raw, dtype=dtype, count=count)
    raw = _read_record_body(fh, length)
    return np.frombuffer(raw, dtype=dtype, count=count)


def _read_record_raw(fh) -> bytes:
    """读取一个 Fortran 风格记录，返回原始字节。"""
    length = _read_record_header(fh)
    if length < 0:
        raise ValueError(f"DCD 记录长度为负: {length}")
    raw = fh.read(length)
    if len(raw) < length:
        raise EOFError("DCD 记录数据意外结束")
    tail = fh.read(4)
    if len(tail) < 4:
        raise EOFError("DCD 记录尾标记意外结束")
    (tail_length,) = struct.unpack("<i", tail)
    if tail_length != length:
        raise ValueError(
            f"DCD 记录长度不匹配: 头 {length} != 尾 {tail_length}（文件可能已损坏）"
        )
    return raw


def _read_dcd_header(fh) -> tuple[np.ndarray, int]:
    """读取 DCD 头记录，返回 ``(icntrl, header_len)``。

    DCD 头是一条 Fortran 记录，其**记录体**布局为::

        4 字节魔数 'CORD'
        80 字节 icntrl（20 个 int32）
        [4 字节尾长度标记]

    记录体总长为 84 字节（CHARMM）或 84+80=164 字节（带 unitcell 的版本）。
    注意记录体**本身不含**开头的长度标记——长度标记在记录体之外，
    因此不能像普通 Fortran 记录那样用 ``_read_record`` 读取，否则会把
    'CORD' 魔数误当作第一个 int32。

    Args:
        fh: 已打开的二进制文件对象（定位在文件开头）。

    Returns:
        ``(icntrl, header_len)``，其中 icntrl 为 20 个 int32 的数组。

    Raises:
        ValueError: 魔数不是 'CORD'，或记录长度不合法。
    """
    record_len = _read_record_header(fh)
    if record_len not in (84, 164):
        raise ValueError(f"不是有效的 DCD 文件（头记录长度 {record_len}，期望 84 或 164）")

    body = fh.read(record_len)
    if len(body) < record_len:
        raise EOFError("DCD 头记录意外结束")

    magic = body[:4]
    if magic != b"CORD":
        raise ValueError(f"DCD 魔数错误：{magic!r}，期望 b'CORD'")

    icntrl = np.frombuffer(body[4:84], dtype="<i4", count=20).copy()

    # 84 字节版本后面跟着 4 字节尾长度标记；164 字节版本则把 unitcell
    # 字段直接放在 icntrl 之后，不再有尾标记。
    if record_len == 84:
        tail = fh.read(4)
        if len(tail) < 4:
            raise EOFError("DCD 头记录尾标记意外结束")
        (tail_len,) = struct.unpack("<i", tail)
        if tail_len != 84:
            raise ValueError(f"DCD 头记录长度不匹配：头 {84} != 尾 {tail_len}")

    return icntrl, record_len


def _sniff_natoms(fh) -> int:
    """在头记录之后嗅探 natoms 记录。

    真实布局（CHARMM / NAMD）::

        [头记录 84 或 164 字节]
        [标题记录]     长度通常是 80 的倍数（244 等），NAMD 会写 4+80+pad
        [natoms 记录]  长度 4，值为原子数
        [自由原子记录] 可选
        [坐标记录...]

    因此不能假定 natoms 记录紧跟头部，必须跳过任意数量的辅助记录，
    直到遇到长度为 4 且取值合理的记录。

    Args:
        fh: 已打开的二进制文件对象。

    Returns:
        原子数；未找到时返回 0。
    """
    for _ in range(64):
        try:
            rec_len = _read_record_header(fh)
        except EOFError:
            return 0
        if rec_len < 0:
            return 0
        if rec_len == 4:
            raw = fh.read(4)
            if len(raw) < 4:
                return 0
            (value,) = struct.unpack("<i", raw)
            fh.read(4)  # 尾长度标记
            if 0 < value <= 10_000_000:
                return value
            continue
        if rec_len == 0:
            fh.read(4)  # 空记录，只有尾标记
            continue
        if rec_len < (1 << 20) and rec_len % 4 == 0:
            fh.read(rec_len)
            fh.read(4)
            continue
        return 0
    return 0


def dcd_header_info(path: str | Path) -> dict:
    """解析 DCD header，返回关键字段（不读坐标）。

    Args:
        path: DCD 文件路径。

    Returns:
        含 ``n_frames`` / ``n_atoms`` / ``nsavc`` / ``delta`` 等字段的字典。
    """
    path = Path(path)
    with path.open("rb") as fh:
        icntrl, _ = _read_dcd_header(fh)
        natoms = _sniff_natoms(fh)

    charmm = int(icntrl[10])

    # DELTA 以 float32 存放在 icntrl 的第 9、10 个槽
    delta = float(np.frombuffer(icntrl[7].tobytes() + icntrl[8].tobytes(), dtype="<f4")[0])
    return {
        "n_frames": int(icntrl[0]),
        "istart": int(icntrl[1]),
        "nsavc": int(icntrl[2]),
        "nsteps": int(icntrl[3]),
        "ndegf": int(icntrl[5]),
        "delta": delta,
        "charmm": charmm,
        "unitcell": int(icntrl[19]),
        "n_atoms": natoms,
    }


def read_dcd(path: str | Path, read_frames: int | None = None) -> tuple[np.ndarray, dict]:
    """读取 DCD 轨迹（纯 Python + numpy 实现）。

    DCD 是 Fortran 风格的定长记录流。文件布局（CHARMM / NAMD 规范）::

        [header 记录]     84 字节 = 21 个 int32
        [辅助记录...]      标题记录（长度是 80 的倍数）等
        [natoms 记录]     1 个 int32
        [自由原子记录]    可选，长度 = nfixed * 4
        重复 n_frames 次:
            [X 记录][Y 记录][Z 记录]    每个 natoms * 4 字节
            [unit cell 记录]            可选，48 字节（6 个 float64）

    不同发行版（NAMD / CHARMM / X-PLOR / VMD）对 header 内字段的偏移并不统一，
    因此本实现**不依赖 header 字段定位 natoms**，而是通过记录长度自描述地嗅探：
    长度为 4 的记录即 natoms 记录，长度为 ``natoms*4`` 或 ``natoms*3*4`` 的记录
    即坐标记录，其余记录按辅助记录跳过。

    Args:
        path: DCD 文件路径。
        read_frames: 只读前 N 帧（None 表示全部）。

    Returns:
        ``(coords, info)``，其中 ``coords`` 形状为 ``(n_frames, n_atoms, 3)``。

    Raises:
        ValueError: 文件格式不合法。
        EOFError: 文件提前结束。
    """
    path = Path(path)
    with path.open("rb") as fh:
        icntrl, _ = _read_dcd_header(fh)

        n_frames_header = int(icntrl[0])
        charmm = int(icntrl[10])

        # --- 嗅探 natoms ----------------------------------------------------
        natoms = _sniff_natoms(fh)
        if natoms <= 0:
            raise ValueError("DCD 文件中未找到合法的 natoms 记录")

        total = n_frames_header if n_frames_header > 0 else 0
        limit = total if read_frames is None else min(total or read_frames, read_frames)

        coords = np.zeros((max(limit, 0), natoms, 3), dtype=np.float32)
        rec_single = natoms * 3 * 4
        rec_axis = natoms * 4
        rec_unitcell = 48  # 6 个 float64 的晶胞记录

        # 真实 DCD（CHARMM/NAMD）每帧的记录序列为：
        #   [unitcell 48B] [X  natoms*4] [Y  natoms*4] [Z  natoms*4]
        # 或（无晶胞时）
        #   [X natoms*4] [Y natoms*4] [Z natoms*4]
        # 另有少数发行版把三个轴合并成一条 [X..Y..Z] 记录。
        # 记录本身没有标签，只能靠长度与出现顺序判定，故用状态机累积。
        #
        # icntrl 字段布局（CHARMM 规范）：
        #   [0] NSET 帧数  [1] ISTART  [2] NSAVC  [3] NSTEP
        #   [7]/[8] DELTA 的 float32 高/低位  [10] CHARM 版本标志
        #   [19] CHARMM 版本号
        # 是否带晶胞记录由 icntrl[19] 的 CHARMM 版本判定：CHARMM 写出的
        # 轨迹每帧都带 48 字节晶胞记录。
        unitcell_flag = int(icntrl[19]) >= 20 or charmm != 0

        n_read = 0
        pending_axes: list[np.ndarray] = []
        saw_unitcell_in_frame = False

        while True:
            if limit and n_read >= limit:
                break
            try:
                length = _read_record_header(fh)
            except EOFError:
                break
            if length < 0 or length > (1 << 26):
                break

            if length == rec_axis:
                # 轴记录：X / Y / Z。
                # 特殊情形：natoms == 12 时 rec_axis == 48 == rec_unitcell，
                # 此时必须靠 icntrl[19] 的晶胞标志与“帧内位置”判定是坐标还是晶胞。
                if (
                    length == rec_unitcell
                    and not pending_axes
                    and unitcell_flag
                    and not saw_unitcell_in_frame
                ):
                    _read_record_body(fh, length)
                    saw_unitcell_in_frame = True
                    continue
                block = np.frombuffer(
                    _read_record_body(fh, length), dtype="<f4", count=natoms
                )
                pending_axes.append(block)
                if len(pending_axes) == 3:
                    if n_read < coords.shape[0]:
                        coords[n_read, :, 0] = pending_axes[0]
                        coords[n_read, :, 1] = pending_axes[1]
                        coords[n_read, :, 2] = pending_axes[2]
                    pending_axes = []
                    n_read += 1
                    saw_unitcell_in_frame = False
                continue

            if length == rec_single:
                # 合并记录：一帧一条，内部按 [x..., y..., z...] 排列
                block = np.frombuffer(
                    _read_record_body(fh, length), dtype="<f4", count=natoms * 3
                )
                if n_read < coords.shape[0]:
                    coords[n_read, :, 0] = block[:natoms]
                    coords[n_read, :, 1] = block[natoms : 2 * natoms]
                    coords[n_read, :, 2] = block[2 * natoms :]
                n_read += 1
                pending_axes = []
                saw_unitcell_in_frame = False
                continue

            if length == rec_unitcell and not pending_axes and not saw_unitcell_in_frame:
                # 晶胞记录：每帧的第一条，6 个 float64，直接跳过
                _read_record_body(fh, length)
                saw_unitcell_in_frame = True
                continue

            if 0 < length < (1 << 24) and length % 4 == 0:
                # 其它辅助记录（标题 / 自由原子列表），跳过
                _read_record_body(fh, length)
                continue

            break

    if n_read == 0:
        raise ValueError("DCD 文件中未读到任何坐标帧")

    coords = coords[:n_read].astype(np.float64)
    delta = float(np.frombuffer(icntrl[7].tobytes() + icntrl[8].tobytes(), dtype="<f4")[0])
    info = {
        "n_frames": n_read,
        "n_atoms": natoms,
        "charmm": charmm,
        "unitcell": int(icntrl[19]),
        "nsavc": int(icntrl[2]),
        "delta": delta,
        "source": str(path),
    }
    return coords, info


# ---------------------------------------------------------------------------
# PDB / XYZ 轨迹
# ---------------------------------------------------------------------------

_THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "HSD": "H", "HSE": "H", "HSP": "H", "CYX": "C",
}


def _pdb_topology_and_frame(
    lines: Sequence[str],
) -> tuple[SimpleTopology, np.ndarray]:
    """从一组 PDB 行解析拓扑与坐标。"""
    atoms: list[Atom] = []
    coords: list[tuple[float, float, float]] = []
    residue_map: dict[tuple[str, int, str], Residue] = {}
    residues: list[Residue] = []

    for line in lines:
        record = line[0:6].strip()
        if record not in ("ATOM", "HETATM"):
            continue
        altloc = line[16]
        if altloc not in (" ", "A", "1"):
            continue
        try:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except ValueError:
            continue
        name = line[12:16].strip()
        resname = line[17:20].strip()
        chain = line[21].strip() or "A"
        try:
            resid = int(line[22:26])
        except ValueError:
            continue
        icode = line[26].strip()
        element = line[76:78].strip()
        if not element:
            element = "".join(ch for ch in name if ch.isalpha())[:2].upper()

        idx = len(atoms)
        atoms.append(Atom(idx, name, resname, resid, chain, element))
        coords.append((x, y, z))

        key = (chain, resid, icode)
        res = residue_map.get(key)
        if res is None:
            res = Residue(resid=resid, resname=resname, chain=chain, index=len(residues))
            residue_map[key] = res
            residues.append(res)
        res.atoms.append(idx)

    if not atoms:
        raise ValueError("PDB 轨迹中没有解析到任何原子")
    topology = SimpleTopology(atoms, residues)
    return topology, np.asarray(coords, dtype=np.float64)


def read_pdb_trajectory(path: str | Path) -> Trajectory:
    """读取多模型（MODEL/ENDMDL）PDB 轨迹。

    单模型 PDB 也可读取，此时轨迹长度为 1 帧。
    """
    path = Path(path)
    text = path.read_text(errors="ignore").splitlines()
    models: list[list[str]] = []
    current: list[str] = []
    in_model = False
    for line in text:
        record = line[0:6].strip()
        if record == "MODEL":
            in_model = True
            current = []
            continue
        if record == "ENDMDL":
            if current:
                models.append(current)
            current = []
            in_model = False
            continue
        if in_model or record in ("ATOM", "HETATM", "TER", "ANISOU"):
            current.append(line)
    if current and not models:
        models.append(current)
    if not models:
        raise ValueError("PDB 文件中没有找到任何原子记录")

    topologies = []
    frames = []
    for model_lines in models:
        topo, xyz = _pdb_topology_and_frame(model_lines)
        topologies.append(topo)
        frames.append(xyz)

    n_atoms = topologies[0].n_atoms
    for topo in topologies[1:]:
        if topo.n_atoms != n_atoms:
            raise ValueError("多模型 PDB 的各帧原子数不一致，无法堆叠为轨迹")

    coords = np.stack(frames, axis=0)
    times = np.arange(coords.shape[0], dtype=np.float64)
    return Trajectory(
        coords=coords,
        topology=topologies[0],
        times=times,
        source=str(path),
        reader="builtin-pdb",
    )


def read_xyz_trajectory(path: str | Path) -> Trajectory:
    """读取多帧 XYZ 轨迹（无拓扑，元素按顺序猜测）。"""
    path = Path(path)
    lines = path.read_text(errors="ignore").splitlines()
    frames: list[np.ndarray] = []
    elements: list[str] = []
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if not line:
            i += 1
            continue
        try:
            natoms = int(line.split()[0])
        except (ValueError, IndexError):
            i += 1
            continue
        i += 1
        i += 1
        xyz = np.zeros((natoms, 3), dtype=np.float64)
        for k in range(natoms):
            parts = lines[i + k].split()
            if len(parts) >= 4:
                xyz[k] = [float(parts[1]), float(parts[2]), float(parts[3])]
                if not elements:
                    elements.append(parts[0])
        i += natoms
        frames.append(xyz)

    if not frames:
        raise ValueError("XYZ 文件中没有解析到任何帧")

    n_atoms = frames[0].shape[0]
    frames = [f for f in frames if f.shape[0] == n_atoms]
    coords = np.stack(frames, axis=0)
    atoms = [
        Atom(i, elements[i] if i < len(elements) else "X", "UNK", 1, "A", elements[i] if i < len(elements) else "X")
        for i in range(n_atoms)
    ]
    res = Residue(resid=1, resname="UNK", chain="A", index=0, atoms=list(range(n_atoms)))
    topology = SimpleTopology(atoms, [res])
    times = np.arange(coords.shape[0], dtype=np.float64)
    return Trajectory(coords=coords, topology=topology, times=times, source=str(path), reader="builtin-xyz")


# ---------------------------------------------------------------------------
# 统一入口
# ---------------------------------------------------------------------------


def _load_with_mdanalysis(path: Path, start: int, stop: int, stride: int) -> Trajectory | None:
    """尝试用 MDAnalysis 读取；不可用时返回 None。"""
    try:
        import MDAnalysis as mda  # type: ignore
    except Exception:  # noqa: BLE001
        return None
    try:
        universe = mda.Universe(str(path))
        if universe.atoms.n_atoms == 0:
            return None
        frames = []
        for ts in universe.trajectory[start:stop:stride]:
            frames.append(ts.positions.astype(np.float64))
        if not frames:
            return None
        coords = np.stack(frames, axis=0)
        atoms = []
        for atom in universe.atoms:
            atoms.append(
                Atom(
                    index=int(atom.ix),
                    name=str(atom.name),
                    resname=str(atom.resname),
                    resid=int(atom.resid),
                    chain=str(atom.chainID),
                    element=str(atom.element),
                )
            )
        residues = []
        for res in universe.residues:
            residues.append(
                Residue(
                    resid=int(res.resid),
                    resname=str(res.resname),
                    chain=str(res.chainID),
                    index=int(res.ix),
                    atoms=[int(a.ix) for a in res.atoms],
                )
            )
        topology = SimpleTopology(atoms, residues)
        times = np.asarray(universe.trajectory[start:stop:stride].ts.time, dtype=np.float64)
        if times.shape[0] != coords.shape[0]:
            times = np.arange(coords.shape[0], dtype=np.float64)
        return Trajectory(
            coords=coords, topology=topology, times=times, source=str(path), reader="mdanalysis"
        )
    except Exception:  # noqa: BLE001
        return None


def load_trajectory(
    path: str | Path,
    start: int = 0,
    stop: int | None = None,
    stride: int = 1,
    prefer_mdanalysis: bool = True,
) -> Trajectory:
    """统一轨迹加载入口。

    优先使用 MDAnalysis（若已安装且可解析），否则回退到内置解析器。

    Args:
        path: 轨迹文件（``.dcd`` / ``.pdb`` / ``.xyz`` / 其它 MDAnalysis 支持的格式）。
        start: 起始帧。
        stop: 结束帧（不含）。
        stride: 帧步长。
        prefer_mdanalysis: 是否优先使用 MDAnalysis。

    Returns:
        :class:`Trajectory` 实例。

    Raises:
        FileNotFoundError: 文件不存在。
        SkillError-like ValueError: 无法解析。
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"轨迹文件不存在: {path}")

    if prefer_mdanalysis:
        traj = _load_with_mdanalysis(path, start, stop or -1, stride)
        if traj is not None:
            return traj

    suffix = path.suffix.lower()
    if suffix == ".dcd":
        coords, info = read_dcd(path)
        traj = _trajectory_from_coords(coords, info, path)
    elif suffix in (".pdb", ".ent"):
        traj = read_pdb_trajectory(path)
    elif suffix == ".xyz":
        traj = read_xyz_trajectory(path)
    else:
        # 未知扩展名：尝试 PDB，再尝试 DCD
        try:
            traj = read_pdb_trajectory(path)
        except Exception:  # noqa: BLE001
            coords, info = read_dcd(path)
            traj = _trajectory_from_coords(coords, info, path)

    if start or stop or stride != 1:
        end = stop if stop is not None else traj.n_frames
        idx = np.arange(start, min(end, traj.n_frames), stride)
        traj = Trajectory(
            coords=traj.coords[idx],
            topology=traj.topology,
            times=traj.times[idx],
            source=traj.source,
            reader=traj.reader,
        )
    return traj


def _trajectory_from_coords(coords: np.ndarray, info: dict, path: Path) -> Trajectory:
    """为没有拓扑的坐标构造占位拓扑。"""
    n_atoms = coords.shape[1]
    atoms = [Atom(i, "C", "UNK", 1, "A", "C") for i in range(n_atoms)]
    res = Residue(resid=1, resname="UNK", chain="A", index=0, atoms=list(range(n_atoms)))
    topology = SimpleTopology(atoms, [res])
    nsavc = int(info.get("nsavc", 1) or 1)
    delta = float(info.get("delta", 1.0) or 1.0)
    times = np.arange(coords.shape[0], dtype=np.float64) * nsavc * delta
    return Trajectory(
        coords=coords, topology=topology, times=times, source=str(path), reader="builtin-dcd"
    )


# ---------------------------------------------------------------------------
# 便捷工具
# ---------------------------------------------------------------------------


def selection_by_resname(topology: SimpleTopology, resnames: Sequence[str]) -> np.ndarray:
    """返回属于指定残基名的残基索引数组。"""
    wanted = {r.upper() for r in resnames}
    return np.array(
        [res.index for res in topology.residues if res.resname.upper() in wanted],
        dtype=int,
    )


def backbone_atom_indices(topology: SimpleTopology) -> np.ndarray:
    """返回主链原子（N、CA、C、O）索引数组。"""
    names = {"N", "CA", "C", "O"}
    mask = topology.atom_mask_by_name(names)
    return np.flatnonzero(mask)


def extract_backbone(topology: SimpleTopology) -> dict:
    """从拓扑中提取逐残基的主链原子索引。

    Returns:
        字典，含 ``residue_index``（残基在拓扑中的序号）、``n`` / ``ca`` / ``c`` / ``o``
        （每个元素是对应原子索引，缺失为 -1）、``chain``、``resname``。
    """
    res_index: list[int] = []
    n_idx: list[int] = []
    ca_idx: list[int] = []
    c_idx: list[int] = []
    o_idx: list[int] = []
    chain: list[str] = []
    resname: list[str] = []

    for res in topology.residues:
        names = {}
        for ai in res.atoms:
            name = topology.atoms[ai].name.strip().upper()
            if name in ("N", "CA", "C", "O") and name not in names:
                names[name] = ai
        if "CA" not in names or "C" not in names:
            continue
        res_index.append(res.index)
        n_idx.append(names.get("N", -1))
        ca_idx.append(names["CA"])
        c_idx.append(names["C"])
        o_idx.append(names.get("O", -1))
        chain.append(res.chain)
        resname.append(res.resname)

    return {
        "residue_index": np.asarray(res_index, dtype=int),
        "n": np.asarray(n_idx, dtype=int),
        "ca": np.asarray(ca_idx, dtype=int),
        "c": np.asarray(c_idx, dtype=int),
        "o": np.asarray(o_idx, dtype=int),
        "chain": chain,
        "resname": resname,
    }


def atom_elements(topology: SimpleTopology) -> list[str]:
    """返回每个原子的元素符号（大写）。"""
    return [str(a.element or a.name).strip().upper() for a in topology.atoms]


def atom_chain_ids(topology: SimpleTopology) -> list[str]:
    """返回每个原子的链 ID。"""
    return [a.chain for a in topology.atoms]


def residue_map_by_atom(topology: SimpleTopology) -> dict[int, tuple[str, int]]:
    """构造 原子索引 -> (链, 残基号) 的映射。"""
    mapping: dict[int, tuple[str, int]] = {}
    for res in topology.residues:
        for ai in res.atoms:
            mapping[int(ai)] = (res.chain, int(res.resid))
    return mapping


def _main(argv: list[str] | None = None) -> int:
    """命令行自检入口：打印轨迹概览。"""
    import argparse

    sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "packages"))
    from biopipe_core import get_logger

    logger = get_logger("biopipe.molecular-dynamics.trajectory_io")
    parser = argparse.ArgumentParser(
        prog="trajectory_io", description="读取轨迹文件并打印概览（DCD / PDB / XYZ）"
    )
    parser.add_argument("--input", required=True, help="轨迹文件路径")
    parser.add_argument("--start", type=int, default=0, help="起始帧")
    parser.add_argument("--stop", type=int, default=None, help="结束帧（不含）")
    parser.add_argument("--stride", type=int, default=1, help="帧步长")
    args = parser.parse_args(argv)

    traj = load_trajectory(args.input, start=args.start, stop=args.stop, stride=args.stride)
    logger.info(
        "轨迹 %s: %d 帧, %d 原子, %d 残基, %d 链, 读取器=%s",
        traj.source,
        traj.n_frames,
        traj.n_atoms,
        traj.topology.n_residues,
        traj.topology.n_chains,
        traj.reader,
    )
    logger.info("坐标范围 x=[%.2f, %.2f] Å", traj.coords[..., 0].min(), traj.coords[..., 0].max())
    return 0


if __name__ == "__main__":
    sys.exit(_main())
