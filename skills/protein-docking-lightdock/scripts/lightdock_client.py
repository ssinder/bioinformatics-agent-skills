#!/usr/bin/env python3
"""LightDock 蛋白质对接的本地执行客户端。

对应老仓库 ``dock_lightdock/scripts/`` 下的 5 个脚本。

**LightDock 与 HDOCK 是完全不同的技术栈**：HDOCK 是网络服务，提交后
由服务端算；LightDock 是**本地程序**（``pip install lightdock``），
对接计算全部发生在你自己的机器上，CPU 密集，可能跑几小时。所以这里的
设计重点不是"提交后取结果"，而是"如何安全地组织一次本地批量计算"。

从老仓库继承并修掉的问题
------------------------
1. **硬编码个人路径**。老脚本把作者的 macOS 桌面路径与 conda 环境里的
   LightDock 可执行文件目录写死在全部 5 个脚本里，换台机器立刻失效。
   本模块改为在 PATH 与 ``<sys.prefix>/bin`` 里自动探测可执行文件，
   也允许显式指定路径。
2. **三份几乎逐字重复的代码**。``dock_completess.py``、``dock_simple.py``、
   ``dock_test_simple.py`` 各自复制了一份 ``run_docking`` /
   ``extract_best_score`` / ``generate_scores_csv``，改一处要改三处。
   本模块只保留一份实现。
3. **把 LigID 当成 glowworm 索引**。见下方"GSO 文件格式"。
4. **字典序排序取"第一个"构象**。见下方"构象文件命名"。
5. **坐标系错配**。见下方"坐标系陷阱"。

GSO 文件格式（最容易读错的地方）
--------------------------------
LightDock 每个 swarm 目录下有 ``gso_<step>.out``，**每行一个 glowworm**。
真实格式（官方 0.9.4 ``lightdock/util/analysis.py`` 里的
``read_lightdock_output`` 为准）::

    (tx, ty, tz, qw, qx, qy, qz[, anm...])  RecID LigID Luciferin Neighbors VisionRange Scoring
    (31.4171143, 1.8570079, -6.3956223, -0.1058407, -0.4849369, 0.5997430, -0.6276482) 0 0 11.25395618 0 4.200 7.52800101

字段依次是：坐标括号、四元数、**受体结构索引 RecID**、**配体结构索引
LigID**、**萤光素值 luciferin**、**邻居数**、**视野半径**、**打���分**。

**关键**：glowworm 的编号（``id_glowworm``）是**第几个坐标行**，**不是**文件
里的任何一列，也**不是**物理行号（注释行不占号）。官方
``read_lightdock_output`` 里写得很清楚::

    raw_lines = [line for line in fin if line[0] != "#"]
    for id_line, line in enumerate(raw_lines):
        ...
        id_glowworm=id_line,
        receptor_id=int(rest[0]),
        ligand_id=int(rest[1]),

``lgd_generate_conformations.py`` 的 ``parse_output_file`` 用的是同一口径
（只在 ``line[0] == "("`` 时递增计数器），所以 ``lightdock_<i>.pdb`` 里的
``i`` 就是 glowworm 编号。两边必须一致，否则按编号挑构象会挑错。

老脚本的正则 ``r'\\([^\\)]+\\)\\s+\\d+\\s+(\\d+)\\s+...'`` 把捕获组 1
命名为 ``glowworm_idx`` 并当成 glowworm 编号，但它捕获的是 **LigID**——
单体对接时 LigID 恒为 ``0``。于是"最佳构象"永远是 glowworm #0，
与真实得分最高的 glowworm 毫无关系。这个 bug 产出的结果看起来完全合理，
极难察觉。

构象文件命名
------------
``lgd_generate_conformations.py`` 输出 ``lightdock_<i>.pdb``，**i 从 0 起**
且**没有补零**。所以按字典序排会得到::

    lightdock_0.pdb, lightdock_1.pdb, lightdock_10.pdb, lightdock_11.pdb, lightdock_2.pdb, ...

``sorted(glob(...))[0]`` 恰好能拿到 ``lightdock_0.pdb``（看起来对），
但任何"取第 N 个"的下标操作全错。老脚本用 ``sorted(...)[0]`` 取"最佳",
在只生成 1 个构象时碰巧可用，一旦生成多个就取错。

坐标系陷阱
----------
``lightdock3_setup.py`` 会把受体与配体**平移到坐标原点**，并另存为
``lightdock_<原文件名>.pdb``。GSO 文件里记录的平移量正是相对这个原点的。

所以生成构象时**必须用 setup 产出的 ``lightdock_receptor.pdb`` /
``lightdock_ligand.pdb``**，而不是原始输入的 ``receptor.pdb`` /
``ligand.pdb``。老脚本传的是原始文件，复合物的两个部分因此对不上——
生成的 PDB 里配体位置是错的。这种错误不会报错，只会让结果"看起来像那么
回事"，是最危险的一类 bug。

官方命令清单（0.9.4，来自 ``bin/`` 目录）
-----------------------------------------
``lightdock3_setup.py``（准备）、``lightdock3.py``（运行）、
``lgd_generate_conformations.py``（生成构象）、``lgd_rank.py``（排名）、
``lgd_top.py``（提取 top 模型）等。这些是通过 setup.py 的 ``scripts=``
声明安装的**裸文件**，不是 console_scripts 入口点，因此不会出现在
``pip show`` 的 entry point 列表里，但确实会落在 ``<env>/bin/`` 下。
"""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

#: LightDock 0.9.4 的准备脚本名（在 <env>/bin/ 下）。
SETUP_SCRIPT = "lightdock3_setup.py"

#: LightDock 0.9.4 的运行脚本名（在 <env>/bin/ 下）。
RUN_SCRIPT = "lightdock3.py"

#: 生成构象的脚本名。
CONFORMATIONS_SCRIPT = "lgd_generate_conformations.py"

#: 排名脚本名。
RANK_SCRIPT = "lgd_rank.py"

#: 提取 top 模型的脚本名。
TOP_SCRIPT = "lgd_top.py"

#: setup 阶段把结构移到原点后另存的文件前缀（``lightdock_%s``）。
LIGHTDOCK_PREFIX = "lightdock_"

#: 单个 GSO 行的解析正则。
#:
#: 坐标括号里可能有 7 个数（刚性）或更多（ANM 开启），所以用
#: ``[^)]*`` 整体匹配括号内容，再按空白切分后面的定长字段。
#:
#: 注意捕获组的含义：
#: - group(1)  坐标部分（丢弃）
#: - group(2)  RecID  受体结构索引
#: - group(3)  LigID  配体结构索引
#: - group(4)  luciferin  萤光素值
#: - group(5)  neighbors  邻居数
#: - group(6)  vision_range  视野半径
#: - group(7)  scoring  打分（DFIRE 等）
#:
#: **glowworm 编号不在这里**——它是行号。
GSO_LINE_RE = re.compile(
    r"^\s*\((?P<coord>[^)]*)\)\s+"
    r"(?P<rec_id>\d+)\s+"
    r"(?P<lig_id>\d+)\s+"
    r"(?P<luciferin>[-+0-9.eE]+)\s+"
    r"(?P<neighbors>\d+)\s+"
    r"(?P<vision_range>[-+0-9.eE]+)\s+"
    r"(?P<scoring>[-+0-9.eE]+)\s*$"
)

#: 兼容老脚本用的宽松正则（只有坐标 + 若干数值），用于行格式有出入时兜底。
#: 注意它**不**区分 RecID / LigID，因此只能用来取打分，不能用来定位 glowworm。
GSO_FALLBACK_RE = re.compile(r"\([^\)]+\)\s+\d+\s+\d+\s+[-+0-9.eE]+\s+\d+\s+[-+0-9.eE]+\s+(?P<scoring>[-+0-9.eE]+)")


class LightDockError(RuntimeError):
    """LightDock 执行或解析失败。"""


@dataclass(frozen=True)
class GsoEntry:
    """GSO 输出里的一行，即一个 glowworm 的最终状态。

    Attributes:
        glowworm: glowworm 编号，即**第几个坐标行**（注释行不占号，从 0 开始）。
        rec_id: 受体结构索引。
        lig_id: 配体结构索引。
        luciferin: 萤光素值，GSO 算法里与打分成正比。
        neighbors: 该 glowworm 视野内的邻居数量。
        vision_range: 当前视野半径。
        scoring: 打分函数原始值（默认 DFIRE）。**越负越好**。
    """

    glowworm: int
    rec_id: int
    lig_id: int
    luciferin: float
    neighbors: int
    vision_range: float
    scoring: float

    def to_row(self) -> dict:
        """转成扁平字典，便于写 CSV。"""
        return {
            "glowworm": self.glowworm,
            "rec_id": self.rec_id,
            "lig_id": self.lig_id,
            "luciferin": self.luciferin,
            "neighbors": self.neighbors,
            "vision_range": self.vision_range,
            "scoring": self.scoring,
        }


@dataclass(frozen=True)
class SwarmResult:
    """一个 swarm 的最终结果。

    Attributes:
        swarm: swarm 编号（对应目录名 ``swarm_<n>``）。
        gso_file: 用于取分的 GSO 文件路径。
        best: 得分最低（最优）的 glowworm 条目。
        entries: 该 swarm 全部 glowworm 条目。
    """

    swarm: int
    gso_file: Path
    best: Optional[GsoEntry]
    entries: list[GsoEntry] = field(default_factory=list)

    @property
    def best_score(self) -> Optional[float]:
        """最优得分，无有效条目时返回 ``None``。"""
        return self.best.scoring if self.best else None


@dataclass(frozen=True)
class DockingJob:
    """一次对接作业的完整记录。

    Attributes:
        job_id: 作业标识。
        receptor: 受体 PDB 路径。
        ligand: 配体 PDB 路径。
        output_dir: 该作业的输出目录。
        status: ``success`` / ``failed`` / ``skipped``。
        best_score: 全局最优打分，越负越好。
        best_swarm: 最优打分所在的 swarm 编号。
        best_glowworm: 最优打分所在的 glowworm 编号。
        message: 失败或跳过的原因。
    """

    job_id: str
    receptor: Path
    ligand: Path
    output_dir: Path
    status: str
    best_score: Optional[float] = None
    best_swarm: Optional[int] = None
    best_glowworm: Optional[int] = None
    message: str = ""

    def to_row(self) -> dict:
        """转成扁平字典，便于写 CSV。"""
        return {
            "job_id": self.job_id,
            "receptor": self.receptor.name,
            "ligand": self.ligand.name,
            "status": self.status,
            "best_score": (
                f"{self.best_score:.4f}" if self.best_score is not None else ""
            ),
            "best_swarm": "" if self.best_swarm is None else self.best_swarm,
            "best_glowworm": (
                "" if self.best_glowworm is None else self.best_glowworm
            ),
            "output_dir": str(self.output_dir),
            "message": self.message,
        }


def parse_gso_line(line: str, glowworm: int) -> Optional[GsoEntry]:
    """解析一行 GSO 输出。

    Args:
        line: GSO 文件中的一行。
        glowworm: 该行的 glowworm 编号，即**它是第几个坐标行**
            （注释行不占号）。

    Returns:
        解析成功返回 :class:`GsoEntry`；行是注释或格式无法识别时返回 ``None``。
    """
    stripped = line.strip()
    if not stripped or stripped.startswith("#"):
        return None

    match = GSO_LINE_RE.match(stripped)
    if match:
        try:
            return GsoEntry(
                glowworm=glowworm,
                rec_id=int(match.group("rec_id")),
                lig_id=int(match.group("lig_id")),
                luciferin=float(match.group("luciferin")),
                neighbors=int(match.group("neighbors")),
                vision_range=float(match.group("vision_range")),
                scoring=float(match.group("scoring")),
            )
        except ValueError:
            return None

    # 兜底：字段布局与预期不符时，至少把打分取出来。
    fallback = GSO_FALLBACK_RE.search(stripped)
    if fallback:
        try:
            return GsoEntry(
                glowworm=glowworm,
                rec_id=0,
                lig_id=0,
                luciferin=0.0,
                neighbors=0,
                vision_range=0.0,
                scoring=float(fallback.group("scoring")),
            )
        except ValueError:
            return None
    return None


def parse_gso_file(path: Path) -> list[GsoEntry]:
    """解析整个 GSO 文件。

    **glowworm 编号是"第几个坐标行"，不是文件里的物理行号。**

    官方 ``read_lightdock_output`` 先剔除注释行再编号::

        raw_lines = [line for line in fin if line[0] != "#"]
        for id_line, line in enumerate(raw_lines):
            ... id_glowworm=id_line ...

    ``lgd_generate_conformations.py`` 的 ``parse_output_file`` 也是同样的
    口径——只在 ``line[0] == "("`` 时递增计数器。两者必须一致，否则
    "最优 glowworm 的编号"与"``lightdock_<i>.pdb`` 里的 i"会差一位，
    于是按编号去挑构象文件时挑错。

    Args:
        path: ``gso_<step>.out`` 文件路径。

    Returns:
        所有成功解析的 glowworm 条目，按编号顺序。

    Raises:
        LightDockError: 文件不存在或无法读取。
    """
    path = Path(path)
    if not path.is_file():
        raise LightDockError(f"GSO 文件不存在: {path}")

    entries: list[GsoEntry] = []
    try:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            glowworm = 0
            for line in handle:
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    # 注释行不占编号——与官方口径一致
                    continue
                entry = parse_gso_line(line, glowworm)
                if entry is not None:
                    entries.append(entry)
                glowworm += 1
    except OSError as exc:
        raise LightDockError(f"读取 GSO 文件失败: {path} —— {exc}") from exc
    return entries


def find_latest_gso(swarm_dir: Path) -> Optional[Path]:
    """找出 swarm 目录里步数最大的 GSO 文件。

    LightDock 每 10 步存一次结果，文件名形如 ``gso_10.out``、``gso_100.out``。
    取**步数最大**的那个才是收敛结果。

    老脚本用 ``sorted(glob(...))[-1]`` 取最后一个——字典序下
    ``gso_100.out`` 排在 ``gso_20.out`` 前面，于是 100 步的结果会被
    20 步的结果覆盖掉。这是一个真实的静默 bug。

    Args:
        swarm_dir: ``swarm_<n>`` 目录。

    Returns:
        步数最大的 GSO 文件；没有则返回 ``None``。
    """
    swarm_dir = Path(swarm_dir)
    if not swarm_dir.is_dir():
        return None

    best: Optional[tuple[int, Path]] = None
    for candidate in swarm_dir.glob("gso_*.out"):
        match = re.search(r"gso_(\d+)\.out$", candidate.name)
        if not match:
            continue
        step = int(match.group(1))
        if best is None or step > best[0]:
            best = (step, candidate)
    return best[1] if best else None


def collect_swarm_results(output_dir: Path) -> list[SwarmResult]:
    """扫描输出目录下的所有 swarm，提取每个 swarm 的最优 glowworm。

    Args:
        output_dir: 对接作业的输出目录（内含 ``swarm_*`` 子目录）。

    Returns:
        每个 swarm 一个 :class:`SwarmResult`，按 swarm 编号排序。
    """
    output_dir = Path(output_dir)
    results: list[SwarmResult] = []

    for swarm_dir in sorted(output_dir.glob("swarm_*")):
        match = re.search(r"swarm_(\d+)$", swarm_dir.name)
        swarm_id = int(match.group(1)) if match else len(results)
        gso_file = find_latest_gso(swarm_dir)
        if gso_file is None:
            continue
        entries = parse_gso_file(gso_file)
        best = min(entries, key=lambda e: e.scoring) if entries else None
        results.append(
            SwarmResult(swarm=swarm_id, gso_file=gso_file, best=best, entries=entries)
        )
    return results


def find_global_best(output_dir: Path) -> Optional[tuple[SwarmResult, GsoEntry]]:
    """找出整个作业目录里得分最优的 glowworm。

    Args:
        output_dir: 对接作业的输出目录。

    Returns:
        ``(swarm 结果, 最优条目)``；没有任何有效条目时返回 ``None``。
    """
    best_pair: Optional[tuple[SwarmResult, GsoEntry]] = None
    for swarm in collect_swarm_results(output_dir):
        if swarm.best is None:
            continue
        if best_pair is None or swarm.best.scoring < best_pair[1].scoring:
            best_pair = (swarm, swarm.best)
    return best_pair


def locate_executable(name: str, explicit: Optional[str] = None) -> Optional[Path]:
    """定位 LightDock 的可执行脚本。

    探测顺序：

    1. 显式传入的路径；
    2. 当前 Python 环境的 ``<prefix>/bin/``（conda / venv 装 lightdock 时
       脚本落在这里，但**不在 PATH 里**——这是老脚本踩坑的原因）；
    3. 系统 ``PATH``。

    Args:
        name: 脚本名，如 ``lightdock3_setup.py``。
        explicit: 用户显式指定的路径。

    Returns:
        找到则返回路径；找不到返回 ``None``。
    """
    if explicit:
        candidate = Path(explicit).expanduser()
        if candidate.is_file():
            return candidate
        # 允许只给目录
        if candidate.is_dir():
            for sub in (candidate / name, candidate / "bin" / name):
                if sub.is_file():
                    return sub
        return None

    # conda / venv 环境：脚本在 <prefix>/bin 下
    prefix = Path(sys.prefix)
    for sub in (prefix / "bin" / name, prefix / "Scripts" / name):
        if sub.is_file():
            return sub

    found = shutil.which(name)
    return Path(found) if found else None


def resolve_executables(
    setup: Optional[str] = None,
    run: Optional[str] = None,
    conformations: Optional[str] = None,
    rank: Optional[str] = None,
    top: Optional[str] = None,
) -> dict[str, Path]:
    """一次性定位全部需要的可执行脚本。

    Args:
        setup: ``lightdock3_setup.py`` 路径。
        run: ``lightdock3.py`` 路径。
        conformations: ``lgd_generate_conformations.py`` 路径。
        rank: ``lgd_rank.py`` 路径。
        top: ``lgd_top.py`` 路径。

    Returns:
        脚本名到路径的映射。

    Raises:
        LightDockError: 缺少 ``setup`` 或 ``run`` 这两个必需脚本。
    """
    mapping = {
        SETUP_SCRIPT: locate_executable(SETUP_SCRIPT, setup),
        RUN_SCRIPT: locate_executable(RUN_SCRIPT, run),
        CONFORMATIONS_SCRIPT: locate_executable(CONFORMATIONS_SCRIPT, conformations),
        RANK_SCRIPT: locate_executable(RANK_SCRIPT, rank),
        TOP_SCRIPT: locate_executable(TOP_SCRIPT, top),
    }
    missing = [name for name in (SETUP_SCRIPT, RUN_SCRIPT) if mapping[name] is None]
    if missing:
        raise LightDockError(
            "找不到 LightDock 可执行脚本: "
            + ", ".join(missing)
            + "。请先 `pip install lightdock==0.9.4`（Windows 官方不支持，"
            "请用 Linux/macOS），或用 --setup-exe/--run-exe 显式指定路径。"
        )
    return {name: path for name, path in mapping.items() if path is not None}


def run_command(
    command: Sequence[str],
    cwd: Path,
    timeout: int = 3600,
    env: Optional[dict[str, str]] = None,
) -> subprocess.CompletedProcess:
    """执行一条外部命令并返回结果。

    Args:
        command: 命令参数列表。
        cwd: 工作目录。
        timeout: 超时秒数。
        env: 追加的环境变量。

    Returns:
        ``subprocess.CompletedProcess``。

    Raises:
        LightDockError: 命令超时或无法启动。
    """
    run_env = None
    if env:
        run_env = os.environ.copy()
        run_env.update(env)

    try:
        return subprocess.run(
            list(command),
            cwd=str(cwd),
            capture_output=True,
            text=True,
            timeout=timeout,
            env=run_env,
        )
    except subprocess.TimeoutExpired as exc:
        raise LightDockError(
            f"命令执行超时（{timeout}s）: {' '.join(command)}"
        ) from exc
    except OSError as exc:
        raise LightDockError(
            f"命令无法启动: {' '.join(command)} —— {exc}"
        ) from exc


def setup_docking(
    receptor: Path,
    ligand: Path,
    output_dir: Path,
    setup_exe: Path,
    swarms: int = 50,
    glowworms: int = 50,
    timeout: int = 600,
    no_hydrogens: bool = True,
    no_oxt: bool = True,
    no_water: bool = True,
) -> Path:
    """执行 LightDock 的准备阶段。

    会把受体与配体复制进 ``output_dir`` 并运行 ``lightdock3_setup.py``。
    准备阶段会把两个结构**平移到坐标原点**并另存为
    ``lightdock_receptor.pdb`` / ``lightdock_ligand.pdb``——后续生成构象
    必须用这两个文件。

    Args:
        receptor: 受体 PDB 路径。
        ligand: 配体 PDB 路径。
        output_dir: 输出目录。
        setup_exe: ``lightdock3_setup.py`` 路径。
        swarms: swarm 数量。
        glowworms: 每个 swarm 的 glowworm 数量。
        timeout: 超时秒数。
        no_hydrogens: 移除氢原子。
        no_oxt: 移除 OXT 原子。
        no_water: 移除水分子。

    Returns:
        准备阶段生成的 ``setup.json`` 路径。

    Raises:
        LightDockError: 准备阶段失败。
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    receptor_copy = output_dir / "receptor.pdb"
    ligand_copy = output_dir / "ligand.pdb"
    shutil.copy2(receptor, receptor_copy)
    shutil.copy2(ligand, ligand_copy)

    command = [
        str(setup_exe),
        str(receptor_copy),
        str(ligand_copy),
        "-s",
        str(swarms),
        "-g",
        str(glowworms),
    ]
    if no_hydrogens:
        command.append("--noh")
    if no_oxt:
        command.append("--noxt")
    if no_water:
        command.append("--now")

    result = run_command(command, cwd=output_dir, timeout=timeout)
    if result.returncode != 0:
        detail = (result.stderr or result.stdout or "").strip()[-2000:]
        raise LightDockError(f"LightDock 准备阶段失败（退出码 {result.returncode}）:\n{detail}")

    setup_json = output_dir / "setup.json"
    if not setup_json.is_file():
        raise LightDockError(
            f"准备阶段没有产出 setup.json，LightDock 可能装坏了: {output_dir}"
        )
    return setup_json


def run_docking(
    setup_json: Path,
    output_dir: Path,
    run_exe: Path,
    steps: int = 50,
    cores: int = 8,
    timeout: int = 86400,
) -> None:
    """执行 LightDock 的模拟阶段。

    Args:
        setup_json: 准备阶段产出的 ``setup.json``。
        output_dir: 输出目录（也是工作目录）。
        run_exe: ``lightdock3.py`` 路径。
        steps: 模拟步数。
        cores: 并行 CPU 核数。
        timeout: 超时秒数，默认一天。

    Raises:
        LightDockError: 模拟阶段失败。
    """
    result = run_command(
        [str(run_exe), "setup.json", str(steps), "-c", str(cores)],
        cwd=Path(output_dir),
        timeout=timeout,
    )
    if result.returncode != 0:
        detail = (result.stderr or result.stdout or "").strip()[-2000:]
        raise LightDockError(
            f"LightDock 模拟阶段失败（退出码 {result.returncode}）:\n{detail}"
        )


def moved_structure_path(original: Path) -> Path:
    """返回 setup 阶段平移后另存的同名结构路径。

    ``lightdock3_setup.py`` 用 ``lightdock_%s`` 前缀保存平移后的结构。
    """
    original = Path(original)
    return original.parent / f"{LIGHTDOCK_PREFIX}{original.name}"


def generate_conformations(
    output_dir: Path,
    gso_file: Path,
    num_glowworms: int = 1,
    conformations_exe: Optional[Path] = None,
    timeout: int = 600,
) -> list[Path]:
    """从 GSO 结果生成复合物 PDB。

    **必须传 setup 阶段产出的平移后结构**（``lightdock_receptor.pdb`` /
    ``lightdock_ligand.pdb``），否则配体与受体会错位。

    Args:
        output_dir: 对接作业的输出目录。
        gso_file: 用于生成构象的 GSO 文件。
        num_glowworms: 生成多少个构象。
        conformations_exe: ``lgd_generate_conformations.py`` 路径。
        timeout: 超时秒数。

    Returns:
        生成的 PDB 文件列表（按 glowworm 编号数值排序）。

    Raises:
        LightDockError: 找不到可执行脚本或生成失败。
    """
    output_dir = Path(output_dir)
    gso_file = Path(gso_file)

    if conformations_exe is None:
        conformations_exe = locate_executable(CONFORMATIONS_SCRIPT)
    if conformations_exe is None:
        raise LightDockError(
            f"找不到 {CONFORMATIONS_SCRIPT}，请用 --conformations-exe 指定路径"
        )

    receptor = moved_structure_path(output_dir / "receptor.pdb")
    ligand = moved_structure_path(output_dir / "ligand.pdb")
    for structure in (receptor, ligand):
        if not structure.is_file():
            raise LightDockError(
                f"缺少平移后的结构文件 {structure.name}——"
                "请确认准备阶段已成功完成（该文件由 lightdock3_setup.py 生成）"
            )

    # 官方用法是在 swarm 目录下运行，参考相对路径 ../lightdock_receptor.pdb
    swarm_dir = gso_file.parent
    command = [
        str(conformations_exe),
        f"../{receptor.name}",
        f"../{ligand.name}",
        gso_file.name,
        str(num_glowworms),
    ]
    result = run_command(command, cwd=swarm_dir, timeout=timeout)
    if result.returncode != 0:
        detail = (result.stderr or result.stdout or "").strip()[-2000:]
        raise LightDockError(
            f"生成构象失败（退出码 {result.returncode}）:\n{detail}"
        )

    return sort_conformation_files(swarm_dir)


def sort_conformation_files(directory: Path) -> list[Path]:
    """按 glowworm 编号**数值**排序构象文件。

    文件名形如 ``lightdock_0.pdb``、``lightdock_10.pdb``，没有补零，
    字典序会把 10 排到 2 前面。

    Args:
        directory: 构象文件所在目录。

    Returns:
        排序后的文件列表。
    """
    directory = Path(directory)

    def key(path: Path) -> tuple[int, str]:
        match = re.search(r"lightdock_(\d+)\.pdb$", path.name)
        if match:
            return (int(match.group(1)), path.name)
        return (10**9, path.name)

    return sorted(directory.glob("lightdock_*.pdb"), key=key)


def read_setup_config(setup_json: Path) -> dict:
    """读取 setup.json 的配置。

    Args:
        setup_json: 配置文件路径。

    Returns:
        解析后的字典；文件不存在或损坏时返回空字典。
    """
    path = Path(setup_json)
    if not path.is_file():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return {}
