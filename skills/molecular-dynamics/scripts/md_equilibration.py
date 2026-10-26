"""分子动力学的分阶段平衡协议（NVT / NPT）、位置限制与链分离监测。

为什么需要这个模块
------------------
复合物（蛋白-蛋白、蛋白-配体、蛋白-DNA）的常规模拟失败模式往往不是
"算不动"，而是**平衡阶段就散架**：链间界面在 300 K 下被热涨落拉开，
或者体系还没从晶体堆积的非物理构象里松弛过来就进入生产期。

标准的应对协议是分阶段的：

1. **能量最小化** —— 消除初始结构中的严重原子重叠；
2. **NVT 平衡** —— 固定盒体积，逐步升温到目标温度，让原子找到局部平衡；
3. **NPT 平衡** —— 开启压力耦合，让密度松弛到真实值；
4. **生产期** —— 采样真实动力学。

其中 NVT 与 NPT 阶段常配**位置限制**（position restraints），把骨架
固定住，只让侧链与溶剂松弛；生产期再撤除限制。

本模块提供这些能力的**纯逻辑部分**：

- :func:`select_restraint_atoms` —— 按选择模式挑出要施加限制的原子
- :func:`build_restraint_force` —— 构造 OpenMM ``CustomExternalForce``
- :func:`chain_separation_report` —— 计算链间质心距离
- :class:`ChainSeparationMonitor` —— 模拟过程中的周期性解离监测
- :func:`stabilize_after_nan` —— NaN 检测与安全回滚

这些函数刻意设计成**不 import OpenMM 就能单测**：只有真正需要
OpenMM 对象时才在函数内部 ``import openmm``。

与老仓库 md_modules/runner.py 的差异（修掉的真实缺陷）
---------------------------------------------------------
1. **NaN 回滚传 None** —— 老代码 ``simulation.context.setPositions(
   modeller.positions if 'modeller' in dir() else None)`` 中 ``modeller``
   从未定义，``dir()`` 返回局部名列表，判断恒为 False，于是向 OpenMM
   传 ``None`` 并抛异常。现在改为显式接收回滚坐标。
2. **Context 固化后加力无效** —— 老代码往 ``simulation.system.addForce()``
   加位置限制，但 Context 已按旧 System 创建，改 System **不会**重建
   Context，力根本没生效。现在要求在创建 Simulation **之前**加力。
3. **Quantity 未转 float** —— 老代码 ``positions[i].x`` 返回带单位的
   ``Quantity``，``np.array`` 得到 object 数组，``.mean(axis=0)`` 直接报错。
   现在统一用 ``value_in_unit(unit.angstroms)``。
4. **链 ID 由 index 猜测** —— 老代码用 ``chr(65 + chain.index)`` 猜链 ID，
   与 PDB 文件里的真实链标识（可能本身就是数字或非字母）不符。
   现在优先用 OpenMM Chain 的 ``id`` 属性。

本模块只依赖 numpy（numpy 缺失时降级为纯 Python 计算）。
"""

from __future__ import annotations

import sys
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

# 允许本模块被单独调用（biopipe_core 是软依赖，取不到时退回标准库 logging）
try:
    sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "packages"))
except (IndexError, OSError):  # pragma: no cover
    pass

try:  # numpy 是可选依赖：缺失时降级到纯 Python
    import numpy as np
except ImportError:  # pragma: no cover - 仅在无 numpy 环境触发
    np = None  # type: ignore[assignment]


__all__ = [
    "RESTRAINT_SELECTIONS",
    "ChainSeparationMonitor",
    "build_restraint_force",
    "chain_separation_report",
    "select_restraint_atoms",
    "stabilize_after_nan",
]

# 支持的位置限制选择模式
RESTRAINT_SELECTIONS = ("backbone", "ca", "heavy", "all")

# 骨架原子名（N / CA / C）—— 位置限制的默认目标
_BACKBONE_NAMES = frozenset({"N", "CA", "C"})


def _require_numpy() -> None:
    """numpy 缺失时抛出带指引的异常。"""
    if np is None:  # pragma: no cover
        raise ImportError(
            "该功能需要 numpy：pip install numpy"
        )


def _positions_to_array(positions: Sequence[Any], unit_module: Any = None) -> Any:
    """把 OpenMM 位置序列转成 (N, 3) 的 numpy 数组。

    老实现直接取 ``positions[i].x``，得到的是带单位的 ``Quantity``，
    拼成的数组 dtype 是 object，后续 ``.mean()`` 会抛类型错误。
    这里统一做单位换算；如果传入的已经是裸 float，直接使用。

    参数
    ----
    positions
        OpenMM 的 ``Vec3`` 列表，或等价的 ``(x, y, z)`` 序列。
    unit_module
        ``openmm.unit`` 模块。提供时会把 Quantity 转成埃。
    """
    _require_numpy()
    rows = []
    for pos in positions:
        if unit_module is not None and hasattr(pos, "value_in_unit"):
            vec = pos.value_in_unit(unit_module.angstroms)
            rows.append((float(vec[0]), float(vec[1]), float(vec[2])))
        elif hasattr(pos, "x"):
            # 裸 Quantity 但没给 unit 模块：退回数值部分
            rows.append((float(pos.x), float(pos.y), float(pos.z)))
        else:
            rows.append((float(pos[0]), float(pos[1]), float(pos[2])))
    return np.array(rows, dtype=float)


def select_restraint_atoms(
    topology: Any,
    selection: str = "backbone",
) -> list[int]:
    """挑选需要施加位置限制的原子索引。

    参数
    ----
    topology
        OpenMM ``Topology`` 对象。
    selection
        选择模式：

        - ``backbone`` —— 主链 N / CA / C（默认，最常用）
        - ``ca`` —— 仅 Cα 原子
        - ``heavy`` —— 所有非氢原子
        - ``all`` —— 所有原子

    返回
    ----
    需要限制的原子索引列表（升序）。
    """
    key = (selection or "backbone").strip().lower()
    if key not in RESTRAINT_SELECTIONS:
        raise ValueError(
            f"未知的位置限制选择 '{selection}'，"
            f"可选：{' / '.join(RESTRAINT_SELECTIONS)}"
        )

    indices: list[int] = []
    for atom in topology.atoms():
        if key == "all":
            indices.append(atom.index)
        elif key == "backbone":
            if atom.name in _BACKBONE_NAMES:
                indices.append(atom.index)
        elif key == "ca":
            if atom.name == "CA":
                indices.append(atom.index)
        elif key == "heavy":
            element = getattr(atom, "element", None)
            symbol = getattr(element, "symbol", None) if element is not None else None
            if symbol is None:
                # 拿不到元素信息时，退回按原子名判断（OpenMM 的原子名不含 H）
                if not atom.name.upper().startswith("H"):
                    indices.append(atom.index)
            elif symbol != "H":
                indices.append(atom.index)
    return indices


def build_restraint_force(
    system: Any,
    atom_indices: Sequence[int],
    reference_positions: Any,
    force_constant: float = 100.0,
) -> Any:
    """构造位置限制力并加到 System 上。

    **必须在创建 Simulation / Context 之前调用** —— Context 一旦创建，
    对 ``system`` 的后续修改不会生效。老仓库的实现正是在这一点上出错：
    它在 NVT 阶段往已存在的 ``simulation.system`` 加力，限制从未真正生效。

    势能形式::

        U = 0.5 * k * ((x-x0)^2 + (y-y0)^2 + (z-z0)^2)

    参数
    ----
    system
        OpenMM ``System`` 对象。
    atom_indices
        要限制的原子索引。
    reference_positions
        参考坐标。元素可以是 ``(x, y, z)`` 三元组，也可以是带单位的
        OpenMM ``Vec3``。
    force_constant
        力常数（kJ/mol/Å²）。默认 100.0。

    返回
    ----
    新建的 ``CustomExternalForce``（已加入 system）。
    """
    from openmm import CustomExternalForce  # 延迟 import，便于无 OpenMM 时单测

    force = CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    force.addPerParticleParameter("k")
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")

    for index in atom_indices:
        ref = reference_positions[index]
        if hasattr(ref, "value_in_unit"):
            # 带单位的 Vec3：转成埃
            from openmm import unit

            vec = ref.value_in_unit(unit.angstroms)
            x0, y0, z0 = float(vec[0]), float(vec[1]), float(vec[2])
        else:
            x0, y0, z0 = float(ref[0]), float(ref[1]), float(ref[2])
        force.addParticle(index, [float(force_constant), x0, y0, z0])

    system.addForce(force)
    return force


@dataclass
class ChainSeparationMonitor:
    """模拟过程中的链分离监测器。

    定期检查各链质心间距，一旦超过阈值即判定复合物正在解离。

    老仓库同名实现的问题：

    - 用 ``chr(65 + chain.index)`` 猜链 ID，与 PDB 真实链标识不符；
    - ``positions[i].x`` 未做单位换算，numpy 数组变成 object dtype；
    - 只在超阈值时打印，缺完整的历史记录与摘要。
    """

    topology: Any
    chain_ids: Optional[Sequence[str]] = None
    threshold: float = 30.0
    report_interval: int = 10000
    initial_distances: Optional[dict[str, float]] = None
    distance_history: list[dict[str, float]] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    def __post_init__(self) -> None:
        self._chain_atom_map: dict[str, list[int]] = {}
        self._chain_labels: dict[str, str] = {}
        for atom in self.topology.atoms():
            chain = atom.residue.chain
            label = self._chain_label(chain)
            if self.chain_ids is None or label in self.chain_ids:
                self._chain_atom_map.setdefault(label, []).append(atom.index)

    def _chain_label(self, chain: Any) -> str:
        """给 OpenMM Chain 生成稳定的标签。

        优先用 ``chain.id``（来自 PDB 的真实链标识），退回到 index 推断。
        老实现**总是**用 index 推断，会把 ``chr(65 + index)`` 当成链 ID，
        与文件里的实际标识脱节。
        """
        chain_id = getattr(chain, "id", None)
        if chain_id:
            return str(chain_id)
        index = getattr(chain, "index", 0)
        return chr(65 + index) if index < 26 else f"Chain{index}"

    def set_initial_positions(self, positions: Sequence[Any], unit_module: Any = None) -> dict[str, float]:
        """记录初始链间距离，作为后续比较的基线。"""
        report = chain_separation_report(
            self.topology,
            positions,
            chain_ids=self.chain_ids,
            unit_module=unit_module,
        )
        self.initial_distances = dict(report["chain_distances"])
        return self.initial_distances

    def check(
        self,
        positions: Sequence[Any],
        unit_module: Any = None,
        step: Optional[int] = None,
    ) -> dict[str, Any]:
        """检查一次链间距离。

        返回
        ----
        含 ``chain_distances`` / ``min_distance`` / ``max_distance`` /
        ``separated`` / ``warnings`` 的字典。
        """
        report = chain_separation_report(
            self.topology,
            positions,
            chain_ids=self.chain_ids,
            unit_module=unit_module,
        )
        distances: dict[str, float] = report["chain_distances"]

        entry = dict(distances)
        if step is not None:
            entry["_step"] = float(step)
        self.distance_history.append(entry)

        current_warnings: list[str] = []
        for pair, distance in distances.items():
            if distance > self.threshold:
                message = f"链 {pair} 间距 {distance:.2f} Å 超过阈值 {self.threshold:.1f} Å"
                current_warnings.append(message)
                self.warnings.append(message)
            if self.initial_distances and pair in self.initial_distances:
                delta = distance - self.initial_distances[pair]
                if abs(delta) > self.threshold:
                    message = (
                        f"链 {pair} 间距相对初始值变化 {delta:+.2f} Å"
                        f"（初始 {self.initial_distances[pair]:.2f} Å）"
                    )
                    current_warnings.append(message)
                    self.warnings.append(message)

        report["separated"] = bool(current_warnings)
        report["warnings"] = current_warnings
        return report

    def get_summary(self) -> dict[str, Any]:
        """获取监测摘要。"""
        return {
            "threshold": self.threshold,
            "report_interval": self.report_interval,
            "initial_distances": self.initial_distances or {},
            "distance_history": self.distance_history,
            "warning_count": len(self.warnings),
            "warnings": self.warnings,
        }


def chain_separation_report(
    topology: Any,
    positions: Sequence[Any],
    chain_ids: Optional[Sequence[str]] = None,
    unit_module: Any = None,
) -> dict[str, Any]:
    """计算各链质心之间的间距。

    参数
    ----
    topology
        OpenMM ``Topology``。
    positions
        原子坐标序列。
    chain_ids
        只统计这些链；``None`` 表示全部。
    unit_module
        ``openmm.unit`` 模块，用于把 Quantity 转成埃。

    返回
    ----
    含 ``chain_centroids`` / ``chain_distances`` / ``min_distance`` /
    ``max_distance`` 的字典。
    """
    _require_numpy()

    labels: dict[int, str] = {}
    atom_indices: dict[str, list[int]] = {}
    for atom in topology.atoms():
        chain = atom.residue.chain
        index = getattr(chain, "index", 0)
        chain_id = getattr(chain, "id", None)
        label = str(chain_id) if chain_id else (
            chr(65 + index) if index < 26 else f"Chain{index}"
        )
        labels[index] = label
        if chain_ids is None or label in chain_ids:
            atom_indices.setdefault(label, []).append(atom.index)

    if not atom_indices:
        return {
            "chain_centroids": {},
            "chain_distances": {},
            "min_distance": 0.0,
            "max_distance": 0.0,
        }

    coords = _positions_to_array(positions, unit_module)

    centroids: dict[str, Any] = {}
    for label, indices in atom_indices.items():
        subset = coords[list(indices)]
        centroids[label] = subset.mean(axis=0)

    labels_sorted = sorted(centroids.keys())
    distances: dict[str, float] = {}
    for i in range(len(labels_sorted)):
        for j in range(i + 1, len(labels_sorted)):
            a, b = labels_sorted[i], labels_sorted[j]
            distances[f"{a}-{b}"] = float(
                np.linalg.norm(centroids[a] - centroids[b])
            )

    values = list(distances.values())
    return {
        "chain_centroids": centroids,
        "chain_distances": distances,
        "min_distance": min(values) if values else 0.0,
        "max_distance": max(values) if values else 0.0,
    }


def _to_float_array(data: Any) -> Any:
    """把坐标数据转成 (N, 3) 或 (N*3,) 的 float 数组。"""
    if hasattr(data, "value_in_unit"):
        from openmm import unit

        values = data.value_in_unit(unit.angstroms)
        return np.array(values, dtype=float)
    return np.array(data, dtype=float)


def stabilize_after_nan(
    positions: Any,
    rollback_positions: Optional[Any] = None,
    tolerance: float = 1e6,
) -> dict[str, Any]:
    """检测坐标中的 NaN / Inf，必要时给出回滚建议。

    老仓库的 NVT 收尾逻辑::

        if np.any(np.isnan(pos_array)):
            simulation.context.setPositions(
                modeller.positions if 'modeller' in dir() else None
            )

    其中 ``modeller`` 从未在函数里定义过，``'modeller' in dir()`` 恒为
    False，于是实际执行的是 ``setPositions(None)`` —— OpenMM 会直接抛
    TypeError，**把一个可恢复的数值问题变成一次崩溃**。

    本函数只做**检测与决策**，不碰 OpenMM 对象，方便单测：

    - ``positions`` —— 被检测的**当前**坐标（模拟末帧）；
    - ``rollback_positions`` —— 回滚目标坐标（通常是输入结构）。

    回滚目标必须自身健全 —— 否则回滚到一个同样坏的结构毫无意义，
    此时会标记 ``rollback_usable=False``。

    返回
    ----
    含 ``has_nan`` / ``n_nan`` / ``n_inf`` / ``exploded`` /
    ``should_rollback`` / ``rollback_usable`` / ``reason`` 的字典。
    """
    _require_numpy()

    if rollback_positions is None:
        raise ValueError("rollback_positions 不能为 None")

    current = _to_float_array(positions)
    rollback = _to_float_array(rollback_positions)

    flat = current.reshape(-1)
    nan_mask = np.isnan(flat)
    inf_mask = np.isinf(flat)
    n_nan = int(nan_mask.sum())
    n_inf = int(inf_mask.sum())

    has_nan = n_nan > 0 or n_inf > 0
    magnitudes = np.abs(flat[~nan_mask & ~inf_mask]) if flat.size else np.array([])
    exploded = bool(magnitudes.size and magnitudes.max() > tolerance)

    # 回滚目标自身是否健全
    rollback_flat = rollback.reshape(-1)
    rollback_usable = bool(
        not np.isnan(rollback_flat).any() and not np.isinf(rollback_flat).any()
    )

    if has_nan or exploded:
        reason_parts = []
        if n_nan:
            reason_parts.append(f"{n_nan} 个 NaN")
        if n_inf:
            reason_parts.append(f"{n_inf} 个 Inf")
        if exploded:
            reason_parts.append(
                f"坐标量级超过 {tolerance:g} Å（疑似数值发散）"
            )
        if not rollback_usable:
            reason_parts.append("回滚目标本身不健全，无法安全回滚")
        reason = "、".join(reason_parts)
    else:
        reason = ""

    return {
        "has_nan": has_nan,
        "n_nan": n_nan,
        "n_inf": n_inf,
        "exploded": exploded,
        "should_rollback": has_nan or exploded,
        "rollback_usable": rollback_usable,
        "reason": reason,
    }


# ---------------------------------------------------------------- 自检 CLI


def _get_logger() -> Any:
    """获取 logger。

    优先用 ``biopipe_core``（仓库统一日志格式）；拿不到时退回标准库
    ``logging``，保证本模块在裸环境下也能用。
    """
    try:
        from biopipe_core import get_logger

        return get_logger("biopipe.molecular-dynamics")
    except Exception:  # noqa: BLE001
        import logging

        return logging.getLogger("biopipe.molecular-dynamics")


def _main(argv: Optional[list[str]] = None) -> int:
    """命令行入口：对一个 PDB 报告链间距离与位置限制建议。

    用途：

    - 在跑长模拟**之前**先看一眼复合物的链间距离，
      判断起始姿态是否合理、阈值该设多少；
    - 检查一个结构里有几链、每链多少原子。

    用法::

        python md_equilibration.py --input complex.pdb [--threshold 30]
        python md_equilibration.py --input complex.pdb --json result.json
    """
    import argparse
    import json

    parser = argparse.ArgumentParser(
        prog="md_equilibration",
        description="检查结构的链间距离与位置限制建议（需要 OpenMM 读 PDB）",
    )
    parser.add_argument("--input", required=True, metavar="PATH", help="输入 PDB 文件")
    parser.add_argument("--threshold", type=float, default=30.0,
                        help="链间距告警阈值（Å，默认 30）")
    parser.add_argument("--json", metavar="PATH",
                        help="把结果写成 JSON 文件（不写则用日志输出）")
    args = parser.parse_args(argv)

    logger = _get_logger()

    structure = Path(args.input)
    if not structure.exists():
        logger.error("输入文件不存在: %s", structure)
        return 1

    try:
        from openmm import app, unit
    except Exception as exc:  # noqa: BLE001
        logger.error(
            "读取 PDB 需要 OpenMM：%s\n"
            "  conda install -c conda-forge openmm\n"
            "  或 pip install openmm",
            exc,
        )
        return 1

    pdb = app.PDBFile(str(structure))
    topology = pdb.topology
    positions = pdb.positions

    report = chain_separation_report(topology, positions, unit_module=unit)
    monitor = ChainSeparationMonitor(topology, threshold=args.threshold)
    monitor.set_initial_positions(positions, unit)

    restraint_counts = {
        mode: len(select_restraint_atoms(topology, mode))
        for mode in RESTRAINT_SELECTIONS
    }

    separated = bool(
        report["chain_distances"] and report["max_distance"] > args.threshold
    )

    result = {
        "file": str(structure),
        "n_atoms": topology.getNumAtoms(),
        "n_residues": topology.getNumResidues(),
        "n_chains": len(report["chain_centroids"]),
        "chain_distances": report["chain_distances"],
        "min_distance": report["min_distance"],
        "max_distance": report["max_distance"],
        "threshold": args.threshold,
        "separated": separated,
        "restraint_atom_counts": restraint_counts,
    }

    if args.json:
        out_path = Path(args.json)
        out_path.write_text(
            json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8"
        )
        logger.info("已写出 %s", out_path)
        return 0

    logger.info("结构: %s", structure)
    logger.info(
        "  原子 %d  残基 %d  链 %d",
        result["n_atoms"], result["n_residues"], result["n_chains"],
    )
    if report["chain_distances"]:
        logger.info("  链间质心距离（阈值 %.4g Å）:", args.threshold)
        for pair, distance in sorted(report["chain_distances"].items()):
            flag = "  <== 超过阈值" if distance > args.threshold else ""
            logger.info("    %s: %.2f Å%s", pair, distance, flag)
    else:
        logger.info("  单链体系，无需链分离监测")
    logger.info("  位置限制可选原子数:")
    for mode, count in restraint_counts.items():
        logger.info("    %-9s %d", mode, count)
    return 0


if __name__ == "__main__":
    sys.exit(_main())
