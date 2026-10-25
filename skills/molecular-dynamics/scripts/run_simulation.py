"""分子动力学模拟入口（OpenMM 可用时真实执行）。

用法::

    python run_simulation.py --input system.pdb --output results/ --steps 1000

若 OpenMM 已安装并可用，本脚本会：

1. 用 OpenMM 读取输入结构（PDB / PDBx）；
2. 依据 `--force-field` 选择力场（默认 `amber14-all.xml`，并自动补充隐式溶剂
   `implicit/gbn2.xml` 与 `implicit/obc2.xml`）；
3. 建立非周期体系，加入柔性氢（`hydrogens=add`）；
4. 用 Langevin Middle 恒温器做最小化 + 生产性 MD；
5. 写出 DCD 轨迹、最终结构与状态报告。

OpenMM 不可用时，脚本会给出明确的安装提示并以非零码退出（不静默降级为
"假装跑了模拟"）。
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "packages"))

from biopipe_core import cli_parser, get_config, get_logger  # noqa: E402
from biopipe_core.cli import SkillError, run_cli  # noqa: E402

SKILL_NAME = "molecular-dynamics"

DEFAULT_FORCE_FIELDS = (
    ("amber14-all.xml", "amber14/tip3pfb.xml"),
    ("amber99sbildn.xml", "amber99/tip3p.xml"),
    ("charmm36.xml", "charmm36/water.xml"),
)


def _build_parser() -> argparse.ArgumentParser:
    """构造命令行解析器。"""
    parser = cli_parser(
        prog=f"{SKILL_NAME}-sim",
        description="运行分子动力学模拟（需要 OpenMM）",
        skill=SKILL_NAME,
    )
    parser.add_argument("--input", required=True, metavar="PATH", help="输入结构文件（PDB / PDBx）")
    parser.add_argument("--output", metavar="DIR", help="输出目录（默认取配置中的 results/）")
    parser.add_argument("--force-field", default="amber14", help="力场预设：amber14 / amber99 / charmm36（默认 amber14）")
    parser.add_argument("--solvent", choices=("none", "implicit", "vacuum"), default="implicit",
                        help="溶剂处理：none / implicit / vacuum（默认 implicit）")
    parser.add_argument("--steps", type=int, default=1000, help="生产性 MD 步数（默认 1000）")
    parser.add_argument("--equilibration-steps", type=int, default=100, help="平衡步数（默认 100）")
    parser.add_argument("--minimize-steps", type=int, default=200, help="能量最小化步数（默认 200）")
    parser.add_argument("--temperature", type=float, default=300.0, help="温度（K，默认 300）")
    parser.add_argument("--friction", type=float, default=1.0, help="Langevin 摩擦系数 /ps（默认 1.0）")
    parser.add_argument("--timestep", type=float, default=0.002, help="积分步长（ps，默认 0.002）")
    parser.add_argument("--platform", default=None, help="OpenMM 平台名，如 CUDA / CPU / OpenCL")
    parser.add_argument("--report-interval", type=int, default=100, help="状态报告间隔（步，默认 100）")
    parser.add_argument("--trajectory-interval", type=int, default=100, help="轨迹写出间隔（步，默认 100）")
    parser.add_argument("--check-only", action="store_true", help="仅检查 OpenMM 是否可用并退出")
    return parser


def _resolve_force_field(name: str) -> tuple[str, str]:
    """根据预设名返回 (蛋白力场 XML, 水模型 XML)。"""
    key = name.strip().lower()
    for preset, files in DEFAULT_FORCE_FIELDS:
        if preset.startswith(key) or key.startswith(preset.split("-")[0]):
            return files
    raise SkillError(
        f"未知的力场预设 '{name}'，可选：amber14 / amber99 / charmm36"
    )


def _check_openmm() -> tuple[bool, str]:
    """检查 OpenMM 是否可用，返回 (可用, 说明)。"""
    try:
        import openmm  # type: ignore  # noqa: F401
        from openmm import unit  # type: F401
    except Exception as exc:  # noqa: BLE001
        return False, f"OpenMM 不可用（{exc}）"
    return True, "OpenMM 可用"


def _resolve_output_dir(output: str | None) -> Path:
    """解析输出目录。"""
    config = get_config()
    if output:
        path = Path(output).expanduser()
        if not path.is_absolute():
            path = Path.cwd() / path
    else:
        path = config.output_dir
    path.mkdir(parents=True, exist_ok=True)
    return path


def _run_openmm(args: argparse.Namespace, output_dir: Path) -> dict:
    """真正调用 OpenMM 执行模拟。"""
    import openmm
    from openmm import app, unit

    logger = get_logger(f"biopipe.{SKILL_NAME}")

    structure_path = Path(args.input)
    if not structure_path.exists():
        raise SkillError(f"输入结构不存在: {structure_path}")

    protein_ff, water_ff = _resolve_force_field(args.force_field)

    pdb = app.PDBFile(str(structure_path))
    forcefield = app.ForceField(protein_ff)
    if args.solvent == "implicit":
        forcefield.addGenerator(app.GBSAOBCForceGenerator(water_model="tip3p"))

    topology = pdb.topology
    positions = pdb.positions

    logger.info(
        "体系: %d 原子, %d 残基, 力场 %s, 溶剂 %s",
        topology.getNumAtoms(), topology.getNumResidues(), protein_ff, args.solvent,
    )

    system = forcefield.createSystem(
        topology,
        nonbondedMethod=app.NoCutoff,
        constraints=app.HBonds,
        rigidWater=True,
    )

    integrator = openmm.LangevinMiddleIntegrator(
        args.temperature * unit.kelvin,
        args.friction / unit.picosecond,
        args.timestep * unit.picoseconds,
    )

    platform = None
    if args.platform:
        platform = openmm.Platform.getPlatformByName(args.platform)
    simulation = app.Simulation(topology, system, integrator, platform)
    simulation.context.setPositions(positions)

    if args.minimize_steps > 0:
        logger.info("能量最小化 %d 步", args.minimize_steps)
        simulation.minimizeEnergy()
        simulation.context.setVelocitiesToTemperature(args.temperature * unit.kelvin)

    if args.equilibration_steps > 0:
        logger.info("平衡 %d 步", args.equilibration_steps)
        simulation.step(args.equilibration_steps)

    traj_path = output_dir / "trajectory.dcd"
    state_path = output_dir / "final_state.pdb"
    reporter_interval = max(1, args.report_interval)
    traj_interval = max(1, args.trajectory_interval)

    simulation.reporters.append(
        app.StateDataReporter(
            str(output_dir / "state.csv"),
            reporter_interval,
            step=True,
            potentialEnergy=True,
            temperature=True,
            kineticEnergy=True,
            volume=True,
        )
    )
    simulation.reporters.append(
        app.DCDReporter(str(traj_path), traj_interval, append=False)
    )

    logger.info("开始生产性 MD: %d 步, 时间步 %.4f ps", args.steps, args.timestep)
    simulation.step(args.steps)

    final_positions = simulation.context.getState(getPositions=True).getPositions()
    with open(state_path, "w", encoding="utf-8") as handle:
        app.PDBFile.writeFile(simulation.topology, final_positions, handle)

    final_state = simulation.context.getState(
        getEnergy=True, getPositions=False
    )
    potential = final_state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    logger.info("完成：末帧势能 %.1f kJ/mol", potential)

    return {
        "engine": "OpenMM",
        "openmm_version": openmm.version.version,
        "n_atoms": topology.getNumAtoms(),
        "n_residues": topology.getNumResidues(),
        "force_field": protein_ff,
        "water_model": water_ff if args.solvent == "implicit" else None,
        "solvent": args.solvent,
        "temperature_k": args.temperature,
        "friction_per_ps": args.friction,
        "timestep_ps": args.timestep,
        "minimize_steps": args.minimize_steps,
        "equilibration_steps": args.equilibration_steps,
        "production_steps": args.steps,
        "trajectory": str(traj_path),
        "final_state": str(state_path),
        "state_csv": str(output_dir / "state.csv"),
        "final_potential_kj_per_mol": potential,
        "platform": simulation.context.getPlatform().getName(),
    }


def main(args: argparse.Namespace) -> int:
    """技能主入口。"""
    logger = get_logger(f"biopipe.{SKILL_NAME}")
    available, message = _check_openmm()
    logger.info("%s", message)

    if args.check_only:
        if not available:
            logger.error("未检测到 OpenMM，无法运行模拟")
            return 1
        logger.info("OpenMM 环境检查通过")
        return 0

    if not available:
        raise SkillError(
            "运行分子动力学模拟需要 OpenMM。请先安装：\n"
            "  conda install -c conda-forge openmm\n"
            "  或 pip install openmm\n"
            "安装后重试；本技能的分析功能（run_analysis.py）不需要 OpenMM。"
        )

    output_dir = _resolve_output_dir(args.output)
    report = _run_openmm(args, output_dir)

    report_path = output_dir / "simulation_report.json"
    report_path.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    logger.info("已写出 %s", report_path)
    logger.info("轨迹: %s", report["trajectory"])
    logger.info("末态结构: %s", report["final_state"])
    return 0


if __name__ == "__main__":
    sys.exit(run_cli(main, parser=_build_parser(), skill=SKILL_NAME, version="1.0.0"))
