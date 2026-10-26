#!/usr/bin/env python3
"""批量运行 LightDock 蛋白质对接。

对应老仓库 ``dock_lightdock/scripts/dock_completess.py`` 与
``dock_simple.py``——这两个脚本除了参数不同，逻辑几乎逐字重复，
本脚本把两者合并成一份。

因为 LightDock 是**本地程序**，一次批量对接可能跑几小时，所以本脚本：

- 每完成一个作业就**增量写一行**到 CSV（老脚本每 500 个作业才写一次，
  中途 Ctrl-C 会丢掉前面全部结果）；
- 支持 ``--resume``，跳过已有有效结果的作业；
- 单个作业失败**不影响**其余作业，老脚本在取不到分数时会中止整个批次。

Examples:
    # 一个受体 vs 一组配体
    python run_docking.py --receptor rec.pdb --ligands ligands/*.pdb \\
        --output-dir results/

    # 多个受体 vs 多个配体（笛卡尔积）
    python run_docking.py --receptors recs/*.pdb --ligands ligands/*.pdb \\
        --output-dir results/

    # 快速试跑：少 swarm、少步数
    python run_docking.py --receptor rec.pdb --ligands a.pdb b.pdb \\
        --swarms 5 --glowworms 10 --steps 20 --output-dir test/

    # 断点续跑
    python run_docking.py --receptor rec.pdb --ligands ligands/*.pdb \\
        --output-dir results/ --resume
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from collections.abc import Sequence
from pathlib import Path
from typing import Optional

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "packages"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from biopipe_core import ensure_dir, get_logger  # noqa: E402
from lightdock_client import (  # noqa: E402
    DockingJob,
    LightDockError,
    find_global_best,
    resolve_executables,
    setup_docking,
)
from lightdock_client import (
    run_docking as lightdock_run,
)

LOGGER = get_logger("biopipe.protein-docking-lightdock")

FIELDNAMES = [
    "job_id",
    "receptor",
    "ligand",
    "status",
    "best_score",
    "best_swarm",
    "best_glowworm",
    "output_dir",
    "message",
]


def build_parser() -> argparse.ArgumentParser:
    """构造命令行参数解析器。"""
    parser = argparse.ArgumentParser(
        prog="run_docking",
        description="批量运行 LightDock 本地蛋白质对接",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "示例:\n"
            "  python run_docking.py --receptor rec.pdb --ligands ligands/*.pdb \\\n"
            "      --output-dir results/\n"
            "  python run_docking.py --receptors recs/*.pdb --ligands ligs/*.pdb \\\n"
            "      --output-dir results/\n"
            "  python run_docking.py --receptor rec.pdb --ligands a.pdb \\\n"
            "      --swarms 5 --glowworms 10 --steps 20 --output-dir test/\n"
            "\n"
            "说明:\n"
            "  · LightDock 是本地程序，计算量取决于 CPU 核数与步数。\n"
            "  · 结果每完成一个作业就增量写入 CSV，中断也不丢。\n"
            "  · 用 --resume 可以跳过已经算完的作业。\n"
        ),
    )

    parser.add_argument("--receptor", help="受体 PDB 文件路径")
    parser.add_argument(
        "--receptors",
        nargs="+",
        help="多个受体 PDB 文件路径（与 --ligands 组成笛卡尔积）",
    )
    parser.add_argument(
        "--ligands",
        nargs="+",
        required=True,
        help="配体 PDB 文件路径列表（可传 glob 展开后的多个路径）",
    )
    parser.add_argument(
        "--output-dir", default="lightdock_results", help="输出目录（默认: lightdock_results）"
    )
    parser.add_argument(
        "--swarms", type=int, default=50, help="swarm 数量（默认: 50）"
    )
    parser.add_argument(
        "--glowworms", type=int, default=50, help="每个 swarm 的 glowworm 数量（默认: 50）"
    )
    parser.add_argument("--steps", type=int, default=50, help="模拟步数（默认: 50）")
    parser.add_argument(
        "--cores",
        type=int,
        default=4,
        help="并行 CPU 核数（默认: 4；LightDock 会按核数开进程）",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="跳过已有有效结果的作业",
    )
    parser.add_argument(
        "--keep-hydrogens",
        action="store_true",
        help="保留氢原子（默认移除）",
    )
    parser.add_argument(
        "--keep-oxt",
        action="store_true",
        help="保留 OXT 原子（默认移除）",
    )
    parser.add_argument(
        "--keep-water",
        action="store_true",
        help="保留水分子（默认移除）",
    )
    parser.add_argument(
        "--setup-timeout",
        type=int,
        default=600,
        help="准备阶段超时秒数（默认: 600）",
    )
    parser.add_argument(
        "--run-timeout",
        type=int,
        default=86400,
        help="模拟阶段超时秒数（默认: 86400，即一天）",
    )
    parser.add_argument(
        "--setup-exe", default="", help="lightdock3_setup.py 的路径（缺省自动探测）"
    )
    parser.add_argument(
        "--run-exe", default="", help="lightdock3.py 的路径（缺省自动探测）"
    )
    return parser


def collect_receptors(args: argparse.Namespace) -> list[Path]:
    """收集受体文件列表。"""
    paths: list[Path] = []
    if args.receptors:
        paths.extend(Path(p) for p in args.receptors)
    if args.receptor:
        paths.append(Path(args.receptor))
    return paths


def expand_ligands(items: Sequence[str]) -> list[Path]:
    """展开配体参数。

    支持三种写法：glob 模式（``ligands/*.pdb``）、目录、以及直接的文件路径。
    """
    results: list[Path] = []
    for item in items:
        candidate = Path(item)
        if candidate.is_dir():
            results.extend(sorted(candidate.glob("*.pdb")))
            continue
        if any(ch in item for ch in "*?["):
            results.extend(sorted(Path().glob(item)))
            continue
        results.append(candidate)
    return results


def write_row(csv_path: Path, json_path: Path, job: DockingJob, append: bool) -> None:
    """把一条作业记录写入 CSV 与 JSON。

    CSV 是**追加**模式且随时 flush——这样 Ctrl-C 也不会丢掉已完成的结果。
    JSON 则整体重写（作业量通常不大，重写代价可接受）。

    Args:
        csv_path: CSV 路径。
        json_path: JSON 路径。
        job: 作业记录。
        append: CSV 是否需要写表头（首次写时为 ``True``）。
    """
    row = job.to_row()
    with csv_path.open("a", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        if append:
            writer.writeheader()
        writer.writerow(row)
        handle.flush()

    existing: list[dict] = []
    if json_path.is_file():
        try:
            existing = json.loads(json_path.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError):
            existing = []
    existing.append(row)
    json_path.write_text(json.dumps(existing, ensure_ascii=False, indent=2), encoding="utf-8")


def main(argv: Optional[list[str]] = None) -> int:
    """命令行主入口。"""
    args = build_parser().parse_args(argv)

    receptors = collect_receptors(args)
    ligands = expand_ligands(args.ligands)

    if not receptors:
        LOGGER.error("必须提供 --receptor 或 --receptors")
        return 2
    if not ligands:
        LOGGER.error("没有找到任何配体文件")
        return 2

    missing = [p for p in receptors + ligands if not p.is_file()]
    if missing:
        for path in missing:
            LOGGER.error("文件不存在: %s", path)
        return 2

    try:
        executables = resolve_executables(
            setup=args.setup_exe or None, run=args.run_exe or None
        )
    except LightDockError as exc:
        LOGGER.error("%s", exc)
        return 3

    LOGGER.info("=" * 60)
    LOGGER.info("LightDock 批量对接")
    LOGGER.info("=" * 60)
    LOGGER.info("受体 %d 个，配体 %d 个，共 %d 个作业", len(receptors), len(ligands), len(receptors) * len(ligands))
    LOGGER.info("参数: swarms=%d glowworms=%d steps=%d cores=%d",
                args.swarms, args.glowworms, args.steps, args.cores)
    LOGGER.info("LightDock 脚本目录: %s", executables)

    output_dir = ensure_dir(args.output_dir)
    csv_path = output_dir / "docking_scores.csv"
    json_path = output_dir / "docking_scores.json"
    first_write = not csv_path.exists()

    stats = {"success": 0, "failed": 0, "skipped": 0}
    counter = 0
    started = time.time()

    for receptor in receptors:
        for ligand in ligands:
            counter += 1
            job_id = f"{receptor.stem}__{ligand.stem}"
            job_dir = output_dir / job_id

            LOGGER.info("-" * 60)
            LOGGER.info("[%d/%d] %s", counter, len(receptors) * len(ligands), job_id)

            if args.resume and job_dir.is_dir():
                existing = find_global_best(job_dir)
                if existing is not None:
                    swarm, entry = existing
                    LOGGER.info("已存在有效结果，跳过（得分 %.4f）", entry.scoring)
                    job = DockingJob(
                        job_id=job_id,
                        receptor=receptor,
                        ligand=ligand,
                        output_dir=job_dir,
                        status="skipped",
                        best_score=entry.scoring,
                        best_swarm=swarm.swarm,
                        best_glowworm=entry.glowworm,
                        message="跳过（--resume，已有结果）",
                    )
                    write_row(csv_path, json_path, job, first_write)
                    first_write = False
                    stats["skipped"] += 1
                    continue

            try:
                setup_docking(
                    receptor=receptor,
                    ligand=ligand,
                    output_dir=job_dir,
                    setup_exe=executables["lightdock3_setup.py"],
                    swarms=args.swarms,
                    glowworms=args.glowworms,
                    timeout=args.setup_timeout,
                    no_hydrogens=not args.keep_hydrogens,
                    no_oxt=not args.keep_oxt,
                    no_water=not args.keep_water,
                )
                LOGGER.info("准备阶段完成，开始模拟（这一步最慢）")

                lightdock_run(
                    setup_json=job_dir / "setup.json",
                    output_dir=job_dir,
                    run_exe=executables["lightdock3.py"],
                    steps=args.steps,
                    cores=args.cores,
                    timeout=args.run_timeout,
                )
            except LightDockError as exc:
                LOGGER.error("作业失败: %s", exc)
                job = DockingJob(
                    job_id=job_id,
                    receptor=receptor,
                    ligand=ligand,
                    output_dir=job_dir,
                    status="failed",
                    message=str(exc)[:500],
                )
                write_row(csv_path, json_path, job, first_write)
                first_write = False
                stats["failed"] += 1
                continue

            best = find_global_best(job_dir)
            if best is None:
                LOGGER.error("模拟结束但没有解析到任何有效得分")
                job = DockingJob(
                    job_id=job_id,
                    receptor=receptor,
                    ligand=ligand,
                    output_dir=job_dir,
                    status="failed",
                    message="无法从 GSO 输出提取得分",
                )
                write_row(csv_path, json_path, job, first_write)
                first_write = False
                stats["failed"] += 1
                continue

            swarm, entry = best
            LOGGER.info(
                "完成：最优得分 %.4f（swarm %d, glowworm %d）",
                entry.scoring,
                swarm.swarm,
                entry.glowworm,
            )
            job = DockingJob(
                job_id=job_id,
                receptor=receptor,
                ligand=ligand,
                output_dir=job_dir,
                status="success",
                best_score=entry.scoring,
                best_swarm=swarm.swarm,
                best_glowworm=entry.glowworm,
            )
            write_row(csv_path, json_path, job, first_write)
            first_write = False
            stats["success"] += 1

    duration = time.time() - started
    LOGGER.info("=" * 60)
    LOGGER.info(
        "全部完成：成功 %d，失败 %d，跳过 %d，耗时 %s",
        stats["success"],
        stats["failed"],
        stats["skipped"],
        f"{duration:.0f}s",
    )
    LOGGER.info("结果文件: %s", csv_path)
    LOGGER.info("=" * 60)
    return 0 if stats["failed"] == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
