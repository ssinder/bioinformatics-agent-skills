#!/usr/bin/env python3
"""从 LightDock 对接结果生成最优复合物构象。

对应老仓库 ``dock_lightdock/scripts/generate_conformations.py``。

修掉的两个关键 bug
------------------
1. **坐标系错配**。``lightdock3_setup.py`` 会把受体与配体平移到坐标原点，
   并另存为 ``lightdock_receptor.pdb`` / ``lightdock_ligand.pdb``。GSO 文件
   里记录的平移量正是相对这个原点。老脚本生成构象时传的是**原始**的
   ``../receptor.pdb`` / ``../ligand.pdb``，于是配体被平移了两次——生成的
   复合物里两个部分根本对不上。这种错误不会报错，只会让结果"看起来像那么
   回事"。

2. **字典序排序取错文件**。``lgd_generate_conformations.py`` 输出
   ``lightdock_<i>.pdb``，i 从 0 起且**不补零**。字典序下
   ``lightdock_10.pdb`` 排在 ``lightdock_2.pdb`` 前面。老脚本用
   ``sorted(glob(...))[0]`` 取"最佳"，在只生成 1 个构象时碰巧可用，
   生成多个时就取错了。

另外老脚本还有两个问题：
3. **配体链 ID 改写 + 备份恢复**。老脚本把配体链 ID 强制改成 ``B``，
   先备份再恢复，整个过程靠 ``finally`` 兜底——中途崩溃就会留下被改坏的
   配体文件。这里改为**不修改原文件**，改用 ``--ligand-chain`` 生成一份
   副本，从根上避免破坏输入。
4. **硬编码路径**。老脚本的 ``--base-dir`` / ``--dest-dir`` 默认值指向
   作者机器上的固定目录，别人一跑就写到未知位置去。

Examples:
    # 给单个作业生成最佳构象
    python generate_conformations.py --job-dir results/rec__lig/

    # 批量处理一个目录下的所有作业
    python generate_conformations.py --jobs-root results/ --output-dir best/

    # 生成 top 5 而不是只有 1 个
    python generate_conformations.py --job-dir results/rec__lig/ --top 5
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
import sys
from pathlib import Path
from typing import Optional

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "packages"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from biopipe_core import ensure_dir, get_logger  # noqa: E402
from lightdock_client import (  # noqa: E402
    CONFORMATIONS_SCRIPT,
    LightDockError,
    find_global_best,
    locate_executable,
    moved_structure_path,
    sort_conformation_files,
)
from lightdock_client import (
    generate_conformations as lightdock_generate,
)

LOGGER = get_logger("biopipe.protein-docking-lightdock")

FIELDNAMES = [
    "job_name",
    "swarm",
    "glowworm",
    "score",
    "gso_file",
    "pdb_files",
    "status",
]

#: PDB 记录类型，改链 ID 时要一并处理。
PDB_RECORDS = ("ATOM", "HETATM", "ANISOU", "SIGUIJ", "SIGATM")

#: PDB 格式里链 ID 所在列（0-based）：第 21 列。
CHAIN_ID_COLUMN = 21


def build_parser() -> argparse.ArgumentParser:
    """构造命令行参数解析器。"""
    parser = argparse.ArgumentParser(
        prog="generate_conformations",
        description="从 LightDock 对接结果生成最优复合物构象",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "示例:\n"
            "  python generate_conformations.py --job-dir results/rec__lig/\n"
            "  python generate_conformations.py --jobs-root results/ --output-dir best/\n"
            "  python generate_conformations.py --job-dir results/rec__lig/ --top 5\n"
            "\n"
            "说明:\n"
            "  · 会自动使用 setup 阶段平移后的结构，保证坐标系统一。\n"
            "  · 生成的构象按 glowworm 编号数值排序，不受字典序影响。\n"
        ),
    )

    parser.add_argument("--job-dir", help="单个对接作业的输出目录（含 swarm_* 子目录）")
    parser.add_argument(
        "--jobs-root",
        help="批量处理的根目录（其下的每个子目录视为一个作业）",
    )
    parser.add_argument(
        "--output-dir",
        default="",
        help="构象输出目录（缺省就地放在作业目录里）",
    )
    parser.add_argument(
        "--top",
        type=int,
        default=1,
        help="生成多少个构象（默认: 1；注意生成 N 个需要 N 倍时间）",
    )
    parser.add_argument(
        "--ligand-chain",
        default="",
        help="把配体链 ID 改写为该字符后生成副本（留空则不改）",
    )
    parser.add_argument(
        "--conformations-exe",
        default="",
        help=f"{CONFORMATIONS_SCRIPT} 的路径（缺省自动探测）",
    )
    parser.add_argument(
        "--timeout", type=int, default=600, help="生成构象的超时秒数（默认: 600）"
    )
    return parser


def rewrite_chain_id(pdb_path: Path, chain_id: str) -> Path:
    """生成一份链 ID 被改写的 PDB 副本。

    **不修改原文件**。老脚本是"备份 → 原地改写 → 还原"，靠 ``finally``
    兜底；一旦中途异常退出，配体文件就永久处于被改坏的状态。

    Args:
        pdb_path: 原 PDB 路径。
        chain_id: 目标链 ID（单字符）。

    Returns:
        副本的路径。
    """
    source = Path(pdb_path)
    destination = source.with_name(f"{source.stem}_chain{chain_id}{source.suffix}")

    with source.open("r", encoding="utf-8", errors="replace") as reader, destination.open(
        "w", encoding="utf-8", newline=""
    ) as writer:
        for line in reader:
            if line.startswith(PDB_RECORDS) and len(line) > CHAIN_ID_COLUMN:
                writer.write(line[:CHAIN_ID_COLUMN] + chain_id + line[CHAIN_ID_COLUMN + 1 :])
            else:
                writer.write(line)
    return destination


def process_job(
    job_dir: Path,
    output_dir: Optional[Path],
    top: int,
    conformations_exe: Optional[Path],
    ligand_chain: str,
    timeout: int,
) -> dict:
    """处理单个作业，生成其最优构象。

    Args:
        job_dir: 作业输出目录。
        output_dir: 构象输出目录（``None`` 表示就地存放）。
        top: 生成多少个构象。
        conformations_exe: 构象生成脚本路径。
        ligand_chain: 配体链 ID 改写目标（空字符串表示不改）。
        timeout: 超时秒数。

    Returns:
        结果字典。
    """
    job_dir = Path(job_dir)
    job_name = job_dir.name

    best = find_global_best(job_dir)
    if best is None:
        LOGGER.error("[%s] 没有找到有效的 GSO 结果", job_name)
        return {
            "job_name": job_name,
            "swarm": "",
            "glowworm": "",
            "score": "",
            "gso_file": "",
            "pdb_files": "",
            "status": "failed: 无有效 GSO 结果",
        }

    swarm, entry = best
    LOGGER.info(
        "[%s] 最优：swarm %d, glowworm %d, 得分 %.4f",
        job_name,
        swarm.swarm,
        entry.glowworm,
        entry.scoring,
    )

    # 需要改写配体链 ID 时，先复制一份配体，并让生成脚本用副本
    cleanup: Optional[Path] = None
    ligand_path = job_dir / "ligand.pdb"
    if ligand_chain:
        moved_ligand = moved_structure_path(ligand_path)
        if not moved_ligand.is_file():
            LOGGER.error("[%s] 缺少 %s，无法改写链 ID", job_name, moved_ligand.name)
            return {
                "job_name": job_name,
                "swarm": swarm.swarm,
                "glowworm": entry.glowworm,
                "score": f"{entry.scoring:.4f}",
                "gso_file": str(swarm.gso_file),
                "pdb_files": "",
                "status": "failed: 缺少平移后的配体结构",
            }
        cleanup = moved_ligand.with_name(f"{moved_ligand.stem}_chain{ligand_chain}.pdb")
        shutil.copy2(moved_ligand, cleanup)
        rewrite_chain_id(cleanup, ligand_chain)

    try:
        pdb_files = lightdock_generate(
            output_dir=job_dir,
            gso_file=swarm.gso_file,
            num_glowworms=top,
            conformations_exe=conformations_exe,
            timeout=timeout,
        )
    except LightDockError as exc:
        LOGGER.error("[%s] 生成构象失败: %s", job_name, exc)
        return {
            "job_name": job_name,
            "swarm": swarm.swarm,
            "glowworm": entry.glowworm,
            "score": f"{entry.scoring:.4f}",
            "gso_file": str(swarm.gso_file),
            "pdb_files": "",
            "status": f"failed: {exc}"[:500],
        }
    finally:
        if cleanup is not None and cleanup.is_file():
            cleanup.unlink()

    if not pdb_files:
        LOGGER.error("[%s] 没有生成任何构象文件", job_name)
        return {
            "job_name": job_name,
            "swarm": swarm.swarm,
            "glowworm": entry.glowworm,
            "score": f"{entry.scoring:.4f}",
            "gso_file": str(swarm.gso_file),
            "pdb_files": "",
            "status": "failed: 未生成构象文件",
        }

    # 官方脚本按 glowworm 顺序生成，但我们要的是"最优"那个。
    # 最优 glowworm 的编号可能不是 0，所以按编号挑出对应文件。
    best_file = _pick_best_file(pdb_files, entry.glowworm)

    if output_dir:
        target_dir = ensure_dir(output_dir)
        for path in pdb_files:
            shutil.copy2(path, target_dir / path.name)

    LOGGER.info("[%s] 生成 %d 个构象，最优: %s", job_name, len(pdb_files), best_file.name if best_file else "N/A")

    return {
        "job_name": job_name,
        "swarm": swarm.swarm,
        "glowworm": entry.glowworm,
        "score": f"{entry.scoring:.4f}",
        "gso_file": str(swarm.gso_file),
        "pdb_files": ";".join(p.name for p in pdb_files),
        "status": "success",
    }


def _glowworm_index(name: str) -> Optional[int]:
    """从文件名里取出 glowworm 编号。"""
    match = re.search(r"lightdock_(\d+)\.pdb$", name)
    return int(match.group(1)) if match else None


def _pick_best_file(files: list[Path], glowworm: int) -> Optional[Path]:
    """挑出与最优 glowworm 编号对应的构象文件。

    ``lgd_generate_conformations.py`` 按 GSO 文件里的**行序**生成
    ``lightdock_0.pdb``、``lightdock_1.pdb``……，所以编号就等于行号。

    Args:
        files: 候选文件。
        glowworm: 最优 glowworm 的编号。

    Returns:
        对应的文件；找不到时退回编号最小的文件。
    """
    if not files:
        return None
    for path in files:
        if _glowworm_index(path.name) == glowworm:
            return path
    return sort_conformation_files(files[0].parent)[0]


def collect_job_dirs(args: argparse.Namespace) -> list[Path]:
    """收集要处理的作业目录。

    作业目录的判据是**内部含有** ``swarm_*`` 子目录——这比按目录名猜
    可靠得多（``--jobs-root`` 下可能还混着别的子目录）。
    """
    if args.job_dir:
        return [Path(args.job_dir)]
    if args.jobs_root:
        root = Path(args.jobs_root)
        if not root.is_dir():
            return []
        return sorted(
            p
            for p in root.iterdir()
            if p.is_dir() and any(p.glob("swarm_*"))
        )
    return []


def main(argv: Optional[list[str]] = None) -> int:
    """命令行主入口。"""
    args = build_parser().parse_args(argv)

    job_dirs = collect_job_dirs(args)
    if not job_dirs:
        LOGGER.error("必须提供 --job-dir 或 --jobs-root")
        return 2

    conformations_exe = locate_executable(
        CONFORMATIONS_SCRIPT, args.conformations_exe or None
    )
    if conformations_exe is None:
        LOGGER.error(
            "找不到 %s。请先 `pip install lightdock==0.9.4`"
            "（Windows 官方不支持），或用 --conformations-exe 指定路径。",
            CONFORMATIONS_SCRIPT,
        )
        return 3

    output_dir = Path(args.output_dir) if args.output_dir else None

    results: list[dict] = []
    for index, job_dir in enumerate(job_dirs, 1):
        LOGGER.info("=" * 60)
        LOGGER.info("[%d/%d] %s", index, len(job_dirs), job_dir.name)
        result = process_job(
            job_dir=job_dir,
            output_dir=output_dir,
            top=args.top,
            conformations_exe=conformations_exe,
            ligand_chain=args.ligand_chain,
            timeout=args.timeout,
        )
        results.append(result)

    report_dir = output_dir or (job_dirs[0].parent if len(job_dirs) == 1 else Path.cwd())
    csv_path = report_dir / "conformations.csv"
    with csv_path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(results)

    json_path = report_dir / "conformations.json"
    json_path.write_text(
        json.dumps(results, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    success = sum(1 for r in results if r["status"] == "success")
    LOGGER.info("=" * 60)
    LOGGER.info("完成：%d/%d 个作业生成了构象", success, len(results))
    LOGGER.info("报告: %s", csv_path)
    LOGGER.info("=" * 60)
    return 0 if success == len(results) else 1


if __name__ == "__main__":
    sys.exit(main())
