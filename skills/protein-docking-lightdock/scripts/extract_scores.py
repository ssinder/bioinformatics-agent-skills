#!/usr/bin/env python3
"""汇总 LightDock 对接结果的分数。

对应老仓库三个脚本里重复了三遍的 ``extract_best_score`` 与
``generate_scores_csv``。

老实现有两个真实缺陷
--------------------
1. **只取最后一个 GSO 文件**（``sorted(glob("gso_*.out"))[-1]``）。LightDock
   每 10 步存一次结果，文件名不补零，所以字典序下 ``gso_100.out`` 排在
   ``gso_20.out`` 前面——100 步的收敛结果会被 20 步的覆盖掉。
2. **glowworm 编号取错字段**。老正则捕获的是 ``LigID``（配体结构索引，
   单体对接时恒为 0），却当成 glowworm 编号使用。

本脚本按**步数数值**取最新 GSO 文件，并按**行号**确定 glowworm 编号。

Examples:
    # 汇总单个作业
    python extract_scores.py --job-dir results/rec__lig/

    # 汇总整棵目录树
    python extract_scores.py --jobs-root results/ --output-dir scores/

    # 只看每个 swarm 的最优（不展开全部 glowworm）
    python extract_scores.py --jobs-root results/ --per-swarm
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Optional

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "packages"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from biopipe_core import ensure_dir, get_logger  # noqa: E402
from lightdock_client import collect_swarm_results, find_global_best  # noqa: E402

LOGGER = get_logger("biopipe.protein-docking-lightdock")

GLOBAL_FIELDS = [
    "job_name",
    "swarm",
    "glowworm",
    "score",
    "luciferin",
    "gso_file",
    "num_glowworms",
]

PER_SWARM_FIELDS = [
    "job_name",
    "swarm",
    "best_score",
    "best_glowworm",
    "num_glowworms",
    "gso_file",
]


def build_parser() -> argparse.ArgumentParser:
    """构造命令行参数解析器。"""
    parser = argparse.ArgumentParser(
        prog="extract_scores",
        description="汇总 LightDock 对接结果的分数",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "示例:\n"
            "  python extract_scores.py --job-dir results/rec__lig/\n"
            "  python extract_scores.py --jobs-root results/ --output-dir scores/\n"
            "  python extract_scores.py --jobs-root results/ --per-swarm\n"
        ),
    )
    parser.add_argument("--job-dir", help="单个对接作业的输出目录")
    parser.add_argument("--jobs-root", help="批量处理的根目录")
    parser.add_argument(
        "--output-dir", default="", help="输出目录（缺省放在作业目录里）"
    )
    parser.add_argument(
        "--per-swarm",
        action="store_true",
        help="每个 swarm 输出一行（而不是每个作业一行）",
    )
    return parser


def collect_job_dirs(args: argparse.Namespace) -> list[Path]:
    """收集要汇总的作业目录（判据：内部含有 swarm_* 子目录）。"""
    if args.job_dir:
        return [Path(args.job_dir)]
    if args.jobs_root:
        root = Path(args.jobs_root)
        if not root.is_dir():
            return []
        return sorted(p for p in root.iterdir() if p.is_dir() and any(p.glob("swarm_*")))
    return []


def summarize_job(job_dir: Path) -> list[dict]:
    """汇总单个作业。

    Args:
        job_dir: 作业输出目录。

    Returns:
        每个 swarm 一行的字典列表。
    """
    rows: list[dict] = []
    for swarm in collect_swarm_results(job_dir):
        if swarm.best is None:
            continue
        rows.append(
            {
                "job_name": job_dir.name,
                "swarm": swarm.swarm,
                "best_score": f"{swarm.best.scoring:.4f}",
                "best_glowworm": swarm.best.glowworm,
                "num_glowworms": len(swarm.entries),
                "gso_file": str(swarm.gso_file),
            }
        )
    return rows


def global_row(job_dir: Path) -> Optional[dict]:
    """取作业的全局最优一行。"""
    best = find_global_best(job_dir)
    if best is None:
        return None
    swarm, entry = best
    return {
        "job_name": job_dir.name,
        "swarm": swarm.swarm,
        "glowworm": entry.glowworm,
        "score": f"{entry.scoring:.4f}",
        "luciferin": f"{entry.luciferin:.4f}",
        "gso_file": str(swarm.gso_file),
        "num_glowworms": len(swarm.entries),
    }


def main(argv: Optional[list[str]] = None) -> int:
    """命令行主入口。"""
    args = build_parser().parse_args(argv)

    job_dirs = collect_job_dirs(args)
    if not job_dirs:
        LOGGER.error("必须提供 --job-dir 或 --jobs-root")
        return 2

    if args.output_dir:
        report_dir = ensure_dir(args.output_dir)
    elif len(job_dirs) == 1:
        report_dir = job_dirs[0]
    else:
        report_dir = ensure_dir(job_dirs[0].parent / "scores")

    if args.per_swarm:
        rows: list[dict] = []
        for job_dir in job_dirs:
            rows.extend(summarize_job(job_dir))
        fieldnames = PER_SWARM_FIELDS
    else:
        rows = []
        for job_dir in job_dirs:
            row = global_row(job_dir)
            if row is not None:
                rows.append(row)
            else:
                LOGGER.warning("[%s] 没有有效得分", job_dir.name)
        fieldnames = GLOBAL_FIELDS

    rows.sort(key=lambda r: float(r.get("score") or r.get("best_score") or 0.0))

    csv_path = report_dir / "docking_scores.csv"
    with csv_path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    json_path = report_dir / "docking_scores.json"
    json_path.write_text(
        json.dumps(rows, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    LOGGER.info("已汇总 %d 行 -> %s", len(rows), csv_path)
    if rows and not args.per_swarm:
        LOGGER.info("最优: %s（%s）", rows[0]["job_name"], rows[0]["score"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
