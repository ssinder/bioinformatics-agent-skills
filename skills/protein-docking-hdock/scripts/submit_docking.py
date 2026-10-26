#!/usr/bin/env python3
"""向 HDOCK Server 提交蛋白质对接作业。

HDOCK 服务端计算通常需要 5–30 分钟，所以本命令**只提交、不等待**。
提交后会立刻得到一个作业 URL，把它保存下来，之后用
``fetch_docking_results.py`` 取结果。

对应老仓库 ``dock_hdock/scripts/upload hdock（two pdb）.py`` 与
``upload hdock （p and d）.py``——那两个脚本用 Selenium 起浏览器填表，
并且用线程池并发提交。HDOCK 是公益服务器，密集提交既容易被限流，也会
拖慢所有人的作业，所以这里默认串行、每次提交之间留出间隔。

Examples:
    # 用 PDB ID:ChainID 提交
    python submit_docking.py --receptor-pdb-id 1CGI:E --ligand-pdb-id 1CGI:I \\
        --job-name my_docking --output-dir results/

    # 用 FASTA 序列提交
    python submit_docking.py --receptor-seq receptor.fasta --ligand-seq ligand.fasta \\
        --output-dir results/

    # 用 PDB 文件提交并指定结合位点
    python submit_docking.py --receptor-pdb rec.pdb --ligand-pdb lig.pdb \\
        --receptor-site 45,79,102 --output-dir results/
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from pathlib import Path
from typing import Optional

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "packages"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from biopipe_core import ensure_dir, get_logger  # noqa: E402
from hdock_client import (  # noqa: E402
    LIGAND_TYPES,
    JobStatus,
    submit_docking_job,
)

LOGGER = get_logger("biopipe.protein-docking-hdock")

FIELDNAMES = [
    "job_name",
    "job_id",
    "url",
    "status",
    "docking_score",
    "confidence_score",
    "error",
]


def build_parser() -> argparse.ArgumentParser:
    """构造命令行参数解析器。"""
    parser = argparse.ArgumentParser(
        prog="submit_docking",
        description="向 HDOCK Server 提交蛋白质对接作业（提交后不等待结果）",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "示例:\n"
            "  python submit_docking.py --receptor-pdb-id 1CGI:E --ligand-pdb-id 1CGI:I \\\n"
            "      --output-dir results/\n"
            "  python submit_docking.py --receptor-seq rec.fasta --ligand-seq lig.fasta \\\n"
            "      --output-dir results/\n"
            "\n"
            "说明:\n"
            "  · 受体与配体各需一种输入：PDB 文件 / PDB ID:ChainID / FASTA 序列。\n"
            "  · FASTA 参数既接受序列文本，也接受 FASTA 文件路径。\n"
            "  · HDOCK 计算需 5-30 分钟，本命令提交完就退出。\n"
            "  · 用 fetch_docking_results.py 取结果。\n"
        ),
    )

    receptor = parser.add_argument_group("受体输入（四选一）")
    receptor.add_argument("--receptor-pdb-id", help="受体 PDB ID:ChainID，如 1CGI:E")
    receptor.add_argument("--receptor-seq", help="受体 FASTA 序列（文本或文件路径）")
    receptor.add_argument("--receptor-pdb", help="受体 PDB 文件路径")

    ligand = parser.add_argument_group("配体输入（四选一）")
    ligand.add_argument("--ligand-pdb-id", help="配体 PDB ID:ChainID，如 1CGI:I")
    ligand.add_argument("--ligand-seq", help="配体 FASTA 序列（文本或文件路径）")
    ligand.add_argument("--ligand-pdb", help="配体 PDB 文件路径")

    parser.add_argument("--job-name", default="", help="作业名（留空则用时间戳）")
    parser.add_argument(
        "--email", default="", help="通知邮箱（留空则服务端不发通知）"
    )
    parser.add_argument(
        "--ligand-type",
        default="protein",
        choices=sorted(LIGAND_TYPES),
        help="配体类型（默认: protein）",
    )
    parser.add_argument(
        "--symmetry", default="", help="对称性，如 C2 / D3；填了走对称对接入口"
    )
    parser.add_argument("--receptor-site", default="", help="受体结合位点残基号，逗号分隔")
    parser.add_argument("--ligand-site", default="", help="配体结合位点残基号，逗号分隔")
    parser.add_argument("--restraint", default="", help="残基距离限制")
    parser.add_argument(
        "--output-dir", default="hdock_results", help="输出目录（默认: hdock_results）"
    )
    parser.add_argument(
        "--timeout", type=int, default=180, help="单次请求超时秒数（默认: 180）"
    )
    return parser


def _read_sequence(value: str) -> str:
    """把 ``--*-seq`` 参数解析成 FASTA 文本。

    参数既可以是序列本身，也可以是 FASTA 文件路径。
    """
    if not value:
        return ""
    stripped = value.strip()
    if stripped.startswith(">") or "\n" in stripped:
        return stripped
    candidate = Path(stripped)
    try:
        if candidate.is_file():
            return candidate.read_text(encoding="utf-8", errors="replace")
    except OSError:  # pragma: no cover
        pass
    return stripped


def main(argv: Optional[list[str]] = None) -> int:
    """命令行主入口。"""
    args = build_parser().parse_args(argv)

    job_name = args.job_name or f"biopipe_{time.strftime('%Y%m%d_%H%M%S')}"

    receptor_seq = _read_sequence(args.receptor_seq)
    ligand_seq = _read_sequence(args.ligand_seq)

    try:
        result = submit_docking_job(
            job_name=job_name,
            email=args.email,
            ligand_type=args.ligand_type,
            symmetry=args.symmetry,
            receptor_pdb_id=args.receptor_pdb_id or "",
            ligand_pdb_id=args.ligand_pdb_id or "",
            receptor_seq=receptor_seq,
            ligand_seq=ligand_seq,
            receptor_pdb_file=args.receptor_pdb,
            ligand_pdb_file=args.ligand_pdb,
            receptor_site=args.receptor_site,
            ligand_site=args.ligand_site,
            restraint=args.restraint,
            timeout=args.timeout,
        )
    except ValueError as exc:
        LOGGER.error("%s", exc)
        return 2

    if result.status == JobStatus.FAILED:
        LOGGER.error("提交失败: %s", "; ".join(result.errors))
        return 1

    output_dir = ensure_dir(Path(args.output_dir))

    csv_path = output_dir / "jobs.csv"
    write_header = not csv_path.exists()
    with open(csv_path, "a", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        if write_header:
            writer.writeheader()
        writer.writerow(result.to_row())
    LOGGER.info("作业记录已写入: %s", csv_path)

    json_path = output_dir / "jobs.json"
    existing: list[dict] = []
    if json_path.exists():
        try:
            existing = json.loads(json_path.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError):
            existing = []
    existing.append(result.to_row())
    json_path.write_text(
        json.dumps(existing, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    LOGGER.info("作业记录已写入: %s", json_path)

    LOGGER.info("=" * 60)
    LOGGER.info("作业已提交（服务端计算中，通常需要 5-30 分钟）")
    LOGGER.info("作业名: %s", result.job_name)
    LOGGER.info("作业 ID: %s", result.job_id)
    LOGGER.info("结果 URL: %s", result.url)
    LOGGER.info("=" * 60)
    LOGGER.info("取结果：python fetch_docking_results.py --url %s", result.url)
    return 0


if __name__ == "__main__":
    sys.exit(main())
