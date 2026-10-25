#!/usr/bin/env python3
"""批量结构分析：遍历目录中的 PDB 文件并汇总统计。

与 :mod:`analyze_structure` 的单文件分析共用同一套算法（纯 Python PDB 解析 +
内置简化 DSSP），本脚本额外提供：

- 递归扫描目录、并发处理（线程池，IO 与纯计算混合场景下有效）
- 长度分布 / 链数分布 / Rg 分布等直方图统计
- 汇总 CSV 与 JSON 报告

用法示例::

    python batch_analyze.py --input ./pdb_data --output ./results --workers 4
    python batch_analyze.py --input ./pdb_data --output ./results --min-residues 50
"""

from __future__ import annotations

import argparse
import json
import statistics
import sys
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "packages"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from analyze_structure import analyze_structure  # noqa: E402
from biopipe_core import SkillError, cli_parser, get_logger, write_csv  # noqa: E402
from biopipe_core.cli import run_cli  # noqa: E402

LOGGER = get_logger("biopipe.protein-structure")
SKILL_NAME = "protein-structure"


def _build_parser() -> argparse.ArgumentParser:
    """构造命令行解析器。"""
    parser = cli_parser(
        prog=f"{SKILL_NAME}-batch",
        description="批量 PDB 结构分析：并行处理目录并输出汇总统计",
        skill=SKILL_NAME,
    )
    parser.add_argument("--input", required=True, metavar="DIR", help="包含 PDB 文件的目录")
    parser.add_argument("--output", metavar="DIR", help="输出目录（默认取配置中的 results/）")
    parser.add_argument("--pattern", default="*.pdb", help="文件匹配模式（默认 *.pdb）")
    parser.add_argument("--recursive", action="store_true", help="递归扫描子目录")
    parser.add_argument("--workers", type=int, default=4, help="并发线程数（默认 4）")
    parser.add_argument("--max-files", type=int, default=0, help="最多处理的文件数（0 = 不限）")
    parser.add_argument("--min-residues", type=int, default=0, help="跳过残基数低于此值的结构")
    parser.add_argument("--dssp-method", choices=("auto", "external", "builtin"), default="auto")
    parser.add_argument("--dssp-executable", metavar="PATH", help="DSSP 可执行文件路径")
    parser.add_argument("--hbond-energy", type=float, default=-0.5, help="内置 DSSP 氢键能量阈值")
    parser.add_argument("--query-uniprot", action="store_true", help="联网查询 UniProt 元数据")
    return parser


def _collect_files(input_dir: Path, pattern: str, recursive: bool) -> list[Path]:
    """收集 PDB 文件。"""
    if not input_dir.is_dir():
        raise SkillError(f"输入目录不存在: {input_dir}")
    if recursive:
        files = sorted(input_dir.rglob(pattern))
    else:
        files = sorted(input_dir.glob(pattern))
    return files


def _length_bucket(length: int) -> str:
    """把残基数分桶。"""
    if length < 100:
        return "<100"
    if length < 200:
        return "100-200"
    if length < 300:
        return "200-300"
    if length < 500:
        return "300-500"
    return ">=500"


def _summarize(results: list[dict]) -> dict:
    """对批量结果做汇总统计。"""
    total_residues = sum(r["structure"]["n_residues"] for r in results)
    total_atoms = sum(r["structure"]["n_atoms"] for r in results)
    n = len(results)

    length_dist = Counter(_length_bucket(r["structure"]["seq_length"]) for r in results)
    chain_dist = Counter(r["structure"]["n_chains"] for r in results)
    method_dist = Counter(r.get("dssp_method", "none") for r in results)

    helix = sum(r.get("ss_statistics", {}).get("helix_count", 0) for r in results)
    sheet = sum(r.get("ss_statistics", {}).get("sheet_count", 0) for r in results)
    coil = sum(r.get("ss_statistics", {}).get("coil_count", 0) for r in results)

    rg_values = [
        r["geometry"]["radius_of_gyration"]["rg"]
        for r in results
        if r["geometry"]["radius_of_gyration"].get("rg")
    ]

    return {
        "n_structures": n,
        "total_residues": total_residues,
        "total_atoms": total_atoms,
        "avg_residues": round(total_residues / n, 2) if n else 0.0,
        "avg_atoms": round(total_atoms / n, 2) if n else 0.0,
        "avg_rg": round(statistics.fmean(rg_values), 3) if rg_values else 0.0,
        "stdev_rg": round(statistics.pstdev(rg_values), 3) if len(rg_values) > 1 else 0.0,
        "length_distribution": dict(length_dist),
        "chains_distribution": {str(k): v for k, v in sorted(chain_dist.items())},
        "dssp_method_distribution": dict(method_dist),
        "helix_pct": round(helix / total_residues * 100, 2) if total_residues else 0.0,
        "sheet_pct": round(sheet / total_residues * 100, 2) if total_residues else 0.0,
        "coil_pct": round(coil / total_residues * 100, 2) if total_residues else 0.0,
    }


def main(args: argparse.Namespace) -> int:
    """入口主函数。"""
    input_dir = Path(args.input).expanduser()
    if not input_dir.is_absolute():
        input_dir = Path.cwd() / input_dir

    config_dir = Path(args.output).expanduser() if args.output else Path.cwd() / "results"
    if not config_dir.is_absolute():
        config_dir = Path.cwd() / config_dir
    config_dir.mkdir(parents=True, exist_ok=True)

    files = _collect_files(input_dir, args.pattern, args.recursive)
    if args.max_files and args.max_files > 0:
        files = files[: args.max_files]
    if not files:
        raise SkillError(f"在 {input_dir} 下没有匹配 {args.pattern} 的文件")
    LOGGER.info("共发现 %d 个文件，%d 个线程", len(files), max(1, args.workers))

    results: list[dict] = []
    failures: list[tuple[str, str]] = []
    workers = max(1, args.workers)

    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {
            pool.submit(
                analyze_structure,
                path,
                args.dssp_method,
                args.dssp_executable,
                args.hbond_energy,
                1,
                args.query_uniprot,
            ): path
            for path in files
        }
        for future in as_completed(futures):
            path = futures[future]
            try:
                result = future.result()
            except Exception as exc:  # noqa: BLE001
                failures.append((path.name, str(exc)))
                continue
            if args.min_residues and result["structure"]["n_residues"] < args.min_residues:
                continue
            results.append(result)

    if not results:
        raise SkillError("没有成功分析任何结构")

    summary = _summarize(results)
    summary["generated_at"] = datetime.now().isoformat(timespec="seconds")
    summary["failed_files"] = [{"file": n, "reason": r} for n, r in failures]

    json_path = config_dir / "batch_summary.json"
    json_path.write_text(
        json.dumps({"summary": summary, "structures": results}, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    LOGGER.info("已写出 %s", json_path)

    rows = []
    for item in results:
        struct = item["structure"]
        rg = item["geometry"]["radius_of_gyration"]
        stats = item.get("ss_statistics", {})
        rows.append(
            {
                "structure_id": item["structure_id"],
                "n_chains": struct["n_chains"],
                "n_residues": struct["n_residues"],
                "n_atoms": struct["n_atoms"],
                "seq_length": struct["seq_length"],
                "rg_angstrom": round(rg.get("rg", 0.0), 3),
                "rg_ca": round(rg["rg_ca"], 3) if rg.get("rg_ca") else "",
                "shape_class": item["geometry"]["shape"]["shape_class"],
                "helix_pct": stats.get("helix_pct", 0.0),
                "sheet_pct": stats.get("sheet_pct", 0.0),
                "coil_pct": stats.get("coil_pct", 0.0),
                "dssp_method": item.get("dssp_method", ""),
            }
        )
    csv_path = config_dir / "batch_results.csv"
    write_csv(rows, csv_path)
    LOGGER.info("已写出 %s", csv_path)

    LOGGER.info(
        "完成：%d 个结构，平均残基数 %.1f，平均 Rg %.2f Å，螺旋 %.1f%%",
        summary["n_structures"],
        summary["avg_residues"],
        summary["avg_rg"],
        summary["helix_pct"],
    )
    if failures:
        LOGGER.warning("失败 %d 个文件", len(failures))
    return 0


if __name__ == "__main__":
    sys.exit(run_cli(main, parser=_build_parser(), skill=SKILL_NAME, version="1.0.0"))
