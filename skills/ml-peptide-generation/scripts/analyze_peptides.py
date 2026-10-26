#!/usr/bin/env python3
"""分析 PepMLM 生成的候选肽链。

对应老仓库 ``design_pepmlm/scripts/filter_analysis.py``。那个脚本把待分析
的文件名**硬编码**成 12 个（``results_10_3.csv`` … ``results_15_10.csv``），
换一组参数就得改源码；这里改成命令行传入。

老脚本还有两个问题：

1. ``df_all = pd.concat(all_valid, ...)`` —— 若所有文件都过滤空了，
   ``all_valid`` 是空列表，``pd.concat([])`` 会抛 ``ValueError: No
   objects to concatenate``。这里在空输入时给出明确提示。
2. 只按文件名解析参数（``results_{len}_{topk}.csv``），文件名不符合约定
   就直接 ``ValueError``。这里改为从 CSV 的 ``config`` / ``peptide_length``
   / ``top_k`` **列**读取，文件名只是兜底。

本脚本**不依赖 torch / transformers**——纯 CSV 处理，因此在没有 GPU 的
机器上也能做参数扫描结果的分析。

Examples:
    python analyze_peptides.py --input results/peptides.csv --output-dir analysis/
    python analyze_peptides.py --input results/peptides.csv --max-ppl 5.0 --min-length 8
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import statistics
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Optional

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "packages"))

from biopipe_core import ensure_dir, get_logger  # noqa: E402

LOGGER = get_logger("biopipe.ml-peptide-generation")

#: 从文件名推断参数的兜底正则，如 ``results_10_3.csv``
_FILENAME_PATTERN = re.compile(r"(\d+)[_-](\d+)")


def read_rows(path: Path) -> list[dict]:
    """读取 CSV 结果文件。"""
    with open(path, encoding="utf-8-sig", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def parse_float(value: object) -> Optional[float]:
    """把单元格解析成浮点数，失败返回 ``None``。"""
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    try:
        result = float(text)
    except ValueError:
        return None
    if result != result:  # NaN
        return None
    return result


def parse_int(value: object) -> Optional[int]:
    """把单元格解析成整数，失败返回 ``None``。"""
    number = parse_float(value)
    if number is None:
        return None
    return int(number)


def is_valid_peptide(row: dict) -> bool:
    """判断一条结果是否有效（不含未知氨基酸）。

    优先读 ``has_unknown_aa`` 列；该列不存在时退回检查序列本身是否只含
    标准氨基酸。老脚本只依赖 ``Has_Unknown_AA`` 一列，列名一变就失效。
    """
    flag = str(row.get("has_unknown_aa", "")).strip().lower()
    if flag in {"true", "1", "yes"}:
        return False
    if flag in {"false", "0", "no"}:
        return True

    binder = str(row.get("binder", "")).strip()
    if not binder:
        return False
    return all(char in "ACDEFGHIKLMNPQRSTVWY" for char in binder.upper())


def infer_config(row: dict, source: Optional[Path] = None) -> tuple[int, int, str]:
    """推断一条结果的长度、top-k 与配置名。

    优先从 CSV 的列读取；列缺失时才从文件名推断。

    Args:
        row: 结果行。
        source: 来源文件路径（用于文件名兜底）。

    Returns:
        ``(peptide_length, top_k, config)``。
    """
    peptide_length = parse_int(row.get("peptide_length"))
    top_k = parse_int(row.get("top_k"))
    config = str(row.get("config", "")).strip()

    if peptide_length is None or top_k is None:
        match = _FILENAME_PATTERN.search(source.stem) if source else None
        if match:
            if peptide_length is None:
                peptide_length = int(match.group(1))
            if top_k is None:
                top_k = int(match.group(2))

    if peptide_length is None:
        peptide_length = parse_int(row.get("binder_length")) or 0
    if top_k is None:
        top_k = 0
    if not config:
        config = f"{peptide_length}aa_topk{top_k}"
    return peptide_length, top_k, config


def filter_peptides(
    rows: Sequence[dict],
    max_ppl: Optional[float] = None,
    min_length: Optional[int] = None,
    max_length: Optional[int] = None,
    drop_invalid: bool = True,
    source: Optional[Path] = None,
) -> list[dict]:
    """筛选有效肽链。

    Args:
        rows: 原始结果行。
        max_ppl: 伪困惑度上限（只保留 ≤ 该值的）。
        min_length: 肽链长度下限。
        max_length: 肽链长度上限。
        drop_invalid: 是否丢弃含非标准氨基酸的肽链。
        source: 来源文件路径（用于文件名兜底推断参数）。

    Returns:
        通过筛选的行（已补上 ``config`` / ``peptide_length`` / ``top_k``）。
    """
    kept: list[dict] = []
    for row in rows:
        if drop_invalid and not is_valid_peptide(row):
            continue

        ppl = parse_float(row.get("pseudo_perplexity"))
        if max_ppl is not None:
            if ppl is None or ppl > max_ppl:
                continue

        length = parse_int(row.get("binder_length"))
        if length is None:
            length = len(str(row.get("binder", "")).strip())
        if min_length is not None and length < min_length:
            continue
        if max_length is not None and length > max_length:
            continue

        peptide_length, top_k, config = infer_config(row, source)
        enriched = dict(row)
        enriched["peptide_length"] = peptide_length
        enriched["top_k"] = top_k
        enriched["config"] = config
        enriched["pseudo_perplexity_value"] = ppl
        kept.append(enriched)
    return kept


def summarize_by_config(rows: Sequence[dict]) -> list[dict]:
    """按参数组合统计效果。

    输出每组配置的：总条数、有效条数、有效率、平均/最小伪困惑度。
    """
    groups: dict[str, list[dict]] = {}
    for row in rows:
        groups.setdefault(str(row.get("config", "unknown")), []).append(row)

    summary: list[dict] = []
    for config, group in groups.items():
        values = [
            row["pseudo_perplexity_value"]
            for row in group
            if row.get("pseudo_perplexity_value") is not None
        ]
        lengths = [row.get("peptide_length") or 0 for row in group]
        top_ks = [row.get("top_k") or 0 for row in group]

        summary.append(
            {
                "config": config,
                "peptide_length": max(lengths) if lengths else 0,
                "top_k": max(top_ks) if top_ks else 0,
                "total": len(group),
                "mean_ppl": statistics.fmean(values) if values else "",
                "median_ppl": statistics.median(values) if values else "",
                "min_ppl": min(values) if values else "",
                "max_ppl": max(values) if values else "",
            }
        )

    summary.sort(key=lambda item: (item["mean_ppl"] == "", item["mean_ppl"]))
    return summary


def build_parser() -> argparse.ArgumentParser:
    """构造命令行参数解析器。"""
    parser = argparse.ArgumentParser(
        prog="analyze_peptides",
        description="分析 PepMLM 生成的候选肽链：过滤、排序、参数统计",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "示例:\n"
            "  python analyze_peptides.py --input results/peptides.csv --output-dir analysis/\n"
            "  python analyze_peptides.py --input results/peptides.csv --max-ppl 5.0\n"
            "\n"
            "说明:\n"
            "  本脚本只用标准库，不需要 torch / transformers。\n"
        ),
    )
    parser.add_argument(
        "--input",
        nargs="+",
        required=True,
        help="待分析的 CSV 结果文件（可多个）",
    )
    parser.add_argument(
        "--output-dir",
        default="peptide_analysis",
        help="输出目录（默认: peptide_analysis）",
    )
    parser.add_argument(
        "--max-ppl",
        type=float,
        default=None,
        help="伪困惑度上限，只保留 ≤ 该值的肽链",
    )
    parser.add_argument(
        "--min-length",
        type=int,
        default=None,
        help="肽链长度下限",
    )
    parser.add_argument(
        "--max-length",
        type=int,
        default=None,
        help="肽链长度上限",
    )
    parser.add_argument(
        "--keep-invalid",
        action="store_true",
        help="保留含 X/B/Z 等非标准氨基酸的肽链",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=20,
        help="终端打印的前 N 条（默认: 20）",
    )
    return parser


def main(argv: Optional[list[str]] = None) -> int:
    """命令行主入口。"""
    args = build_parser().parse_args(argv)

    all_rows: list[dict] = []
    for item in args.input:
        path = Path(item)
        if not path.is_file():
            LOGGER.error("文件不存在: %s", path)
            return 2
        rows = read_rows(path)
        LOGGER.info("读取 %s: %d 条", path, len(rows))
        all_rows.extend(rows)

    if not all_rows:
        LOGGER.error("没有读到任何结果行")
        return 1

    kept = filter_peptides(
        all_rows,
        max_ppl=args.max_ppl,
        min_length=args.min_length,
        max_length=args.max_length,
        drop_invalid=not args.keep_invalid,
    )

    if not kept:
        LOGGER.error(
            "筛选后没有任何肽链（输入 %d 条）。若用了 --max-ppl，"
            "试着放宽阈值；或加 --keep-invalid 保留含未知氨基酸的结果。",
            len(all_rows),
        )
        return 1

    kept.sort(key=lambda row: (row["pseudo_perplexity_value"] is None,
                               row["pseudo_perplexity_value"] if row["pseudo_perplexity_value"] is not None else 0.0))

    output_dir = ensure_dir(Path(args.output_dir))

    filtered_path = output_dir / "filtered_peptides.csv"
    fieldnames = list(kept[0].keys())
    with open(filtered_path, "w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(kept)
    LOGGER.info("筛选结果已写入: %s（%d 条）", filtered_path, len(kept))

    summary = summarize_by_config(kept)
    summary_path = output_dir / "parameter_analysis.csv"
    summary_fields = [
        "config", "peptide_length", "top_k", "total",
        "mean_ppl", "median_ppl", "min_ppl", "max_ppl",
    ]
    with open(summary_path, "w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(summary)
    LOGGER.info("参数统计已写入: %s", summary_path)

    json_path = output_dir / "parameter_analysis.json"
    json_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    LOGGER.info("参数统计已写入: %s", json_path)

    print("=" * 64)
    print(f"筛选结果: 共 {len(kept)} 条有效肽链（输入 {len(all_rows)} 条）")
    print("=" * 64)
    print(f"\nTop {min(args.top_n, len(kept))} 最佳肽链（按伪困惑度排序）:")
    for row in kept[: args.top_n]:
        print(
            f"  {row.get('config', ''):18s} {row.get('binder', ''):20s} "
            f"PPL={row.get('pseudo_perplexity_value')}"
        )

    print("\n" + "=" * 64)
    print("参数分析（按平均伪困惑度排序）:")
    print("=" * 64)
    for item in summary:
        print(
            f"  {item['config']:18s} total={item['total']:4d} "
            f"mean_ppl={item['mean_ppl']} min_ppl={item['min_ppl']}"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
