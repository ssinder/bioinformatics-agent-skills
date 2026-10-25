#!/usr/bin/env python3
"""引用格式导出脚本。

把 ``search_literature.py`` 产出的统一 schema CSV 导出为多种参考文献格式：
RIS（BibTeX 之外的 EndNote/Zotero/Mendeley 通用格式）、BibTeX、
GB/T 7714-2015 顺序编码制、中文纯文本编号、JSON。

用法示例::

    python export_citations.py --input results/literature_results.csv \\
        --format bibtex --output results/references.bib

    python export_citations.py --input results/literature_results.csv \\
        --format ris gbt7714 plain --output-dir results/citations
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "packages"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from biopipe_core import (  # noqa: E402
    SkillError,
    cli_parser,
    get_config,
    get_logger,
    read_csv_safe,
    resolve_path,
    run_cli,
)
from literature_core import parse_record_rows, write_citations  # noqa: E402

LOGGER = get_logger("biopipe.literature-search")

EXTENSIONS = {
    "ris": ".ris",
    "bibtex": ".bib",
    "gbt7714": ".txt",
    "plain": ".txt",
    "json": ".json",
}

FORMAT_ALIASES = {
    "bib": "bibtex",
    "gbt": "gbt7714",
    "gb/t7714": "gbt7714",
    "txt": "plain",
    "text": "plain",
    "bibtex": "bibtex",
    "ris": "ris",
    "json": "json",
}


def build_parser() -> argparse.ArgumentParser:
    """构建命令行解析器。"""
    parser = cli_parser(
        prog="literature-search-export",
        description="把文献检索结果导出为 RIS / BibTeX / GB-T 7714 等引用格式",
    )
    parser.add_argument("--input", required=True, help="输入 CSV（统一 schema）")
    parser.add_argument(
        "--format",
        nargs="+",
        default=["gbt7714"],
        choices=sorted(set(EXTENSIONS) | set(FORMAT_ALIASES)),
        help="导出的引用格式，可多选",
    )
    parser.add_argument("--output", help="输出文件路径（单格式时有效）")
    parser.add_argument("--output-dir", help="输出目录（多格式时按格式分文件写出）")
    parser.add_argument("--encoding", default="utf-8-sig", help="输出编码")
    parser.add_argument("--sort-by", default="input", choices=["input", "citations", "year"])
    parser.add_argument("--limit", type=int, default=0, help="最多导出多少条")
    return parser


def main(args: argparse.Namespace) -> int:
    """主流程。"""
    config = get_config(args.config)
    input_path = resolve_path(args.input)
    if not input_path.exists():
        raise SkillError(f"输入文件不存在: {input_path}")

    rows = read_csv_safe(input_path)
    if not rows:
        LOGGER.warning("输入文件没有记录")
        return 0

    records = parse_record_rows(rows)
    LOGGER.info("载入 %d 条记录", len(records))

    if args.sort_by == "citations":
        records.sort(key=lambda r: int(r.get("citations", 0)), reverse=True)
    elif args.sort_by == "year":
        records.sort(key=lambda r: (r.get("year") or 0), reverse=True)
    if args.limit > 0:
        records = records[: args.limit]

    formats = [FORMAT_ALIASES.get(f.lower(), f.lower()) for f in args.format]
    formats = list(dict.fromkeys(formats))

    outputs: list[str] = []
    used_names: set[str] = set()
    for fmt in formats:
        if args.output and len(formats) == 1:
            out_path = resolve_path(args.output)
        else:
            base_dir = resolve_path(args.output_dir) if args.output_dir else config.output_dir
            ext = EXTENSIONS.get(fmt, ".txt")
            stem = "references"
            # 不同格式可能共用同一扩展名（如 gbt7714 与 plain 都是 .txt），避免互相覆盖
            candidate = f"{stem}{ext}"
            counter = 2
            while candidate.lower() in used_names:
                candidate = f"{stem}_{fmt}{ext}"
                counter += 1
            used_names.add(candidate.lower())
            out_path = base_dir / candidate
        out_path.parent.mkdir(parents=True, exist_ok=True)
        write_citations(records, fmt, str(out_path), encoding=args.encoding)
        LOGGER.info("已导出 %s 格式: %s", fmt, out_path)
        outputs.append(str(out_path))

    LOGGER.info("导出完成，共 %d 个文件", len(outputs))
    return 0


if __name__ == "__main__":
    sys.exit(
        run_cli(
            main,
            parser=build_parser(),
            skill="literature-search",
        )
    )
