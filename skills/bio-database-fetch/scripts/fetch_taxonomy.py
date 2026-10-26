#!/usr/bin/env python3
"""按 NCBI Taxonomy ID 批量获取物种分类信息。

对应老仓库 ``crawl_database/scripts/02_ncbi_taxonomy_example.py`` 与
``03_ncbi_taxonomy_batch.py``。老脚本用 Biopython 的 ``Entrez`` 接口，
单线程查询后把结果拼成一行行文本；这里直接调 E-utilities 的
``efetch.fcgi``（Taxonomy 的 ID 本身就是主键，不需要先走 esearch），
**一次请求可取 100 个 ID**，输出结构化的 CSV/JSON。

相对老脚本的改进：

1. **email 是必填且真实的**。NCBI 要求调用方标识身份，老脚本准备了
   5 个 ``userN@example.com`` 轮流冒充——既违规，也会因为无法联系到
   调用方而被更严格地限速。这里从 ``--email`` 或 ``BIOPIPE_EMAIL``
   环境变量读取，缺失时直接报错。
2. **不再多线程写同一个文件**。老脚本用 ``ThreadPoolExecutor``，每个
   线程各自 ``open(output_tsv, "a")`` 追加一行，两个线程的写入会交错，
   行可能被撕裂成半行。这里在内存里收集结果，一次性写出。
3. **完整保留所有分类层级**。老脚本只挑 8 个标准 rank，会漏掉
   ``domain``（Eukaryota）这类信息量很大的层级。这里把
   ``LineageEx`` 里的**全部** rank 都保留下来。
4. **lineage 包含自身**。``LineageEx`` 只含祖先、不含自身，所以
   species 层级必须从顶层 ``<Rank>``/``<ScientificName>`` 补进去。

Examples:
    # 从文件读取 Taxonomy ID
    python fetch_taxonomy.py --tax-ids taxids.txt --output-dir results/

    # 直接给 ID
    python fetch_taxonomy.py --tax-ids 9606 10090 --output-dir results/

    # email 也可从环境变量读
    export BIOPIPE_EMAIL="you@your-institution.edu"
    python fetch_taxonomy.py --tax-ids 9606 --output-dir results/
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Optional

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "packages"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from biopipe_core import ensure_dir, get_logger  # noqa: E402
from database_client import ELEMENTARY_RANKS, fetch_taxonomy  # noqa: E402

LOGGER = get_logger("biopipe.bio-database-fetch")

#: 输出 CSV 的列顺序。``by_rank`` 展开成固定的几列，其余层级放在
#: ``full_lineage`` 与 ``lineage_json`` 里。
FIELDNAMES = [
    "tax_id",
    "scientific_name",
    "common_name",
    "rank",
    "superkingdom",
    "domain",
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
    "lineage_depth",
    "full_lineage",
]

#: 会展开成独立列的层级。
#:
#: 注意 ``superkingdom`` 在现代 NCBI 分类学里**已经不存在**了——2018 年
#: 之后 "cellular organisms" 的 rank 改成了 ``cellular root``，顶层大分类
#: 改用 ``domain``（Eukaryota / Bacteria / Archaea）。这里仍然保留
#: ``superkingdom`` 一列是为了兼容更老的数据快照；``domain`` 才是现在
#: 真正有值的顶层列。
FLATTENED_RANKS = (
    "superkingdom",
    "domain",
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
)


def build_parser() -> argparse.ArgumentParser:
    """构造命令行参数解析器。"""
    parser = argparse.ArgumentParser(
        prog="fetch_taxonomy",
        description="从 NCBI Taxonomy 批量获取物种分类信息（需 email）",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "示例:\n"
            "  python fetch_taxonomy.py --tax-ids taxids.txt --output-dir results/\n"
            "  python fetch_taxonomy.py --tax-ids 9606 10090 --output-dir results/\n"
            "  export BIOPIPE_EMAIL=you@example.edu\n"
            "  python fetch_taxonomy.py --tax-ids 9606 --output-dir results/\n"
            "\n"
            "说明:\n"
            "  · email 必须提供：NCBI 要求调用方标识身份。\n"
            "  · 一次请求可取 100 个 Taxonomy ID。\n"
            "  · lineage 包含祖先与自身，保留全部 rank（含 domain）。\n"
        ),
    )
    parser.add_argument(
        "--tax-ids",
        nargs="+",
        required=True,
        help="Taxonomy ID 列表，或包含 ID 的文本文件路径（每行一个）",
    )
    parser.add_argument(
        "--output-dir",
        default="taxonomy_results",
        help="结果输出目录（默认: taxonomy_results）",
    )
    parser.add_argument(
        "--email",
        default=os.environ.get("BIOPIPE_EMAIL", ""),
        help="联系邮箱（NCBI 要求）；默认取环境变量 BIOPIPE_EMAIL",
    )
    parser.add_argument(
        "--api-key",
        default=os.environ.get("NCBI_API_KEY", ""),
        help="NCBI API key；限速从 3 req/s 提到 10 req/s",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=60,
        help="单次请求超时秒数（默认: 60）",
    )
    return parser


def read_id_file(path: Path) -> list[str]:
    """从文本文件读取 Taxonomy ID 列表。

    自动跳过表头行（如果首行不是纯数字）、注释行（``#``）与空行。
    """
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    ids: list[str] = []
    for index, raw in enumerate(lines):
        stripped = raw.strip()
        if not stripped or stripped.startswith("#"):
            continue
        candidate = stripped.split()[0].split(",")[0]
        if not candidate.isdigit():
            if index == 0:
                LOGGER.info("跳过表头行: %s", stripped)
            continue
        if candidate not in ids:
            ids.append(candidate)
    return ids


def collect_tax_ids(items: Sequence[str]) -> list[str]:
    """把 ``--tax-ids`` 参数解析成 ID 列表。"""
    ids: list[str] = []
    for item in items:
        try:
            is_file = Path(item).is_file()
        except OSError:  # pragma: no cover
            is_file = False
        if is_file:
            LOGGER.info("从文件读取 Taxonomy ID: %s", item)
            ids.extend(read_id_file(Path(item)))
        else:
            stripped = item.strip()
            if stripped and stripped not in ids:
                ids.append(stripped)
    return ids


def build_rows(records: Sequence[dict]) -> list[dict]:
    """把解析出的 taxon 记录展开成扁平的输出行。

    ``by_rank`` 里保存的是 ``{rank: name}``；这里把常用层级摊平成独立列，
    同时保留一条完整的 lineage 字符串。
    """
    rows: list[dict] = []
    for record in records:
        by_rank = record.get("by_rank") or {}
        lineage = record.get("lineage") or []

        row: dict = {
            "tax_id": record.get("tax_id", ""),
            "scientific_name": record.get("scientific_name", ""),
            "common_name": record.get("common_name", ""),
            "rank": record.get("rank", ""),
            "lineage_depth": record.get("lineage_depth", len(lineage) + 1),
            "full_lineage": " > ".join(
                item.get("name", "") for item in lineage if item.get("name")
            ),
        }
        for rank_name in FLATTENED_RANKS:
            row[rank_name] = by_rank.get(rank_name, "")
        # 自身不在 LineageEx 里（那里面只有祖先），所以单独补一列
        scientific_name = record.get("scientific_name", "")
        if scientific_name and scientific_name not in row["full_lineage"]:
            row["full_lineage"] = (
                f"{row['full_lineage']} > {scientific_name}" if row["full_lineage"] else scientific_name
            )
        rows.append(row)
    return rows


def summarize_rank_coverage(rows: Sequence[dict]) -> dict[str, int]:
    """统计各 rank 的覆盖率，用于提示哪些层级缺失。"""
    coverage: dict[str, int] = {}
    for row in rows:
        for rank_name in ELEMENTARY_RANKS:
            if row.get(rank_name):
                coverage[rank_name] = coverage.get(rank_name, 0) + 1
    return coverage


def main(argv: Optional[list[str]] = None) -> int:
    """命令行主入口。"""
    parser = build_parser()
    args = parser.parse_args(argv)

    if not args.email:
        LOGGER.error(
            "必须提供 email —— NCBI 要求调用方标识身份。"
            "用 --email 传入，或设置环境变量 BIOPIPE_EMAIL"
        )
        return 2

    tax_ids = collect_tax_ids(args.tax_ids)
    if not tax_ids:
        LOGGER.error("没有可用的 Taxonomy ID")
        return 2

    LOGGER.info("共 %d 个 Taxonomy ID，开始查询 NCBI", len(tax_ids))

    try:
        records = fetch_taxonomy(
            tax_ids,
            email=args.email,
            api_key=args.api_key or None,
            timeout=args.timeout,
        )
    except ValueError as exc:
        LOGGER.error("%s", exc)
        return 2

    if not records:
        LOGGER.error("没有拿到任何结果")
        return 1

    rows = build_rows(records)

    output_dir = ensure_dir(Path(args.output_dir))

    csv_path = output_dir / "taxonomy.csv"
    with open(csv_path, "w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    LOGGER.info("结果已写入: %s", csv_path)

    json_path = output_dir / "taxonomy.json"
    json_path.write_text(
        json.dumps(records, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    LOGGER.info("结果已写入: %s", json_path)

    coverage = summarize_rank_coverage(rows)
    LOGGER.info("共 %d 条记录，输出目录: %s", len(rows), output_dir.resolve())
    LOGGER.info("层级覆盖: %s", ", ".join(f"{k}={v}" for k, v in coverage.items()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
