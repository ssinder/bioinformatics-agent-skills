#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""多源学术文献检索与筛选入口。

功能链路：多源检索（PubMed / OpenAlex / Semantic Scholar）-> 归一化到统一 schema
-> 去重（DOI 优先 + 标题模糊匹配）-> 相关性打分与筛选 -> 摘要富集 -> 导出 CSV/JSON。

所有外部 API 调用都通过共享层 ``http_session`` 完成（内置指数退避重试与限速）。
主题、关键词、年份范围、文献源全部参数化，不写死任何具体疾病主题。

用法示例::

    # PubMed 检索 2020 年以来的抗菌肽文献
    python search_literature.py --topic "antimicrobial peptide" --source pubmed \\
        --date-range 2020: --output results/amp.csv

    # 三源联合检索 + 摘要富集 + 相关性过滤
    python search_literature.py --topic "protein structure prediction" \\
        --keywords "AlphaFold" "molecular dynamics" \\
        --source pubmed openalex semantic-scholar \\
        --date-range 2018:2024 --min-citations 20 --enrich-abstracts \\
        --output results/prediction.csv

    # 只处理已有 CSV（不联网），重新打分并去重
    python search_literature.py --input results/raw.csv --keywords "antifungal" \\
        --min-relevance 3 --output results/filtered.csv
"""

from __future__ import annotations

import argparse
import os
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
    write_csv,
)

from literature_core import (  # noqa: E402
    CSV_FIELDS,
    SOURCE_LABELS,
    LiteratureClient,
    build_query,
    deduplicate,
    enrich_abstracts,
    filter_records,
    parse_date_range,
    parse_record_rows,
    score_relevance,
)

LOGGER = get_logger("biopipe.literature-search")


def build_parser() -> argparse.ArgumentParser:
    """构建命令行解析器（在共享层公共参数之上追加业务参数）。"""
    parser = cli_parser(
        prog="literature-search",
        description="多源学术文献检索、去重、筛选与摘要富集",
    )
    parser.add_argument("--topic", help="自由文本检索主题，如 'antimicrobial peptide'")
    parser.add_argument(
        "--keywords",
        nargs="+",
        default=[],
        help="附加关键词（作为 AND 条件收窄检索，同时用于相关性打分）",
    )
    parser.add_argument(
        "--source",
        nargs="+",
        default=["pubmed"],
        choices=sorted(SOURCE_LABELS),
        help="检索的数据源，可多选",
    )
    parser.add_argument("--date-range", help="年份范围，形如 2018:2024 / 2020: / :2020")
    parser.add_argument("--max-results", type=int, default=50, help="每个数据源每个查询的最大返回条数")
    parser.add_argument(
        "--input",
        help="已有 CSV 输入（提供时不联网，直接对本地数据去重/打分/筛选）",
    )
    parser.add_argument("--output", help="输出 CSV 路径")
    parser.add_argument("--json-output", help="额外输出完整 JSON 记录路径")
    parser.add_argument("--min-citations", type=int, default=0, help="最低引用数")
    parser.add_argument(
        "--min-relevance",
        type=float,
        default=0.0,
        help="最低相关性得分（见 references/methodology.md 评分公式）",
    )
    parser.add_argument("--min-abstract-length", type=int, default=0, help="最短摘要长度")
    parser.add_argument(
        "--include-terms",
        nargs="+",
        default=[],
        help="必须包含的词（标题/摘要/期刊/关键词，命中任一即保留）",
    )
    parser.add_argument(
        "--exclude-terms",
        nargs="+",
        default=[],
        help="必须排除的词（命中任一即剔除）",
    )
    parser.add_argument("--title-similarity", type=float, default=0.92, help="标题模糊去重阈值")
    parser.add_argument("--enrich-abstracts", action="store_true", help="为缺失摘要的文献补全摘要")
    parser.add_argument("--min-abstract-for-enrich", type=int, default=60, help="视为需要富化的摘要长度阈值")
    parser.add_argument("--email", help="联系邮箱（OpenAlex 礼貌池 + NCBI 建议提供）")
    parser.add_argument("--ncbi-api-key", help="NCBI E-utilities API Key（环境变量 NCBI_API_KEY 优先）")
    parser.add_argument("--s2-api-key", help="Semantic Scholar API Key（环境变量 S2_API_KEY 优先）")
    parser.add_argument("--min-interval", type=float, default=0.34, help="请求最小间隔秒数")
    parser.add_argument(
        "--sort-by",
        default="relevance",
        choices=["relevance", "citations", "year"],
        help="输出排序依据",
    )
    parser.add_argument("--limit", type=int, default=0, help="最多输出多少条（0 表示不限）")
    return parser


def _run_search(
    client: LiteratureClient,
    args: argparse.Namespace,
    query: str,
    start_year: int | None,
    end_year: int | None,
) -> list[dict]:
    """按选定数据源执行检索并汇总记录。"""
    records: list[dict] = []
    for source in args.source:
        try:
            if source == "pubmed":
                found = client.search_pubmed(query, args.max_results, start_year, end_year)
            elif source == "openalex":
                found = client.search_openalex(query, args.max_results, start_year, end_year)
            else:
                found = client.search_semantic_scholar(query, args.max_results, start_year, end_year)
        except Exception as exc:  # noqa: BLE001 - 单源失败不应中断整体检索
            LOGGER.warning("数据源 %s 检索失败: %s", source, exc)
            continue
        records.extend(found)
    return records


def main(args: argparse.Namespace) -> int:
    """主流程。"""
    config = get_config(args.config)
    logger = LOGGER

    start_year, end_year = parse_date_range(args.date_range)
    keywords = [k.strip() for k in (args.keywords or []) if k.strip()]

    # ---- 1. 载入或检索 -------------------------------------------------
    if args.input:
        input_path = resolve_path(args.input)
        if not input_path.exists():
            raise SkillError(f"输入文件不存在: {input_path}")
        logger.info("读取本地文件: %s", input_path)
        records = parse_record_rows(read_csv_safe(input_path))
        logger.info("本地载入 %d 条记录", len(records))
    else:
        query = build_query(args.topic, keywords)
        if not query:
            raise SkillError("必须提供 --topic、--keywords 或 --input 之一")
        logger.info("检索式: %s", query)
        client = LiteratureClient(
            email=args.email,
            ncbi_api_key=args.ncbi_api_key or os.environ.get("NCBI_API_KEY"),
            s2_api_key=args.s2_api_key or os.environ.get("S2_API_KEY"),
            min_interval=args.min_interval,
        )
        try:
            records = _run_search(client, args, query, start_year, end_year)
            if args.enrich_abstracts:
                enrich_abstracts(
                    records,
                    client,
                    min_length=args.min_abstract_for_enrich,
                )
        finally:
            client.close()

    if not records:
        logger.warning("没有检索到任何记录")
        return 0

    # ---- 2. 打分 -------------------------------------------------------
    for record in records:
        record["relevance"] = score_relevance(record, keywords)

    # ---- 3. 去重 -------------------------------------------------------
    deduped, removed = deduplicate(records, title_threshold=args.title_similarity)
    logger.info("去重完成: %d -> %d 条（移除 %d）", len(records), len(deduped), removed)

    # ---- 4. 筛选 -------------------------------------------------------
    selected = filter_records(
        deduped,
        min_relevance=args.min_relevance,
        min_citations=args.min_citations,
        min_abstract_length=args.min_abstract_length,
        include_terms=args.include_terms,
        exclude_terms=args.exclude_terms,
        start_year=start_year,
        end_year=end_year,
    )
    logger.info("筛选后保留 %d 条", len(selected))

    # ---- 5. 排序与截断 -------------------------------------------------
    if args.sort_by == "citations":
        selected.sort(key=lambda r: int(r.get("citations", 0)), reverse=True)
    elif args.sort_by == "year":
        selected.sort(key=lambda r: (r.get("year") or 0), reverse=True)
    else:
        selected.sort(
            key=lambda r: (float(r.get("relevance", 0)), int(r.get("citations", 0))),
            reverse=True,
        )
    if args.limit > 0:
        selected = selected[: args.limit]

    # ---- 6. 输出 -------------------------------------------------------
    if not selected:
        logger.warning("筛选后没有记录，未写出结果文件")
        return 0

    output_path = (
        resolve_path(args.output)
        if args.output
        else config.output_dir / "literature_results.csv"
    )
    rows = []
    for record in selected:
        row = {field: record.get(field, "") for field in CSV_FIELDS}
        row["year"] = record.get("year") or ""
        row["citations"] = record.get("citations", 0)
        row["relevance"] = record.get("relevance", 0.0)
        rows.append(row)
    write_csv(rows, output_path, fieldnames=CSV_FIELDS)
    logger.info("已写出 %d 条记录: %s", len(rows), output_path)

    if args.json_output:
        import json

        json_path = resolve_path(args.json_output)
        json_path.parent.mkdir(parents=True, exist_ok=True)
        with open(json_path, "w", encoding="utf-8") as handle:
            json.dump(selected, handle, ensure_ascii=False, indent=2)
        logger.info("已写出 JSON: %s", json_path)

    return 0


if __name__ == "__main__":
    sys.exit(
        run_cli(
            main,
            parser=build_parser(),
            skill="literature-search",
        )
    )
