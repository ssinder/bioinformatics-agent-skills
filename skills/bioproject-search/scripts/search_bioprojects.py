#!/usr/bin/env python3
"""批量检索 NCBI BioProject 元数据。

核心改进（相对老仓库 ``crawl_bioproject/``）：

1. **用 NCBI E-utilities 取代 HTML 爬虫**。老脚本用 requests +
   BeautifulSoup 抓网页，靠 ``nth-child(5)`` 绝对 XPath 定位描述，
   页面改版即失效。E-utilities 返回结构化 XML。
2. **不再写死 Windows 桌面路径**。老脚本的 ``base_dir`` 指向
   ``C:\\Users\\DELL\\Desktop\\trae\\bioproject``，换台机器首行日志就
   FileNotFoundError。
3. **不再无条件跳过前三行**。老脚本用 ``all_lines[3:]`` 假设文件
   前三行是示例，配置与数据混在一起。
4. **不伪造 NCBI email**。老脚本随机生成 ``userN@example.com``——Entrez
   的 email 是身份标识，伪造既违规也会被 NCBI 限速封禁。
5. **删掉死代码**。老脚本的 ``check_has_description_link()`` 定义了
   从未被调用，且内部绕过了重试与 UA 逻辑，另有一套请求路径。
6. **不再每条记录重写整个 CSV**（O(n²) 磁盘写入），改为先收集后
   一次性写出。

Examples:
    # 从文件读取 accession
    python search_bioprojects.py --ids prj_ids.txt --output-dir results/

    # 直接给 accession
    python search_bioprojects.py --ids PRJNA172563 PRJNA212358 --output-dir results/

    # 只要描述的关键词筛选
    python search_bioprojects.py --ids prj_ids.txt --output-dir results/ \
        --filter "CRISPR"
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Optional

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "packages"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from biopipe_core import ensure_dir, get_logger  # noqa: E402
from bioproject_client import search_bioprojects  # noqa: E402

LOGGER = get_logger("biopipe.bioproject-search")

# BioProject accession 的典型形态：PRJNA/PRJEB/PRJDB + 数字
_ACCESSION_RE = re.compile(r"^PRJ(?:NA|EB|DB)\d+$", re.IGNORECASE)

# 常见的导航/页脚噪声词。老脚本的 invalid_keywords 列表与真实描述语义
# 不符，容易误杀有效描述——这里只在"描述过短"时才启用这些词做兜底判断。
_NOISE_MARKERS = (
    "navigate across",
    "additional projects",
    "learn more",
    "click here",
)


def read_id_file(path: Path) -> list[str]:
    """从文本文件读取 BioProject accession 列表。

    与老脚本的"无条件跳过前三行"不同，这里按 accession 格式自动判断
    首行是否为表头；注释行（``#``）与空行一律忽略。

    Args:
        path: 文件路径，每行一个 accession。

    Returns:
        去重后的 accession 列表（保持原顺序）。
    """
    if not path.is_file():
        raise FileNotFoundError(f"accession 列表文件不存在: {path}")

    lines: list[str] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        stripped = raw.strip()
        if not stripped or stripped.startswith("#"):
            continue
        # 兼容 "PRJNA123  some note" 这种带注释的写法
        candidate = stripped.split()[0]
        if candidate and candidate.lower() not in {x.lower() for x in lines}:
            lines.append(candidate)
    return lines


def is_valid_accession(text: str) -> bool:
    """判断字符串是否像 BioProject accession。"""
    return bool(_ACCESSION_RE.match(text.strip()))


def _looks_like_noise(description: str) -> bool:
    """判断描述是否只是页面导航噪声而非真实描述。

    只有描述**很短**且含导航词时才判定为噪声——老脚本对任意长度
    的描述都套用噪声词列表，会误杀正文里恰好提到这些词的合法描述。
    """
    if len(description) > 200:
        return False
    lowered = description.lower()
    return any(marker in lowered for marker in _NOISE_MARKERS)


def filter_records(records: Sequence[dict], keyword: Optional[str]) -> list[dict]:
    """按关键词过滤记录（不区分大小写，匹配标题或描述）。

    Args:
        records: 记录列表。
        keyword: 关键词；None 或空串表示不过滤。

    Returns:
        保留的记录列表。
    """
    if not keyword:
        return list(records)

    needle = keyword.lower()
    kept: list[dict] = []
    for record in records:
        haystack = f"{record.get('title', '')} {record.get('description', '')} {record.get('name', '')}"
        if needle in haystack.lower():
            kept.append(record)
    return kept


def write_outputs(records: Sequence[dict], output_dir: Path) -> None:
    """把记录写成 CSV 与 JSON。"""
    output_dir.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "accession",
        "title",
        "name",
        "description",
        "organism",
        "taxid",
        "study_type",
        "submitter",
        "registration_date",
        "release_date",
        "last_update_date",
        "uid",
    ]

    csv_path = output_dir / "bioproject_data.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for record in records:
            writer.writerow(record)

    json_path = output_dir / "bioproject_data.json"
    json_path.write_text(
        json.dumps(list(records), ensure_ascii=False, indent=2), encoding="utf-8"
    )

    LOGGER.info("结果已写入: %s / %s", csv_path, json_path)


def build_parser() -> argparse.ArgumentParser:
    """构造命令行参数解析器。"""
    parser = argparse.ArgumentParser(
        prog="search_bioprojects",
        description="批量检索 NCBI BioProject 元数据（标题 / 描述 / 物种等）",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "示例:\n"
            "  python search_bioprojects.py --ids prj_ids.txt --output-dir results/ --email you@example.org\n"
            "  python search_bioprojects.py --ids PRJNA172563 --output-dir results/ --email you@example.org\n"
            "\n"
            "说明:\n"
            "  数据来自 NCBI E-utilities（结构化 XML），不解析网页，因此不受网页改版影响。\n"
            "  NCBI 要求提供联系邮箱；可用 --email 传入或设置环境变量 BIOPIPE_EMAIL。\n"
        ),
    )
    parser.add_argument(
        "--ids",
        nargs="+",
        required=True,
        help="BioProject accession 列表，或包含 accession 的文本文件路径",
    )
    parser.add_argument(
        "--output-dir",
        default="bioproject_results",
        help="结果输出目录（默认: bioproject_results）",
    )
    parser.add_argument(
        "--email",
        default=os.environ.get("BIOPIPE_EMAIL", ""),
        help="联系邮箱（NCBI 要求；也可用环境变量 BIOPIPE_EMAIL）",
    )
    parser.add_argument(
        "--api-key",
        default=os.environ.get("BIOPIPE_API_KEY", ""),
        help="NCBI API key（可选；也可用环境变量 BIOPIPE_API_KEY）",
    )
    parser.add_argument(
        "--filter",
        dest="keyword",
        default=None,
        help="只保留标题/描述中包含该关键词的记录",
    )
    parser.add_argument(
        "--drop-noisy",
        action="store_true",
        help="丢弃疑似页面导航噪声的极短描述",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=60,
        help="单次请求超时秒数（默认: 60）",
    )
    return parser


def _collect_accessions(items: Sequence[str]) -> list[str]:
    """把 --ids 参数解析成 accession 列表。"""
    accessions: list[str] = []
    for item in items:
        try:
            is_file = Path(item).is_file()
        except OSError:  # pragma: no cover
            is_file = False
        if is_file:
            LOGGER.info("从文件读取 accession: %s", item)
            accessions.extend(read_id_file(Path(item)))
        else:
            accessions.append(item.strip())
    return accessions


def main(argv: Optional[list[str]] = None) -> int:
    """命令行主入口。"""
    args = build_parser().parse_args(argv)

    if not args.email:
        LOGGER.error(
            "必须提供 --email 或设置环境变量 BIOPIPE_EMAIL —— "
            "NCBI 要求调用方标识身份（老脚本伪造 userN@example.com 会被限速封禁）"
        )
        return 2

    accessions = _collect_accessions(args.ids)
    if not accessions:
        LOGGER.error("没有可用的 BioProject accession")
        return 2

    # 过滤掉明显不是 accession 的行（但保留它们用于提示）
    valid = [a for a in accessions if is_valid_accession(a)]
    invalid = [a for a in accessions if a not in valid]
    if invalid:
        LOGGER.warning("以下条目不像 BioProject accession，已跳过: %s", ", ".join(invalid[:10]))
    if not valid:
        LOGGER.error("没有有效的 BioProject accession")
        return 2

    LOGGER.info("开始查询 %d 个 BioProject", len(valid))
    records = search_bioprojects(
        valid,
        email=args.email,
        api_key=args.api_key or None,
        timeout=args.timeout,
    )

    if not records:
        LOGGER.error("没有查到任何 BioProject 记录")
        return 1

    if args.drop_noisy:
        before = len(records)
        records = [r for r in records if not _looks_like_noise(r.get("description", ""))]
        if before != len(records):
            LOGGER.info("丢弃 %d 条疑似导航噪声记录", before - len(records))

    records = filter_records(records, args.keyword)
    if args.keyword:
        LOGGER.info("按关键词 '%s' 过滤后剩 %d 条", args.keyword, len(records))

    if not records:
        LOGGER.error("过滤后没有剩余记录")
        return 1

    output_dir = ensure_dir(args.output_dir)
    write_outputs(records, output_dir)

    LOGGER.info("共 %d 条记录，输出目录: %s", len(records), output_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
