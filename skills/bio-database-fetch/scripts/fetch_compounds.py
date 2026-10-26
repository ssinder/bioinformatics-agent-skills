#!/usr/bin/env python3
"""从 PubChem 获取化合物属性、同义词与 CAS 号。

对应老仓库 ``crawl_database/scripts/04_pubchem_compound_info.py`` 与
``05_pubchem_cas_retriever.py``。老脚本用 ``pubchempy`` 逐个查 CID，
输出是一行行的文本；这里直接调 PUG REST，**一次请求可取 100 个 CID**，
输出结构化的 CSV/JSON。

相对老脚本的改进：

1. **批量查询**。老脚本每个 CID 一次请求，1000 个化合物就是 1000 次
   请求；这里 100 个一批，降到 1/100。
2. **CAS 号取对了地方**。老脚本从 ``xrefs/RegistryID`` 取"CAS"，但那个
   端点返回的是几百个不同数据库的注册号，真正的 CAS 只是其中符合数字
   格式的少数几个。这里从 ``synonyms`` 里按 CAS 的格式筛选。
3. **能发现"没找到"**。PubChem 对不存在的 CID 返回 HTTP 200，内容是
   ``{"CID": 999999999}``——只有 CID 没有属性。老脚本靠 ``if compounds:``
   判断，而 ``pubchempy`` 对不存在的 CID 可能返回空列表也可能抛异常。
   这里显式给出 ``found`` 标记。
4. **属性名兼容**。PubChem 把 ``CanonicalSMILES`` 改名成了 ``SMILES``，
   老脚本按老名字取会拿到空值。

Examples:
    # 从文件读取 CID
    python fetch_compounds.py --cids cids.txt --output-dir results/

    # 直接给 CID
    python fetch_compounds.py --cids 2244 3672 --output-dir results/

    # 只要基本属性，不要几百个同义词
    python fetch_compounds.py --cids cids.txt --output-dir results/ --no-synonyms
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
from database_client import (  # noqa: E402
    COMPOUND_PROPERTIES,
    extract_cas_numbers,
    fetch_compounds,
    fetch_registry_ids,
    fetch_synonyms,
)

LOGGER = get_logger("biopipe.bio-database-fetch")

#: 输出 CSV 的列顺序
FIELDNAMES = [
    "CID",
    "found",
    "IUPACName",
    "MolecularFormula",
    "MolecularWeight",
    "ExactMass",
    "XLogP",
    "TPSA",
    "Charge",
    "SMILES",
    "ConnectivitySMILES",
    "InChI",
    "InChIKey",
    "HBondDonorCount",
    "HBondAcceptorCount",
    "RotatableBondCount",
    "HeavyAtomCount",
    "CAS",
    "synonym_count",
]


def build_parser() -> argparse.ArgumentParser:
    """构造命令行参数解析器。"""
    parser = argparse.ArgumentParser(
        prog="fetch_compounds",
        description="从 PubChem PUG REST 批量获取化合物属性、同义词与 CAS 号",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "示例:\n"
            "  python fetch_compounds.py --cids cids.txt --output-dir results/\n"
            "  python fetch_compounds.py --cids 2244 3672 --output-dir results/\n"
            "\n"
            "说明:\n"
            "  · 一次请求可取 100 个 CID，比逐个查询快两个数量级。\n"
            "  · CAS 号从同义词里按 CAS 的数字格式筛选（不是从 RegistryID）。\n"
            "  · 不存在的 CID 会带上 found=False 标记，不会静默丢弃。\n"
        ),
    )
    parser.add_argument(
        "--cids",
        nargs="+",
        required=True,
        help="CID 列表，或包含 CID 的文本文件路径（每行一个）",
    )
    parser.add_argument(
        "--output-dir",
        default="compound_results",
        help="结果输出目录（默认: compound_results）",
    )
    parser.add_argument(
        "--no-synonyms",
        action="store_true",
        help="不同义词（每个化合物可能有几百个同义词，输出会很大）",
    )
    parser.add_argument(
        "--no-cas",
        action="store_true",
        help="不提取 CAS 号",
    )
    parser.add_argument(
        "--max-synonyms",
        type=int,
        default=20,
        help="每条记录最多保留多少个同义词（默认: 20）",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=60,
        help="单次请求超时秒数（默认: 60）",
    )
    return parser


def read_cid_file(path: Path) -> list[str]:
    """从文本文件读取 CID 列表。

    自动跳过表头行（如果首行不是纯数字）、注释行（``#``）与空行。
    """
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    cids: list[str] = []
    for index, raw in enumerate(lines):
        stripped = raw.strip()
        if not stripped or stripped.startswith("#"):
            continue
        candidate = stripped.split()[0].split(",")[0]
        if not candidate.isdigit():
            # 首行可能是表头；其它位置的非法行只跳过不报错
            if index == 0:
                LOGGER.info("跳过表头行: %s", stripped)
            continue
        if candidate not in cids:
            cids.append(candidate)
    return cids


def collect_cids(items: Sequence[str]) -> list[str]:
    """把 ``--cids`` 参数解析成 CID 列表。"""
    cids: list[str] = []
    for item in items:
        try:
            is_file = Path(item).is_file()
        except OSError:  # pragma: no cover
            is_file = False
        if is_file:
            LOGGER.info("从文件读取 CID: %s", item)
            cids.extend(read_cid_file(Path(item)))
        else:
            stripped = item.strip()
            if stripped and stripped not in cids:
                cids.append(stripped)
    return cids


def build_rows(
    records: Sequence[dict],
    synonyms: dict[str, list[str]],
    registry: dict[str, list[str]],
    max_synonyms: int,
) -> list[dict]:
    """把各端点的结果合成输出行。"""
    rows: list[dict] = []
    for record in records:
        cid = str(record.get("CID", ""))
        synonym_list = synonyms.get(cid, [])
        row = dict(record)
        row["CID"] = cid

        # CAS 号：优先从同义词里找（同义词更全），再补 RegistryID 里的
        cas_numbers = extract_cas_numbers(synonym_list)
        for extra in extract_cas_numbers(registry.get(cid, [])):
            if extra not in cas_numbers:
                cas_numbers.append(extra)
        row["CAS"] = "; ".join(cas_numbers)
        row["synonym_count"] = len(synonym_list)

        if max_synonyms > 0:
            row["synonyms"] = "; ".join(synonym_list[:max_synonyms])

        rows.append(row)
    return rows


def main(argv: Optional[list[str]] = None) -> int:
    """命令行主入口。"""
    args = build_parser().parse_args(argv)

    cids = collect_cids(args.cids)
    if not cids:
        LOGGER.error("没有可用的 CID")
        return 2

    LOGGER.info("共 %d 个 CID，开始查询 PubChem", len(cids))

    records = fetch_compounds(cids, properties=COMPOUND_PROPERTIES, timeout=args.timeout)
    if not records:
        LOGGER.error("没有拿到任何结果")
        return 1

    synonyms: dict[str, list[str]] = {}
    registry: dict[str, list[str]] = {}
    if not args.no_synonyms:
        synonyms = fetch_synonyms(cids, timeout=args.timeout)
    if not args.no_cas:
        registry = fetch_registry_ids(cids, timeout=args.timeout)

    rows = build_rows(records, synonyms, registry, args.max_synonyms)

    output_dir = ensure_dir(Path(args.output_dir))

    csv_path = output_dir / "compounds.csv"
    fieldnames = list(FIELDNAMES)
    if rows and "synonyms" in rows[0]:
        fieldnames.append("synonyms")
    with open(csv_path, "w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    LOGGER.info("结果已写入: %s", csv_path)

    json_path = output_dir / "compounds.json"
    json_path.write_text(json.dumps(rows, ensure_ascii=False, indent=2), encoding="utf-8")
    LOGGER.info("结果已写入: %s", json_path)

    found = sum(1 for row in rows if row.get("found"))
    with_cas = sum(1 for row in rows if row.get("CAS"))
    LOGGER.info("共 %d 条记录（找到 %d，未找到 %d，含 CAS %d），输出目录: %s",
                len(rows), found, len(rows) - found, with_cas, output_dir.resolve())
    return 0


if __name__ == "__main__":
    sys.exit(main())
