#!/usr/bin/env python3
"""UniProt 元数据查询客户端（基于共享层 ``http_session``）。

仅在显式调用时才会发起网络请求；``requests`` 缺失时优雅降级为 ``None``。

接口：

- :func:`fetch_uniprot_entry` —— 按 accession 拉取单条 UniProtKB 记录
- :func:`search_uniprot` —— 批量搜索（``accession:X OR accession:Y``）
- :func:`parse_uniprot_entry` —— 把 UniProt JSON 扁平化为表格行

Usage:
    from uniprot_client import fetch_uniprot_entries
    rows = fetch_uniprot_entries(["P69905"], session=None)
"""

from __future__ import annotations

import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "packages"))

from biopipe_core import get_logger, http_session  # noqa: E402

__all__ = ["fetch_uniprot_entries", "parse_uniprot_entry", "UNIPROT_BASE_URL"]

UNIPROT_BASE_URL = "https://rest.uniprot.org"

LOGGER = get_logger("biopipe.protein-structure")

# 需要从 UniProt 请求的字段
QUERY_FIELDS = (
    "accession",
    "id",
    "protein_name",
    "gene_names",
    "organism_name",
    "organism_id",
    "length",
    "mass",
    "entry_type",
    "reviewed",
)


def parse_uniprot_entry(entry: dict) -> dict:
    """把 UniProtKB JSON 记录扁平化为一行表格数据。

    提取字段：accession、entry_name、entry_type、protein_name、gene_names、
    organism_name、organism_id、length、mass、function、keywords、go_terms。
    """
    record: dict[str, Any] = {
        "accession": entry.get("primaryAccession", ""),
        "entry_name": entry.get("uniProtkbId", ""),
        "entry_type": entry.get("entryType", ""),
        "protein_name": "",
        "gene_names": "",
        "organism_name": "",
        "organism_id": "",
        "length": "",
        "mass": "",
        "function": "",
        "keywords": "",
        "go_terms": "",
    }

    description = entry.get("proteinDescription") or {}
    recommended = description.get("recommendedName") or {}
    full_name = recommended.get("fullName") or {}
    record["protein_name"] = full_name.get("value", "")

    genes = entry.get("genes") or []
    names = []
    for gene in genes:
        value = (gene.get("geneName") or {}).get("value")
        if value:
            names.append(value)
        for synonym in gene.get("synonyms") or []:
            svalue = synonym.get("value")
            if svalue:
                names.append(svalue)
    record["gene_names"] = ", ".join(names)

    organism = entry.get("organism") or {}
    record["organism_name"] = organism.get("scientificName", "")
    record["organism_id"] = organism.get("taxonId", "")

    sequence = entry.get("sequence") or {}
    record["length"] = sequence.get("length", "")
    record["mass"] = sequence.get("molWeight", "")

    for comment in entry.get("comments") or []:
        if comment.get("commentType") == "FUNCTION":
            texts = comment.get("texts") or []
            if texts:
                record["function"] = texts[0].get("value", "")
            break

    keywords = entry.get("keywords") or []
    record["keywords"] = "; ".join(k.get("name", "") for k in keywords if k.get("name"))

    go_terms = entry.get("uniProtKBCrossReferences") or []
    go_list = []
    for xref in go_terms:
        if xref.get("database") == "GO":
            props = {p.get("key"): p.get("value") for p in xref.get("properties") or []}
            term = props.get("GoTerm", "")
            aspect = props.get("GoAspect", "")
            if term:
                go_list.append(f"{term} ({aspect})")
    record["go_terms"] = "; ".join(go_list)

    return record


def fetch_uniprot_entries(
    accessions: Sequence[str],
    session=None,
    batch_size: int = 100,
) -> list[dict]:
    """按 accession 批量拉取 UniProtKB 元数据。

    使用 UniProt 搜索接口，一次最多提交 ``batch_size`` 个 accession，
    查询串形如 ``accession:P69905 OR accession:P02340``。

    Args:
        accessions: UniProt accession 列表。
        session: 可复用的 ``requests.Session``，None 时用共享层创建带重试的会话。
        batch_size: 单次请求的 accession 数量上限。

    Returns:
        扁平化后的记录列表（每个字典一行）。网络不可用时返回空列表。
    """
    if not accessions:
        return []

    own_session = session is None
    if own_session:
        session = http_session(min_interval=0.34, retries=4)
    if session is None:
        LOGGER.warning("requests 不可用，跳过 UniProt 查询")
        return []

    results: list[dict] = []
    try:
        for start in range(0, len(accessions), batch_size):
            chunk = [a for a in accessions[start : start + batch_size] if a]
            if not chunk:
                continue
            query = " OR ".join(f"accession:{acc}" for acc in chunk)
            url = f"{UNIPROT_BASE_URL}/uniprotkb/search"
            params = {
                "query": query,
                "format": "json",
                "fields": ",".join(QUERY_FIELDS),
                "size": str(len(chunk)),
            }
            LOGGER.info("查询 UniProt: %d 个 accession", len(chunk))
            try:
                response = session.get(url, params=params, timeout=60)
                response.raise_for_status()
                payload = response.json()
            except Exception as exc:  # noqa: BLE001
                LOGGER.warning("UniProt 查询失败（%d 个 accession）: %s", len(chunk), exc)
                continue

            for entry in payload.get("results", []):
                results.append(parse_uniprot_entry(entry))
    finally:
        if own_session:
            session.close()

    LOGGER.info("UniProt 返回 %d 条记录", len(results))
    return results


def _main(argv: list[str] | None = None) -> int:
    """命令行入口：查询 UniProt 并打印结果。"""
    import argparse
    import json

    parser = argparse.ArgumentParser(
        prog="uniprot_client", description="UniProt 元数据查询（命令行自检工具）"
    )
    parser.add_argument("accessions", nargs="+", help="UniProt accession（可多个）")
    args = parser.parse_args(argv)

    rows = fetch_uniprot_entries(args.accessions)
    if not rows:
        LOGGER.warning("没有查到记录")
        return 1
    for row in rows:
        LOGGER.info(
            "%s | %s | %s | %s aa | %s",
            row.get("accession", ""),
            row.get("protein_name", ""),
            row.get("organism_name", ""),
            row.get("length", ""),
            row.get("gene_names", ""),
        )
    LOGGER.debug("完整记录: %s", json.dumps(rows, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(_main())
