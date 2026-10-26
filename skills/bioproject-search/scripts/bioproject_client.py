#!/usr/bin/env python3
"""NCBI BioProject 检索客户端（基于 NCBI E-utilities）。

老仓库的 ``crawl_bioproject/`` 用 requests + BeautifulSoup 爬 NCBI
BioProject 网页取描述，依赖脆弱的 ``nth-child(5)`` 绝对 XPath。本模块改用
NCBI 官方的 E-utilities 接口，直接拿到结构化的 XML 记录。

三个端点的分工
--------------
BioProject 的元数据分散在三个 E-utilities 端点里，**没有任何一个端点
能拿到全部字段**。这是本模块设计上最关键的一点——只调一个端点会静默
丢失一半的字段，而且不会报错：

======================  ==========================================
端点                   独有字段
======================  ==========================================
``esearch``             UID（accession → 内部 ID 的映射）
``esummary``           TaxId、Organism_Name、Project_Type、
                       Submitter_Organization、Registration_Date
``efetch``             ProjectReleaseDate、外部链接（dbGaP 等）
======================  ==========================================

所以这里用 **esummary + efetch 双源合并**：efetch 拿不到的字段
（taxid、organism、submitter）由 esummary 补齐，esummary 里没有的
字段（release date）由 efetch 补齐。

为什么不能只用 esummary
------------------------
``esummary`` 的 ``Project_Description`` 有时比 ``efetch`` 短，而且
``Registration_Date`` 的格式是 ``2013/07/16 00:00``（斜杠分隔、
无时区），与 ``efetch`` 的 ISO-8601 格式不一致。统一用 efetch 的
描述文本、esummary 的日期，并在合并时做归一化。

E-utilities 使用要求
--------------------
NCBI 要求所有调用方提供 ``email`` 与 ``api_key``（可选但强烈建议）。
没有 email 会被限速，频繁调用可能收到 NCBI 的封禁。因此本模块的
``email`` 参数是**必填**，并且不提供任何伪造的占位邮箱——老代码里
``userN@example.com`` 那种随机伪造的邮箱既违规又无用。

Usage:
    from bioproject_client import search_bioprojects
    rows = search_bioprojects(["PRJNA172563"], email="you@example.org")
"""

from __future__ import annotations

import re
import sys
import xml.etree.ElementTree as ET
from collections.abc import Sequence
from pathlib import Path
from typing import Any, Optional

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "packages"))

from biopipe_core import get_logger, http_session  # noqa: E402

__all__ = [
    "search_bioprojects",
    "parse_bioproject_xml",
    "parse_esummary_xml",
    "merge_records",
    "build_esearch_term",
    "normalize_date",
    "EUTILS_BASE",
    "BIOPROJECT_DB",
]

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
BIOPROJECT_DB = "bioproject"

LOGGER = get_logger("biopipe.bioproject-search")

# esummary 的 Registration_Date 形如 "2013/07/16 00:00"，还可能是
# "1/01/01 00:00"（占位日期，公共项目里很常见——NCBI 用它表示
# "日期不详"）。这两种都需要归一化。
_DATE_SLASH_RE = re.compile(r"^\s*(\d{4})/(\d{1,2})/(\d{1,2})\s")


def normalize_date(raw: str) -> str:
    """把各种日期写法归一化成 ``YYYY-MM-DD``。

    NCBI 的三个端点用了三套日期格式：

    - ``efetch`` 的 ``ProjectReleaseDate``：``2016-12-15T00:00:00Z``
    - ``esummary`` 的 ``Registration_Date``：``2013/07/16 00:00``
    - ``efetch`` 的 ``Submission/@submitted``：``2013-07-16``

    归一化成 ISO 日期后，下游做时间排序/比较才不会被字符串格式坑到。

    Args:
        raw: 原始日期字符串。

    Returns:
        ``YYYY-MM-DD`` 形式的日期；无法解析时返回空串。
    """
    if not raw:
        return ""

    text = raw.strip()

    # ISO-8601：2016-12-15T00:00:00Z
    m = re.match(r"^(\d{4})-(\d{2})-(\d{2})", text)
    if m:
        return f"{m.group(1)}-{m.group(2)}-{m.group(3)}"

    # 斜杠格式：2013/07/16 00:00
    m = _DATE_SLASH_RE.match(text)
    if m:
        return f"{m.group(1)}-{int(m.group(2)):02d}-{int(m.group(3)):02d}"

    return ""


def build_esearch_term(accessions: Sequence[str]) -> str:
    """构造 E-utilities 的 accession 检索式。

    每个 accession 用 ``OR`` 连接，并加上 ``[Project Accession]`` 字段限定符
    ——不加限定符时，``PRJNA172563`` 这种字符串可能被当作自由词匹配到
    无关记录。

    Args:
        accessions: BioProject accession 列表（PRJNA/PRJEB/PRJDB 等）。

    Returns:
        检索式字符串；accessions 为空时返回空串。
    """
    cleaned = [
        acc.strip()
        for acc in accessions
        if acc and acc.strip()
    ]
    if not cleaned:
        return ""
    return " OR ".join(f"{acc}[Project Accession]" for acc in cleaned)


def _esearch_to_uids(accessions: Sequence[str], session, email: str, api_key: Optional[str]) -> list[str]:
    """把 accession 列表解析成 BioProject 的内部 UID 列表。"""
    term = build_esearch_term(accessions)
    if not term:
        return []

    params: dict[str, Any] = {
        "db": BIOPROJECT_DB,
        "term": term,
        "retmode": "json",
        "retmax": str(len(accessions) * 2 + 10),
    }
    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key

    response = session.get(f"{EUTILS_BASE}/esearch.fcgi", params=params, timeout=60)
    response.raise_for_status()
    payload = response.json()
    result = payload.get("esearchresult", {})
    return [str(uid) for uid in result.get("idlist", [])]


def _parse_organism(doc: ET.Element) -> tuple[str, str]:
    """从 efetch 记录里取物种名与 taxid。

    这里有个**极易踩的坑**：``<Organism>`` 元素的 ``species`` 属性里存的
    是 **taxid**，不是物种名。真实数据长这样::

        <Organism species="28901" taxID="90370">
          <OrganismName>Salmonella enterica subsp. enterica serovar Typhi</OrganismName>

    ``species="28901"`` 是 Salmonella enterica 的物种 taxid，
    ``taxID="90370"`` 是亚种/血清型的 taxid，而真正的名称在
    ``<OrganismName>`` 子元素里。把 ``species`` 当名称用，会得到
    "organism=28901" 这种明显荒谬的结果——而且不会报错，只是值错了。

    因此这里只从子元素 ``OrganismName`` 取名称；taxid 则依次尝试
    ``TaxonomyID`` / ``taxID`` / ``taxid``（不同提交批次用词不一致）。

    Args:
        doc: DocumentSummary 元素。

    Returns:
        ``(organism, taxid)`` 二元组。
    """
    organism = ""
    taxid = ""
    for node in doc.iter("Organism"):
        name_node = node.find("OrganismName")
        if name_node is not None and name_node.text:
            organism = name_node.text.strip()
        elif node.get("OrganismName"):
            organism = node.get("OrganismName").strip()

        tid = (
            node.get("TaxonomyID")
            or node.get("taxID")
            or node.get("taxid")
            or ""
        ).strip()
        if tid:
            taxid = tid
        if organism and taxid:
            break
    return organism, taxid


def _parse_document(doc: ET.Element) -> dict:
    """把一个 efetch 的 DocumentSummary 元素解析成字典。"""
    uid = doc.get("uid", "")

    def _text(path: str) -> str:
        node = doc.find(path)
        if node is None or node.text is None:
            return ""
        return node.text.strip()

    # accession 在 ArchiveID/@accession 属性上
    accession = ""
    archive = doc.find("Project/ProjectID/ArchiveID")
    if archive is not None:
        accession = (archive.get("accession") or "").strip()

    organism, taxid = _parse_organism(doc)

    # 描述可能是多段的 <Description> 文本
    description = ""
    desc_node = doc.find("Project/ProjectDescr/Description")
    if desc_node is not None and desc_node.text:
        description = desc_node.text.strip()
    else:
        # 兜底：拼接所有 Description 片段
        parts = [
            (n.text or "").strip()
            for n in doc.findall("Project/ProjectDescr/Description")
        ]
        description = " ".join(p for p in parts if p)

    # 提交日期在 <Submission submitted="2013-07-16"> 的**属性**上，
    # 不是子元素——只看元素节点会永远拿到空串。
    submission = doc.find("Submission")
    submitted = ""
    last_update = ""
    if submission is not None:
        submitted = normalize_date(submission.get("submitted") or "")
        last_update = normalize_date(submission.get("last_update") or "")

    # 发布日期：ProjectDescr 下的 ProjectReleaseDate（部分记录才有）
    release_date = normalize_date(_text("Project/ProjectDescr/ProjectReleaseDate"))

    # 提交机构在 Submission/Description/Organization/Name
    submitter = _text("Submission/Description/Organization/Name")

    return {
        "uid": uid,
        "accession": accession,
        "title": _text("Project/ProjectDescr/Title"),
        "name": _text("Project/ProjectDescr/Name"),
        "description": description,
        "organism": organism,
        "taxid": taxid,
        "study_type": _text("Project/ProjectType/ProjectTypeSubmission/Method"),
        "submitter": submitter,
        "registration_date": submitted,
        "release_date": release_date,
        "last_update_date": last_update,
    }


def parse_bioproject_xml(xml_text: str) -> list[dict]:
    """把 efetch 返回的 XML 解析成结构化记录列表。

    记录里提取：uid、accession、title、name、description、organism、
    taxid、study_type、submitter、registration_date、release_date。

    Args:
        xml_text: efetch 返回的 XML 文本。

    Returns:
        每个 BioProject 一条记录的列表。
    """
    if not xml_text.strip():
        return []

    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError as exc:
        LOGGER.warning("BioProject XML 解析失败: %s", exc)
        return []

    records: list[dict] = []
    for doc in root.iter("DocumentSummary"):
        record = _parse_document(doc)
        if record:
            records.append(record)
    return records


def parse_esummary_xml(xml_text: str) -> list[dict]:
    """把 esummary 返回的 XML 解析成结构化记录列表。

    esummary 提供 efetch 拿不到的结构化字段：TaxId、Organism_Name、
    Project_Type、Submitter_Organization、Registration_Date。

    Args:
        xml_text: esummary 返回的 XML 文本。

    Returns:
        每个 BioProject 一条记录的列表（字段与 efetch 解析结果同名，
        以便合并）。
    """
    if not xml_text.strip():
        return []

    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError as exc:
        LOGGER.warning("BioProject esummary XML 解析失败: %s", exc)
        return []

    records: list[dict] = []
    for doc in root.iter("DocumentSummary"):
        uid = doc.get("uid", "")

        # 默认参数绑定当前的 doc，避免闭包捕获循环变量
        def _elem_text(name: str, _doc: ET.Element = doc) -> str:
            node = _doc.find(name)
            if node is None or node.text is None:
                return ""
            return node.text.strip()

        submitter = _elem_text("Submitter_Organization")
        if not submitter:
            # Submitter_Organization_List 是 <string>…</string> 的列表
            node = doc.find("Submitter_Organization_List")
            if node is not None:
                parts = [
                    (n.text or "").strip()
                    for n in node.findall("string")
                ]
                submitter = ", ".join(p for p in parts if p)

        records.append(
            {
                "uid": uid,
                "accession": _elem_text("Project_Acc"),
                "title": _elem_text("Project_Title"),
                "name": _elem_text("Project_Name"),
                "description": _elem_text("Project_Description"),
                "organism": _elem_text("Organism_Name"),
                "taxid": _elem_text("TaxId"),
                "study_type": _elem_text("Project_Type"),
                "submitter": submitter,
                "registration_date": normalize_date(_elem_text("Registration_Date")),
                "release_date": "",
                "last_update_date": "",
            }
        )
    return records


def merge_records(efetch_records: Sequence[dict], esummary_records: Sequence[dict]) -> list[dict]:
    """把 esummary 与 efetch 的记录按 UID 合并。

    合并规则：**efetch 优先**，只在 efetch 的字段为空时才用 esummary 的值
    补齐。原因有二：

    1. efetch 的描述更完整（esummary 的 ``Project_Description`` 偶尔会被
       截断）；
    2. efetch 的 ``ProjectReleaseDate`` 是 ISO 格式，而 esummary 的
       ``Registration_Date`` 是斜杠格式，已在解析时各自归一化。

    Args:
        efetch_records: efetch 解析出的记录。
        esummary_records: esummary 解析出的记录。

    Returns:
        合并后的记录列表，以 efetch 的顺序为准。
    """
    summary_by_uid = {rec.get("uid", ""): rec for rec in esummary_records}

    merged: list[dict] = []
    for record in efetch_records:
        uid = record.get("uid", "")
        fallback = summary_by_uid.get(uid, {})
        combined = dict(record)
        for key, value in fallback.items():
            if not combined.get(key):
                combined[key] = value
        merged.append(combined)

    # efetch 里没有、但 esummary 有的记录（理论上不该出现，但要兜底）
    seen = {rec.get("uid", "") for rec in efetch_records}
    for uid, record in summary_by_uid.items():
        if uid and uid not in seen:
            merged.append(dict(record))

    return merged


def search_bioprojects(
    accessions: Sequence[str],
    email: str,
    api_key: Optional[str] = None,
    session=None,
    timeout: int = 60,
) -> list[dict]:
    """按 accession 批量查询 BioProject 记录。

    流程：``esearch`` 把 accession 转成内部 UID，然后**同时**用
    ``esummary`` 与 ``efetch`` 拉详情并合并——两个端点的字段互补，
    只用一个会静默丢失另一半元数据（详见模块文档）。

    Args:
        accessions: BioProject accession 列表。
        email: 联系邮箱——**必填**，NCBI 要求调用方标识身份。
        api_key: NCBI API key（可选，可把限速从 3 req/s 提到 10 req/s）。
        session: 可复用的 ``requests.Session``。
        timeout: 单次请求超时秒数。

    Returns:
        每个 BioProject 一条记录的列表。
    """
    if not accessions:
        return []
    if not email:
        raise ValueError(
            "必须提供 email —— NCBI 要求调用方标识身份，"
            "可用 --email 传入或设置环境变量 BIOPIPE_EMAIL"
        )

    own_session = session is None
    if own_session:
        # NCBI 无 key 时限速 3 req/s，有 key 时 10 req/s
        session = http_session(min_interval=0.34, retries=4)
    if session is None:
        LOGGER.warning("requests 不可用，无法查询 NCBI")
        return []

    records: list[dict] = []
    try:
        uids = _esearch_to_uids(accessions, session, email, api_key)
        if not uids:
            LOGGER.warning("没有匹配到 BioProject 记录")
            return []

        LOGGER.info("匹配到 %d 个 BioProject，拉取详情", len(uids))

        common: dict[str, Any] = {"db": BIOPROJECT_DB}
        if email:
            common["email"] = email
        if api_key:
            common["api_key"] = api_key

        # esummary：结构化字段（taxid / organism / submitter / 注册日期）
        summary_params = {**common, "id": ",".join(uids), "retmode": "xml"}
        summary_resp = session.get(
            f"{EUTILS_BASE}/esummary.fcgi", params=summary_params, timeout=timeout
        )
        summary_resp.raise_for_status()
        esummary_records = parse_esummary_xml(summary_resp.text)

        # efetch：完整描述、发布日期、外部链接
        fetch_params = {**common, "id": ",".join(uids), "rettype": "xml", "retmode": "xml"}
        fetch_resp = session.get(
            f"{EUTILS_BASE}/efetch.fcgi", params=fetch_params, timeout=timeout
        )
        fetch_resp.raise_for_status()
        efetch_records = parse_bioproject_xml(fetch_resp.text)

        records = merge_records(efetch_records, esummary_records)
    except ValueError:
        raise
    except Exception as exc:  # noqa: BLE001
        LOGGER.warning("BioProject 查询失败: %s", exc)
    finally:
        if own_session:
            session.close()

    LOGGER.info("返回 %d 条 BioProject 记录", len(records))
    return records


def _main(argv: Optional[list[str]] = None) -> int:
    """命令行入口：查询 BioProject（自检工具）。"""
    import argparse
    import json
    import os

    parser = argparse.ArgumentParser(
        prog="bioproject_client", description="查询 NCBI BioProject 记录（命令行自检工具）"
    )
    parser.add_argument("accessions", nargs="+", help="BioProject accession（PRJNA…）")
    parser.add_argument(
        "--email",
        default=os.environ.get("BIOPIPE_EMAIL", ""),
        help="联系邮箱（NCBI 要求；也可用环境变量 BIOPIPE_EMAIL）",
    )
    parser.add_argument("--api-key", default=os.environ.get("BIOPIPE_API_KEY", ""), help="NCBI API key")
    parser.add_argument("--json", action="store_true", help="以 JSON 输出完整记录")
    args = parser.parse_args(argv)

    if not args.email:
        LOGGER.error("必须提供 --email 或设置环境变量 BIOPIPE_EMAIL")
        return 2

    rows = search_bioprojects(args.accessions, email=args.email, api_key=args.api_key or None)
    if not rows:
        LOGGER.warning("没有查到任何记录")
        return 1

    if args.json:
        print(json.dumps(rows, ensure_ascii=False, indent=2))
    else:
        for row in rows:
            LOGGER.info(
                "%s | %s | %s | %s",
                row.get("accession", ""),
                row.get("organism", ""),
                row.get("title", "")[:60],
                (row.get("description", "") or "")[:80],
            )
    return 0


if __name__ == "__main__":
    sys.exit(_main())
