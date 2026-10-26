#!/usr/bin/env python3
"""PubChem 与 NCBI Taxonomy 的官方 API 客户端。

对应老仓库 ``crawl_database/scripts/`` 下的 6 个脚本：

- ``01_hit2lead_compound_search.py`` —— 爬公司内部 hit2lead 网站。
  **已删除**：内部商业系统，公开仓库里既跑不通也不该存在。
- ``02_ncbi_taxonomy_example.py`` / ``03_ncbi_taxonomy_batch.py`` ——
  用 Biopython 的 ``Entrez`` 查物种分类。
- ``04_pubchem_compound_info.py`` —— 用 ``pubchempy`` 查化合物信息。
- ``05_pubchem_cas_retriever.py`` —— 查 CAS 号。
- ``06_sid_to_iupac.py`` —— 用 IDEXCHANGE 查 SID → IUPAC 名。

本模块把其中 5 个脚本能公开的部分合并成两个纯 HTTP 客户端，直接调
官方 REST 接口，不依赖 ``pubchempy`` / ``Biopython``，也不需要浏览器。

老脚本的缺陷
------------
1. **伪造 NCBI email**。``03_ncbi_taxonomy_batch.py`` 准备了 5 个
   ``userN@example.com`` 并在每次请求时随机挑一个——Entrez 的 email 是
   身份标识，伪造既违规也会被限速封禁。
2. **多线程并发写同一个文件**。老脚本用 ``ThreadPoolExecutor(max_workers=2)``，
   每个线程各自 ``open(output_tsv, "a")`` 追加一行。两个线程的写入会
   交错，行可能被撕裂成半行。
3. **硬编码个人路径**。``C:\\Users\\DELL\\Desktop\\...`` 遍布全部脚本。
4. **CAS 号取错了地方**。老脚本从 ``xrefs/RegistryID`` 取"CAS"，但那个
   端点返回的是**几百个不同数据库的注册号**（``CHEBI:15365``、
   ``DB00945``、``LGCAMP0133.00-11``…），真正的 CAS 号混在里面。
   正确做法是从 ``synonyms`` 里按 CAS 的数字格式筛选。
5. **``CanonicalSMILES`` 字段名已变**。PubChem 现在把这个属性返回为
   ``SMILES``（或 ``ConnectivitySMILES``），按老名字取会得到空值。

PubChem 的三个 API 陷阱
-----------------------
1. **请求 ``CanonicalSMILES``，返回的键叫 ``SMILES``**。要拿到
   ``ConnectivitySMILES`` 就得请求 ``ConnectivitySMILES``。
2. **不存在的 CID 返回 HTTP 200**，内容是 ``{"CID": 999999999}``——
   只有 CID 字段、没有任何属性。必须检查属性是否齐全，不能只看状态码。
3. **``RegistryID`` 不是 CAS 号**，见上文第 4 条。

NCBI Taxonomy 的三个陷阱
-------------------------
1. **``LineageEx`` 不包含自身**，只含祖先。所以 species 层级要从顶层的
   ``<Rank>``/``<ScientificName>`` 取，不能只遍历 LineageEx。
2. **``LineageEx`` 里的 rank 远不止 8 个**。除了标准的
   superkingdom/kingdom/phylum/class/order/family/genus/species，还有
   domain、subphylum、superclass、clade、suborder、infraorder、parvorder、
   superfamily、subfamily 等几十种。老脚本只挑 8 个标准 rank，会漏掉
   domain（Eukaryota）这类信息量很大的层级。
3. **不能用 ``root.iter("Taxon")`` 遍历**。它会连``LineageEx`` 里的祖先
   节点一起收进来：查 3 个物种会得到 75 条记录，其中 72 条是
   "cellular organisms"、"Eukaryota" 这类祖先物种。必须只取根节点的
   **直接子** ``Taxon``（``root.findall("Taxon")``）。
"""

from __future__ import annotations

import re
import sys
import xml.etree.ElementTree as ET
from collections.abc import Sequence
from pathlib import Path
from typing import Any, Optional, Union

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "packages"))

from biopipe_core import get_logger, http_session  # noqa: E402

__all__ = [
    "PUG_REST_BASE",
    "EUTILS_BASE",
    "COMPOUND_PROPERTIES",
    "CAS_PATTERN",
    "fetch_compounds",
    "fetch_synonyms",
    "fetch_registry_ids",
    "extract_cas_numbers",
    "fetch_taxonomy",
    "parse_taxonomy_xml",
    "build_taxonomy_term",
    "ELEMENTARY_RANKS",
]

LOGGER = get_logger("biopipe.bio-database-fetch")

PUG_REST_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TAXONOMY_DB = "taxonomy"

#: 默认拉取的化合物属性。``ConnectivitySMILES`` 与 ``SMILES`` 都列出——
#: PubChem 会把 ``CanonicalSMILES`` 改名成 ``SMILES``，两个都要覆盖。
COMPOUND_PROPERTIES = (
    "MolecularFormula",
    "MolecularWeight",
    "SMILES",
    "ConnectivitySMILES",
    "InChI",
    "InChIKey",
    "IUPACName",
    "XLogP",
    "ExactMass",
    "TPSA",
    "HBondDonorCount",
    "HBondAcceptorCount",
    "RotatableBondCount",
    "HeavyAtomCount",
    "Charge",
)

#: CAS 号的标准格式：2-7 位数字 - 2 位数字 - 1 位数字，如 ``50-78-2``
CAS_PATTERN = re.compile(r"^\d{2,7}-\d{2}-\d$")

#: 分类学里最常用的"标准"层级，按从大到小排列。
#:
#: ``superkingdom`` 是历史遗留的 rank——现代 NCBI 已改用 ``cellular root``
#: 与 ``domain``，所以这个层级在当前数据里通常是空的。保留它是为了兼容
#: 更老的数据快照。
ELEMENTARY_RANKS = (
    "cellular root",
    "superkingdom",
    "domain",
    "kingdom",
    "subkingdom",
    "phylum",
    "subphylum",
    "superclass",
    "class",
    "subclass",
    "order",
    "suborder",
    "infraorder",
    "parvorder",
    "superfamily",
    "family",
    "subfamily",
    "tribe",
    "genus",
    "subgenus",
    "species",
    "subspecies",
)


def _json_get(url: str, session, timeout: int = 60) -> Any:
    """发一个 GET 请求并解析 JSON。"""
    response = session.get(url, timeout=timeout)
    response.raise_for_status()
    return response.json()


def fetch_compounds(
    cids: Sequence[Union[int, str]],
    properties: Sequence[str] = COMPOUND_PROPERTIES,
    session=None,
    timeout: int = 60,
    chunk_size: int = 100,
) -> list[dict]:
    """按 CID 批量获取化合物属性。

    PubChem 的 GET 接口单次最多接受 100 个 CID，这里自动分块。

    **注意**：不存在的 CID 不会报错——PubChem 返回 HTTP 200，内容是
    ``{"CID": 999999999}``，只有 CID 字段没有属性。因此每条记录都带
    ``found`` 标记，调用方应当检查它。

    Args:
        cids: CID 列表。
        properties: 要拉取的属性名。
        session: 可复用的 ``requests.Session``。
        timeout: 单次请求超时秒数。
        chunk_size: 单次请求的 CID 数量上限。

    Returns:
        每个 CID 一条记录的列表（含 ``found`` 标记）。
    """
    if not cids:
        return []

    own_session = session is None
    if own_session:
        # PubChem 的限速是 5 req/s、400 req/min
        session = http_session(min_interval=0.25, retries=4)
    if session is None:
        LOGGER.warning("requests 不可用，无法查询 PubChem")
        return []

    normalized = [str(cid).strip() for cid in cids if str(cid).strip()]
    property_list = ",".join(properties)
    records: list[dict] = []

    try:
        for start in range(0, len(normalized), chunk_size):
            chunk = normalized[start : start + chunk_size]
            url = (
                f"{PUG_REST_BASE}/compound/cid/{','.join(chunk)}"
                f"/property/{property_list}/JSON"
            )
            payload = _json_get(url, session, timeout)
            table = payload.get("PropertyTable", {}).get("Properties", [])

            by_cid = {str(item.get("CID")): item for item in table}
            for cid in chunk:
                item = by_cid.get(cid)
                if item is None:
                    records.append({"CID": cid, "found": False})
                    continue
                # PubChem 对不存在的 CID 只返回 CID 字段——据此判定是否找到
                found = len(item) > 1
                record = dict(item)
                # CID 统一成字符串：PubChem 的 JSON 里 CID 是数字，而请求侧
                # 始终是字符串。不统一的话，同一列里会混着 int 和 str，
                # 下游按 CID 归并时会漏掉一半记录。
                record["CID"] = cid
                record["found"] = found
                records.append(record)
    except Exception as exc:  # noqa: BLE001
        LOGGER.warning("PubChem 查询失败: %s", exc)
    finally:
        if own_session:
            session.close()

    return records


def fetch_synonyms(
    cids: Sequence[Union[int, str]],
    session=None,
    timeout: int = 60,
    chunk_size: int = 100,
) -> dict[str, list[str]]:
    """获取化合物的同义词列表。

    Args:
        cids: CID 列表。
        session: 可复用的 ``requests.Session``。
        timeout: 单次请求超时秒数。
        chunk_size: 单次请求的 CID 数量上限。

    Returns:
        ``{CID: [同义词, ...]}`` 字典；没有同义词的 CID 不在字典里。
    """
    if not cids:
        return {}

    own_session = session is None
    if own_session:
        session = http_session(min_interval=0.25, retries=4)
    if session is None:
        return {}

    normalized = [str(cid).strip() for cid in cids if str(cid).strip()]
    result: dict[str, list[str]] = {}

    try:
        for start in range(0, len(normalized), chunk_size):
            chunk = normalized[start : start + chunk_size]
            url = (
                f"{PUG_REST_BASE}/compound/cid/{','.join(chunk)}"
                f"/synonyms/JSON"
            )
            payload = _json_get(url, session, timeout)
            for item in payload.get("InformationList", {}).get("Information", []):
                cid = str(item.get("CID", ""))
                synonyms = item.get("Synonym") or []
                if cid:
                    result[cid] = [str(s) for s in synonyms]
    except Exception as exc:  # noqa: BLE001
        LOGGER.warning("PubChem 同义词查询失败: %s", exc)
    finally:
        if own_session:
            session.close()

    return result


def fetch_registry_ids(
    cids: Sequence[Union[int, str]],
    session=None,
    timeout: int = 60,
    chunk_size: int = 100,
) -> dict[str, list[str]]:
    """获取化合物的 RegistryID 列表。

    **注意：这不是 CAS 号列表。** RegistryID 里混着几百个不同数据库的
    注册号（``CHEBI:15365``、``DB00945``、``LGCAMP0133.00-11``…），
    CAS 号只是其中符合数字格式的少数几个。要拿 CAS 请用
    :func:`extract_cas_numbers` 配合 :func:`fetch_synonyms`。

    Args:
        cids: CID 列表。
        session: 可复用的 ``requests.Session``。
        timeout: 单次请求超时秒数。
        chunk_size: 单次请求的 CID 数量上限。

    Returns:
        ``{CID: [RegistryID, ...]}`` 字典。
    """
    if not cids:
        return {}

    own_session = session is None
    if own_session:
        session = http_session(min_interval=0.25, retries=4)
    if session is None:
        return {}

    normalized = [str(cid).strip() for cid in cids if str(cid).strip()]
    result: dict[str, list[str]] = {}

    try:
        for start in range(0, len(normalized), chunk_size):
            chunk = normalized[start : start + chunk_size]
            url = (
                f"{PUG_REST_BASE}/compound/cid/{','.join(chunk)}"
                f"/xrefs/RegistryID/JSON"
            )
            payload = _json_get(url, session, timeout)
            for item in payload.get("InformationList", {}).get("Information", []):
                cid = str(item.get("CID", ""))
                ids = item.get("RegistryID") or []
                if cid:
                    result[cid] = [str(s) for s in ids]
    except Exception as exc:  # noqa: BLE001
        LOGGER.warning("PubChem RegistryID 查询失败: %s", exc)
    finally:
        if own_session:
            session.close()

    return result


def extract_cas_numbers(values: Sequence[str]) -> list[str]:
    """从同义词或 RegistryID 列表里筛出 CAS 号。

    CAS 号的格式是固定的：``NNNNN-NN-N``（末段校验位，2-7 位数字前缀）。
    按这个模式筛选能把 ``CHEBI:15365``、``DB00945``、``50-78-2`` 这类
    混在一起的条目分开。

    一个化合物可能有多个 CAS 号（例如同一物质的不同盐形式），因此返回
    列表而非单个值。

    Args:
        values: 待筛选的字符串列表。

    Returns:
        去重后的 CAS 号列表（保持原始顺序）。
    """
    seen: set[str] = set()
    cas_numbers: list[str] = []
    for value in values:
        candidate = str(value).strip()
        if CAS_PATTERN.match(candidate) and candidate not in seen:
            seen.add(candidate)
            cas_numbers.append(candidate)
    return cas_numbers


def build_taxonomy_term(tax_ids: Sequence[Union[int, str]]) -> str:
    """构造 Taxonomy 的检索式。

    Args:
        tax_ids: NCBI Taxonomy ID 列表。

    Returns:
        检索式字符串；为空时返回空串。
    """
    cleaned = [str(tid).strip() for tid in tax_ids if str(tid).strip()]
    if not cleaned:
        return ""
    return " OR ".join(f"{tid}[Taxonomy ID]" for tid in cleaned)


def parse_taxonomy_xml(xml_text: str) -> list[dict]:
    """把 efetch 返回的 Taxonomy XML 解析成结构化记录。

    每条记录包含 taxid、学名、俗名、各级分类名，以及完整 lineage。

    Args:
        xml_text: efetch 返回的 XML 文本。

    Returns:
        每个 taxon 一条记录的列表。
    """
    if not xml_text.strip():
        return []

    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError as exc:
        LOGGER.warning("Taxonomy XML 解析失败: %s", exc)
        return []

    records: list[dict] = []
    # 只取根节点下的**直接子** Taxon。用 root.iter("Taxon") 会把每个
    # 记录 LineageEx 里的祖先节点也当成独立记录——查 3 个物种会返回
    # 75 条，其中 72 条是 "cellular organisms"、"Eukaryota" 这类祖先。
    for taxon in root.findall("Taxon"):
        record = _parse_taxon(taxon)
        if record:
            records.append(record)
    return records


def _parse_taxon(taxon: ET.Element) -> dict:
    """解析单个 Taxon 元素。"""
    tax_id = _find_text(taxon, "TaxId")
    scientific_name = _find_text(taxon, "ScientificName")
    rank = _find_text(taxon, "Rank")

    # 俗名（GenbankCommonName）
    common_name = ""
    genbank = taxon.find("OtherNames/GenbankCommonName")
    if genbank is not None and genbank.text:
        common_name = genbank.text.strip()

    # LineageEx 只含**祖先**，不含自身。所以自身的 rank（通常是 species）
    # 必须从顶层元素取，不能只遍历 LineageEx。
    lineage: list[dict] = []
    lineage_ex = taxon.find("LineageEx")
    if lineage_ex is not None:
        for entry in lineage_ex.findall("Taxon"):
            lineage.append(
                {
                    "tax_id": _find_text(entry, "TaxId"),
                    "rank": _find_text(entry, "Rank").lower(),
                    "name": _find_text(entry, "ScientificName"),
                }
            )

    # 按 rank 归类。LineageEx 里有几十种 rank（含 clade、superclass 等），
    # 全部保留——只挑 8 个标准 rank 会漏掉 domain 这类信息量大的层级。
    by_rank: dict[str, str] = {}
    for item in lineage:
        if item["rank"] and item["rank"] not in by_rank:
            by_rank[item["rank"]] = item["name"]

    # 自身也计入 by_rank（species 常常只在这里出现）
    if rank and rank.lower() not in by_rank:
        by_rank[rank.lower()] = scientific_name

    record = {
        "tax_id": tax_id,
        "scientific_name": scientific_name,
        "common_name": common_name,
        "rank": rank.lower(),
        "by_rank": by_rank,
        "lineage": lineage,
        # lineage 深度 = 祖先数 + 1（自身）。自身不在 LineageEx 里，
        # 所以要加上，否则 depth 会少算一级。
        "lineage_depth": len(lineage) + 1 if lineage else 1,
    }
    return record


def _find_text(element: ET.Element, tag: str) -> str:
    """在元素（含子元素）中查找第一个匹配标签的文本。"""
    node = element.find(tag)
    if node is None or node.text is None:
        return ""
    return node.text.strip()


def fetch_taxonomy(
    tax_ids: Sequence[Union[int, str]],
    email: str,
    api_key: Optional[str] = None,
    session=None,
    timeout: int = 60,
    chunk_size: int = 100,
) -> list[dict]:
    """按 Taxonomy ID 批量获取分类信息。

    流程：``efetch`` 直接接受逗号分隔的 Taxonomy ID，不需要先走
    esearch（这与 BioProject 不同——Taxonomy 的 ID 本身就是主键）。

    **email 是必填的**：NCBI 要求调用方标识身份。老脚本用 5 个随机
    生成的 ``userN@example.com`` 轮流冒充，既违规也拿不到更高的限速。

    Args:
        tax_ids: Taxonomy ID 列表。
        email: 联系邮箱（NCBI 要求）。
        api_key: NCBI API key（可选，限速从 3 req/s 提到 10 req/s）。
        session: 可复用的 ``requests.Session``。
        timeout: 单次请求超时秒数。
        chunk_size: 单次请求的 ID 数量上限。

    Returns:
        每个 taxon 一条记录的列表。
    """
    if not tax_ids:
        return []
    if not email:
        raise ValueError(
            "必须提供 email —— NCBI 要求调用方标识身份，"
            "可用 --email 传入或设置环境变量 BIOPIPE_EMAIL"
        )

    own_session = session is None
    if own_session:
        session = http_session(min_interval=0.34, retries=4)
    if session is None:
        LOGGER.warning("requests 不可用，无法查询 NCBI")
        return []

    normalized = [str(tid).strip() for tid in tax_ids if str(tid).strip()]
    records: list[dict] = []

    try:
        for start in range(0, len(normalized), chunk_size):
            chunk = normalized[start : start + chunk_size]
            params: dict[str, Any] = {
                "db": TAXONOMY_DB,
                "id": ",".join(chunk),
                "retmode": "xml",
                "email": email,
            }
            if api_key:
                params["api_key"] = api_key

            response = session.get(
                f"{EUTILS_BASE}/efetch.fcgi", params=params, timeout=timeout
            )
            response.raise_for_status()
            records.extend(parse_taxonomy_xml(response.text))
    except ValueError:
        raise
    except Exception as exc:  # noqa: BLE001
        LOGGER.warning("NCBI Taxonomy 查询失败: %s", exc)
    finally:
        if own_session:
            session.close()

    return records
