#!/usr/bin/env python3
"""文献检索核心库：多源 API 客户端、归一化、去重、相关性打分、摘要富集。

本模块被 ``search_literature.py`` 与 ``export_citations.py`` 复用，
所有网络访问都通过共享层 :func:`biopipe_core.http_session` 完成（自带重试与限速）。

数据源与端点
------------
- PubMed   : ``https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi``
             + ``esummary.fcgi`` + ``efetch.fcgi``
- OpenAlex : ``https://api.openalex.org/works``
- Semantic Scholar : ``https://api.semanticscholar.org/graph/v1/paper/search``

统一记录 schema（``dict``）
---------------------------
``title`` ``authors`` ``year`` ``doi`` ``journal`` ``abstract`` ``citations``
``source`` ``pmid`` ``url`` ``type`` ``oa_status`` ``query`` ``relevance``
"""

from __future__ import annotations

import difflib
import json
import re
import xml.etree.ElementTree as ET
from collections.abc import Iterable, Sequence
from typing import Any

from biopipe_core import SkillError, get_logger, http_session

LOGGER = get_logger("biopipe.literature-search")

# ---------------------------------------------------------------- API 端点

PUBMED_ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
PUBMED_ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
PUBMED_EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
OPENALEX_WORKS_URL = "https://api.openalex.org/works"
OPENALEX_WORK_URL = "https://api.openalex.org/works"
S2_SEARCH_URL = "https://api.semanticscholar.org/graph/v1/paper/search"
S2_PAPER_URL = "https://api.semanticscholar.org/graph/v1/paper/DOI:"

SOURCE_LABELS = {
    "pubmed": "PubMed",
    "openalex": "OpenAlex",
    "semantic-scholar": "Semantic Scholar",
}

# 统一输出列
CSV_FIELDS = [
    "title",
    "authors",
    "year",
    "doi",
    "journal",
    "abstract",
    "citations",
    "source",
    "pmid",
    "url",
    "type",
    "oa_status",
    "relevance",
    "query",
]


# ---------------------------------------------------------------- 工具函数


def normalize_doi(doi: str | None) -> str:
    """把各种写法的 DOI 归一化为裸小写 DOI。

    ``https://doi.org/10.1234/ABC`` -> ``10.1234/abc``
    """
    if not doi:
        return ""
    text = str(doi).strip().lower()
    for prefix in ("https://doi.org/", "http://doi.org/", "https://dx.doi.org/", "doi:"):
        if text.startswith(prefix):
            text = text[len(prefix) :]
            break
    text = text.strip().strip(".").strip()
    return text if text.startswith("10.") else ""


def normalize_title(title: str | None) -> str:
    """标题归一化：转小写、去掉标点与多余空白，用于去重与模糊匹配。"""
    if not title:
        return ""
    text = str(title).lower()
    text = re.sub(r"[^a-z0-9\u4e00-\u9fff]+", " ", text)
    return re.sub(r"\s+", " ", text).strip()


def clean_text(text: str | None) -> str:
    """压缩空白并去掉 HTML 残留。"""
    if not text:
        return ""
    text = re.sub(r"<[^>]+>", " ", str(text))
    text = re.sub(r"&[a-zA-Z]+;", " ", text)
    return re.sub(r"\s+", " ", text).strip()


def parse_year(value: Any) -> int | None:
    """从多种字段形态中解析出 4 位年份。"""
    if value is None:
        return None
    if isinstance(value, int):
        return value if 1000 <= value <= 2999 else None
    text = str(value).strip()
    match = re.search(r"(1[6-9]\d{2}|20\d{2}|21\d{2})", text)
    if match:
        return int(match.group(1))
    return None


def parse_date_range(spec: str | None) -> tuple[int | None, int | None]:
    """解析 ``--date-range`` 参数。

    支持 ``START:END``（如 ``2018:2024``）、``START:``（至今）、``:END``。
    """
    if not spec:
        return None, None
    text = spec.strip()
    if ":" not in text:
        raise SkillError(f"--date-range 格式应为 START:END，收到 {spec!r}")
    start_text, end_text = text.split(":", 1)
    start = parse_year(start_text) if start_text.strip() else None
    end = parse_year(end_text) if end_text.strip() else None
    if start and end and start > end:
        raise SkillError(f"--date-range 起始年份 {start} 大于结束年份 {end}")
    return start, end


def build_query(topic: str | None, keywords: Sequence[str] | None) -> str:
    """把主题与关键词拼成一条检索式。

    ``--topic`` 视为自由文本主检索式；``--keywords`` 作为附加 AND 条件，
    使检索可以收窄到具体子方向。
    """
    parts: list[str] = []
    if topic:
        parts.append(f"({topic})")
    if keywords:
        joined = " OR ".join(f'"{kw}"' for kw in keywords if kw.strip())
        if joined:
            parts.append(f"({joined})")
    return " ".join(parts) if parts else ""


def _author_list(raw: Any) -> list[str]:
    """把各数据源的作者字段统一成字符串列表。"""
    authors: list[str] = []
    if not raw:
        return authors
    if isinstance(raw, str):
        return [a.strip() for a in re.split(r"[;,]| and ", raw) if a.strip()]
    if isinstance(raw, dict):
        raw = raw.get("authors") or raw.get("author") or []
    if isinstance(raw, Sequence):
        for item in raw:
            if isinstance(item, str):
                if item.strip():
                    authors.append(item.strip())
            elif isinstance(item, dict):
                name = (
                    item.get("name")
                    or item.get("display_name")
                    or item.get("fullName")
                    or (item.get("author") or {}).get("display_name")
                    or ""
                )
                if name:
                    authors.append(str(name).strip())
    return authors


# ---------------------------------------------------------------- 记录构造


def make_record(
    *,
    title: str,
    authors: Any = None,
    year: Any = None,
    doi: str = "",
    journal: str = "",
    abstract: str = "",
    citations: Any = 0,
    source: str = "",
    pmid: str = "",
    url: str = "",
    doc_type: str = "",
    oa_status: str = "",
    query: str = "",
) -> dict[str, Any]:
    """构造统一 schema 的一条文献记录。"""
    doi_norm = normalize_doi(doi)
    title_clean = clean_text(title)
    return {
        "title": title_clean,
        "authors": "; ".join(_author_list(authors)),
        "year": parse_year(year),
        "doi": doi_norm,
        "journal": clean_text(journal),
        "abstract": clean_text(abstract),
        "citations": int(citations or 0),
        "source": source,
        "pmid": str(pmid or ""),
        "url": url or (f"https://doi.org/{doi_norm}" if doi_norm else ""),
        "type": doc_type,
        "oa_status": oa_status,
        "query": query,
        "relevance": 0.0,
    }


# ---------------------------------------------------------------- API 客户端


class LiteratureClient:
    """封装三个公开文献 API 的查询逻辑（带重试、限速与 API Key 支持）。"""

    def __init__(
        self,
        email: str | None = None,
        ncbi_api_key: str | None = None,
        s2_api_key: str | None = None,
        min_interval: float = 0.34,
    ) -> None:
        self.email = email
        self.ncbi_api_key = ncbi_api_key or ""
        self.s2_api_key = s2_api_key or ""
        # Semantic Scholar 无 key 时限速更严（约 1 req/s）
        self._sessions: dict[str, Any] = {}
        self._min_interval = min_interval
        self._s2_min_interval = min_interval if self.s2_api_key else max(min_interval, 1.1)

    # -- session 管理 ----------------------------------------------------
    def _session(self, source: str):
        if source not in self._sessions:
            headers = {}
            if source == "semantic-scholar" and self.s2_api_key:
                headers["x-api-key"] = self.s2_api_key
            self._sessions[source] = http_session(
                min_interval=self._s2_min_interval if source == "semantic-scholar" else self._min_interval,
                extra_headers=headers or None,
            )
        return self._sessions[source]

    def close(self) -> None:
        """关闭所有会话。"""
        for session in self._sessions.values():
            if session is not None:
                session.close()
        self._sessions.clear()

    def _get(self, source: str, url: str, params: dict[str, Any]) -> Any:
        session = self._session(source)
        if session is None:
            raise SkillError("requests 未安装，无法访问外部 API（pip install requests）")
        response = session.get(url, params=params, timeout=60)
        response.raise_for_status()
        return response.json()

    # -- PubMed ----------------------------------------------------------
    def search_pubmed(
        self,
        query: str,
        max_results: int = 50,
        start_year: int | None = None,
        end_year: int | None = None,
    ) -> list[dict[str, Any]]:
        """PubMed E-utilities 检索：esearch 取 PMID，再 esummary 取元数据。"""
        term = query
        if start_year or end_year:
            low = start_year or 1800
            high = end_year or 3000
            term = f"({query}) AND (\"{low}\"[PDAT] : \"{high}\"[PDAT])"

        params: dict[str, Any] = {
            "db": "pubmed",
            "term": term,
            "retmax": min(int(max_results), 10000),
            "retmode": "json",
            "sort": "relevance",
        }
        if self.ncbi_api_key:
            params["api_key"] = self.ncbi_api_key
            params["tool"] = "biopipe"
        if self.email:
            params["email"] = self.email

        payload = self._get("pubmed", PUBMED_ESEARCH_URL, params)
        id_list = payload.get("esearchresult", {}).get("idlist", []) or []
        if not id_list:
            LOGGER.info("PubMed 未命中: %s", query)
            return []

        LOGGER.info("PubMed 命中 %d 条: %s", len(id_list), query)
        return self._pubmed_summaries(id_list, query)

    def _pubmed_summaries(self, pmids: list[str], query: str) -> list[dict[str, Any]]:
        """批量 esummary：每次最多 200 个 PMID。"""
        records: list[dict[str, Any]] = []
        for start in range(0, len(pmids), 200):
            chunk = pmids[start : start + 200]
            params: dict[str, Any] = {
                "db": "pubmed",
                "id": ",".join(chunk),
                "retmode": "json",
            }
            if self.ncbi_api_key:
                params["api_key"] = self.ncbi_api_key
            payload = self._get("pubmed", PUBMED_ESUMMARY_URL, params)
            result = payload.get("result", {}) or {}
            for pmid in chunk:
                doc = result.get(pmid)
                if not isinstance(doc, dict):
                    continue
                doi = ""
                for aid in doc.get("articleids", []) or []:
                    if aid.get("idtype") == "doi":
                        doi = aid.get("value") or aid.get("id") or ""
                        break
                pubtype = doc.get("pubtype", []) or []
                records.append(
                    make_record(
                        title=doc.get("title", ""),
                        authors=doc.get("authors", []),
                        year=doc.get("pubdate", ""),
                        doi=doi,
                        journal=doc.get("source", ""),
                        abstract="",  # esummary 不含摘要，由 enrich 阶段补全
                        citations=doc.get("pmcrefcount", 0) or 0,
                        source=SOURCE_LABELS["pubmed"],
                        pmid=pmid,
                        url=f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                        doc_type="; ".join(pubtype),
                        oa_status=("Y" if doc.get("pubstatus", "") == "pubmed" else ""),
                        query=query,
                    )
                )
        return records

    def fetch_pubmed_abstracts(self, pmids: Sequence[str]) -> dict[str, str]:
        """用 efetch(abstract, xml) 批量取摘要，返回 ``pmid -> abstract``。"""
        abstracts: dict[str, str] = {}
        pmids = [p for p in pmids if p]
        if not pmids:
            return abstracts
        for start in range(0, len(pmids), 100):
            chunk = pmids[start : start + 100]
            params: dict[str, Any] = {
                "db": "pubmed",
                "id": ",".join(chunk),
                "rettype": "abstract",
                "retmode": "xml",
            }
            if self.ncbi_api_key:
                params["api_key"] = self.ncbi_api_key
            session = self._session("pubmed")
            if session is None:
                raise SkillError("requests 未安装，无法访问 PubMed")
            response = session.get(PUBMED_EFETCH_URL, params=params, timeout=60)
            response.raise_for_status()
            abstracts.update(self._parse_pubmed_abstracts(response.text))
        return abstracts

    @staticmethod
    def _parse_pubmed_abstracts(xml_text: str) -> dict[str, str]:
        """解析 efetch XML，抽取每个 PubMedArticle 的摘要与 DOI。"""
        result: dict[str, str] = {}
        try:
            root = ET.fromstring(xml_text)
        except ET.ParseError as exc:
            LOGGER.warning("PubMed 摘要 XML 解析失败: %s", exc)
            return result
        for article in root.iter("PubMedArticle"):
            pmid_node = article.find(".//MedlineCitation/PMID")
            if pmid_node is None or not pmid_node.text:
                continue
            pmid = pmid_node.text.strip()
            parts: list[str] = []
            for node in article.findall(".//Abstract/AbstractText"):
                label = node.get("Label")
                text = clean_text("".join(node.itertext()))
                if not text:
                    continue
                parts.append(f"{label}: {text}" if label else text)
            result[pmid] = " ".join(parts)
        return result

    # -- OpenAlex --------------------------------------------------------
    def search_openalex(
        self,
        query: str,
        max_results: int = 50,
        start_year: int | None = None,
        end_year: int | None = None,
        min_citations: int = 0,
    ) -> list[dict[str, Any]]:
        """OpenAlex works 检索，返回归一化记录。"""
        filters: list[str] = []
        if start_year or end_year:
            low = start_year or 1800
            high = end_year or 3000
            filters.append(f"from_publication_date:{low}-01-01,to_publication_date:{high}-12-31")
        if min_citations > 0:
            filters.append(f"cited_by_count:>{min_citations}")
        filters.append("is_retracted:false")

        params: dict[str, Any] = {
            "search": query,
            "per-page": min(int(max_results), 200),
            "select": ",".join(
                [
                    "id",
                    "doi",
                    "title",
                    "display_name",
                    "publication_year",
                    "publication_date",
                    "cited_by_count",
                    "type",
                    "authorships",
                    "primary_location",
                    "open_access",
                    "abstract_inverted_index",
                    "best_oa_location",
                ]
            ),
        }
        if filters:
            params["filter"] = ",".join(filters)
        if self.email:
            params["mailto"] = self.email

        payload = self._get("openalex", OPENALEX_WORKS_URL, params)
        results = payload.get("results", []) or []
        LOGGER.info("OpenAlex 命中 %d 条: %s", len(results), query)

        records: list[dict[str, Any]] = []
        for work in results:
            location = work.get("primary_location") or {}
            source_info = location.get("source") or {}
            oa = work.get("open_access") or {}
            records.append(
                make_record(
                    title=work.get("title") or work.get("display_name") or "",
                    authors=work.get("authorships", []),
                    year=work.get("publication_year"),
                    doi=work.get("doi", ""),
                    journal=source_info.get("display_name", ""),
                    abstract=reconstruct_abstract(work.get("abstract_inverted_index")),
                    citations=work.get("cited_by_count", 0),
                    source=SOURCE_LABELS["openalex"],
                    url=work.get("id", "") or (location.get("landing_page_url", "")),
                    doc_type=work.get("type", ""),
                    oa_status=oa.get("oa_status", "") or "",
                    query=query,
                )
            )
        return records

    def fetch_openalex_abstract(self, doi: str) -> str:
        """按 DOI 从 OpenAlex 单条取摘要（用于摘要富集）。"""
        doi_norm = normalize_doi(doi)
        if not doi_norm:
            return ""
        params: dict[str, Any] = {}
        if self.email:
            params["mailto"] = self.email
        try:
            payload = self._get(
                "openalex",
                f"{OPENALEX_WORK_URL}/https://doi.org/{doi_norm}",
                params,
            )
        except Exception as exc:  # noqa: BLE001 - 摘要富集为尽力而为
            LOGGER.debug("OpenAlex 摘要获取失败 (%s): %s", doi_norm, exc)
            return ""
        return reconstruct_abstract(payload.get("abstract_inverted_index"))

    # -- Semantic Scholar ------------------------------------------------
    def search_semantic_scholar(
        self,
        query: str,
        max_results: int = 50,
        start_year: int | None = None,
        end_year: int | None = None,
    ) -> list[dict[str, Any]]:
        """Semantic Scholar graph search，返回归一化记录。"""
        params: dict[str, Any] = {
            "query": query,
            "limit": min(int(max_results), 100),
            "fields": ",".join(
                [
                    "paperId",
                    "title",
                    "authors",
                    "year",
                    "venue",
                    "publicationVenue",
                    "publicationTypes",
                    "citationCount",
                    "isOpenAccess",
                    "openAccessPdf",
                    "abstract",
                    "externalIds",
                    "url",
                ]
            ),
        }
        if start_year or end_year:
            low = start_year or 1000
            high = end_year or 3000
            params["year"] = f"{low}-{high}"

        payload = self._get("semantic-scholar", S2_SEARCH_URL, params)
        papers = payload.get("data", []) or []
        LOGGER.info("Semantic Scholar 命中 %d 条: %s", len(papers), query)

        records: list[dict[str, Any]] = []
        for paper in papers:
            ext = paper.get("externalIds") or {}
            venue = paper.get("venue") or (paper.get("publicationVenue") or {}).get("name", "")
            oa = "Y" if paper.get("isOpenAccess") else "N"
            records.append(
                make_record(
                    title=paper.get("title", ""),
                    authors=paper.get("authors", []),
                    year=paper.get("year"),
                    doi=ext.get("DOI", ""),
                    journal=venue,
                    abstract=paper.get("abstract", ""),
                    citations=paper.get("citationCount", 0),
                    source=SOURCE_LABELS["semantic-scholar"],
                    pmid=str(ext.get("PMID", "") or ""),
                    url=paper.get("url", "") or paper.get("paperId", ""),
                    doc_type="; ".join(paper.get("publicationTypes", []) or []),
                    oa_status=oa,
                    query=query,
                )
            )
        return records

    def fetch_s2_abstract(self, doi: str) -> str:
        """按 DOI 从 Semantic Scholar 取摘要（无摘要时回退到 tldr）。"""
        doi_norm = normalize_doi(doi)
        if not doi_norm:
            return ""
        try:
            payload = self._get(
                "semantic-scholar",
                f"{S2_PAPER_URL}{doi_norm}",
                {"fields": "abstract,tldr"},
            )
        except Exception as exc:  # noqa: BLE001
            LOGGER.debug("S2 摘要获取失败 (%s): %s", doi_norm, exc)
            return ""
        abstract = clean_text(payload.get("abstract", ""))
        if abstract:
            return abstract
        tldr = payload.get("tldr") or {}
        return clean_text(tldr.get("text", ""))


def reconstruct_abstract(inverted_index: Any) -> str:
    """从 OpenAlex 的倒排索引重建摘要文本。

    ``{"word": [pos0, pos1, ...]}`` -> 正常句子。
    """
    if not inverted_index or not isinstance(inverted_index, dict):
        return ""
    positions: dict[int, str] = {}
    max_pos = -1
    for word, word_positions in inverted_index.items():
        if not isinstance(word_positions, list):
            continue
        for pos in word_positions:
            try:
                pos_int = int(pos)
            except (TypeError, ValueError):
                continue
            positions[pos_int] = str(word)
            max_pos = max(max_pos, pos_int)
    if max_pos < 0:
        return ""
    return clean_text(" ".join(positions.get(i, "") for i in range(max_pos + 1)))


# ---------------------------------------------------------------- 去重


def deduplicate(
    records: Sequence[dict[str, Any]],
    title_threshold: float = 0.92,
) -> tuple[list[dict[str, Any]], int]:
    """两级去重：DOI 精确匹配优先，其次标题模糊匹配。

    重要约束：**两条记录都带 DOI 且 DOI 不同时，判定为不同文献**。
    学术出版中不同期刊常刊出同名文章（如多篇都叫 "Antimicrobial Peptides"
    的综述），仅凭标题相似就合并会造成严重的文献丢失。
    标题模糊匹配只在以下情况生效：

    - 至少有一条记录缺少 DOI；
    - 或两条记录的 DOI 相同（此时已由第一级规则覆盖）。

    Returns:
        ``(去重后记录列表, 被移除的条数)``
    """
    kept: list[dict[str, Any]] = []
    kept_dois: list[str] = []
    kept_titles: list[str] = []
    removed = 0

    for record in records:
        doi = str(record.get("doi", "") or "").strip().lower()
        norm_title = normalize_title(record.get("title", ""))

        # ---- 第一级：DOI 精确匹配 ----
        if doi and doi in kept_dois:
            removed += 1
            _merge_into(kept, kept_dois.index(doi), record)
            continue

        # ---- 第二级：标题模糊匹配（仅在 DOI 不足以区分时才生效）----
        duplicate_index = -1
        for idx, existing_title in enumerate(kept_titles):
            if not existing_title or not norm_title:
                continue
            existing_doi = kept_dois[idx]
            if doi and existing_doi and doi != existing_doi:
                # 双方都有不同 DOI —— 视为不同文献，跳过标题比较
                continue
            if existing_title == norm_title:
                duplicate_index = idx
                break
            if difflib.SequenceMatcher(None, existing_title, norm_title).ratio() >= title_threshold:
                duplicate_index = idx
                break

        if duplicate_index >= 0:
            removed += 1
            _merge_into(kept, duplicate_index, record)
            continue

        kept.append(record)
        kept_dois.append(doi)
        kept_titles.append(norm_title)

    return kept, removed


def _merge_into(kept: list[dict[str, Any]], index: int, record: dict[str, Any]) -> None:
    """把被去重掉的记录信息合并到已有记录（补摘要、累计引用、标记多来源）。"""
    if index < 0 or index >= len(kept):
        return
    existing = kept[index]
    if not existing.get("abstract") and record.get("abstract"):
        existing["abstract"] = record["abstract"]
    existing["citations"] = max(int(existing.get("citations", 0)), int(record.get("citations", 0)))
    sources = {existing.get("source", ""), record.get("source", "")} - {""}
    if len(sources) > 1:
        existing["source"] = "+".join(sorted(sources))
    # 补齐缺失的标识符
    for key in ("pmid", "url"):
        if not existing.get(key) and record.get(key):
            existing[key] = record[key]


# ---------------------------------------------------------------- 相关性打分


def score_relevance(
    record: dict[str, Any],
    keywords: Sequence[str],
    title_weight: float = 3.0,
    abstract_weight: float = 1.0,
) -> float:
    """相关性打分。

    - 关键词命中标题权重更高
    - 引用数做对数饱和加成（避免老论文碾压新论文）
    - 近 5 年论文小幅加成
    """
    title = str(record.get("title", "")).lower()
    abstract = str(record.get("abstract", "")).lower()
    score = 0.0

    for kw in keywords:
        needle = kw.lower().strip()
        if not needle:
            continue
        if needle in title:
            score += title_weight
        elif needle in abstract:
            score += abstract_weight

    citations = int(record.get("citations", 0) or 0)
    if citations > 0:
        score += min(2.0, (citations ** 0.25) * 0.6)

    year = record.get("year")
    if isinstance(year, int) and year >= 2021:
        score += 0.5

    return round(score, 3)


def filter_records(
    records: Sequence[dict[str, Any]],
    *,
    min_relevance: float = 0.0,
    min_citations: int = 0,
    min_abstract_length: int = 0,
    include_terms: Sequence[str] = (),
    exclude_terms: Sequence[str] = (),
    start_year: int | None = None,
    end_year: int | None = None,
) -> list[dict[str, Any]]:
    """按阈值与包含/排除词筛选文献。"""
    include_lower = [t.lower() for t in include_terms if t.strip()]
    exclude_lower = [t.lower() for t in exclude_terms if t.strip()]

    selected: list[dict[str, Any]] = []
    for record in records:
        if float(record.get("relevance", 0.0)) < min_relevance:
            continue
        if int(record.get("citations", 0) or 0) < min_citations:
            continue
        if min_abstract_length and len(str(record.get("abstract", ""))) < min_abstract_length:
            continue
        year = record.get("year")
        if start_year and isinstance(year, int) and year < start_year:
            continue
        if end_year and isinstance(year, int) and year > end_year:
            continue

        haystack = " ".join(
            [
                str(record.get("title", "")),
                str(record.get("abstract", "")),
                str(record.get("journal", "")),
                str(record.get("keywords_text", "")),
            ]
        ).lower()

        if include_lower and not any(term in haystack for term in include_lower):
            continue
        if exclude_lower and any(term in haystack for term in exclude_lower):
            continue
        selected.append(record)
    return selected


# ---------------------------------------------------------------- 摘要富集


def enrich_abstracts(
    records: list[dict[str, Any]],
    client: LiteratureClient,
    min_length: int = 60,
) -> int:
    """为缺摘要的记录补全摘要。

    依次尝试：PubMed（按 PMID）-> OpenAlex（按 DOI）-> Semantic Scholar（按 DOI）。

    Returns:
        被成功补全摘要的条数。
    """
    missing = [
        r
        for r in records
        if len(str(r.get("abstract", ""))) < min_length
    ]
    if not missing:
        return 0

    LOGGER.info("需要富化摘要: %d 条", len(missing))

    pmids = [r["pmid"] for r in missing if r.get("pmid")]
    if pmids:
        try:
            abstracts = client.fetch_pubmed_abstracts(pmids)
        except Exception as exc:  # noqa: BLE001
            LOGGER.warning("PubMed 摘要批量获取失败: %s", exc)
            abstracts = {}
        for record in missing:
            text = abstracts.get(record.get("pmid", ""))
            if text and len(text) >= min_length:
                record["abstract"] = text

    still_missing = [r for r in records if len(str(r.get("abstract", ""))) < min_length]
    for record in still_missing:
        doi = record.get("doi", "")
        if not doi:
            continue
        text = client.fetch_openalex_abstract(doi)
        if text and len(text) >= min_length:
            record["abstract"] = text
            continue
        text = client.fetch_s2_abstract(doi)
        if text and len(text) >= min_length:
            record["abstract"] = text

    filled = sum(1 for r in records if len(str(r.get("abstract", ""))) >= min_length)
    LOGGER.info("摘要富集完成，当前有摘要记录 %d/%d", filled, len(records))
    return filled


# ---------------------------------------------------------------- 引用格式化


def format_gbt7714(record: dict[str, Any]) -> str:
    """生成 GB/T 7714-2015 顺序编码制参考文献。

    格式：``主要责任者. 题名[J]. 刊名, 年, 卷(期): 页码. DOI.``
    """
    authors = [a.strip() for a in str(record.get("authors", "")).split(";") if a.strip()]
    if not authors:
        author_str = "佚名"
    elif len(authors) <= 3:
        author_str = ", ".join(authors)
    else:
        author_str = ", ".join(authors[:3]) + ", 等"

    title = (record.get("title") or "无题名").rstrip(".")
    journal = record.get("journal") or "未知刊名"
    year = record.get("year") or ""
    volume = record.get("volume") or ""
    issue = record.get("issue") or ""
    pages = record.get("pages") or ""
    doi = record.get("doi") or ""

    citation = f"{author_str}. {title}[J]. {journal}"
    if year:
        citation += f", {year}"
    if volume:
        citation += f", {volume}"
    if issue:
        citation += f"({issue})"
    if pages:
        citation += f": {pages}"
    citation += "."
    if doi:
        citation += f" DOI: {doi}."
    return citation


def format_ris(record: dict[str, Any], index: int = 1) -> str:
    """生成 RIS 格式（EndNote / Zotero / Mendeley 可直接导入）。"""
    authors = [a.strip() for a in str(record.get("authors", "")).split(";") if a.strip()]
    lines = ["TY  - JOUR"]
    for author in authors:
        parts = author.rsplit(" ", 1)
        if len(parts) == 2:
            lines.append(f"AU  - {parts[1]}, {parts[0]}")
        else:
            lines.append(f"AU  - {author}")
    lines.append(f"TI  - {record.get('title', '')}")
    if record.get("journal"):
        lines.append(f"JO  - {record['journal']}")
    if record.get("year"):
        lines.append(f"PY  - {record['year']}")
    if record.get("volume"):
        lines.append(f"VL  - {record['volume']}")
    if record.get("issue"):
        lines.append(f"IS  - {record['issue']}")
    if record.get("pages"):
        lines.append(f"SP  - {record['pages']}")
    if record.get("doi"):
        lines.append(f"DO  - {record['doi']}")
    if record.get("abstract"):
        lines.append(f"AB  - {record['abstract']}")
    if record.get("url"):
        lines.append(f"UR  - {record['url']}")
    lines.append(f"AN  - {index}")
    lines.append("ER  - ")
    return "\n".join(lines)


def _bibtex_key(record: dict[str, Any], index: int) -> str:
    """生成稳定的 BibTeX 引用键：首位作者姓 + 年份 + 标题首词。"""
    authors = [a.strip() for a in str(record.get("authors", "")).split(";") if a.strip()]
    surname = "anon"
    if authors:
        parts = authors[0].rsplit(" ", 1)
        surname = parts[-1] if parts else authors[0]
        surname = re.sub(r"[^A-Za-z]", "", surname) or "anon"
    year = record.get("year") or ""
    words = re.findall(r"[A-Za-z]+", str(record.get("title", "")))
    word = words[0].lower() if words else "ref"
    return f"{surname.lower()}{year}{word}{index}"


def format_bibtex(record: dict[str, Any], index: int = 1) -> str:
    """生成 BibTeX 条目（``@article{key, ...}``）。"""
    key = _bibtex_key(record, index)
    authors = [a.strip() for a in str(record.get("authors", "")).split(";") if a.strip()]
    fields = [
        ("author", " and ".join(authors) if authors else ""),
        ("title", record.get("title", "")),
        ("journal", record.get("journal", "")),
        ("year", str(record.get("year", ""))),
        ("volume", record.get("volume", "")),
        ("number", record.get("issue", "")),
        ("pages", record.get("pages", "")),
        ("doi", record.get("doi", "")),
        ("url", record.get("url", "")),
        ("abstract", record.get("abstract", "")),
    ]
    lines = [f"@article{{{key},"]
    for name, value in fields:
        if value:
            value = value.replace("{", "").replace("}", "")
            lines.append(f"  {name} = {{{value}}},")
    lines[-1] = lines[-1].rstrip(",")
    lines.append("}")
    return "\n".join(lines)


def format_plain(record: dict[str, Any], index: int = 1) -> str:
    """纯文本编号引用（中文综述常用）。"""
    return f"[{index}] {format_gbt7714(record)}"


def write_citations(
    records: Sequence[dict[str, Any]],
    fmt: str,
    path: str,
    encoding: str = "utf-8",
) -> str:
    """把记录写成指定引用格式，返回写出的文件路径。"""
    fmt = fmt.lower()
    if fmt == "ris":
        content = "\n\n".join(format_ris(r, i) for i, r in enumerate(records, 1)) + "\n"
    elif fmt == "bibtex":
        content = "\n\n".join(format_bibtex(r, i) for i, r in enumerate(records, 1)) + "\n"
    elif fmt in ("gbt7714", "gbt"):
        content = "\n".join(format_gbt7714(r) for r in records) + "\n"
    elif fmt in ("plain", "txt"):
        content = "\n".join(format_plain(r, i) for i, r in enumerate(records, 1)) + "\n"
    elif fmt == "json":
        content = json.dumps(list(records), ensure_ascii=False, indent=2) + "\n"
    else:
        raise SkillError(f"不支持的引用格式: {fmt}（可选 ris/bibtex/gbt7714/plain/json）")

    with open(path, "w", encoding=encoding, newline="\n") as handle:
        handle.write(content)
    return path


def parse_record_rows(rows: Iterable[dict[str, str]]) -> list[dict[str, Any]]:
    """把 CSV 行（字符串值）转成统一 schema 记录，供导出脚本复用。"""
    records: list[dict[str, Any]] = []
    for row in rows:
        record = make_record(
            title=row.get("title", ""),
            authors=row.get("authors", ""),
            year=row.get("year", ""),
            doi=row.get("doi", ""),
            journal=row.get("journal", ""),
            abstract=row.get("abstract", ""),
            citations=row.get("citations", 0),
            source=row.get("source", ""),
            pmid=row.get("pmid", ""),
            url=row.get("url", ""),
            doc_type=row.get("type", ""),
            oa_status=row.get("oa_status", ""),
            query=row.get("query", ""),
        )
        if row.get("volume"):
            record["volume"] = row["volume"]
        if row.get("issue"):
            record["issue"] = row["issue"]
        if row.get("pages"):
            record["pages"] = row["pages"]
        records.append(record)
    return records
