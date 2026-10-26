"""TEI XML 解析：把 GROBID 输出的 TEI 文档转成结构化数据。

TEI（Text Encoding Initiative）是 GROBID 的标准输出格式。命名空间为
``http://www.tei-c.org/ns/1.0``。

本模块提供
----------

- :func:`parse_tei` —— 解析成统一的 ``dict`` 结构
- :func:`extract_sections` —— 抽取章节树
- :func:`extract_keywords` —— 抽取关键词
- :func:`extract_first_sentences` —— 取正文各段首句（用于"引言精读"）
- :func:`is_valid_title` —— 判断 TEI 标题是否真的是论文标题
- :func:`to_markdown` / :func:`to_html` —— 输出转换

修掉的老仓库问题
-----------------
1. **关键词被空格切开** —— 老代码 ``keywords_text.split()`` 把
   ``"tumor necrosis factor"`` 拆成三个独立关键词。TEI 的关键词
   元素下通常有多个 ``<term>`` 子元素，应该按子元素取；退而求其次
   才按逗号/分号切分。
2. **标题误判** —— 老代码用"标题是否**包含** INTRODUCTION / METHODS /
   ABSTRACT"来判断标题是否有效。正常论文标题
   ``"Abstract Methods for Targeting X"`` 会被误杀。正确做法是判断
   标题**等于**章节名（GROBID 把章节名误填进 title 元素时，
   取值通常与某个 ``<head>`` 完全一致）。
3. **章节名未做空白归一** —— PDF 抽取常把 ``"INTRODUC TI ON"`` 这类
   带异常空格的词当标题，导致匹配失败。
4. **只取 ``abstract`` 下的第一个 ``<p>``** —— 多段摘要会丢内容。

本模块只依赖标准库 ``xml.etree.ElementTree``。
"""

from __future__ import annotations

import argparse
import html as html_lib
import json
import re
import sys
import unicodedata
import xml.etree.ElementTree as ET
from collections.abc import Iterable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional, Union

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "packages"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from biopipe_core import cli_parser, get_logger, run_cli  # noqa: E402

LOGGER = get_logger("biopipe.paper-fulltext-extract")


__all__ = [
    "TEI_NS",
    "Section",
    "ParsedPaper",
    "element_text",
    "extract_first_sentences",
    "extract_keywords",
    "extract_sections",
    "is_valid_title",
    "normalize_heading",
    "parse_tei",
    "to_html",
    "to_markdown",
]

TEI_NS = {"tei": "http://www.tei-c.org/ns/1.0"}

#: ``xml:id`` 属性在 ElementTree 里的真实键名。
#:
#: ``xml`` 是 XML 规范保留的前缀，ElementTree 会把它展开成
#: ``{http://www.w3.org/XML/1998/namespace}id``。
#: 直接 ``elem.get("xml:id")`` **永远取不到值**（返回 None 而不报错），
#: 这是个非常隐蔽的坑 —— 图表 id、参考文献 id 全部静默变成空字符串。
XML_ID_ATTR = "{http://www.w3.org/XML/1998/namespace}id"

# GROBID 偶尔把章节名误填进 title 元素。这些是**章节名**，
# 用来和 <head> 做精确比对，而不是用"包含"来判断。
_SECTION_HEADINGS = frozenset({
    "abstract", "introduction", "background", "results", "discussion",
    "conclusion", "conclusions", "methods", "materials and methods",
    "experimental procedures", "experimental procedure", "acknowledgments",
    "acknowledgements", "references", "supplementary information",
    "data availability", "author contributions", "competing interests",
    "ethics statement", "supplementary material",
})

#: 常见学术章节名（归一化后的小写形式）。
#:
#: 这份字典有两个用途：
#:
#: 1. **断字修复** —— PDF 抽取会把 ``"INTRODUC TI ON"`` 这种
#:    被空格劈开的词送过来。判断它是不是断字，靠的是
#:    "去掉空格后是不是一个已知章节名"，而不是猜正则。
#: 2. **章节名归一** —— ``"Materials  and   Methods"`` 直接命中。
#:
#: 字典比 :data:`_SECTION_HEADINGS` 更宽，因为断字可能出现在
#: 任意章节名上；而标题有效性判断要更保守（宁可漏判也不误杀）。
_KNOWN_HEADINGS = frozenset(
    _SECTION_HEADINGS
    | {
        "abstract and introduction",
        "acknowledgement", "acknowledgment",
        "appendix", "appendices",
        "author contributions", "authors contributions",
        "availability of data and materials",
        "bibliography",
        "competing interests", "conflict of interest", "conflicts of interest",
        "conclusion and future work", "concluding remarks",
        "contributions", "data and code availability",
        "declarations", "ethics approval", "ethics statement",
        "experimental section",
        "future work", "introduction and background",
        "materials and method", "method", "methodology",
        "nomenclature", "outlook",
        "related work", "related works",
        "results and discussion",
        "supplementary material", "supplementary materials",
        "table of contents",
    }
)

_WS_RE = re.compile(r"\s+")
# PDF 抽取常见的断字（字间距被误当空格）
_BROKEN_WORD_RE = re.compile(r"\b([A-Z])\s+([A-Z]{2,})\b")


def normalize_heading(text: str) -> str:
    """归一化章节标题：压缩空白、修掉 PDF 断字、统一大小写比较。

    PDF 抽取会把 ``"INTRODUC TI ON"`` 这样的词原样送过来 ——
    断字发生在渲染层，成因是字间距被当成了空格。

    断字修复靠**字典**而不是猜正则：把标题里的空格全部去掉，
    看是否正好等于某个已知章节名。``"INTRODUC TI ON"`` 去空格后是
    ``"introduction"`` —— 命中，于是还原成 ``"introduction"``。
    这种做法对 ``"DISCUSS I ON"``、``"MATER I ALS AND METHODS"``
    都成立，而 ``"Materials  and   Methods"``（本来就是多个词）
    去空格后是 ``"materialsandmethods"``，不命中任何单词，
    保持原样只压缩空白。

    例::

        >>> normalize_heading("INTRODUC TI ON")
        'introduction'
        >>> normalize_heading("DISCUSS I ON")
        'discussion'
        >>> normalize_heading("Materials  and   Methods")
        'materials and methods'
    """
    if not text:
        return ""
    text = _WS_RE.sub(" ", text).strip()
    lowered = text.lower()

    # 1) 原样命中已知章节名
    if lowered in _KNOWN_HEADINGS:
        return lowered

    # 2) 断字：去空格后命中已知章节名
    collapsed = lowered.replace(" ", "")
    if collapsed != lowered:
        for canonical in _KNOWN_HEADINGS:
            if canonical.replace(" ", "") == collapsed:
                return canonical

    # 3) 通用断字启发式：单字母 + 大写块（如 "A BCD"）
    text = _BROKEN_WORD_RE.sub(r"\1\2", text)
    text = _WS_RE.sub(" ", text).strip()
    return text.lower()


def _element_text(element: Optional[ET.Element]) -> str:
    """取元素的全部文本（含嵌套元素的文本）。"""
    if element is None:
        return ""
    parts: list[str] = []
    if element.text:
        parts.append(element.text)
    for child in element:
        parts.append(_element_text(child))
        if child.tail:
            parts.append(child.tail)
    return _WS_RE.sub(" ", "".join(parts)).strip()


#: 公开别名 —— 同包内的其他模块需要这个辅助函数
element_text = _element_text


def _iter_text(element: Optional[ET.Element]) -> Iterable[str]:
    """递归产出所有非空文本片段（保持段落级粒度）。"""
    if element is None:
        return
    if element.text and element.text.strip():
        yield element.text.strip()
    for child in element:
        yield from _iter_text(child)
        if child.tail and child.tail.strip():
            yield child.tail.strip()


@dataclass
class Section:
    """一个章节。"""

    title: str
    paragraphs: list[str] = field(default_factory=list)
    level: int = 0
    children: list[Section] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "title": self.title,
            "level": self.level,
            "paragraphs": self.paragraphs,
            "children": [child.to_dict() for child in self.children],
        }


@dataclass
class ParsedPaper:
    """一篇文献的结构化表示。"""

    title: str = ""
    authors: list[str] = field(default_factory=list)
    affiliations: list[str] = field(default_factory=list)
    abstract: str = ""
    keywords: list[str] = field(default_factory=list)
    doi: str = ""
    year: str = ""
    journal: str = ""
    sections: list[Section] = field(default_factory=list)
    references: list[dict[str, str]] = field(default_factory=list)
    figures: list[dict[str, Any]] = field(default_factory=list)
    source_file: str = ""

    def to_dict(self) -> dict[str, Any]:
        return {
            "title": self.title,
            "authors": self.authors,
            "affiliations": self.affiliations,
            "abstract": self.abstract,
            "keywords": self.keywords,
            "doi": self.doi,
            "year": self.year,
            "journal": self.journal,
            "source_file": self.source_file,
            "sections": [section.to_dict() for section in self.sections],
            "references": self.references,
            "figures": self.figures,
        }

    def iter_sections(self) -> Iterable[Section]:
        """深度优先遍历所有章节。"""
        stack = list(reversed(self.sections))
        while stack:
            node = stack.pop()
            yield node
            stack.extend(reversed(node.children))

    def full_text(self) -> str:
        """拼接全部章节正文。"""
        chunks: list[str] = []
        for section in self.iter_sections():
            chunks.extend(section.paragraphs)
        return "\n\n".join(chunks)


# ---------------------------------------------------------------- 标题


def _collect_headings(root: ET.Element) -> set[str]:
    """收集文档里出现过的所有章节名（归一化后）。"""
    headings: set[str] = set()
    for head in root.iter(f"{{{TEI_NS['tei']}}}head"):
        text = _element_text(head)
        if text:
            headings.add(normalize_heading(text))
    return headings


def is_valid_title(title: str, headings: Optional[set[str]] = None) -> bool:
    """判断一个字符串是否**真的**是论文标题。

    背景：GROBID 在某些 PDF 上会把章节名误填进 ``<title type="main">``，
    比如标题变成 ``"INTRODUCTION"``。老仓库用"是否包含 INTRODUCTION /
    METHODS / ABSTRACT"来判断 —— 这会误杀正常标题，例如
    ``"Abstract Methods for Targeting X"``。

    正确判据：标题**归一化后等于**某个章节名，才认为它是被误填的。

    参数
    ----
    title
        待判断的标题。
    headings
        文档中出现过的章节名集合（归一化后）。传入时启用精确比对；
        不传时只做基础的长度与形状检查，以及与**已知章节名词典**
        的比对 —— 后者是必须的：``is_valid_title("INTRODUCTION")``
        在不传 ``headings`` 时同样要判为无效，否则调用方一旦漏传
        就会把章节名当成论文标题。
    """
    if not title or not title.strip():
        return False

    normalized = normalize_heading(title)

    # 过短 / 明显不是标题
    if len(title.strip()) < 10:
        return False
    # 含竖线且整体很短 —— 多半是 PDF 侧栏文字被拼进来了
    if "|" in title and len(title) < 20:
        return False
    # 全部是数字或标点
    if not re.search(r"[A-Za-z一-鿿]", title):
        return False

    # 与已知章节名完全相同 —— 无论有没有传 headings 都要判无效。
    # 注意这里是**精确相等**而不是"包含"：
    # "Abstract Methods for Targeting X" 包含 "abstract" 但不是章节名，
    # 用包含匹配会把它误杀（老仓库就是这么错的）。
    if normalized in _SECTION_HEADINGS:
        return False

    if headings:
        # 归一化后与文档中真实出现的章节名完全相同
        if normalized in headings:
            return False

    return True


# ---------------------------------------------------------------- 抽取


def extract_keywords(root: ET.Element) -> list[str]:
    """抽取关键词。

    TEI 里关键词有两种常见形态：

    1. ``<keywords><term>keyword1</term><term>keyword2</term></keywords>``
       —— 每个 ``<term>`` 一个关键词（GROBID 的标准输出）；
    2. ``<keywords>keyword1, keyword2</keywords>``
       —— 纯文本，用逗号或分号分隔。

    老仓库一律用 ``.split()`` 按**空格**切分，
    ``"tumor necrosis factor"`` 会被拆成三个词。
    这里优先按 ``<term>`` 子元素取，其次按逗号/分号切分，
    最后才按空格切分（并去重、保留原顺序）。
    """
    for keywords_elem in root.iter(f"{{{TEI_NS['tei']}}}keywords"):
        terms: list[str] = []

        # 形态 1：<term> 子元素
        for term in keywords_elem.findall(f".//{{{TEI_NS['tei']}}}term"):
            text = _element_text(term)
            if text:
                terms.append(text)

        if terms:
            return _dedupe(terms)

        # 形态 2：纯文本，用逗号/分号分隔
        text = _element_text(keywords_elem)
        if not text:
            continue
        parts = re.split(r"[;,，；、]", text)
        cleaned = [p.strip() for p in parts if p.strip()]
        if len(cleaned) > 1:
            return _dedupe(cleaned)

        # 形态 3：空格分隔（最后兜底）
        return _dedupe([p for p in text.split() if p])

    return []


def _dedupe(items: list[str]) -> list[str]:
    """去重但保留原顺序。"""
    seen: set[str] = set()
    out: list[str] = []
    for item in items:
        key = item.strip().lower()
        if not key or key in seen:
            continue
        seen.add(key)
        out.append(item.strip())
    return out


def extract_sections(root: ET.Element) -> list[Section]:
    """抽取正文章节树。

    GROBID 的正文结构是 ``<body><div><head>...</head><p>...</p></div></body>``，
    ``div`` 可以嵌套。返回按文档顺序排列的顶层 ``Section`` 列表。

    少数 TEI 变体不使用 ``<div>``，而是把章节平铺在 ``<body>`` 下
    （``<head>`` 与 ``<p>`` 互为兄弟）。这种情况有专门的兜底分支。
    """
    body = root.find(f".//{{{TEI_NS['tei']}}}body")
    if body is None:
        return []

    sections = _parse_divs(body, level=0)
    if sections:
        return sections

    # 兜底：没有 <div>，章节直接平铺在 <body> 下
    flat: list[Section] = []
    current: Optional[Section] = None
    for child in body:
        tag = child.tag.rsplit("}", 1)[-1]
        if tag == "head":
            current = Section(title=_element_text(child), level=0)
            flat.append(current)
        elif tag == "p":
            text = _element_text(child)
            if text:
                if current is None:
                    current = Section(title="", level=0)
                    flat.append(current)
                current.paragraphs.append(text)
    return flat


def _parse_divs(parent: ET.Element, level: int) -> list[Section]:
    """递归解析 ``<div>`` 层级。

    注意：末尾的兜底分支只在**顶层** ``<body>`` 上生效。有些 TEI 变体
    不使用 ``<div>``，而是把章节直接平铺在 ``<body>`` 下（``<head>``
    + ``<p>`` 兄弟）。如果在任意层级都做这个兜底，会把每个 div 自己的
    ``<head>`` 重复解析一次 —— 表现为章节树里每个标题都出现两遍
    （一次在本层，一次在下一层的空兜底里）。
    """
    sections: list[Section] = []

    for div in parent.findall(f"{{{TEI_NS['tei']}}}div"):
        head = div.find(f"{{{TEI_NS['tei']}}}head")
        title = _element_text(head) if head is not None else ""

        paragraphs: list[str] = []
        for p in div.findall(f"{{{TEI_NS['tei']}}}p"):
            text = _element_text(p)
            if text:
                paragraphs.append(text)
        # 图表引用段落（如 <figElem>）也当作正文的一部分
        for elem in div.findall(f"{{{TEI_NS['tei']}}}figElem"):
            text = _element_text(elem)
            if text:
                paragraphs.append(text)

        section = Section(title=title, paragraphs=paragraphs, level=level)
        section.children = _parse_divs(div, level + 1)
        sections.append(section)

    return sections


def extract_first_sentences(
    root: ET.Element,
    section_keyword: str = "introduction",
    max_sentences: int = 3,
) -> list[str]:
    """取出指定章节（默认 Introduction）各段落的**首句**。

    用途：快速读一批文献的引言 —— 每个引言段落的第一句通常就点明了
    这篇文章要解决什么问题，比读完整段更高效。

    匹配是**归一化后精确相等**，所以 ``"Introduction"``、
    ``"INTRODUCTION"``、``"INTRODUC TI ON"`` 都能匹配上。
    """
    target = normalize_heading(section_keyword)
    found: list[str] = []

    for section in _iter_section_nodes(root):
        if normalize_heading(section.title) != target:
            continue
        for paragraph in section.paragraphs:
            sentence = _first_sentence(paragraph)
            if sentence:
                found.append(sentence)
            if len(found) >= max_sentences:
                return found
    return found


def _iter_section_nodes(root: ET.Element) -> Iterable[Section]:
    """遍历文档里的所有 ``<div>``（含嵌套）。"""
    for div in root.iter(f"{{{TEI_NS['tei']}}}div"):
        head = div.find(f"{{{TEI_NS['tei']}}}head")
        title = _element_text(head) if head is not None else ""
        paragraphs = [
            text
            for text in (
                _element_text(p)
                for p in div.findall(f"{{{TEI_NS['tei']}}}p")
            )
            if text
        ]
        yield Section(title=title, paragraphs=paragraphs)


_ABBREV_RE = re.compile(
    r"\b(?:e\.g|i\.e|et al|vs|cf|Fig|Fig\.|Eq|Eq\.|No|No\.|ref|refs|approx)\.$",
    re.IGNORECASE,
)


def _first_sentence(text: str) -> str:
    """取段落首句。

    不能简单按 ``"."`` 切 —— ``"e.g."``、``"et al."``、``"Fig. 3"``
    都会造成误切。这里用一个小启发式：找到句末标点后，若下一个字符
    是大写或空白，才认为是句子边界。
    """
    text = text.strip()
    if not text:
        return ""

    for match in re.finditer(r"[.!?](?=\s|$)", text):
        end = match.end()
        candidate = text[:end]
        # 排除缩写结尾
        if _ABBREV_RE.search(candidate):
            continue
        # 排除 "3. " 这样的编号后紧跟句号的情况
        if re.search(r"\d\.$", candidate) and end < len(text) and text[end:end + 1].isdigit():
            continue
        return candidate.strip()

    return text


def _extract_authors(root: ET.Element) -> list[str]:
    """抽取作者列表。"""
    authors: list[str] = []
    for pers_name in root.iter(f"{{{TEI_NS['tei']}}}persName"):
        forenames = [
            _element_text(node)
            for node in pers_name.findall(f"{{{TEI_NS['tei']}}}forename")
        ]
        surname = _element_text(pers_name.find(f"{{{TEI_NS['tei']}}}surname"))
        parts = [p for p in forenames + ([surname] if surname else []) if p]
        if parts:
            authors.append(" ".join(parts))
    return _dedupe(authors)


def _extract_affiliations(root: ET.Element) -> list[str]:
    """抽取机构列表。"""
    out: list[str] = []
    for aff in root.iter(f"{{{TEI_NS['tei']}}}affiliation"):
        text = _element_text(aff)
        if text:
            out.append(text)
    return _dedupe(out)


def _extract_idno(root: ET.Element, kind: str) -> str:
    """抽取指定类型的标识符（DOI / PMID 等）。"""
    for idno in root.iter(f"{{{TEI_NS['tei']}}}idno"):
        if idno.get("type") == kind:
            return _element_text(idno)
    return ""


def _extract_biblio(root: ET.Element, field_name: str) -> str:
    """从 ``<biblScope>`` 里取某个字段。"""
    for scope in root.iter(f"{{{TEI_NS['tei']}}}biblScope"):
        if scope.get("unit") == field_name:
            return _element_text(scope)
    return ""


def _extract_references(root: ET.Element) -> list[dict[str, str]]:
    """抽取参考文献列表。

    **范围很重要**：只取 ``<back>`` / ``<listBibl>`` 下的
    ``<biblStruct>``。

    ``<sourceDesc>`` 里也有 ``<biblStruct>``，但那是**当前这篇文献
    自己的元数据**（期刊名、DOI），不是参考文献。遍历全文会把
    它也数进去，导致参考文献数量多 1 —— 而且排第一。
    """
    refs: list[dict[str, str]] = []

    back = root.find(f".//{{{TEI_NS['tei']}}}back")
    if back is None:
        return refs

    for bibl in back.iter(f"{{{TEI_NS['tei']}}}biblStruct"):
        entry = {
            "id": bibl.get(XML_ID_ATTR, "") or "",
            "raw": _element_text(bibl),
        }
        title_elem = bibl.find(f".//{{{TEI_NS['tei']}}}title")
        if title_elem is not None:
            entry["title"] = _element_text(title_elem)
        doi = _extract_idno(bibl, "DOI")
        if doi:
            entry["doi"] = doi
        refs.append(entry)
    return refs


def _extract_figures(root: ET.Element) -> list[dict[str, Any]]:
    """抽取图表信息。"""
    figures: list[dict[str, Any]] = []
    for fig in root.iter(f"{{{TEI_NS['tei']}}}figure"):
        xml_id = fig.get(XML_ID_ATTR, "") or ""
        label = _element_text(fig.find(f"{{{TEI_NS['tei']}}}label"))
        caption = ""
        caption_elem = fig.find(f"{{{TEI_NS['tei']}}}figDesc")
        if caption_elem is not None:
            caption = _element_text(caption_elem)
        graphics: list[str] = []
        for graphic in fig.iter(f"{{{TEI_NS['tei']}}}graphic"):
            url = graphic.get("url", "")
            if url:
                graphics.append(url)
        figures.append({
            "id": xml_id,
            "label": label,
            "caption": caption,
            "graphics": graphics,
        })
    return figures


# ---------------------------------------------------------------- 主解析


def _candidate_titles(root: ET.Element) -> list[str]:
    """列出候选的论文标题。

    TEI 里 ``<title>`` 元素有好几处：``<titleStmt><title type="main">``
    是论文标题；``<monogr><title>`` 是**期刊名**；``<biblStruct><analytic>
    <title>`` 是文献条目的标题；``<figure><head>`` 之类更是无关。

    候选只从 ``<titleStmt>`` 下取 —— 这样才不会被期刊名或参考文献标题
    污染。老仓库的回退逻辑是"遍历所有 title 元素"，结果在标题被
    GROBID 误填成章节名时，会回退到 **期刊名**，比回退到文件名更糟。
    """
    candidates: list[str] = []
    for title_stmt in root.iter(f"{{{TEI_NS['tei']}}}titleStmt"):
        for title_elem in title_stmt.findall(f"{{{TEI_NS['tei']}}}title"):
            text = _element_text(title_elem)
            if text:
                candidates.append(text)
    return candidates


def parse_tei(
    xml_path: Union[str, Path],
    fallback_title: str = "",
) -> ParsedPaper:
    """解析一个 TEI XML 文件。

    参数
    ----
    xml_path
        TEI XML 文件路径。
    fallback_title
        当 TEI 里的标题无效时使用的兜底标题（通常传文件名）。

    返回
    ----
    :class:`ParsedPaper`。
    """
    path = Path(xml_path)
    tree = ET.parse(path)
    root = tree.getroot()

    headings = _collect_headings(root)

    title = ""
    title_elem = root.find(f'.//{{{TEI_NS["tei"]}}}title[@type="main"]')
    if title_elem is not None:
        title = _element_text(title_elem)

    if not is_valid_title(title, headings):
        # 退化到 titleStmt 下的其他 title（**不**扫描全文，避免拿到期刊名）
        for candidate in _candidate_titles(root):
            if is_valid_title(candidate, headings):
                title = candidate
                break

    if not is_valid_title(title, headings):
        title = fallback_title or path.stem

    # 摘要：拼接所有段落，而不是只取第一段
    abstract_parts: list[str] = []
    for abstract_elem in root.iter(f"{{{TEI_NS['tei']}}}abstract"):
        for p in abstract_elem.findall(f".//{{{TEI_NS['tei']}}}p"):
            text = _element_text(p)
            if text:
                abstract_parts.append(text)
        if not abstract_parts:
            text = _element_text(abstract_elem)
            if text:
                abstract_parts.append(text)

    journal = ""
    for monogr in root.iter(f"{{{TEI_NS['tei']}}}monogr"):
        title_elem = monogr.find(f"{{{TEI_NS['tei']}}}title")
        if title_elem is not None:
            journal = _element_text(title_elem)
            break

    return ParsedPaper(
        title=title,
        authors=_extract_authors(root),
        affiliations=_extract_affiliations(root),
        abstract="\n\n".join(abstract_parts),
        keywords=extract_keywords(root),
        doi=_extract_idno(root, "DOI"),
        year=_extract_biblio(root, "year"),
        journal=journal,
        sections=extract_sections(root),
        references=_extract_references(root),
        figures=_extract_figures(root),
        source_file=str(path),
    )


# ---------------------------------------------------------------- 输出


def to_markdown(paper: ParsedPaper) -> str:
    """把解析结果转成 Markdown。"""
    lines: list[str] = []

    lines.append(f"# {paper.title}")
    lines.append("")

    if paper.authors:
        lines.append(f"**作者**: {', '.join(paper.authors)}")
        lines.append("")
    meta_parts = []
    if paper.journal:
        meta_parts.append(f"**期刊**: {paper.journal}")
    if paper.year:
        meta_parts.append(f"**年份**: {paper.year}")
    if paper.doi:
        meta_parts.append(f"**DOI**: {paper.doi}")
    if meta_parts:
        lines.append(" · ".join(meta_parts))
        lines.append("")
    if paper.keywords:
        lines.append(f"**关键词**: {', '.join(paper.keywords)}")
        lines.append("")

    if paper.abstract:
        lines.append("## 摘要")
        lines.append("")
        lines.append(paper.abstract)
        lines.append("")

    def emit(sections: list[Section]) -> None:
        for section in sections:
            if section.title:
                level = min(2 + section.level, 6)
                lines.append(f"{'#' * level} {section.title}")
                lines.append("")
            for paragraph in section.paragraphs:
                lines.append(paragraph)
                lines.append("")
            emit(section.children)

    lines.append("## 正文")
    lines.append("")
    emit(paper.sections)

    if paper.references:
        lines.append("## 参考文献")
        lines.append("")
        for index, ref in enumerate(paper.references, start=1):
            raw = ref.get("raw") or ref.get("title") or ""
            lines.append(f"{index}. {raw}")
        lines.append("")

    return "\n".join(lines).rstrip() + "\n"


def to_html(paper: ParsedPaper, image_base: str = "") -> str:
    """把解析结果转成 HTML。

    参数
    ----
    image_base
        图像目录的 URL 前缀。TEI 里图形是相对路径（如 ``figures/fig_1.png``），
        渲染时需要拼上前缀才能显示。
    """
    parts: list[str] = []
    parts.append("<!DOCTYPE html>")
    parts.append('<html lang="zh">')
    parts.append("<head>")
    parts.append('<meta charset="utf-8">')
    parts.append(f"<title>{html_lib.escape(paper.title)}</title>")
    parts.append("<style>")
    parts.append(
        "body{font-family:system-ui,-apple-system,'Segoe UI',sans-serif;"
        "max-width:900px;margin:2rem auto;padding:0 1rem;line-height:1.7;color:#222}"
    )
    parts.append("h1{font-size:1.7rem;line-height:1.3}")
    parts.append("h2{font-size:1.35rem;margin-top:2rem;border-bottom:1px solid #eee;padding-bottom:.3rem}")
    parts.append("h3{font-size:1.15rem;margin-top:1.5rem}")
    parts.append(".meta{color:#666;font-size:.92rem}")
    parts.append(".keywords{background:#f6f8fa;padding:.6rem .9rem;border-radius:6px;font-size:.92rem}")
    parts.append("figure{margin:1.5rem 0;text-align:center}")
    parts.append("figure img{max-width:100%;border:1px solid #e5e5e5;border-radius:4px}")
    parts.append("figcaption{font-size:.88rem;color:#555;margin-top:.5rem;text-align:left}")
    parts.append("</style>")
    parts.append("</head>")
    parts.append("<body>")

    parts.append(f"<h1>{html_lib.escape(paper.title)}</h1>")

    if paper.authors:
        authors = ", ".join(html_lib.escape(a) for a in paper.authors)
        parts.append(f'<p class="meta">{authors}</p>')

    meta_bits = []
    if paper.journal:
        meta_bits.append(f"期刊：{html_lib.escape(paper.journal)}")
    if paper.year:
        meta_bits.append(f"年份：{html_lib.escape(paper.year)}")
    if paper.doi:
        doi_link = f"https://doi.org/{html_lib.escape(paper.doi)}"
        meta_bits.append(f'DOI：<a href="{doi_link}">{html_lib.escape(paper.doi)}</a>')
    if meta_bits:
        parts.append(f'<p class="meta">{" · ".join(meta_bits)}</p>')

    if paper.keywords:
        kws = ", ".join(html_lib.escape(k) for k in paper.keywords)
        parts.append(f'<p class="keywords"><strong>关键词</strong>：{kws}</p>')

    if paper.abstract:
        parts.append("<h2>摘要</h2>")
        for para in paper.abstract.split("\n\n"):
            parts.append(f"<p>{html_lib.escape(para)}</p>")

    figure_by_id = {f["id"]: f for f in paper.figures if f.get("id")}

    def emit(sections: list[Section]) -> None:
        for section in sections:
            level = min(2 + section.level, 6)
            if section.title:
                parts.append(f"<h{level}>{html_lib.escape(section.title)}</h{level}>")
            for paragraph in section.paragraphs:
                parts.append(f"<p>{html_lib.escape(paragraph)}</p>")
                # 段落里的图引用（如 <figRef>）—— TEI 已在 p 里内联，
                # 这里只处理独立成段的图表说明
                for fig in figure_by_id.values():
                    if fig.get("label") and fig["label"] in paragraph:
                        pass
            emit(section.children)

    parts.append("<h2>正文</h2>")
    emit(paper.sections)

    if paper.figures:
        parts.append("<h2>图表</h2>")
        for fig in paper.figures:
            label = html_lib.escape(fig.get("label") or "")
            caption = html_lib.escape(fig.get("caption") or "")
            parts.append("<figure>")
            for graphic in fig.get("graphics", []):
                src = f"{image_base.rstrip('/')}/{graphic.lstrip('/')}" if image_base else graphic
                parts.append(f'<img src="{html_lib.escape(src)}" alt="{label}">')
            if label or caption:
                parts.append(f"<figcaption>{label} {caption}</figcaption>")
            parts.append("</figure>")

    if paper.references:
        parts.append("<h2>参考文献</h2>")
        parts.append("<ol>")
        for ref in paper.references:
            raw = html_lib.escape(ref.get("raw") or ref.get("title") or "")
            parts.append(f"<li>{raw}</li>")
        parts.append("</ol>")

    parts.append("</body>")
    parts.append("</html>")
    return "\n".join(parts) + "\n"


def _build_parser() -> argparse.ArgumentParser:
    """构造命令行解析器。"""
    parser = cli_parser(
        prog="tei-parser",
        description="解析 TEI XML（GROBID 输出）为 JSON / Markdown / HTML",
        skill="paper-fulltext-extract",
    )
    parser.add_argument("--input", required=True, metavar="PATH", help="TEI XML 文件")
    parser.add_argument("--json", metavar="PATH", help="写出 JSON")
    parser.add_argument("--markdown", metavar="PATH", help="写出 Markdown")
    parser.add_argument("--html", metavar="PATH", help="写出 HTML")
    parser.add_argument("--image-base", default="", help="HTML 中图像的 URL 前缀")
    return parser


def _main(args: argparse.Namespace) -> int:
    """命令行入口：把 TEI XML 解析成 JSON / Markdown / HTML。"""
    paper = parse_tei(args.input)

    LOGGER.info("标题: %s", paper.title)
    LOGGER.info(
        "作者 %d 位，章节 %d 个，关键词 %d 个，参考文献 %d 条",
        len(paper.authors), len(list(paper.iter_sections())),
        len(paper.keywords), len(paper.references),
    )

    if args.json:
        Path(args.json).write_text(
            json.dumps(paper.to_dict(), indent=2, ensure_ascii=False), encoding="utf-8"
        )
        LOGGER.info("已写出 %s", args.json)
    if args.markdown:
        Path(args.markdown).write_text(to_markdown(paper), encoding="utf-8")
        LOGGER.info("已写出 %s", args.markdown)
    if args.html:
        Path(args.html).write_text(
            to_html(paper, image_base=args.image_base), encoding="utf-8"
        )
        LOGGER.info("已写出 %s", args.html)

    if not (args.json or args.markdown or args.html):
        # 无输出目标时，把 JSON 打到 stdout —— 供管道消费
        print(json.dumps(paper.to_dict(), indent=2, ensure_ascii=False))

    return 0


if __name__ == "__main__":
    sys.exit(run_cli(_main, parser=_build_parser(), skill="paper-fulltext-extract", version="1.0.0"))
