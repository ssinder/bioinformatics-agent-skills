"""TEI 文献语料的全文检索。

为一批 TEI XML 文献建立索引，然后按**章节标题**、**正文关键词**、
**论文标题**检索，并导出 CSV 元数据表。

从老仓库 ``paper_extract/scripts/xml_literature_searcher.py`` 改造而来，
保留其核心能力（索引 / 三种检索 / 上下文 / CSV 导出），并修掉：

1. **关键词被空格切开** —— 见 :func:`tei_parser.extract_keywords`。
2. **标题误判** —— 见 :func:`tei_parser.is_valid_title`。
3. **章节名未归一** —— 检索 ``"introduction"`` 匹配不到
   ``"INTRODUC TI ON"``，现在用归一化后的名字比对。
4. **单个文件解析失败会中断整批** —— 现在跳过并记录。

本模块只依赖标准库。
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
import xml.etree.ElementTree as ET
from collections.abc import Iterable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional, Union

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "packages"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from biopipe_core import cli_parser, get_logger, run_cli  # noqa: E402
from tei_parser import (  # noqa: E402
    TEI_NS,
    element_text,
    extract_keywords,
    is_valid_title,
    normalize_heading,
)

LOGGER = get_logger("biopipe.paper-fulltext-extract")


__all__ = [
    "CorpusIndex",
    "SearchHit",
    "build_index",
]


@dataclass
class SearchHit:
    """一条检索命中。"""

    file: str
    title: str
    kind: str = "paragraph"
    section_title: str = ""
    paragraph_index: int = 0
    match_text: str = ""
    level: int = 0

    def to_dict(self) -> dict[str, Any]:
        return {
            "file": self.file,
            "title": self.title,
            "kind": self.kind,
            "section_title": self.section_title,
            "paragraph_index": self.paragraph_index,
            "match_text": self.match_text,
            "level": self.level,
        }


@dataclass
class CorpusIndex:
    """一批 TEI 文献的索引。

    结构::

        {
          "files": {
            "<文件名>": {
              "title": ..., "abstract": ..., "keywords": [...],
              "sections": [{"title":..., "level":..., "paragraphs":[...]}]
            }
          }
        }
    """

    index: dict[str, Any] = field(default_factory=lambda: {"files": {}})
    errors: list[str] = field(default_factory=list)

    # ------------------------------------------------------------ 构建

    def build(
        self,
        xml_dir: Union[str, Path],
        pattern: str = "*.tei.xml",
        select_titles: Optional[Iterable[str]] = None,
    ) -> CorpusIndex:
        """扫描目录并建立索引。"""
        directory = Path(xml_dir)
        if not directory.is_dir():
            raise FileNotFoundError(f"目录不存在: {directory}")

        selected = {t.strip() for t in select_titles or [] if t.strip()}
        files = sorted(directory.glob(pattern))
        # 也接受裸 .xml（不是 .tei.xml）
        files.extend(sorted(p for p in directory.glob("*.xml") if p not in files))

        self.index = {"files": {}}
        self.errors = []

        for xml_file in files:
            try:
                self._index_file(xml_file)
            except ET.ParseError as exc:
                self.errors.append(f"{xml_file.name}: XML 解析失败 — {exc}")
                continue
            except Exception as exc:  # noqa: BLE001
                self.errors.append(f"{xml_file.name}: {type(exc).__name__} — {exc}")
                continue

            if selected:
                entry = self.index["files"].get(xml_file.name)
                if entry and entry["title"] not in selected:
                    del self.index["files"][xml_file.name]

        return self

    def _index_file(self, xml_file: Path) -> None:
        """索引单个 TEI 文件。"""
        root = ET.parse(xml_file).getroot()

        # 收集章节名，用于标题有效性判断
        headings = {
            normalize_heading(element_text(head))
            for head in root.iter(f"{{{TEI_NS['tei']}}}head")
        }

        title = ""
        title_elem = root.find(f'.//{{{TEI_NS["tei"]}}}title[@type="main"]')
        if title_elem is not None:
            title = element_text(title_elem)
        if not is_valid_title(title, headings):
            # 只在 <titleStmt> 下找候选 —— 全局扫描会命中 <monogr><title>
            # （期刊名），回退到期刊名比回退到文件名更糟
            for title_stmt in root.iter(f"{{{TEI_NS['tei']}}}titleStmt"):
                for candidate in title_stmt.findall(f"{{{TEI_NS['tei']}}}title"):
                    text = element_text(candidate)
                    if is_valid_title(text, headings):
                        title = text
                        break
                if is_valid_title(title, headings):
                    break
        if not is_valid_title(title, headings):
            title = xml_file.stem.replace(".pdf.tei", "").replace(".tei", "")

        abstract_parts: list[str] = []
        for abstract_elem in root.iter(f"{{{TEI_NS['tei']}}}abstract"):
            for p in abstract_elem.findall(f".//{{{TEI_NS['tei']}}}p"):
                text = element_text(p)
                if text:
                    abstract_parts.append(text)

        entry = {
            "title": title,
            "abstract": "\n\n".join(abstract_parts),
            "keywords": extract_keywords(root),
            "sections": [],
        }
        entry["sections"] = _collect_sections(root)
        self.index["files"][xml_file.name] = entry

    # ------------------------------------------------------------ 检索

    def search_sections(
        self, keyword: str, case_sensitive: bool = False
    ) -> list[SearchHit]:
        """按章节标题检索。"""
        hits: list[SearchHit] = []
        needle = keyword if case_sensitive else keyword.lower()
        for file_name, info in self.index["files"].items():
            for section in info["sections"]:
                title = section["title"]
                haystack = title if case_sensitive else title.lower()
                if needle in haystack:
                    hits.append(SearchHit(
                        file=file_name,
                        title=info["title"],
                        kind="section",
                        section_title=title,
                        level=section.get("level", 0),
                    ))
        return hits

    def search_paragraphs(
        self, keyword: str, case_sensitive: bool = False, context: int = 100
    ) -> list[SearchHit]:
        """按正文关键词检索（返回上下文片段）。"""
        hits: list[SearchHit] = []
        needle = keyword if case_sensitive else keyword.lower()
        for file_name, info in self.index["files"].items():
            for section in info["sections"]:
                for index, paragraph in enumerate(section.get("paragraphs", [])):
                    haystack = paragraph if case_sensitive else paragraph.lower()
                    if needle in haystack:
                        hits.append(SearchHit(
                            file=file_name,
                            title=info["title"],
                            kind="paragraph",
                            section_title=section["title"],
                            paragraph_index=index,
                            match_text=_context(paragraph, keyword, context, case_sensitive),
                        ))
        return hits

    def search_titles(
        self, keyword: str, case_sensitive: bool = False
    ) -> list[SearchHit]:
        """按论文标题检索。"""
        hits: list[SearchHit] = []
        needle = keyword if case_sensitive else keyword.lower()
        for file_name, info in self.index["files"].items():
            title = info["title"]
            haystack = title if case_sensitive else title.lower()
            if needle in haystack:
                hits.append(SearchHit(
                    file=file_name, title=title, kind="title"
                ))
        return hits

    def search(
        self,
        keyword: str,
        scope: str = "paragraph",
        case_sensitive: bool = False,
        context: int = 100,
    ) -> list[SearchHit]:
        """统一检索入口。

        scope 取值：``paragraph``（正文）/ ``section``（章节标题）/
        ``title``（论文标题）/ ``all``（三者都要）。
        """
        if scope == "paragraph":
            return self.search_paragraphs(keyword, case_sensitive, context)
        if scope == "section":
            return self.search_sections(keyword, case_sensitive)
        if scope == "title":
            return self.search_titles(keyword, case_sensitive)
        if scope == "all":
            return (
                self.search_paragraphs(keyword, case_sensitive, context)
                + self.search_sections(keyword, case_sensitive)
                + self.search_titles(keyword, case_sensitive)
            )
        raise ValueError(f"未知检索范围: {scope!r}（可用：paragraph/section/title/all）")

    # ------------------------------------------------------------ 结构

    def structure(self, file_name: str) -> Optional[dict[str, Any]]:
        """查看某个文件的章节结构。"""
        return self.index["files"].get(file_name)

    def list_files(self) -> list[str]:
        """列出所有已索引的文件。"""
        return sorted(self.index["files"].keys())

    def titles(self) -> list[str]:
        """列出所有已索引的论文标题。"""
        return [info["title"] for info in self.index["files"].values()]

    # ------------------------------------------------------------ 导出

    def to_csv(self, csv_path: Union[str, Path]) -> Path:
        """导出文献元数据 CSV。"""
        path = Path(csv_path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8-sig", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow([
                "file", "title", "keywords", "n_sections",
                "n_paragraphs", "abstract_chars",
            ])
            for file_name in sorted(self.index["files"]):
                info = self.index["files"][file_name]
                sections = info.get("sections", [])
                n_paragraphs = sum(len(s.get("paragraphs", [])) for s in sections)
                writer.writerow([
                    file_name,
                    info.get("title", ""),
                    "; ".join(info.get("keywords", [])),
                    len(sections),
                    n_paragraphs,
                    len(info.get("abstract", "")),
                ])
        return path

    def save(self, json_path: Union[str, Path]) -> Path:
        """把索引存成 JSON。"""
        path = Path(json_path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(
            json.dumps(self.index, ensure_ascii=False, indent=2), encoding="utf-8"
        )
        return path

    @classmethod
    def load(cls, json_path: Union[str, Path]) -> CorpusIndex:
        """从 JSON 载入索引。"""
        path = Path(json_path)
        data = json.loads(path.read_text(encoding="utf-8"))
        return cls(index=data)

    @property
    def n_files(self) -> int:
        return len(self.index.get("files", {}))


def _collect_sections(root: ET.Element) -> list[dict[str, Any]]:
    """收集章节结构（扁平列表，带 level）。"""
    sections: list[dict[str, Any]] = []

    def walk(parent: ET.Element, level: int) -> None:
        for div in parent.findall(f"{{{TEI_NS['tei']}}}div"):
            head = div.find(f"{{{TEI_NS['tei']}}}head")
            title = element_text(head) if head is not None else ""
            paragraphs = [
                text for text in (
                    element_text(p)
                    for p in div.findall(f"{{{TEI_NS['tei']}}}p")
                ) if text
            ]
            sections.append({
                "title": title,
                "level": level,
                "paragraphs": paragraphs,
            })
            walk(div, level + 1)

    body = root.find(f".//{{{TEI_NS['tei']}}}body")
    if body is not None:
        walk(body, 0)
    return sections


def _context(text: str, keyword: str, length: int, case_sensitive: bool) -> str:
    """截取命中位置前后各若干字符作为上下文。"""
    haystack = text if case_sensitive else text.lower()
    needle = keyword if case_sensitive else keyword.lower()
    position = haystack.find(needle)
    if position < 0:
        return text[: length * 2]

    half = length // 2
    start = max(0, position - half)
    end = min(len(text), position + len(keyword) + half)
    snippet = text[start:end].strip()
    if start > 0:
        snippet = "..." + snippet
    if end < len(text):
        snippet = snippet + "..."
    return snippet


def build_index(
    xml_dir: Union[str, Path],
    pattern: str = "*.tei.xml",
) -> CorpusIndex:
    """便捷函数：为目录建立索引。"""
    return CorpusIndex().build(xml_dir, pattern)


def _build_parser() -> argparse.ArgumentParser:
    """构造命令行解析器。"""
    parser = cli_parser(
        prog="search-corpus",
        description="为 TEI 文献语料建立索引并做全文检索",
        skill="paper-fulltext-extract",
    )
    parser.add_argument("--xml-dir", metavar="DIR", help="TEI XML 目录（建立索引）")
    parser.add_argument("--index-file", metavar="PATH", help="索引 JSON 路径")
    parser.add_argument("--load-index", metavar="PATH", help="从已有索引 JSON 载入")
    parser.add_argument("--keyword", help="检索词")
    parser.add_argument("--scope", default="paragraph",
                        choices=("paragraph", "section", "title", "all"),
                        help="检索范围（默认 paragraph）")
    parser.add_argument("--case-sensitive", action="store_true", help="区分大小写")
    parser.add_argument("--context", type=int, default=100, help="上下文长度（字符）")
    parser.add_argument("--export-csv", metavar="PATH", help="导出文献元数据 CSV")
    parser.add_argument("--list-files", action="store_true", help="列出已索引文件")
    parser.add_argument("--show-structure", metavar="FILE", help="显示某个文件的章节结构")
    parser.add_argument("--json", metavar="PATH", help="把检索结果写成 JSON")
    return parser


def _main(args: argparse.Namespace) -> int:
    """命令行入口：建立索引并检索。"""
    if args.load_index:
        corpus = CorpusIndex.load(args.load_index)
        LOGGER.info("已载入索引：%d 个文件", corpus.n_files)
    elif args.xml_dir:
        corpus = CorpusIndex().build(args.xml_dir)
        LOGGER.info("已建立索引：%d 个文件", corpus.n_files)
        for error in corpus.errors:
            LOGGER.warning("  %s", error)
        if args.index_file:
            corpus.save(args.index_file)
            LOGGER.info("索引已写入 %s", args.index_file)
    else:
        LOGGER.error("需要 --xml-dir 或 --load-index")
        return 1

    if args.export_csv:
        path = corpus.to_csv(args.export_csv)
        LOGGER.info("元数据已导出到 %s", path)

    if args.list_files:
        for name in corpus.list_files():
            LOGGER.info("  %s", name)

    if args.show_structure:
        structure = corpus.structure(args.show_structure)
        if not structure:
            LOGGER.error("未找到文件: %s", args.show_structure)
            return 1
        for section in structure.get("sections", []):
            indent = "  " * section.get("level", 0)
            LOGGER.info(
                "%s%s（%d 段）",
                indent, section.get("title", "(无标题)"),
                len(section.get("paragraphs", [])),
            )

    if args.keyword:
        hits = corpus.search(
            args.keyword,
            scope=args.scope,
            case_sensitive=args.case_sensitive,
            context=args.context,
        )
        LOGGER.info("命中 %d 条", len(hits))
        for hit in hits[:50]:
            LOGGER.info(
                "  [%s] %s :: %s",
                hit.kind, hit.title[:60], hit.match_text[:80] or hit.section_title,
            )
        if len(hits) > 50:
            LOGGER.info("  ...（另有 %d 条未显示）", len(hits) - 50)
        if args.json:
            Path(args.json).write_text(
                json.dumps([h.to_dict() for h in hits], ensure_ascii=False, indent=2),
                encoding="utf-8",
            )
            LOGGER.info("结果已写入 %s", args.json)

    return 0


if __name__ == "__main__":
    sys.exit(run_cli(_main, parser=_build_parser(), skill="paper-fulltext-extract", version="1.0.0"))
