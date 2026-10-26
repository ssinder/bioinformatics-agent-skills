"""文献全文抽取流水线：PDF → GROBID → TEI XML → JSON / Markdown / HTML。

这是技能的主入口，把各步骤串成一条流水线：

.. code-block:: text

    PDF 目录
       │  grobid_client（GROBID 服务，本地部署）
       ▼
    TEI XML 目录
       │  tei_parser
       ▼
    JSON  ──► Markdown
       └─► HTML（带图注）
       │  search_corpus
       ▼
    全文索引 + 检索结果 CSV

单独跑某一步也可以 —— 每个 ``--from`` 决定起点。

用法::

    # 全流程
    python extract_fulltext.py --pdf-dir ./pdfs --output ./out

    # 只从 TEI XML 开始
    python extract_fulltext.py --tei-dir ./tei --output ./out

    # 只提取引言首句
    python extract_fulltext.py --tei-dir ./tei --output ./out --intro-only
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "packages"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from typing import Optional  # noqa: E402

from biopipe_core import cli_parser, get_logger, run_cli  # noqa: E402
from biopipe_core.cli import SkillError  # noqa: E402
from grobid_client import GrobidClient  # noqa: E402
from search_corpus import CorpusIndex  # noqa: E402
from tei_parser import (  # noqa: E402
    ParsedPaper,
    extract_first_sentences,
    parse_tei,
    to_html,
    to_markdown,
)

SKILL_NAME = "paper-fulltext-extract"

# 断字修复的正则（与 tei_parser 保持一致的思路）
_TITLE_FALLBACK_RE = re.compile(r"\.pdf\.tei$|\.tei$|\.pdf$")


def _build_parser() -> argparse.ArgumentParser:
    """构造命令行解析器。"""
    parser = cli_parser(
        prog=f"{SKILL_NAME}",
        description="文献全文抽取流水线：PDF → GROBID → TEI XML → JSON/Markdown/HTML",
        skill=SKILL_NAME,
    )
    parser.add_argument("--pdf-dir", metavar="DIR", help="PDF 输入目录（触发 GROBID 步骤）")
    parser.add_argument("--tei-dir", metavar="DIR", help="TEI XML 目录（跳过 GROBID 步骤）")
    parser.add_argument("--output", required=True, metavar="DIR", help="输出目录")

    parser.add_argument("--grobid-server", default="http://127.0.0.1:8070",
                        help="GROBID 服务地址（默认 http://127.0.0.1:8070）")
    parser.add_argument("--flavor", choices=("precision", "speed"), default="precision",
                        help="GROBID 抽取精度（默认 precision）")
    parser.add_argument("--no-consolidate-header", action="store_true",
                        help="关闭 consolidateHeader")
    parser.add_argument("--no-consolidate-citations", action="store_true",
                        help="关闭 consolidateCitations")
    parser.add_argument("--pause", type=float, default=0.0,
                        help="GROBID 每篇之间的等待秒数")

    parser.add_argument("--formats", default="json,markdown",
                        help="输出格式，逗号分隔：json / markdown / html（默认 json,markdown）")
    parser.add_argument("--image-base", default="images",
                        help="HTML 中图像的相对目录前缀（默认 images）")
    parser.add_argument("--intro-only", action="store_true",
                        help="只提取引言首句，不产出完整文档")
    parser.add_argument("--intro-section", default="introduction",
                        help="引言章节名（默认 introduction）")
    parser.add_argument("--max-intro-sentences", type=int, default=3,
                        help="每篇最多提取的引言首句数（默认 3）")
    parser.add_argument("--build-index", action="store_true",
                        help="建立全文检索索引")
    parser.add_argument("--search", metavar="KEYWORD", help="检索词（需配合 --build-index）")
    parser.add_argument("--search-scope", default="paragraph",
                        choices=("paragraph", "section", "title", "all"),
                        help="检索范围（默认 paragraph）")
    parser.add_argument("--skip-grobid-check", action="store_true",
                        help="跳过 GROBID 可用性检查（服务已确认在线时用）")
    return parser


def _discover_tei_files(tei_dir: Path) -> list[Path]:
    """找出目录下所有 TEI XML 文件。"""
    files = sorted(tei_dir.glob("*.tei.xml"))
    files.extend(p for p in sorted(tei_dir.glob("*.xml")) if p not in files)
    return files


def _title_from_filename(path: Path) -> str:
    """从文件名推出一个可用的标题（作为兜底）。"""
    stem = path.name
    for suffix in (".pdf.tei.xml", ".tei.xml", ".xml", ".pdf"):
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    return stem.replace("_", " ").replace("-", " ").strip() or path.stem


def _sanitize_stem(stem: str) -> str:
    """把文件名净化成可安全用作输出文件名的 stem。"""
    cleaned = re.sub(r'[<>:"/\\|?*\x00-\x1f]', "_", stem).strip(" .")
    return cleaned or "paper"


def _write_outputs(
    paper: ParsedPaper,
    output_dir: Path,
    formats: set[str],
    image_base: str,
) -> dict[str, str]:
    """按指定格式写出文件。

    输出文件名必须从**文件名**（而不是完整路径）推导。
    :attr:`ParsedPaper.source_file` 存的是 TEI 文件的绝对路径，
    直接拿它拼 ``output_dir / "json" / f"{source_file}.json"`` 会得到
    ``out/json/C:/Users/.../paper.json`` 这种带盘符和目录层级的路径 ——
    在 Windows 上更糟，``Path`` 会把它当成锚定路径，文件直接写到
    输出目录**之外**去。
    """
    written: dict[str, str] = {}
    stem = Path(paper.source_file).name
    for suffix in (".pdf.tei.xml", ".tei.xml", ".xml", ".pdf"):
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    stem = _sanitize_stem(stem)

    if "json" in formats:
        path = output_dir / "json" / f"{stem}.json"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(
            json.dumps(paper.to_dict(), ensure_ascii=False, indent=2), encoding="utf-8"
        )
        written["json"] = str(path)

    if "markdown" in formats:
        path = output_dir / "markdown" / f"{stem}.md"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(to_markdown(paper), encoding="utf-8")
        written["markdown"] = str(path)

    if "html" in formats:
        path = output_dir / "html" / f"{stem}.html"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(to_html(paper, image_base=image_base), encoding="utf-8")
        written["html"] = str(path)

    return written


def _extract_intro_csv(
    tei_files: list[Path],
    output_path: Path,
    section_keyword: str,
    max_sentences: int,
) -> int:
    """提取引言首句并写 CSV，返回有效条目数。"""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    rows = 0
    with output_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["file", "title", "sentence_index", "first_sentence"])
        for xml_file in tei_files:
            try:
                import xml.etree.ElementTree as ET

                root = ET.parse(xml_file).getroot()
                sentences = extract_first_sentences(
                    root, section_keyword=section_keyword,
                    max_sentences=max_sentences,
                )
            except Exception:  # noqa: BLE001
                continue
            if not sentences:
                continue
            title = _title_from_filename(xml_file)
            for index, sentence in enumerate(sentences, start=1):
                writer.writerow([xml_file.name, title, index, sentence])
                rows += 1
    return rows


def main(args: argparse.Namespace) -> int:
    """流水线主入口。"""
    logger = get_logger(f"biopipe.{SKILL_NAME}")

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    if not args.pdf_dir and not args.tei_dir:
        raise SkillError("需要提供 --pdf-dir 或 --tei-dir 之一")

    # ---------------------------------------------------------- GROBID
    if args.pdf_dir:
        pdf_dir = Path(args.pdf_dir)
        if not pdf_dir.is_dir():
            raise SkillError(f"PDF 目录不存在: {pdf_dir}")

        client = GrobidClient(
            server=args.grobid_server,
            flavor=args.flavor,
            consolidate_header=not args.no_consolidate_header,
            consolidate_citations=not args.no_consolidate_citations,
        )

        if not args.skip_grobid_check:
            available, message = client.is_available()
            logger.info("%s", message)
            if not available:
                logger.error(
                    "GROBID 服务不可用。可以：\n"
                    "  1) 启动服务：docker run -it -p 8070:8070 kermitt/grobid\n"
                    "  2) 或跳过此步，改用 --tei-dir 指向已有的 TEI XML"
                )
                return 1

        tei_dir = output_dir / "tei"
        logger.info("GROBID 处理 %s → %s", pdf_dir, tei_dir)

        def progress(result, index: int, total: int) -> None:
            marker = {"ok": "✓", "skipped": "-", "failed": "✗"}.get(result.status, "?")
            logger.info("[%d/%d] %s %s", index, total, marker, Path(result.pdf_path).name)
            if result.error:
                logger.warning("    %s", result.error)

        results = client.process_directory(
            pdf_dir, tei_dir, pause=args.pause, on_progress=progress
        )
        stats = client.summary(results)
        logger.info(
            "GROBID 完成：共 %d，成功 %d，跳过 %d，失败 %d",
            stats["total"], stats["succeeded"], stats["skipped"], stats["failed"],
        )
        if stats["failed"] and not args.tei_dir:
            logger.warning("有 %d 篇失败，后续只处理成功的部分", stats["failed"])
    else:
        tei_dir = Path(args.tei_dir)

    if not tei_dir.is_dir():
        raise SkillError(f"TEI XML 目录不存在: {tei_dir}")

    tei_files = _discover_tei_files(tei_dir)
    if not tei_files:
        raise SkillError(f"目录中没有找到 TEI XML 文件: {tei_dir}")
    logger.info("找到 %d 个 TEI XML 文件", len(tei_files))

    # ---------------------------------------------------------- 引言提取
    if args.intro_only:
        csv_path = output_dir / "intro_sentences.csv"
        rows = _extract_intro_csv(
            tei_files, csv_path, args.intro_section, args.max_intro_sentences
        )
        logger.info("引言首句已写入 %s（%d 条）", csv_path, rows)
        return 0

    # ---------------------------------------------------------- 解析输出
    formats = {f.strip().lower() for f in args.formats.split(",") if f.strip()}
    valid_formats = {"json", "markdown", "html"}
    unknown = formats - valid_formats
    if unknown:
        raise SkillError(
            f"未知的输出格式: {sorted(unknown)}（可用：{sorted(valid_formats)}）"
        )

    parsed: list[ParsedPaper] = []
    failures = 0
    for xml_file in tei_files:
        try:
            paper = parse_tei(xml_file, fallback_title=_title_from_filename(xml_file))
        except Exception as exc:  # noqa: BLE001
            logger.warning("解析 %s 失败: %s", xml_file.name, exc)
            failures += 1
            continue
        parsed.append(paper)
        _write_outputs(paper, output_dir, formats, args.image_base)

    logger.info("解析完成：成功 %d 篇，失败 %d 篇", len(parsed), failures)

    # ---------------------------------------------------------- 索引检索
    if args.build_index or args.search:
        corpus = CorpusIndex().build(tei_dir)
        index_path = output_dir / "corpus_index.json"
        corpus.save(index_path)
        logger.info("索引已建立：%d 个文件 → %s", corpus.n_files, index_path)
        for error in corpus.errors[:10]:
            logger.warning("  %s", error)

        meta_csv = output_dir / "corpus_metadata.csv"
        corpus.to_csv(meta_csv)
        logger.info("元数据 CSV → %s", meta_csv)

        if args.search:
            hits = corpus.search(args.search, scope=args.search_scope)
            logger.info("检索 %r 命中 %d 条", args.search, len(hits))
            hits_path = output_dir / "search_hits.json"
            hits_path.write_text(
                json.dumps([h.to_dict() for h in hits], ensure_ascii=False, indent=2),
                encoding="utf-8",
            )
            logger.info("检索结果 → %s", hits_path)
            for hit in hits[:20]:
                logger.info(
                    "  [%s] %s :: %s",
                    hit.kind, hit.title[:60],
                    hit.match_text[:80] or hit.section_title,
                )

    return 0


if __name__ == "__main__":
    from biopipe_core.cli import run_cli

    sys.exit(run_cli(main, parser=_build_parser(), skill=SKILL_NAME, version="1.0.0"))
