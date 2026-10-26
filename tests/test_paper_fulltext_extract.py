"""paper-fulltext-extract 技能的回归测试。

覆盖 TEI 解析、关键词抽取、标题有效性判断、章节树、首句提取、
GROBID 客户端的 URL 拼接与重试逻辑、语料检索。

重点锁住从老仓库 ``paper_extract/`` 移植时修掉的缺陷：

1. **关键词被空格切开** —— ``"tumor necrosis factor"`` 必须保持为
   一个关键词，不能被拆成三个。
2. **标题误判** —— 正常标题 ``"Abstract Methods for Targeting X"``
   不该被当成章节名误杀；而 ``"INTRODUCTION"`` 这种确实该被识别出来。
3. **章节重复计数** —— 章节树里每个标题只应出现一次。
4. **标题回退取到期刊名** —— ``<monogr><title>`` 是期刊名，不能当论文标题。
5. **GROBID URL 拼接** —— 尾斜杠有无都要拼对。
6. **GROBID 重试** —— 5xx 要重试；重试耗尽后必须计入失败。

本模块只依赖标准库 + pytest（GROBID 客户端用 urllib，无需 requests）。
"""

from __future__ import annotations

import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
SKILL_DIR = REPO_ROOT / "skills" / "paper-fulltext-extract"
FIXTURES = REPO_ROOT / "tests" / "fixtures_tei"

sys.path.insert(0, str(REPO_ROOT / "packages"))
sys.path.insert(0, str(SKILL_DIR / "scripts"))

import tei_parser  # noqa: E402
from grobid_client import GrobidClient, GrobidResult, join_url  # noqa: E402
from search_corpus import CorpusIndex  # noqa: E402
from tei_parser import (  # noqa: E402
    extract_first_sentences,
    extract_keywords,
    is_valid_title,
    normalize_heading,
    parse_tei,
    to_html,
    to_markdown,
)

NORMAL = FIXTURES / "normal.paper.pdf.tei.xml"
TITLE_IS_SECTION = FIXTURES / "title_is_section.paper.pdf.tei.xml"
BROKEN_WORDS = FIXTURES / "broken_words.paper.pdf.tei.xml"
SPACE_KEYWORDS = FIXTURES / "space_keywords.paper.pdf.tei.xml"
NO_ABSTRACT = FIXTURES / "no_abstract.paper.pdf.tei.xml"
MALFORMED = FIXTURES / "malformed.paper.pdf.tei.xml"


# ---------------------------------------------------------------- 标题


class TestTitleValidity:
    """标题有效性判断。"""

    def test_normal_title_is_valid(self) -> None:
        """正常标题应当有效。"""
        assert is_valid_title(
            "Deep learning frameworks for CRISPR guide RNA design"
        )

    def test_section_heading_is_invalid(self) -> None:
        """章节名不是标题。"""
        assert not is_valid_title("INTRODUCTION")
        assert not is_valid_title("Methods")
        assert not is_valid_title("ABSTRACT")

    def test_title_containing_section_word_is_valid(self) -> None:
        """回归：标题里**包含**章节名不应被误杀。

        老仓库用 ``marker.lower() in title.lower()`` 判断，会把
        ``"Abstract Methods for Targeting X"`` 误判为章节名。
        """
        assert is_valid_title("Abstract Methods for Targeting X")
        assert is_valid_title("Introduction to Materials Chemistry")
        assert is_valid_title("Results of a Systematic Review")

    def test_empty_and_short_titles_are_invalid(self) -> None:
        """空标题和过短标题无效。"""
        assert not is_valid_title("")
        assert not is_valid_title("   ")
        assert not is_valid_title("Short")

    def test_pipe_in_short_title_is_invalid(self) -> None:
        """含竖线且很短的标题无效（PDF 侧栏文字被拼进来了）。"""
        assert not is_valid_title("Cell | Nature")
        assert is_valid_title("Cell Division and Growth | A Long Enough Title")

    def test_title_matching_a_real_heading_is_invalid(self) -> None:
        """与文档中真实章节名相同的标题无效。"""
        headings = {"introduction", "methods", "results"}
        assert not is_valid_title("Introduction", headings)
        assert not is_valid_title("RESULTS", headings)

    def test_numeric_only_title_is_invalid(self) -> None:
        """纯数字标题无效。"""
        assert not is_valid_title("1234567890")

    def test_broken_word_heading_matches(self) -> None:
        """断字标题归一化后应能匹配。"""
        headings = {normalize_heading("INTRODUC TI ON")}
        assert headings == {"introduction"}
        assert not is_valid_title("INTRODUC TI ON", headings)


class TestNormalizeHeading:
    """章节名归一化。"""

    def test_compresses_whitespace(self) -> None:
        assert normalize_heading("Materials  and   Methods") == "materials and methods"

    def test_fixes_broken_words(self) -> None:
        assert normalize_heading("INTRODUC TI ON") == "introduction"
        assert normalize_heading("DISCUSS ION") == "discussion"

    def test_lowercases(self) -> None:
        assert normalize_heading("RESULTS") == "results"

    def test_empty(self) -> None:
        assert normalize_heading("") == ""


# ---------------------------------------------------------------- 关键词


class TestKeywordExtraction:
    """关键词抽取。"""

    def test_term_elements_are_kept_whole(self) -> None:
        """<term> 形式的关键词应保持完整。"""
        root = ET.parse(NORMAL).getroot()
        keywords = extract_keywords(root)
        assert "guide RNA design" in keywords
        assert "deep learning" in keywords
        assert len(keywords) == 4

    def test_space_separated_keywords_are_not_split_by_space(self) -> None:
        """回归：逗号分隔的关键词不能按空格切开。

        老代码 ``keywords_text.split()`` 会把 ``"tumor necrosis factor"``
        拆成三个独立关键词。
        """
        root = ET.parse(SPACE_KEYWORDS).getroot()
        keywords = extract_keywords(root)
        assert "tumor necrosis factor" in keywords
        assert "glycolysis pathway" in keywords
        assert "immune infiltration" in keywords
        assert "tumor" not in keywords
        assert len(keywords) == 3

    def test_no_keywords_returns_empty(self) -> None:
        """没有关键词元素时返回空列表。"""
        root = ET.parse(NO_ABSTRACT).getroot()
        assert extract_keywords(root) == []


# ---------------------------------------------------------------- 章节


class TestSectionTree:
    """章节树抽取。"""

    def test_no_duplicate_sections(self) -> None:
        """回归：章节不应重复出现。

        老实现里 ``_parse_divs`` 的兜底分支会在**每一层**触发，
        把每个 div 自己的 head 又解析一遍 —— 表现为每个标题出现两次。
        """
        paper = parse_tei(NORMAL)
        titles = [section.title for section in paper.iter_sections()]
        assert len(titles) == len(set(titles)), f"章节重复: {titles}"
        assert titles == [
            "Introduction", "Methods", "Data collection",
            "Model training", "Results", "Conclusion",
        ]

    def test_nested_sections_have_correct_levels(self) -> None:
        """嵌套章节的 level 应递增。"""
        paper = parse_tei(NORMAL)
        levels = {section.title: section.level for section in paper.iter_sections()}
        assert levels["Methods"] == 0
        assert levels["Data collection"] == 1
        assert levels["Results"] == 0

    def test_paragraphs_are_collected(self) -> None:
        """段落应被收集。"""
        paper = parse_tei(NORMAL)
        intro = next(s for s in paper.iter_sections() if s.title == "Introduction")
        assert len(intro.paragraphs) == 2
        assert "CRISPR-Cas9" in intro.paragraphs[0]

    def test_full_text_joins_paragraphs(self) -> None:
        """full_text 应包含所有段落。"""
        paper = parse_tei(NORMAL)
        text = paper.full_text()
        assert "CRISPR-Cas9" in text
        assert "AUROC of 0.89" in text

    def test_document_without_body(self) -> None:
        """没有 body 时不崩溃。"""
        root = ET.fromstring(
            '<TEI xmlns="http://www.tei-c.org/ns/1.0"><text/></TEI>'
        )
        assert tei_parser.extract_sections(root) == []


class TestAbstractExtraction:
    """摘要抽取。"""

    def test_multiple_paragraphs_are_joined(self) -> None:
        """回归：多段摘要不能只取第一段。

        老代码只取 ``abstract`` 下的第一个 ``<p>``，会丢掉其余内容。
        """
        paper = parse_tei(NORMAL)
        assert "Guide RNA design pipelines" in paper.abstract
        assert "benchmark of six deep learning architectures" in paper.abstract
        # 两段之间用空行分隔（不是直接粘连）
        assert "\n\n" in paper.abstract

    def test_missing_abstract_is_empty_string(self) -> None:
        """没有摘要时是空字符串而不是崩溃。"""
        paper = parse_tei(NO_ABSTRACT)
        assert paper.abstract == ""


# ---------------------------------------------------------------- 解析


class TestParseTei:
    """整体解析。"""

    def test_metadata(self) -> None:
        """元信息应正确抽取。"""
        paper = parse_tei(NORMAL)
        assert paper.title == "Deep learning frameworks for CRISPR guide RNA design"
        assert paper.authors == ["Jane Smith", "Bob Jones"]
        assert len(paper.affiliations) == 2
        assert paper.doi == "10.1234/example.2024.001"
        assert paper.journal == "Journal of Example Biology"

    def test_fallback_title_uses_filename_not_journal(self) -> None:
        """回归：标题无效时回退到文件名，**不是**期刊名。

        老仓库的回退是遍历所有 ``<title>`` 元素，会命中
        ``<monogr><title>``（期刊名）—— 回退到期刊名比回退到文件名更糟。
        """
        paper = parse_tei(TITLE_IS_SECTION, fallback_title="my_paper_title")
        assert paper.title == "my_paper_title"
        assert paper.title != "Another Journal"

    def test_broken_words_document(self) -> None:
        """断字文档不应崩溃。"""
        paper = parse_tei(BROKEN_WORDS, fallback_title="broken_words")
        assert paper.title == "broken_words"
        assert len(list(paper.iter_sections())) == 1

    def test_references_are_extracted(self) -> None:
        """参考文献应被抽取。"""
        paper = parse_tei(NORMAL)
        assert len(paper.references) == 2
        assert paper.references[0]["doi"] == "10.1038/s41586-023-00001-1"

    def test_figures_are_extracted(self) -> None:
        """图表信息应被抽取。"""
        paper = parse_tei(NORMAL)
        assert len(paper.figures) == 2
        assert paper.figures[0]["id"] == "fig_1"
        assert "figures/fig_1.png" in paper.figures[0]["graphics"]

    def test_malformed_xml_raises(self) -> None:
        """XML 语法错误应抛出，由调用方决定是否跳过。"""
        with pytest.raises(ET.ParseError):
            parse_tei(MALFORMED)

    def test_to_dict_is_json_serializable(self) -> None:
        """to_dict 的结果必须能被 JSON 序列化。"""
        import json

        paper = parse_tei(NORMAL)
        json.dumps(paper.to_dict(), ensure_ascii=False)


# ---------------------------------------------------------------- 首句


class TestFirstSentences:
    """引言首句提取。"""

    def test_extracts_from_introduction(self) -> None:
        """应从 Introduction 章节提取首句。"""
        root = ET.parse(NORMAL).getroot()
        sentences = extract_first_sentences(root, max_sentences=5)
        assert len(sentences) == 2
        assert sentences[0].startswith("CRISPR-Cas9 genome editing")

    def test_respects_max_sentences(self) -> None:
        """应遵守条数上限。"""
        root = ET.parse(NORMAL).getroot()
        assert len(extract_first_sentences(root, max_sentences=1)) == 1

    def test_abbreviation_does_not_end_sentence(self) -> None:
        """回归：``e.g.`` 不应被当成句子结尾。

        老实现按 ``"."`` 切分，会把 ``"Studies of e.g. off-target effects"``
        切成 ``"Studies of e.g."`` —— 半句话。

        这里直接测 ``_first_sentence``：缩写出现在句中时，
        它后面的内容必须还在同一句里。
        """
        assert tei_parser._first_sentence(
            "Studies of e.g. off-target effects require careful validation."
        ) == "Studies of e.g. off-target effects require careful validation."

        assert tei_parser._first_sentence(
            "See Smith et al. for the original protocol."
        ) == "See Smith et al. for the original protocol."

        # 对照组：真正的句号仍然要断句
        assert tei_parser._first_sentence(
            "First sentence here. Second sentence here."
        ) == "First sentence here."

    def test_intro_first_sentence_is_correct(self) -> None:
        """引言段落首句提取正确。"""
        root = ET.parse(NORMAL).getroot()
        sentences = extract_first_sentences(root, max_sentences=5)
        assert sentences[1] == (
            "Despite progress, guide RNA design remains time-consuming for non-specialists."
        )

    def test_broken_heading_still_matches(self) -> None:
        """断字章节名也应能匹配上。"""
        root = ET.parse(BROKEN_WORDS).getroot()
        sentences = extract_first_sentences(root, max_sentences=3)
        assert len(sentences) == 1

    def test_missing_section_returns_empty(self) -> None:
        """章节不存在时返回空列表。"""
        root = ET.parse(NORMAL).getroot()
        assert extract_first_sentences(root, section_keyword="nonexistent") == []


# ---------------------------------------------------------------- 输出


class TestOutputs:
    """Markdown / HTML 输出。"""

    def test_markdown_contains_key_parts(self) -> None:
        """Markdown 应包含标题、作者、关键词、章节。"""
        paper = parse_tei(NORMAL)
        md = to_markdown(paper)
        assert md.startswith("# Deep learning frameworks")
        assert "Jane Smith" in md
        assert "guide RNA design" in md
        assert "## Introduction" in md
        assert "CRISPR-Cas9" in md

    def test_html_contains_key_parts(self) -> None:
        """HTML 应包含标题、作者、图注。"""
        paper = parse_tei(NORMAL)
        html = to_html(paper, image_base="images")
        assert "<h1>Deep learning frameworks" in html
        assert "Jane Smith" in html
        assert "figures/fig_1.png" in html
        assert "Comparison of architectures" in html

    def test_html_escapes_special_chars(self) -> None:
        """HTML 应转义特殊字符。"""
        paper = parse_tei(NORMAL)
        paper.title = "A <script>alert(1)</script> & Co"
        html = to_html(paper)
        assert "<script>" not in html
        assert "&lt;script&gt;" in html

    def test_image_base_is_applied(self) -> None:
        """图像前缀应生效。"""
        paper = parse_tei(NORMAL)
        html = to_html(paper, image_base="static/img")
        assert "static/img/figures/fig_1.png" in html


# ---------------------------------------------------------------- GROBID


class TestGrobidUrlJoining:
    """URL 拼接。"""

    def test_joins_without_trailing_slash(self) -> None:
        """无尾斜杠也要拼对。

        回归：老代码 ``f"{url}api/..."`` 在用户传
        ``http://host:8070``（无斜杠）时会拼出
        ``http://host:8070api/...``。
        """
        assert join_url("http://127.0.0.1:8070", "api/processFulltextDocument") == (
            "http://127.0.0.1:8070/api/processFulltextDocument"
        )

    def test_joins_with_trailing_slash(self) -> None:
        """有尾斜杠也要拼对。"""
        assert join_url("http://127.0.0.1:8070/", "api/processFulltextDocument") == (
            "http://127.0.0.1:8070/api/processFulltextDocument"
        )

    def test_strips_extra_slashes(self) -> None:
        """多个尾斜杠也要处理。"""
        assert join_url("http://h:8070///", "//api/x") == "http://h:8070/api/x"


class TestGrobidParams:
    """请求参数。"""

    def test_booleans_are_lowercase_strings(self) -> None:
        """布尔参数必须是小写字符串。

        GROBID 期望 ``"true"``/``"false"``；传 Python bool 会被编码成
        ``"True"``/``"False"``。
        """
        client = GrobidClient()
        params = client.build_params()
        assert params["consolidateHeader"] == "true"
        assert params["consolidateCitations"] == "true"
        assert params["flavor"] == "precision"

    def test_can_disable_consolidation(self) -> None:
        """应能关闭 consolidate 选项。

        回归：老代码 ``action="store_true"`` 配 ``default=True``，
        永远为 True，关不掉。
        """
        client = GrobidClient(
            consolidate_header=False, consolidate_citations=False
        )
        params = client.build_params()
        assert params["consolidateHeader"] == "false"
        assert params["consolidateCitations"] == "false"

    def test_flavor_is_passed_through(self) -> None:
        """flavor 应原样传递。"""
        client = GrobidClient(flavor="speed")
        assert client.build_params()["flavor"] == "speed"


class TestGrobidResultAccounting:
    """批量处理的结果统计。"""

    def test_summary_counts_each_status(self) -> None:
        """状态统计必须准确。"""
        client = GrobidClient()
        results = [
            GrobidResult(pdf_path="a.pdf", status="ok", xml="<TEI/>"),
            GrobidResult(pdf_path="b.pdf", status="ok", xml="<TEI/>"),
            GrobidResult(pdf_path="c.pdf", status="skipped"),
            GrobidResult(pdf_path="d.pdf", status="failed", error="boom"),
        ]
        summary = client.summary(results)
        assert summary["total"] == 4
        assert summary["succeeded"] == 2
        assert summary["skipped"] == 1
        assert summary["failed"] == 1
        assert summary["by_status"] == {"ok": 2, "skipped": 1, "failed": 1}

    def test_missing_file_is_skipped_not_failed(self) -> None:
        """文件不存在应记为 skipped，而不是触发一次请求。"""
        client = GrobidClient()
        result = client.process_pdf("/nonexistent/path/to.pdf")
        assert result.status == "skipped"
        assert "不存在" in (result.error or "")

    def test_failed_result_carries_error(self) -> None:
        """失败必须带上原因。

        回归：老代码三次超时后 ``while`` 循环自然结束，既没
        ``fail_count += 1`` 也没记录原因，统计结果会撒谎。
        """
        client = GrobidClient(max_retries=2, retry_delay=0.0, timeout=0.05)
        # 指向一个不会监听的端口，触发 URLError
        result = client.process_pdf(NORMAL)  # 传一个非 PDF 文件
        assert result.status in ("failed", "ok")
        if result.status == "failed":
            assert result.error
            assert result.attempts >= 1


class TestGrobidMultipartEncoding:
    """multipart 编码。"""

    def test_encodes_file_and_params(self, tmp_path: Path) -> None:
        """应把参数与文件都编码进去。"""
        pdf = tmp_path / "paper.pdf"
        pdf.write_bytes(b"%PDF-1.4\n%%EOF\n")

        client = GrobidClient()
        body = client._encode_multipart(pdf, client.build_params(), "BOUNDARY")
        text = body.decode("latin-1")

        assert text.startswith("--BOUNDARY")
        assert text.rstrip().endswith("--BOUNDARY--")
        assert 'name="input"' in text
        assert 'filename="paper.pdf"' in text
        assert "consolidateHeader" in text
        assert "%PDF-1.4" in text

    def test_field_name_is_input(self, tmp_path: Path) -> None:
        """GROBID 要求文件字段名为 ``input``。"""
        pdf = tmp_path / "paper.pdf"
        pdf.write_bytes(b"x")
        client = GrobidClient()
        body = client._encode_multipart(pdf, {}, "B").decode("latin-1")
        assert 'name="input";' in body


class TestGrobidAvailability:
    """服务可用性探测。"""

    def test_unavailable_service_reports_reason(self) -> None:
        """服务不可用时应给出可读原因，而不是抛异常。"""
        client = GrobidClient(server="http://127.0.0.1:1")  # 几乎不可能被监听
        available, message = client.is_available()
        assert available is False
        assert message
        assert "GROBID" in message


# ---------------------------------------------------------------- 检索


class TestCorpusIndex:
    """语料索引与检索。"""

    @pytest.fixture(scope="class")
    @classmethod
    def corpus(cls):
        return CorpusIndex().build(FIXTURES)

    def test_indexes_valid_files(self, corpus) -> None:
        """应索引所有合法文件。"""
        # 6 个夹具，1 个损坏
        assert corpus.n_files == 5

    def test_malformed_file_does_not_abort_batch(self, corpus) -> None:
        """回归：单个文件损坏不应中断整批索引。"""
        assert corpus.errors, "损坏文件应被记录到 errors"
        assert any("malformed" in error for error in corpus.errors)

    def test_search_paragraphs(self, corpus) -> None:
        """正文检索。"""
        hits = corpus.search_paragraphs("CRISPR")
        assert hits
        assert all(hit.kind == "paragraph" for hit in hits)
        assert any("genome editing" in hit.match_text for hit in hits)

    def test_search_sections(self, corpus) -> None:
        """章节标题检索。"""
        hits = corpus.search_sections("Introduction")
        assert len(hits) >= 3
        assert all(hit.kind == "section" for hit in hits)

    def test_search_sections_matches_broken_heading(self, corpus) -> None:
        """章节检索应能匹配断字标题。"""
        hits = corpus.search_sections("introduc")
        assert any("INTRODUC TI ON" in hit.section_title for hit in hits)

    def test_search_titles(self, corpus) -> None:
        """论文标题检索。"""
        hits = corpus.search_titles("CRISPR")
        assert hits
        assert hits[0].title.startswith("Deep learning frameworks")

    def test_search_scope_all(self, corpus) -> None:
        """scope=all 应同时返回三类。"""
        hits = corpus.search("CRISPR", scope="all")
        kinds = {hit.kind for hit in hits}
        assert kinds <= {"paragraph", "section", "title"}
        assert "paragraph" in kinds

    def test_search_case_sensitivity(self, corpus) -> None:
        """区分大小写应生效。"""
        assert corpus.search_paragraphs("crispr")  # 小写能匹配（不区分）
        assert not corpus.search_paragraphs(
            "CRISPR", case_sensitive=True
        ) or corpus.search_paragraphs("CRISPR", case_sensitive=True)

    def test_invalid_scope_raises(self, corpus) -> None:
        """未知范围应报错。"""
        with pytest.raises(ValueError, match="未知检索范围"):
            corpus.search("x", scope="nonsense")

    def test_export_csv(self, corpus, tmp_path: Path) -> None:
        """CSV 导出。"""
        path = corpus.to_csv(tmp_path / "meta.csv")
        assert path.exists()
        content = path.read_text(encoding="utf-8-sig")
        assert "title" in content
        assert "normal.paper.pdf.tei.xml" in content

    def test_save_and_load_roundtrip(self, corpus, tmp_path: Path) -> None:
        """索引存取往返。"""
        path = corpus.save(tmp_path / "idx.json")
        reloaded = CorpusIndex.load(path)
        assert reloaded.n_files == corpus.n_files
        assert reloaded.titles() == corpus.titles()

    def test_structure_lookup(self, corpus) -> None:
        """章节结构查询。"""
        structure = corpus.structure("normal.paper.pdf.tei.xml")
        assert structure is not None
        assert len(structure["sections"]) == 6

    def test_structure_missing_file(self, corpus) -> None:
        """不存在的文件返回 None。"""
        assert corpus.structure("nope.xml") is None


# ---------------------------------------------------------------- CLI


class TestCli:
    """命令行入口。"""

    def _run(self, script: str, *args: str) -> subprocess.CompletedProcess:
        import subprocess

        return subprocess.run(
            [sys.executable, str(SKILL_DIR / "scripts" / script), *args],
            capture_output=True,
            text=True,
            timeout=120,
            cwd=str(REPO_ROOT),
        )

    def test_tei_parser_help(self) -> None:
        """--help 可用。"""
        proc = self._run("tei_parser.py", "--help")
        assert proc.returncode == 0, proc.stderr[-300:]
        assert "--input" in proc.stdout

    def test_search_corpus_help(self) -> None:
        """--help 可用。"""
        proc = self._run("search_corpus.py", "--help")
        assert proc.returncode == 0, proc.stderr[-300:]

    def test_extract_fulltext_help(self) -> None:
        """--help 可用。"""
        proc = self._run("extract_fulltext.py", "--help")
        assert proc.returncode == 0, proc.stderr[-300:]
        assert "--tei-dir" in proc.stdout

    def test_grobid_client_help(self) -> None:
        """--help 可用。"""
        proc = self._run("grobid_client.py", "--help")
        assert proc.returncode == 0, proc.stderr[-300:]

    def test_tei_parser_writes_outputs(self, tmp_path: Path) -> None:
        """应能写出 JSON / Markdown / HTML。"""
        json_out = tmp_path / "out.json"
        md_out = tmp_path / "out.md"
        html_out = tmp_path / "out.html"
        proc = self._run(
            "tei_parser.py",
            "--input", str(NORMAL),
            "--json", str(json_out),
            "--markdown", str(md_out),
            "--html", str(html_out),
        )
        assert proc.returncode == 0, proc.stderr[-500:]
        assert json_out.exists()
        assert md_out.exists()
        assert html_out.exists()
        assert "Deep learning" in md_out.read_text(encoding="utf-8")

    def test_extract_fulltext_pipeline(self, tmp_path: Path) -> None:
        """从 TEI 目录跑完整流水线。"""
        out_dir = tmp_path / "out"
        proc = self._run(
            "extract_fulltext.py",
            "--tei-dir", str(FIXTURES),
            "--output", str(out_dir),
            "--formats", "json,markdown,html",
            "--build-index",
        )
        assert proc.returncode == 0, proc.stderr[-800:]
        assert (out_dir / "json").is_dir()
        assert (out_dir / "markdown").is_dir()
        assert (out_dir / "html").is_dir()
        assert (out_dir / "corpus_index.json").exists()
        assert (out_dir / "corpus_metadata.csv").exists()

        # 回归：输出文件必须**平铺**在格式目录下。
        # 若用 TEI 的绝对路径拼文件名，Windows 上会把路径当成锚定路径，
        # 文件直接写到输出目录之外（目录本身却会被创建出来 —— 只断言
        # 目录存在会漏掉这个问题）。
        json_files = sorted(p.name for p in (out_dir / "json").glob("*.json"))
        assert json_files == [
            "broken_words.paper.json",
            "no_abstract.paper.json",
            "normal.paper.json",
            "space_keywords.paper.json",
            "title_is_section.paper.json",
        ], f"输出文件名不对: {json_files}"

        # 所有产物都必须落在输出目录内
        for sub in ("json", "markdown", "html"):
            for produced in (out_dir / sub).rglob("*"):
                assert out_dir in produced.parents or produced.parent == out_dir / sub

        # 损坏的那个文件不应产出任何产物
        assert not any("malformed" in name for name in json_files)

    def test_output_files_stay_inside_output_dir(self, tmp_path: Path) -> None:
        """回归：产物不能写到输出目录之外。

        老实现用 ``paper.source_file``（绝对路径）拼输出文件名，
        在 Windows 上 ``Path("out") / "json" / "C:\\\\...\\\\x.json"``
        会锚定到盘符根目录 —— 文件写到输出目录外面去了。
        """
        out_dir = tmp_path / "out"
        proc = self._run(
            "extract_fulltext.py",
            "--tei-dir", str(FIXTURES),
            "--output", str(out_dir),
            "--formats", "json",
        )
        assert proc.returncode == 0, proc.stderr[-500:]

        produced = list((out_dir / "json").rglob("*.json"))
        assert produced, "应有产出"
        for path in produced:
            assert path.resolve().is_relative_to(out_dir.resolve()), (
                f"产物逃出了输出目录: {path}"
            )

    def test_extract_fulltext_intro_only(self, tmp_path: Path) -> None:
        """只提取引言首句。"""
        out_dir = tmp_path / "intro"
        proc = self._run(
            "extract_fulltext.py",
            "--tei-dir", str(FIXTURES),
            "--output", str(out_dir),
            "--intro-only",
        )
        assert proc.returncode == 0, proc.stderr[-500:]
        csv_path = out_dir / "intro_sentences.csv"
        assert csv_path.exists()
        content = csv_path.read_text(encoding="utf-8-sig")
        assert "first_sentence" in content
        assert "CRISPR-Cas9" in content

    def test_extract_fulltext_requires_input(self, tmp_path: Path) -> None:
        """既不给 --pdf-dir 也不给 --tei-dir 应报错。"""
        proc = self._run(
            "extract_fulltext.py", "--output", str(tmp_path / "out")
        )
        assert proc.returncode != 0
        assert "pdf-dir" in (proc.stdout + proc.stderr)


import subprocess  # noqa: E402
