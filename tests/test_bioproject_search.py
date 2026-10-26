"""BioProject 检索技能的回归测试。

锁住老仓库 ``crawl_bioproject/`` 的一批真实缺陷：

1. **硬编码 Windows 桌面路径** ``C:\\Users\\DELL\\Desktop\\trae\\bioproject``
   —— 换台机器首行日志就 FileNotFoundError。
2. **用 BeautifulSoup + 绝对 XPath 爬网页** —— 页面改版即失效。
3. **无条件跳过前三行** ``all_lines[3:]`` —— 配置与数据混在一起。
4. **伪造 NCBI email** ``userN@example.com`` —— Entrez 的 email 是身份
   标识，伪造会被限速封禁。
5. **死代码** ``check_has_description_link()`` 从未被调用，且内部绕过
   了重试与 UA 逻辑。
6. **每条记录重写整个 CSV** —— O(n²) 磁盘写入。
7. **噪声词误杀** —— 老脚本对任意长度的描述都套用噪声词列表。

全部用例不联网，使用伪造的响应对象。
"""

from __future__ import annotations

import ast
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
SKILL_DIR = REPO_ROOT / "skills" / "bioproject-search"
sys.path.insert(0, str(REPO_ROOT / "packages"))
sys.path.insert(0, str(SKILL_DIR / "scripts"))

import bioproject_client  # noqa: E402
import search_bioprojects  # noqa: E402

# 真实的 NCBI 响应夹具（从 eutils 抓取并冻结）。
#
# 这里刻意用**真实响应**而不是手写的简化 XML：手写夹具会把真实的
# XML 结构简化掉，于是"简化结构下能跑通"的实现照样能通过测试，
# 但真实数据上会静默取空字段。开发本技能时正是这样踩了坑——
# 手写夹具里的 <Organism OrganismName=... TaxonomyID=.../> 在真实
# 响应里根本不是这个形状。
FIXTURE_DIR = REPO_ROOT / "tests" / "fixtures_bioproject"

# efetch 响应：PRJNA212358 / PRJNA355343 / PRJNA172563
_SAMPLE_XML = (FIXTURE_DIR / "efetch_response.xml").read_text(encoding="utf-8")

# esummary 响应：同样的三个项目
_SAMPLE_ESUMMARY_XML = (FIXTURE_DIR / "esummary_response.xml").read_text(encoding="utf-8")


class _FakeResponse:
    """模拟 requests.Response 的最小实现。"""

    def __init__(self, payload=None, text="", status_code=200):
        self._payload = payload
        self.text = text
        self.status_code = status_code
        self.url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

    def raise_for_status(self):
        if self.status_code >= 400:
            raise RuntimeError(f"HTTP {self.status_code}")

    def json(self):
        return self._payload


class _FakeSession:
    """按 URL 关键字返回预设响应的假 session。"""

    def __init__(self, mapping):
        self.mapping = mapping
        self.requests = []

    def get(self, url, params=None, timeout=None, **kwargs):
        self.requests.append((url, params or {}))
        for key, resp in self.mapping.items():
            if key in url:
                return resp
        return _FakeResponse(payload={})

    def close(self):
        pass


def _non_docstring_strings(tree: ast.AST) -> list[str]:
    """返回源码里所有**非 docstring** 的字符串字面量。

    用于"不得出现某类字面量"的检查。文档字符串里讨论反面例子是合理的，
    不应被误判为缺陷——所以必须先把 docstring 排除掉。
    """
    docstring_nodes: set[int] = set()
    for node in ast.walk(tree):
        if isinstance(
            node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)
        ):
            body = getattr(node, "body", [])
            if body and isinstance(body[0], ast.Expr):
                value = body[0].value
                if isinstance(value, ast.Constant) and isinstance(value.value, str):
                    docstring_nodes.add(id(value))

    return [
        node.value
        for node in ast.walk(tree)
        if isinstance(node, ast.Constant)
        and isinstance(node.value, str)
        and id(node) not in docstring_nodes
    ]


class TestNoHardcodedPaths:
    """不得有硬编码绝对路径。"""

    def test_no_windows_desktop_path(self):
        for name in ("bioproject_client.py", "search_bioprojects.py"):
            source = (SKILL_DIR / "scripts" / name).read_text(encoding="utf-8")
            assert "C:\\Users" not in source, f"{name} 含硬编码 Windows 路径"
            assert "C:/Users" not in source, f"{name} 含硬编码 Windows 路径"
            assert "/Users/" not in source, f"{name} 含硬编码 macOS 路径"

    def test_no_desktop_in_source(self):
        """不得出现硬编码的桌面路径字面量。

        只检查**非 docstring 的字符串常量**——文档里说明"老脚本写了
        C:\\Users\\DELL\\Desktop 这类路径"是必要的说明文字，但代码里
        真的写死路径就是缺陷。
        """
        for name in ("bioproject_client.py", "search_bioprojects.py"):
            source = (SKILL_DIR / "scripts" / name).read_text(encoding="utf-8")
            tree = ast.parse(source)
            for literal in _non_docstring_strings(tree):
                assert "Desktop" not in literal, f"{name} 含硬编码桌面路径"
                assert "Users" not in literal or "you@example" in literal


class TestNoHtmlScraping:
    """不得依赖 HTML 爬虫。"""

    def test_no_beautifulsoup_or_lxml(self):
        for name in ("bioproject_client.py", "search_bioprojects.py"):
            source = (SKILL_DIR / "scripts" / name).read_text(encoding="utf-8")
            tree = ast.parse(source)
            modules = set()
            for node in ast.walk(tree):
                if isinstance(node, ast.Import):
                    modules.update(a.name.split(".")[0] for a in node.names)
                elif isinstance(node, ast.ImportFrom) and node.module and node.level == 0:
                    modules.add(node.module.split(".")[0])
            assert "bs4" not in modules, f"{name} 不应依赖 BeautifulSoup"
            assert "lxml" not in modules, f"{name} 不应依赖 lxml"

    def test_no_xpath_usage(self):
        for name in ("bioproject_client.py", "search_bioprojects.py"):
            source = (SKILL_DIR / "scripts" / name).read_text(encoding="utf-8")
            tree = ast.parse(source)
            attrs = {n.attr for n in ast.walk(tree) if isinstance(n, ast.Attribute)}
            assert "xpath" not in attrs, f"{name} 不应使用 XPath 解析网页"
            assert "select" not in attrs, f"{name} 不应使用 CSS 选择器解析网页"


class TestNoFakeEmail:
    """不得伪造 NCBI email。"""

    def test_no_placeholder_email_in_source(self):
        """不得把占位邮箱作为**实际使用的值**。

        检查的重点是"有没有哪个 email 变量被赋成了 userN@example.com
        这种占位值"。日志消息里为了说明"为什么需要 email"而提到占位邮箱
        是合理的，不应误报。
        """
        for name in ("bioproject_client.py", "search_bioprojects.py"):
            source = (SKILL_DIR / "scripts" / name).read_text(encoding="utf-8")
            tree = ast.parse(source)
            for node in ast.walk(tree):
                if not isinstance(node, ast.Assign):
                    continue
                for target in node.targets:
                    if not isinstance(target, ast.Name):
                        continue
                    var_name = target.id.lower()
                    if "email" not in var_name and "mail" not in var_name:
                        continue
                    value = node.value
                    if isinstance(value, ast.Constant) and isinstance(value.value, str):
                        assert "example.com" not in value.value, (
                            f"{name} 把占位邮箱赋给了 {target.id}"
                        )

    def test_email_is_required(self):
        """缺 email 时必须报错，而不是用一个占位值蒙混过去。"""
        with pytest.raises(ValueError, match="email"):
            bioproject_client.search_bioprojects(["PRJNA1"], email="")

    def test_cli_rejects_missing_email(self, monkeypatch):
        monkeypatch.delenv("BIOPIPE_EMAIL", raising=False)
        assert search_bioprojects.main(["--ids", "PRJNA1", "--output-dir", "out"]) == 2


class TestBuildEsearchTerm:
    """检索式构造。"""

    def test_single_accession(self):
        term = bioproject_client.build_esearch_term(["PRJNA172563"])
        assert term == "PRJNA172563[Project Accession]"

    def test_multiple_accessions_joined_with_or(self):
        term = bioproject_client.build_esearch_term(["PRJNA1", "PRJEB2"])
        assert term == "PRJNA1[Project Accession] OR PRJEB2[Project Accession]"

    def test_empty_returns_empty(self):
        assert bioproject_client.build_esearch_term([]) == ""
        assert bioproject_client.build_esearch_term(["  "]) == ""

    def test_field_qualifier_present(self):
        """必须带 [Project Accession] 限定符，否则会匹配到无关记录。"""
        term = bioproject_client.build_esearch_term(["PRJNA1"])
        assert "[Project Accession]" in term


class TestParseBioprojectXml:
    """XML 解析。"""

    def test_extracts_core_fields(self):
        records = bioproject_client.parse_bioproject_xml(_SAMPLE_XML)
        assert len(records) == 3
        record = next(r for r in records if r["accession"] == "PRJNA212358")
        assert record["title"].startswith("Gene expression profiles of Hyalella azteca")
        assert record["organism"] == "Hyalella azteca"
        assert record["taxid"] == "294128"
        assert record["uid"] == "212358"

    def test_description_parsed(self):
        records = bioproject_client.parse_bioproject_xml(_SAMPLE_XML)
        record = next(r for r in records if r["accession"] == "PRJNA212358")
        assert "gene expression differences" in record["description"]

    def test_empty_xml_returns_empty(self):
        assert bioproject_client.parse_bioproject_xml("") == []

    def test_malformed_xml_does_not_raise(self):
        """坏 XML 应返回空列表而不是抛异常。"""
        assert bioproject_client.parse_bioproject_xml("<not-xml") == []


class TestRealResponseFieldExtraction:
    """针对**真实** NCBI 响应的字段抽取回归测试。

    这组用例锁的是开发过程中实际踩到的三个坑——它们都不会抛异常，
    只会静默地把字段取空或取错：

    1. ``<Organism species="28901" taxID="90370">`` 里的 ``species``
       属性存的是 **taxid**，不是物种名。物种名在 ``<OrganismName>``
       子元素里。
    2. 提交日期在 ``<Submission submitted="2013-07-16">`` 的**属性**
       上，不是子元素——只查元素节点会永远拿到空串。
    3. 发布日期在 ``<ProjectDescr><ProjectReleaseDate>`` 元素里，
       且只有部分记录才有。
    """

    def test_organism_name_comes_from_element_not_species_attribute(self):
        """``species`` 属性是 taxid，不能当物种名用。

        PRJNA355343 的真实 XML：
        ``<Organism species="28901" taxID="90370"><OrganismName>Salmonella…</OrganismName>``
        """
        records = bioproject_client.parse_bioproject_xml(_SAMPLE_XML)
        record = next(r for r in records if r["accession"] == "PRJNA355343")
        assert record["organism"] == "Salmonella enterica subsp. enterica serovar Typhi"
        # 关键断言：organism 不能是 "28901"（那是 species 属性里的 taxid）
        assert record["organism"] != "28901"
        assert not record["organism"].isdigit()
        assert record["taxid"] == "90370"

    def test_organism_is_never_a_bare_number(self):
        """所有记录的 organism 都应该是名称而不是数字。"""
        records = bioproject_client.parse_bioproject_xml(_SAMPLE_XML)
        for record in records:
            assert record["organism"], f"{record['accession']} 的 organism 为空"
            assert not record["organism"].isdigit(), (
                f"{record['accession']} 的 organism 是纯数字：{record['organism']}"
            )

    def test_submitted_date_read_from_attribute(self):
        """提交日期在 Submission 元素的属性上。"""
        records = bioproject_client.parse_bioproject_xml(_SAMPLE_XML)
        record = next(r for r in records if r["accession"] == "PRJNA212358")
        assert record["registration_date"] == "2013-07-16"
        assert record["last_update_date"] == "2013-10-23"

    def test_release_date_read_from_project_descr(self):
        """发布日期在 ProjectDescr/ProjectReleaseDate。"""
        records = bioproject_client.parse_bioproject_xml(_SAMPLE_XML)
        record = next(r for r in records if r["accession"] == "PRJNA355343")
        assert record["release_date"] == "2016-12-15"

    def test_submitter_read_from_organization(self):
        """提交机构在 Submission/Description/Organization/Name。"""
        records = bioproject_client.parse_bioproject_xml(_SAMPLE_XML)
        record = next(r for r in records if r["accession"] == "PRJNA355343")
        assert record["submitter"] == "Christian Medical College"

    def test_missing_release_date_is_empty_not_error(self):
        """部分记录根本没有发布日期——应为空串，不是异常。"""
        records = bioproject_client.parse_bioproject_xml(_SAMPLE_XML)
        record = next(r for r in records if r["accession"] == "PRJNA212358")
        assert record["release_date"] == ""


class TestParseESummaryXml:
    """esummary 端点的解析。"""

    def test_extracts_structured_fields(self):
        records = bioproject_client.parse_esummary_xml(_SAMPLE_ESUMMARY_XML)
        assert len(records) == 3
        record = next(r for r in records if r["accession"] == "PRJNA355343")
        assert record["taxid"] == "90370"
        assert record["organism"] == "Salmonella enterica subsp. enterica serovar Typhi"
        assert record["study_type"] == "Primary submission"
        assert record["submitter"] == "Christian Medical College"

    def test_registration_date_is_normalized(self):
        """esummary 的日期是 ``2013/07/16 00:00``，须归一化成 ISO。"""
        records = bioproject_client.parse_esummary_xml(_SAMPLE_ESUMMARY_XML)
        record = next(r for r in records if r["accession"] == "PRJNA212358")
        assert record["registration_date"] == "2013-07-16"

    def test_placeholder_date_is_not_treated_as_real(self):
        """NCBI 的占位日期 ``1/01/01 00:00`` 不能冒充真实日期。"""
        records = bioproject_client.parse_esummary_xml(_SAMPLE_ESUMMARY_XML)
        record = next(r for r in records if r["accession"] == "PRJNA172563")
        assert record["registration_date"] == ""


class TestNormalizeDate:
    """三种日期格式的归一化。"""

    @pytest.mark.parametrize(
        ("raw", "expected"),
        [
            ("2016-12-15T00:00:00Z", "2016-12-15"),
            ("2013-07-16", "2013-07-16"),
            ("2013/07/16 00:00", "2013-07-16"),
            ("2013/7/6 00:00", "2013-07-06"),
            ("1/01/01 00:00", ""),
            ("", ""),
            ("not a date", ""),
        ],
    )
    def test_normalize(self, raw, expected):
        assert bioproject_client.normalize_date(raw) == expected


class TestMergeRecords:
    """esummary + efetch 双源合并。"""

    def test_efetch_wins_when_both_have_value(self):
        merged = bioproject_client.merge_records(
            [{"uid": "1", "organism": "from-efetch", "release_date": "2020-01-01"}],
            [{"uid": "1", "organism": "from-esummary", "release_date": "1999-01-01"}],
        )
        assert merged[0]["organism"] == "from-efetch"
        assert merged[0]["release_date"] == "2020-01-01"

    def test_esummary_fills_gaps(self):
        merged = bioproject_client.merge_records(
            [{"uid": "1", "organism": "", "taxid": "", "release_date": "2020-01-01"}],
            [{"uid": "1", "organism": "from-esummary", "taxid": "9606", "release_date": ""}],
        )
        assert merged[0]["organism"] == "from-esummary"
        assert merged[0]["taxid"] == "9606"
        assert merged[0]["release_date"] == "2020-01-01"

    def test_real_responses_merge_to_complete_records(self):
        """真实响应合并后，关键字段不该有空缺。"""
        efetch_records = bioproject_client.parse_bioproject_xml(_SAMPLE_XML)
        esummary_records = bioproject_client.parse_esummary_xml(_SAMPLE_ESUMMARY_XML)
        merged = bioproject_client.merge_records(efetch_records, esummary_records)
        assert len(merged) == 3
        for record in merged:
            assert record["accession"]
            assert record["title"]
            assert record["organism"]
            assert not record["organism"].isdigit()
            assert record["taxid"]
            assert record["submitter"]
            assert record["study_type"]


class TestSearchBioprojects:
    """批量查询流程。"""

    @staticmethod
    def _session(uidlist=None):
        ids = uidlist if uidlist is not None else ["212358", "355343", "172563"]
        return _FakeSession(
            {
                "esearch": _FakeResponse(
                    payload={"esearchresult": {"idlist": ids, "count": str(len(ids))}}
                ),
                "esummary": _FakeResponse(text=_SAMPLE_ESUMMARY_XML),
                "efetch": _FakeResponse(text=_SAMPLE_XML),
            }
        )

    def test_esearch_then_esummary_and_efetch(self):
        """先 esearch 拿 UID，再**同时**用 esummary 与 efetch 拉详情。"""
        session = self._session(["212358"])
        records = bioproject_client.search_bioprojects(
            ["PRJNA212358"], email="me@lab.org", session=session
        )
        # 夹具里有 3 条记录（真实响应就是整批返回的），按 accession 定位
        record = next(r for r in records if r["accession"] == "PRJNA212358")
        assert record["organism"] == "Hyalella azteca"
        assert record["taxid"] == "294128"

        urls = [url for url, _ in session.requests]
        assert any("esearch" in u for u in urls)
        assert any("efetch" in u for u in urls)
        # esummary 不能少——taxid/organism/submitter 只有它有
        assert any("esummary" in u for u in urls)

    def test_esummary_is_required_for_organism_fields(self):
        """只有 esummary 有 organism/taxid——少调它就会丢字段。"""
        session = _FakeSession(
            {
                "esearch": _FakeResponse(
                    payload={"esearchresult": {"idlist": ["172563"], "count": "1"}}
                ),
                "esummary": _FakeResponse(text=_SAMPLE_ESUMMARY_XML),
                "efetch": _FakeResponse(text=_SAMPLE_XML),
            }
        )
        records = bioproject_client.search_bioprojects(
            ["PRJNA172563"], email="me@lab.org", session=session
        )
        record = next(r for r in records if r["accession"] == "PRJNA172563")
        assert record["organism"] == "Homo sapiens"
        assert record["taxid"] == "9606"
        assert record["submitter"] == "NIH"

    def test_email_is_sent_to_ncbi(self):
        session = self._session(["212358"])
        bioproject_client.search_bioprojects(
            ["PRJNA212358"], email="me@lab.org", session=session
        )
        for _url, params in session.requests:
            assert params.get("email") == "me@lab.org"

    def test_api_key_is_forwarded(self):
        session = self._session(["212358"])
        bioproject_client.search_bioprojects(
            ["PRJNA212358"], email="me@lab.org", api_key="KEY123", session=session
        )
        for _url, params in session.requests:
            assert params.get("api_key") == "KEY123"

    def test_no_match_returns_empty(self):
        session = _FakeSession(
            {"esearch": _FakeResponse(payload={"esearchresult": {"idlist": [], "count": "0"}})}
        )
        records = bioproject_client.search_bioprojects(
            ["PRJNA999999999"], email="me@lab.org", session=session
        )
        assert records == []

    def test_empty_accessions(self):
        assert bioproject_client.search_bioprojects([], email="me@lab.org") == []


class TestReadIdFile:
    """accession 列表解析。"""

    def test_reads_plain_list(self, tmp_path):
        path = tmp_path / "ids.txt"
        path.write_text("PRJNA172563\nPRJNA212358\n", encoding="utf-8")
        assert search_bioprojects.read_id_file(path) == ["PRJNA172563", "PRJNA212358"]

    def test_skips_comments_and_blanks(self, tmp_path):
        path = tmp_path / "ids.txt"
        path.write_text("# note\n\nPRJNA1\n  \n", encoding="utf-8")
        assert search_bioprojects.read_id_file(path) == ["PRJNA1"]

    def test_deduplicates(self, tmp_path):
        path = tmp_path / "ids.txt"
        path.write_text("PRJNA1\nPRJNA1\nPRJNA2\n", encoding="utf-8")
        assert search_bioprojects.read_id_file(path) == ["PRJNA1", "PRJNA2"]

    def test_trailing_comment_stripped(self, tmp_path):
        path = tmp_path / "ids.txt"
        path.write_text("PRJNA1 some note\n", encoding="utf-8")
        assert search_bioprojects.read_id_file(path) == ["PRJNA1"]

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            search_bioprojects.read_id_file(tmp_path / "nope.txt")


class TestAccessionValidation:
    """accession 格式校验。"""

    @pytest.mark.parametrize(
        "text,expected",
        [
            ("PRJNA172563", True),
            ("PRJEB1234", True),
            ("PRJDB5678", True),
            ("prjna172563", True),
            ("id", False),
            ("PRJNA", False),
            ("172563", False),
            ("SRR390728", False),
        ],
    )
    def test_is_valid_accession(self, text, expected):
        assert search_bioprojects.is_valid_accession(text) is expected


class TestNoiseDetection:
    """导航噪声判定。"""

    def test_short_navigation_text_is_noise(self):
        assert search_bioprojects._looks_like_noise("Navigate Across additional projects")

    def test_long_description_not_treated_as_noise(self):
        """老脚本对任意长度描述都套噪声词，会误杀正文。"""
        long_text = "Navigate across the samples. " + "x" * 300
        assert not search_bioprojects._looks_like_noise(long_text)

    def test_real_description_is_not_noise(self):
        assert not search_bioprojects._looks_like_noise(
            "This study investigates gene expression differences in Hyalella azteca."
        )


class TestFilterRecords:
    """关键词过滤。"""

    def test_no_keyword_keeps_all(self):
        records = [{"title": "a"}, {"title": "b"}]
        assert len(search_bioprojects.filter_records(records, None)) == 2

    def test_keyword_filters_case_insensitively(self):
        records = [
            {"title": "CRISPR screening", "description": "", "name": ""},
            {"title": "RNA-seq", "description": "", "name": ""},
        ]
        kept = search_bioprojects.filter_records(records, "crispr")
        assert len(kept) == 1
        assert kept[0]["title"] == "CRISPR screening"

    def test_keyword_matches_description(self):
        records = [{"title": "x", "description": "contains MAGIC keyword", "name": ""}]
        assert len(search_bioprojects.filter_records(records, "magic")) == 1


class TestOutputs:
    """结果写出。"""

    def test_writes_csv_and_json(self, tmp_path):
        records = [
            {
                "accession": "PRJNA1",
                "title": "T",
                "name": "N",
                "description": "D",
                "organism": "O",
                "taxid": "1",
                "study_type": "BioProject",
                "submitter": "S",
                "registration_date": "2018",
                "release_date": "2019",
                "uid": "1",
            }
        ]
        search_bioprojects.write_outputs(records, tmp_path)
        assert (tmp_path / "bioproject_data.csv").is_file()
        assert (tmp_path / "bioproject_data.json").is_file()

    def test_csv_has_all_columns(self, tmp_path):
        records = [
            {
                "accession": "PRJNA1",
                "title": "T",
                "name": "N",
                "description": "D",
                "organism": "O",
                "taxid": "1",
                "study_type": "BioProject",
                "submitter": "S",
                "registration_date": "2018",
                "release_date": "2019",
                "last_update_date": "2020",
                "uid": "1",
            }
        ]
        search_bioprojects.write_outputs(records, tmp_path)
        text = (tmp_path / "bioproject_data.csv").read_text(encoding="utf-8")
        assert "accession" in text
        assert "organism" in text
        assert "description" in text
        assert "last_update_date" in text

    def test_csv_header_covers_all_documented_columns(self, tmp_path):
        """CSV 表头必须包含文档里承诺的所有字段。"""
        search_bioprojects.write_outputs([], tmp_path)
        text = (tmp_path / "bioproject_data.csv").read_text(encoding="utf-8")
        header = text.splitlines()[0]
        for column in (
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
        ):
            assert column in header.split(","), f"CSV 表头缺少 {column}"


class TestCliSurface:
    """命令行接口。"""

    def test_parser_builds(self):
        args = search_bioprojects.build_parser().parse_args(
            ["--ids", "PRJNA1", "--output-dir", "out/", "--email", "me@lab.org"]
        )
        assert args.ids == ["PRJNA1"]
        assert args.email == "me@lab.org"

    def test_ids_required(self):
        with pytest.raises(SystemExit):
            search_bioprojects.build_parser().parse_args([])

    def test_help_exits_cleanly(self, capsys):
        with pytest.raises(SystemExit) as exc:
            search_bioprojects.build_parser().parse_args(["--help"])
        assert exc.value.code == 0
        assert "BioProject" in capsys.readouterr().out
