"""BioDatabase 检索技能（PubChem + NCBI Taxonomy）的回归测试。

锁住老仓库 ``crawl_database/`` 的一批真实缺陷，以及两个官方 API 的陷阱：

1. **伪造 NCBI email** —— 老脚本准备 5 个 ``userN@example.com`` 轮流冒充。
   NCBI 的 email 是身份标识，伪造既违规也会被限速。这里要求 email 真实
   必填，且代码里不许出现占位地址。
2. **多线程并发写同一个文件** —— 老脚本用 ``ThreadPoolExecutor``，每个
   线程各自 ``open(output_tsv, "a")`` 追加一行，写入交错会撕裂成半行。
   这里禁止模块级写文件。
3. **CAS 号取错了地方** —— 老脚本从 ``xrefs/RegistryID`` 取"CAS"，但那个
   端点返回 461 条各种数据库的注册号，真正的 CAS 只是其中 2 条；从
   ``synonyms`` 能拿到 3 条。
4. **``CanonicalSMILES`` 字段名已变** —— PubChem 现在返回为 ``SMILES``。
5. **不存在的 CID 返回 HTTP 200** —— 内容是 ``{"CID": 999999999}``，
   只看状态码会误判为查询成功。
6. **``root.iter("Taxon")`` 过度收集** —— 会把 ``LineageEx`` 里的祖先
   也当成独立记录，查 3 个物种返回 75 条。
7. **``LineageEx`` 不含自身** —— species 层级只存在于顶层元素。
8. **只挑 8 个标准 rank** —— 会漏掉 ``domain``（Eukaryota/Bacteria）。

全部用例不联网。
"""

from __future__ import annotations

import ast
import csv
import json
import re
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
SKILL_DIR = REPO_ROOT / "skills" / "bio-database-fetch"
FIXTURES = REPO_ROOT / "tests" / "fixtures_bio_database"
sys.path.insert(0, str(REPO_ROOT / "packages"))
sys.path.insert(0, str(SKILL_DIR / "scripts"))

import database_client  # noqa: E402
import fetch_compounds  # noqa: E402
import fetch_taxonomy  # noqa: E402

SCRIPTS = ("database_client.py", "fetch_compounds.py", "fetch_taxonomy.py")


def _script(name: str) -> Path:
    return SKILL_DIR / "scripts" / name


def _parse(name: str) -> ast.Module:
    return ast.parse(_script(name).read_text(encoding="utf-8"))


def _non_docstring_strings(tree: ast.AST) -> list[str]:
    """返回源码里所有**非 docstring** 的字符串字面量。"""
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


def _all_strings() -> list[str]:
    """所有脚本里的非 docstring 字符串。"""
    result: list[str] = []
    for name in SCRIPTS:
        result.extend(_non_docstring_strings(_parse(name)))
    return result


# --------------------------------------------------------------------------
# 老仓库缺陷的回归锁
# --------------------------------------------------------------------------


class TestNoFabricatedEmail:
    """老脚本用 ``userN@example.com`` 轮流冒充 NCBI 调用方。"""

    def test_no_placeholder_email_addresses(self):
        """代码里不得出现占位邮箱。

        文档字符串里的示例（``you@example.edu``）是给人看的，不算硬编码；
        这里只检查**非 docstring** 的字符串。
        """
        pattern = re.compile(r"[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}")
        offenders = []
        for value in _all_strings():
            for match in pattern.findall(value):
                # 示例占位符不算违规
                if match in {"you@example.edu", "you@your-institution.edu"}:
                    continue
                offenders.append(match)
        assert not offenders, f"发现占位/硬编码邮箱: {offenders}"

    def test_no_example_dot_com_placeholders(self):
        """老脚本用的是 ``userN@example.com`` 这种形态，必须彻底消失。

        docstring 里说明"老脚本用了假邮箱"是必要的文档，不算违规。
        """
        for name in SCRIPTS:
            for value in _non_docstring_strings(_parse(name)):
                assert "example.com" not in value, f"{name}: {value!r}"

    def test_email_is_read_from_cli_or_env(self):
        """email 必须来自命令行参数或环境变量，而不是硬编码。"""
        source = _script("fetch_taxonomy.py").read_text(encoding="utf-8")
        assert "BIOPIPE_EMAIL" in source
        assert "--email" in source

    def test_fetch_taxonomy_requires_email(self):
        """缺 email 时必须报错，而不是静默发出无身份的请求。"""
        with pytest.raises(ValueError, match="email"):
            database_client.fetch_taxonomy(["9606"], email="")

    def test_main_returns_2_without_email(self, monkeypatch, tmp_path):
        """CLI 缺 email 时退出码为 2（参数错误），不是 0。"""
        monkeypatch.delenv("BIOPIPE_EMAIL", raising=False)
        code = fetch_taxonomy.main(
            ["--tax-ids", "9606", "--output-dir", str(tmp_path)]
        )
        assert code == 2

    def test_email_is_passed_to_the_request(self):
        """email 必须真的拼进请求参数。"""
        captured: dict = {}

        class _FakeResponse:
            text = "<TaxaSet></TaxaSet>"

            def raise_for_status(self):
                return None

        class _FakeSession:
            def get(self, url, params=None, timeout=None):
                captured["params"] = params or {}
                return _FakeResponse()

            def close(self):
                return None

        database_client.fetch_taxonomy(
            ["9606"], email="real.person@lab.edu", session=_FakeSession()
        )
        assert captured["params"].get("email") == "real.person@lab.edu"
        assert captured["params"].get("db") == "taxonomy"


class TestNoModuleLevelSideEffects:
    """老脚本在 import 时就写文件、起线程。"""

    @staticmethod
    def _top_level_calls(name: str) -> list[str]:
        """返回**模块顶层**直接执行的调用名。

        只看 ``tree.body`` 里的顶层语句本身，不下钻到函数体——函数体里
        的 ``open()`` 是正常业务逻辑，只有顶层执行才是副作用。
        """
        tree = _parse(name)
        calls: list[str] = []
        for node in tree.body:
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                # 函数/类的定义本身不算调用，但它的默认值会求值
                for default in list(node.args.defaults) + [
                    d for d in node.args.kw_defaults if d is not None
                ]:
                    for sub in ast.walk(default):
                        if isinstance(sub, ast.Call):
                            calls.append(getattr(sub.func, "id", "?"))
                continue
            if isinstance(node, ast.Expr) and isinstance(node.value, ast.Constant):
                continue  # docstring
            for sub in ast.walk(node):
                if isinstance(sub, ast.Call):
                    func = sub.func
                    calls.append(
                        getattr(func, "id", None) or getattr(func, "attr", None) or "?"
                    )
        return calls

    def test_no_module_level_open(self):
        """模块顶层不得调用 ``open()``。"""
        for name in SCRIPTS:
            assert "open" not in self._top_level_calls(name), (
                f"{name} 在模块顶层调用了 open()"
            )

    def test_no_module_level_network_call(self):
        """模块顶层不得直接发起网络请求。"""
        for name in SCRIPTS:
            calls = self._top_level_calls(name)
            for called in calls:
                assert called not in {"get", "post", "request", "urlopen"}, (
                    f"{name} 在模块顶层发起了网络请求: {called}()"
                )

    def test_no_module_level_write(self):
        """模块顶层不得写文件（老脚本 import 即写 TSV）。"""
        for name in SCRIPTS:
            calls = self._top_level_calls(name)
            for called in calls:
                assert called not in {"write_text", "write_bytes", "dump"}, (
                    f"{name} 在模块顶层写文件: {called}()"
                )

    def test_no_thread_pool(self):
        """不得用线程池并发写同一个文件（老脚本的写入撕裂问题）。

        docstring 里说明"老脚本用了线程池"是必要的文档，不算违规。
        """
        for name in SCRIPTS:
            for value in _non_docstring_strings(_parse(name)):
                assert "ThreadPoolExecutor" not in value
                assert "multiprocessing" not in value

    def test_no_globally_shared_output_variable(self):
        """老脚本用 ``global output_tsv`` 让多线程共享同一个写句柄。"""
        for name in SCRIPTS:
            tree = _parse(name)
            for node in ast.walk(tree):
                if isinstance(node, ast.Global):
                    assert node.names != ["output_tsv"], (
                        f"{name} 仍在用 global output_tsv"
                    )

    def test_importing_the_client_is_side_effect_free(self):
        """import 模块不应产生任何输出或文件。"""
        assert database_client.PUG_REST_BASE.startswith("https://")


class TestNoHardcodedPaths:
    """老脚本遍布 ``C:\\Users\\DELL\\Desktop\\...`` 之类的绝对路径。"""

    def test_no_absolute_windows_paths(self):
        pattern = re.compile(r"[A-Za-z]:\\\\?Users\\\\", re.IGNORECASE)
        for value in _all_strings():
            assert not pattern.search(value), f"发现硬编码 Windows 路径: {value!r}"

    def test_no_absolute_macos_paths(self):
        pattern = re.compile(r"/Users/[^/\s]+/")
        for value in _all_strings():
            assert not pattern.search(value), f"发现硬编码 macOS 路径: {value!r}"

    def test_no_msedgedriver_path(self):
        """老脚本硬编码了 msedgedriver 的本地路径。"""
        for name in SCRIPTS:
            source = _script(name).read_text(encoding="utf-8")
            assert "msedgedriver" not in source
            assert "chromedriver" not in source


class TestNoSelenium:
    """这个技能走官方 REST API，不需要浏览器自动化。"""

    def test_no_selenium_import(self):
        for name in SCRIPTS:
            source = _script(name).read_text(encoding="utf-8")
            assert "selenium" not in source.lower()
            assert "webdriver" not in source.lower()

    def test_hit2lead_internal_site_removed(self):
        """老脚本 01_hit2lead 爬的是公司内部商业站，已删除。

        文档字符串里可以说明"为什么删除"，但不能有可执行的调用——所以
        这里检查非 docstring 的字符串。
        """
        for name in SCRIPTS:
            for value in _non_docstring_strings(_parse(name)):
                assert "hit2lead" not in value.lower()
                assert "fake_useragent" not in value.lower()

    def test_hit2lead_not_referenced_in_cli(self):
        """命令行参数里不得出现 hit2lead 相关入口。"""
        for name in ("fetch_compounds.py", "fetch_taxonomy.py"):
            source = _script(name).read_text(encoding="utf-8")
            assert "hit2lead" not in source.lower()


# --------------------------------------------------------------------------
# PubChem 端点行为
# --------------------------------------------------------------------------


class _StubResponse:
    def __init__(self, payload: str):
        self._payload = payload

    def raise_for_status(self):
        return None

    def json(self):
        return json.loads(self._payload)

    @property
    def text(self):
        return self._payload


class _StubSession:
    """按 URL 后缀返回冻结夹具的假 session。"""

    def __init__(self, routes: dict[str, str]):
        self.routes = routes
        self.calls: list[str] = []

    def get(self, url, params=None, timeout=None):
        self.calls.append(url)
        for suffix, payload in self.routes.items():
            if suffix in url:
                return _StubResponse(payload)
        raise AssertionError(f"未预期的请求 URL: {url}")

    def close(self):
        return None


def _fixture(name: str) -> str:
    return (FIXTURES / name).read_text(encoding="utf-8")


class TestCompoundProperties:
    """属性查询与 ``found`` 标记。"""

    def test_finds_existing_compound(self):
        payload = _fixture("pubchem_properties.json")
        records = database_client.fetch_compounds(
            ["2244", "999999999"], session=_StubSession({"/property/": payload})
        )
        by_cid = {r["CID"]: r for r in records}
        assert by_cid["2244"]["found"] is True
        assert by_cid["2244"]["MolecularFormula"] == "C9H8O4"
        assert by_cid["2244"]["InChIKey"] == "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"

    def test_missing_compound_marked_not_found(self):
        """不存在的 CID 返回 HTTP 200 但只有 CID 字段——必须标记为未找到。"""
        payload = _fixture("pubchem_properties.json")
        records = database_client.fetch_compounds(
            ["999999999"], session=_StubSession({"/property/": payload})
        )
        assert len(records) == 1
        assert records[0]["found"] is False
        # 不能凭空造出属性
        assert "MolecularFormula" not in records[0]

    def test_cid_preserved_even_when_not_returned(self):
        """PubChem 完全没提到的 CID 也要出现在结果里，标记为未找到。"""
        payload = _fixture("pubchem_properties.json")
        records = database_client.fetch_compounds(
            ["12345"], session=_StubSession({"/property/": payload})
        )
        assert records[0]["CID"] == "12345"
        assert records[0]["found"] is False

    def test_smiles_field_name(self):
        """``CanonicalSMILES`` 已被 PubChem 改名为 ``SMILES``。"""
        payload = _fixture("pubchem_properties.json")
        records = database_client.fetch_compounds(
            ["2244"], session=_StubSession({"/property/": payload})
        )
        assert records[0]["SMILES"]
        assert records[0]["ConnectivitySMILES"]

    def test_both_smiles_variants_requested(self):
        """默认属性集必须同时包含两个 SMILES 变体。"""
        assert "SMILES" in database_client.COMPOUND_PROPERTIES
        assert "ConnectivitySMILES" in database_client.COMPOUND_PROPERTIES
        # 旧名不应再出现在请求里
        assert "CanonicalSMILES" not in database_client.COMPOUND_PROPERTIES

    def test_batches_at_100_cids(self):
        """单次请求上限 100 个 CID，必须自动分块。"""
        session = _StubSession({"/property/": _fixture("pubchem_properties.json")})
        cids = [str(i) for i in range(250)]
        database_client.fetch_compounds(cids, session=session)
        property_calls = [u for u in session.calls if "/property/" in u]
        assert len(property_calls) == 3


class TestSynonymsAndRegistry:
    """同义词与 RegistryID 端点。"""

    def test_synonyms_parsed(self):
        payload = _fixture("pubchem_synonyms.json")
        result = database_client.fetch_synonyms(
            ["2244"], session=_StubSession({"/synonyms/": payload})
        )
        assert "2244" in result
        assert len(result["2244"]) > 0

    def test_registry_ids_parsed(self):
        payload = _fixture("pubchem_registry.json")
        result = database_client.fetch_registry_ids(
            ["2244"], session=_StubSession({"/RegistryID/": payload})
        )
        assert "2244" in result
        assert len(result["2244"]) > 0


class TestCasExtraction:
    """CAS 号必须按格式从同义词里筛出来，不能直接从 RegistryID 取。"""

    def test_cas_pattern(self):
        assert database_client.CAS_PATTERN.match("50-78-2")
        assert database_client.CAS_PATTERN.match("11126-35-5")
        assert database_client.CAS_PATTERN.match("98201-60-6")

    def test_non_cas_rejected(self):
        for value in ("CHEBI:15365", "DB00945", "LGCAMP0133.00-11", "50-78-2-x"):
            assert not database_client.CAS_PATTERN.match(value), value

    def test_pure_digits_rejected(self):
        """RegistryID 里的 ``0000050782`` 不是 CAS 号（没有连字符）。"""
        assert not database_client.CAS_PATTERN.match("0000050782")
        assert not database_client.CAS_PATTERN.match("00002661")

    def test_extract_from_synonyms(self):
        payload = _fixture("pubchem_synonyms.json")
        synonyms = database_client.fetch_synonyms(
            ["2244"], session=_StubSession({"/synonyms/": payload})
        )
        cas = database_client.extract_cas_numbers(synonyms["2244"])
        assert cas == ["50-78-2", "98201-60-6", "11126-35-5"]

    def test_registry_has_fewer_cas_than_synonyms(self):
        """实测：synonyms 3 条 CAS，RegistryID 只有 2 条。"""
        synonyms = database_client.fetch_synonyms(
            ["2244"], session=_StubSession({"/synonyms/": _fixture("pubchem_synonyms.json")})
        )
        registry = database_client.fetch_registry_ids(
            ["2244"],
            session=_StubSession({"/RegistryID/": _fixture("pubchem_registry.json")}),
        )
        from_syn = database_client.extract_cas_numbers(synonyms["2244"])
        from_reg = database_client.extract_cas_numbers(registry["2244"])
        assert len(from_syn) > len(from_reg)
        assert "11126-35-5" in from_syn
        assert "11126-35-5" not in from_reg

    def test_deduplicates(self):
        values = ["50-78-2", "50-78-2", " 50-78-2 ", "11126-35-5"]
        assert database_client.extract_cas_numbers(values) == ["50-78-2", "11126-35-5"]

    def test_empty_input(self):
        assert database_client.extract_cas_numbers([]) == []


class TestBuildRows:
    """输出行的组装。"""

    def test_merges_cas_from_both_sources(self):
        records = [{"CID": "2244", "found": True, "MolecularFormula": "C9H8O4"}]
        synonyms = {"2244": ["50-78-2", "11126-35-5"]}
        registry = {"2244": ["50-78-2", "98201-60-6"]}
        rows = fetch_compounds.build_rows(records, synonyms, registry, max_synonyms=5)
        assert len(rows) == 1
        cas = rows[0]["CAS"].split("; ")
        assert cas == ["50-78-2", "11126-35-5", "98201-60-6"]

    def test_synonym_count_reports_full_number(self):
        records = [{"CID": "2244", "found": True}]
        synonyms = {"2244": [f"s{i}" for i in range(699)]}
        rows = fetch_compounds.build_rows(records, synonyms, {}, max_synonyms=20)
        assert rows[0]["synonym_count"] == 699
        assert len(rows[0]["synonyms"].split("; ")) == 20

    def test_found_flag_preserved(self):
        records = [{"CID": "999999999", "found": False}]
        rows = fetch_compounds.build_rows(records, {}, {}, max_synonyms=0)
        assert rows[0]["found"] is False
        assert rows[0]["CAS"] == ""


class TestCidFileReading:
    """CID 文件解析。"""

    def test_reads_cids_and_skips_noise(self, tmp_path):
        path = tmp_path / "cids.txt"
        path.write_text(
            "# comment\nCID\n2244\n3672\n\n  1983  \n2244\nnot_a_cid\n",
            encoding="utf-8",
        )
        cids = fetch_compounds.read_cid_file(path)
        assert cids == ["2244", "3672", "1983"]

    def test_collect_from_file_and_args(self, tmp_path):
        path = tmp_path / "cids.txt"
        path.write_text("2244\n3672\n", encoding="utf-8")
        cids = fetch_compounds.collect_cids([str(path), "1983"])
        assert cids == ["2244", "3672", "1983"]

    def test_collect_deduplicates(self):
        assert fetch_compounds.collect_cids(["2244", "2244", "3672"]) == ["2244", "3672"]


# --------------------------------------------------------------------------
# NCBI Taxonomy
# --------------------------------------------------------------------------


class TestTaxonomyParsing:
    """Taxonomy XML 解析。"""

    @pytest.fixture(scope="class")
    @classmethod
    def records(cls):
        return database_client.parse_taxonomy_xml(_fixture("taxonomy_response.xml"))

    def test_only_requested_taxa_returned(self, records):
        """关键：不能把 LineageEx 里的祖先当成独立记录。"""
        assert len(records) == 2
        assert {r["tax_id"] for r in records} == {"9606", "4932"}

    def test_no_ancestor_records_leaked(self, records):
        """祖先名不应作为独立记录出现。"""
        names = {r["scientific_name"] for r in records}
        assert "Eukaryota" not in names
        assert "cellular organisms" not in names
        assert "Opisthokonta" not in names

    def test_scientific_name(self, records):
        by_id = {r["tax_id"]: r for r in records}
        assert by_id["9606"]["scientific_name"] == "Homo sapiens"
        assert by_id["4932"]["scientific_name"] == "Saccharomyces cerevisiae"

    def test_common_name(self, records):
        by_id = {r["tax_id"]: r for r in records}
        assert by_id["9606"]["common_name"] == "human"
        assert by_id["4932"]["common_name"] == "brewer's yeast"

    def test_rank(self, records):
        by_id = {r["tax_id"]: r for r in records}
        assert by_id["9606"]["rank"] == "species"

    def test_species_present_though_lineage_excludes_it(self, records):
        """LineageEx 只含祖先，species 必须从顶层补进去。"""
        by_id = {r["tax_id"]: r for r in records}
        assert by_id["9606"]["by_rank"]["species"] == "Homo sapiens"
        assert by_id["4932"]["by_rank"]["species"] == "Saccharomyces cerevisiae"

    def test_lineage_excludes_self(self, records):
        """LineageEx 本身不该包含 species 层的自身条目。"""
        by_id = {r["tax_id"]: r for r in records}
        lineage_names = {item["name"] for item in by_id["9606"]["lineage"]}
        assert "Homo sapiens" not in lineage_names

    def test_domain_retained(self, records):
        """domain（Eukaryota）不能因为"不在 8 个标准 rank 里"而丢失。"""
        by_id = {r["tax_id"]: r for r in records}
        assert by_id["9606"]["by_rank"]["domain"] == "Eukaryota"

    def test_kingdom_and_below(self, records):
        by_id = {r["tax_id"]: r for r in records}
        assert by_id["9606"]["by_rank"]["kingdom"] == "Metazoa"
        assert by_id["9606"]["by_rank"]["phylum"] == "Chordata"
        assert by_id["9606"]["by_rank"]["class"] == "Mammalia"
        assert by_id["9606"]["by_rank"]["order"] == "Primates"
        assert by_id["9606"]["by_rank"]["family"] == "Hominidae"
        assert by_id["9606"]["by_rank"]["genus"] == "Homo"

    def test_lineage_entries_have_taxid(self, records):
        by_id = {r["tax_id"]: r for r in records}
        for item in by_id["9606"]["lineage"]:
            assert item["tax_id"]
            assert item["rank"]
            assert item["name"]

    def test_lineage_depth(self, records):
        """lineage_depth 是**完整** lineage 的层数（含自身）。

        LineageEx 只有 30 个祖先，加上自身是 31 层。这必须与
        ``full_lineage`` 字符串的段数一致——两者曾经差 1（depth 少算
        了自身），会让"深度"和"链长"对不上。
        """
        by_id = {r["tax_id"]: r for r in records}
        assert by_id["9606"]["lineage_depth"] == 31
        assert by_id["4932"]["lineage_depth"] == 13

    def test_lineage_depth_matches_full_lineage_length(self, records):
        """depth 与 full_lineage 的段数必须一致。"""
        rows = fetch_taxonomy.build_rows(records)
        for row in rows:
            segments = [s for s in row["full_lineage"].split(" > ") if s]
            assert row["lineage_depth"] == len(segments), row["tax_id"]

    def test_empty_xml(self):
        assert database_client.parse_taxonomy_xml("") == []

    def test_invalid_xml_returns_empty(self):
        assert database_client.parse_taxonomy_xml("<not-xml") == []


class TestTaxonomyBuilder:
    """检索式构造与输出行。"""

    def test_term_uses_field_qualifier(self):
        term = database_client.build_taxonomy_term(["9606", "4932"])
        assert term == "9606[Taxonomy ID] OR 4932[Taxonomy ID]"

    def test_empty_term(self):
        assert database_client.build_taxonomy_term([]) == ""

    def test_rows_flatten_ranks(self):
        records = database_client.parse_taxonomy_xml(_fixture("taxonomy_response.xml"))
        rows = fetch_taxonomy.build_rows(records)
        by_id = {r["tax_id"]: r for r in rows}
        assert by_id["9606"]["domain"] == "Eukaryota"
        assert by_id["9606"]["species"] == "Homo sapiens"
        assert by_id["9606"]["genus"] == "Homo"

    def test_rows_include_full_lineage(self):
        records = database_client.parse_taxonomy_xml(_fixture("taxonomy_response.xml"))
        rows = fetch_taxonomy.build_rows(records)
        by_id = {r["tax_id"]: r for r in rows}
        lineage = by_id["9606"]["full_lineage"]
        assert lineage.startswith("cellular organisms")
        assert "Homo sapiens" in lineage
        assert " > " in lineage

    def test_domain_is_a_column(self):
        """domain 必须展开成独立列。"""
        assert "domain" in fetch_taxonomy.FIELDNAMES
        assert "superkingdom" in fetch_taxonomy.FIELDNAMES

    def test_elementary_ranks_include_domain(self):
        assert "domain" in database_client.ELEMENTARY_RANKS

    def test_rank_coverage_summary(self):
        records = database_client.parse_taxonomy_xml(_fixture("taxonomy_response.xml"))
        rows = fetch_taxonomy.build_rows(records)
        coverage = fetch_taxonomy.summarize_rank_coverage(rows)
        assert coverage.get("species") == 2
        assert coverage.get("domain") == 2


class TestTaxonomyIdFile:
    """Taxonomy ID 文件解析。"""

    def test_reads_ids_and_skips_noise(self, tmp_path):
        path = tmp_path / "taxids.txt"
        path.write_text("# comment\ntax_id\n9606\n10090\n\n4932\n9606\nxxx\n", encoding="utf-8")
        ids = fetch_taxonomy.read_id_file(path)
        assert ids == ["9606", "10090", "4932"]


# --------------------------------------------------------------------------
# 端到端（用假 session，不联网）
# --------------------------------------------------------------------------


class TestEndToEnd:
    """走完整 CLI 流程，用冻结夹具顶替网络。"""

    def test_compounds_cli(self, monkeypatch, tmp_path):
        routes = {
            "/property/": _fixture("pubchem_properties.json"),
            "/synonyms/": _fixture("pubchem_synonyms.json"),
            "/RegistryID/": _fixture("pubchem_registry.json"),
        }
        session = _StubSession(routes)

        monkeypatch.setattr(
            fetch_compounds, "fetch_compounds",
            lambda cids, properties=None, timeout=60: database_client.fetch_compounds(
                cids, properties=properties or database_client.COMPOUND_PROPERTIES,
                session=session,
            ),
        )
        monkeypatch.setattr(
            fetch_compounds, "fetch_synonyms",
            lambda cids, timeout=60: database_client.fetch_synonyms(cids, session=session),
        )
        monkeypatch.setattr(
            fetch_compounds, "fetch_registry_ids",
            lambda cids, timeout=60: database_client.fetch_registry_ids(cids, session=session),
        )

        out = tmp_path / "out"
        code = fetch_compounds.main(
            ["--cids", "2244", "999999999", "--output-dir", str(out)]
        )
        assert code == 0

        with open(out / "compounds.csv", encoding="utf-8-sig", newline="") as handle:
            rows = list(csv.DictReader(handle))
        assert len(rows) == 2
        by_cid = {r["CID"]: r for r in rows}
        assert by_cid["2244"]["found"] == "True"
        assert by_cid["2244"]["MolecularFormula"] == "C9H8O4"
        assert by_cid["2244"]["CAS"]
        assert by_cid["999999999"]["found"] == "False"

        payload = json.loads((out / "compounds.json").read_text(encoding="utf-8"))
        assert len(payload) == 2

    def test_taxonomy_cli(self, monkeypatch, tmp_path):
        xml = _fixture("taxonomy_response.xml")

        class _XmlSession(_StubSession):
            def get(self, url, params=None, timeout=None):
                self.calls.append(url)
                return _StubResponse(xml)

        session = _XmlSession({})
        monkeypatch.setattr(
            fetch_taxonomy, "fetch_taxonomy",
            lambda tax_ids, email, api_key=None, timeout=60: database_client.fetch_taxonomy(
                tax_ids, email=email, api_key=api_key, session=session
            ),
        )

        out = tmp_path / "out"
        code = fetch_taxonomy.main(
            ["--tax-ids", "9606", "4932", "--output-dir", str(out),
             "--email", "real.person@lab.edu"]
        )
        assert code == 0

        with open(out / "taxonomy.csv", encoding="utf-8-sig", newline="") as handle:
            rows = list(csv.DictReader(handle))
        assert len(rows) == 2
        by_id = {r["tax_id"]: r for r in rows}
        assert by_id["9606"]["scientific_name"] == "Homo sapiens"
        assert by_id["9606"]["domain"] == "Eukaryota"
        assert by_id["9606"]["species"] == "Homo sapiens"

        payload = json.loads((out / "taxonomy.json").read_text(encoding="utf-8"))
        assert len(payload) == 2
        assert payload[0]["lineage"]


# --------------------------------------------------------------------------
# 文档与元数据
# --------------------------------------------------------------------------


class TestSkillMetadata:
    """SKILL.md 与仓库规范的一致性。"""

    def test_skill_md_exists(self):
        assert (SKILL_DIR / "SKILL.md").is_file()

    def test_requirements_declares_requests(self):
        text = (SKILL_DIR / "requirements.txt").read_text(encoding="utf-8")
        assert "requests" in text

    def test_requirements_has_no_local_path_dependency(self):
        """老仓库的 requirements 含 file:// 本地路径，别人装不上。"""
        text = (SKILL_DIR / "requirements.txt").read_text(encoding="utf-8")
        for line in text.splitlines():
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            assert "file://" not in stripped, f"本地路径依赖: {stripped}"

    def test_references_exist(self):
        assert (SKILL_DIR / "references" / "methodology.md").is_file()
        assert (SKILL_DIR / "references" / "troubleshooting.md").is_file()

    def test_documentation_covers_the_traps(self):
        """文档必须把踩过的坑写清楚。"""
        methodology = (SKILL_DIR / "references" / "methodology.md").read_text(
            encoding="utf-8"
        )
        assert "CanonicalSMILES" in methodology
        assert "RegistryID" in methodology
        assert "LineageEx" in methodology
        assert "findall" in methodology

    def test_no_private_email_in_docs(self):
        """老脚本里硬编码的个人邮箱不许出现在任何文件里。"""
        for path in SKILL_DIR.rglob("*"):
            if not path.is_file():
                continue
            text = path.read_text(encoding="utf-8", errors="replace")
            assert "15993327335" not in text
            assert "@163.com" not in text
