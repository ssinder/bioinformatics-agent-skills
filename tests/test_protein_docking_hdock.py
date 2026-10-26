"""HDOCK 对接技能的回归测试。

锁住老仓库 ``dock_hdock/scripts/`` 的 6 个 selenium 脚本带来的真实问题，
以及 HDOCK 提交协议里最容易踩的坑：

1. **提交到了首页而不是 runhdock.php** —— 首页 form 的 ``action`` 是空的，
   真实目标由 JS 在提交前改写。POST 到首页会返回**字节级相同的首页**
   （HTTP 200、无错误、无作业），脚本会误报"提交成功"。
2. **硬编码个人路径与 driver** —— ``C:\\Users\\DELL\\Desktop\\protein\\
   msedgedriver.exe``、``C:\\Users\\DELL\\Desktop\\hdock\\protein website.txt``。
3. **硬编码个人邮箱** —— ``lisenlikun@foxmail.com`` 被写死在提交代码里，
   别人用这个技能会往作者的邮箱发通知。
4. **多线程并发提交** —— 老脚本用 ``threading.Thread`` 起 N 个浏览器实例
   同时提交，公益服务器扛不住。
5. **"排队中"被当成"失败"** —— 老脚本取不到分数就写 N/A 并**中止整个
   批次**，一个未就绪的作业会让后面所有作业都不被检查。
6. **XPath 取分** —— ``//th[text()='Docking Score']/following-sibling::td``
   在页面改版时静默返回 None，不报错。
7. **只取排名第 1 的分数** —— 结果页是排名矩阵（行=指标、列=排名 1..10），
   只取第一格丢掉了全部对比信息。
8. **Confidence Score 要在页面上等** —— HDOCK 公开了闭式公式，可以直接算。

全部用例不联网。
"""

from __future__ import annotations

import ast
import math
import re
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
SKILL_DIR = REPO_ROOT / "skills" / "protein-docking-hdock"
FIXTURES = REPO_ROOT / "tests" / "fixtures_hdock"
sys.path.insert(0, str(REPO_ROOT / "packages"))
sys.path.insert(0, str(SKILL_DIR / "scripts"))

import fetch_docking_results  # noqa: E402
import hdock_client  # noqa: E402
import submit_docking  # noqa: E402

SCRIPTS = ("hdock_client.py", "submit_docking.py", "fetch_docking_results.py")


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
    result: list[str] = []
    for name in SCRIPTS:
        result.extend(_non_docstring_strings(_parse(name)))
    return result


def _fixture(name: str) -> str:
    return (FIXTURES / name).read_text(encoding="utf-8")


# --------------------------------------------------------------------------
# 提交协议
# --------------------------------------------------------------------------


class TestSubmitEndpoint:
    """提交地址——最容易踩的坑。"""

    def test_submit_url_is_runhdock(self):
        assert hdock_client.HDOCK_SUBMIT_URL.endswith("runhdock.php")

    def test_symmetric_submit_url(self):
        assert hdock_client.HDOCK_SYMMETRIC_SUBMIT_URL.endswith("runhsymdock.php")

    def test_not_posting_to_homepage(self):
        """不能把 form 直接提交到首页——那会返回字节级相同的首页。"""
        assert hdock_client.HDOCK_SUBMIT_URL != hdock_client.HDOCK_BASE_URL
        assert hdock_client.HDOCK_SYMMETRIC_SUBMIT_URL != hdock_client.HDOCK_BASE_URL

    def test_spellcheck(self):
        """``rundock.php`` 是 404——名字里有个 h。"""
        assert "rundock.php" not in hdock_client.HDOCK_SUBMIT_URL
        assert "rundock.php" not in hdock_client.HDOCK_SYMMETRIC_SUBMIT_URL

    def test_posts_to_submit_endpoint(self, monkeypatch):
        """submit_docking_job 必须 POST 到 runhdock.php 而不是首页。"""
        posted: dict = {}

        class _Response:
            url = "http://hdock.phys.hust.edu.cn/data/abc123456789/"
            text = ""

            def raise_for_status(self):
                return None

        class _Session:
            def get(self, url, timeout=None):
                return _Response()

            def post(self, url, data=None, files=None, timeout=None, allow_redirects=True):
                posted["url"] = url
                posted["data"] = data
                return _Response()

            def close(self):
                return None

        hdock_client.submit_docking_job(
            receptor_pdb_id="1CGI:E",
            ligand_pdb_id="1CGI:I",
            session=_Session(),
        )
        assert posted["url"] == hdock_client.HDOCK_SUBMIT_URL
        assert "runhdock.php" in posted["url"]

    def test_symmetric_uses_other_endpoint(self):
        posted: dict = {}

        class _Response:
            url = "http://hdock.phys.hust.edu.cn/data/abc123456789/"
            text = ""

            def raise_for_status(self):
                return None

        class _Session:
            def get(self, url, timeout=None):
                return _Response()

            def post(self, url, data=None, files=None, timeout=None, allow_redirects=True):
                posted["url"] = url
                return _Response()

            def close(self):
                return None

        hdock_client.submit_docking_job(
            receptor_pdb_id="1CGI:E",
            ligand_pdb_id="1CGI:I",
            symmetry="C2",
            session=_Session(),
        )
        assert posted["url"] == hdock_client.HDOCK_SYMMETRIC_SUBMIT_URL


class TestFormFields:
    """表单字段。"""

    def test_all_fields_present(self):
        fields = hdock_client.build_form_fields()
        for name in (
            "fastaseq1", "fastaseq2", "pdbid1", "pdbid2",
            "email", "jobname", "ligtyp", "symmetry",
            "sitenum1", "sitenum2", "restrnum", "upload",
        ):
            assert name in fields, name

    def test_receptor_and_ligand_fields(self):
        fields = hdock_client.build_form_fields(
            receptor_pdb_id="1CGI:E", ligand_pdb_id="1CGI:I"
        )
        assert fields["pdbid1"] == "1CGI:E"
        assert fields["pdbid2"] == "1CGI:I"

    def test_sequence_fields(self):
        fields = hdock_client.build_form_fields(
            receptor_seq=">r\nAAAA", ligand_seq=">l\nCCCC"
        )
        assert fields["fastaseq1"] == ">r\nAAAA"
        assert fields["fastaseq2"] == ">l\nCCCC"

    def test_file_fields_default_empty(self):
        fields = hdock_client.build_form_fields()
        assert fields["pdbfile1"] == ""
        assert fields["fastafile1"] == ""

    def test_submit_value(self):
        assert hdock_client.build_form_fields()["upload"] == "Submit"

    def test_ligand_types(self):
        assert "protein" in hdock_client.LIGAND_TYPES
        assert "rna" in hdock_client.LIGAND_TYPES
        assert "dna" in hdock_client.LIGAND_TYPES


class TestInputValidation:
    """受体与配体各需一种输入。"""

    def test_missing_receptor_raises(self):
        with pytest.raises(ValueError, match="受体"):
            hdock_client.submit_docking_job(ligand_pdb_id="1CGI:I")

    def test_missing_ligand_raises(self):
        with pytest.raises(ValueError, match="配体"):
            hdock_client.submit_docking_job(receptor_pdb_id="1CGI:E")

    def test_empty_raises(self):
        with pytest.raises(ValueError):
            hdock_client.submit_docking_job()

    def test_cli_returns_2_on_missing_input(self, tmp_path):
        code = submit_docking.main(
            ["--receptor-pdb-id", "1CGI:E", "--output-dir", str(tmp_path)]
        )
        assert code == 2


# --------------------------------------------------------------------------
# 老仓库缺陷的回归锁
# --------------------------------------------------------------------------


class TestNoHardcodedPersonalPaths:
    """老脚本硬编码了个人桌面路径与 driver。"""

    def test_no_windows_user_paths(self):
        pattern = re.compile(r"[A-Za-z]:\\\\?Users\\\\", re.IGNORECASE)
        for value in _all_strings():
            assert not pattern.search(value), f"硬编码 Windows 路径: {value!r}"

    def test_no_macos_user_paths(self):
        pattern = re.compile(r"/Users/[^/\s]+/")
        for value in _all_strings():
            assert not pattern.search(value), f"硬编码 macOS 路径: {value!r}"

    def test_no_driver_references(self):
        """docstring 里说明"老脚本硬编码了 driver"是必要文档，不算违规。"""
        for name in SCRIPTS:
            for value in _non_docstring_strings(_parse(name)):
                assert "msedgedriver" not in value
                assert "chromedriver" not in value
                assert "webdriver_manager" not in value


class TestNoHardcodedEmail:
    """老脚本把作者的个人邮箱写死在提交代码里。"""

    def test_no_personal_email(self):
        pattern = re.compile(r"[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}")
        for value in _all_strings():
            for match in pattern.findall(value):
                assert match in {"you@example.edu", "you@your-institution.edu"}, match

    def test_foxmail_not_present(self):
        for name in SCRIPTS:
            for value in _non_docstring_strings(_parse(name)):
                assert "foxmail" not in value.lower()

    def test_email_is_a_parameter(self):
        assert "email" in _script("submit_docking.py").read_text(encoding="utf-8")


class TestNoSelenium:
    """这个技能走纯 HTTP，不需要浏览器。"""

    def test_no_selenium_import(self):
        """docstring 里说明"老脚本用 selenium"是必要文档，不算违规。"""
        for name in SCRIPTS:
            for value in _non_docstring_strings(_parse(name)):
                assert "selenium" not in value.lower()
                assert "webdriver" not in value.lower()

    def test_requirements_do_not_need_selenium(self):
        text = (SKILL_DIR / "requirements.txt").read_text(encoding="utf-8")
        for line in text.splitlines():
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            assert "selenium" not in stripped.lower()


class TestNoConcurrentSubmission:
    """老脚本用线程池并发提交，公益服务器扛不住。"""

    def test_no_thread_pool(self):
        for name in SCRIPTS:
            for value in _non_docstring_strings(_parse(name)):
                assert "ThreadPoolExecutor" not in value
                assert "threading.Thread" not in value

    def test_no_threading_import(self):
        for name in SCRIPTS:
            tree = _parse(name)
            for node in ast.walk(tree):
                if isinstance(node, ast.Import):
                    for alias in node.names:
                        assert alias.name != "threading", f"{name} 导入了 threading"

    def test_no_global_output_variable(self):
        for name in SCRIPTS:
            tree = _parse(name)
            for node in ast.walk(tree):
                if isinstance(node, ast.Global):
                    assert node.names != ["output_tsv"]


class TestNoModuleLevelSideEffects:
    """模块顶层不得写文件或发请求。"""

    @staticmethod
    def _top_level_calls(name: str) -> list[str]:
        tree = _parse(name)
        calls: list[str] = []
        for node in tree.body:
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                continue
            if isinstance(node, ast.Expr) and isinstance(node.value, ast.Constant):
                continue
            for sub in ast.walk(node):
                if isinstance(sub, ast.Call):
                    func = sub.func
                    calls.append(
                        getattr(func, "id", None) or getattr(func, "attr", None) or "?"
                    )
        return calls

    def test_no_module_level_open(self):
        for name in SCRIPTS:
            assert "open" not in self._top_level_calls(name)

    def test_no_module_level_network(self):
        for name in SCRIPTS:
            for called in self._top_level_calls(name):
                assert called not in {"get", "post", "request", "urlopen"}


# --------------------------------------------------------------------------
# 结果页解析
# --------------------------------------------------------------------------


class TestParseFinishedPage:
    """已完成的结果页。"""

    @pytest.fixture(scope="class")
    @classmethod
    def result(cls):
        html = _fixture("result_example.html")
        return hdock_client.parse_result_page(
            html, job_name="example", url="http://hdock.phys.hust.edu.cn/data/example/"
        )

    def test_status_done(self, result):
        assert result.status == hdock_client.JobStatus.DONE

    def test_docking_score(self, result):
        assert result.docking_score == "-444.95"

    def test_docking_score_is_numeric(self, result):
        assert float(result.docking_score) < 0

    def test_confidence_score_computed(self, result):
        """页面没给 Confidence Score 时，用官方公式算。"""
        expected = 1.0 / (1.0 + math.exp(0.02 * (-444.95 + 150)))
        assert abs(float(result.confidence_score) - expected) < 1e-4

    def test_confidence_is_high_for_strong_score(self, result):
        assert float(result.confidence_score) > 0.9

    def test_rankings_preserve_full_row(self, result):
        """结果页是排名矩阵，整行都要保留。"""
        scores = result.rankings["Docking Score"]
        assert len(scores) == 10
        assert scores[0] == "-444.95"
        assert float(scores[1]) > float(scores[0])  # 排名 2 没那么负

    def test_rmsd_row_parsed(self, result):
        assert any("rmsd" in key.lower() for key in result.rankings)

    def test_rankings_are_descending(self, result):
        values = [float(v) for v in result.rankings["Docking Score"]]
        assert values == sorted(values)

    def test_job_id_extracted(self, result):
        assert result.job_id == "example"

    def test_fields_contains_score(self, result):
        assert result.fields["Docking Score"] == "-444.95"


class TestParseQueuedPage:
    """排队中的结果页——必须与"失败"区分开。"""

    @pytest.fixture(scope="class")
    @classmethod
    def result(cls):
        html = _fixture("result_queued.html")
        return hdock_client.parse_result_page(
            html, job_name="queued", url="http://hdock.phys.hust.edu.cn/data/6adf1c9e1eb6e/"
        )

    def test_status_is_queued_not_failed(self, result):
        assert result.status == hdock_client.JobStatus.QUEUED

    def test_queued_is_not_done(self, result):
        assert result.status != hdock_client.JobStatus.DONE

    def test_queued_has_no_score(self, result):
        assert not result.docking_score

    def test_queued_has_no_errors(self, result):
        assert not result.errors

    def test_job_id_from_url(self, result):
        assert result.job_id == "6adf1c9e1eb6e"


class TestParseErrorPage:
    """失败的结果页。"""

    def test_error_detected(self):
        html = (
            "<html><body><h2>Your HDOCK results for job <i>x</i></h2>"
            "<h2>ERROR: No structure found for the <font color=red>Receptor</font> molecule.</h2>"
            "</body></html>"
        )
        result = hdock_client.parse_result_page(html, url="http://hdock.phys.hust.edu.cn/data/x/")
        assert result.status == hdock_client.JobStatus.FAILED
        assert result.errors
        assert "Receptor" in result.errors[0]

    def test_error_has_no_score(self):
        html = "<h2>ERROR: something went wrong</h2>"
        result = hdock_client.parse_result_page(html)
        assert not result.docking_score

    def test_empty_html_is_unknown(self):
        result = hdock_client.parse_result_page("")
        assert result.status == hdock_client.JobStatus.UNKNOWN


class TestConfidenceScoreFormula:
    """官方闭式公式。"""

    def test_known_value(self):
        assert abs(
            hdock_client.compute_confidence_score(-444.95) - 0.99727
        ) < 1e-4

    def test_matches_documented_formula(self):
        for score in (-500.0, -300.0, -200.0, -100.0, 0.0):
            expected = 1.0 / (1.0 + math.exp(0.02 * (score + 150)))
            assert abs(
                hdock_client.compute_confidence_score(score) - expected
            ) < 1e-12

    def test_monotonic(self):
        """分数越负（越好），置信度越高。"""
        assert hdock_client.compute_confidence_score(-400) > hdock_client.compute_confidence_score(-200)

    def test_threshold_semantics(self):
        """文档说 >0.7 表示很可能结合。"""
        assert hdock_client.compute_confidence_score(-200) > 0.7
        assert hdock_client.compute_confidence_score(-100) < 0.7

    def test_invalid_input(self):
        assert hdock_client.compute_confidence_score("abc") == 0.0
        assert hdock_client.compute_confidence_score("") == 0.0

    def test_range(self):
        for score in (-1000.0, -200.0, 100.0):
            value = hdock_client.compute_confidence_score(score)
            assert 0.0 <= value <= 1.0


class TestResultRow:
    """CSV 输出行。"""

    def test_row_has_all_fields(self):
        result = hdock_client.DockingResult(
            job_name="x", job_id="y", url="z", status=hdock_client.JobStatus.DONE,
            docking_score="-444.95", confidence_score="0.9973",
        )
        row = result.to_row()
        for key in (
            "job_name", "job_id", "url", "status",
            "docking_score", "confidence_score", "error",
        ):
            assert key in row

    def test_computed_confidence_property(self):
        result = hdock_client.DockingResult(
            job_name="x", status=hdock_client.JobStatus.DONE, docking_score="-444.95"
        )
        assert abs(result.confidence_score_computed - 0.99727) < 1e-4


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------


class TestSequenceArgument:
    """``--*-seq`` 既接受文本也接受文件路径。"""

    def test_reads_file(self, tmp_path):
        path = tmp_path / "rec.fasta"
        path.write_text(">receptor\nMSTAVL\n", encoding="utf-8")
        assert submit_docking._read_sequence(str(path)).strip() == ">receptor\nMSTAVL"

    def test_passes_through_inline_fasta(self):
        assert submit_docking._read_sequence(">r\nAAAA").startswith(">r")

    def test_empty(self):
        assert submit_docking._read_sequence("") == ""


class TestFetchCli:
    """取结果命令。"""

    def test_requires_url_or_jobs_file(self, tmp_path):
        code = fetch_docking_results.main(["--output-dir", str(tmp_path)])
        assert code == 2

    def test_missing_jobs_file(self, tmp_path):
        code = fetch_docking_results.main(
            ["--jobs-file", str(tmp_path / "nope.csv"), "--output-dir", str(tmp_path)]
        )
        assert code == 2

    def test_reports_queued_state(self, monkeypatch, tmp_path):
        """排队中的作业必须报 queued，而不是当成失败。"""
        html = _fixture("result_queued.html")
        monkeypatch.setattr(
            fetch_docking_results, "fetch_result_page", lambda url, timeout=60: html
        )
        out = tmp_path / "out"
        code = fetch_docking_results.main(
            ["--url", "http://hdock.phys.hust.edu.cn/data/6adf1c9e1eb6e/",
             "--output-dir", str(out)]
        )
        assert code == 0
        content = (out / "results.csv").read_text(encoding="utf-8-sig")
        assert "queued" in content

    def test_reports_done_state(self, monkeypatch, tmp_path):
        html = _fixture("result_example.html")
        monkeypatch.setattr(
            fetch_docking_results, "fetch_result_page", lambda url, timeout=60: html
        )
        out = tmp_path / "out"
        code = fetch_docking_results.main(
            ["--url", "http://hdock.phys.hust.edu.cn/data/example/",
             "--output-dir", str(out)]
        )
        assert code == 0
        content = (out / "results.csv").read_text(encoding="utf-8-sig")
        assert "done" in content
        assert "-444.95" in content


# --------------------------------------------------------------------------
# 文档与元数据
# --------------------------------------------------------------------------


class TestSkillMetadata:
    def test_skill_md_exists(self):
        assert (SKILL_DIR / "SKILL.md").is_file()

    def test_references_exist(self):
        assert (SKILL_DIR / "references" / "methodology.md").is_file()
        assert (SKILL_DIR / "references" / "troubleshooting.md").is_file()

    def test_requirements_declares_requests(self):
        text = (SKILL_DIR / "requirements.txt").read_text(encoding="utf-8")
        assert "requests" in text

    def test_no_local_path_dependency(self):
        text = (SKILL_DIR / "requirements.txt").read_text(encoding="utf-8")
        for line in text.splitlines():
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            assert "file://" not in stripped

    def test_documentation_mentions_the_endpoint_gotcha(self):
        methodology = (SKILL_DIR / "references" / "methodology.md").read_text(
            encoding="utf-8"
        )
        assert "runhdock.php" in methodology
        assert "runhsymdock.php" in methodology

    def test_documentation_mentions_the_confidence_formula(self):
        methodology = (SKILL_DIR / "references" / "methodology.md").read_text(
            encoding="utf-8"
        )
        assert "Confidence_score" in methodology

    def test_skill_md_mentions_queued_state(self):
        text = (SKILL_DIR / "SKILL.md").read_text(encoding="utf-8")
        assert "queued" in text
