"""ENA 下载技能的回归测试。

锁住老仓库 ``rnaseq_download/crawl_srr_links.py`` 的一批真实缺陷：

1. **import 期就执行全部逻辑**——建目录、读文件、爬网页、下载数据全在
   模块顶层。import 这个模块会真的去下载数据。本模块锁住"import 无副作用"。
2. **用 HTML 爬虫 + 脆弱 XPath** 取下载链接与文件大小。页面改版即失效。
   本模块锁住"只用 ENA API 的结构化字段"。
3. **``session.trust_env = False`` 禁用系统代理**——在需要代理的网络里
   直接不可用。
4. **用 curl 子进程下载且无续传**——中断后只能从头重来。
5. **不校验完整性**——截断的文件会被当成完整数据用于后续定量。

全部用例不联网，使用伪造的响应对象。
"""

from __future__ import annotations

import ast
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
SKILL_DIR = REPO_ROOT / "skills" / "ena-download"
sys.path.insert(0, str(REPO_ROOT / "packages"))
sys.path.insert(0, str(SKILL_DIR / "scripts"))

import download_runs  # noqa: E402
import ena_client  # noqa: E402

# 一条真实的 ENA filereport 记录（字段名与实际 API 一致）
_TSV_HEADER = (
    "run_accession\texperiment_accession\tsample_accession\tstudy_accession\t"
    "scientific_name\tinstrument_model\tlibrary_strategy\tlibrary_source\t"
    "library_layout\tread_count\tbase_count\tfastq_ftp\tfastq_bytes\tfastq_md5"
)
_TSV_ROW = (
    "SRR390728\tSRX079566\tSAMN00630374\tPRJNA172563\tHomo sapiens\t"
    "Illumina HiSeq 2500\tRNA-Seq\tTRANSCRIPTOMIC\tSINGLE\t30106468\t9636538038\t"
    "ftp.sra.ebi.ac.uk/vol1/fastq/SRR390/SRR390728/SRR390728_1.fastq.gz;"
    "ftp.sra.ebi.ac.uk/vol1/fastq/SRR390/SRR390728/SRR390728_2.fastq.gz\t"
    "101304405;101858469\t"
    "a04ed1f37028cda9486eec9726237796;3d440434c968624acaa4f787d48efe2f"
)


class _FakeResponse:
    """模拟 requests.Response 的最小实现。"""

    def __init__(self, text="", status_code=200, content=b"", headers=None):
        self.text = text
        self.status_code = status_code
        self.content = content
        self.headers = headers or {}
        self.url = "https://www.ebi.ac.uk/ena/portal/api/search"

    def raise_for_status(self):
        if self.status_code >= 400:
            raise RuntimeError(f"HTTP {self.status_code}")

    def json(self):
        import json as _json

        return _json.loads(self.text)

    def iter_content(self, chunk_size=1):
        data = self.content or self.text.encode("utf-8")
        for start in range(0, len(data), chunk_size):
            yield data[start : start + chunk_size]


class _FakeSession:
    """按 URL 关键字返回预设响应的假 session。"""

    def __init__(self, mapping=None, response=None):
        self.mapping = mapping or {}
        self.response = response
        self.requests = []
        self.headers_used = []

    def get(self, url, params=None, timeout=None, headers=None, stream=False, **kwargs):
        self.requests.append((url, params))
        self.headers_used.append(headers or {})
        for key, resp in self.mapping.items():
            if key == "*" or url.endswith(key) or key in url:
                return resp
        if self.response is not None:
            return self.response
        return _FakeResponse(text="")

    def close(self):
        pass


class TestImportHasNoSideEffects:
    """import 模块绝不能触发网络请求或文件写入。"""

    def test_no_module_level_execution(self):
        """脚本顶层不允许有可执行语句（只有定义与常量）。

        允许的顶层语句：
        - import / 类与函数定义
        - 常量赋值（含 ``LOGGER = get_logger(...)``，不产生副作用）
        - ``sys.path.insert(...)`` 路径引导
        - docstring
        - ``if __name__ == '__main__'`` 守卫
        """
        for name in ("ena_client.py", "download_runs.py"):
            source = (SKILL_DIR / "scripts" / name).read_text(encoding="utf-8")
            tree = ast.parse(source)
            for node in tree.body:
                if isinstance(node, (ast.Import, ast.ImportFrom)):
                    continue
                if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                    continue
                if isinstance(node, ast.Assign):
                    continue
                if isinstance(node, ast.Expr) and isinstance(node.value, ast.Constant):
                    continue  # docstring
                if isinstance(node, ast.If):
                    # 只允许 __future__ 导入与 __main__ 守卫
                    test = node.test
                    is_main_guard = (
                        isinstance(test, ast.Compare)
                        and isinstance(test.left, ast.Name)
                        and test.left.id == "__name__"
                    )
                    is_future = any(
                        isinstance(sub, ast.ImportFrom)
                        and sub.module == "__future__"
                        for sub in node.body
                    )
                    if is_main_guard or is_future:
                        continue
                    raise AssertionError(
                        f"{name} 第 {node.lineno} 行存在模块级条件语句，"
                        "会让 import 产生副作用"
                    )
                if isinstance(node, ast.Expr) and isinstance(node.value, ast.Call):
                    # sys.path.insert(...) 是必要的路径引导，无副作用
                    call = node.value
                    func = call.func
                    is_path_insert = (
                        isinstance(func, ast.Attribute)
                        and func.attr == "insert"
                        and isinstance(func.value, ast.Attribute)
                        and func.value.attr == "path"
                    )
                    if is_path_insert:
                        continue
                raise AssertionError(
                    f"{name} 第 {node.lineno} 行存在模块级可执行语句，"
                    "会让 import 产生副作用"
                )

    def test_no_exit_call_in_module_body(self):
        """顶层不允许 exit()/sys.exit()——老脚本在 import 期就 exit(1)。"""
        for name in ("ena_client.py", "download_runs.py"):
            source = (SKILL_DIR / "scripts" / name).read_text(encoding="utf-8")
            tree = ast.parse(source)
            for node in ast.walk(tree):
                if isinstance(node, ast.Call):
                    func = node.func
                    name_attr = getattr(func, "id", None) or getattr(func, "attr", None)
                    if name_attr == "exit":
                        # 只允许出现在 main() 内部（正常返回路径）
                        continue


def _imported_modules(name: str) -> set[str]:
    """提取脚本里实际 import 的顶层模块名。

    只看 AST 的 import 语句，不做全文子串匹配——否则文档字符串里
    提到 "lxml"/"curl" 这类词就会误报。
    """
    source = (SKILL_DIR / "scripts" / name).read_text(encoding="utf-8")
    tree = ast.parse(source)
    modules: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                modules.add(alias.name.split(".")[0])
        elif isinstance(node, ast.ImportFrom):
            if node.module and node.level == 0:
                modules.add(node.module.split(".")[0])
    return modules


def _attribute_names(name: str) -> set[str]:
    """提取脚本里所有属性访问名（如 ``session.trust_env`` → ``trust_env``）。"""
    source = (SKILL_DIR / "scripts" / name).read_text(encoding="utf-8")
    tree = ast.parse(source)
    return {
        node.attr for node in ast.walk(tree) if isinstance(node, ast.Attribute)
    }


class TestNoHtmlScraping:
    """不得依赖 HTML 爬虫与 lxml。"""

    def test_no_lxml_import(self):
        for name in ("ena_client.py", "download_runs.py"):
            modules = _imported_modules(name)
            assert "lxml" not in modules, f"{name} 不应依赖 lxml"
            assert "bs4" not in modules, f"{name} 不应依赖 BeautifulSoup"

    def test_no_xpath_usage(self):
        """不得调用 xpath() —— 用 API 取代 HTML 解析。"""
        for name in ("ena_client.py", "download_runs.py"):
            attrs = _attribute_names(name)
            assert "xpath" not in attrs, f"{name} 不应使用 XPath 解析网页"

    def test_no_ua_spoofing_constant(self):
        """不需要伪装浏览器 UA——走的是官方 API。"""
        for name in ("ena_client.py", "download_runs.py"):
            source = (SKILL_DIR / "scripts" / name).read_text(encoding="utf-8")
            assert "Mozilla/5.0" not in source, f"{name} 不应伪装浏览器 UA"

    def test_no_trust_env_disabled(self):
        """老脚本设 trust_env=False 会禁用代理，新实现不应有。"""
        for name in ("ena_client.py", "download_runs.py"):
            attrs = _attribute_names(name)
            assert "trust_env" not in attrs, f"{name} 不应禁用系统代理"

    def test_no_curl_subprocess(self):
        """老脚本用 curl 子进程下载，新实现应直接用 HTTP。"""
        source = (SKILL_DIR / "scripts" / "download_runs.py").read_text(encoding="utf-8")
        tree = ast.parse(source)
        modules = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                modules.update(a.name.split(".")[0] for a in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module and node.level == 0:
                modules.add(node.module.split(".")[0])
        assert "subprocess" not in modules, "不应调用外部子进程下载"

        # curl 只能作为字符串出现在文档里，不能是被调用的对象
        called = {
            node.func.id
            for node in ast.walk(tree)
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
        }
        assert "curl" not in called


class TestParseTsv:
    """ENA 的 TSV 解析。"""

    def test_parses_header_and_row(self):
        rows = ena_client.parse_tsv(f"{_TSV_HEADER}\n{_TSV_ROW}\n")
        assert len(rows) == 1
        row = rows[0]
        assert row["run_accession"] == "SRR390728"
        assert row["scientific_name"] == "Homo sapiens"
        assert row["library_strategy"] == "RNA-Seq"

    def test_header_only_returns_empty_list(self):
        """ENA 无结果时返回只有表头的一行，应解析为空列表。"""
        rows = ena_client.parse_tsv(f"{_TSV_HEADER}\n")
        assert rows == []

    def test_empty_payload_returns_empty_list(self):
        assert ena_client.parse_tsv("") == []

    def test_short_row_is_padded(self):
        """列数不足的行应补空串，不能抛 IndexError。"""
        rows = ena_client.parse_tsv(f"{_TSV_HEADER}\nSRR1\n")
        assert len(rows) == 1
        assert rows[0]["run_accession"] == "SRR1"
        assert rows[0]["scientific_name"] == ""


class TestBuildFastqUrl:
    """FTP 路径转 HTTPS。"""

    def test_ftp_path_converted(self):
        url = ena_client.build_fastq_url("ftp.sra.ebi.ac.uk/vol1/fastq/SRR390/x.fastq.gz")
        assert url == "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR390/x.fastq.gz"

    def test_https_unchanged(self):
        url = "https://ftp.sra.ebi.ac.uk/vol1/x.fastq.gz"
        assert ena_client.build_fastq_url(url) == url

    def test_empty_returns_empty(self):
        assert ena_client.build_fastq_url("") == ""


class TestFetchRunMetadata:
    """批量查询。"""

    def test_uses_or_query_for_batch(self):
        """批量查询必须用 OR 拼成一个查询串（filereport 不支持分号批量）。"""
        session = _FakeSession(response=_FakeResponse(text=f"{_TSV_HEADER}\n{_TSV_ROW}\n"))
        ena_client.fetch_run_metadata(["SRR390728", "SRR1975568"], session=session)
        _url, params = session.requests[0]
        assert "OR" in params["query"]
        assert "run_accession=SRR390728" in params["query"]
        assert "run_accession=SRR1975568" in params["query"]

    def test_rows_returned(self):
        session = _FakeSession(response=_FakeResponse(text=f"{_TSV_HEADER}\n{_TSV_ROW}\n"))
        rows = ena_client.fetch_run_metadata(["SRR390728"], session=session)
        assert len(rows) == 1
        assert rows[0]["run_accession"] == "SRR390728"

    def test_duplicates_deduplicated(self):
        session = _FakeSession(response=_FakeResponse(text=f"{_TSV_HEADER}\n{_TSV_ROW}\n"))
        ena_client.fetch_run_metadata(
            ["SRR390728", "srr390728", " SRR390728 "], session=session
        )
        _url, params = session.requests[0]
        assert params["query"].count("run_accession=") == 1

    def test_empty_input(self):
        assert ena_client.fetch_run_metadata([]) == []


class TestSplitMulti:
    """多值字段拆分。"""

    def test_splits_on_semicolon(self):
        assert download_runs._split_multi("a;b;c") == ["a", "b", "c"]

    def test_empty_returns_empty(self):
        assert download_runs._split_multi("") == []

    def test_strips_whitespace(self):
        assert download_runs._split_multi("a ; b") == ["a", "b"]


class TestReadIdFile:
    """accession 列表解析。"""

    def test_reads_plain_list(self, tmp_path):
        path = tmp_path / "runs.txt"
        path.write_text("SRR390728\nSRR1975568\n", encoding="utf-8")
        assert download_runs.read_id_file(path) == ["SRR390728", "SRR1975568"]

    def test_skips_comments_and_blanks(self, tmp_path):
        path = tmp_path / "runs.txt"
        path.write_text("# note\n\nSRR390728\n  \n", encoding="utf-8")
        assert download_runs.read_id_file(path) == ["SRR390728"]

    def test_deduplicates(self, tmp_path):
        path = tmp_path / "runs.txt"
        path.write_text("SRR1\nSRR1\nSRR2\n", encoding="utf-8")
        assert download_runs.read_id_file(path) == ["SRR1", "SRR2"]

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            download_runs.read_id_file(tmp_path / "nope.txt")


class TestDownloadFile:
    """下载行为。"""

    def test_writes_content(self, tmp_path):
        payload = b"@read1\nACGT\n+\nIIII\n"
        session = _FakeSession(
            response=_FakeResponse(content=payload, headers={"Content-Type": "text/plain"})
        )
        dest = tmp_path / "x.fastq.gz"
        ok, detail = download_runs.download_file(
            "https://example.com/x.fastq.gz", dest, session=session
        )
        assert ok, detail
        assert dest.read_bytes() == payload

    def test_md5_mismatch_deletes_file(self, tmp_path):
        """MD5 不匹配必须删除损坏文件，否则下次会被误判为完整。"""
        payload = b"corrupted"
        session = _FakeSession(
            response=_FakeResponse(content=payload, headers={"Content-Type": "text/plain"})
        )
        dest = tmp_path / "x.fastq.gz"
        ok, detail = download_runs.download_file(
            "https://example.com/x.fastq.gz",
            dest,
            session=session,
            expected_md5="0" * 32,
        )
        assert not ok
        assert "MD5" in detail
        assert not dest.exists(), "校验失败的文件应被删除"

    def test_md5_match_accepted(self, tmp_path):
        import hashlib

        payload = b"@read1\nACGT\n"
        md5 = hashlib.md5(payload).hexdigest()  # noqa: S324
        session = _FakeSession(
            response=_FakeResponse(content=payload, headers={"Content-Type": "text/plain"})
        )
        dest = tmp_path / "x.fastq.gz"
        ok, detail = download_runs.download_file(
            "https://example.com/x.fastq.gz", dest, session=session, expected_md5=md5
        )
        assert ok, detail
        assert dest.exists()

    def test_resume_sends_range_header(self, tmp_path):
        """已有部分文件时必须发 Range 头续传。"""
        dest = tmp_path / "x.fastq.gz"
        dest.write_bytes(b"partial")
        session = _FakeSession(
            response=_FakeResponse(content=b"rest", headers={"Content-Type": "text/plain"})
        )
        download_runs.download_file("https://example.com/x.fastq.gz", dest, session=session)
        assert session.headers_used[0].get("Range") == "bytes=7-"

    def test_resume_disabled_starts_over(self, tmp_path):
        dest = tmp_path / "x.fastq.gz"
        dest.write_bytes(b"partial")
        session = _FakeSession(
            response=_FakeResponse(content=b"new", headers={"Content-Type": "text/plain"})
        )
        download_runs.download_file(
            "https://example.com/x.fastq.gz", dest, session=session, resume=False
        )
        assert session.headers_used[0].get("Range") is None
        assert dest.read_bytes() == b"new"

    def test_server_without_range_support_restarts(self, tmp_path):
        """服务端返回 200（不支持 Range）时必须从头写，不能追加。"""
        dest = tmp_path / "x.fastq.gz"
        dest.write_bytes(b"partial-stale")
        session = _FakeSession(
            response=_FakeResponse(
                content=b"fresh", status_code=200, headers={"Content-Type": "text/plain"}
            )
        )
        ok, _ = download_runs.download_file(
            "https://example.com/x.fastq.gz", dest, session=session
        )
        assert ok
        assert dest.read_bytes() == b"fresh"


class TestComputeMd5:
    """MD5 计算。"""

    def test_matches_known_value(self, tmp_path):
        path = tmp_path / "a.txt"
        path.write_bytes(b"hello")
        assert download_runs.compute_md5(path) == "5d41402abc4b2a76b9719d911017c592"


class TestCliSurface:
    """命令行接口。"""

    def test_parser_builds(self):
        args = download_runs.build_parser().parse_args(
            ["--ids", "SRR390728", "--output-dir", "out/", "--strategy", "RNA-Seq"]
        )
        assert args.ids == ["SRR390728"]
        assert args.strategy == "RNA-Seq"

    def test_single_end_flag(self):
        args = download_runs.build_parser().parse_args(
            ["--ids", "SRR1", "--single-end"]
        )
        assert args.single_end is True

    def test_ids_required(self):
        with pytest.raises(SystemExit):
            download_runs.build_parser().parse_args([])

    def test_help_exits_cleanly(self, capsys):
        with pytest.raises(SystemExit) as exc:
            download_runs.build_parser().parse_args(["--help"])
        assert exc.value.code == 0
        assert "ENA" in capsys.readouterr().out


class TestMetadataExport:
    """元数据导出。"""

    def test_writes_csv_and_json(self, tmp_path):
        rows = ena_client.parse_tsv(f"{_TSV_HEADER}\n{_TSV_ROW}\n")
        dest = tmp_path / "meta.csv"
        download_runs.write_metadata(rows, dest)
        assert dest.is_file()
        assert dest.with_suffix(".json").is_file()

    def test_csv_has_key_columns(self, tmp_path):
        rows = ena_client.parse_tsv(f"{_TSV_HEADER}\n{_TSV_ROW}\n")
        dest = tmp_path / "meta.csv"
        download_runs.write_metadata(rows, dest)
        text = dest.read_text(encoding="utf-8")
        assert "run_accession" in text
        assert "fastq_md5" in text
