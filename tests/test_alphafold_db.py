"""AlphaFold DB 技能的回归测试。

这些用例锁住一批**曾经真实存在过**的缺陷：

1. 老仓库 ``afdb_download`` 把模型版本号 ``model_v6`` 写死在 URL 里。
   AlphaFold DB 重训模型后版本号会递增，写死必然在下一版 404。
   本模块用 AST 静态检查确保代码里不存在写死的版本号。
2. 老仓库把 API 响应按 dict 处理，但 ``/api/prediction/{acc}`` 对单个
   accession 返回的是**数组**。本模块锁住三种响应形态的处理。
3. 老仓库用 ``response.text`` + 文本模式写文件，Windows 下会做换行转换，
   破坏 PDB 的定宽列记录。本模块锁住二进制写入。
4. 老仓库无条件跳过 ID 文件第一行，纯 ID 文件会静默丢失首个 ID。
5. 老仓库把 HTML 错误页当 PDB 存下。本模块锁住 Content-Type 校验。

全部用例不联网，使用伪造的响应对象。
"""

from __future__ import annotations

import ast
import re
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
SKILL_DIR = REPO_ROOT / "skills" / "alphafold-db"
sys.path.insert(0, str(REPO_ROOT / "packages"))
sys.path.insert(0, str(SKILL_DIR / "scripts"))

import alphafold_client  # noqa: E402
import download_structures  # noqa: E402


class _FakeResponse:
    """模拟 requests.Response 的最小实现。"""

    def __init__(self, payload=None, status_code=200, content=b"", headers=None):
        self._payload = payload
        self.status_code = status_code
        self.content = content
        self.headers = headers or {}
        self.text = content.decode("utf-8", errors="replace") if content else ""
        self.url = "https://alphafold.ebi.ac.uk/files/AF-P69905-F1-model_v6.pdb"

    def raise_for_status(self):
        if self.status_code >= 400:
            raise RuntimeError(f"HTTP {self.status_code}")

    def json(self):
        return self._payload


class _FakeSession:
    """按 URL 关键字返回预设响应的假 session。

    键如果是 ``*`` 则匹配任意 URL；否则按 "URL 以该键结尾" 匹配，
    用于模拟 ``/api/prediction/{accession}`` 这类路径。
    """

    def __init__(self, mapping):
        self.mapping = mapping
        self.requests = []

    def get(self, url, timeout=None, **kwargs):
        self.requests.append(url)
        for key, response in self.mapping.items():
            if key == "*" or url.endswith(f"/{key}") or url.endswith(key):
                return response
        return _FakeResponse(status_code=404)

    def close(self):
        pass


# 一条真实的 AlphaFold DB 记录（字段名与实际 API 一致）
_SAMPLE_ENTRY = {
    "entryId": "AF-P69905-F1",
    "uniprotAccession": "P69905",
    "uniprotId": "HBA_HUMAN",
    "uniprotDescription": "Hemoglobin subunit alpha",
    "gene": "HBA1",
    "organismScientificName": "Homo sapiens",
    "taxId": "9606",
    "chainId": "A",
    "modelCreatedDate": "2025-08-01T00:00:00Z",
    "latestVersion": 6,
    "globalMetricValue": 98.06,
    "fractionPlddtVeryLow": 0.0,
    "fractionPlddtLow": 0.0,
    "fractionPlddtConfident": 0.0146,
    "fractionPlddtVeryHigh": 0.9854,
    "sequence": "M" * 142,
    "pdbUrl": "https://alphafold.ebi.ac.uk/files/AF-P69905-F1-model_v6.pdb",
    "cifUrl": "https://alphafold.ebi.ac.uk/files/AF-P69905-F1-model_v6.cif",
    "isReviewed": True,
    "isComplex": False,
}


class TestResponseShapeHandling:
    """API 的三种响应形态都必须被正确处理。"""

    def test_list_response_is_normalised(self):
        """单个 accession 的正常响应是数组，不能当成 dict 处理。"""
        entries = list(alphafold_client._iter_entries([_SAMPLE_ENTRY]))
        assert len(entries) == 1
        assert entries[0]["uniprotAccession"] == "P69905"

    def test_dict_response_is_accepted(self):
        """少数情况下直接返回对象，也应被接受。"""
        entries = list(alphafold_client._iter_entries(_SAMPLE_ENTRY))
        assert len(entries) == 1

    def test_error_object_is_not_treated_as_record(self):
        """``{"error": ...}`` 是错误信息，不能当成元数据记录。"""
        payload = {
            "error": "Invalid identifier format. Please use a UniProt accession "
            "or a supported AlphaFold DB ID."
        }
        assert list(alphafold_client._iter_entries(payload)) == []

    def test_error_object_with_accession_is_record(self):
        """带 accession 的 dict 仍应被当作记录（防御性）。"""
        entries = list(alphafold_client._iter_entries(_SAMPLE_ENTRY))
        assert entries and entries[0]["uniprotAccession"] == "P69905"


class TestParsePrediction:
    """元数据扁平化。"""

    def test_core_fields_extracted(self):
        record = alphafold_client.parse_prediction(_SAMPLE_ENTRY)
        assert record["accession"] == "P69905"
        assert record["gene"] == "HBA1"
        assert record["global_plddt"] == 98.06
        assert record["model_version"] == 6
        assert record["pdb_url"].endswith(".pdb")

    def test_sequence_length_computed(self):
        record = alphafold_client.parse_prediction(_SAMPLE_ENTRY)
        assert record["sequence_length"] == 142

    def test_missing_fields_become_empty_not_keyerror(self):
        """字段缺失时返回空串而不是抛 KeyError。"""
        record = alphafold_client.parse_prediction({})
        assert record["accession"] == ""
        assert record["global_plddt"] == ""

    def test_zero_values_are_preserved_not_treated_as_missing(self):
        """0.0 是有意义的值，不能被当成缺失转成空串。

        AlphaFold 的置信度占比字段里 0.0 完全合法（某个区间的残基
        占比恰好为 0）。用 ``or ""`` 会把它变成空字符串，导致下游
        求和校验失败。
        """
        entry = dict(_SAMPLE_ENTRY)
        entry["fractionPlddtVeryLow"] = 0.0
        entry["fractionPlddtLow"] = 0.0
        record = alphafold_client.parse_prediction(entry)
        assert record["fraction_plddt_very_low"] == 0.0
        assert record["fraction_plddt_low"] == 0.0
        assert record["fraction_plddt_confident"] == 0.0146

    def test_false_values_are_preserved(self):
        """False 也是有效值（如 isComplex=False）。"""
        entry = dict(_SAMPLE_ENTRY)
        entry["isComplex"] = False
        record = alphafold_client.parse_prediction(entry)
        assert record["is_complex"] is False


class TestFetchPredictions:
    """批量查询逻辑。"""

    def test_returns_records_for_valid_accessions(self):
        session = _FakeSession({"P69905": _FakeResponse(payload=[_SAMPLE_ENTRY])})
        rows = alphafold_client.fetch_predictions(["P69905"], session=session)
        assert len(rows) == 1
        assert rows[0]["accession"] == "P69905"

    def test_missing_accession_is_skipped_not_fatal(self):
        session = _FakeSession({})
        rows = alphafold_client.fetch_predictions(["NOPE00"], session=session)
        assert rows == []

    def test_duplicate_accessions_are_deduplicated(self):
        session = _FakeSession({"P69905": _FakeResponse(payload=[_SAMPLE_ENTRY])})
        rows = alphafold_client.fetch_predictions(
            ["P69905", "p69905", " P69905 "], session=session
        )
        assert len(rows) == 1
        assert len(session.requests) == 1

    def test_empty_input_returns_empty(self):
        assert alphafold_client.fetch_predictions([]) == []

    def test_accession_falls_back_when_record_lacks_it(self):
        """记录里没有 accession 时用请求的 ID 补上。"""
        stripped = {k: v for k, v in _SAMPLE_ENTRY.items() if k != "uniprotAccession"}
        session = _FakeSession({"P69905": _FakeResponse(payload=[stripped])})
        rows = alphafold_client.fetch_predictions(["P69905"], session=session)
        assert rows[0]["accession"] == "P69905"


class TestNoHardcodedModelVersion:
    """模型版本号绝不能写死——这是老脚本最致命的缺陷。"""

    def test_no_hardcoded_model_v6_in_source(self):
        """源码里不允许出现写死的 ``model_vN`` 版本号字面量。

        只检查**非 docstring** 的字符串常量：文档字符串里讨论
        ``model_v6`` 是必要的说明文字，但代码里拼 URL 时写死版本号
        才是缺陷。字段名 ``model_version`` 之类的普通标识不受影响。
        """
        source = (SKILL_DIR / "scripts" / "alphafold_client.py").read_text(
            encoding="utf-8"
        )
        tree = ast.parse(source)

        # 收集所有 docstring 节点，豁免它们
        docstring_nodes = set()
        for node in ast.walk(tree):
            if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
                body = getattr(node, "body", [])
                if body and isinstance(body[0], ast.Expr):
                    value = body[0].value
                    if isinstance(value, ast.Constant) and isinstance(value.value, str):
                        docstring_nodes.add(id(value))

        # 匹配 model_v + 数字 的版本号形态（model_v6 / model_v10）
        version_pattern = re.compile(r"model_v\d")

        offenders = []
        for node in ast.walk(tree):
            if not isinstance(node, ast.Constant) or not isinstance(node.value, str):
                continue
            if id(node) in docstring_nodes:
                continue
            if version_pattern.search(node.value):
                offenders.append((node.lineno, node.value))

        assert not offenders, f"发现写死的模型版本号: {offenders}"

    def test_fallback_versions_are_not_absolute(self):
        """兜底版本列表只能是候选，不能是唯一固定值。"""
        versions = alphafold_client._FALLBACK_MODEL_VERSIONS
        assert len(versions) > 1, "兜底版本应提供多个候选，避免单一版本失效"

    def test_build_pdb_url_does_not_hardcode_single_version(self):
        """不指定版本时也不能返回写死的单一版本 URL。"""
        url = alphafold_client.build_pdb_url("P69905")
        # URL 里必须含 model_ 前缀，且版本来自候选列表
        assert "model_" in url
        assert any(f"model_{v}" in url for v in alphafold_client._FALLBACK_MODEL_VERSIONS)


class TestReadIdFile:
    """ID 文件解析——老脚本在这里会静默丢 ID。"""

    def test_header_line_is_skipped(self, tmp_path):
        path = tmp_path / "ids.txt"
        path.write_text("id\nP69905\nP04637\n", encoding="utf-8")
        assert download_structures.read_id_file(path) == ["P69905", "P04637"]

    def test_pure_id_file_keeps_first_entry(self, tmp_path):
        """纯 ID 文件（无表头）绝不能丢掉第一个 ID。"""
        path = tmp_path / "ids.txt"
        path.write_text("P69905\nP04637\n", encoding="utf-8")
        assert download_structures.read_id_file(path) == ["P69905", "P04637"]

    def test_comments_and_blank_lines_ignored(self, tmp_path):
        path = tmp_path / "ids.txt"
        path.write_text("# a comment\n\nP69905\n  \nP04637\n", encoding="utf-8")
        assert download_structures.read_id_file(path) == ["P69905", "P04637"]

    def test_duplicates_removed_preserving_order(self, tmp_path):
        path = tmp_path / "ids.txt"
        path.write_text("P69905\nP04637\nP69905\n", encoding="utf-8")
        assert download_structures.read_id_file(path) == ["P69905", "P04637"]

    def test_comma_annotation_stripped(self, tmp_path):
        path = tmp_path / "ids.txt"
        path.write_text("P69905, hemoglobin alpha\n", encoding="utf-8")
        assert download_structures.read_id_file(path) == ["P69905"]

    def test_empty_file_returns_empty_list(self, tmp_path):
        path = tmp_path / "ids.txt"
        path.write_text("", encoding="utf-8")
        assert download_structures.read_id_file(path) == []

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            download_structures.read_id_file(tmp_path / "nope.txt")


class TestAccessionHeuristic:
    """表头自动识别。"""

    @pytest.mark.parametrize(
        "text,expected",
        [
            ("P69905", True),
            ("Q9Y6K9", True),
            ("A0A024RBG1", True),
            ("id", False),
            ("uniprot_id", False),
            ("Accession", False),
            ("P1234", False),  # 太短
        ],
    )
    def test_looks_like_accession(self, text, expected):
        assert download_structures._looks_like_accession(text) is expected


class TestDownloadStructure:
    """下载行为。"""

    def test_binary_write_preserves_bytes(self, tmp_path):
        """必须按字节写入——文本模式会破坏 PDB 定宽列。"""
        payload = b"ATOM      1  N   MET A   1      1.0\nEND\n"
        session = _FakeSession(
            {
                "*": _FakeResponse(
                    content=payload, headers={"Content-Type": "text/plain"}
                )
            }
        )
        dest = tmp_path / "P69905.pdb"
        ok, detail = download_structures.download_structure(
            "https://example.com/x.pdb", dest, session=session
        )
        assert ok, detail
        assert dest.read_bytes() == payload

    def test_html_content_type_rejected(self, tmp_path):
        """HTML 错误页不能被当成 PDB 存下来。"""
        session = _FakeSession(
            {
                "*": _FakeResponse(
                    content=b"<!DOCTYPE html><html></html>",
                    headers={"Content-Type": "text/html"},
                )
            }
        )
        dest = tmp_path / "bad.pdb"
        ok, detail = download_structures.download_structure(
            "https://example.com/x.pdb", dest, session=session
        )
        assert not ok
        assert "内容类型" in detail or "Content-Type" in detail
        assert not dest.exists(), "被拒绝的响应不应留下文件"

    def test_empty_body_rejected(self, tmp_path):
        session = _FakeSession(
            {"*": _FakeResponse(content=b"", headers={"Content-Type": "text/plain"})}
        )
        dest = tmp_path / "empty.pdb"
        ok, _ = download_structures.download_structure(
            "https://example.com/x.pdb", dest, session=session
        )
        assert not ok

    def test_octet_stream_accepted(self, tmp_path):
        """部分 CDN 返回 octet-stream，应被接受。"""
        session = _FakeSession(
            {
                "*": _FakeResponse(
                    content=b"ATOM\nEND\n",
                    headers={"Content-Type": "application/octet-stream"},
                )
            }
        )
        dest = tmp_path / "ok.pdb"
        ok, _ = download_structures.download_structure(
            "https://example.com/x.pdb", dest, session=session
        )
        assert ok


class TestTargetName:
    """输出文件名规则。"""

    def test_default_is_accession_pdb(self):
        name = download_structures._target_name(
            "P69905", keep_original=False, pdb_url="https://x/AF-P69905-F1-model_v6.pdb"
        )
        assert name == "P69905.pdb"

    def test_keep_original_uses_url_basename(self):
        name = download_structures._target_name(
            "P69905", keep_original=True, pdb_url="https://x/AF-P69905-F1-model_v6.pdb"
        )
        assert name == "AF-P69905-F1-model_v6.pdb"

    def test_keep_original_falls_back_when_url_empty(self):
        name = download_structures._target_name("P69905", keep_original=True, pdb_url="")
        assert name == "P69905.pdb"


class TestMetadataExport:
    """元数据导出。"""

    def test_writes_csv_and_json(self, tmp_path):
        records = [alphafold_client.parse_prediction(_SAMPLE_ENTRY)]
        dest = tmp_path / "alphafold_metadata.csv"
        download_structures.write_metadata(records, dest)
        assert dest.is_file()
        assert dest.with_suffix(".json").is_file()

        import json

        payload = json.loads(dest.with_suffix(".json").read_text(encoding="utf-8"))
        assert payload[0]["accession"] == "P69905"

    def test_csv_contains_plddt_columns(self, tmp_path):
        records = [alphafold_client.parse_prediction(_SAMPLE_ENTRY)]
        dest = tmp_path / "meta.csv"
        download_structures.write_metadata(records, dest)
        text = dest.read_text(encoding="utf-8")
        assert "global_plddt" in text
        assert "fraction_plddt_confident" in text


class TestCliSurface:
    """命令行接口必须存在且可解析 --help。"""

    def test_parser_builds(self):
        parser = download_structures.build_parser()
        args = parser.parse_args(
            ["--ids", "P69905", "--output-dir", "out/", "--min-plddt", "70"]
        )
        assert args.ids == ["P69905"]
        assert args.min_plddt == 70
        assert args.output_dir == "out/"

    def test_metadata_only_flag(self):
        args = download_structures.build_parser().parse_args(
            ["--ids", "P69905", "--metadata-only"]
        )
        assert args.metadata_only is True

    def test_ids_is_required(self):
        with pytest.raises(SystemExit):
            download_structures.build_parser().parse_args([])

    def test_help_exits_cleanly(self, capsys):
        with pytest.raises(SystemExit) as exc:
            download_structures.build_parser().parse_args(["--help"])
        assert exc.value.code == 0
        assert "AlphaFold" in capsys.readouterr().out
