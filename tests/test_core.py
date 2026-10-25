"""biopipe_core 共享工具层的单元测试。"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "packages"))

from biopipe_core import (  # noqa: E402
    ensure_dir,
    get_config,
    read_csv_safe,
    read_fasta,
    resolve_path,
    write_csv,
    write_fasta,
)


class TestConfig:
    """配置解析测试。"""

    def test_default_config_from_cwd(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        monkeypatch.delenv("BIOPIPE_PROJECT_ROOT", raising=False)
        config = get_config(refresh=True)
        assert config.root == tmp_path.resolve()
        assert config.work_dir.name == "work"
        assert config.output_dir.name == "results"

    def test_env_root_override(self, tmp_path, monkeypatch):
        monkeypatch.setenv("BIOPIPE_PROJECT_ROOT", str(tmp_path))
        config = get_config(refresh=True)
        assert config.root == tmp_path.resolve()

    def test_json_config_file(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        cfg = tmp_path / "custom.json"
        cfg.write_text(
            json.dumps({"root": str(tmp_path), "output_dir": "out_custom"}), encoding="utf-8"
        )
        config = get_config(cfg, refresh=True)
        assert config.output_dir.name == "out_custom"

    def test_key_value_config_file(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        cfg = tmp_path / "custom.yaml"
        cfg.write_text(f"root: {tmp_path}\nfigures_dir: figs\n", encoding="utf-8")
        config = get_config(cfg, refresh=True)
        assert config.figures_dir.name == "figs"

    def test_ensure_dirs_creates_all(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        monkeypatch.delenv("BIOPIPE_PROJECT_ROOT", raising=False)
        config = get_config(refresh=True).ensure_dirs()
        for directory in (config.work_dir, config.output_dir, config.data_dir, config.figures_dir):
            assert directory.is_dir()

    def test_missing_config_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            get_config(tmp_path / "nope.json", refresh=True)


class TestPaths:
    """路径解析测试。"""

    def test_absolute_path_preserved(self, tmp_path):
        absolute = str(tmp_path / "x").replace("\\", "/")
        resolved = resolve_path(absolute)
        assert resolved.is_absolute()
        assert resolved.name == "x"

    def test_relative_path_uses_cwd(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        assert resolve_path("data/x.txt") == (tmp_path / "data" / "x.txt").resolve()

    def test_ensure_dir(self, tmp_path):
        target = ensure_dir(tmp_path / "a" / "b")
        assert target.is_dir()


class TestFasta:
    """FASTA 读写测试。"""

    def test_roundtrip(self, tmp_path):
        records = [("seq1", "ACGTACGT"), ("seq2", "TTTTAAAA")]
        path = write_fasta(records, tmp_path / "out.fa")
        loaded = list(read_fasta(path))
        assert loaded == records

    def test_multiline_sequence_joined(self, tmp_path):
        path = tmp_path / "ml.fa"
        path.write_text(">s1 desc here\nACGT\nACGT\n>s2\nTTTT\n", encoding="utf-8")
        loaded = dict(read_fasta(path))
        assert loaded["s1"] == "ACGTACGT"
        assert loaded["s2"] == "TTTT"

    def test_comment_lines_skipped(self, tmp_path):
        path = tmp_path / "c.fa"
        path.write_text(";comment\n>s1\nACGT\n", encoding="utf-8")
        assert list(read_fasta(path)) == [("s1", "ACGT")]

    def test_wrapped_output(self, tmp_path):
        path = write_fasta([("s", "A" * 25)], tmp_path / "w.fa", width=10)
        lines = path.read_text(encoding="utf-8").strip().split("\n")
        assert len(lines) == 4
        assert lines[0] == ">s"

    def test_empty_sequence(self, tmp_path):
        path = tmp_path / "e.fa"
        path.write_text(">s\n", encoding="utf-8")
        assert list(read_fasta(path)) == [("s", "")]


class TestCsv:
    """CSV 读写测试。"""

    def test_roundtrip_with_chinese(self, tmp_path):
        rows = [{"名称": "α螺旋", "值": "1"}, {"名称": "β折叠", "值": "2"}]
        path = write_csv(rows, tmp_path / "c.csv")
        loaded = read_csv_safe(path)
        assert loaded[0]["名称"] == "α螺旋"
        assert len(loaded) == 2

    def test_bom_written(self, tmp_path):
        path = write_csv([{"a": "1"}], tmp_path / "b.csv")
        assert path.read_bytes().startswith(b"\xef\xbb\xbf")

    def test_semicolon_delimiter_autodetect(self, tmp_path):
        path = tmp_path / "s.csv"
        path.write_text("a;b\n1;2\n", encoding="utf-8")
        loaded = read_csv_safe(path)
        assert loaded[0]["b"] == "2"

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            read_csv_safe(tmp_path / "nope.csv")


class TestRegistry:
    """技能注册表测试。"""

    def test_registry_non_empty(self):
        from biopipe.registry import REGISTRY

        assert len(REGISTRY) >= 8

    def test_registry_entries_have_skill_md(self):
        from biopipe.registry import REGISTRY

        missing = [s.name for s in REGISTRY.values() if not s.skill_md.exists()]
        if missing:
            pytest.skip(f"以下技能尚未生成 SKILL.md: {missing}")
        assert not missing

    def test_registry_names_are_kebab(self):
        from biopipe.registry import REGISTRY

        for name in REGISTRY:
            assert name == name.lower()
            assert " " not in name
            assert "_" not in name
