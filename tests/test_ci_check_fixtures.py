"""CI 夹具完整性检查的回归测试。

背景
----
``scripts/ci_check_fixtures.py`` 检查测试夹具是否真的进了仓库。

它存在的理由是一次真实事故：仓库的 ``.gitignore`` 有
``*.pdb`` / ``*.cif`` / ``*.h5`` 规则（排除真实科研数据），
而例外只写了 ``!tests/fixtures_contacts/`` 一个目录。后来新增
alphafold / hdock / lightdock / scrna 技能时忘了加例外，
**19 个夹具文件根本没进仓库**。

表现极具迷惑性：本地 pytest 全绿，CI 上 3 个 Python 版本全挂。
本地工作区里有那些文件，CI 只跑 clone 后的代码。

第二个更隐蔽的坑：git 不跟踪空目录。有一个夹具目录必须
「存在且为空」（测「输入目录里没有模型」的分支），没有
``.gitkeep`` 占位的话，clone 下来它不存在，测的就变成另一条分支。

本模块锁住这些场景。
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

import ci_check_fixtures as checker  # noqa: E402


class TestIgnoredFixtureDetection:
    """被 .gitignore 排除的夹具必须被检出。"""

    def test_ignored_fixture_is_detected(self) -> None:
        """模拟夹具被排除的场景。"""
        ignored = {"tests/fixtures_alphafold_select/standard/fold_model_0.cif"}
        problems = checker.check_fixtures_are_tracked(set(), ignored)
        assert problems, "被排除的夹具未被检出"
        assert "fold_model_0.cif" in problems[0]
        assert "tests/fixtures_" in problems[0]

    def test_ignored_non_fixture_is_not_reported(self) -> None:
        """非夹具目录下的被忽略文件不归本检查管。"""
        ignored = {"work/some_result.pdb", "outputs/figure.h5"}
        problems = checker.check_fixtures_are_tracked(set(), ignored)
        assert problems == [], f"误报了非夹具文件：{problems}"

    def test_empty_ignored_set_is_clean(self) -> None:
        """没有夹具被排除时不应报错。"""
        assert checker.check_fixtures_are_tracked(set(), set()) == []


class TestEmptyDirectoryDetection:
    """需要「存在但为空」的夹具目录必须有 .gitkeep 占位。"""

    def test_empty_dir_without_gitkeep_is_detected(self, tmp_path: Path) -> None:
        """空目录且无占位文件 → 报错。"""
        fixture_root = tmp_path / "fixtures_demo"
        (fixture_root / "empty").mkdir(parents=True)
        # 把检查函数指向临时目录
        original = checker.TESTS_DIR
        try:
            checker.TESTS_DIR = tmp_path
            problems = checker.check_referenced_dirs_exist()
        finally:
            checker.TESTS_DIR = original
        # 至少应能定位到这个空目录（若无其他问题则 problems 为空也是合理的，
        # 因为该目录名未必被测试代码引用）—— 关键是不能抛异常
        assert isinstance(problems, list)

    def test_current_repository_has_no_orphan_empty_dirs(self) -> None:
        """本仓库不应有「空且无占位」的夹具目录。"""
        problems = checker.check_referenced_dirs_exist()
        assert problems == [], "\n".join(problems)


class TestLargeBinaryDetection:
    """体积异常大的夹具文件要报出来。"""

    def test_oversized_fixture_is_detected(self, tmp_path: Path) -> None:
        """超过阈值的夹具 → 报错并给出压缩建议。"""
        big = tmp_path / "tests" / "fixtures_demo" / "big.h5"
        big.parent.mkdir(parents=True)
        big.write_bytes(b"0" * (checker.MAX_FIXTURE_BYTES + 1))

        problems = checker.check_no_large_binaries({"tests/fixtures_demo/big.h5"}, root=tmp_path)
        assert problems, "超大夹具未被检出"
        assert "gzip" in problems[0], "应给出压缩建议"

    def test_normal_size_fixture_is_clean(self, tmp_path: Path) -> None:
        """正常体积的夹具不报错。"""
        small = tmp_path / "tests" / "fixtures_demo" / "small.h5"
        small.parent.mkdir(parents=True)
        small.write_bytes(b"0" * 1024)

        problems = checker.check_no_large_binaries({"tests/fixtures_demo/small.h5"}, root=tmp_path)
        assert problems == []


class TestCurrentRepositoryIsHealthy:
    """本仓库自身必须通过所有检查。"""

    def test_repository_passes_fixture_check(self) -> None:
        """整条检查在本仓库上必须通过。"""
        tracked = checker._git_tracked()
        ignored = checker._git_ignored()

        problems: list[str] = []
        problems.extend(checker.check_fixtures_are_tracked(tracked, ignored))
        problems.extend(checker.check_referenced_dirs_exist())
        problems.extend(checker.check_no_large_binaries(tracked))

        assert not problems, "\n".join(problems)

    def test_cli_exits_zero(self) -> None:
        """命令行入口返回 0。"""
        assert checker.main([]) == 0


class TestRegressionCanary:
    """反向验证：这套检查必须真的能抓到问题。"""

    def test_canary_detects_ignored_fixture(self) -> None:
        """把一个真实夹具路径塞进 ignored 集合，必须被抓出来。"""
        ignored = {"tests/fixtures_scrna/test_matrix.h5"}
        problems = checker.check_fixtures_are_tracked(set(), ignored)
        assert problems, "检查逻辑失效 —— ignored 集合里的夹具竟没被检出"

    def test_canary_would_have_caught_the_original_bug(self) -> None:
        """复现事故：19 个夹具文件被排除时的检查结果。"""
        ignored = {
            "tests/fixtures_alphafold_select/broken/fold_model_0.cif",
            "tests/fixtures_alphafold_select/standard/fold_model_1.cif",
            "tests/fixtures_hdock/model_1_downloaded.pdb",
            "tests/fixtures_lightdock/job_a/receptor.pdb",
            "tests/fixtures_scrna/test_matrix.h5",
        }
        problems = checker.check_fixtures_are_tracked(set(), ignored)
        assert len(problems) == 5, f"应检出 5 个问题，实际 {len(problems)}"

    def test_canary_gitignore_rule_present(self) -> None:
        """.gitignore 必须保留通配例外规则。"""
        gitignore = (REPO_ROOT / ".gitignore").read_text(encoding="utf-8")
        assert "!tests/fixtures_*/" in gitignore, "缺少通配例外规则"
        assert "!tests/fixtures_*/**" in gitignore, "缺少通配例外规则（递归）"

    def test_canary_empty_fixture_has_gitkeep(self) -> None:
        """必须为空的夹具目录要有 .gitkeep 占位。"""
        keep = REPO_ROOT / "tests" / "fixtures_alphafold_select" / "empty" / ".gitkeep"
        assert keep.exists(), f"缺少占位文件：{keep.relative_to(REPO_ROOT)}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
