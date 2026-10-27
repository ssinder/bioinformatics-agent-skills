"""仓库一致性检查器的回归测试。

背景
----
``scripts/ci_consistency_check.py`` 负责守住几处"手工维护的计数"
（README 徽章上的技能数、技能清单表格、REGISTRY 与磁盘的对应关系等）。
这些数字一旦漂移，CI 不会报错，只有读者会发现文档在撒谎。

本模块锁定其中最容易复发的一项：**文档里声明的技能总数必须与实际一致**。
``docs/skills.md`` 曾长期写着「共 8 个技能」而实际已有 17 个 —— 它是
早期 ``build_catalog.py`` 的产物，被 ``docs/skills/index.md`` 取代后
忘了删除，没有任何检查会报错。

本模块只依赖标准库 + pytest。
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "scripts" / "ci_consistency_check.py"


def _load_module():
    """以文件路径加载被测脚本（脚本名不是合法模块名）。"""
    spec = importlib.util.spec_from_file_location("ci_consistency_check", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules["ci_consistency_check"] = module
    spec.loader.exec_module(module)
    return module


mod = _load_module()


class TestSkillCountDetection:
    """文档里的技能总数声明必须被识别出来。"""

    def test_detects_stale_chinese_count(self, tmp_path: Path) -> None:
        """过时计数（中文写法）应被检出。"""
        doc = tmp_path / "stale.md"
        doc.write_text("本项目共 **8** 个技能,按分类分组。\n", encoding="utf-8")

        problems: list[str] = []
        mod._check_doc_skill_counts(problems, 17, paths=[doc])

        assert problems, "过时的技能数声明应被检出"
        assert "8" in problems[0]
        assert "17" in problems[0]

    def test_detects_stale_plain_chinese_count(self, tmp_path: Path) -> None:
        """没有加粗标记的中文计数也应被检出。"""
        doc = tmp_path / "stale.md"
        doc.write_text("本项目共 8 个技能。\n", encoding="utf-8")

        problems: list[str] = []
        mod._check_doc_skill_counts(problems, 17, paths=[doc])

        assert problems, "无加粗标记的过时计数也应被检出"

    def test_detects_stale_english_count(self, tmp_path: Path) -> None:
        """英文写法的过时计数同样应被检出。"""
        doc = tmp_path / "stale.md"
        doc.write_text("We provide 8 skills in total.\n", encoding="utf-8")

        problems: list[str] = []
        mod._check_doc_skill_counts(problems, 17, paths=[doc])

        assert problems, "英文的过时计数也应被检出"

    def test_correct_count_passes(self, tmp_path: Path) -> None:
        """与实际一致的计数不应报错。"""
        doc = tmp_path / "ok.md"
        doc.write_text("本项目共 **17** 个技能。\n", encoding="utf-8")

        problems: list[str] = []
        mod._check_doc_skill_counts(problems, 17, paths=[doc])

        assert not problems, f"正确计数不应报错，实际报了: {problems}"


class TestNoFalsePositives:
    """举例中的数字不应被误判为技能总数。"""

    def test_example_counts_are_not_flagged(self, tmp_path: Path) -> None:
        """「这 3 个技能应先做质控」是举例，不是总数声明。"""
        doc = tmp_path / "example.md"
        doc.write_text(
            "这 3 个技能应先做质控再进入下游分析。\n\n"
            "对于 5 skills 场景，先确认输入格式。\n",
            encoding="utf-8",
        )

        problems: list[str] = []
        mod._check_doc_skill_counts(problems, 17, paths=[doc])

        assert not problems, f"举例中的数字不应被误判，实际报了: {problems}"

    def test_skill_listing_numbers_are_not_flagged(self, tmp_path: Path) -> None:
        """技能表格里的行号/序号不应被误判。"""
        doc = tmp_path / "table.md"
        doc.write_text(
            "| Skill | Scripts |\n|---|---|\n"
            "| `bulk-rnaseq` | 16 scripts |\n",
            encoding="utf-8",
        )

        problems: list[str] = []
        mod._check_doc_skill_counts(problems, 17, paths=[doc])

        assert not problems, f"表格中的脚本数不应被误判，实际报了: {problems}"


class TestContextHeuristic:
    """上下文判定以数字位置为基准，而非整个匹配的位置。"""

    def test_context_uses_number_position(self) -> None:
        """「共 **8** 个技能」中，匹配起点在「共」上，不能以它为基准。

        这是本检查最初的 bug：以匹配起点为基准时，前向窗口里找不到
        「共」，判定永远为假，检查形同虚设。
        """
        text = "本项目共 **8** 个技能,按分类分组。"
        # 「共」的位置
        assert mod._is_skill_count_context(text, text.index("8")) is True
        # 「共」的位置（不是数字位置）
        assert mod._is_skill_count_context(text, text.index("共")) is False

    def test_context_requires_total_qualifier(self) -> None:
        """没有总数量措辞时不应判为总数语境。"""
        text = "在某些流程中，8 个技能可以并行执行。"
        assert mod._is_skill_count_context(text, text.index("8")) is False

    def test_context_requires_skill_word_nearby(self) -> None:
        """数字后面必须紧跟「个技能 / skills」。"""
        text = "本项目共 8 个模块。"
        assert mod._is_skill_count_context(text, text.index("8")) is False


class TestCurrentRepositoryIsClean:
    """当前仓库自身必须通过这项检查。"""

    def test_repository_has_no_stale_counts(self) -> None:
        problems: list[str] = []
        mod._check_doc_skill_counts(problems, len(mod._skill_dirs()))
        assert not problems, f"仓库内存在过时的技能数声明: {problems}"

    def test_script_reports_success_on_current_repo(self) -> None:
        """整脚本在当前仓库上必须返回 0。"""
        assert mod.main() == 0
