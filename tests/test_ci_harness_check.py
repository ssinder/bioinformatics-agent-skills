"""跨宿主兼容性检查器的测试。

这组测试守护的是一个**静默失败**：技能装上去了，但宿主不加载它。
没有报错、没有日志，只有一句"这个技能好像不灵"。
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "scripts" / "ci_harness_check.py"

sys.path.insert(0, str(REPO_ROOT / "scripts"))


def _run(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(SCRIPT), *args],
        capture_output=True,
        text=True,
        timeout=120,
    )


class TestHarnessCheck:
    """检查器本身的行为。"""

    def test_repo_passes(self) -> None:
        """仓库当前状态必须通过检查。"""
        proc = _run()
        assert proc.returncode == 0, f"检查未通过:\n{proc.stdout}\n{proc.stderr}"

    def test_reports_skill_count(self) -> None:
        """输出应说明检查了多少个技能。"""
        proc = _run()
        assert "17 个技能" in proc.stdout

    def test_json_output_is_machine_readable(self) -> None:
        """--json 应输出可解析的 JSON，供 CI 消费。"""
        proc = _run("--json")
        assert proc.returncode == 0, proc.stdout
        data = json.loads(proc.stdout)
        assert data["total"] == 17
        assert data["problems"] == []


class TestDescriptionsAreHostSafe:
    """每个技能的描述都要在各宿主的限制内。"""

    def test_no_description_exceeds_host_limit(self) -> None:
        """description 不得超过 1,024 字符。

        超限时部分宿主（copilot-cli）会**整个技能不加载**——
        不是截断显示，而是直接消失。
        """
        import yaml

        limit = 1024
        offenders = []
        for skill_md in sorted((REPO_ROOT / "skills").glob("*/SKILL.md")):
            text = skill_md.read_text(encoding="utf-8")
            fm = yaml.safe_load(text.split("---", 2)[1])
            desc = str(fm.get("description", ""))
            if len(desc) > limit:
                offenders.append((skill_md.parent.name, len(desc)))
        assert not offenders, f"description 超限: {offenders}"

    def test_trigger_clause_comes_first(self) -> None:
        """触发用例必须在描述开头。

        宿主为了控制上下文预算会从**尾部**截断描述，
        把触发语放前面才不会被截掉。
        """
        import yaml

        offenders = []
        for skill_md in sorted((REPO_ROOT / "skills").glob("*/SKILL.md")):
            text = skill_md.read_text(encoding="utf-8")
            fm = yaml.safe_load(text.split("---", 2)[1])
            desc = str(fm.get("description", "")).strip()
            if not desc.lower().startswith(("use when", "this skill should be used when")):
                offenders.append(skill_md.parent.name)
        assert not offenders, f"触发语不在开头: {offenders}"

    def test_frontmatter_only_uses_spec_fields(self) -> None:
        """frontmatter 只能含规范定义的字段。

        claude.ai 上传与 Skills API 只接受这六个字段，
        多出的字段会让打包直接失败。
        """
        import yaml

        allowed = {
            "name",
            "description",
            "license",
            "compatibility",
            "metadata",
            "allowed-tools",
        }
        offenders = []
        for skill_md in sorted((REPO_ROOT / "skills").glob("*/SKILL.md")):
            text = skill_md.read_text(encoding="utf-8")
            fm = yaml.safe_load(text.split("---", 2)[1]) or {}
            extra = set(fm.keys()) - allowed
            if extra:
                offenders.append((skill_md.parent.name, sorted(extra)))
        assert not offenders, f"frontmatter 含规范外字段: {offenders}"

    def test_metadata_values_are_strings(self) -> None:
        """metadata 必须是 string→string 映射（规范要求）。"""
        import yaml

        offenders = []
        for skill_md in sorted((REPO_ROOT / "skills").glob("*/SKILL.md")):
            text = skill_md.read_text(encoding="utf-8")
            fm = yaml.safe_load(text.split("---", 2)[1]) or {}
            metadata = fm.get("metadata")
            if not isinstance(metadata, dict):
                continue
            for key, value in metadata.items():
                if not isinstance(value, str):
                    offenders.append((skill_md.parent.name, key, type(value).__name__))
        assert not offenders, f"metadata 值不是字符串: {offenders}"


class TestListingBudget:
    """技能列表的常驻成本。"""

    def test_total_description_budget(self) -> None:
        """全部技能的 name+description 必须在 8,000 字符预算内。

        超出后宿主会从最少使用的技能开始丢弃描述，
        那些技能只剩裸名字而永远无法被触发。
        """
        import yaml

        total = 0
        for skill_md in sorted((REPO_ROOT / "skills").glob("*/SKILL.md")):
            text = skill_md.read_text(encoding="utf-8")
            fm = yaml.safe_load(text.split("---", 2)[1]) or {}
            total += len(str(fm.get("name", ""))) + len(str(fm.get("description", "")))
        assert total <= 8000, f"技能列表成本 {total} 超出 8,000 字符预算"

    def test_budget_headroom(self) -> None:
        """应留有余量给用户自己安装的其他技能。

        硬预算是 8,000 字符；软限设在 7,000（87.5%）——
        全量 17 个技能装齐后仍有约 1,000 字符余量，
        而按 skill-guides 推荐只装 2–3 个时仅占约 15%。
        """
        import yaml

        total = 0
        for skill_md in sorted((REPO_ROOT / "skills").glob("*/SKILL.md")):
            text = skill_md.read_text(encoding="utf-8")
            fm = yaml.safe_load(text.split("---", 2)[1]) or {}
            total += len(str(fm.get("name", ""))) + len(str(fm.get("description", "")))
        assert total <= 7000, f"技能列表成本 {total} 超出建议上限 7,000 字符"
