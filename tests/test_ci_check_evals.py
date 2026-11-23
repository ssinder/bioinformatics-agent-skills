"""评测用例体系的测试。

评测套件的价值在于**能不能挡住问题**，所以这里既检查现状合规，
也验证检查器本身能抓住问题（用故意写坏的临时文件）。
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "scripts" / "ci_check_evals.py"
SCHEMA = REPO_ROOT / "docs" / "evals.schema.json"

sys.path.insert(0, str(REPO_ROOT / "scripts"))


def _run(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(SCRIPT), *args],
        capture_output=True,
        text=True,
        timeout=120,
    )


class TestEvalFiles:
    """评测用例文件本身。"""

    def test_every_skill_has_evals(self) -> None:
        """每个技能都必须有评测用例。"""
        for skill_dir in sorted((REPO_ROOT / "skills").iterdir()):
            if not skill_dir.is_dir() or skill_dir.name.startswith((".", "_")):
                continue
            evals = skill_dir / "evals" / "evals.json"
            assert evals.is_file(), f"{skill_dir.name} 缺少 evals/evals.json"

    def test_checker_passes(self) -> None:
        """检查器当前应通过。"""
        proc = _run("--require-all")
        assert proc.returncode == 0, f"检查未通过:\n{proc.stdout}"

    def test_all_skills_covered(self) -> None:
        """报告应覆盖全部技能。"""
        proc = _run("--json")
        assert proc.returncode == 0, proc.stdout
        data = json.loads(proc.stdout)
        assert data["skills_with_evals"] == data["total_skills"]

    def test_minimum_case_count(self) -> None:
        """每个技能至少 3 个用例（官方建议的起步数量）。"""
        for skill_dir in sorted((REPO_ROOT / "skills").iterdir()):
            if not skill_dir.is_dir():
                continue
            data = json.loads(
                (skill_dir / "evals" / "evals.json").read_text(encoding="utf-8")
            )
            assert len(data["evals"]) >= 3, f"{skill_dir.name} 用例过少"


class TestNegativeCases:
    """近似负例是防止技能互相抢触发的关键。"""

    def test_every_skill_has_near_miss(self) -> None:
        """每个技能都必须至少有一个 should_not_trigger 近似负例。"""
        for skill_dir in sorted((REPO_ROOT / "skills").iterdir()):
            if not skill_dir.is_dir() or skill_dir.name.startswith((".", "_")):
                continue
            data = json.loads(
                (skill_dir / "evals" / "evals.json").read_text(encoding="utf-8")
            )
            negatives = [
                case
                for case in data["evals"]
                for assertion in case.get("assertions", [])
                if assertion.get("type") == "should_not_trigger"
            ]
            assert negatives, f"{skill_dir.name} 没有近似负例"

    def test_near_miss_points_at_a_sibling_skill(self) -> None:
        """近似负例应当指向某个兄弟技能，而不是空泛地说"不该触发"。

        否则代理拿到负例请求时不知道该转给谁，负例就失去意义。
        """
        skill_names = {
            d.name for d in (REPO_ROOT / "skills").iterdir() if d.is_dir()
        }
        for skill_dir in sorted((REPO_ROOT / "skills").iterdir()):
            if not skill_dir.is_dir() or skill_dir.name.startswith((".", "_")):
                continue
            data = json.loads(
                (skill_dir / "evals" / "evals.json").read_text(encoding="utf-8")
            )
            for case in data["evals"]:
                if not any(
                    a.get("type") == "should_not_trigger"
                    for a in case.get("assertions", [])
                ):
                    continue
                expected = str(case.get("expected_output", "")).lower()
                mentioned = {n for n in skill_names if n in expected}
                assert mentioned or "skill" in expected, (
                    f"{skill_dir.name} 用例 {case.get('id')} 的负例既没指向具体技能，"
                    f"也没说明应转向何处"
                )


class TestCaseQuality:
    """用例本身的质量。"""

    def test_prompts_are_specific(self) -> None:
        """prompt 应当包含具体场景，而不是抽象需求。"""
        for skill_dir in sorted((REPO_ROOT / "skills").iterdir()):
            if not skill_dir.is_dir() or skill_dir.name.startswith((".", "_")):
                continue
            data = json.loads(
                (skill_dir / "evals" / "evals.json").read_text(encoding="utf-8")
            )
            for case in data["evals"]:
                prompt = str(case.get("prompt", ""))
                assert len(prompt) >= 30, (
                    f"{skill_dir.name} 用例 {case.get('id')} 的 prompt 过短"
                )

    def test_expected_outputs_are_checkable(self) -> None:
        """expected_output 要具体到可核查，不能是空话。"""
        vague = ("give the right answer", "正确回答", "输出正确结果")
        for skill_dir in sorted((REPO_ROOT / "skills").iterdir()):
            if not skill_dir.is_dir() or skill_dir.name.startswith((".", "_")):
                continue
            data = json.loads(
                (skill_dir / "evals" / "evals.json").read_text(encoding="utf-8")
            )
            for case in data["evals"]:
                expected = str(case.get("expected_output", ""))
                assert len(expected) >= 40, (
                    f"{skill_dir.name} 用例 {case.get('id')} 的 expected_output 过短"
                )
                lowered = expected.lower()
                for phrase in vague:
                    assert phrase.lower() not in lowered, (
                        f"{skill_dir.name} 用例 {case.get('id')} 的 expected_output 过于笼统"
                    )

    def test_skill_name_matches_directory(self) -> None:
        """skill_name 必须与目录名一致。"""
        for skill_dir in sorted((REPO_ROOT / "skills").iterdir()):
            if not skill_dir.is_dir() or skill_dir.name.startswith((".", "_")):
                continue
            data = json.loads(
                (skill_dir / "evals" / "evals.json").read_text(encoding="utf-8")
            )
            assert data["skill_name"] == skill_dir.name


class TestSchema:
    """schema 文件本身。"""

    def test_schema_exists_and_is_valid_json(self) -> None:
        """schema 必须是合法 JSON。"""
        assert SCHEMA.is_file()
        data = json.loads(SCHEMA.read_text(encoding="utf-8"))
        assert data.get("type") == "object"
        assert "skill_name" in data.get("properties", {})
        assert "evals" in data.get("properties", {})

    def test_schema_declares_both_assertion_types(self) -> None:
        """schema 必须同时声明正例与负例断言类型。"""
        data = json.loads(SCHEMA.read_text(encoding="utf-8"))
        props = data["properties"]["evals"]["items"]["properties"]
        assertion_schema = props["assertions"]["items"]["properties"]
        enum = assertion_schema["type"]["enum"]
        assert "should_trigger" in enum
        assert "should_not_trigger" in enum
