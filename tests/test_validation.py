"""技能规范校验器的自测试。

用一个临时构造的"坏技能"验证校验器确实能抓出问题。
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from validate_skills import validate_skill  # noqa: E402


def _write(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def _good_skill(root: Path, name: str = "good-skill") -> Path:
    skill_dir = root / name
    _write(
        skill_dir / "SKILL.md",
        f"""---
name: {name}
description: This skill should be used when a user needs to analyse bioinformatics data with reproducible pipelines and validated scientific methods.
license: MIT
metadata:
  version: "1.0"
  last-reviewed: "2026-10-25"
  skill-author: sinder
---

# Good Skill

## Overview
A skill that satisfies every validation rule.

## When to Use This Skill
- When analysing bio data

## Quick Start
Run it.

## Workflow
1. Step one

## Parameters
| Name | Description |

## Outputs
| File | Description |

## Validation
Check outputs.

## References
See references/methodology.md
""",
    )
    _write(
        skill_dir / "scripts" / "run.py",
        '''#!/usr/bin/env python3
"""Test script."""
import argparse


def main(args):
    print(args.input)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input")
    main(parser.parse_args())
''',
    )
    _write(skill_dir / "references" / "methodology.md", "# Methodology\n" + "content " * 50)
    _write(skill_dir / "requirements.txt", "requests>=2.28\n")
    return skill_dir


class TestValidator:
    """校验器行为测试。"""

    def test_good_skill_passes(self, tmp_path):
        skill_dir = _good_skill(tmp_path)
        result = validate_skill(skill_dir)
        assert result.ok, [str(i) for i in result.issues]

    def test_missing_skill_md_fails(self, tmp_path):
        skill_dir = tmp_path / "no-md"
        skill_dir.mkdir()
        result = validate_skill(skill_dir)
        assert not result.ok
        assert any(i.code == "no-skill-md" for i in result.errors)

    def test_name_mismatch_detected(self, tmp_path):
        skill_dir = _good_skill(tmp_path)
        text = (skill_dir / "SKILL.md").read_text(encoding="utf-8")
        text = text.replace("name: good-skill", "name: other-name")
        (skill_dir / "SKILL.md").write_text(text, encoding="utf-8")
        result = validate_skill(skill_dir)
        assert not result.ok
        assert any(i.code == "name-mismatch" for i in result.errors)

    def test_second_person_description_rejected(self, tmp_path):
        skill_dir = _good_skill(tmp_path)
        text = (skill_dir / "SKILL.md").read_text(encoding="utf-8")
        text = text.replace(
            "This skill should be used when a user needs to analyse bioinformatics data with reproducible pipelines and validated scientific methods.",
            "Use this skill when you need to analyse bioinformatics data with reproducible pipelines and validated scientific methods.",
        )
        (skill_dir / "SKILL.md").write_text(text, encoding="utf-8")
        result = validate_skill(skill_dir)
        assert not result.ok
        assert any(i.code == "desc-second-person" for i in result.errors)

    def test_missing_section_detected(self, tmp_path):
        skill_dir = _good_skill(tmp_path)
        text = (skill_dir / "SKILL.md").read_text(encoding="utf-8")
        text = text.replace("## Validation\nCheck outputs.\n", "")
        (skill_dir / "SKILL.md").write_text(text, encoding="utf-8")
        result = validate_skill(skill_dir)
        assert not result.ok
        assert any(i.code == "missing-section" for i in result.errors)

    def test_hardcoded_path_detected(self, tmp_path):
        skill_dir = _good_skill(tmp_path)
        script = skill_dir / "scripts" / "run.py"
        script.write_text(
            'BASE = "/Users/jung/Desktop/trae/tool/data"\n'
            "import argparse\n"
            "def main(a):\n    print(BASE)\n"
            'if __name__ == "__main__":\n    main(None)\n',
            encoding="utf-8",
        )
        result = validate_skill(skill_dir)
        assert not result.ok
        assert any(i.code == "script-hardcoded-path" for i in result.errors)

    def test_syntax_error_detected(self, tmp_path):
        skill_dir = _good_skill(tmp_path)
        (skill_dir / "scripts" / "run.py").write_text("def broken(:\n", encoding="utf-8")
        result = validate_skill(skill_dir)
        assert not result.ok
        assert any(i.code == "script-syntax" for i in result.errors)

    def test_missing_main_detected(self, tmp_path):
        skill_dir = _good_skill(tmp_path)
        (skill_dir / "scripts" / "run.py").write_text(
            "import argparse\ndef main(a):\n    pass\n", encoding="utf-8"
        )
        result = validate_skill(skill_dir)
        assert not result.ok
        assert any(i.code == "script-no-main" for i in result.errors)


class TestRealSkills:
    """仓库内真实技能必须通过校验。"""

    def test_all_repo_skills_valid(self):
        from validate_skills import SKILLS_DIR, validate_skill

        skills = [p for p in SKILLS_DIR.iterdir() if p.is_dir() and not p.name.startswith(".")]
        if not skills:
            pytest.skip("no skills present")
        for skill_dir in skills:
            result = validate_skill(skill_dir)
            assert result.ok, f"{skill_dir.name}: {[str(i) for i in result.errors]}"
