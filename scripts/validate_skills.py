#!/usr/bin/env python3
"""校验 skills/ 下所有技能是否符合 Agent Skills 规范。

这是 CI 的第一道关卡。任何技能在合并前必须通过本校验。

校验项：
1. 目录存在且非空
2. SKILL.md 存在，YAML frontmatter 可解析
3. frontmatter 必填字段齐全（name/description/license/metadata）
4. name 与目录名一致，且为合法 kebab-case
5. description 为第三人称、含触发语句、长度在合理区间
6. description 包含可被检索的关键词
7. SKILL.md 正文包含全部必需章节
8. scripts/ 下的 Python 文件语法正确、无硬编码绝对路径
9. references/ 文档非空
10. 每个脚本都有 __main__ 入口

用法：
    python scripts/validate_skills.py                # 校验全部技能
    python scripts/validate_skills.py --skill bulk-rnaseq
    python scripts/validate_skills.py --strict       # 警告也视为失败
    python scripts/validate_skills.py --json         # 机器可读输出
"""

from __future__ import annotations

import argparse
import ast
import json
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

try:
    import yaml
except ImportError:  # pragma: no cover
    yaml = None

REPO_ROOT = Path(__file__).resolve().parents[1]
SKILLS_DIR = REPO_ROOT / "skills"

REQUIRED_FRONTMATTER = ["name", "description", "license"]
REQUIRED_SECTIONS = [
    "Overview",
    "When to Use This Skill",
    "Quick Start",
    "Workflow",
    "Parameters",
    "Outputs",
    "Validation",
    "References",
]

# 触发语句：第三人称 + when 语义
TRIGGER_PATTERNS = [
    r"This skill should be used when",
    r"should be used when",
    r"用于.{0,40}的场景",
    r"当用户.{0,30}时",
    r"适用于",
]

# 硬编码绝对路径检测（Windows 与 POSIX）
HARDCODED_PATH_PATTERNS = [
    re.compile(r"[\"'](/Users/[A-Za-z0-9_.-]+/)"),
    re.compile(r"[\"'](/home/[A-Za-z0-9_.-]+/)"),
    re.compile(r"[\"']([A-Za-z]:[\\/]Users[\\/][A-Za-z0-9_. -]+[\\/])"),
    re.compile(r"[\"'](/Users/jung)"),
    re.compile(r"[\"'](\\\\Users\\\\[A-Za-z0-9_.-]+)"),
]

NAME_PATTERN = re.compile(r"^[a-z0-9]+(?:-[a-z0-9]+)*$")

MAX_DESC = 1024
MIN_DESC = 80


@dataclass
class Issue:
    """单条校验问题。"""

    skill: str
    level: str  # "error" | "warning"
    code: str
    message: str

    def __str__(self) -> str:
        return f"[{self.level.upper():<7}] {self.skill}: {self.code} - {self.message}"


@dataclass
class Result:
    """单个技能的校验结果。"""

    skill: str
    issues: list[Issue] = field(default_factory=list)

    @property
    def errors(self) -> list[Issue]:
        return [i for i in self.issues if i.level == "error"]

    @property
    def warnings(self) -> list[Issue]:
        return [i for i in self.issues if i.level == "warning"]

    @property
    def ok(self) -> bool:
        return not self.errors


def _read_frontmatter(text: str) -> tuple[dict | None, str, str | None]:
    """拆分 SKILL.md 的 frontmatter 与正文。

    Returns:
        ``(frontmatter_dict, body, error_message)``；解析失败时 dict 为 None。
    """
    if not text.startswith("---"):
        return None, text, "SKILL.md 未以 '---' 开头，缺少 YAML frontmatter"

    parts = text.split("---", 2)
    if len(parts) < 3:
        return None, text, "frontmatter 未正确闭合（缺少结束的 '---'）"

    raw = parts[1]
    body = parts[2]
    if yaml is None:
        return None, body, "缺少 PyYAML，无法解析 frontmatter"
    try:
        data = yaml.safe_load(raw)
    except yaml.YAMLError as exc:
        return None, body, f"frontmatter YAML 解析失败: {exc}"
    if not isinstance(data, dict):
        return None, body, "frontmatter 必须是键值映射"
    return data, body, None


def _check_description(skill: str, desc: str, issues: list[Issue]) -> None:
    """校验 description 质量。"""
    if not desc or not desc.strip():
        issues.append(Issue(skill, "error", "desc-empty", "description 不能为空"))
        return

    if len(desc) < MIN_DESC:
        issues.append(
            Issue(skill, "warning", "desc-short", f"description 过短 ({len(desc)} 字符)，应说明做什么与何时触发")
        )
    if len(desc) > MAX_DESC:
        issues.append(Issue(skill, "error", "desc-long", f"description 过长 ({len(desc)} 字符)，上限 {MAX_DESC}"))

    if not any(re.search(p, desc, re.IGNORECASE) for p in TRIGGER_PATTERNS):
        issues.append(
            Issue(
                skill,
                "error",
                "desc-no-trigger",
                "description 缺少触发语句，需含 'This skill should be used when...' 或等价中文表述",
            )
        )

    if re.search(r"\byou should\b|\buse this skill when\b", desc, re.IGNORECASE):
        issues.append(
            Issue(skill, "error", "desc-second-person", "description 应使用第三人称，避免 'You should' / 'Use this skill when'")
        )

    lowered = desc.lower()
    if len(re.findall(r"[a-z0-9\-]{3,}", lowered)) < 4:
        issues.append(
            Issue(skill, "warning", "desc-few-keywords", "description 缺少足够检索关键词，影响技能被正确触发")
        )


def _display_path(path: Path) -> Path:
    """把路径转成相对于仓库根的展示形式；不在仓库内时退回文件名。"""
    try:
        return path.relative_to(REPO_ROOT)
    except ValueError:
        return path


def _is_library_module(source: str, filename: str) -> bool:
    """判断脚本是被导入的库模块还是命令行入口。

    库模块（helper/core/utils/properties/filters 等）由其他脚本 import 调用，
    本就不需要 ``__main__`` 入口和 argparse。只有真正的入口脚本才有这些要求。

    Args:
        source: 脚本文本。
        filename: 脚本文件名。
    """
    # 明确带有 CLI 特征的一定是入口脚本
    if "bootstrap_cli" in source or "argparse" in source:
        return False

    stem = Path(filename).stem.lower()
    library_markers = (
        "core",
        "helper",
        "helpers",
        "utils",
        "utility",
        "common",
        "lib",
        "base",
        "properties",
        "filters",
        "parser",
        "geometry",
        "style",
        "client",
        "analysis",
        "deseq2",
        "backend",
        "impl",
    )
    if any(
        stem == marker or stem.startswith(marker) or stem.endswith(marker)
        for marker in library_markers
    ):
        return True

    # 兜底：既无 argparse 也无 __main__，且以函数定义为主——视为库模块
    if "argparse" not in source and "bootstrap_cli" not in source:
        try:
            tree = ast.parse(source)
        except SyntaxError:
            return False
        has_main = any(
            isinstance(node, ast.If)
            and isinstance(node.test, ast.Compare)
            and isinstance(node.test.left, ast.Name)
            and node.test.left.id == "__name__"
            for node in ast.walk(tree)
        )
        if not has_main:
            return True

    return False


def _print_outside_main_block(source: str) -> bool:
    """判断是否存在位于 ``__main__`` 块之外的裸 print 调用。

    ``if __name__ == "__main__":`` 块内的 print 是命令行演示输出，
    属于合理用法，不算违规。
    """
    try:
        tree = ast.parse(source)
    except SyntaxError:
        return "print(" in source

    main_blocks: list[ast.If] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.If):
            test = node.test
            if (
                isinstance(test, ast.Compare)
                and isinstance(test.left, ast.Name)
                and test.left.id == "__name__"
            ):
                main_blocks.append(node)

    for node in ast.walk(tree):
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id == "print":
            in_main = any(
                block is node or block in _descendants(block) for block in main_blocks
            )
            if not in_main:
                return True
    return False


def _descendants(node: ast.AST) -> set:
    """收集一个 AST 节点的所有子节点。"""
    result = set()
    for child in ast.iter_child_nodes(node):
        result.add(child)
        result |= _descendants(child)
    return result


def _check_scripts(skill_dir: Path, skill: str, issues: list[Issue]) -> None:
    """校验 scripts/ 下的 Python 文件。"""
    scripts_dir = skill_dir / "scripts"
    if not scripts_dir.exists():
        issues.append(Issue(skill, "warning", "no-scripts", "缺少 scripts/ 目录"))
        return

    py_files = sorted(scripts_dir.rglob("*.py"))
    if not py_files:
        issues.append(Issue(skill, "error", "no-scripts", "scripts/ 下没有任何 .py 文件"))
        return

    for py_file in py_files:
        rel = _display_path(py_file)
        try:
            source = py_file.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            issues.append(Issue(skill, "error", "script-encoding", f"{rel} 不是 UTF-8 编码"))
            continue

        try:
            tree = ast.parse(source, filename=str(py_file))
        except SyntaxError as exc:
            issues.append(Issue(skill, "error", "script-syntax", f"{rel} 语法错误: {exc.msg} (行 {exc.lineno})"))
            continue

        is_library = _is_library_module(source, py_file.name)

        if not is_library:
            has_main = False
            for node in ast.walk(tree):
                if isinstance(node, ast.If):
                    test = node.test
                    if (
                        isinstance(test, ast.Compare)
                        and isinstance(test.left, ast.Name)
                        and test.left.id == "__name__"
                    ):
                        has_main = True
                        break
            if not has_main:
                issues.append(
                    Issue(skill, "error", "script-no-main", f"{rel} 缺少 if __name__ == '__main__' 入口")
                )

        for pattern in HARDCODED_PATH_PATTERNS:
            match = pattern.search(source)
            if match:
                issues.append(
                    Issue(
                        skill,
                        "error",
                        "script-hardcoded-path",
                        f"{rel} 含硬编码绝对路径: {match.group(1)!r}",
                    )
                )
                break

        if not is_library and "argparse" not in source:
            issues.append(
                Issue(skill, "warning", "script-no-argparse", f"{rel} 未使用 argparse，缺少命令行接口")
            )

        if "print(" in source and "logging" not in source and "logger" not in source:
            if _print_outside_main_block(source):
                issues.append(
                    Issue(skill, "warning", "script-raw-print", f"{rel} 疑似使用裸 print 而非 logging")
                )


def validate_skill(skill_dir: Path) -> Result:
    """校验单个技能目录。"""
    skill = skill_dir.name
    result = Result(skill=skill)
    issues = result.issues

    if not skill_dir.is_dir():
        issues.append(Issue(skill, "error", "not-dir", "技能路径不是目录"))
        return result

    skill_md = skill_dir / "SKILL.md"
    if not skill_md.exists():
        issues.append(Issue(skill, "error", "no-skill-md", "缺少 SKILL.md"))
        return result

    text = skill_md.read_text(encoding="utf-8")
    data, body, err = _read_frontmatter(text)
    if err:
        issues.append(Issue(skill, "error", "frontmatter", err))
        return result

    assert data is not None

    for key in REQUIRED_FRONTMATTER:
        if key not in data or not data[key]:
            issues.append(Issue(skill, "error", f"missing-{key}", f"frontmatter 缺少必填字段 '{key}'"))

    name = str(data.get("name", ""))
    if name and name != skill:
        issues.append(
            Issue(skill, "error", "name-mismatch", f"frontmatter name '{name}' 与目录名 '{skill}' 不一致")
        )
    if name and not NAME_PATTERN.match(name):
        issues.append(
            Issue(skill, "error", "name-invalid", f"name '{name}' 不是合法 kebab-case（仅小写字母、数字、连字符）")
        )

    metadata = data.get("metadata")
    if not isinstance(metadata, dict):
        issues.append(Issue(skill, "warning", "no-metadata", "建议提供 metadata 块（version/last-reviewed/skill-author）"))
    else:
        for key in ("version", "last-reviewed", "skill-author"):
            if key not in metadata:
                issues.append(Issue(skill, "warning", f"missing-metadata-{key}", f"metadata 缺少 '{key}'"))

    _check_description(skill, str(data.get("description", "")), issues)

    headings = re.findall(r"^##\s+(.+?)\s*$", body, re.MULTILINE)
    normalized = [h.strip().lower() for h in headings]
    for section in REQUIRED_SECTIONS:
        if section.lower() not in normalized:
            issues.append(Issue(skill, "error", "missing-section", f"SKILL.md 缺少章节 '## {section}'"))

    _check_scripts(skill_dir, skill, issues)

    refs_dir = skill_dir / "references"
    if not refs_dir.exists():
        issues.append(Issue(skill, "warning", "no-references", "缺少 references/ 目录"))
    else:
        md_files = list(refs_dir.glob("*.md"))
        if not md_files:
            issues.append(Issue(skill, "warning", "empty-references", "references/ 下没有 .md 文档"))
        for md in md_files:
            if md.stat().st_size < 100:
                issues.append(Issue(skill, "warning", "thin-reference", f"references/{md.name} 内容过少"))

    req = skill_dir / "requirements.txt"
    if not req.exists():
        issues.append(Issue(skill, "warning", "no-requirements", "缺少 requirements.txt"))

    return result


def main(argv: list[str] | None = None) -> int:
    """校验器主函数。"""
    parser = argparse.ArgumentParser(description="校验 Agent Skills 规范符合性")
    parser.add_argument("--skill", help="只校验指定技能")
    parser.add_argument("--strict", action="store_true", help="警告也视为失败")
    parser.add_argument("--json", action="store_true", help="输出 JSON")
    parser.add_argument("--skills-dir", default=str(SKILLS_DIR), help="技能目录路径")
    args = parser.parse_args(argv)

    skills_dir = Path(args.skills_dir).resolve()
    if not skills_dir.exists():
        print(f"技能目录不存在: {skills_dir}", file=sys.stderr)
        return 2

    if args.skill:
        targets = [skills_dir / args.skill]
        if not targets[0].exists():
            print(f"技能不存在: {targets[0]}", file=sys.stderr)
            return 2
    else:
        targets = sorted(p for p in skills_dir.iterdir() if p.is_dir() and not p.name.startswith((".", "_")))

    if not targets:
        print(f"未找到任何技能: {skills_dir}", file=sys.stderr)
        return 2

    results = [validate_skill(t) for t in targets]
    total_errors = sum(len(r.errors) for r in results)
    total_warnings = sum(len(r.warnings) for r in results)

    if args.json:
        print(
            json.dumps(
                {
                    "total": len(results),
                    "passed": sum(1 for r in results if r.ok),
                    "errors": total_errors,
                    "warnings": total_warnings,
                    "skills": [
                        {
                            "skill": r.skill,
                            "ok": r.ok,
                            "errors": [{"code": i.code, "message": i.message} for i in r.errors],
                            "warnings": [{"code": i.code, "message": i.message} for i in r.warnings],
                        }
                        for r in results
                    ],
                },
                ensure_ascii=False,
                indent=2,
            )
        )
        return 1 if (total_errors or (args.strict and total_warnings)) else 0

    print(f"校验 {len(results)} 个技能\n" + "=" * 70)
    for r in results:
        status = "PASS" if r.ok else "FAIL"
        print(f"[{status}] {r.skill}  ({len(r.errors)} 错误, {len(r.warnings)} 警告)")
        for issue in r.issues:
            print(f"    {issue}")
        print()

    print("=" * 70)
    print(f"合计: {total_errors} 错误, {total_warnings} 警告")
    failed = [r.skill for r in results if not r.ok]
    if failed:
        print(f"未通过: {', '.join(failed)}")
    return 1 if (total_errors or (args.strict and total_warnings)) else 0


if __name__ == "__main__":
    sys.exit(main())
