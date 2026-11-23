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

sys.path.insert(0, str(Path(__file__).resolve().parent))
from skill_frontmatter import FrontmatterError, split_frontmatter  # noqa: E402

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

# 触发语句：第三人称 + when 语义。
#
# 两种等价写法都接受：
#   * "This skill should be used when ..."（仓库早期的写法）
#   * "Use when ..."（Agent Skills 官方 optimizing-descriptions 指南推荐的
#     祈使句式——描述被当作给代理的指令，且官方明确要求把关键用例放在最前面）
TRIGGER_PATTERNS = [
    r"This skill should be used when",
    r"\bshould be used when\b",
    r"\bUse when\b",
    r"\bUse this skill when\b",
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

# 技能列表的共享字符预算。Claude Code 在启动时把所有技能的
# name + description 装进上下文；超出预算后会**从最少被调用的技能开始
# 丢弃 description**，那些技能只剩一个裸名字，代理便无从匹配用户请求，
# 也就永远不会触发它们。这是最隐蔽的一类失败——没有任何报错。
#
# 官方默认值是「上下文窗口的 1%，fallback 到 8,000 字符」。
# 这里用一个保守上限做门禁：留出余量给用户自己安装的其他技能。
DESC_BUDGET = 8000
# 安全余量：目标上限，比硬预算低，留给下游宿主的额外开销与用户
# 自己安装的其他技能。设为 7,000（硬预算的 87.5%）——全量 17 个技能
# 装齐时占 84.7%，仍有 1,000 字符余量；而按 skill-guides 推荐只装
# 2–3 个时仅占约 15%。
#
# 不再往下压：描述的触发准确率比省下的字符更值钱。
DESC_BUDGET_SOFT = 7000


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

    刻意**不**用 ``text.split("---", 2)``：Markdown 表格的分隔行
    （``|---|---|``）里含有连续的 ``-``，会被误判为 frontmatter 的结束
    标记，把大半个正文当成 YAML 解析。统一走
    ``skill_frontmatter.split_frontmatter``，它只认独占一行的 ``---``。

    Returns:
        ``(frontmatter_dict, body, error_message)``；解析失败时 dict 为 None。
    """
    try:
        raw, body = split_frontmatter(text)
    except FrontmatterError as exc:
        return None, text, str(exc)

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


def _check_edam(skill: str, metadata: dict) -> list[Issue]:
    """校验 metadata 里的 EDAM 本体标注。

    EDAM（https://edamontology.org）是生信工具元数据的标准本体。
    标注本身不强制，但**一旦写了就必须合法**——错误的 ID 比没有 ID
    更糟：它会让工具在基于本体的检索里"存在但查不到"，或者被挂到
    错误的语义分支下。

    **键的形态**：Agent Skills 规范规定 ``metadata`` 是
    "a map from string keys to string values"（字符串键 → 字符串值），
    因此不能嵌套映射或列表。所以 EDAM 标注用扁平的
    ``edam-<axis>`` 键，值是逗号分隔的 ID 字符串。

    允许的轴与 ID 前缀：
      topic / operation / input-data / output-data → 对应 EDAM 四个子本体
      input-format / output-format → format_ 前缀
    """
    issues: list[Issue] = []

    axis_prefix = {
        "edam-topic": "topic_",
        "edam-operation": "operation_",
        "edam-input-data": "data_",
        "edam-output-data": "data_",
        "edam-input-format": "format_",
        "edam-output-format": "format_",
    }

    found_axes: list[str] = []
    for key, value in metadata.items():
        if not key.startswith("edam-"):
            continue

        if key not in axis_prefix:
            issues.append(
                Issue(
                    skill,
                    "error",
                    "edam-unknown-axis",
                    f"metadata.{key} 不是已知的 EDAM 轴，允许的键为 {sorted(axis_prefix)}",
                )
            )
            continue

        found_axes.append(key)

        if not isinstance(value, str):
            issues.append(
                Issue(
                    skill,
                    "error",
                    "edam-value-not-string",
                    f"metadata.{key} 的值必须是字符串"
                    f"（规范要求 metadata 是字符串到字符串的映射）",
                )
            )
            continue

        ids = [item.strip() for item in value.split(",") if item.strip()]
        if not ids:
            issues.append(
                Issue(skill, "warning", "edam-axis-empty", f"metadata.{key} 没有列出任何 ID")
            )
            continue

        prefix = axis_prefix[key]
        for edam_id in ids:
            if not edam_id.startswith(prefix):
                issues.append(
                    Issue(
                        skill,
                        "error",
                        "edam-id-prefix",
                        f"metadata.{key} 的 '{edam_id}' 前缀应为 '{prefix}'",
                    )
                )

    if found_axes and not (
        "edam-topic" in found_axes or "edam-operation" in found_axes
    ):
        issues.append(
            Issue(
                skill,
                "warning",
                "edam-missing-function",
                "EDAM 标注至少应包含 edam-topic 或 edam-operation，否则无法表达工具的功能",
            )
        )

    return issues


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
        issues.extend(_check_edam(skill, metadata))

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


def _collect_listing_cost(results: list[Result]) -> tuple[int, list[tuple[str, int]]]:
    """统计所有技能的 name + description 总字符数。

    这是技能列表在每次会话启动时的常驻成本。超出宿主的共享预算后，
    描述会被丢弃，技能便不再能被触发——因此必须作为门禁。
    """
    per_skill: list[tuple[str, int]] = []
    total = 0
    for r in results:
        fm = _read_frontmatter((SKILLS_DIR / r.skill / "SKILL.md").read_text(encoding="utf-8"))[0]
        if fm is None:
            continue
        try:
            data = yaml.safe_load(fm) or {}
        except Exception:
            continue
        cost = len(str(data.get("name", ""))) + len(str(data.get("description", "")))
        per_skill.append((r.skill, cost))
        total += cost
    return total, per_skill


def _check_listing_budget(results: list[Result]) -> list[Issue]:
    """检查技能列表的常驻成本是否在预算内。"""
    issues: list[Issue] = []
    total, per_skill = _collect_listing_cost(results)
    if total == 0:
        return issues

    if total > DESC_BUDGET:
        issues.append(
            Issue(
                "ALL",
                "error",
                "listing-budget-exceeded",
                f"全部技能的 name+description 合计 {total} 字符，"
                f"超出宿主共享预算 {DESC_BUDGET} 字符。超出的技能会被丢弃描述、"
                f"只剩裸名字而无法被触发。请压缩描述或改用 skillOverrides 降级低频技能。",
            )
        )
    elif total > DESC_BUDGET_SOFT:
        issues.append(
            Issue(
                "ALL",
                "warning",
                "listing-budget-tight",
                f"全部技能的 name+description 合计 {total} 字符，"
                f"已超过建议上限 {DESC_BUDGET_SOFT} 字符（硬预算 {DESC_BUDGET}）。"
                f"请预留余量给用户自己安装的其他技能。",
            )
        )

    # 单条描述接近硬上限时提示——超限会让技能在部分宿主上整个消失
    for skill, cost in per_skill:
        if cost > MAX_DESC:
            issues.append(
                Issue(
                    skill,
                    "error",
                    "listing-entry-too-long",
                    f"name+description 合计 {cost} 字符，超过单条上限 {MAX_DESC}。"
                    f"部分宿主（如 copilot-cli）会直接不加载超限技能。",
                )
            )
    return issues


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

    # 技能列表的常驻成本是跨技能的全局约束，单独校验一次。
    # 只在检查全量技能时生效（--skill 单项检查无法评估全局预算）。
    if not args.skill:
        budget_issues = _check_listing_budget(results)
        if budget_issues:
            results.append(Result(skill="__listing__", issues=budget_issues))

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
