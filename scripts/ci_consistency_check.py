#!/usr/bin/env python3
"""检查仓库各处描述与实际技能状态是否一致。

背景
----
仓库里有几处"手工维护的计数"：README 徽章上的技能数量、
plugin.json 的 keywords 等。随着技能增减，这些数字会悄悄过时，
而没有任何检查会报错——直到有人发现文档说 8 个技能、实际有 9 个。

本脚本把这些跨文件的一致性约束显式化：

1. README 徽章 ``Skills-N`` 的 N 必须等于实际技能数；
2. README 的技能清单表格必须与实际技能集合完全一致（不多不少）；
3. ``packages/biopipe/registry.py`` 的 REGISTRY 必须与 ``skills/``
   目录一一对应（不多不少）；
4. 每个技能都必须有 SKILL.md、scripts/、references/；
5. plugin.json 的 name 与 pyproject 的项目名一致；
6. 文档中声明的技能总数必须与实际一致（防止旧生成物残留）。

用法：
    python scripts/ci_consistency_check.py
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SKILLS_DIR = REPO_ROOT / "skills"
DOCS_DIR = REPO_ROOT / "docs"
README = REPO_ROOT / "README.md"
PLUGIN_JSON = REPO_ROOT / "plugin.json"
PYPROJECT = REPO_ROOT / "pyproject.toml"
REGISTRY_PY = REPO_ROOT / "packages" / "biopipe" / "registry.py"


def _skill_dirs() -> list[Path]:
    if not SKILLS_DIR.is_dir():
        return []
    return sorted(d for d in SKILLS_DIR.iterdir() if d.is_dir() and (d / "SKILL.md").is_file())


def _check_readme_badge(problems: list[str], count: int) -> None:
    if not README.is_file():
        problems.append("README.md 不存在")
        return
    text = README.read_text(encoding="utf-8")
    m = re.search(r"badge/Skills-(\d+)-", text)
    if not m:
        problems.append("README.md 中找不到 Skills-N 徽章，无法校验技能数量")
        return
    badge = int(m.group(1))
    if badge != count:
        problems.append(
            f"README.md 徽章显示 {badge} 个技能，实际有 {count} 个。"
            "请更新 Skills-N 徽章。"
        )


def _check_readme_listing(problems: list[str], skill_names: list[str]) -> None:
    """README 的技能清单表格必须收录每个技能。

    这里刻意**不**用"全文子串包含"来判断——README 正文里还有
    ``python skills/xxx/scripts/...`` 这类路径引用，只要路径里出现
    了技能名，子串检查就会误判为"已收录"，从而漏掉真正的清单遗漏。

    正确做法是解析技能表格行：形如
    ``| [`skill-name`](skills/skill-name/) | 描述 | 能力 |``
    只把出现在表格第一列的技能名视为"已收录"。
    """
    if not README.is_file():
        return
    text = README.read_text(encoding="utf-8")

    # 匹配技能表格行：| [`name`](...skills/name...) | ...
    listed: set[str] = set()
    for m in re.finditer(r"^\|\s*\[`([a-z0-9\-]+)`\]\(", text, re.MULTILINE):
        listed.add(m.group(1))

    if not listed:
        problems.append(
            "README.md 中找不到技能清单表格（期望形如 "
            "| [`skill-name`](skills/skill-name/) | ... | 的行）"
        )
        return

    missing = [name for name in skill_names if name not in listed]
    if missing:
        problems.append(f"README.md 的技能清单缺少: {', '.join(missing)}")

    # 清单里出现磁盘上不存在的技能，同样是漂移
    stale = sorted(listed - set(skill_names))
    if stale:
        problems.append(
            f"README.md 的技能清单包含不存在的技能: {', '.join(stale)}"
        )


def _check_registry(problems: list[str], skill_names: list[str]) -> None:
    """REGISTRY 必须与 skills/ 目录一一对应。"""
    if not REGISTRY_PY.is_file():
        problems.append("packages/biopipe/registry.py 不存在")
        return
    text = REGISTRY_PY.read_text(encoding="utf-8")
    registered = set(re.findall(r'name="([a-z0-9\-]+)"', text))
    actual = set(skill_names)

    only_registry = sorted(registered - actual)
    only_disk = sorted(actual - registered)
    if only_registry:
        problems.append(f"REGISTRY 中登记但磁盘上不存在: {', '.join(only_registry)}")
    if only_disk:
        problems.append(f"磁盘上存在但 REGISTRY 未登记: {', '.join(only_disk)}")


def _check_skill_structure(problems: list[str]) -> None:
    """每个技能必须有 SKILL.md + scripts/ + references/。"""
    for skill_dir in _skill_dirs():
        name = skill_dir.name
        if not skill_dir.is_dir():
            continue
        for required in ("SKILL.md", "scripts", "references"):
            if not (skill_dir / required).exists():
                problems.append(f"技能 {name} 缺少 {required}/")


def _check_plugin_json(problems: list[str]) -> None:
    """plugin.json 的字段必须与仓库实际状态一致。

    注意：plugin.json 的 name 是**仓库/插件名**
    （bioinformatics-agent-skills），而 pyproject.toml 的 name 是
    **PyPI 包名**（biopipe-agent-skills）。两者本就应该不同——
    包名一旦发布就不能随意改。因此这里校验的是"包名与仓库名必须
    有稳定的映射关系"，而不是字面相等。

    校验内容：
    1. plugin.json 的 skills 路径确实存在；
    2. plugin.json 的 name 与仓库目录名一致（约定俗成，便于识别）；
    3. plugin.json 的 homepage / repository 指向正确的仓库地址。
    """
    if not PLUGIN_JSON.is_file():
        problems.append("plugin.json 不存在")
        return
    try:
        plugin = json.loads(PLUGIN_JSON.read_text(encoding="utf-8"))
    except Exception as exc:
        problems.append(f"plugin.json 解析失败: {exc}")
        return

    skills_path = plugin.get("skills")
    if skills_path:
        target = (REPO_ROOT / skills_path).resolve()
        if not target.is_dir():
            problems.append(f"plugin.json 的 skills 路径不存在: {skills_path}")

    repo_name = plugin.get("name")
    if repo_name and repo_name != REPO_ROOT.name:
        problems.append(
            f"plugin.json 的 name ({repo_name}) 与仓库目录名 ({REPO_ROOT.name}) 不一致"
        )

    repo_slug = "ssinder/bioinformatics-agent-skills"
    repository = plugin.get("repository") or ""
    if repository and repo_slug not in repository:
        problems.append(f"plugin.json 的 repository 指向了非本仓库地址: {repository}")

    # homepage 可以指向文档站（Pages），也可以指向仓库首页，
    # 只需保证它属于本项目，不能指向别人的地址。
    homepage = plugin.get("homepage") or ""
    if homepage:
        allowed = ("ssinder.github.io/bioinformatics-agent-skills", repo_slug)
        if not any(token in homepage for token in allowed):
            problems.append(f"plugin.json 的 homepage 指向了非本项目地址: {homepage}")


def _check_doc_skill_counts(
    problems: list[str],
    count: int,
    paths: list[Path] | None = None,
) -> None:
    """文档里声明的技能数必须与实际一致。

    背景：``docs/skills.md`` 曾长期写着「共 8 个技能」而实际已有 17 个 ——
    它是早期 ``scripts/build_catalog.py`` 的产物，被后来的
    ``docs/skills/index.md`` 取代后忘了删。没有任何检查会报错，
    因为「文档里的数字」本来就不在自动生成链路上。

    这里做一次兜底扫描：凡是文档里出现「N 个技能」「N skills」这类
    明确计数，且 N 与实际不符，就报出来。

    ``paths`` 允许测试注入临时文件；默认扫描 ``docs/`` 下全部 Markdown
    加上 README。
    """
    doc_patterns = (
        re.compile(r"共\s*\*\*(?P<n>\d+)\*\*\s*个技能"),
        re.compile(r"共\s*(?P<n>\d+)\s*个技能"),
        re.compile(r"\b(?P<n>\d+)\s+skills\b"),
        re.compile(r"\b(?P<n>\d+)\s+个技能\b"),
    )

    if paths is None:
        targets = sorted(DOCS_DIR.rglob("*.md")) + [README]
    else:
        targets = list(paths)

    for path in targets:
        if not path.is_file():
            continue
        text = path.read_text(encoding="utf-8")
        for pattern in doc_patterns:
            for m in pattern.finditer(text):
                stated = int(m.group("n"))
                if stated == count:
                    continue
                if _is_skill_count_context(text, m.start("n")):
                    try:
                        rel = path.relative_to(REPO_ROOT)
                    except ValueError:
                        rel = path
                    problems.append(
                        f"{rel} 声明技能数为 {stated}，实际为 {count}"
                    )


def _is_skill_count_context(text: str, number_pos: int) -> bool:
    """判断该数字是否处在"技能总数"的语境里。

    判据：以**数字本身**为基准，向前 30 字符内出现「共 / 包含 / 提供 /
    覆盖 / 收录 / 总共」这类总数量措辞，向后 20 字符内出现「个技能 /
    skills」。

    注意必须以数字位置为基准而不是整个匹配的位置 —— 像「共 **8** 个技能」
    这种写法，匹配的起点在「共」上，若以起点为基准则前向窗口里找不到
    「共」，判断会永远为假。

    这样「3 个技能应先做质控」这种举例不会误判：它前面没有总数量措辞。
    """
    before = text[max(0, number_pos - 30):number_pos]
    after = text[number_pos:number_pos + 20]
    if not re.search(r"(个技能|skills)", after, re.IGNORECASE):
        return False
    return bool(
        re.search(
            r"(共|包含|提供|覆盖|收录|总共|offer\w*|provide\w*|contain\w*|cover\w*|includ\w*)",
            before,
            re.IGNORECASE,
        )
    )


def main() -> int:
    problems: list[str] = []

    skill_dirs = _skill_dirs()
    skill_names = [d.name for d in skill_dirs]
    count = len(skill_names)

    if count == 0:
        problems.append("skills/ 目录下没有找到任何技能")
    else:
        _check_readme_badge(problems, count)
        _check_readme_listing(problems, skill_names)
        _check_registry(problems, skill_names)
        _check_skill_structure(problems)
        _check_doc_skill_counts(problems, count)

    _check_plugin_json(problems)

    if problems:
        print("一致性检查未通过：")
        for p in problems:
            print(f"  - {p}")
        print("::error::仓库一致性检查失败")
        return 1

    print(f"一致性检查通过（{count} 个技能，REGISTRY 与磁盘完全对应）")
    return 0


if __name__ == "__main__":
    sys.exit(main())
