#!/usr/bin/env python3
"""批量改写技能 frontmatter 里的 description 字段。

为什么需要这个脚本
------------------
SKILL.md 的 frontmatter 里，``description`` 可能是两种写法：

* 单行：``description: "文本"``
* 折叠块：``description: >-`` 后跟若干缩进行

折叠块是本仓库多数技能采用的写法（长描述用多行更好读），但也
意味着**简单的单行正则替换会静默失配**——不报错，却什么都没改。
这个脚本用 YAML 解析器读取 frontmatter，按块边界切分并重写，避免
这类"看起来跑了、其实没生效"的失败。

用法：
    python scripts/set_skill_description.py --from-file descriptions.json
    python scripts/set_skill_description.py --from-file descriptions.json --check
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from skill_frontmatter import split_frontmatter  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[1]
SKILLS_DIR = REPO_ROOT / "skills"

try:
    import yaml
except ImportError:  # pragma: no cover
    yaml = None


def _split_frontmatter(text: str) -> tuple[str, str]:
    """把文件拆成 (frontmatter 原文, 正文)。

    走共用的健壮实现——不能用 ``split("---", 2)``，Markdown 表格的
    分隔行里含有连续的 ``-``，会被误判为 frontmatter 结束标记。
    """
    return split_frontmatter(text)


def _replace_description_block(fm: str, new_desc: str) -> str:
    """在 frontmatter 文本中替换 description 字段（兼容单行与折叠块）。"""
    lines = fm.split("\n")
    out: list[str] = []
    i = 0
    replaced = False
    while i < len(lines):
        line = lines[i]
        m = re.match(r"^description\s*:", line)
        if m and not replaced:
            out.append(f'description: "{new_desc}"')
            i += 1
            # 跳过原折叠块的续行（缩进的行）
            while i < len(lines) and lines[i].strip() and lines[i][:1] in (" ", "\t"):
                i += 1
            replaced = True
            continue
        out.append(line)
        i += 1
    if not replaced:
        raise ValueError("frontmatter 中找不到 description 字段")
    return "\n".join(out)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="批量改写技能 description")
    parser.add_argument("--from-file", required=True, help="JSON 文件：{技能名: 描述}")
    parser.add_argument("--check", action="store_true", help="只检查不写入")
    args = parser.parse_args(argv)

    if yaml is None:
        print("缺少 PyYAML，无法解析 frontmatter", file=sys.stderr)
        return 2

    data = json.loads(Path(args.from_file).read_text(encoding="utf-8"))

    problems: list[str] = []
    changed = 0
    total = 0

    for name, desc in data.items():
        skill_md = SKILLS_DIR / name / "SKILL.md"
        if not skill_md.is_file():
            problems.append(f"技能 {name} 不存在")
            continue

        text = skill_md.read_text(encoding="utf-8")
        fm, body = _split_frontmatter(text)

        parsed = yaml.safe_load(fm)
        current = str((parsed or {}).get("description", ""))
        total += len(str((parsed or {}).get("name", name))) + len(current)

        if current == desc:
            continue

        if args.check:
            problems.append(f"{name}: description 与目标不一致（当前 {len(current)}c → 目标 {len(desc)}c）")
            continue

        esc = desc.replace("\\", "\\\\").replace('"', '\\"')
        new_fm = _replace_description_block(fm, esc)
        skill_md.write_text("---" + new_fm + "---" + body, encoding="utf-8")
        changed += 1

    if args.check:
        if problems:
            print("description 检查未通过：")
            for p in problems:
                print(f"  - {p}")
            return 1
        print(f"description 检查通过：{len(data)} 个技能均与目标一致")
        return 0

    # 改写后重新统计全部技能的预算占用
    grand_total = 0
    for skill_dir in sorted(SKILLS_DIR.iterdir()):
        if not skill_dir.is_dir() or skill_dir.name.startswith((".", "_")):
            continue
        skill_md = skill_dir / "SKILL.md"
        if not skill_md.is_file():
            continue
        fm, _ = _split_frontmatter(skill_md.read_text(encoding="utf-8"))
        parsed = yaml.safe_load(fm) or {}
        grand_total += len(str(parsed.get("name", ""))) + len(str(parsed.get("description", "")))

    print(f"改写 {changed} 个技能；全库 name+description 合计 {grand_total} 字符")
    if problems:
        for p in problems:
            print(f"  ! {p}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
