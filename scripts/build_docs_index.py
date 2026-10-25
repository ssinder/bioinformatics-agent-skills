#!/usr/bin/env python3
"""生成技能索引页 docs/skills/index.md，并把同一份表格内嵌到首页。

从每个技能的 SKILL.md frontmatter 读取 name / description，
生成一个自动维护的技能总览表。

为什么用脚本生成
----------------
手工维护的索引一定会和代码漂移——加了新技能却忘了更新文档表，
或者改了 description 却没同步索引。这里从 SKILL.md 直接读取，
保证索引与实际技能一致。首页与索引页共用同一份数据，两处不会不一致。

用法：
    python scripts/build_docs_index.py            # 生成
    python scripts/build_docs_index.py --check    # 只校验是否最新（CI 用）
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SKILLS_DIR = REPO_ROOT / "skills"
OUTPUT = REPO_ROOT / "docs" / "skills" / "index.md"
HOME = REPO_ROOT / "docs" / "index.md"

# 技能源码在仓库根的 skills/ 下，不在 docs/ 下，
# 因此站点上的链接必须指向 GitHub 上的实际文件，否则会 404。
REPO_BLOB = "https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills"

_MARKER_START = "<!-- SKILLS_TABLE_START -->"
_MARKER_END = "<!-- SKILLS_TABLE_END -->"

_HEADER = """---
layout: default
title: 技能总览
---

# 技能总览 / Skills Overview

> 本页由 `scripts/build_docs_index.py` 自动生成，请勿手工编辑。

下表列出本仓库提供的全部 Agent Skills。

| 技能 / Skill | 说明 / Description |
|---|---|
"""

_FOOTER = """
## 安装方式

将 `skills/` 目录复制到宿主的技能目录即可被发现：

```bash
git clone https://github.com/ssinder/bioinformatics-agent-skills.git
cp -r bioinformatics-agent-skills/skills/* ~/.claude/skills/
```

详见 [快速上手](/getting-started.html)。
"""


def _parse_frontmatter(path: Path) -> dict[str, str]:
    """解析 SKILL.md 的 YAML frontmatter（只取顶层 key: value，够用即可）。"""
    text = path.read_text(encoding="utf-8")
    if not text.startswith("---"):
        return {}
    end = text.find("\n---", 3)
    if end == -1:
        return {}
    block = text[3:end]

    data: dict[str, str] = {}
    current_key = None
    current_val: list[str] = []

    for raw in block.splitlines():
        if not raw.strip():
            continue
        # 顶层 key: value
        m = re.match(r"^([A-Za-z_][A-Za-z0-9_\-]*):\s*(.*)$", raw)
        if m and not raw.startswith((" ", "\t", "-")):
            if current_key:
                data[current_key] = " ".join(current_val).strip().strip('"').strip("'")
            current_key = m.group(1)
            current_val = [m.group(2)] if m.group(2).strip() else []
            continue
        if current_key:
            current_val.append(raw.strip().lstrip("-").strip())

    if current_key:
        data[current_key] = " ".join(current_val).strip().strip('"').strip("'")

    return data


def _clean_description(raw: str) -> str:
    """整理 description 为适合表格的一行文本。"""
    text = raw.replace("\n", " ").strip()
    # 去掉 YAML 折叠标记残留（如 ">-"、">"）
    text = re.sub(r"^>[-+]?\s*", "", text)
    text = text.strip('"').strip("'").strip()
    text = re.sub(r"\s+", " ", text)
    # 表格里不能出现裸的竖线
    text = text.replace("|", "\\|")
    return text


def _skill_rows() -> list[str]:
    """构建技能表格行。"""
    rows: list[str] = []
    if not SKILLS_DIR.is_dir():
        return rows

    for skill_dir in sorted(d for d in SKILLS_DIR.iterdir() if d.is_dir()):
        skill_md = skill_dir / "SKILL.md"
        if not skill_md.is_file():
            continue
        meta = _parse_frontmatter(skill_md)
        name = meta.get("name", skill_dir.name)
        desc = _clean_description(meta.get("description", ""))
        link = f"{REPO_BLOB}/{skill_dir.name}/SKILL.md"
        rows.append(f"| [`{name}`]({link}) | {desc} |")
    return rows


def build_index() -> str:
    """生成 docs/skills/index.md 的完整内容。"""
    return _HEADER + "\n".join(_skill_rows()) + "\n" + _FOOTER


def _apply_to_home(home_text: str, table: str) -> str:
    """把表格写入首页的标记区之间。"""
    if _MARKER_START not in home_text or _MARKER_END not in home_text:
        return home_text
    pre = home_text.split(_MARKER_START)[0]
    post = home_text.split(_MARKER_END)[1]
    return f"{pre}{_MARKER_START}\n{table}\n{_MARKER_END}{post}"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="只校验，不写入")
    args = parser.parse_args()

    content = build_index()
    table = "\n".join(_skill_rows())

    home_current = HOME.read_text(encoding="utf-8") if HOME.is_file() else ""
    home_updated = _apply_to_home(home_current, table)

    if args.check:
        problems: list[str] = []
        if not OUTPUT.is_file():
            problems.append(f"{OUTPUT.relative_to(REPO_ROOT)} 不存在")
        elif OUTPUT.read_text(encoding="utf-8") != content:
            problems.append(f"{OUTPUT.relative_to(REPO_ROOT)} 与 SKILL.md 不同步")
        if home_updated != home_current:
            problems.append("docs/index.md 中的技能表格与 SKILL.md 不同步")

        if problems:
            for p in problems:
                print(f"::error::{p}，请运行 python scripts/build_docs_index.py")
            return 1
        print("技能索引页与 SKILL.md 保持同步")
        return 0

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(content, encoding="utf-8")
    if home_updated != home_current:
        HOME.write_text(home_updated, encoding="utf-8")
        print("已更新 docs/index.md 中的技能表格")

    print(f"已生成 {OUTPUT.relative_to(REPO_ROOT)}（{len(_skill_rows())} 个技能）")
    return 0


if __name__ == "__main__":
    sys.exit(main())
