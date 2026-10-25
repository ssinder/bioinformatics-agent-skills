#!/usr/bin/env python3
"""生成技能索引页 docs/skills/index.md。

从每个技能的 SKILL.md frontmatter 读取 name / description，
生成一个自动维护的技能总览表。

为什么用脚本生成
----------------
手工维护的索引一定会和代码漂移——加了新技能却忘了更新文档表，
或者改了 description 却没同步索引。这里从 SKILL.md 直接读取，
保证索引与实际技能一致。

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

_HEADER = """# 技能总览 / Skills Overview

> 本页由 `scripts/build_docs_index.py` 自动生成，请勿手工编辑。

下表列出本仓库提供的全部 Agent Skills。点击技能名可查看该技能的
详细文档与源码。

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

详见 [快速上手](../getting-started.md)。
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


def build_index() -> str:
    rows: list[str] = []
    skill_dirs = sorted(d for d in SKILLS_DIR.iterdir() if d.is_dir()) if SKILLS_DIR.is_dir() else []

    for skill_dir in skill_dirs:
        skill_md = skill_dir / "SKILL.md"
        if not skill_md.is_file():
            continue
        meta = _parse_frontmatter(skill_md)
        name = meta.get("name", skill_dir.name)
        desc = _clean_description(meta.get("description", ""))
        rows.append(f"| [`{name}`](./{skill_dir.name}/SKILL.md) | {desc} |")

    return _HEADER + "\n".join(rows) + "\n" + _FOOTER


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="只校验，不写入")
    args = parser.parse_args()

    content = build_index()

    if args.check:
        if not OUTPUT.is_file():
            print(f"::error::{OUTPUT.relative_to(REPO_ROOT)} 不存在，请运行 scripts/build_docs_index.py")
            return 1
        current = OUTPUT.read_text(encoding="utf-8")
        if current != content:
            print(
                f"::error::{OUTPUT.relative_to(REPO_ROOT)} 与 SKILL.md 不同步，"
                "请运行 python scripts/build_docs_index.py"
            )
            return 1
        print("技能索引页与 SKILL.md 保持同步")
        return 0

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(content, encoding="utf-8")
    skill_count = content.count("\n| [")
    print(f"已生成 {OUTPUT.relative_to(REPO_ROOT)}（{skill_count} 个技能）")
    return 0


if __name__ == "__main__":
    sys.exit(main())
