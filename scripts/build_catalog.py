#!/usr/bin/env python3
"""从各技能的 SKILL.md 生成技能目录文档。

用法：
    python scripts/build_catalog.py             # 输出到 stdout
    python scripts/build_catalog.py -o docs/skills.md
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SKILLS_DIR = REPO_ROOT / "skills"

sys.path.insert(0, str(REPO_ROOT / "packages"))

HEADER = """# 技能目录 / Skills Catalog

<!-- 本文件由 scripts/build_catalog.py 自动生成，请勿手工编辑 -->

共 **{count}** 个技能,按分类分组。

"""


def _parse_frontmatter(text: str) -> dict:
    """解析 SKILL.md 的 YAML frontmatter。"""
    if not text.startswith("---"):
        return {}
    parts = text.split("---", 2)
    if len(parts) < 3:
        return {}
    try:
        import yaml

        data = yaml.safe_load(parts[1])
        return data if isinstance(data, dict) else {}
    except Exception:  # noqa: BLE001
        return {}


def main(argv: list[str] | None = None) -> int:
    """生成目录文档。"""
    parser = argparse.ArgumentParser(description="生成技能目录")
    parser.add_argument("-o", "--output", help="输出文件路径，默认打印到 stdout")
    args = parser.parse_args(argv)

    skills_dir = SKILLS_DIR
    if not skills_dir.exists():
        print(f"技能目录不存在: {skills_dir}", file=sys.stderr)
        return 1

    try:
        from biopipe.registry import REGISTRY
    except ImportError:
        REGISTRY = {}

    entries = []
    for skill_dir in sorted(p for p in skills_dir.iterdir() if p.is_dir()):
        skill_md = skill_dir / "SKILL.md"
        if not skill_md.exists():
            continue
        data = _parse_frontmatter(skill_md.read_text(encoding="utf-8"))
        spec = REGISTRY.get(skill_dir.name)
        entries.append(
            {
                "name": skill_dir.name,
                "description": data.get("description", ""),
                "category": spec.category if spec else "其他",
                "scripts": len(list((skill_dir / "scripts").glob("*.py"))) if (skill_dir / "scripts").exists() else 0,
                "refs": len(list((skill_dir / "references").glob("*.md"))) if (skill_dir / "references").exists() else 0,
            }
        )

    categories: dict[str, list[dict]] = {}
    for entry in entries:
        categories.setdefault(entry["category"], []).append(entry)

    lines = [HEADER.format(count=len(entries))]
    for category in sorted(categories):
        lines.append(f"## {category}\n")
        lines.append("| 技能 | 描述 | 脚本 | 参考文档 |")
        lines.append("|------|------|------|---------|")
        for entry in categories[category]:
            desc = entry["description"].split(".")[0].strip()
            if len(desc) > 110:
                desc = desc[:107] + "..."
            lines.append(
                f"| [`{entry['name']}`](../skills/{entry['name']}/) | {desc} | "
                f"{entry['scripts']} | {entry['refs']} |"
            )
        lines.append("")

    output = "\n".join(lines)
    if args.output:
        out = Path(args.output)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(output, encoding="utf-8")
        print(f"已生成: {out}（{len(entries)} 个技能）")
    else:
        print(output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
