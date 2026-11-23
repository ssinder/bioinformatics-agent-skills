#!/usr/bin/env python3
"""导出 bio.tools 风格的工具注册素材。

为什么需要这个
--------------
bio.tools（https://bio.tools）是 ELIXIR 维护的生信工具注册表，目前收录
约 3.3 万个工具。它的条目 ID 会被 Bioconda、BioContainers、Galaxy、
Debian Med 复用——也就是说，**注册进 bio.tools 就是注册进一整条分发链**。

每个技能的 ``metadata.edam`` 已经按 EDAM 本体标注好了，这里把它们汇总成
一份符合 bio.tools ``biotoolsSchema`` 的 JSON，便于人工提交或半自动提交。

用法：
    python scripts/export_biotools.py                    # 输出到 stdout
    python scripts/export_biotools.py --output docs/biotools.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from skill_frontmatter import read_frontmatter  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[1]
SKILLS_DIR = REPO_ROOT / "skills"

try:
    import yaml
except ImportError:  # pragma: no cover
    yaml = None


def _skill_entries() -> list[dict]:
    entries: list[dict] = []
    for skill_dir in sorted(SKILLS_DIR.iterdir()):
        if not skill_dir.is_dir() or skill_dir.name.startswith((".", "_")):
            continue
        skill_md = skill_dir / "SKILL.md"
        if not skill_md.is_file():
            continue

        fm = read_frontmatter(skill_md)
        metadata = fm.get("metadata") or {}

        # EDAM 标注以扁平的 edam-<axis> 字符串键存放
        # （规范要求 metadata 是 string→string 映射）。
        edam: dict[str, list[str]] = {}
        for key, value in metadata.items():
            if key.startswith("edam-") and isinstance(value, str):
                axis = key[len("edam-"):].replace("-", "_")
                edam[axis] = [v.strip() for v in value.split(",") if v.strip()]

        description = str(fm.get("description", "")).split(" Use when", 1)[0]
        # 取第一句作为简短描述
        short = description.strip().split(". ")[0].rstrip(".") or description[:120]

        entry = {
            "name": skill_dir.name,
            "description": short,
            "fullDescription": str(fm.get("description", "")),
            "homepage": "https://github.com/ssinder/bioinformatics-agent-skills",
            "repository": "https://github.com/ssinder/bioinformatics-agent-skills",
            "documentation": f"https://github.com/ssinder/bioinformatics-agent-skills/tree/main/skills/{skill_dir.name}",
            "toolType": "Command line tool",
            "topic": edam.get("topic", []),
            "function": {
                "operation": edam.get("operation", []),
                "inputData": edam.get("input_data", []),
                "outputData": edam.get("output_data", []),
                "format": edam.get("output_format", []),
            },
            "language": "Python",
            "license": str(fm.get("license", "MIT")),
            "publication": [],
            "contributors": [
                {
                    "name": str(metadata.get("skill-author", "sinder")),
                    "affiliation": "Independent researcher",
                }
            ],
            "keywords": _keywords(fm),
            "operatingSystem": ["Linux", "Mac OS", "Windows"],
            "dependencies": _dependencies(skill_dir),
        }
        entries.append(entry)
    return entries


def _keywords(fm: dict) -> list[str]:
    """从 description 里抽取关键词（简单的停用词过滤）。"""
    import re

    stop = {
        "use", "when", "this", "the", "and", "for", "with", "from", "that",
        "should", "are", "its", "has", "been", "was", "were", "not", "but",
        "you", "can", "will", "skill", "using", "used", "into", "via", "per",
    }
    text = str(fm.get("description", ""))
    words = re.findall(r"[A-Za-z][A-Za-z0-9\-]{2,}", text.lower())
    seen: list[str] = []
    for w in words:
        if w in stop or w in seen:
            continue
        seen.append(w)
    return seen[:15]


def _dependencies(skill_dir: Path) -> list[dict]:
    """读取 requirements.txt，转成 bio.tools 风格的依赖列表。"""
    req = skill_dir / "requirements.txt"
    if not req.is_file():
        return []
    deps = []
    for line in req.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        deps.append({"name": line, "packageManager": "pip"})
    return deps


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="导出 bio.tools 注册素材")
    parser.add_argument("--output", help="输出文件（默认 stdout）")
    args = parser.parse_args(argv)

    if yaml is None:
        print("缺少 PyYAML，无法解析 frontmatter", file=sys.stderr)
        return 2

    entries = _skill_entries()
    if not entries:
        print("未找到任何技能", file=sys.stderr)
        return 2

    payload = {
        "schema": "biotoolsSchema",
        "tools": entries,
    }
    text = json.dumps(payload, ensure_ascii=False, indent=2)

    if args.output:
        out = Path(args.output)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(text + "\n", encoding="utf-8")
        print(f"已导出 {len(entries)} 个工具条目到 {out}")
    else:
        print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
