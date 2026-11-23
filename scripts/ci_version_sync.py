#!/usr/bin/env python3
"""校验版本号在多处手工维护的位置之间保持一致。

背景
----
本仓库的版本号分散在四个地方，且此前没有任何检查约束它们：

1. ``pyproject.toml`` 的 ``[project] version`` —— PyPI 包的实际版本；
2. ``plugin.json`` 的 ``version`` —— Agent Plugins 清单的版本；
3. ``CITATION.cff`` 的 ``version`` —— 引用元数据的版本；
4. 各技能 ``SKILL.md`` frontmatter 的 ``metadata.version`` —— 单个技能的版本。

这四处都是**手工维护**的。一旦有人只改了其中一处，发布出去的
wheel、插件清单和引用信息就会互相矛盾——用户 ``pip install`` 到的
版本与仓库页面上显示的版本对不上，而这类不一致很难被察觉。

同时校验 ``CITATION.cff`` 的 ``date-released`` 不晚于今天，避免
遗留旧日期。

用法：
    python scripts/ci_version_sync.py
    python scripts/ci_version_sync.py --set 1.2.0   # 统一改写为指定版本
"""

from __future__ import annotations

import argparse
import datetime as dt
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SKILLS_DIR = REPO_ROOT / "skills"

PYPROJECT = REPO_ROOT / "pyproject.toml"
PLUGIN_JSON = REPO_ROOT / "plugin.json"
CITATION_CFF = REPO_ROOT / "CITATION.cff"

# 技能自身的版本允许落后于仓库版本——技能是独立演进的，
# 只要求格式合法即可（仓库版本与技能版本本来就不是同一件事）。
SKILL_VERSION_RE = re.compile(r"^\d+\.\d+(\.\d+)?$")


def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def _pyproject_version(text: str) -> str | None:
    m = re.search(r"^\[project\]\s*$(.*?)(?=^\[|\Z)", text, re.MULTILINE | re.DOTALL)
    if not m:
        return None
    v = re.search(r'^version\s*=\s*"([^"]+)"', m.group(1), re.MULTILINE)
    return v.group(1) if v else None


def _plugin_version(text: str) -> str | None:
    m = re.search(r'"version"\s*:\s*"([^"]+)"', text)
    return m.group(1) if m else None


def _citation_version(text: str) -> str | None:
    m = re.search(r"^version:\s*\"?([^\"\n]+)\"?\s*$", text, re.MULTILINE)
    return m.group(1).strip() if m else None


def _skill_versions() -> list[tuple[str, str | None]]:
    out: list[tuple[str, str | None]] = []
    for skill_dir in sorted(SKILLS_DIR.iterdir()) if SKILLS_DIR.is_dir() else []:
        if not skill_dir.is_dir() or skill_dir.name.startswith((".", "_")):
            continue
        skill_md = skill_dir / "SKILL.md"
        if not skill_md.is_file():
            continue
        text = _read(skill_md)
        m = re.search(r"^\s+version:\s*\"?([^\"\n]+)\"?\s*$", text, re.MULTILINE)
        out.append((skill_dir.name, m.group(1).strip() if m else None))
    return out


def _rewrite(path: Path, pattern: str, replacement: str) -> bool:
    """按正则改写文件中的版本字段，返回是否发生了改动。"""
    text = _read(path)
    new_text, n = re.subn(pattern, replacement, text, count=1, flags=re.MULTILINE)
    if n and new_text != text:
        path.write_text(new_text, encoding="utf-8")
        return True
    return False


def _set_version(version: str) -> int:
    changed = False
    changed |= _rewrite(
        PYPROJECT, r'^version\s*=\s*"[^"]+"', f'version = "{version}"'
    )
    changed |= _rewrite(
        PLUGIN_JSON, r'"version"\s*:\s*"[^"]+"', f'"version": "{version}"'
    )
    changed |= _rewrite(
        CITATION_CFF, r'^version:.*$', f'version: {version}'
    )
    print(f"已将仓库版本统一设为 {version}" + ("（有文件被改动）" if changed else "（无需改动）"))
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="校验版本号在各文件之间一致")
    parser.add_argument("--set", dest="set_version", help="把各处版本统一改写为指定值")
    args = parser.parse_args(argv)

    if args.set_version:
        return _set_version(args.set_version)

    problems: list[str] = []

    for path in (PYPROJECT, PLUGIN_JSON, CITATION_CFF):
        if not path.is_file():
            problems.append(f"缺少文件: {path.relative_to(REPO_ROOT)}")

    if problems:
        for p in problems:
            print(f"  - {p}")
        print("::error::版本同步检查失败")
        return 1

    pyproject_v = _pyproject_version(_read(PYPROJECT))
    plugin_v = _plugin_version(_read(PLUGIN_JSON))
    citation_v = _citation_version(_read(CITATION_CFF))

    if not pyproject_v:
        problems.append("pyproject.toml 中解析不到 [project].version")
    if not plugin_v:
        problems.append("plugin.json 中解析不到 version 字段")
    if not citation_v:
        problems.append("CITATION.cff 中解析不到 version 字段")

    # 三处必须完全一致
    versions = {
        "pyproject.toml": pyproject_v,
        "plugin.json": plugin_v,
        "CITATION.cff": citation_v,
    }
    distinct = {v for v in versions.values() if v}
    if len(distinct) > 1:
        detail = ", ".join(f"{k}={v}" for k, v in versions.items())
        problems.append(
            f"版本号在多个文件之间不一致: {detail}。"
            "请用 `python scripts/ci_version_sync.py --set <version>` 统一。"
        )

    # 技能自身的版本只需格式合法
    for name, version in _skill_versions():
        if version is None:
            problems.append(f"技能 {name} 的 SKILL.md 缺少 metadata.version")
        elif not SKILL_VERSION_RE.match(version):
            problems.append(f"技能 {name} 的 metadata.version '{version}' 不是合法版本号")

    # CITATION.cff 的发布日期不应是未来日期
    m = re.search(r'^date-released:\s*"?([\d\-]+)"?', _read(CITATION_CFF), re.MULTILINE)
    if m:
        try:
            released = dt.date.fromisoformat(m.group(1))
            if released > dt.date.today():
                problems.append(
                    f"CITATION.cff 的 date-released ({released}) 晚于今天，"
                    "请在正式发布后再更新。"
                )
        except ValueError:
            problems.append(f"CITATION.cff 的 date-released '{m.group(1)}' 不是合法日期")

    if problems:
        print("版本同步检查未通过：")
        for p in problems:
            print(f"  - {p}")
        print("::error::版本同步检查失败")
        return 1

    print(f"版本同步检查通过：仓库版本 {pyproject_v}，{len(_skill_versions())} 个技能版本格式合法")
    return 0


if __name__ == "__main__":
    sys.exit(main())
