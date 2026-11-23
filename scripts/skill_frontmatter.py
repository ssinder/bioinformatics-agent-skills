#!/usr/bin/env python3
"""SKILL.md frontmatter 的健壮解析。

为什么需要这个模块
------------------
最直觉的写法是 ``text.split("---", 2)``，但它是**错的**：

    ---
    name: bulk-rnaseq
    ...
    ---

    ## Parameters

    | Parameter | Default |
    |---|---|        <-- 这一行的 '---' 会被误当成 frontmatter 结束标记
    | `--input` | required |

``split("---", 2)`` 会在表格分隔行处切开，把大半个正文当成 YAML 去解析，
于是抛出一个和真实原因毫无关系的错误（"could not find expected ':'"），
或者更糟——静默地把 description 解析成空字符串，让人以为描述丢了。

正确做法：只把**独占一行**的 ``---`` 当作分隔符。本模块提供统一实现，
供校验器、导出脚本、描述改写脚本与测试共用。

用法：
    from skill_frontmatter import split_frontmatter, read_frontmatter

    fm_text, body = split_frontmatter(path)
    data = read_frontmatter(path)
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

try:
    import yaml
except ImportError:  # pragma: no cover
    yaml = None

__all__ = ["split_frontmatter", "read_frontmatter", "FrontmatterError"]


class FrontmatterError(ValueError):
    """frontmatter 缺失或格式错误。"""


def split_frontmatter(text: str) -> tuple[str, str]:
    """把文本拆成 (frontmatter 原文, 正文)。

    Raises:
        FrontmatterError: 未以 ``---`` 开头，或找不到结束标记。
    """
    lines = text.split("\n")

    if not lines or lines[0].strip() != "---":
        raise FrontmatterError("文件未以 '---' 开头，缺少 frontmatter")

    for i in range(1, len(lines)):
        if lines[i].strip() == "---":
            return "\n".join(lines[1:i]), "\n".join(lines[i + 1:])

    raise FrontmatterError("frontmatter 未正确闭合（缺少独占一行的 '---'）")


def read_frontmatter(path: Path) -> dict[str, Any]:
    """读取并解析 SKILL.md 的 frontmatter，返回字典。

    解析失败时返回空字典而不是抛异常——调用方通常是在做检查，
    空的 frontmatter 会被后续的"缺少必填字段"检查捕获，报错信息
    反而更准确。
    """
    if yaml is None:
        return {}
    try:
        fm_text, _ = split_frontmatter(path.read_text(encoding="utf-8"))
    except (FrontmatterError, OSError, UnicodeDecodeError):
        return {}
    try:
        data = yaml.safe_load(fm_text)
    except Exception:
        return {}
    return data if isinstance(data, dict) else {}
