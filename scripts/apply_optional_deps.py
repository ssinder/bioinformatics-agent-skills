#!/usr/bin/env python3
"""把技能脚本的顶层重型依赖导入改为可降级导入。

背景：技能脚本在模块顶层 ``import numpy`` 会导致 ``--help`` 在未安装
科学计算栈的环境中直接崩溃。正确做法是顶层用 ``optional_import`` 取到
``None``（缺失时不抛异常），再在真正执行计算的入口处显式检查。

本脚本用 AST 精确定位，保证只改动目标 import 语句与紧随其后的插入点：

1. 把顶层 ``import numpy as np`` 替换为
   ``np = optional_import("numpy")`` + 一个 ``_require_np()`` 辅助函数；
2. 不改动任何函数体 —— 函数里引用的仍是模块级 ``np``，
   缺失时为 ``None``，由入口的依赖检查统一给出清晰报错。

用法：
    python scripts/apply_optional_deps.py --check          # 只报告
    python scripts/apply_optional_deps.py <file> [<file>...]
"""

from __future__ import annotations

import argparse
import ast
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]

#: 顶层 import 的模块名 -> (别名, 用途提示)
TARGETS = {
    "numpy": ("np", "数值计算基础库"),
    "pandas": ("pd", "表格数据处理"),
    "scipy.sparse": ("sp", "稀疏矩阵运算"),
}


def _find_import(tree: ast.Module, module: str) -> ast.stmt | None:
    """在顶层找到 import 指定模块的语句节点。"""
    for node in tree.body:
        if isinstance(node, ast.Import):
            for alias in node.names:
                if alias.name == module:
                    return node
        elif isinstance(node, ast.ImportFrom):
            if node.module == module and node.level == 0:
                return node
    return None


def _build_block(alias: str, module: str, hint: str) -> str:
    """构造替换用的代码块。"""
    req = f"_require_{alias}"
    return (
        f"{alias} = optional_import({module!r})\n"
        f"{alias}_HINT = {hint!r}\n"
        f"\n"
        f"\n"
        f"def {req}():\n"
        f'    """确保 {module} 可用，否则抛出带安装指引的异常。"""\n'
        f"    if {alias} is None:\n"
        f"        from biopipe_core.optional_deps import MissingDependencyError\n"
        f"\n"
        f"        raise MissingDependencyError({module!r}, hint={alias}_HINT)\n"
        f"    return {alias}\n"
        f"\n"
        f"\n"
    )


#: 插入到第一个替换块之前的导入语句
IMPORT_LINE = "from biopipe_core.optional_deps import optional_import  # noqa: E402\n\n"


def transform(path: Path, check_only: bool = False) -> list[str]:
    """转换单个文件，返回发生的改动描述列表。"""
    src = path.read_text(encoding="utf-8")
    tree = ast.parse(src)
    lines = src.splitlines(keepends=True)

    # 收集替换任务：(起始行号0based, 结束行号0based不含, 文本)
    edits: list[tuple[int, int, str]] = []
    notes: list[str] = []

    for module, (alias, hint) in TARGETS.items():
        node = _find_import(tree, module)
        if node is None:
            continue
        # 确认该 import 用的别名与预期一致
        actual = alias
        if isinstance(node, ast.Import):
            for a in node.names:
                if a.name == module:
                    actual = a.asname or module.split(".")[0]
        elif isinstance(node, ast.ImportFrom):
            for a in node.names:
                actual = a.asname or a.name
        if actual != alias:
            notes.append(f"{path.name}: {module} 别名为 {actual}，跳过")
            continue

        block = _build_block(alias, module, hint)
        edits.append((node.lineno - 1, node.end_lineno, block))
        notes.append(f"{path.name}: {module} -> optional_import")

    if not edits:
        return []

    # 从后往前替换，避免行号漂移
    for start, end, text in sorted(edits, key=lambda e: -e[0]):
        lines[start:end] = [text]

    # 在第一个替换块之前插入 optional_import 的导入语句
    if edits:
        first_start = min(e[0] for e in edits)
        # 避免重复插入
        if not lines[max(first_start - 3, 0) : first_start] or IMPORT_LINE.strip() not in "".join(
            lines[max(first_start - 3, 0) : first_start]
        ):
            lines[first_start:first_start] = [IMPORT_LINE]

    new_src = "".join(lines)
    if new_src != src:
        try:
            ast.parse(new_src)
        except SyntaxError as exc:
            raise SystemExit(f"{path}: 转换后语法错误: {exc}") from exc
        if not check_only:
            path.write_text(new_src, encoding="utf-8")

    return notes


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="把顶层重型依赖改为可降级导入")
    parser.add_argument("paths", nargs="*", help="要处理的 Python 文件")
    parser.add_argument("--check", action="store_true", help="只报告不修改")
    args = parser.parse_args(argv)

    targets = [Path(p) for p in args.paths]
    if not targets:
        targets = sorted(REPO_ROOT.glob("skills/*/scripts/*.py"))

    total = 0
    for target in targets:
        if not target.exists():
            continue
        notes = transform(target, check_only=args.check)
        for note in notes:
            print(note)
        total += len(notes)

    verb = "需要修改" if args.check else "已修改"
    print(f"\n共 {total} 处依赖导入{verb}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
