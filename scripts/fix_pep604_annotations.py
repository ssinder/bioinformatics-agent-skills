#!/usr/bin/env python3
"""把 PEP 604 联合类型注解改写为 Python 3.9 兼容的 typing 写法。

问题
----
`from __future__ import annotations` 只让**注解不求值**，但
``typing.get_type_hints()`` 会**强制求值**。在 Python 3.9 上求值
``str | None`` 会抛：

    TypeError: unsupported operand type(s) for |: 'type' and 'type'

因为 PEP 604 的 ``X | Y`` 运算符是 Python 3.10 才引入的。

修复
----
把 ``X | None`` 改写为 ``Optional[X]``，``X | Y`` 改写为 ``Union[X, Y]``，
并自动补上 ``typing`` 中需要的名字导入。

实现要点
--------
用 ``ast.parse`` 精确定位注解节点的位置（``node.lineno`` / ``node.end_lineno``
/ ``col_offset``），而不是用正则——正则无法可靠处理嵌套括号、多行注解、
字符串字面量里的 ``|``。

用法：
    python scripts/fix_pep604_annotations.py --check
    python scripts/fix_pep604_annotations.py --apply
"""

from __future__ import annotations

import argparse
import ast
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]

TARGET_DIRS = ("packages", "skills")

# 已被 future import 覆盖的模块不需要处理（future import 只在模块顶层有效，
# 而我们的目标注解都在模块顶层函数/方法签名里）。
NEEDED_TYPING_NAMES = ("Optional", "Union")


def _is_pep604_binop(node: ast.AST) -> bool:
    """判断节点是否是 PEP 604 的 BinOp(|)。"""
    return isinstance(node, ast.BinOp) and isinstance(node.op, ast.BitOr)


def _collect_pep604_annotations(tree: ast.AST) -> list[ast.BinOp]:
    """收集所有出现在注解位置的 PEP 604 BinOp 节点。

    只看注解位置（arg.annotation、returns、AnnAssign.annotation），
    避免误伤普通表达式里的按位或运算。
    """
    found: list[ast.BinOp] = []

    def visit_annotation(node: ast.AST | None) -> None:
        if node is None:
            return
        # 情况 A：注解本身就是 PEP 604 联合 —— 收集这个 BinOp 节点
        if _is_pep604_binop(node):
            found.append(node)
            return
        # 情况 B：注解是带下标的复合类型（list[...] / tuple[...] /
        # Optional[...] / Union[...]），内部可能藏着 PEP 604 联合。
        # 收集这个**最外层** Subscript，由 _render 递归改写内部。
        # 这里必须收集最外层而不是内部 BinOp——否则 `tuple[str, int | None]`
        # 只会改写成 `tuple[str, Optional[int]]` 却保留了错误的括号结构。
        if isinstance(node, ast.Subscript):
            found.append(node)
            return
        # 情况 C：其他复合注解
        for sub in ast.walk(node):
            if _is_pep604_binop(sub):
                found.append(sub)
                break

    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            args = node.args
            for arg in (
                list(args.posonlyargs)
                + list(args.args)
                + list(args.kwonlyargs)
            ):
                visit_annotation(arg.annotation)
            if args.vararg is not None:
                visit_annotation(args.vararg.annotation)
            if args.kwarg is not None:
                visit_annotation(args.kwarg.annotation)
            visit_annotation(node.returns)
        elif isinstance(node, ast.AnnAssign):
            visit_annotation(node.annotation)
        elif isinstance(node, ast.arg):
            visit_annotation(node.annotation)

    # 去重（同一个节点可能被 visit 多次）
    unique: dict[int, ast.BinOp] = {}
    for node in found:
        unique[id(node)] = node
    return list(unique.values())


def _render(node: ast.AST, source_lines: list[str]) -> str:
    """把 AST 节点还原为源码文本。

    关键点：**只处理 PEP 604 联合类型那一层**，其余部分原样保留原始源码文本。

    做法是用节点的行列位置从源文件里切片，而不是递归重建 AST。
    早期版本用递归渲染 Subscript/Attribute，结果在
    ``os.PathLike[str] | None`` 这种「下标 + 联合」组合下生成了
    ``os.PathLike[str]].PathLike[str]]]`` 这样的垃圾输出。
    源码切片天然保真，不存在这个问题。
    """
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.BitOr):
        # 递归渲染左右子树：内层可能还有嵌套的联合类型
        # （如 ``Optional[str | os.PathLike[str]]``）
        left = _render(node.left, source_lines)
        right = _render(node.right, source_lines)
        # None 优先转成 Optional
        if isinstance(node.right, ast.Constant) and node.right.value is None:
            return f"Optional[{left}]"
        if isinstance(node.left, ast.Constant) and node.left.value is None:
            return f"Optional[{right}]"
        return f"Union[{left}, {right}]"

    if isinstance(node, ast.Subscript):
        # 处理「外层已包 typing、内层仍是 PEP 604」的情况：
        #   Optional[str | os.PathLike[str]] -> Optional[Union[str, os.PathLike[str]]]
        #   tuple[str, int | None]          -> tuple[str, Optional[int]]
        slice_node = node.slice
        if isinstance(slice_node, ast.Index):  # Python 3.8 兼容
            slice_node = slice_node.value

        # slice 可能是 Tuple（tuple[str, int | None]）或单个表达式
        if isinstance(slice_node, ast.Tuple):
            inner_has_pep604 = any(
                _is_pep604_binop(s) for elt in slice_node.elts for s in ast.walk(elt)
            )
            if inner_has_pep604:
                base = _slice_source(node.value, source_lines)
                parts = [_render(elt, source_lines) for elt in slice_node.elts]
                return f"{base}[{', '.join(parts)}]"
            return _slice_source(node, source_lines)

        inner_has_pep604 = any(_is_pep604_binop(s) for s in ast.walk(slice_node))
        if inner_has_pep604:
            base = _slice_source(node.value, source_lines)
            return f"{base}[{_render(slice_node, source_lines)}]"
        return _slice_source(node, source_lines)

    if isinstance(node, ast.Tuple):
        # tuple 的每个元素都可能含 PEP 604 联合，逐个递归渲染
        parts = [_render(elt, source_lines) for elt in node.elts]
        return f"({', '.join(parts)})"

    if isinstance(node, ast.List):
        parts = [_render(elt, source_lines) for elt in node.elts]
        return f"[{', '.join(parts)}]"

    # 非联合节点：直接用源码切片
    return _slice_source(node, source_lines)


def _slice_source(node: ast.AST, source_lines: list[str]) -> str:
    """按节点的行列位置从源文件中切出原始文本。"""
    start_line = node.lineno
    end_line = getattr(node, "end_lineno", start_line)
    start_col = node.col_offset
    end_col = getattr(node, "end_col_offset", None)

    if start_line == end_line:
        line = source_lines[start_line - 1]
        if end_col is None:
            return line[start_col:].rstrip()
        return line[start_col:end_col]

    # 多行：拼接各行
    chunks = [source_lines[start_line - 1][start_col:]]
    for i in range(start_line + 1, end_line):
        chunks.append(source_lines[i - 1])
    if end_col is not None:
        chunks.append(source_lines[end_line - 1][:end_col])
    return "".join(chunks).strip()


def _needs_typing_import(source: str) -> bool:
    """检查是否已有 typing 的 Optional / Union 导入。"""
    try:
        tree = ast.parse(source)
    except SyntaxError:
        return False
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and node.module == "typing":
            for alias in node.names:
                if alias.name in NEEDED_TYPING_NAMES:
                    return True
    return False


def _ensure_typing_import(source: str, names: set[str]) -> str:
    """确保 typing 中需要的名字被导入。

    注意要检查**所有** typing 导入语句，而不是只看第一个——
    早期版本用 ast.walk 找到第一个就 break，会漏掉后续的导入语句，
    导致 "NameError: name 'Union' is not defined"。
    """
    if not names:
        return source

    try:
        tree = ast.parse(source)
    except SyntaxError:
        return source

    existing: set[str] = set()
    import_nodes: list[ast.ImportFrom] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and node.module == "typing":
            existing |= {alias.name for alias in node.names}
            import_nodes.append(node)

    missing = sorted(names - existing)
    if not missing:
        return source

    if import_nodes:
        # 已有 from typing import ...：扩展第一个可用的单行导入；
        # 若是多行形式（带括号），则往里插入名字。
        for stmt in import_nodes:
            lines = source.splitlines(keepends=True)
            if stmt.end_lineno != stmt.lineno:
                # 多行导入：逐个名字插入到括号内第一行之后
                insert_at = stmt.lineno  # 0-based，插到 end_lineno-1 行之前
                for name in reversed(missing):
                    lines.insert(insert_at, f"    {name},\n")
                return "".join(lines)
            old_line = lines[stmt.lineno - 1]
            indent = old_line[: len(old_line) - len(old_line.lstrip())]
            new_names = sorted(existing | set(names))
            lines[stmt.lineno - 1] = f"{indent}from typing import {', '.join(new_names)}\n"
            return "".join(lines)

    # 没有 typing 导入：插到最后一个顶层 import 之后
    tree = ast.parse(source)
    last_import_end = 0
    for node in tree.body:
        if isinstance(node, (ast.Import, ast.ImportFrom)):
            last_import_end = node.end_lineno
    lines = source.splitlines(keepends=True)
    new_import = f"from typing import {', '.join(sorted(names))}\n"
    lines.insert(last_import_end, new_import)
    return "".join(lines)


def _fix_file(path: Path, apply_changes: bool) -> tuple[bool, str]:
    """处理单个文件，返回 (是否有改动, 说明)。"""
    source = path.read_text(encoding="utf-8")
    try:
        tree = ast.parse(source)
    except SyntaxError as exc:
        return False, f"跳过（语法错误）: {exc}"

    nodes = _collect_pep604_annotations(tree)
    if not nodes:
        return False, ""

    source_lines = source.splitlines(keepends=True)
    # 按位置倒序处理，避免前面的替换影响后面的行号
    edits: list[tuple[int, int, int, int, str]] = []  # (start_line, start_col, end_line, end_col, text)

    # 去掉被其他节点包含的节点（嵌套时外层已覆盖内层）。
    # 判定标准：区间**严格重叠**才算包含；同一行上并列的注解区间互不重叠，
    # 必须判为"非包含"，否则会漏改第二个。
    positions = [(n.lineno, n.col_offset, n.end_lineno, n.end_col_offset) for n in nodes]
    included = set()
    for i, outer in enumerate(positions):
        for j, inner in enumerate(positions):
            if i == j:
                continue
            o_start, o_end = (outer[0], outer[1]), (outer[2], outer[3])
            i_start, i_end = (inner[0], inner[1]), (inner[2], inner[3])
            # outer 严格包含 inner（起点更早且终点更晚，或同起点但终点更长）
            if o_start <= i_start and o_end >= i_end and (
                o_start < i_start or o_end > i_end
            ):
                included.add(j)
                break

    for idx, node in enumerate(nodes):
        if idx in included:
            continue
        rendered = _render(node, source_lines)
        edits.append((node.lineno, node.col_offset, node.end_lineno, node.end_col_offset, rendered))

    # 编辑顺序：按 (行号, 列号) **倒序**。
    # 这样同一行上的多处替换从右往左做，右侧替换不会影响左侧的列偏移。
    edits.sort(key=lambda e: (e[0], e[1]), reverse=True)

    new_lines = list(source_lines)
    changed = False
    for start_line, start_col, end_line, end_col, text in edits:
        if start_line == end_line:
            line = new_lines[start_line - 1]
            new_lines[start_line - 1] = line[:start_col] + text + line[end_col:]
        else:
            # 多行注解：整体替换
            head = new_lines[start_line - 1][:start_col]
            tail = new_lines[end_line - 1][end_col:]
            new_lines[start_line - 1 : end_line] = [head + text + tail]
        changed = True

    new_source = "".join(new_lines)

    # 确保 typing 导入齐全
    used_names: set[str] = set()
    for _, _, _, _, text in edits:
        if text.startswith("Optional["):
            used_names.add("Optional")
        elif text.startswith("Union["):
            used_names.add("Union")
    new_source = _ensure_typing_import(new_source, used_names)

    if changed and apply_changes:
        path.write_text(new_source, encoding="utf-8")

    rel = path.relative_to(REPO_ROOT)
    return changed, f"{rel}: 改写 {len(edits)} 处 PEP 604 注解"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--check", action="store_true", help="只报告，不修改")
    mode.add_argument("--apply", action="store_true", help="执行改写")
    args = parser.parse_args()
    apply_changes = args.apply

    files: list[Path] = []
    for target in TARGET_DIRS:
        base = REPO_ROOT / target
        if not base.is_dir():
            continue
        files.extend(sorted(base.rglob("*.py")))

    total = 0
    for path in files:
        changed, message = _fix_file(path, apply_changes)
        if changed:
            total += 1
            if message:
                print(message)

    print(f"\n共 {total} 个文件包含 PEP 604 注解" + ("（已改写）" if apply_changes else "（检查模式，未修改）"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
