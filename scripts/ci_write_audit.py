#!/usr/bin/env python3
"""检查技能脚本是否写到预期目录之外。

目的
----
技能脚本允许写文件（结果表、图片、中间产物），但**不应该**写到用户项目
目录之外的系统位置（如 ``/etc``、``/usr/lib``、``~/.ssh``）。
本脚本列出所有"写文件"的调用点，并标出那些路径变量名看起来像写到
非输出目录的情况，供人工复核。

为什么用 Python 而不是 shell grep
----------------------------------
早先的 workflow 用 shell 正则做这件事：

    PATTERN='open\\s*\\([^)]*["'"']w["'"']|Path\\([^)]*\\)\\.write_'

这段在 ``sh``（dash）里因为单双引号混用而直接语法错误
（``syntax error near unexpected token ')'``）。
正则在 shell 里转义链条太长、跨平台不可控，Python 脚本更可靠。

用法：
    python scripts/ci_write_audit.py            # 报告（默认，不失败）
    python scripts/ci_write_audit.py --strict   # 发现可疑写法就失败
"""

from __future__ import annotations

import argparse
import ast
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]

# 这些路径变量名/目录名表示"预期的输出位置"，命中即视为合规。
SAFE_NAME_HINTS = (
    "output",
    "result",
    "work",
    "data",
    "figure",
    "tmp",
    "temp",
    "cache",
    "out",
    "dir",
    "path",
    "file",
)

# 明确禁止写入的系统路径前缀。
FORBIDDEN_PREFIXES = (
    "/etc",
    "/usr/lib",
    "/usr/bin",
    "/bin",
    "/sbin",
    "/boot",
    "/proc",
    "/sys",
    "/dev",
    "~/.ssh",
    "~/.aws",
    "~/.config/gcloud",
    "C:\\Windows",
    "C:\\Program Files",
)


def _iter_python_files():
    for base in ("skills", "packages"):
        root = REPO_ROOT / base
        if not root.is_dir():
            continue
        for path in sorted(root.rglob("*.py")):  # noqa: UP028
            yield path


def _target_names(node: ast.AST) -> list[str]:
    """提取写文件调用的第一个参数（目标路径）对应的源码文本。"""
    names: list[str] = []

    if isinstance(node, ast.Call):
        func = node.func
        # open(path, "w") / Path(path).write_text(...)
        if isinstance(func, ast.Name) and func.id == "open" and node.args:
            names.append(_expr_text(node.args[0]))
        elif isinstance(func, ast.Attribute) and func.attr in (
            "write_text",
            "write_bytes",
            "write_csv",
            "write_table",
            "to_csv",
            "savefig",
            "save",
            "writelines",
            "dump",
        ):
            # Path(...).write_text(...)：目标是 func.value
            names.append(_expr_text(func.value))
        elif isinstance(func, ast.Attribute) and func.attr in ("write_text", "write_bytes"):
            names.append(_expr_text(func.value))
    return names


def _expr_text(node: ast.AST) -> str:
    """把表达式渲染成便于人工阅读的短文本。"""
    if isinstance(node, ast.Name):
        return node.id
    if isinstance(node, ast.Attribute):
        return f"{_expr_text(node.value)}.{node.attr}"
    if isinstance(node, ast.Call):
        return f"{_expr_text(node.func)}(...)"
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        return repr(node.value)
    if isinstance(node, ast.Subscript):
        return f"{_expr_text(node.value)}[...]"
    try:
        import ast as _ast

        return _ast.unparse(node)[:60]
    except Exception:
        return "<expr>"


def _is_bare_filename(target: str) -> bool:
    """判断目标是否只是一个文件名（不含目录分隔符）。

    这类写操作的完整路径通常由上层 ``output_dir`` / ``results_dir`` 等变量
    拼接而成，单独看文件名无法判断风险，不应报为可疑。
    """
    if not target:
        return True
    # 已包含路径分隔符或路径拼接痕迹的，不是纯文件名
    if any(sep in target for sep in ("/", "\\", "...")):
        return False
    # 形如 a.b.c 的属性访问也不是纯文件名
    if "." in target:
        return False
    return True


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--strict",
        action="store_true",
        help="发现可疑写入路径时以非零退出码结束",
    )
    args = parser.parse_args()

    total_writes = 0
    suspicious: list[str] = []

    for path in _iter_python_files():
        source = path.read_text(encoding="utf-8")
        try:
            tree = ast.parse(source)
        except SyntaxError:
            continue

        rel = path.relative_to(REPO_ROOT)
        for node in ast.walk(tree):
            for target in _target_names(node):
                total_writes += 1
                lowered = target.lower()
                if any(prefix.lower() in lowered for prefix in FORBIDDEN_PREFIXES):
                    suspicious.append(
                        f"{rel}:{getattr(node, 'lineno', '?')} 写入受限系统路径: {target}"
                    )
                    continue
                # 纯文件名（如 `fig`、`table`、`json`）说明完整路径由上层的
                # output_dir 之类变量拼接而来，不构成风险，跳过。
                if _is_bare_filename(target):
                    continue
                if not any(hint in lowered for hint in SAFE_NAME_HINTS):
                    suspicious.append(
                        f"{rel}:{getattr(node, 'lineno', '?')} 目标名不含输出目录关键词，请复核: {target}"
                    )

    print(f"共发现 {total_writes} 处写文件调用")

    if suspicious:
        print(f"\n以下 {len(suspicious)} 处建议人工复核：")
        for item in suspicious:
            print(f"  - {item}")
    else:
        print("所有写文件调用的目标均指向预期输出目录")

    if args.strict and suspicious:
        print("::error::存在需要复核的写入路径")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
