"""biopipe 统一命令行入口。

用法::

    biopipe --help                 # 列出全部技能
    biopipe --list                 # 以表格列出全部技能
    biopipe bulk-rnaseq --help     # 查看技能用法
    biopipe bulk-rnaseq --run ...  # 执行技能
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from typing import Optional

from .registry import REGISTRY, get_skill

__all__ = ["main"]

_CATEGORIES: dict[str, list[str]] = {}
for _spec in REGISTRY.values():
    _CATEGORIES.setdefault(_spec.category, []).append(_spec.name)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="biopipe",
        description="科研技能统一命令行入口 —— RNA-seq、单细胞、结构分析、分子动力学、文献挖掘",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "示例:\n"
            "  biopipe --list                      列出全部技能\n"
            "  biopipe bulk-rnaseq --help          查看某个技能的用法\n"
            "  biopipe scientific-figures --run    执行技能脚本\n"
        ),
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="以表格形式列出全部已注册技能",
    )
    parser.add_argument(
        "--version",
        action="store_true",
        help="显示版本号",
    )
    parser.add_argument(
        "skill",
        nargs="?",
        help="技能名称",
    )
    parser.add_argument(
        "args",
        nargs=argparse.REMAINDER,
        help="传递给技能脚本的参数",
    )
    return parser


def _print_skill_list() -> None:
    print("可用技能 (可用技能总览)\n")
    for category in sorted(_CATEGORIES):
        print(f"【{category}】")
        for name in _CATEGORIES[category]:
            spec = REGISTRY[name]
            print(f"  {name:<24} {spec.summary}")
        print()
    print("用法: biopipe <技能名> --help")


def main(argv: Optional[list[str]] = None) -> int:
    """biopipe 命令主函数。"""
    parser = _build_parser()
    args = parser.parse_args(argv)

    if args.version:
        from . import __version__

        print(f"biopipe {__version__}")
        return 0

    if args.list or not args.skill:
        _print_skill_list()
        return 0

    spec = get_skill(args.skill)
    if spec is None:
        print(f"未知技能: {args.skill}", file=sys.stderr)
        print(f"可用技能: {', '.join(sorted(REGISTRY))}", file=sys.stderr)
        return 2

    entry = spec.entry_path
    if not entry.exists():
        print(f"技能入口脚本不存在: {entry}", file=sys.stderr)
        return 1

    cmd = [sys.executable, str(entry), *args.args]
    try:
        return subprocess.call(cmd)
    except KeyboardInterrupt:
        return 130


if __name__ == "__main__":
    sys.exit(main())
