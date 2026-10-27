"""Unified command-line entry point for every skill.

Usage::

    biopipe --help                 # list every skill
    biopipe --list                 # print every skill as a table
    biopipe bulk-rnaseq --help     # show the usage of one skill
    biopipe bulk-rnaseq --run ...  # run a skill
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
        description=(
            "Unified command-line entry point for research skills — "
            "RNA-seq, single-cell, structure analysis, molecular dynamics, literature mining"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  biopipe --list                      list every skill\n"
            "  biopipe bulk-rnaseq --help          show the usage of one skill\n"
            "  biopipe scientific-figures --run    run a skill script\n"
        ),
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="list every registered skill as a table",
    )
    parser.add_argument(
        "--version",
        action="store_true",
        help="show the version number",
    )
    parser.add_argument(
        "skill",
        nargs="?",
        help="skill name",
    )
    parser.add_argument(
        "args",
        nargs=argparse.REMAINDER,
        help="arguments forwarded to the skill script",
    )
    return parser


def _print_skill_list() -> None:
    print("Available skills (overview)\n")
    for category in sorted(_CATEGORIES):
        print(f"[{category}]")
        for name in _CATEGORIES[category]:
            spec = REGISTRY[name]
            print(f"  {name:<24} {spec.summary}")
        print()
    print("Usage: biopipe <skill-name> --help")


def main(argv: Optional[list[str]] = None) -> int:
    """Main entry point of the biopipe command."""
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
        print(f"Unknown skill: {args.skill}", file=sys.stderr)
        print(f"Available skills: {', '.join(sorted(REGISTRY))}", file=sys.stderr)
        return 2

    entry = spec.entry_path
    if not entry.exists():
        print(f"Skill entry script not found: {entry}", file=sys.stderr)
        return 1

    cmd = [sys.executable, str(entry), *args.args]
    try:
        return subprocess.call(cmd)
    except KeyboardInterrupt:
        return 130


if __name__ == "__main__":
    sys.exit(main())
