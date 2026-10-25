"""统一的命令行入口骨架。

每个技能脚本都通过 :func:`cli_parser` 与 :func:`run_cli` 获得一致的
命令行体验：``--input`` / ``--output`` / ``--config`` / ``--log-file`` /
``--verbose`` / ``--quiet``，以及统一的异常处理与退出码。
"""

from __future__ import annotations

import argparse
import logging
import sys
from typing import Any, Callable, Sequence

__all__ = ["cli_parser", "run_cli", "SkillError"]


class SkillError(RuntimeError):
    """技能执行过程中可预期的错误，向用户展示友好信息并以非零码退出。"""


def cli_parser(
    prog: str | None = None,
    description: str = "",
    skill: str | None = None,
) -> argparse.ArgumentParser:
    """创建带公共参数的解析器。

    Args:
        prog: 程序名，通常传 ``技能名``。
        description: 命令描述。
        skill: 技能名，用于日志标注与元数据。

    Returns:
        配置好的 ``ArgumentParser``。
    """
    parser = argparse.ArgumentParser(
        prog=prog,
        description=description,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--config",
        metavar="PATH",
        help="项目配置文件路径（JSON 或 key=value 格式）",
    )
    parser.add_argument(
        "--log-file",
        metavar="PATH",
        help="日志输出文件路径",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="输出 DEBUG 级别日志",
    )
    parser.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="仅输出警告与错误",
    )
    parser.add_argument(
        "--version",
        action="store_true",
        help="显示版本信息",
    )
    return parser


def run_cli(
    main_func: Callable[[argparse.Namespace], Any],
    argv: Sequence[str] | None = None,
    parser: argparse.ArgumentParser | None = None,
    skill: str = "biopipe",
    version: str = "1.0.0",
) -> int:
    """执行技能主函数并统一处理异常与退出码。

    Args:
        main_func: 接收 ``argparse.Namespace`` 并返回退出码或 None 的函数。
        argv: 命令行参数列表，默认取 ``sys.argv[1:]``。
        parser: 已创建的解析器，为 None 时用 :func:`cli_parser` 创建。
        skill: 技能名。
        version: 版本号。

    Returns:
        进程退出码（0 表示成功）。
    """
    parser = parser or cli_parser(prog=skill, description=skill, skill=skill)
    args = parser.parse_args(argv)

    if getattr(args, "version", False):
        print(f"{skill} {version}")
        return 0

    level = logging.DEBUG if args.verbose else logging.INFO
    if args.quiet:
        level = logging.WARNING

    from .logging_utils import setup_logging

    logger = setup_logging(f"biopipe.{skill}", level=level, log_file=args.log_file, quiet=args.quiet)

    try:
        result = main_func(args)
    except SkillError as exc:
        logger.error("%s", exc)
        return 1
    except KeyboardInterrupt:
        logger.warning("已中断")
        return 130
    except FileNotFoundError as exc:
        logger.error("文件不存在: %s", exc)
        return 1
    except Exception as exc:  # noqa: BLE001
        logger.exception("执行失败: %s", exc)
        return 1

    return int(result) if isinstance(result, int) else 0


def bootstrap_cli(
    description: str,
    skill: str,
    main_func: Callable[[argparse.Namespace], Any],
    version: str = "1.0.0",
):
    """在 ``__main__`` 中使用的标准入口。

    用法::

        if __name__ == "__main__":
            sys.exit(bootstrap_cli("描述", "skill-name", main))
    """
    parser = cli_parser(prog=skill, description=description, skill=skill)
    return run_cli(main_func, parser=parser, skill=skill, version=version)
