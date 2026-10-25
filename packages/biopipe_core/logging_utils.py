"""统一日志工具。

替代散落在 217 个脚本中的 3408 处 ``print()`` 调用，提供：

- 控制台 + 文件双通道输出
- 统一的等级前缀与时间戳
- 技能名自动标注
- 可配置的静默模式（``--quiet``），便于 Agent 解析输出
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import Optional, Union

__all__ = ["get_logger", "setup_logging"]

_FORMAT = "%(asctime)s | %(levelname)-7s | %(name)s | %(message)s"
_DATEFMT = "%Y-%m-%d %H:%M:%S"

_CONFIGURED: set[str] = set()


def setup_logging(
    name: str,
    level: int = logging.INFO,
    log_file: Optional[Union[str, Path]] = None,
    quiet: bool = False,
) -> logging.Logger:
    """配置并返回一个 logger。

    Args:
        name: logger 名称，通常为 ``biopipe.<skill>``。
        level: 日志等级。
        log_file: 可选的日志文件路径，父目录会自动创建。
        quiet: 为 True 时控制台不输出 INFO 及以下内容（错误仍显示）。

    Returns:
        配置好的 ``logging.Logger`` 实例。
    """
    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.propagate = False

    if name in _CONFIGURED:
        return logger

    console = logging.StreamHandler(stream=sys.stderr)
    console.setLevel(logging.ERROR if quiet else level)
    console.setFormatter(logging.Formatter(_FORMAT, datefmt=_DATEFMT))
    logger.addHandler(console)

    if log_file:
        path = Path(log_file).expanduser()
        path.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(path, encoding="utf-8")
        file_handler.setLevel(level)
        file_handler.setFormatter(logging.Formatter(_FORMAT, datefmt=_DATEFMT))
        logger.addHandler(file_handler)

    _CONFIGURED.add(name)
    return logger


def get_logger(name: str = "biopipe", level: int = logging.INFO, **kwargs) -> logging.Logger:
    """获取 logger 的便捷封装。

    Args:
        name: logger 名称。
        level: 日志等级，默认 INFO。
        **kwargs: 透传给 :func:`setup_logging`（如 ``log_file``、``quiet``）。
    """
    return setup_logging(name, level=level, **kwargs)
