"""biopipe_core - 共享工具层

为所有生物信息学技能提供统一的路径配置、日志、HTTP 重试与命令行入口。
所有技能脚本都应通过本模块获取路径与日志，禁止在脚本内硬编码绝对路径。

Usage:
    from biopipe_core import get_logger, resolve_path, http_session, cli_parser
"""

from .config import (
    ProjectConfig,
    get_config,
    resolve_path,
    ensure_dir,
    get_project_root,
)
from .logging_utils import get_logger, setup_logging
from .http import http_session, HttpError, get_json, post_json
from .cli import cli_parser, run_cli, SkillError
from .io_utils import (
    read_fasta,
    write_fasta,
    read_csv_safe,
    write_csv,
    ensure_unicode_csv,
)

__version__ = "1.0.0"

__all__ = [
    "ProjectConfig",
    "get_config",
    "resolve_path",
    "ensure_dir",
    "get_project_root",
    "get_logger",
    "setup_logging",
    "http_session",
    "HttpError",
    "get_json",
    "post_json",
    "cli_parser",
    "run_cli",
    "SkillError",
    "read_fasta",
    "write_fasta",
    "read_csv_safe",
    "write_csv",
    "ensure_unicode_csv",
    "__version__",
]
