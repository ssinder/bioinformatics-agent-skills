"""路径与项目配置管理。

所有技能共享同一套路径解析规则：

1. 显式参数（命令行传入的绝对/相对路径）优先
2. 环境变量 ``BIOPIPE_PROJECT_ROOT`` 次之
3. 当前工作目录（用户实际的数据目录）再次
4. 仓库根目录（仅用于读取随包资源）兜底

这样技能脚本既能被命令行直接调用，也能被 Agent 在任意工作目录中调用，
不再出现 ``/Users/jung/Desktop/trae/tool/...`` 这类硬编码路径。
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional, Union

__all__ = [
    "ProjectConfig",
    "get_config",
    "resolve_path",
    "ensure_dir",
    "get_project_root",
]


def get_project_root() -> Path:
    """返回本仓库根目录（packages/biopipe_core/config.py 的上上级）。"""
    return Path(__file__).resolve().parents[2]


@dataclass
class ProjectConfig:
    """一次技能运行期间共享的配置对象。

    Attributes:
        root: 数据根目录，所有相对路径都基于它解析。
        work_dir: 工作目录（中间产物）。
        output_dir: 结果目录（最终产物）。
        data_dir: 原始数据目录。
        figures_dir: 图片目录。
        config_file: 显式指定的配置文件路径（可选）。
    """

    root: Path
    work_dir: Path
    output_dir: Path
    data_dir: Path
    figures_dir: Path
    config_file: Optional[Path] = None
    _extra: dict[str, Any] = field(default_factory=dict)

    @classmethod
    def from_env(cls, config_file: Optional[Union[str, os.PathLike[str]]] = None) -> ProjectConfig:
        """根据环境变量与可选配置文件构建配置。

        支持的配置文件格式为 JSON 或简单的 ``key=value`` 行，键为
        ``root`` / ``work_dir`` / ``output_dir`` / ``data_dir`` / ``figures_dir``。
        """
        env_root = os.environ.get("BIOPIPE_PROJECT_ROOT")
        root = Path(env_root).expanduser().resolve() if env_root else Path.cwd().resolve()

        extra: dict[str, Any] = {}
        cfg_path: Optional[Path] = None
        if config_file:
            cfg_path = Path(config_file).expanduser().resolve()
            if not cfg_path.exists():
                raise FileNotFoundError(f"配置文件不存在: {cfg_path}")
            extra = _parse_config_file(cfg_path)
            root = Path(extra.get("root", root)).expanduser().resolve()
        else:
            default_cfg = root / "biopipe.yaml"
            if default_cfg.exists():
                extra = _parse_config_file(default_cfg)
                cfg_path = default_cfg
                root = Path(extra.get("root", root)).expanduser().resolve()

        def _sub(key: str, default: str) -> Path:
            value = extra.get(key)
            path = Path(value).expanduser() if value else root / default
            if not path.is_absolute():
                path = root / path
            return path.resolve()

        return cls(
            root=root,
            work_dir=_sub("work_dir", "work"),
            output_dir=_sub("output_dir", "results"),
            data_dir=_sub("data_dir", "data"),
            figures_dir=_sub("figures_dir", "figures"),
            config_file=cfg_path,
            _extra=extra,
        )

    def ensure_dirs(self) -> ProjectConfig:
        """创建所有输出目录，返回自身以便链式调用。"""
        for directory in (self.work_dir, self.output_dir, self.data_dir, self.figures_dir):
            directory.mkdir(parents=True, exist_ok=True)
        return self

    def as_dict(self) -> dict[str, Any]:
        """返回可序列化的配置摘要，用于写入运行元数据。"""
        return {
            "root": str(self.root),
            "work_dir": str(self.work_dir),
            "output_dir": str(self.output_dir),
            "data_dir": str(self.data_dir),
            "figures_dir": str(self.figures_dir),
            "config_file": str(self.config_file) if self.config_file else None,
        }


_CONFIG_CACHE: dict[str, ProjectConfig] = {}


def get_config(config_file: Optional[Union[str, os.PathLike[str]]] = None, refresh: bool = False) -> ProjectConfig:
    """获取（并缓存）当前项目配置。

    Args:
        config_file: 可选的配置文件路径。
        refresh: 为 True 时忽略缓存重新读取。
    """
    key = str(Path(config_file).resolve()) if config_file else "<default>"
    if refresh or key not in _CONFIG_CACHE:
        _CONFIG_CACHE[key] = ProjectConfig.from_env(config_file)
    return _CONFIG_CACHE[key]


def resolve_path(path: Union[str, os.PathLike[str]], base: Optional[Path] = None) -> Path:
    """把用户提供的路径解析为绝对路径。

    Args:
        path: 待解析路径。
        base: 相对路径的基准目录，默认为当前工作目录。

    Returns:
        解析后的绝对路径（不要求存在）。
    """
    p = Path(path).expanduser()
    if p.is_absolute():
        return p.resolve()
    root = base if base is not None else Path.cwd()
    return (root / p).resolve()


def ensure_dir(path: Union[str, os.PathLike[str]]) -> Path:
    """确保目录存在（递归创建）并返回其绝对路径。"""
    p = resolve_path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def _parse_config_file(path: Path) -> dict[str, Any]:
    """解析 JSON 或 ``key=value`` 格式的配置文件。"""
    text = path.read_text(encoding="utf-8")
    if path.suffix.lower() == ".json":
        import json

        try:
            data = json.loads(text)
        except json.JSONDecodeError as exc:
            raise ValueError(f"配置文件 JSON 解析失败: {path}: {exc}") from exc
        return {str(k): v for k, v in data.items()}

    result: dict[str, Any] = {}
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#") or ":" not in line and "=" not in line:
            continue
        if "=" in line and (":" not in line or line.index("=") < line.index(":")):
            key, value = line.split("=", 1)
        else:
            key, value = line.split(":", 1)
        result[key.strip()] = value.strip().strip("\"'")
    return result
