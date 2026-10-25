"""可选依赖的优雅降级支持。

背景：许多技能依赖 numpy / scipy / pandas / matplotlib 等科学计算栈，
但这些库体积大且并非所有用户都安装。技能脚本的 ``--help`` 属于纯元信息
查询，不应因为缺少科学计算栈而无法执行。

本模块提供两种能力：

1. :func:`optional_import` —— 尝试导入模块，失败返回 ``None`` 而不是抛异常；
2. :func:`require` —— 在真正需要该依赖的位置调用，若依赖缺失则抛出
   带明确安装指引的 :class:`MissingDependencyError`。

用法::

    from biopipe_core.optional_deps import optional_import, require

    np = optional_import("numpy")
    if np is None:
        raise MissingDependencyError("numpy", "pip install numpy")

    # 或者在关键路径上直接：
    np = require("numpy", hint="差异表达分析需要 numpy")
"""

from __future__ import annotations

import importlib
from types import ModuleType

__all__ = [
    "MissingDependencyError",
    "optional_import",
    "require",
    "is_available",
    "install_hint",
]


class MissingDependencyError(ImportError):
    """必需的第三方依赖缺失。"""

    def __init__(self, package: str, hint: str = "", extra: str = "") -> None:
        self.package = package
        self.hint = hint
        self.extra = extra
        message = f"缺少依赖: {package}"
        if hint:
            message += f"（{hint}）"
        message += f"\n请先安装: pip install {package}"
        if extra:
            message += f"  或: pip install 'biopipe-agent-skills[{extra}]'"
        super().__init__(message)


#: 常见科学计算依赖的安装提示
install_hint = {
    "numpy": "数值计算基础库，几乎所有技能都需要",
    "scipy": "科学计算库（统计检验、稀疏矩阵、优化）",
    "pandas": "表格数据处理",
    "matplotlib": "绘图",
    "seaborn": "统计图形",
    "scanpy": "单细胞分析主库",
    "anndata": "单细胞数据容器",
    "leidenalg": "Leiden 聚类算法",
    "MDAnalysis": "分子动力学轨迹读取（权威参考实现）",
    "biopython": "生物序列处理",
    "h5py": "HDF5 文件读写",
    "feedparser": "RSS/Atom 解析",
}


def optional_import(name: str) -> ModuleType | None:
    """尝试导入模块，失败时返回 ``None`` 而不抛异常。

    Args:
        name: 模块名，如 ``"numpy"``、``"scipy.sparse"``。

    Returns:
        导入成功返回模块对象，失败返回 ``None``。
    """
    try:
        return importlib.import_module(name)
    except ImportError:
        return None


def require(name: str, hint: str = "", extra: str = "") -> ModuleType:
    """导入模块，缺失时抛出带安装指引的异常。

    Args:
        name: 模块名。
        hint: 说明该依赖在当前上下文中用于什么，便于用户理解。
        extra: 可选依赖组名（对应 ``pyproject.toml`` 中的 extras）。

    Returns:
        导入成功的模块对象。

    Raises:
        MissingDependencyError: 依赖缺失时抛出。
    """
    try:
        return importlib.import_module(name)
    except ImportError as exc:
        raise MissingDependencyError(name, hint=hint, extra=extra) from exc


def is_available(name: str) -> bool:
    """检查模块是否可用。"""
    try:
        importlib.import_module(name)
    except ImportError:
        return False
    return True
