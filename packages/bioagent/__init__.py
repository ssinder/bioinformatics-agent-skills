"""bioagent - 科研技能统一命令行入口。

一个命令调度全部技能::

    bioagent --help              # 列出所有技能
    bioagent rnaseq --help       # 查看某个技能的用法
    bioagent scrna run --help    # 查看某个技能子命令的用法
"""

from .registry import REGISTRY, get_skill, list_skills, skill_names

__version__ = "1.0.0"

__all__ = ["REGISTRY", "get_skill", "list_skills", "skill_names", "__version__"]
