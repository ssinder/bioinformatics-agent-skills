"""构建钩子：把 ``skills/`` 打进 wheel。

为什么需要这个
--------------
本仓库有两类用户：

* **把技能装进 AI 宿主**（Claude Code、Codex、Cursor…）——用
  ``npx skills add`` 或直接复制 ``skills/`` 目录，不需要 pip；
* **在命令行里用**——``pip install biopipe-agent-skills`` 后跑
  ``bioagent bulk-rnaseq --run ...``。

第二类必须能在**装好的包里**找到技能脚本。而 ``skills/`` 在源码仓库里
位于 ``packages/`` 之外，setuptools 默认不会把它放进 wheel。结果是
``pip install`` 之后 ``bioagent --list`` 能列出 17 个技能，却一个都跑不了
——一个只有列表没有功能的空壳。

因此这里在 ``build_py`` 阶段把 ``skills/`` 复制进
``build/lib/bioagent/skills/``，让它作为 package data 随 wheel 分发。

装好后，``bioagent.registry.skills_root()`` 会优先找到这个位置；
源码检出（editable install）则仍走仓库里的 ``skills/``。

项目元数据全部在 ``pyproject.toml`` 里，这里只做构建期的文件搬运。
"""

from __future__ import annotations

import shutil
from pathlib import Path

from setuptools import setup
from setuptools.command.build_py import build_py as _build_py

REPO_ROOT = Path(__file__).resolve().parent
SKILLS_DIR = REPO_ROOT / "skills"
PACKAGE_DIR = REPO_ROOT / "packages" / "bioagent"

# 复制时忽略的目录：版本控制与构建缓存不需要进 wheel。
IGNORE_DIRS = {"__pycache__", ".git", ".pytest_cache", ".mypy_cache", ".ruff_cache"}
IGNORE_SUFFIXES = {".pyc", ".pyo"}


class BuildPyWithSkills(_build_py):
    """在标准 build_py 之后，把 skills/ 复制进包内。"""

    def run(self) -> None:
        super().run()
        self._copy_skills()

    def _copy_skills(self) -> None:
        if not SKILLS_DIR.is_dir():
            return

        target_root = Path(self.build_lib) / "bioagent" / "skills"
        if target_root.exists():
            shutil.rmtree(target_root)

        copied = 0
        for src in sorted(SKILLS_DIR.iterdir()):
            if not src.is_dir() or src.name.startswith((".", "_")):
                continue
            dst = target_root / src.name
            shutil.copytree(
                src,
                dst,
                ignore=shutil.ignore_patterns(
                    *IGNORE_DIRS,
                    *[f"*{suffix}" for suffix in IGNORE_SUFFIXES],
                ),
            )
            copied += 1

        if copied:
            print(f"[setup.py] packaged {copied} skills into bioagent/skills/")


setup(cmdclass={"build_py": BuildPyWithSkills})
