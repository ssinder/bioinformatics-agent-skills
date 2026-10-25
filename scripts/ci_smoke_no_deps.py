#!/usr/bin/env python3
"""无科学计算栈环境下的技能入口冒烟测试。

目的
----
`--help` 属于**元信息查询**：用户（或 Agent）只想知道这个技能有哪些参数、
怎么调用。它不应该因为缺少 numpy / pandas / scipy 这类科学计算栈而崩溃。

这条约束的意义在于：技能被 Agent 以「先看用法，再决定是否装依赖」的方式调用。
如果 `--help` 就崩了，Agent 还没机会告诉用户"请先 pip install numpy"。

做法
----
用 ``sys.meta_path`` 插入一个 import blocker，把 numpy / pandas / scipy /
matplotlib / seaborn 全部屏蔽掉，然后对 ``skills/*/scripts/*.py`` 全量跑
``--help``。任何脚本因此崩溃都说明它存在顶层硬依赖，需要改为可选依赖降级
（见 ``packages/biopipe_core/optional_deps.py``）。

用法：
    python scripts/ci_smoke_no_deps.py
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]

# 被屏蔽的科学计算栈
BLOCKED_PACKAGES = ("numpy", "pandas", "scipy", "matplotlib", "seaborn")

# 在子进程里执行的引导代码：先装 blocker，再运行目标脚本。
_BLOCKER_TEMPLATE = '''
import sys


class _Blocker:
    BLOCKED = {blocked!r}

    def find_module(self, name, path=None):
        return self if name.split(".")[0] in self.BLOCKED else None

    def load_module(self, name):
        raise ImportError("blocked for no-deps smoke test: " + name)


sys.meta_path.insert(0, _Blocker())
'''


def _iter_entry_scripts():
    """遍历所有技能脚本目录下的 .py 文件。"""
    skills_dir = REPO_ROOT / "skills"
    if not skills_dir.is_dir():
        return
    for skill_dir in sorted(skills_dir.iterdir()):
        scripts_dir = skill_dir / "scripts"
        if not scripts_dir.is_dir():
            continue
        for script in sorted(scripts_dir.glob("*.py")):
            yield skill_dir.name, script


def main() -> int:
    blocker_code = _BLOCKER_TEMPLATE.format(blocked=list(BLOCKED_PACKAGES))
    packages_dir = REPO_ROOT / "packages"

    env = dict(os.environ)
    # PYTHONPATH 里加上 packages，保证 biopipe_core 可被导入
    existing = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = os.pathsep.join(
        p for p in (str(packages_dir), existing) if p
    )

    failures: list[str] = []
    checked = 0

    for skill_name, script in _iter_entry_scripts():
        checked += 1
        # 用 runpy 运行目标脚本，同时注入 blocker
        bootstrap = (
            f"{blocker_code}\n"
            f"import runpy, sys\n"
            f"sys.argv = [{str(script)!r}, '--help']\n"
            f"runpy.run_path({str(script)!r}, run_name='__main__')\n"
        )
        try:
            proc = subprocess.run(
                [sys.executable, "-c", bootstrap],
                capture_output=True,
                text=True,
                timeout=120,
                cwd=str(REPO_ROOT),
                env=env,
            )
        except subprocess.TimeoutExpired:
            failures.append(f"{skill_name}/{script.name}: --help 超时（120s）")
            continue

        if proc.returncode != 0:
            tail = (proc.stderr or proc.stdout or "")[-500:]
            failures.append(f"{skill_name}/{script.name}: 退出码 {proc.returncode}\n{tail}")

    print(f"已检查 {checked} 个技能脚本（屏蔽 {'/'.join(BLOCKED_PACKAGES)}）")
    if failures:
        print("\n以下脚本在无科学计算栈环境下 --help 失败：")
        print("\n".join(failures))
        return 1

    print("所有脚本在无科学计算栈环境下均能正常响应 --help")
    return 0


if __name__ == "__main__":
    sys.exit(main())
