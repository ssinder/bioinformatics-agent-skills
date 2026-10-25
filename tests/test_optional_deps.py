"""可选依赖降级能力的回归测试。

背景：技能脚本曾在模块顶层 ``import numpy``，导致在未安装科学计算栈的
环境中连 ``--help`` 都无法执行。修复后顶层改为 ``optional_import``，
缺失时返回 ``None``，由入口的依赖守卫给出明确安装指引。

本模块验证两件事：

1. :mod:`biopipe_core.optional_deps` 的行为正确；
2. **每个技能脚本在屏蔽 numpy / pandas / scipy / matplotlib 后，
   仍能响应 ``--help``** —— 这是 CI 中新增的冒烟步骤的本地等价物。
"""

from __future__ import annotations

import os
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "packages"))

BLOCKED_MODULES = ("numpy", "pandas", "scipy", "matplotlib", "seaborn")

_SITECUSTOMIZE = '''\
import sys


class _Blocker:
    BLOCKED = {blocked!r}

    def find_module(self, name, path=None):
        return self if name.split(".")[0] in self.BLOCKED else None

    def load_module(self, name):
        raise ImportError("blocked for test: " + name)


sys.meta_path.insert(0, _Blocker())
'''


@pytest.fixture(scope="module")
def blocked_env():
    """构造一个屏蔽科学计算栈的子进程环境（模块级 fixture，兼容 Python 3.9）。"""
    tmpdir = Path(tempfile.mkdtemp(prefix="biopipe_nodeps_"))
    (tmpdir / "sitecustomize.py").write_text(
        _SITECUSTOMIZE.format(blocked=set(BLOCKED_MODULES)), encoding="utf-8"
    )
    env = dict(os.environ)
    existing = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = str(tmpdir) + (os.pathsep + existing if existing else "")
    return env


class TestOptionalDepsHelper:
    """optional_deps 模块的基本行为。"""

    def test_optional_import_returns_module(self) -> None:
        from biopipe_core.optional_deps import optional_import

        mod = optional_import("json")
        assert mod is not None
        assert mod.dumps({"a": 1}) == '{"a": 1}'

    def test_optional_import_returns_none_for_missing(self) -> None:
        from biopipe_core.optional_deps import optional_import

        assert optional_import("definitely_not_a_real_module_xyz") is None

    def test_require_raises_with_hint(self) -> None:
        from biopipe_core.optional_deps import MissingDependencyError, require

        with pytest.raises(MissingDependencyError) as excinfo:
            require("definitely_not_a_real_module_xyz", hint="测试提示", extra="all")
        message = str(excinfo.value)
        assert "definitely_not_a_real_module_xyz" in message
        assert "测试提示" in message
        assert "pip install" in message

    def test_is_available(self) -> None:
        from biopipe_core.optional_deps import is_available

        assert is_available("json") is True
        assert is_available("definitely_not_a_real_module_xyz") is False

    def test_core_exports_optional_deps_api(self) -> None:
        import biopipe_core

        for name in ("optional_import", "require", "is_available", "MissingDependencyError"):
            assert name in biopipe_core.__all__
            assert hasattr(biopipe_core, name)


class TestSkillScriptsDegradeGracefully:
    """技能脚本在缺科学计算栈时仍须响应 --help。"""

    def _entry_scripts(self) -> list[Path]:
        return sorted(REPO_ROOT.glob("skills/*/scripts/*.py"))

    def test_has_entry_scripts(self) -> None:
        assert len(self._entry_scripts()) >= 20, "技能脚本数量异常少"

    @pytest.mark.parametrize(
        "script",
        [
            "skills/bulk-rnaseq/scripts/run_pipeline.py",
            "skills/bulk-rnaseq/scripts/build_counts_matrix.py",
            "skills/single-cell-rnaseq/scripts/run_pipeline.py",
            "skills/molecular-dynamics/scripts/run_analysis.py",
        ],
        ids=lambda p: p.split("/")[-1],
    )
    def test_help_works_without_scientific_stack(self, script: str, blocked_env) -> None:
        """关键回归：--help 是元信息查询，不应因缺 numpy/pandas/scipy 而崩溃。"""
        entry = REPO_ROOT / script
        assert entry.exists(), f"{script} 不存在"

        proc = subprocess.run(
            [sys.executable, str(entry), "--help"],
            capture_output=True,
            text=True,
            timeout=120,
            env=blocked_env,
            cwd=str(REPO_ROOT),
        )
        assert proc.returncode == 0, (
            f"{script} 在无科学计算栈下 --help 失败 (exit {proc.returncode}):\n"
            f"{proc.stderr[-800:]}"
        )
        assert "usage" in proc.stdout.lower() or "用法" in proc.stdout, (
            f"{script} 的 --help 未输出用法说明"
        )

    def test_all_entry_scripts_help_without_scientific_stack(self, blocked_env) -> None:
        """全量回归：所有技能脚本都不能有顶层硬依赖。

        背景：`secondary_structure.py` 曾因未把脚本目录加入 sys.path 而在
        独立运行时报 `No module named 'pdb_parser'`。该缺陷恰好落在
        上面那个参数化列表之外，因此本用例做全量兜底。
        """
        failures: list[str] = []
        scripts = self._entry_scripts()

        for entry in scripts:
            proc = subprocess.run(
                [sys.executable, str(entry), "--help"],
                capture_output=True,
                text=True,
                timeout=120,
                env=blocked_env,
                cwd=str(REPO_ROOT),
            )
            if proc.returncode != 0:
                rel = entry.relative_to(REPO_ROOT)
                failures.append(f"{rel}: 退出码 {proc.returncode}\n{(proc.stderr or '')[-400:]}")

        assert not failures, "以下脚本在无科学计算栈下 --help 失败：\n" + "\n".join(failures)
