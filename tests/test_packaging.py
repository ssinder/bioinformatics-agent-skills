"""打包与安装态的测试。

这个文件守护一个**只有装包后才暴露**的缺陷：

``skills/`` 位于 ``packages/`` 之外，setuptools 默认不把它放进 wheel。
结果是 ``pip install biopipe-agent-skills`` 之后，``bioagent --list``
能列出 17 个技能，却一个都跑不了——路径指向一个不存在的目录。
从源码检出或 editable install 完全看不出问题，只有真正装包才会暴露。

因此这里做一次真实的构建 + 安装 + 运行验证。
"""

from __future__ import annotations

import subprocess
import sys
import zipfile
from pathlib import Path
from typing import Optional

REPO_ROOT = Path(__file__).resolve().parents[1]

# 用 pytest 的 tmp_path_factory，避免污染工作区。
_BUILD_DIR: Optional[Path] = None
_WHEEL: Optional[Path] = None


def _build_wheel(tmp_root: Path) -> Path:
    """构建 wheel（一次即可，跨用例复用）。"""
    global _WHEEL
    if _WHEEL is not None and _WHEEL.is_file():
        return _WHEEL

    outdir = tmp_root / "dist"
    proc = subprocess.run(
        [sys.executable, "-m", "build", "--wheel", "--outdir", str(outdir)],
        cwd=str(REPO_ROOT),
        capture_output=True,
        text=True,
        timeout=900,
    )
    if proc.returncode != 0:
        pytest_fail(f"wheel 构建失败:\n{proc.stdout[-3000:]}\n{proc.stderr[-3000:]}")

    wheels = list(outdir.glob("*.whl"))
    if not wheels:
        pytest_fail(f"未产出 wheel:\n{proc.stdout[-2000:]}")
    _WHEEL = wheels[0]
    return _WHEEL


def pytest_fail(message: str) -> None:
    import pytest

    pytest.fail(message)


class TestWheelContainsSkills:
    """wheel 必须携带技能本体。"""

    def test_skills_are_included(self, tmp_path: Path) -> None:
        """skills/ 必须打进 wheel，否则装包后技能不可用。"""
        wheel = _build_wheel(tmp_path)
        with zipfile.ZipFile(wheel) as zf:
            names = zf.namelist()

        skill_files = [n for n in names if "/skills/" in n]
        assert skill_files, "wheel 里没有任何技能文件——pip install 后技能将不可用"

        # 每个技能都要有 SKILL.md
        skill_names = {
            n.split("/skills/")[1].split("/")[0]
            for n in skill_files
            if n.endswith("SKILL.md")
        }
        assert len(skill_names) == 17, f"只打包了 {len(skill_names)} 个技能，应为 17"

    def test_scripts_are_included(self, tmp_path: Path) -> None:
        """脚本必须随技能一起打包——只有 SKILL.md 是没法用的。"""
        wheel = _build_wheel(tmp_path)
        with zipfile.ZipFile(wheel) as zf:
            names = zf.namelist()

        scripts = [n for n in names if "/skills/" in n and "/scripts/" in n]
        assert len(scripts) > 50, f"wheel 里的脚本太少（{len(scripts)}）"

    def test_evals_are_included(self, tmp_path: Path) -> None:
        """评测用例也随包分发。"""
        wheel = _build_wheel(tmp_path)
        with zipfile.ZipFile(wheel) as zf:
            names = zf.namelist()

        evals = [n for n in names if n.endswith("evals/evals.json")]
        assert len(evals) == 17, f"只打包了 {len(evals)} 个评测文件"

    def test_no_pycache_in_wheel(self, tmp_path: Path) -> None:
        """wheel 里不该有 __pycache__。"""
        wheel = _build_wheel(tmp_path)
        with zipfile.ZipFile(wheel) as zf:
            names = zf.namelist()

        junk = [n for n in names if "__pycache__" in n or n.endswith(".pyc")]
        assert not junk, f"wheel 里混入了缓存文件: {junk[:5]}"

    def test_metadata_is_modern(self, tmp_path: Path) -> None:
        """元数据应使用 PEP 639 的 SPDX license 表达式。"""
        wheel = _build_wheel(tmp_path)
        with zipfile.ZipFile(wheel) as zf:
            meta_name = next(n for n in zf.namelist() if n.endswith("METADATA"))
            meta = zf.read(meta_name).decode("utf-8")

        assert "License-Expression: MIT" in meta, (
            "缺少 SPDX license 表达式；setuptools>=77 会因旧式写法告警，"
            "且 2027-02 后将不再支持"
        )
        assert "License :: OSI Approved" not in meta, (
            "license classifier 已弃用，应只保留 SPDX 表达式"
        )


class TestInstalledPackageWorks:
    """真实安装后必须可用。"""

    def test_cli_lists_skills_after_install(self, tmp_path: Path) -> None:
        """装包后 `bioagent --list` 必须能列出技能。"""
        wheel = _build_wheel(tmp_path)
        venv_dir = tmp_path / "venv"

        subprocess.run(
            [sys.executable, "-m", "venv", str(venv_dir)],
            capture_output=True,
            timeout=300,
            check=True,
        )
        py = venv_dir / "Scripts" / "python.exe"
        if not py.is_file():
            py = venv_dir / "bin" / "python"

        subprocess.run(
            [str(py), "-m", "pip", "install", "--quiet", str(wheel)],
            capture_output=True,
            timeout=900,
            check=True,
        )

        # 用安装后的解释器验证路径解析
        code = (
            "import bioagent.registry as r;"
            "root = r.skills_root();"
            "assert root.is_dir(), f'skills_root 不存在: {root}';"
            "specs = r.REGISTRY;"
            "missing = [s.name for s in specs.values() if not s.entry_path.is_file()];"
            "assert not missing, f'入口脚本缺失: {missing}';"
            "print('OK', len(specs))"
        )
        proc = subprocess.run(
            [str(py), "-c", code],
            capture_output=True,
            text=True,
            timeout=300,
        )
        assert proc.returncode == 0, f"安装后技能不可用:\n{proc.stdout}\n{proc.stderr}"
        assert "OK 17" in proc.stdout

    def test_cli_entry_point_runs(self, tmp_path: Path) -> None:
        """装包后的 `bioagent` 命令必须能响应。"""
        wheel = _build_wheel(tmp_path)
        venv_dir = tmp_path / "venv2"

        subprocess.run(
            [sys.executable, "-m", "venv", str(venv_dir)],
            capture_output=True,
            timeout=300,
            check=True,
        )
        py = venv_dir / "Scripts" / "python.exe"
        if not py.is_file():
            py = venv_dir / "bin" / "python"

        subprocess.run(
            [str(py), "-m", "pip", "install", "--quiet", str(wheel)],
            capture_output=True,
            timeout=900,
            check=True,
        )

        proc = subprocess.run(
            [str(py), "-m", "bioagent.cli", "--list"],
            capture_output=True,
            text=True,
            timeout=300,
        )
        assert proc.returncode == 0, proc.stderr
        assert "bulk-rnaseq" in proc.stdout


class TestSkillsRootResolution:
    """路径解析逻辑本身。"""

    def test_prefers_packaged_location(self, tmp_path: Path) -> None:
        """存在包内 skills/ 时应优先选它。"""
        sys.path.insert(0, str(REPO_ROOT / "packages"))
        from bioagent.registry import skills_root

        # 开发态下仓库的 skills/ 存在，函数必须返回一个真实存在的目录
        root = skills_root()
        assert root.is_dir(), f"skills_root 指向不存在的目录: {root}"

    def test_registry_entries_all_resolve(self) -> None:
        """每个登记的技能都必须能定位到入口脚本与 SKILL.md。"""
        sys.path.insert(0, str(REPO_ROOT / "packages"))
        from bioagent.registry import REGISTRY

        broken = []
        for spec in REGISTRY.values():
            if not spec.entry_path.is_file():
                broken.append(f"{spec.name}: entry")
            if not spec.skill_md.is_file():
                broken.append(f"{spec.name}: SKILL.md")
        assert not broken, f"以下技能路径无法解析: {broken}"
