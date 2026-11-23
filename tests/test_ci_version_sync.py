"""版本同步检查器的测试。

版本号分散在 ``pyproject.toml``、``plugin.json``、``CITATION.cff`` 三处，
都是手工维护的。一旦只改了其中一处，发布出去的 wheel、插件清单和
引用元数据就会互相矛盾——而这类不一致很难被察觉：用户 ``pip install``
到的版本与仓库页面上显示的版本对不上。

这里既检查现状一致，也验证检查器本身能抓住漂移。
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "scripts" / "ci_version_sync.py"

sys.path.insert(0, str(REPO_ROOT / "scripts"))

# tomllib 是 Python 3.11+ 才有的标准库；3.9 / 3.10 用 tomli。
# 本仓库声明支持 3.9，因此测试也必须能在 3.9 上跑。
try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover
    import tomli as tomllib


def _run(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(SCRIPT), *args],
        capture_output=True,
        text=True,
        timeout=120,
    )


class TestCurrentState:
    """当前仓库状态。"""

    def test_versions_are_consistent(self) -> None:
        """三处版本必须一致。"""
        proc = _run()
        assert proc.returncode == 0, f"版本不一致:\n{proc.stdout}\n{proc.stderr}"

    def test_reports_the_version(self) -> None:
        """输出应包含当前版本号。"""
        proc = _run()
        assert "1.0.0" in proc.stdout

    def test_all_skills_have_valid_versions(self) -> None:
        """每个技能的 metadata.version 都应合法。"""
        proc = _run()
        assert "17 个技能版本格式合法" in proc.stdout


class TestConsistencyAcrossFiles:
    """三处版本确实指向同一个值。"""

    def test_pyproject_and_plugin_match(self) -> None:
        """pyproject.toml 与 plugin.json 版本一致。"""
        pyproject = tomllib.loads(
            (REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8")
        )
        plugin = json.loads((REPO_ROOT / "plugin.json").read_text(encoding="utf-8"))
        assert pyproject["project"]["version"] == plugin["version"]

    def test_citation_matches(self) -> None:
        """CITATION.cff 版本与包版本一致。"""
        import re

        pyproject = tomllib.loads(
            (REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8")
        )
        citation = (REPO_ROOT / "CITATION.cff").read_text(encoding="utf-8")
        m = re.search(r'^version:\s*"?([^"\n]+)"?\s*$', citation, re.MULTILINE)
        assert m, "CITATION.cff 中找不到 version 字段"
        assert m.group(1).strip() == pyproject["project"]["version"]

    def test_release_date_is_not_in_the_future(self) -> None:
        """CITATION.cff 的发布日期不得晚于今天。"""
        import datetime as dt
        import re

        citation = (REPO_ROOT / "CITATION.cff").read_text(encoding="utf-8")
        m = re.search(r'^date-released:\s*"?([\d\-]+)"?', citation, re.MULTILINE)
        assert m, "CITATION.cff 中找不到 date-released"
        released = dt.date.fromisoformat(m.group(1))
        assert released <= dt.date.today(), f"发布日期 {released} 晚于今天"


class TestCheckerCatchesDrift:
    """检查器本身必须能抓住漂移。

    用临时目录模拟，而不是改动真实文件。
    """

    def _run_against(self, tmp_root: Path, mutate) -> subprocess.CompletedProcess:
        """复制最小仓库结构到临时目录，跑检查器。"""
        import shutil

        (tmp_root / "scripts").mkdir(parents=True, exist_ok=True)
        shutil.copy(SCRIPT, tmp_root / "scripts" / SCRIPT.name)
        shutil.copy(REPO_ROOT / "pyproject.toml", tmp_root / "pyproject.toml")
        shutil.copy(REPO_ROOT / "plugin.json", tmp_root / "plugin.json")
        shutil.copy(REPO_ROOT / "CITATION.cff", tmp_root / "CITATION.cff")
        shutil.copytree(REPO_ROOT / "skills", tmp_root / "skills")
        mutate(tmp_root)
        return subprocess.run(
            [sys.executable, str(tmp_root / "scripts" / SCRIPT.name)],
            capture_output=True,
            text=True,
            timeout=180,
            cwd=str(tmp_root),
        )

    def test_detects_plugin_json_drift(self, tmp_path: Path) -> None:
        """plugin.json 版本被改动 → 必须报错。"""

        def mutate(root: Path) -> None:
            path = root / "plugin.json"
            data = json.loads(path.read_text(encoding="utf-8"))
            data["version"] = "0.9.0"
            path.write_text(json.dumps(data, indent=2), encoding="utf-8")

        proc = self._run_against(tmp_path, mutate)
        assert proc.returncode == 1, "版本漂移必须被拦截"
        assert "不一致" in proc.stdout

    def test_detects_pyproject_drift(self, tmp_path: Path) -> None:
        """pyproject.toml 版本被改动 → 必须报错。"""

        def mutate(root: Path) -> None:
            path = root / "pyproject.toml"
            text = path.read_text(encoding="utf-8")
            text = text.replace(
                'version = "1.0.0"', 'version = "0.8.1"', 1
            )
            path.write_text(text, encoding="utf-8")

        proc = self._run_against(tmp_path, mutate)
        assert proc.returncode == 1, "版本漂移必须被拦截"

    def test_detects_invalid_skill_version(self, tmp_path: Path) -> None:
        """技能的 metadata.version 非法 → 必须报错。"""

        def mutate(root: Path) -> None:
            path = root / "skills" / "bulk-rnaseq" / "SKILL.md"
            text = path.read_text(encoding="utf-8")
            text = text.replace('version: "1.0"', 'version: "not-a-version"', 1)
            path.write_text(text, encoding="utf-8")

        proc = self._run_against(tmp_path, mutate)
        assert proc.returncode == 1, "非法技能版本必须被拦截"
        assert "不是合法版本号" in proc.stdout

    def test_set_mode_normalizes_versions(self, tmp_path: Path) -> None:
        """--set 应把三处统一改写为指定版本。"""

        def mutate(root: Path) -> None:
            path = root / "plugin.json"
            data = json.loads(path.read_text(encoding="utf-8"))
            data["version"] = "0.9.0"
            path.write_text(json.dumps(data, indent=2), encoding="utf-8")

        self._run_against(tmp_path, mutate)

        proc = subprocess.run(
            [sys.executable, str(tmp_path / "scripts" / SCRIPT.name), "--set", "2.0.0"],
            capture_output=True,
            text=True,
            timeout=180,
            cwd=str(tmp_path),
        )
        assert proc.returncode == 0, proc.stdout

        plugin = json.loads((tmp_path / "plugin.json").read_text(encoding="utf-8"))
        assert plugin["version"] == "2.0.0"

        import re

        pyproject = tomllib.loads(
            (tmp_path / "pyproject.toml").read_text(encoding="utf-8")
        )
        assert pyproject["project"]["version"] == "2.0.0"

        citation = (tmp_path / "CITATION.cff").read_text(encoding="utf-8")
        m = re.search(r'^version:\s*"?([^"\n]+)"?', citation, re.MULTILINE)
        assert m.group(1).strip() == "2.0.0"
