"""技能脚本 ``sys.path`` 引导的回归测试。

背景
----
技能脚本有两种常见启动方式：

1. ``python skills/<skill>/scripts/foo.py --help`` —— 此时 Python 会把
   **脚本所在目录**自动插入 ``sys.path[0]``，所以 ``from sibling import x``
   这类同目录导入可以正常工作。
2. ``runpy.run_path(...)`` —— 例如 ``scripts/ci_smoke_no_deps.py`` 就是这么
   跑的。**此时不会**把脚本目录加入 ``sys.path``，同目录导入会直接抛
   ``ModuleNotFoundError: No module named 'sibling'``。

第 2 种方式在本仓库不是假想：``ci_smoke_no_deps.py`` 是 CI 的一环，
它要确保「Agent 先看用法、再决定是否装依赖」这条链路是通的。若某个脚本
的 ``--help`` 在这种启动方式下崩溃，Agent 根本没机会告诉用户「请先
``pip install numpy``」。

约定
----
**任何导入同级模块的脚本，都必须自己把脚本目录加入 ``sys.path``**，
不能依赖解释器的隐式行为。

仓库里此前已有同类事故的记录，见
``tests/test_optional_deps.py::test_all_entry_scripts_help_without_scientific_stack``
的注释：``secondary_structure.py`` 曾因未把脚本目录加入 ``sys.path``
而在独立运行时报 ``No module named 'pdb_parser'``。

本模块锁住这条约定，并用「反向验证」确保用例真的能抓到问题：
把 ``sys.path`` 引导删掉，用例必须失败。

本模块只依赖标准库 + pytest，不依赖 numpy/pandas。
"""

from __future__ import annotations

import ast
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
SKILLS_DIR = REPO_ROOT / "skills"

# 与 scripts/ci_smoke_no_deps.py 保持一致的被屏蔽包
BLOCKED_PACKAGES = ("numpy", "pandas", "scipy", "matplotlib", "seaborn")

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


def _iter_scripts() -> list[Path]:
    """遍历所有技能脚本（排除 __init__ 与下划线开头的模块）。"""
    out: list[Path] = []
    for skill_dir in sorted(SKILLS_DIR.iterdir()):
        scripts_dir = skill_dir / "scripts"
        if not scripts_dir.is_dir():
            continue
        for py in sorted(scripts_dir.glob("*.py")):
            if py.name.startswith("_"):
                continue
            out.append(py)
    return out


def _sibling_imports(script: Path) -> set[str]:
    """返回该脚本导入的**同级**模块名集合。

    同级 = 与脚本同目录、且不是脚本自身。
    """
    scripts_dir = script.parent
    local_modules = {p.stem for p in scripts_dir.glob("*.py")}
    tree = ast.parse(script.read_text(encoding="utf-8"))

    siblings: set[str] = set()
    for node in ast.walk(tree):
        names: list[str] = []
        if isinstance(node, ast.ImportFrom) and node.module and node.level == 0:
            names.append(node.module.split(".")[0])
        elif isinstance(node, ast.Import):
            names.extend(alias.name.split(".")[0] for alias in node.names)
        for name in names:
            if name in local_modules and name != script.stem:
                siblings.add(name)
    return siblings


def _bootstrap_lines(script: Path) -> list[tuple[int, str]]:
    """找出脚本中所有把「自身目录」加入 sys.path 的语句。

    返回 (行号, 源码行) 列表。

    匹配要点：必须用 ``.parent`` 而**不是** ``.parents[N]``。
    ``parents[3] / "packages"`` 指向的是仓库的 packages 目录，
    与脚本自身目录无关，不能算作「引导自身目录」。
    """
    lines = script.read_text(encoding="utf-8").splitlines()
    hits: list[tuple[int, str]] = []
    for idx, line in enumerate(lines, start=1):
        if "sys.path.insert" not in line:
            continue
        # 排除 parents[N] 形式（那是仓库级 packages 目录）
        if re.search(r"\.parents\s*\[", line):
            continue
        if ".parent" in line and "__file__" in line:
            hits.append((idx, line))
    return hits


def _inserts_own_dir(script: Path) -> bool:
    """脚本是否把自身所在目录加入了 sys.path。"""
    return bool(_bootstrap_lines(script))


class TestSysPathBootstrapConvention:
    """任何导入同级模块的脚本都必须自己引导 sys.path。"""

    @pytest.mark.parametrize("script", _iter_scripts(), ids=lambda p: p.name)
    def test_sibling_import_scripts_bootstrap_own_dir(self, script: Path) -> None:
        """导入同级模块的脚本必须把脚本目录加入 sys.path。

        背景：``python foo.py`` 会自动把脚本目录放进 ``sys.path[0]``，
        但 ``runpy.run_path("foo.py")`` 不会。后者是 CI 冒烟测试与
        部分 Agent 工具链的实际启动方式。

        没有这行引导，脚本在第二种方式下会抛
        ``ModuleNotFoundError: No module named '<sibling>'``，
        连 ``--help`` 都打不出来。
        """
        siblings = _sibling_imports(script)
        if not siblings:
            pytest.skip(f"{script.name} 不导入同级模块，无需引导")

        assert _inserts_own_dir(script), (
            f"{script.relative_to(REPO_ROOT)} 导入了同级模块 "
            f"{sorted(siblings)}，但没有把脚本目录加入 sys.path。"
            f"在 runpy.run_path 方式下会抛 ModuleNotFoundError。"
            f"应补上：sys.path.insert(0, str(Path(__file__).resolve().parent))"
        )

    @pytest.mark.parametrize("script", _iter_scripts(), ids=lambda p: p.name)
    def test_bootstrap_precedes_sibling_import(self, script: Path) -> None:
        """sys.path 引导必须发生在同级模块导入**之前**。

        背景：``geometry.py`` 与 ``interaction_network.py`` 曾把
        ``sys.path.insert(...)`` 放在 ``if __name__ == "__main__":`` 块里，
        而同级模块导入放在 ``_main()`` 函数体内。执行顺序上是
        ``__main__`` 块先跑（设好 sys.path）再进 ``_main()``（导入），
        恰好能工作——但这依赖「导入被推迟到函数调用时」这一隐式前提，
        一旦有人把导入提到模块顶层就会立刻炸。

        这里锁住的更强约束是：**模块级导入语句执行时，sys.path 必须已就绪**。
        """
        siblings = _sibling_imports(script)
        if not siblings:
            pytest.skip(f"{script.name} 不导入同级模块，无需引导")

        bootstrap_lines = _bootstrap_lines(script)
        if not bootstrap_lines:
            pytest.skip(f"{script.name} 已由其他用例覆盖（缺少引导）")

        # 找同级模块的**模块级** import 行号
        tree = ast.parse(script.read_text(encoding="utf-8"))
        module_level_imports: list[int] = []
        for node in tree.body:  # 只看模块顶层，不进函数体
            if isinstance(node, ast.ImportFrom) and node.module and node.level == 0:
                if node.module.split(".")[0] in siblings:
                    module_level_imports.append(node.lineno)
            elif isinstance(node, ast.Import):
                for alias in node.names:
                    if alias.name.split(".")[0] in siblings:
                        module_level_imports.append(node.lineno)

        if not module_level_imports:
            # 导入全在函数体内（惰性导入），__main__ 块的引导够用
            pytest.skip(f"{script.name} 的同级导入在函数体内（惰性导入），顺序无风险")

        first_bootstrap = min(line for line, _ in bootstrap_lines)
        first_import = min(module_level_imports)
        assert first_bootstrap < first_import, (
            f"{script.relative_to(REPO_ROOT)}: sys.path 引导在第 {first_bootstrap} 行，"
            f"但同级模块导入在第 {first_import} 行 —— 引导必须在前。"
        )


def _has_cli(script: Path) -> bool:
    """脚本是否定义了命令行入口（argparse + __main__）。"""
    src = script.read_text(encoding="utf-8")
    return "argparse" in src and "__main__" in src


class TestRunpyExecution:
    """用 runpy 真实执行脚本，确认 ``--help`` 不因缺依赖而崩溃。"""

    @pytest.mark.parametrize(
        "script",
        [p for p in _iter_scripts() if _sibling_imports(p)],
        ids=lambda p: p.name,
    )
    def test_help_works_under_runpy(self, script: Path) -> None:
        """``runpy.run_path`` 方式下 ``--help`` 必须成功。

        这与 ``scripts/ci_smoke_no_deps.py`` 的执行方式一致，
        但这里只跑「导入同级模块」的脚本，聚焦路径引导这一类失败。

        注意：部分脚本（如 ``peptide_filters.py``、``scalgorithms.py``）
        是**纯库模块**，没有 CLI。它们能正常 import 即视为通过；
        对有 CLI 的脚本，额外要求输出 usage。
        """
        env = dict(os.environ)
        env["PYTHONPATH"] = os.pathsep.join(
            p for p in (str(REPO_ROOT / "packages"), env.get("PYTHONPATH", "")) if p
        )
        bootstrap = (
            _BLOCKER_TEMPLATE.format(blocked=list(BLOCKED_PACKAGES))
            + "\nimport runpy, sys\n"
            + f"sys.argv = [{str(script)!r}, '--help']\n"
            + f"runpy.run_path({str(script)!r}, run_name='__main__')\n"
        )
        proc = subprocess.run(
            [sys.executable, "-c", bootstrap],
            capture_output=True,
            text=True,
            timeout=120,
            cwd=str(REPO_ROOT),
            env=env,
        )
        assert proc.returncode == 0, (
            f"{script.relative_to(REPO_ROOT)} 在 runpy 方式下 --help 失败：\n"
            f"{(proc.stderr or proc.stdout or '')[-500:]}"
        )
        if _has_cli(script):
            assert "usage" in (proc.stdout or "").lower(), (
                f"{script.relative_to(REPO_ROOT)} 的 --help 未输出用法说明"
            )


class TestRegressionCanary:
    """反向验证：这套用例必须真的能抓到问题。

    如果有人把某个脚本的 ``sys.path`` 引导删掉，对应用例必须失败。
    这里用一个**临时构造的最小脚本**来验证检测逻辑本身有效，
    而不去篡改仓库里的真实文件。
    """

    def test_detector_flags_missing_bootstrap(self, tmp_path: Path) -> None:
        """构造一个缺引导的脚本，检测逻辑必须报警。"""
        scripts_dir = tmp_path / "scripts"
        scripts_dir.mkdir()
        sibling = scripts_dir / "sibling_mod.py"
        sibling.write_text("VALUE = 1\n", encoding="utf-8")
        offender = scripts_dir / "offender.py"
        offender.write_text(
            "from __future__ import annotations\n\n"
            "import sys\n"
            "from pathlib import Path\n\n"
            'sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "packages"))\n\n'
            "from sibling_mod import VALUE  # noqa: E402\n\n"
            "print(VALUE)\n",
            encoding="utf-8",
        )

        assert _sibling_imports(offender) == {"sibling_mod"}, "检测逻辑失效：没识别出同级导入"
        assert not _inserts_own_dir(offender), "检测逻辑失效：误报了已有引导"

    def test_detector_accepts_correct_bootstrap(self, tmp_path: Path) -> None:
        """正确引导的脚本必须通过检测。"""
        scripts_dir = tmp_path / "scripts"
        scripts_dir.mkdir()
        (scripts_dir / "sibling_mod.py").write_text("VALUE = 1\n", encoding="utf-8")
        good = scripts_dir / "good.py"
        good.write_text(
            "from __future__ import annotations\n\n"
            "import sys\n"
            "from pathlib import Path\n\n"
            'sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "packages"))\n'
            "sys.path.insert(0, str(Path(__file__).resolve().parent))\n\n"
            "from sibling_mod import VALUE  # noqa: E402\n\n"
            "print(VALUE)\n",
            encoding="utf-8",
        )

        assert _sibling_imports(good) == {"sibling_mod"}
        assert _inserts_own_dir(good), "检测逻辑失效：没认出正确的引导写法"
