#!/usr/bin/env python3
"""跨平台依赖漏洞审计，供 CI 调用。

为什么需要这个脚本？
------------------
`pip-audit` 默认审计**整个 Python 环境**。在 CI 中我们用 `pip install -e .`
把项目自身装成 editable 发行版，而项目并不发布到 PyPI，于是 pip-audit 会报：

    ERROR: Dependency not found on PyPI and could not be audited: biopipe-agent-skills

配合 `--strict`（任何"无法审计"的依赖都视为致命）时，这会让流水线直接失败，
即使第三方依赖本身完全没有漏洞。同理，任何第三方库以 editable 方式安装
（常见于 CI 缓存、monorepo 本地依赖）也会踩到同一个坑。

解决思路：**审计依赖清单，而不是审计环境本身**。

  1. 冻结实际安装的第三方发行版（`pip freeze --exclude-editable`），
     覆盖包括传递依赖在内的完整依赖树；
  2. 额外审计 `pyproject.toml` 中声明的核心依赖，确保"声明的"和"装上的"都干净。

这样既保持严格（任何真实漏洞都会导致非零退出码），又不会因为项目自身
无法在 PyPI 上检索而误报。

用法：
    python scripts/ci_dependency_audit.py
"""

from __future__ import annotations

import re
import subprocess
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]

# pyproject.toml 中真正需要审计的核心依赖。
# 可选依赖组（dev/science 等）单独审计，避免把 CI 工具链混进核心结论里。
CORE_REQUIREMENTS = (
    "requests>=2.28",
    "pyyaml>=6.0",
)

# 解析 pyproject.toml 里 [project].dependencies 的轻量正则。
# 刻意不引入 tomllib：Python 3.9 没有 tomllib，3.11+ 才有。
_DEP_LINE_RE = re.compile(r"^\s*(?P<dep>[A-Za-z0-9_.\-\[\]]+\s*[><=!~]+\s*[^\s,#]+)\s*,?\s*$")


def _run(cmd: list[str]) -> tuple[int, str, str]:
    """执行命令并返回 (returncode, stdout, stderr)。"""
    proc = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8", errors="replace")
    return proc.returncode, proc.stdout or "", proc.stderr or ""


def _parse_pyproject_dependencies(pyproject: Path) -> list[str]:
    """从 pyproject.toml 的 [project].dependencies 解析核心依赖列表。

    用正则而非 tomllib，是为了兼容 Python 3.9（3.11 才引入 tomllib）。
    """
    if not pyproject.is_file():
        return []

    text = pyproject.read_text(encoding="utf-8")
    # 定位 [project] 段内的 dependencies = [...] 列表
    match = re.search(r"^dependencies\s*=\s*\[(.*?)\]", text, re.DOTALL | re.MULTILINE)
    if not match:
        return []

    deps: list[str] = []
    for line in match.group(1).splitlines():
        line = line.split("#")[0].strip()
        if not line:
            continue
        dep_match = _DEP_LINE_RE.match(line)
        if dep_match:
            deps.append(dep_match.group("dep").replace(" ", ""))
    return deps


def _audit_requirements_file(reqs: list[str], label: str) -> tuple[bool, str]:
    """用 pip-audit 审计给定的依赖清单。

    返回 (是否通过, 完整输出)。
    """
    if not reqs:
        return True, f"{label}: 无依赖需要审计"

    with tempfile.NamedTemporaryFile("w", suffix=".txt", encoding="utf-8", delete=False) as fh:
        fh.write("\n".join(reqs) + "\n")
        req_path = fh.name

    try:
        code, out, err = _run([
            sys.executable, "-m", "pip_audit",
            "--progress-spinner", "off",
            "--strict",
            "-r", req_path,
        ])
    finally:
        Path(req_path).unlink(missing_ok=True)

    combined = (out + err).strip()
    if code == 0:
        return True, f"{label}: OK\n{combined}"
    return False, f"{label}: FAILED (exit={code})\n{combined}"


def _audit_environment() -> tuple[bool, str]:
    """审计环境中实际安装的第三方发行版（排除 editable，即项目自身）。"""
    code, out, err = _run([sys.executable, "-m", "pip", "freeze", "--exclude-editable"])
    if code != 0:
        return False, f"pip freeze 失败: {err.strip()}"

    # 过滤掉 pip 自身的工具链行（如 pip==26.2.2）——pip 是 CI 工具链而非项目依赖，
    # 且其版本受 runner 镜像控制，无法由本项目声明。
    frozen = [
        line.strip()
        for line in out.splitlines()
        if line.strip() and not line.strip().lower().startswith("pip==")
    ]
    return _audit_requirements_file(frozen, "已安装依赖（pip freeze --exclude-editable）")


def main() -> int:
    if _run([sys.executable, "-m", "pip_audit", "--version"])[0] != 0:
        print("::error::pip-audit 未安装，请先运行: pip install pip-audit")
        return 1

    results: list[tuple[bool, str]] = []

    # 1) 核心依赖：从 pyproject.toml 解析
    core = _parse_pyproject_dependencies(REPO_ROOT / "pyproject.toml") or list(CORE_REQUIREMENTS)
    results.append(_audit_requirements_file(core, "核心依赖（pyproject.toml 声明）"))

    # 2) 环境中实际安装的完整依赖树
    results.append(_audit_environment())

    print("=" * 70)
    ok = True
    for passed, message in results:
        status = "PASS" if passed else "FAIL"
        print(f"[{status}] {message}")
        print("-" * 70)
        ok = ok and passed

    if not ok:
        print("::error::依赖审计未通过，请升级受影响依赖后重试")
        return 1

    print("::notice::依赖审计通过，未发现已知漏洞")
    return 0


if __name__ == "__main__":
    sys.exit(main())
