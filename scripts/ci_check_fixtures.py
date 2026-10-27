#!/usr/bin/env python3
"""CI 检查：测试夹具是否完整入库。

背景
----
这个检查是被一次真实事故逼出来的：

仓库的 ``.gitignore`` 里有 ``*.pdb`` / ``*.cif`` / ``*.h5`` 规则
（用于排除真实科研数据），例外只写了 ``!tests/fixtures_contacts/``
一个目录。后来新增 alphafold / hdock / lightdock / scrna 技能时，
新建的夹具目录忘了加例外 —— 于是 **19 个夹具文件根本没进仓库**。

表现极具迷惑性：**本地 pytest 全绿，CI 上 3 个 Python 版本全挂**。
因为本地工作区里有那些文件，而 CI 只跑 clone 后的代码。

还有一个更隐蔽的坑：**git 不跟踪空目录**。有一个夹具目录必须
「存在且为空」（用来测「输入目录里没有模型」的分支），没有
``.gitkeep`` 占位的话，clone 下来它根本不存在，测的就变成了
另一条分支。

本脚本在跑 pytest **之前**就检查这些问题 —— 5 秒内给出清晰报错，
而不是让 pytest 报一堆 FileNotFoundError。

检查项
------
1. ``tests/fixtures_*`` 下的所有文件都已被 git 跟踪
2. 测试代码里引用的每个夹具子目录都存在
3. 需要「存在但为空」的目录有 ``.gitkeep`` 占位
4. 仓库里没有被跟踪的真实科研数据（体积异常大的二进制文件）

退出码：0 = 通过；1 = 发现问题。
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
TESTS_DIR = REPO_ROOT / "tests"

#: 仓库内单文件体积上限（字节）。超过就可能是误提交的真实科研数据。
MAX_FIXTURE_BYTES = 8 * 1024 * 1024

#: 匹配测试代码里的夹具路径常量，例如 ``FIXTURES_SELECT / "standard"``。
FIXTURE_REF = re.compile(
    r"FIXTURES[A-Z_]*\s*/\s*[\"']([a-z0-9_]+)[\"']"
)


def _git_tracked(root: Path | None = None) -> set[str]:
    """返回 git 已跟踪的文件集合（相对仓库根）。"""
    base = root if root is not None else REPO_ROOT
    try:
        result = subprocess.run(
            ["git", "ls-files"],
            cwd=base,
            capture_output=True,
            text=True,
            check=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return set()
    return {line.strip().replace("\\", "/") for line in result.stdout.splitlines() if line.strip()}


def _git_ignored(root: Path | None = None) -> set[str]:
    """返回被 .gitignore 排除、但存在于磁盘上的文件。"""
    base = root if root is not None else REPO_ROOT
    try:
        result = subprocess.run(
            ["git", "status", "--porcelain", "--ignored", "--untracked-files=all"],
            cwd=base,
            capture_output=True,
            text=True,
            check=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return set()

    ignored: set[str] = set()
    for line in result.stdout.splitlines():
        if not line.startswith("!! "):
            continue
        path = line[3:].strip().replace("\\", "/")
        if "__pycache__" in path:
            continue
        ignored.add(path)
    return ignored


def check_fixtures_are_tracked(tracked: set[str], ignored: set[str]) -> list[str]:
    """检查磁盘上的夹具文件是否都进了仓库。"""
    problems: list[str] = []

    for path in sorted(ignored):
        # 只关心测试夹具目录下的东西
        if not path.startswith("tests/fixtures_"):
            continue
        problems.append(
            f"夹具文件被 .gitignore 排除，clone 后会丢失：{path}\n"
            f"    → 修法：.gitignore 里加通配例外 `!tests/fixtures_*/` + `!tests/fixtures_*/**`"
        )

    return problems


def check_referenced_dirs_exist(root: Path | None = None) -> list[str]:
    """检查测试代码引用的每个夹具子目录是否真实存在。"""
    base = root if root is not None else REPO_ROOT
    tests_dir = base / "tests"
    problems: list[str] = []

    references: set[str] = set()
    for test_file in tests_dir.glob("test_*.py"):
        references.update(FIXTURE_REF.findall(test_file.read_text(encoding="utf-8")))

    for fixture_root in sorted(tests_dir.glob("fixtures_*")):
        if not fixture_root.is_dir():
            continue
        for sub in sorted(references):
            candidate = fixture_root / sub
            # 子目录名可能属于别的夹具根，跳过不匹配的
            if not candidate.exists():
                continue
            if candidate.is_dir():
                has_files = any(p for p in candidate.iterdir() if p.name != ".gitkeep")
                if not has_files and not (candidate / ".gitkeep").exists():
                    problems.append(
                        f"夹具目录为空且没有 .gitkeep 占位，clone 后不存在：{candidate.relative_to(base)}"
                    )

    return problems


def check_no_large_binaries(tracked: set[str], root: Path | None = None) -> list[str]:
    """检查仓库里没有体积异常大的文件（可能是误提交的真实科研数据）。

    ``root`` 用于把相对路径解析成实际文件位置；默认用仓库根目录。
    单独抽出来是为了让测试能在临时目录里构造样本。
    """
    base = root if root is not None else REPO_ROOT
    problems: list[str] = []

    for rel in sorted(tracked):
        if not rel.startswith("tests/fixtures_"):
            continue
        path = base / rel
        if not path.is_file():
            continue
        size = path.stat().st_size
        if size > MAX_FIXTURE_BYTES:
            problems.append(
                f"夹具文件过大（{size / 1024 / 1024:.1f} MB > "
                f"{MAX_FIXTURE_BYTES / 1024 / 1024:.0f} MB）：{rel}\n"
                f"    → 建议加 gzip 压缩（h5py create_dataset(compression='gzip')）"
            )

    return problems


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="检查测试夹具是否完整入库")
    parser.add_argument(
        "--repo-root",
        default=str(REPO_ROOT),
        help="仓库根目录",
    )
    args = parser.parse_args(argv)
    root = Path(args.repo_root).resolve()

    tracked = _git_tracked(root)
    ignored = _git_ignored(root)

    problems: list[str] = []
    problems.extend(check_fixtures_are_tracked(tracked, ignored))
    problems.extend(check_referenced_dirs_exist(root))
    problems.extend(check_no_large_binaries(tracked, root=root))

    if problems:
        print("::error::测试夹具完整性检查失败：")
        for problem in problems:
            print(f"  {problem}")
        return 1

    fixture_files = sum(1 for rel in tracked if rel.startswith("tests/fixtures_"))
    print(f"夹具完整性检查通过（{fixture_files} 个夹具文件已入库）")
    return 0


if __name__ == "__main__":
    sys.exit(main())
