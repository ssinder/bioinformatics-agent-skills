#!/usr/bin/env python3
"""一键发布：创建/更新 GitHub 仓库、设置元数据、推送代码。

用法：
    python scripts/publish.py --dry-run          # 只检查，不做任何写操作
    python scripts/publish.py                    # 创建仓库并推送
    python scripts/publish.py --repo ssinder/bioinformatics-agent-skills

需要 GitHub token，优先级：
    1. 环境变量 GITHUB_TOKEN
    2. Git credential manager（git credential fill）
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import urllib.error
import urllib.request
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
API = "https://api.github.com"

DEFAULT_REPO = "ssinder/bioinformatics-agent-skills"
DESCRIPTION = (
    "Agent Skills for computational biology — RNA-seq, single-cell, protein structure, "
    "molecular dynamics, peptide design, literature mining and publication figures. "
    "Turn any AI coding agent into a computational biology assistant."
)
HOMEPAGE = "https://ssinder.github.io/bioinformatics-agent-skills/"
TOPICS = [
    "agent-skills",
    "bioinformatics",
    "computational-biology",
    "rna-seq",
    "single-cell",
    "molecular-dynamics",
    "protein-structure",
    "peptide-design",
    "literature-mining",
    "scientific-visualization",
    "claude-skills",
    "claude-code",
]


def get_token() -> str | None:
    """获取 GitHub token。"""
    import os

    token = os.environ.get("GITHUB_TOKEN")
    if token:
        return token

    try:
        proc = subprocess.run(
            ["git", "credential", "fill"],
            input="protocol=https\nhost=github.com\n\n",
            capture_output=True,
            text=True,
            timeout=30,
        )
        for line in proc.stdout.splitlines():
            if line.startswith("password="):
                return line.split("=", 1)[1].strip()
    except Exception:  # noqa: BLE001
        pass
    return None


def api_call(method: str, path: str, token: str, payload: dict | None = None) -> tuple[int, dict]:
    """调用 GitHub API，返回 (状态码, 响应体)。"""
    url = f"{API}{path}"
    data = json.dumps(payload).encode() if payload is not None else None
    req = urllib.request.Request(url, data=data, method=method)
    req.add_header("Authorization", f"Bearer {token}")
    req.add_header("Accept", "application/vnd.github+json")
    req.add_header("X-GitHub-Api-Version", "2022-11-28")
    if data:
        req.add_header("Content-Type", "application/json")

    try:
        with urllib.request.urlopen(req, timeout=60) as resp:
            body = resp.read().decode()
            return resp.status, json.loads(body) if body else {}
    except urllib.error.HTTPError as exc:
        body = exc.read().decode()
        try:
            return exc.code, json.loads(body)
        except json.JSONDecodeError:
            return exc.code, {"message": body[:500]}
    except urllib.error.URLError as exc:
        return 0, {"message": str(exc)}


def check_prereqs() -> bool:
    """检查前置条件。"""
    ok = True
    for tool in ("git",):
        if subprocess.run(["which", tool], capture_output=True).returncode != 0:
            print(f"[FAIL] 未找到 {tool}")
            ok = False
    if not (REPO_ROOT / "skills").exists():
        print("[FAIL] 缺少 skills/ 目录")
        ok = False
    if not (REPO_ROOT / "README.md").exists():
        print("[FAIL] 缺少 README.md")
        ok = False
    skill_dirs = [p for p in (REPO_ROOT / "skills").iterdir() if p.is_dir()]
    print(f"[INFO] 发现 {len(skill_dirs)} 个技能目录")
    for d in sorted(skill_dirs):
        has_md = (d / "SKILL.md").exists()
        print(f"       {'OK ' if has_md else 'MISS'} {d.name}/SKILL.md")
    return ok


def main(argv: list[str] | None = None) -> int:
    """发布主函数。"""
    parser = argparse.ArgumentParser(description="发布 BioinfoAgent Skills 到 GitHub")
    parser.add_argument("--repo", default=DEFAULT_REPO, help="目标仓库 owner/name")
    parser.add_argument("--dry-run", action="store_true", help="只检查不执行")
    parser.add_argument("--no-push", action="store_true", help="只更新元数据不推送代码")
    args = parser.parse_args(argv)

    print("=" * 70)
    print("BioinfoAgent Skills 发布工具")
    print("=" * 70)

    if not check_prereqs():
        print("\n前置条件检查未通过")
        return 1

    token = get_token()
    if not token:
        print("[FAIL] 未找到 GitHub token，请设置 GITHUB_TOKEN 或配置 git credential")
        return 1
    print(f"[OK] 已获取 token (长度 {len(token)})")

    owner, repo = args.repo.split("/", 1)

    if args.dry_run:
        print("\n[DRY RUN] 将执行以下操作:")
        print(f"  1. 创建/确认仓库 {args.repo}")
        print("  2. 设置 description、homepage、topics")
        print("  3. 推送 main 分支")
        return 0

    print(f"\n[1/4] 检查仓库 {args.repo}")
    status, data = api_call("GET", f"/repos/{owner}/{repo}", token)
    if status == 200:
        print(f"      仓库已存在（当前 {data.get('stargazers_count', 0)} stars）")
    elif status == 404:
        print("      仓库不存在，创建中...")
        status, data = api_call(
            "POST",
            "/user/repos",
            token,
            {
                "name": repo,
                "description": DESCRIPTION,
                "homepage": HOMEPAGE,
                "private": False,
                "auto_init": False,
                "has_issues": True,
                "has_wiki": True,
                "has_projects": False,
                "license_template": "mit",
            },
        )
        if status not in (200, 201):
            print(f"      [FAIL] 创建失败: {data.get('message')}")
            return 1
        print(f"      [OK] 创建成功: {data.get('html_url')}")
    else:
        print(f"      [FAIL] 查询失败: {data.get('message')}")
        return 1

    print("\n[2/4] 更新仓库元数据")
    api_call(
        "PATCH",
        f"/repos/{owner}/{repo}",
        token,
        {"description": DESCRIPTION, "homepage": HOMEPAGE, "has_wiki": True, "topics": TOPICS},
    )
    print("      [OK] description / homepage / topics 已设置")

    print("\n[3/4] 设置分支保护与 Pages")
    api_call(
        "PUT",
        f"/repos/{owner}/{repo}/topics",
        token,
        {"names": TOPICS},
    )
    print("      [OK] topics 已同步")

    if args.no_push:
        print("\n[4/4] 跳过推送 (--no-push)")
        return 0

    print("\n[4/4] 推送代码")
    remote_url = f"https://github.com/{args.repo}.git"

    def run(cmd: list[str]) -> tuple[int, str]:
        proc = subprocess.run(cmd, cwd=REPO_ROOT, capture_output=True, text=True, timeout=600)
        return proc.returncode, (proc.stdout + proc.stderr)

    rc, out = run(["git", "rev-parse", "--git-dir"])
    if rc != 0:
        run(["git", "init", "-b", "main"])

    rc, out = run(["git", "remote", "remove", "origin"])
    if rc != 0:
        pass
    rc, out = run(["git", "remote", "add", "origin", remote_url])
    if rc != 0:
        print(f"      [FAIL] 设置 remote 失败: {out[:300]}")
        return 1

    rc, out = run(["git", "add", "-A"])
    if rc != 0:
        print(f"      [FAIL] git add 失败: {out[:300]}")
        return 1

    rc, out = run(["git", "status", "--porcelain"])
    if not out.strip():
        print("      没有待提交变更")
    else:
        rc, out = run(["git", "commit", "-m", "feat: 改造为 Agent Skills 标准的科研技能库"])
        if rc != 0:
            print(f"      [FAIL] commit 失败: {out[:300]}")
            return 1
        print("      [OK] commit 完成")

    rc, out = run(["git", "push", "-u", "origin", "main"])
    if rc != 0:
        print(f"      [FAIL] push 失败: {out[:500]}")
        print("      若为认证问题,先执行: git config credential.helper store")
        return 1
    print("      [OK] 推送成功")

    print("\n" + "=" * 70)
    print(f"发布完成: https://github.com/{args.repo}")
    print("=" * 70)
    return 0


if __name__ == "__main__":
    sys.exit(main())
