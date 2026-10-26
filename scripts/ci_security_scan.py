#!/usr/bin/env python3
"""CI 安全扫描：检测源码中硬编码的凭据。

以 Python 实现而非 shell grep，原因有三：
1. 避免 shell 层的引号转义地狱（正则里同时含单引号和双引号）；
2. 跨平台行为一致（Windows / Linux 输出完全相同）；
3. 规则可单元测试。

检测的凭据模式覆盖：api_key、apikey、secret、secret_key、password、passwd、
pwd、access_token、auth_token、private_key 等常见命名，赋值为长度 >= 8 的
字面量字符串时告警。

豁免规则（不视为硬编码）：
- 取值来自 os.environ / os.getenv / config 对象；
- 取值是明显的占位符（如 "changeme"、"your-api-key-here"、"xxx"）；
- 取值长度 < 8（过于简单的常量，误报率高）。

退出码：0 = 通过；1 = 发现疑似硬编码凭据。
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]

# 形如  API_KEY = "....."  /  password = '.....'
#
# 关于前缀的处理（重要）
# ----------------------
# 这里刻意**不要求**变量名从行首开始，而是允许任意前缀
# （``APP_SECRET``、``FEISHU_APP_SECRET``、``_MY_TOKEN`` 等）。
#
# 早期实现用的是 ``\b(secret)\b``，而 ``\b`` 要求两侧是词边界 ——
# 但下划线与字母**同属单词字符**，所以 ``APP_SECRET`` 里 ``P`` 与 ``S``
# 之间不存在词边界，整个模式匹配不上。这导致老仓库里泄露的
# ``APP_SECRET = "prOd..."`` **完全逃过了扫描** ——
# 而那正是本次改造中发现的真实安全问题。
#
# 现在改用 ``(?<![A-Za-z0-9])`` 作为前缀断言：只阻止前面紧邻字母/数字，
# 但允许下划线与大写字母组合出的常见命名风格。
#
# 另外三点：
# 1. 变量名用 ``(?P<name>...)`` **整体捕获**，报告里给出的是
#    ``FEISHU_APP_SECRET`` 这样的完整名字，而不是最后一个片段
#    ``app_secret`` —— 后者对定位毫无帮助。
# 2. 关键词后带负向前瞻 ``(?![A-Za-z0-9_])``：变量名必须**以**关键词
#    结尾。``TOKENIZER`` / ``SECRETARY`` 这类只是碰巧含有关键词的
#    普通变量名不算命中，避免误报。
# 3. 赋值处允许**类型注解**：``APP_SECRET: str = "..."``。
#    老实现的 ``\s*[:=]\s*`` 会先撞上注解里的冒号，导致整行漏报。
CREDENTIAL_ASSIGNMENT = re.compile(
    r"""(?ix)
    (?<![A-Za-z0-9])
    (?P<name>
        (?: [A-Za-z_][A-Za-z0-9_]*? )?
        (?:
            api[_-]?key | apikey | secret[_-]?key | secret |
            client[_-]?secret | password | passwd | pwd |
            access[_-]?token | auth[_-]?token | api[_-]?token |
            refresh[_-]?token | private[_-]?key |
            app[_-]?secret | app[_-]?id | bot[_-]?token | token
        )
        (?![A-Za-z0-9_])
    )
    \s*
    (?:
        :\s*[A-Za-z_][\w\[\], .|]*\s*=\s*     # 带类型注解的赋值
      | [:=]\s*                                # 普通赋值 / 字典键
    )
    (?P<quote>["'])
    (?P<value>[^"'\n]{8,})
    (?P=quote)
    """
)

# 明显是占位符而非真实凭据
PLACEHOLDER_VALUES = {
    "changeme",
    "change_me",
    "your-api-key",
    "your_api_key",
    "your-api-key-here",
    "your-secret-here",
    "your-token-here",
    "xxxxxxxx",
    "xxxxxxxxxxxx",
    "todo",
    "none",
    "null",
    "example",
    "placeholder",
    "dummy",
    "fake",
    "test-token",
    "insert-key-here",
    "replace-me",
}

# 占位符子串（子串匹配，覆盖 xxx-your-token-here 之类的组合写法）
PLACEHOLDER_TOKENS = (
    "placeholder",
    "your-",
    "your_",
    "changeme",
    "change_me",
    "example",
    "dummy",
    "fake",
    "example_key",
    "todo",
    "replace",
    "insert",
    "xxxx",
    "abcdef123456",
)

# 明显不是凭据的上下文（这些取值来自环境或配置对象）
SAFE_CONTEXTS = (
    "os.environ",
    "os.getenv",
    "environ.get",
    "config.get",
    "settings.",
    "self.",
    "os.path",
)


def _is_placeholder(value: str) -> bool:
    """判断字面量是否为占位符而非真实凭据。"""
    lowered = value.strip().lower()
    if lowered in PLACEHOLDER_VALUES:
        return True
    # 含有占位类词汇（子串匹配，覆盖 your-api-key-here / my-placeholder-token 等）
    if any(token in lowered for token in PLACEHOLDER_TOKENS):
        return True
    # 全由同一字符重复构成（如 "xxxxxxxx"）
    if len(set(lowered)) <= 2 and len(lowered) >= 8:
        return True
    # 含明显的模板标记
    if any(marker in lowered for marker in ("<", ">", "{", "}", "$", "%", "{{", "}}")):
        return True
    return False


def _line_is_safe(line: str, match: re.Match) -> bool:
    """判断该行是否属于安全上下文。"""
    # 检查赋值左侧到匹配点之前的部分
    prefix = line[: match.start()]
    if any(ctx in prefix for ctx in SAFE_CONTEXTS):
        return True
    # 检查整行
    if any(ctx in line for ctx in ("os.environ", "os.getenv", "environ.get")):
        return True
    return False


def scan_paths(targets: list[Path]) -> list[tuple[str, int, str, str]]:
    """扫描给定目录/文件，返回 (文件, 行号, 变量名, 值) 列表。"""
    findings: list[tuple[str, int, str, str]] = []
    for target in targets:
        if not target.exists():
            continue
        files = [target] if target.is_file() else sorted(target.rglob("*.py"))
        for file in files:
            try:
                content = file.read_text(encoding="utf-8", errors="ignore")
            except OSError:
                continue
            for lineno, line in enumerate(content.splitlines(), start=1):
                if line.lstrip().startswith("#"):
                    continue
                for match in CREDENTIAL_ASSIGNMENT.finditer(line):
                    value = match.group("value")
                    if _is_placeholder(value):
                        continue
                    if _line_is_safe(line, match):
                        continue
                    name = match.group("name")
                    try:
                        rel = file.relative_to(REPO_ROOT)
                    except ValueError:
                        rel = file
                    findings.append((str(rel), lineno, name, value))
    return findings


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="检测源码中的硬编码凭据")
    parser.add_argument(
        "targets",
        nargs="*",
        default=["skills", "packages", "scripts"],
        help="要扫描的目录或文件",
    )
    args = parser.parse_args(argv)

    targets = [Path(t) for t in (args.targets or ["skills", "packages", "scripts"])]
    findings = scan_paths(targets)

    if findings:
        print("::error::检测到疑似硬编码凭据：")
        for path, lineno, name, value in findings:
            masked = value[:3] + "*" * max(len(value) - 6, 1) + value[-2:] if len(value) > 8 else "***"
            print(f"  {path}:{lineno}: {name} = {masked}")
        return 1

    print("未检测到硬编码凭据")
    return 0


if __name__ == "__main__":
    sys.exit(main())
