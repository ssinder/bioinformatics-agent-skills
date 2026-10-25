#!/usr/bin/env python3
"""静态校验文档站配置（无需 Ruby 环境）。

背景
----
本机没有 Ruby/Jekyll，文档站只能靠 CI 构建。而 Jekyll 构建失败时
错误信息往往晦涩（Liquid 语法、插件白名单、baseurl 拼接……）。
本脚本在提交阶段就把可静态发现的问题挑出来。

检查项
------
1. ``docs/_config.yml`` 存在且可被 YAML 解析。
2. 声明的插件都在 **GitHub Pages 官方白名单**内。
   非白名单插件会导致线上构建失败（本地 jekyll build 反而能通过），
   这是最常见也最隐蔽的一类问题。
3. 每个 Markdown 文档的 front matter 可解析（``---`` 成对、YAML 合法）。
4. 站内链接指向的文件确实存在（``/xxx.html`` → ``docs/xxx.md``）。
5. 不残留非白名单插件特有的 Liquid 标签
   （``include_relative`` 是 Jekyll 内置的，不在此列）。

用法：
    python scripts/ci_docs_check.py
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
DOCS_DIR = REPO_ROOT / "docs"

# GitHub Pages 官方支持的插件白名单
#
# 来源：https://docs.github.com/en/pages/using-github-pages/using-custom-actions-with-github-pages
# 的 "Pages dependencies" 章节。注意这份列表**不含** jekyll-feed 与
# jekyll-seo-tag —— 它们不在 Pages 自带依赖里，但作为普通 gem 依赖
# 可以通过 Gemfile 正常安装并启用，因此这里额外允许。
ALLOWED_PLUGINS = {
    # Pages 内置依赖
    "jekyll-relative-links",
    "jekyll-optional-front-matter",
    "jekyll-readme-index",
    "jekyll-default-layout",
    "jekyll-titles-from-headings",
    "jekyll-github-metadata",
    # 主题
    "jekyll-theme-primer",
    "jekyll-theme-cayman",
    "jekyll-theme-minimal",
    "jekyll-theme-modern",
    "jekyll-theme-hacker",
    # 可通过 Gemfile 安装的普通 gem
    "jekyll-feed",
    "jekyll-seo-tag",
    "jekyll-sitemap",
    "jekyll-gist",
    "jekyll-paginate",
    "jekyll-theme-jekyll-theme-primer",
    "jekyll-theme-jekyll-theme-cayman",
}

# 非白名单插件才会提供的 Liquid 标签
PLUGIN_ONLY_TAGS = ("readme_index",)


def _fail(problems: list[str], message: str) -> None:
    problems.append(message)


def _check_config(problems: list[str]) -> None:
    config_path = DOCS_DIR / "_config.yml"
    if not config_path.is_file():
        _fail(problems, "docs/_config.yml 不存在")
        return

    try:
        import yaml
    except ImportError:
        print("::warning::未安装 PyYAML，跳过 _config.yml 解析检查")
        yaml = None  # type: ignore[assignment]

    if yaml is not None:
        try:
            config = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
        except Exception as exc:
            _fail(problems, f"docs/_config.yml 解析失败: {exc}")
            return

        plugins = config.get("plugins") or []
        if isinstance(plugins, str):
            plugins = [plugins]
        for plugin in plugins:
            if plugin not in ALLOWED_PLUGINS:
                _fail(
                    problems,
                    f"docs/_config.yml 声明了非白名单插件 '{plugin}'。"
                    "GitHub Pages 只支持官方白名单内的插件，线上构建会失败。",
                )


def _check_front_matter(problems: list[str]) -> None:
    """检查每个 Markdown 的 front matter 是否成对且 YAML 合法。"""
    try:
        import yaml
    except ImportError:
        print("::warning::未安装 PyYAML，跳过 front matter 解析检查")
        return

    for md in sorted(DOCS_DIR.rglob("*.md")):
        text = md.read_text(encoding="utf-8")
        if not text.startswith("---"):
            continue  # 无 front matter 的 Markdown 是合法的
        end = text.find("\n---", 3)
        if end == -1:
            _fail(problems, f"{md.relative_to(REPO_ROOT)}: front matter 起始 --- 没有闭合")
            continue
        block = text[3:end]
        try:
            yaml.safe_load(block)
        except Exception as exc:
            _fail(problems, f"{md.relative_to(REPO_ROOT)}: front matter YAML 非法: {exc}")


def _check_links(problems: list[str]) -> None:
    """检查站内链接指向的文件是否存在。"""
    # {{ '/xxx.html' | relative_url }} 形式
    pattern = re.compile(r"\{\{\s*['\"](/[^'\"]+?)['\"]\s*\|\s*relative_url\s*\}\}")
    # 纯 Markdown 相对链接
    md_pattern = re.compile(r"\]\((?!https?://|#)([^)\s]+?)\.md\)")

    for md in sorted(DOCS_DIR.rglob("*.md")):
        text = md.read_text(encoding="utf-8")
        rel = md.relative_to(REPO_ROOT)

        for match in pattern.finditer(text):
            url = match.group(1)
            target = DOCS_DIR / (url.lstrip("/").replace(".html", ".md"))
            if not target.exists():
                _fail(problems, f"{rel}: 站内链接 {url} 指向不存在的文件")

        for match in md_pattern.finditer(text):
            link = match.group(1)
            target = (md.parent / f"{link}.md").resolve()
            if not target.exists():
                _fail(problems, f"{rel}: 相对链接 {link}.md 指向不存在的文件")


def _check_plugin_tags(problems: list[str]) -> None:
    """检查是否使用了非白名单插件才提供的 Liquid 标签。"""
    for md in sorted(DOCS_DIR.rglob("*.md")):
        text = md.read_text(encoding="utf-8")
        rel = md.relative_to(REPO_ROOT)
        for tag in PLUGIN_ONLY_TAGS:
            if re.search(r"\{%-?\s*" + tag + r"\b", text):
                _fail(
                    problems,
                    f"{rel}: 使用了 '{tag}' 标签，它由非白名单插件提供，线上构建会失败",
                )


def main() -> int:
    problems: list[str] = []

    if not DOCS_DIR.is_dir():
        print("::warning::docs/ 目录不存在，跳过文档站检查")
        return 0

    _check_config(problems)
    _check_front_matter(problems)
    _check_links(problems)
    _check_plugin_tags(problems)

    if problems:
        print("文档站检查未通过：")
        for p in problems:
            print(f"  - {p}")
        print("::error::文档站静态检查失败")
        return 1

    md_count = len(list(DOCS_DIR.rglob("*.md")))
    print(f"文档站静态检查通过（{md_count} 个 Markdown 文件）")
    return 0


if __name__ == "__main__":
    sys.exit(main())
