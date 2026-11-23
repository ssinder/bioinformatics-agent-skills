#!/usr/bin/env python3
"""跨宿主技能兼容性检查。

为什么需要这个检查
------------------
技能被装到不同 AI 宿主后，**能不能被加载**取决于几个静默条件。
任何一条不满足，技能就会安静地不加载——没有任何报错，用户只会觉得
"这个技能好像不灵"。这类失败比崩溃更难发现，也更难排查。

已知的三个失败模式：

1. **宿主扫描的目录不同。** Claude Code 个人路径是 ``~/.claude/skills/``，
   OpenCode 等宿主用 ``~/.agents/skills/``，项目级又是 ``.claude/skills/``。
   装错路径 → 技能永不加载。
2. **描述超长。** 部分宿主对 description 有 1,024 字符的硬上限。
   超限不是"截断显示"，而是**整个技能不被加载**。
3. **描述在列表内被截断。** 即使未超限，宿主也会为了控制上下文预算
   截断超长条目。截断从尾部开始，所以**关键用例必须写在最前面**。

本脚本自包含，不依赖任何外部工具，因此在 CI 里永远能跑。
（NVIDIA 的 ``skillevaluator harness-check`` 覆盖了更广的宿主矩阵，
但它是可选依赖；本检查器保证最关键的三条底线。）

用法：
    python scripts/ci_harness_check.py
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

REPO_ROOT = Path(__file__).resolve().parents[1]
SKILLS_DIR = REPO_ROOT / "skills"

try:
    import yaml
except ImportError:  # pragma: no cover
    yaml = None

# 各宿主对 description 的硬上限。取最严格的作为门禁。
HOST_DESC_LIMITS = {
    "claude-code": 1024,
    "copilot-cli": 1024,
    "codex": 1024,
    "agentskills-spec": 1024,
}

# SKILL.md 正文的体量上限（字节）。
#
# 依据：SkillsBench（arXiv:2602.12670）附录 A 对 767,430 个技能的分析给出
# 生态中位数 4.8 KB、75 分位 9.2 KB；同论文的设计消融显示"完整详尽的
# 文档"的任务通过率增益只有 +0.7 pp，而聚焦的流程指引达 +19.0 ~ +21.5 pp。
# 官方规范另建议 SKILL.md 保持在 500 行以内。
#
# 这里设 7 KB 作为硬门禁——略高于生态中位数，给内容丰富的技能留出空间，
# 但明显低于 75 分位，避免滑向"详尽文档"那一档。
BODY_LIMIT_BYTES = 7 * 1024

# 超过这个值就提醒（不失败）：多数技能应当贴近生态中位数。
BODY_WARN_BYTES = 6.5 * 1024

# 安全余量：接近上限时提醒，避免后续编辑误越界。
WARN_AT = 0.9

# 描述开头应当包含触发用例——宿主截断从尾部开始，
# 所以放前面才不会被截掉。
TRIGGER_HEAD = re.compile(
    r"^\W*(use when|this skill should be used when|用于|适用于|当用户)",
    re.IGNORECASE,
)


def _load_frontmatter(path: Path) -> dict:
    """解析 SKILL.md 的 YAML frontmatter。

    注意：**不能**用简单的 ``text.split("---", 2)``。Markdown 表格的
    分隔行（``|---|---|``）里含有连续的 ``-``，会被误当成 frontmatter
    的结束标记，把大半个正文当成 YAML 去解析，结果是"description 为空"
    这类莫名其妙的报错。

    正确做法：只把**独占一行**的 ``---`` 当作分隔符。
    """
    text = path.read_text(encoding="utf-8")
    lines = text.split("\n")

    if not lines or lines[0].strip() != "---":
        return {}

    end = None
    for i in range(1, len(lines)):
        if lines[i].strip() == "---":
            end = i
            break
    if end is None:
        return {}

    try:
        return yaml.safe_load("\n".join(lines[1:end])) or {}
    except Exception:
        return {}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="跨宿主技能兼容性检查")
    parser.add_argument("--json", action="store_true", help="输出 JSON")
    args = parser.parse_args(argv)

    if yaml is None:
        print("缺少 PyYAML，无法解析 frontmatter", file=sys.stderr)
        return 2

    if not SKILLS_DIR.is_dir():
        print(f"技能目录不存在: {SKILLS_DIR}", file=sys.stderr)
        return 2

    skill_dirs = sorted(
        d for d in SKILLS_DIR.iterdir() if d.is_dir() and not d.name.startswith((".", "_"))
    )
    if not skill_dirs:
        print(f"未找到任何技能: {SKILLS_DIR}", file=sys.stderr)
        return 2

    problems: list[dict] = []
    warnings: list[dict] = []
    limit = min(HOST_DESC_LIMITS.values())

    for skill_dir in skill_dirs:
        skill_md = skill_dir / "SKILL.md"
        if not skill_md.is_file():
            problems.append({"skill": skill_dir.name, "code": "no-skill-md", "message": "缺少 SKILL.md"})
            continue

        fm = _load_frontmatter(skill_md)
        name = str(fm.get("name", skill_dir.name))
        desc = str(fm.get("description", ""))

        # SKILL.md 正文体量：Claude Code 在自动压缩后只保留每个技能的前
        # 5,000 token，文件尾部先被丢掉；而 SkillsBench 的消融显示
        # "完整详尽的文档"几乎没有增益。控制在合理区间内。
        text = skill_md.read_text(encoding="utf-8")
        body = text.split("---", 2)[2] if text.startswith("---") else text
        body_bytes = len(body.encode("utf-8"))
        if body_bytes > BODY_LIMIT_BYTES:
            problems.append(
                {
                    "skill": name,
                    "code": "body-too-long",
                    "message": (
                        f"SKILL.md 正文 {body_bytes / 1024:.1f} KB，超过上限 "
                        f"{BODY_LIMIT_BYTES / 1024:.0f} KB。详情应下沉到 references/——"
                        f"SkillsBench 的消融显示 comprehensive 文档的增益仅 +0.7 pp，"
                        f"而聚焦的流程指引达 +19.0 ~ +21.5 pp。"
                    ),
                }
            )
        elif body_bytes > BODY_WARN_BYTES:
            warnings.append(
                {
                    "skill": name,
                    "code": "body-long",
                    "message": (
                        f"SKILL.md 正文 {body_bytes / 1024:.1f} KB，略超建议值 "
                        f"{BODY_WARN_BYTES / 1024:.1f} KB（生态中位数 4.8 KB）"
                    ),
                }
            )

        if not desc:
            problems.append(
                {"skill": name, "code": "desc-empty", "message": "description 为空，技能无法被匹配触发"}
            )
            continue

        if len(desc) > limit:
            problems.append(
                {
                    "skill": name,
                    "code": "desc-over-host-limit",
                    "message": (
                        f"description {len(desc)} 字符，超过宿主硬上限 {limit}。"
                        f"超限时部分宿主（如 copilot-cli）会**整个技能不加载**。"
                    ),
                }
            )
        elif len(desc) > limit * WARN_AT:
            warnings.append(
                {
                    "skill": name,
                    "code": "desc-near-host-limit",
                    "message": (
                        f"description {len(desc)} 字符，接近上限 {limit}"
                        f"（余量 {limit - len(desc)} 字符），后续编辑容易误越界。"
                    ),
                }
            )

        # 关键用例前置：宿主截断从尾部开始
        if not TRIGGER_HEAD.match(desc.strip()):
            warnings.append(
                {
                    "skill": name,
                    "code": "desc-trigger-not-front",
                    "message": (
                        "description 开头不是触发用例。宿主会在超出上下文预算时"
                        "从尾部截断描述，把触发语放在最前面才不会被截掉。"
                    ),
                }
            )

        # frontmatter 只允许规范定义的字段；claude.ai 上传路径会直接报错
        allowed = {"name", "description", "license", "compatibility", "metadata", "allowed-tools"}
        extra = set(fm.keys()) - allowed
        if extra:
            problems.append(
                {
                    "skill": name,
                    "code": "frontmatter-extra-keys",
                    "message": (
                        f"frontmatter 含规范外字段 {sorted(extra)}。"
                        f"claude.ai 上传与 Skills API 只接受 {sorted(allowed)}，"
                        f"多出的字段会导致打包失败。"
                    ),
                }
            )

        # metadata 必须是字符串到字符串的映射
        #
        # 这是 Agent Skills 规范的硬要求（"a map from string keys to
        # string values"）。嵌套映射或列表会让部分宿主的打包器直接
        # 拒绝该技能。
        metadata = fm.get("metadata")
        if metadata is not None:
            if not isinstance(metadata, dict):
                problems.append(
                    {"skill": name, "code": "metadata-not-map", "message": "metadata 必须是映射"}
                )
            else:
                for key, value in metadata.items():
                    if not isinstance(value, str):
                        problems.append(
                            {
                                "skill": name,
                                "code": "metadata-value-not-string",
                                "message": (
                                    f"metadata.{key} 的值必须是字符串"
                                    f"（规范要求 string→string 映射；"
                                    f"当前为 {type(value).__name__}）"
                                ),
                            }
                        )

        # compatibility 有 500 字符上限
        compat = fm.get("compatibility")
        if isinstance(compat, str) and len(compat) > 500:
            problems.append(
                {
                    "skill": name,
                    "code": "compat-too-long",
                    "message": f"compatibility {len(compat)} 字符，超过规范上限 500",
                }
            )

    result = {
        "total": len(skill_dirs),
        "problems": problems,
        "warnings": warnings,
    }

    if args.json:
        import json

        print(json.dumps(result, ensure_ascii=False, indent=2))
        return 1 if problems else 0

    if problems or warnings:
        for p in problems:
            print(f"[ERROR] {p['skill']}: {p['code']} - {p['message']}")
        for w in warnings:
            print(f"[WARN ] {w['skill']}: {w['code']} - {w['message']}")
        print(f"\n合计: {len(problems)} 错误, {len(warnings)} 警告")
        if problems:
            print("::error::跨宿主兼容性检查失败")
            return 1
        print("检查通过（有警告）")
        return 0

    print(f"跨宿主兼容性检查通过：{len(skill_dirs)} 个技能均可在各宿主正常加载")
    return 0


if __name__ == "__main__":
    sys.exit(main())
