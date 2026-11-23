#!/usr/bin/env python3
"""校验每个技能的评测用例文件。

依据 Agent Skills 规范对评测的约定（``docs/skill-creation/evaluating-skills.mdx``）：
每个技能目录放 ``evals/evals.json``，包含若干测试用例。

这里检查的是**结构与覆盖度**，不跑模型：

1. 文件存在且是合法 JSON；
2. ``skill_name`` 与技能目录名一致；
3. 至少 3 个用例（官方建议的起步数量）；
4. 每个用例都有 ``id`` / ``prompt`` / ``expected_output``；
5. prompt 不是空泛的占位（要有具体场景）；
6. 至少一个 ``should_trigger`` 正例；
7. **至少一个 ``should_not_trigger`` 近似负例** —— 近似负例是防止技能互相
   抢触发的关键，缺了它就测不出"该不该用这个技能"；
8. ``expected_output`` 足够具体（不能是"给出正确答案"这类空话）。

用法：
    python scripts/ci_check_evals.py
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SKILLS_DIR = REPO_ROOT / "skills"

MIN_CASES = 3

# expected_output 太短说明写得太笼统，测不出东西
MIN_EXPECTED_LEN = 40

# prompt 至少要有一定长度，说明场景是具体的
MIN_PROMPT_LEN = 30


def _skill_dirs() -> list[Path]:
    if not SKILLS_DIR.is_dir():
        return []
    return sorted(
        d for d in SKILLS_DIR.iterdir() if d.is_dir() and not d.name.startswith((".", "_"))
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="校验技能评测用例")
    parser.add_argument("--json", action="store_true", help="输出 JSON")
    parser.add_argument("--require-all", action="store_true", help="要求每个技能都有评测文件")
    args = parser.parse_args(argv)

    problems: list[str] = []
    warnings: list[str] = []
    stats: list[dict] = []

    skill_dirs = _skill_dirs()
    if not skill_dirs:
        print(f"未找到任何技能: {SKILLS_DIR}", file=sys.stderr)
        return 2

    for skill_dir in skill_dirs:
        name = skill_dir.name
        evals_path = skill_dir / "evals" / "evals.json"

        if not evals_path.is_file():
            message = f"{name}: 缺少 evals/evals.json"
            (problems if args.require_all else warnings).append(message)
            continue

        try:
            data = json.loads(evals_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError as exc:
            problems.append(f"{name}: evals/evals.json 不是合法 JSON（{exc}）")
            continue

        if not isinstance(data, dict):
            problems.append(f"{name}: evals/evals.json 顶层必须是对象")
            continue

        if data.get("skill_name") != name:
            problems.append(
                f"{name}: skill_name ({data.get('skill_name')!r}) 与目录名不一致"
            )

        cases = data.get("evals")
        if not isinstance(cases, list):
            problems.append(f"{name}: evals 必须是列表")
            continue

        if len(cases) < MIN_CASES:
            problems.append(f"{name}: 只有 {len(cases)} 个用例，至少需要 {MIN_CASES} 个")

        seen_ids: set = set()
        positives = 0
        negatives = 0

        for index, case in enumerate(cases, 1):
            if not isinstance(case, dict):
                problems.append(f"{name}: 第 {index} 个用例不是对象")
                continue

            case_id = case.get("id")
            if case_id is None:
                problems.append(f"{name}: 第 {index} 个用例缺少 id")
            elif case_id in seen_ids:
                problems.append(f"{name}: 用例 id {case_id!r} 重复")
            else:
                seen_ids.add(case_id)

            prompt = str(case.get("prompt", ""))
            if not prompt.strip():
                problems.append(f"{name}: 用例 {case_id} 缺少 prompt")
            elif len(prompt) < MIN_PROMPT_LEN:
                warnings.append(
                    f"{name}: 用例 {case_id} 的 prompt 过短（{len(prompt)} 字符），"
                    f"可能不够具体"
                )

            expected = str(case.get("expected_output", ""))
            if not expected.strip():
                problems.append(f"{name}: 用例 {case_id} 缺少 expected_output")
            elif len(expected) < MIN_EXPECTED_LEN:
                warnings.append(
                    f"{name}: 用例 {case_id} 的 expected_output 过短"
                    f"（{len(expected)} 字符），可能测不出行为差异"
                )

            assertions = case.get("assertions", [])
            if not isinstance(assertions, list):
                problems.append(f"{name}: 用例 {case_id} 的 assertions 必须是列表")
                continue

            for assertion in assertions:
                if not isinstance(assertion, dict):
                    problems.append(f"{name}: 用例 {case_id} 的断言不是对象")
                    continue
                kind = assertion.get("type")
                if kind == "should_trigger" and assertion.get("value") is not False:
                    positives += 1
                elif kind == "should_not_trigger":
                    negatives += 1
                    if assertion.get("value") is not True:
                        warnings.append(
                            f"{name}: 用例 {case_id} 的 should_not_trigger 应显式写 value: true"
                        )
                elif kind in ("output_contains", "behavior"):
                    if not str(assertion.get("value", "")).strip():
                        problems.append(
                            f"{name}: 用例 {case_id} 的 {kind} 断言缺少 value"
                        )
                else:
                    problems.append(
                        f"{name}: 用例 {case_id} 含未知断言类型 {kind!r}"
                    )

        if positives == 0:
            problems.append(f"{name}: 没有任何 should_trigger 正例")
        if negatives == 0:
            problems.append(
                f"{name}: 没有任何 should_not_trigger 近似负例——"
                f"无法验证该技能不会抢走相邻技能的请求"
            )

        stats.append(
            {
                "skill": name,
                "cases": len(cases),
                "positive": positives,
                "negative": negatives,
            }
        )

    result = {
        "skills_with_evals": len(stats),
        "total_skills": len(skill_dirs),
        "stats": stats,
        "problems": problems,
        "warnings": warnings,
    }

    if args.json:
        print(json.dumps(result, ensure_ascii=False, indent=2))
        return 1 if problems else 0

    for p in problems:
        print(f"[ERROR] {p}")
    for w in warnings:
        print(f"[WARN ] {w}")

    total_cases = sum(s["cases"] for s in stats)
    print(f"\n评测用例检查：{len(stats)}/{len(skill_dirs)} 个技能，共 {total_cases} 个用例")
    if problems:
        print("::error::评测用例检查失败")
        return 1
    print("评测用例检查通过")
    return 0


if __name__ == "__main__":
    sys.exit(main())
