#!/usr/bin/env python3
"""生成一个符合规范的技能骨架。

用法：
    python scripts/init_skill.py my-new-skill
    python scripts/init_skill.py my-new-skill --category 转录组学 --force

生成结构：
    skills/my-new-skill/
    ├── SKILL.md
    ├── scripts/run.py
    ├── references/methodology.md
    ├── references/troubleshooting.md
    ├── requirements.txt
    └── README.md
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SKILLS_DIR = REPO_ROOT / "skills"

SKILL_MD_TEMPLATE = """---
name: {name}
description: This skill should be used when a user needs to {summary}, including analysis of {data_desc} with reproducible and validated scientific methods. 用于{description_cn}的场景。
license: MIT
compatibility: Requires Python 3.9+ with numpy and pandas.
metadata:
  version: "1.0"
  last-reviewed: "{today}"
  skill-author: {author}
---

# {title}

## Overview

{summary}.

本技能封装了 {data_desc} 的分析流程,提供可复现的计算步骤、参数选择的科学依据,
以及结果正确性的验证方法。

## When to Use This Skill

Use this skill when the user wants to:
- 对 {data_desc} 进行 {summary_cn}
- 批量处理多个样本并保持参数一致
- 需要了解参数选择背后的科学依据,而不仅仅是拿到结果
- 需要把分析过程记录成可复现的流程

对于**其他**类型的数据,请使用对应的专门技能。技能名以 `biopipe --list` 输出为准。

## Quick Start

```bash
export PYTHONPATH="$PWD/packages"

# 查看完整用法
python skills/{name}/scripts/{entry} --help

# 基本用法
python skills/{name}/scripts/{entry} \\
    --input {example_input} \\
    --output {example_output}
```

## Workflow

```mermaid
flowchart TD
    A["准备输入数据"] --> B["质控与检查"]
    B --> C["核心分析步骤"]
    C --> D["结果汇总"]
    D --> E["验证与出图"]
    E --> F["输出最终报告"]
```

1. **准备输入数据** —— 确认数据格式与质量,记录来源与版本
2. **质控与检查** —— 在进入核心步骤前完成基础检查,记录关键指标
3. **核心分析步骤** —— 按标准流程执行,固定随机种子
4. **结果汇总** —— 输出结构化的中间结果与统计量
5. **验证与出图** —— 按 Validation 章节的方法确认结果可靠
6. **输出最终报告** —— 汇总所有产物并记录运行参数

## Parameters

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `--input` | Path | 必填 | 输入文件路径 |
| `--output` | Path | `results/` | 输出目录 |
| `--config` | Path | 无 | 项目配置文件(JSON 或 key=value) |
| `--log-file` | Path | 无 | 日志输出文件 |
| `-v/--verbose` | flag | off | 输出 DEBUG 级日志 |
| `-q/--quiet` | flag | off | 仅输出警告与错误 |

## Outputs

| 文件 | 说明 |
|------|------|
| `results/{name}_report.txt` | 运行摘要与关键统计量 |
| `results/{name}_results.csv` | 结构化结果表 |
| `figures/` | 输出图表 |
| `logs/` | 完整运行日志 |

## Validation

结果可靠性检查清单:

- [ ] 输入样本数量与预期一致,无缺失
- [ ] 关键指标落在合理区间内(记录本次运行的实测值)
- [ ] 结果行数与输入记录数匹配
- [ ] 随机性来源已固定(记录随机种子)
- [ ] 异常值已标记而非静默丢弃

**注意**:本技能输出的是辅助分析结果,需要结合领域知识与实验设计做最终判断。

## References

- [方法学原理](references/methodology.md) —— 参数选择的依据与科学背景
- [常见问题与排错](references/troubleshooting.md) —— 故障现象、根因与解决方案
"""

SCRIPT_TEMPLATE = '''#!/usr/bin/env python3
"""{title} —— 入口脚本。

用法示例：
    python skills/{name}/scripts/{entry} --help
    python skills/{name}/scripts/{entry} --input example.txt --output results/
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "packages"))

from biopipe_core import bootstrap_cli, get_config, get_logger


def main(args) -> int:
    """技能主入口。

    Args:
        args: 解析后的命令行参数。

    Returns:
        退出码,0 表示成功。
    """
    logger = get_logger("biopipe.{name}")
    config = get_config(args.config)

    logger.info("技能启动: {name}")
    logger.info("项目根目录: %s", config.root)
    logger.info("输出目录: %s", config.output_dir)

    # TODO: 实现核心逻辑
    logger.info("核心逻辑待实现")
    return 0


if __name__ == "__main__":
    sys.exit(bootstrap_cli("{title}", "{name}", main))
'''

METHODOLOGY_TEMPLATE = """# {title} —— 方法学原理

本文档说明本技能各步骤的**科学依据**与**参数选择理由**,供研究者判断结果是否合理。

## 概述

{summary}。

分析流程遵循领域内的标准实践,关键决策点如下。

## 关键参数与依据

### 参数一

- **取值范围**:合理区间与含义
- **选择依据**:科学原因
- **敏感性与检查方法**:如何验证该选择是否影响结论

### 参数二

- **取值范围**:合理区间与含义
- **选择依据**:科学原因
- **敏感性与检查方法**:如何验证该选择是否影响结论

## 方法局限

诚实说明本方法的边界:

- 什么情况下结果不可靠
- 哪些混杂因素未被考虑
- 需要哪些额外验证才能支撑结论

## 参考文献

- 领域内的标准方法学文献
"""

TROUBLESHOOTING_TEMPLATE = """# {title} —— 常见问题与排错

## 输入与数据问题

### 症状:读取输入文件失败

**可能原因**

1. 文件路径不存在或拼写错误
2. 文件编码不是 UTF-8
3. 文件格式与预期不符(表头列名变化)

**解决方案**

1. 用 `--verbose` 运行,查看实际解析的路径
2. 转换文件编码:`iconv -f gbk -t utf-8 input.csv > input_utf8.csv`
3. 检查表头列名是否与文档一致

## 计算与结果问题

### 症状:结果为空或行数异常

**排查步骤**

1. 确认输入记录数
2. 检查是否有记录被质控步骤过滤
3. 查看日志中过滤比例,过滤率过高通常说明参数不合适

### 症状:数值结果明显异常

**可能原因**

1. 参数超出合理范围
2. 输入数据格式错误
3. 版本差异导致的算法行为变化

**解决方案**

1. 对照 `references/methodology.md` 中的参数表
2. 用最小样例数据验证
3. 记录环境版本,与参考结果对比

## 性能问题

### 症状:运行时间过长或内存占用高

**优化建议**

1. 分批处理
2. 降低默认的采样密度
3. 检查是否有不必要的中间结果落盘

## 如何提问

遇到未收录的问题时,请在 Issue 中提供:

- 完整的命令行与 `--verbose` 日志
- 输入数据的样例(脱敏后)
- 期望结果与实际结果的差异
- Python 与关键依赖版本
"""

REQUIREMENTS_TEMPLATE = """# {title} 依赖
requests>=2.28
numpy>=1.23
pandas>=1.5
"""

README_TEMPLATE = """# {title}

{summary}。

## 快速开始

```bash
export PYTHONPATH="$PWD/packages"
python skills/{name}/scripts/{entry} --help
```

完整说明见 [SKILL.md](SKILL.md)。

- 方法学原理:[references/methodology.md](references/methodology.md)
- 排错手册:[references/troubleshooting.md](references/troubleshooting.md)
"""


def main(argv: list[str] | None = None) -> int:
    """生成技能骨架。"""
    parser = argparse.ArgumentParser(description="生成符合 Agent Skills 规范的技能骨架")
    parser.add_argument("name", help="技能名(kebab-case)")
    parser.add_argument("--summary", default="perform domain-specific analysis", help="技能功能简述")
    parser.add_argument("--summary-cn", default="领域特定分析", help="技能功能简述(中文)")
    parser.add_argument("--description-cn", default="领域特定分析", help="description 中的中文表述")
    parser.add_argument("--data-desc", default="domain data", help="处理的数据类型描述")
    parser.add_argument("--title", default=None, help="技能标题")
    parser.add_argument("--category", default="其他", help="所属分类")
    parser.add_argument("--author", default="sinder", help="作者")
    parser.add_argument("--entry", default="run.py", help="入口脚本文件名")
    parser.add_argument("--force", action="store_true", help="覆盖已存在的目录")
    args = parser.parse_args(argv)

    if not re.match(r"^[a-z0-9]+(?:-[a-z0-9]+)*$", args.name):
        print(f"技能名不合法: {args.name}(仅允许小写字母、数字与连字符)", file=sys.stderr)
        return 2

    skill_dir = SKILLS_DIR / args.name
    if skill_dir.exists():
        if not args.force:
            print(f"目录已存在: {skill_dir}(使用 --force 覆盖)", file=sys.stderr)
            return 1
        import shutil

        shutil.rmtree(skill_dir)

    title = args.title or args.name.replace("-", " ").title()
    from datetime import date

    today = date.today().isoformat()

    (skill_dir / "scripts").mkdir(parents=True)
    (skill_dir / "references").mkdir(parents=True)

    ctx = {
        "name": args.name,
        "title": title,
        "summary": args.summary,
        "summary_cn": args.summary_cn,
        "description_cn": args.description_cn,
        "data_desc": args.data_desc,
        "author": args.author,
        "entry": args.entry,
        "today": today,
        "example_input": "input/data.csv",
        "example_output": "results/",
    }

    (skill_dir / "SKILL.md").write_text(SKILL_MD_TEMPLATE.format(**ctx), encoding="utf-8")
    (skill_dir / "scripts" / args.entry).write_text(SCRIPT_TEMPLATE.format(**ctx), encoding="utf-8")
    (skill_dir / "references" / "methodology.md").write_text(
        METHODOLOGY_TEMPLATE.format(**ctx), encoding="utf-8"
    )
    (skill_dir / "references" / "troubleshooting.md").write_text(
        TROUBLESHOOTING_TEMPLATE.format(**ctx), encoding="utf-8"
    )
    (skill_dir / "requirements.txt").write_text(REQUIREMENTS_TEMPLATE.format(**ctx), encoding="utf-8")
    (skill_dir / "README.md").write_text(README_TEMPLATE.format(**ctx), encoding="utf-8")

    print(f"技能骨架已生成: {skill_dir}")
    print(f"分类: {args.category}")
    print("\n下一步:")
    print(f"  1. 编辑 {skill_dir / 'SKILL.md'} 的 description(决定何时被触发)")
    print(f"  2. 实现 {skill_dir / 'scripts' / args.entry} 的核心逻辑")
    print("  3. 编写 references/ 下的两份参考文档")
    print(f"  4. 校验: python scripts/validate_skills.py --skill {args.name} --strict")
    return 0


if __name__ == "__main__":
    sys.exit(main())
