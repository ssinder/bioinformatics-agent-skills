# 技能改造规范契约

所有技能子 agent 必须严格遵守本文件。改造后的技能目录结构固定为：

```
skills/<skill-name>/
├── SKILL.md              # 必需，Agent Skills 标准格式
├── scripts/              # 可执行脚本，路径参数化 + logging 化
│   └── <entry>.py
├── references/           # 参考文档（可选但推荐）
│   └── *.md
├── requirements.txt      # 依赖（可选）
└── README.md             # 简短人读说明（可选）
```

## SKILL.md 格式（严格遵守）

YAML frontmatter 必须包含：

```yaml
---
name: <skill-name>          # 小写连字符，与目录名一致
description: <第三人称描述，必须包含 "This skill should be used when..." 句式，说明何时触发>
license: MIT
compatibility: <运行环境要求，如 Python 3.9+ 及关键依赖>
metadata:
  version: "1.0"
  last-reviewed: "2026-10-25"
  skill-author: sinder
---
```

description 要求：
- 第三人称，写成 "This skill should be used when..." / "用于...的场景"
- 必须同时包含 **做什么** 和 **何时触发**（含触发关键词，英文关键词必须出现，因为用户可能用英文提问）
- 长度控制在 100–400 字符
- 禁止出现 "Use this skill when" 这类第二人称表述

正文必须包含以下章节（用 `##` 二级标题）：

1. `## Overview` —— 这个技能解决什么问题，2–4 句
2. `## When to Use This Skill` —— 具体使用场景列表（用 `-`）
3. `## Quick Start` —— 可直接复制运行的命令示例
4. `## Workflow` —— 分步骤说明（编号列表或 mermaid 流程图）
5. `## Parameters` —— 关键参数说明表
6. `## Outputs` —— 产出文件说明表
7. `## Validation` —— 如何验证结果正确性（质控要点）
8. `## References` —— 指向 references/ 下的文档

## 脚本规范（严格遵守）

1. **导入共享层**：脚本开头加入
   ```python
   sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "packages"))
   from biopipe_core import get_logger, get_config, SkillError
   ```
   （`parents[3]` 是因为脚本在 `skills/<name>/scripts/x.py`，往上 3 层到仓库根）

2. **禁止硬编码绝对路径**：所有输入输出路径来自命令行参数，默认值基于 `get_config()`

3. **禁止裸 print**：用 `logger = get_logger("biopipe.<skill-name>")`，输出用 `logger.info(...)`

4. **统一 CLI**：用 `biopipe_core.cli.bootstrap_cli(description, skill_name, main)` 作为入口，
   `main(args)` 返回退出码。必须有 `--input`、`--output` 等业务参数。

5. **文件头**：每个脚本加 shebang + 编码声明 + 模块 docstring（含用法示例）

6. **必须有 `if __name__ == "__main__":` 入口**

7. **脚本必须能在没有输入文件时 `--help` 正常退出**（不依赖真实数据）

8. **不依赖重量级第三方库的模块要在导入时优雅降级**：
   例如 scanpy/matplotlib 缺失时给出清晰报错而不是 ImportError 堆栈

## 内容迁移来源

原仓库 `Bioinformatics_pipeline/<对应目录>/scripts/*.py` 是改造来源。
**保留原脚本的科学计算逻辑**，但要：
- 去掉硬编码路径，改为参数
- 去掉裸 print，改为 logger
- 补上 argparse
- 保留原有核心算法与注释

原目录位置：`C:/Users/DELL/WorkBuddy/2026-10-25-20-38-30/Bioinformatics_pipeline/`

## 禁止事项

- 不要修改原仓库 `Bioinformatics_pipeline/`（只读参考）
- 不要创建 `README.md` 以外的临时文件
- 不要在 SKILL.md 中写具体的人名路径、QQ、邮箱等隐私信息
- 不要使用 emoji（除 README 顶层外，SKILL.md 内不用）
- 不要在脚本里执行网络请求的默认行为（必须显式传参才联网）

## 交付自检

完成前必须自检并在返回中说明：
```bash
cd <repo>
python -c "import yaml,sys; d=yaml.safe_load(open('skills/<name>/SKILL.md').read().split('---')[1]); print(d['name'], '|', d['description'][:80])"
python skills/<name>/scripts/<entry>.py --help
```
