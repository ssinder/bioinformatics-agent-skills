# 贡献指南 / Contributing Guide

感谢你考虑为 BioinfoAgent Skills 贡献技能。本项目遵循 [Agent Skills 开放标准](https://agentskills.io/)，
欢迎任何领域的科研工作者提交技能。

English version available at the bottom of this document.

---

## 贡献什么

我们欢迎三类贡献：

| 类型 | 说明 | 适合谁 |
|------|------|--------|
| **新技能 (New Skill)** | 提交一个全新的领域技能 | 有特定分析流程的科研工作者 |
| **改进技能 (Improve Skill)** | 修复 bug、补充参数、改进文档 | 所有用户 |
| **基础设施 (Infrastructure)** | 改进共享工具层、CI、文档 | 想参与工程建设的开发者 |

## 快速开始贡献一个新技能

### 1. Fork 并创建分支

```bash
git clone https://github.com/ssinder/bioinformatics-agent-skills.git
cd bioinformatics-agent-skills
git checkout -b feature/my-new-skill
```

### 2. 生成技能骨架

```bash
python scripts/init_skill.py my-new-skill
```

这会创建 `skills/my-new-skill/` 目录，含 `SKILL.md` 模板、`scripts/`、`references/`、`assets/`。

### 3. 编写 SKILL.md

**这是最重要的一步。** SKILL.md 决定 AI 何时触发你的技能。

```yaml
---
name: my-new-skill
description: This skill should be used when a user needs to perform X analysis on Y data, including A and B steps. 用于对 Y 数据进行 X 分析的场景。
license: MIT
compatibility: Requires Python 3.9+ with numpy and pandas.
metadata:
  version: "1.0"
  last-reviewed: "2026-10-25"
  skill-author: 你的名字
---
```

description 的三条铁律：
1. **第三人称** —— 写 "This skill should be used when..."，不要写 "You should use this..."
2. **说清何时触发** —— 包含用户可能说的关键词（中英文都要有，用户可能用英文提问）
3. **具体** —— 列出数据格式、工具名、分析步骤，不要写空话

### 4. 编写脚本

所有脚本必须遵守项目约定：

```python
#!/usr/bin/env python3
"""技能说明与用法示例。"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "packages"))

from biopipe_core import get_logger, get_config, bootstrap_cli


def main(args) -> int:
    logger = get_logger("biopipe.my-new-skill")
    config = get_config(args.config)
    logger.info("处理输入: %s", args.input)
    return 0


if __name__ == "__main__":
    sys.exit(bootstrap_cli("描述", "my-new-skill", main))
```

**硬性规则**：
- ❌ 禁止硬编码绝对路径（`/Users/xxx/...`、`C:/Users/xxx/...`）
- ❌ 禁止裸 `print()`，一律用 `logger`
- ✅ 必须用 `argparse` 提供命令行接口
- ✅ 必须有 `if __name__ == "__main__":` 入口
- ✅ 重量级依赖（scanpy、OpenMM 等）要优雅降级，给出清晰报错

### 5. 本地校验

```bash
python scripts/validate_skills.py --skill my-new-skill
python scripts/validate_skills.py --skill my-new-skill --strict   # 警告也算失败
```

必须零错误零警告才能合并。

### 6. 提交 PR

```bash
git add .
git commit -m "feat: add my-new-skill for X analysis"
git push origin feature/my-new-skill
```

打开 Pull Request。CI 会自动运行规范校验、语法检查、依赖审计与冒烟测试。

## 质量要求

| 要求 | 说明 |
|------|------|
| **科学正确性** | 算法必须真实可用，不能是伪代码占位 |
| **可复现** | 固定随机种子，记录参数与版本 |
| **有降级路径** | 核心算法应内置纯 Python 实现，避免强依赖重库 |
| **文档完整** | 必须提供 `references/methodology.md`（方法学原理）与 `references/troubleshooting.md`（排错） |
| **有验证环节** | SKILL.md 的 Validation 章节说明如何验证结果正确 |

## 需要帮助？

- 技能规范疑问 → 提 Issue，标签 `question`
- 不知道怎么写 description → 提 Issue 附上你的分析场景，我帮你写
- 想要更多示例 → 见 `docs/skill-guides/`

---

# English

## Contributing a New Skill

We welcome contributions from researchers and developers. Three kinds:

- **New skills** — a skill for a domain you work in
- **Improvements** — bug fixes, extra parameters, better docs
- **Infrastructure** — shared tooling, CI, documentation

### Quick start

```bash
git checkout -b feature/my-new-skill
python scripts/init_skill.py my-new-skill
# edit skills/my-new-skill/SKILL.md and scripts/
python scripts/validate_skills.py --skill my-new-skill --strict
git commit -m "feat: add my-new-skill"
git push origin feature/my-new-skill
```

### SKILL.md rules

1. Third person: "This skill should be used when..." — never "You should..."
2. State clearly **when** to trigger, with keywords users would actually type (include both English and Chinese)
3. Be specific: name the data formats, tools, and steps

### Script rules

- No hardcoded absolute paths
- Use `logger` instead of bare `print()`
- Use `argparse` and provide a `__main__` entry point
- Heavy dependencies (scanpy, OpenMM, ...) must degrade gracefully
- Import shared helpers from `biopipe_core`

### Quality bar

Zero errors and zero warnings from `scripts/validate_skills.py --strict` before merge.
Scientific code must be genuinely runnable, not a stub.
