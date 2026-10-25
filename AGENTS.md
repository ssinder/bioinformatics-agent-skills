# AGENTS.md

本文件为在本仓库中工作的 AI 编码助手提供上下文与约定。

## 这是什么项目

`bioinformatics-agent-skills` 是一个 **Agent Skills 技能库**,为计算生物学领域提供
可被 AI 助手自动发现与触发的标准化技能包。每个技能是一个独立目录,含 `SKILL.md`、
`scripts/`、`references/`。

## 目录约定

```
skills/<skill-name>/       技能包,一个技能一个目录,kebab-case 命名
packages/biopipe_core/     共享工具层（配置/日志/HTTP/CLI/IO）
packages/biopipe/          统一命令行入口与技能注册表
scripts/                   仓库级工具（校验器、骨架生成器）
tests/                     pytest 测试
docs/                      文档与技能目录
```

## 开发约定

### 1. 新增技能

必须走这个流程,不要手工造目录:

```bash
python scripts/init_skill.py <name> --summary "..." --data-desc "..."
```

然后填写真实的 SKILL.md 与脚本,最后校验:

```bash
python scripts/validate_skills.py --skill <name> --strict
```

### 2. 脚本规范

- 共享能力一律从 `biopipe_core` 导入,不要重复实现路径解析、日志、HTTP 重试
- 禁止硬编码绝对路径;输入输出路径来自 `argparse`
- 禁止裸 `print()`,用 `get_logger("biopipe.<skill>")`
- 入口统一用 `bootstrap_cli(description, skill_name, main)`
- 重量级依赖(scanpy/OpenMM/MDAnalysis 等)必须优雅降级,核心算法内置纯 Python 实现
- 每个脚本都要能被 `--help` 正常响应,且不依赖真实数据

### 3. SKILL.md 规范

- `name` 必须与目录名一致,kebab-case
- `description` 第三人称,含触发语句,中英文关键词都要有
- 必须包含 Overview / When to Use / Quick Start / Workflow / Parameters / Outputs / Validation / References 八个章节

完整规范见 [docs/SKILL_SPEC.md](docs/SKILL_SPEC.md)。

### 4. 注册新技能

新增技能后,必须同步更新 `packages/biopipe/registry.py` 的 `REGISTRY`,
否则 `biopipe <name>` 命令不可用,CI 的集成测试会失败。

## 提交前自检

```bash
python -m pytest tests -q                      # 测试
python scripts/validate_skills.py --strict      # 规范校验
python -m compileall -q skills packages scripts  # 语法
```

三项全绿再提交。

## 科研正确性底线

这是科研工具,代码必须**真实可运行**:

- 不要写 `pass  # TODO` 冒充实现
- 不要返回硬编码的假数据
- 算法可以简化,但要说明简化了什么、以及何时不适用
- 结果必须能被 `SKILL.md` 的 Validation 章节检查

## 提交信息规范

```
feat: add <skill-name> for <purpose>
fix: <what was broken>
docs: <what changed>
refactor: <what moved>
ci: <what changed>
```
