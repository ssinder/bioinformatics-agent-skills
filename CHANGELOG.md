# Changelog

本文件记录本项目的显著变更。格式参考
[Keep a Changelog](https://keepachangelog.com/zh-CN/1.1.0/)，
版本号遵循 [语义化版本](https://semver.org/lang/zh-CN/)。

## [Unreleased]

### 计划中
- 评测套件（`evals/evals.json`）与 `BENCHMARK.md`
- EDAM 本体标注与 bio.tools 注册

## [1.0.0] — 2026-11-22

首个正式发布。

### 新增
- **17 个 Agent Skills**，覆盖转录组、单细胞、结构生物学、分子动力学、
  蛋白对接、多肽设计与文献挖掘，目录见
  [技能清单](https://ssinder.github.io/bioinformatics-agent-skills/skills/)
- **统一命令行入口** `bioagent`（原 `biopipe`；PyPI 上的 `biopipe` 包名
  与控制台脚本已被第三方占用，主命令改名，旧名保留一个版本的兼容别名）
- **双后端交叉验证**：`bulk-rnaseq` 的自研 NumPy 实现可与 Bioconductor
  DESeq2 对照；`single-cell-rnaseq` 的内置实现可与 scanpy 对照
- **纯 Python 降级路径**：缺少 numpy / pandas / scipy / matplotlib 时，
  技能脚本仍可响应 `--help` 并给出明确安装指引
- **可复现溯源文档** `docs/migration-map.md`，记录 217 个研究脚本的
  完整去向与改写理由

### 变更
- 技能描述改为**触发用例前置**的祈使句式（`Use when ...`），
  并保留中英双语触发词；全库 name+description 合计从 8,954 字符
  压缩到 7,478 字符，控制在宿主技能列表共享预算（8,000 字符）之内
- CLI 命令由 `biopipe` 更名为 `bioagent`
- Python 包由 `biopipe` / `biopipe_core` 更名为 `bioagent` / `bioagent_core`

### 修复
- 校验器在缺少 PyYAML 时不再把「依赖缺失」误报为「内容不合规」，
  而是给出明确的安装指引
- 一致性检查不再依赖本地检出目录名，只克隆到任意别名目录都能通过

### 工程
- 7 条 CI 流水线：规范校验、脚本语法、单元测试、无依赖冒烟、
  安全扫描、PR 增量校验、文档站构建
- 发布工作流增加 PyPI Trusted Publishing（OIDC）与 PEP 740 证明
- 新增版本一致性校验：`plugin.json` / `pyproject.toml` / `CITATION.cff`
  三处版本必须相同
- 新增技能列表预算门禁：超出 8,000 字符即失败，防止描述被宿主丢弃

[Unreleased]: https://github.com/ssinder/bioinformatics-agent-skills/compare/v1.0.0...HEAD
[1.0.0]: https://github.com/ssinder/bioinformatics-agent-skills/releases/tag/v1.0.0
