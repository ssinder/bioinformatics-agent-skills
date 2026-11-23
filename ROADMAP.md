# Roadmap

本文件记录**已公开的计划**。它不是承诺书——优先级会随社区反馈调整，
但任何方向的变更都会在此处留下记录。

时间标记用相对季度表述（Q1 = 1–3 月）。当前版本见
[Releases](https://github.com/ssinder/bioinformatics-agent-skills/releases)。

---

## 已完成

- **17 个技能**覆盖转录组、单细胞、结构生物学、分子动力学、蛋白对接、
  多肽设计与文献挖掘
- **纯 Python 核心实现**：不装 scanpy / OpenMM / MDAnalysis 也能跑完主要流程
- **双后端交叉验证**：`bulk-rnaseq` 提供自研 NumPy 实现与 Bioconductor
  DESeq2 双后端，可对照结果；`single-cell-rnaseq` 提供 scanpy 与内置
  numpy 双后端
- **7 条 CI 流水线**：规范校验、脚本语法、单元测试、无依赖冒烟、
  安全扫描、PR 增量校验、文档站构建
- **可复现溯源**：`docs/migration-map.md` 记录 217 个研究脚本的完整去向
- **Agent Skills 规范符合**：全部技能的 frontmatter 符合 agentskills.io 规范

## 进行中

- **评测体系**：按官方 `evals/evals.json` 约定为每个技能建立
  with-skill / without-skill 配对评测，并发布 `BENCHMARK.md`
- **标准与元数据**：EDAM 本体标注、bio.tools 注册、每技能 license 声明

## 计划中

### 近期（下一个季度）

- **跨宿主兼容性矩阵**：把 `skillevaluator harness-check` 纳入 CI，
  覆盖 Claude Code / Codex / OpenCode / Copilot CLI 等宿主
- **BUSCO 与 PhyKIT 技能**：补齐基因组组装质控与系统发育工具链
  （这两个方向在现有生物信息代理评测中被反复标注为缺失）
- **选择式安装指南**：按研究场景给出「每个工作流装 2–3 个技能」的推荐组合

### 中期

- **论文**：面向 AI 编码代理的可复现生物信息学工作流技能库，
  重点写双后端一致性与评测方法
- **跨学科「科学判断」技能**：把现有方法学里的判读约定（背景集选择、
  效应量与显著性之别、图表截断的误导性等）提炼为可跨领域复用的技能
- **容器化参考环境**：为需要重型依赖的技能提供可选的固定版本环境

### 长期

- **工作流编排**：让技能可被组合进 Nextflow / Snakemake 流程
- **上游注册**：Bioconda 配方、WorkflowHub 包装、期刊配套代码可用性声明

---

## 不打算做的事

明确记录下来，避免反复讨论：

- **不做通用 SaaS 或托管服务**。本项目是给研究者的工具集，
  不引入账号体系与云端依赖。
- **不追求技能数量**。广度不等于价值——同等深度下 17 个可用技能
  好过 200 个薄封装。
- **不引入付费墙**。核心能力永久免费且可离线运行。

---

## 参与路线图讨论

路线图的优先级由社区反馈驱动。请到
[Discussions](https://github.com/ssinder/bioinformatics-agent-skills/discussions)
投票或提出调整建议；重大方向变更会先在此处公告再实施。
