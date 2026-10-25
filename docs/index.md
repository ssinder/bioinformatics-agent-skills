---
layout: default
title: BioinfoAgent Skills
---

# BioinfoAgent Skills

面向计算生物学的**标准 Agent Skills 工具库** —— 遵循 [Agent Skills 开放标准](https://agentskills.io/)，
把 RNA-seq、单细胞、结构生物学、分子动力学、多肽设计与文献挖掘等分析流程，
封装成任何 AI Agent 都能即插即用的技能。

<div class="hero">
  <p class="tagline">把可复现的科研分析流程，变成 Agent 可以直接调用的技能</p>
</div>

## 快速开始

```bash
git clone https://github.com/ssinder/bioinformatics-agent-skills.git
cp -r bioinformatics-agent-skills/skills/* ~/.claude/skills/
```

重启宿主后技能即被自动发现。完整说明见 [快速上手](getting-started.md)。

## 为什么是 Agent Skills

| 特性 | 说明 |
|---|---|
| **开放标准** | `SKILL.md` + `scripts/` + `references/`，遵循 [agentskills.io](https://agentskills.io/) 规范，兼容 Claude Code、Cursor、Codex 等宿主 |
| **优雅降级** | `--help` 属于元信息查询，不会因缺少 numpy / pandas 等科学计算栈而崩溃；真正用到该功能时才提示安装 |
| **零重型依赖** | 核心安装仅需 `requests` + `pyyaml`，其余按需通过 extras 安装 |
| **跨版本兼容** | CI 覆盖 Python 3.9 / 3.11 / 3.13 |
| **质量保障** | 6 条流水线：规范校验、语法与 lint、单测、依赖与凭据安全扫描、PR 技能扫描、发布 |

## 技能总览

{% include_relative skills/index.md %}

## 与普通脚本仓库的区别

普通脚本仓库要求使用者自己读文档、拼命令、判断参数。
Agent Skills 把这些知识**写进了技能描述里** —— Agent 读一眼 `SKILL.md`
就知道该用哪个技能、怎么调、要注意什么。

```yaml
---
name: bulk-rnaseq
description: >
  This skill should be used when running bulk RNA-seq differential
  expression analysis: FastQC/trimming, Salmon quantification, counts
  matrix construction, DESeq2 statistics, GO/KEGG enrichment and
  publication figures. 用于 bulk 转录组（RNA-seq）差异表达分析全流程。
---
```

这段描述就是 Agent 的路由依据。

## 文档

- [快速上手](getting-started.md) —— 安装、配置、首次使用
- [技能总览](skills/index.md) —— 全部技能及其适用场景
- [技能规范](SKILL_SPEC.md) —— 如何编写一个新技能
- [贡献指南](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/CONTRIBUTING.md)

## 引用

如果你在研究中使用了本项目，请引用：

```bibtex
@misc{bioinfoagentskills2026,
  title        = {BioinfoAgent Skills: An Open Standard Agent Skills Library for Computational Biology},
  author       = {ssinder},
  year         = {2026},
  howpublished = {\url{https://github.com/ssinder/bioinformatics-agent-skills}},
  note         = {Accessed: 2026-10-25}
}
```

## 许可证

MIT，见 [LICENSE](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/LICENSE)。
