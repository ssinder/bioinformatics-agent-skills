<div align="center">

# BioinfoAgent Skills

**把任意 AI 编程助手变成计算生物学家的科研搭档**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python](https://img.shields.io/badge/Python-3.9%2B-blue.svg)](https://www.python.org/)
[![Skills](https://img.shields.io/badge/Skills-12-brightgreen.svg)](#-available-skills)
[![Agent Skills](https://img.shields.io/badge/Standard-Agent_Skills-blueviolet.svg)](https://agentskills.io/)
[![CI](https://github.com/ssinder/bioinformatics-agent-skills/actions/workflows/skill-spec-validation.yml/badge.svg)](https://github.com/ssinder/bioinformatics-agent-skills/actions/workflows/skill-spec-validation.yml)
[![Tests](https://github.com/ssinder/bioinformatics-agent-skills/actions/workflows/tests.yml/badge.svg)](https://github.com/ssinder/bioinformatics-agent-skills/actions/workflows/tests.yml)
[![Security](https://github.com/ssinder/bioinformatics-agent-skills/actions/workflows/security-scan.yml/badge.svg)](https://github.com/ssinder/bioinformatics-agent-skills/actions/workflows/security-scan.yml)
[![Pages](https://img.shields.io/badge/Docs-GitHub_Pages-blue.svg)](https://ssinder.github.io/bioinformatics-agent-skills/)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo-blue.svg)](https://doi.org/10.5281/zenodo)

📖 **[在线文档 / Documentation](https://ssinder.github.io/bioinformatics-agent-skills/)** ——
快速上手、技能总览、方法学参考与故障排查

[English](#english) | 简体中文

</div>

---

## 这是什么 / What is this

一套面向计算生物学的 **Agent Skills 技能库**。每个技能都是一个标准化的 `SKILL.md` 包,
让 Claude Code、Cursor、Codex 等任意支持 [Agent Skills 开放标准](https://agentskills.io/) 的
AI 助手获得**领域专业知识、可复现流程和验证检查** —— 而不只是会调用 API 的通才。

覆盖转录组学、单细胞、蛋白结构、分子动力学、抗菌肽设计、ML 肽链生成、文献挖掘、NCBI BioProject 元数据与出版级图表。

> 转型自 [`ssinder/Bioinformatics_pipeline`](https://github.com/ssinder/Bioinformatics_pipeline)——
> 217 个真实科研脚本,重构为符合开放标准的可发现、可触发、可复现技能库。

---

## 为什么需要它 / Why

大多数 AI 助手能写 `scanpy` 代码,但不知道:

- 你的单细胞数据里有多少线粒体基因比例该剔除才合理
- DESeq2 的 `padj` 为什么突然全为 NA —— 通常是 `design` 设错了
- RMSF 曲线的第一个尖峰对应哪个柔性 Loop,该拿它去做什么实验
- 检索到的 3000 篇文献里,哪些真正是「靶点发现」而非「功能验证」
- 审稿人要求的 ROC 曲线,坐标轴范围该截到多少才不被质疑

**这些是本项目承载的东西。** 每个技能的 `references/methodology.md` 记录参数选择的科学依据,
`references/troubleshooting.md` 记录常见故障的根因,`SKILL.md` 的 Validation 章节给出结果的自检方法。

---

## 快速开始 / Quick Start

### 方式一:作为 Claude Code 技能安装

```bash
git clone https://github.com/ssinder/bioinformatics-agent-skills.git
cp -r bioinformatics-agent-skills/skills/* ~/.claude/skills/
```

重启 Claude Code 后,技能会被自动发现。当你说「帮我做一次 RNA-seq 差异表达分析」,
它会自行加载 `bulk-rnaseq` 技能。

### 方式二:作为 Python 包使用

```bash
pip install -e .
biopipe --list                      # 列出全部技能
biopipe bulk-rnaseq --help         # 查看某个技能用法
```

### 方式三:直接调用脚本

```bash
export PYTHONPATH="$PWD/packages"

# 蛋白质理化性质计算
python skills/peptide-design/scripts/analyze_peptides.py \
    --input peptides.fasta --output results/peptide_stats.csv

# 绘制 ROC 曲线
python skills/scientific-figures/scripts/make_figures.py roc \
    --input metrics.csv --output figures/roc.png

# 文献检索
python skills/literature-search/scripts/search_literature.py \
    --topic "antimicrobial peptide delivery" --source pubmed --max-results 50
```

---

## 可用技能 / Available Skills

### 🧬 转录组学 / Transcriptomics

| 技能 | 描述 | 关键能力 |
|------|------|----------|
| [`bulk-rnaseq`](skills/bulk-rnaseq/) | RNA-seq 定量与差异表达 | FastQC 质控 · Salmon 定量 · DESeq2 差异分析 · GO/KEGG 富集 |
| [`ena-download`](skills/ena-download/) | ENA 原始数据下载 | Portal API 取 FASTQ · 断点续传 · MD5 校验 · 元数据导出 |
| [`single-cell-rnaseq`](skills/single-cell-rnaseq/) | 单细胞 RNA-seq 标准流程 | 质控过滤 · CP10K 归一化 · HVG · PCA · Leiden 聚类 · UMAP · 细胞注释 |

### 🔬 结构生物学 / Structural Biology

| 技能 | 描述 | 关键能力 |
|------|------|----------|
| [`protein-structure`](skills/protein-structure/) | PDB 结构解析与二级结构 | 自带 PDB 解析器 · DSSP 二级结构 · 几何计算 · UniProt 元数据 |
| [`protein-interactions`](skills/protein-interactions/) | 蛋白质-配体相互作用 | PLIP XML 解析 · 残基级抽取 · 交互网络 · PyMOL 标注 |
| [`molecular-dynamics`](skills/molecular-dynamics/) | 分子动力学模拟与轨迹分析 | RMSD/RMSF/PCA · 氢键 · Ramachandran · 轨迹聚类 · 内置 DCD 解析器 |
| [`alphafold-db`](skills/alphafold-db/) | AlphaFold DB 结构下载 | API 解析下载链接 · 批量下载 · pLDDT 置信度筛选 · 元数据导出 |

### 🧪 蛋白工程 / Protein Engineering

| 技能 | 描述 | 关键能力 |
|------|------|----------|
| [`peptide-design`](skills/peptide-design/) | 抗菌肽设计与筛选 | 理化性质计算 · 疏水性/pI · 稳定性 · MMseqs2 聚类去冗余 |
| [`ml-peptide-generation`](skills/ml-peptide-generation/) | ML 肽链生成与筛选 | PepMLM 掩码语言模型 · 非法氨基酸剔除 · 伪困惑度排序 · 种子可复现 |

### 📚 文献挖掘 / Literature Mining

| 技能 | 描述 | 关键能力 |
|------|------|----------|
| [`literature-search`](skills/literature-search/) | 多源学术文献检索 | PubMed · OpenAlex · Semantic Scholar · DOI 去重 · RIS/BibTeX 导出 |
| [`bioproject-search`](skills/bioproject-search/) | NCBI BioProject 元数据检索 | E-utilities API · 标题/描述/物种 · 关键词筛选 · CSV/JSON 导出 |

### 📊 数据可视化 / Scientific Visualization

| 技能 | 描述 | 关键能力 |
|------|------|----------|
| [`scientific-figures`](skills/scientific-figures/) | 顶刊级图表生成 | ROC · 混淆矩阵 · 热图 · 火山图 · SHAP 摘要 · 统一出版样式 |

---

## 设计原则 / Design Principles

**自包含** —— 核心算法内置纯 Python 实现。scanpy 没装?照样能跑完单细胞流程。
OpenMM 没装?照样能解析 DCD 轨迹并做 RMSD/PCA。依赖可选,能力不打折。

**可复现** —— 固定随机种子、记录参数、锁定版本、输出运行元数据。

**可验证** —— 每个技能都定义了自己的质控环节,结果对不对有明确的检查方法。

**无硬编码** —— 所有路径参数化,所有输出走标准目录约定,没有一条 `/Users/xxx/...` 写死在代码里。

**零门槛** —— `pip install -e .` 即可使用,无需配置云服务、API Key 或数据库。

---

## 项目结构 / Layout

```
bioinformatics-agent-skills/
├── skills/                    # 12 个独立可安装的技能
│   ├── bulk-rnaseq/
│   │   ├── SKILL.md           # 触发描述 + 工作流 + 验证方法
│   │   ├── scripts/           # 可执行脚本
│   │   └── references/         # 方法学原理 + 排错手册
│   └── ...
├── packages/
│   ├── biopipe_core/          # 共享工具层：配置/日志/HTTP重试/CLI/IO
│   └── biopipe/               # 统一命令行入口与技能注册表
├── scripts/
│   ├── validate_skills.py     # 技能规范校验器（CI 核心）
│   └── init_skill.py          # 技能骨架生成器
├── tests/                     # 单元测试与校验器自测
├── docs/                      # 使用指南与技能目录
└── .github/workflows/         # 5 条 CI 流水线
```

---

## 开发 / Development

```bash
git clone https://github.com/ssinder/bioinformatics-agent-skills.git
cd bioinformatics-agent-skills
pip install -e ".[dev]"

python -m pytest tests -v              # 跑测试
python scripts/validate_skills.py      # 校验所有技能
python scripts/validate_skills.py --strict   # 警告也算失败
```

新增技能请阅读 [CONTRIBUTING.md](CONTRIBUTING.md)。

### CI 流水线

| 工作流 | 作用 |
|--------|------|
| `skill-spec-validation` | 校验所有技能符合 Agent Skills 规范 |
| `script-syntax` | 多 Python 版本语法检查 + Ruff 静态检查 |
| `tests` | 单元测试 + 全部技能 `--help` 冒烟测试 |
| `security-scan` | 依赖漏洞审计 + 硬编码凭据扫描 |
| `pr-skill-scan` | PR 中变更技能的定向校验与报告 |

---

## 引用 / Citation

如果你在研究中使用了本项目,请引用:

```bibtex
@software{bioinfoagent2026,
  title  = {BioinfoAgent Skills: Agent Skills for Computational Biology},
  author = {sinder},
  year   = {2026},
  url    = {https://github.com/ssinder/bioinformatics-agent-skills}
}
```

详见 [CITATION.cff](CITATION.cff)。

---

## 贡献 / Contributing

欢迎提交新技能、修复 bug 或改进文档。详见 [CONTRIBUTING.md](CONTRIBUTING.md)。
提交新技能请先读 [docs/SKILL_SPEC.md](docs/SKILL_SPEC.md)。

参与本项目即表示同意 [行为准则](CODE_OF_CONDUCT.md)。

## 许可证 / License

MIT © sinder。详见 [LICENSE](LICENSE)。

## 免责声明 / Disclaimer

本项目提供分析与计算工具,其输出需研究者自行验证,不构成医学、临床或监管建议。

---

<div align="center">

**[给这个项目点个 star ★](https://github.com/ssinder/bioinformatics-agent-skills)** —— 让更多研究者发现这套可复现的科研技能。

</div>

---

# English

**BioinfoAgent Skills** — a collection of Agent Skills for computational biology, built on the
open [Agent Skills standard](https://agentskills.io/). Compatible with Claude Code, Cursor, Codex,
and any host that supports the standard.

Each skill packages a `SKILL.md` trigger description, runnable scripts, and two reference docs:
`methodology.md` (why each parameter is set the way it is) and `troubleshooting.md` (what breaks
and why). Heavy dependencies degrade gracefully — the core algorithms are implemented in pure
Python so the skills work even in minimal environments.

**Domains covered:** bulk RNA-seq, single-cell RNA-seq, protein structure analysis, protein-ligand
interactions, molecular dynamics, antimicrobial peptide design, ML peptide generation, literature
mining, NCBI BioProject metadata retrieval, and publication-quality figures.

```bash
pip install -e .
biopipe --list
```

Contributions welcome — see [CONTRIBUTING.md](CONTRIBUTING.md).

MIT licensed.
