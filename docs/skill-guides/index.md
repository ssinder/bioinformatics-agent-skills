---
layout: default
title: 工作流指南
---

# 工作流指南 / Skill Guides

每个技能都是独立可用的。但真实的科研流程往往需要**两三个技能串联**，
而一次性装下全部 17 个技能反而会拖慢触发准确率——宿主需要为每个技能
保留一段描述，多了以后必然超预算，被丢弃的那些技能就再也匹配不上了。

本页面给出**按研究场景的推荐组合**。每个组合控制在 2–3 个技能。

## 怎么用这些组合

1. 只安装你当前项目需要的技能
2. 组合内的技能可以按顺序串联，也可以按需单独调用
3. 换研究主题时再换组合，不必长期挂着全部技能

```bash
# 例：只装转录组基础流程的三个技能
npx skills add ssinder/bioinformatics-agent-skills \
  --skill ena-download --skill bulk-rnaseq --skill scientific-figures
```

---

## 组合一览

| 场景 | 技能组合 | 说明 |
|---|---|---|
| [从原始数据到差异表达](#1-从原始数据到差异表达) | `ena-download` → `bulk-rnaseq` → `scientific-figures` | 最常见的 bulk RNA-seq 全流程 |
| [单细胞图谱构建](#2-单细胞图谱构建) | `single-cell-rnaseq` → `scientific-figures` | 质控、聚类、注释与出图 |
| [蛋白结构解析与相互作用](#3-蛋白结构解析与相互作用) | `protein-structure` → `protein-interactions` → `pymol-structure-visualization` | 从 PDB 到相互作用网络与可视化 |
| [对接与构象分析](#4-对接与构象分析) | `protein-docking-hdock` → `molecular-dynamics` | 远程对接 + 本地动力学验证 |
| [多肽设计与筛选](#5-多肽设计与筛选) | `peptide-design` → `ml-peptide-generation` | 性质计算 → 机器学习生成 |
| [文献调研闭环](#6-文献调研闭环) | `literature-search` → `paper-fulltext-extract` | 检索去重 → 全文抽取 |
| [化合物与物种信息查询](#7-化合物与物种信息查询) | `bio-database-fetch` | 单技能即可 |

---

## 1. 从原始数据到差异表达

**场景**：拿到一个 BioProject 号，想复现或重新分析它的表达差异。

```bash
npx skills add ssinder/bioinformatics-agent-skills \
  --skill ena-download --skill bulk-rnaseq --skill scientific-figures
```

**串联方式**：

1. `ena-download` 按 SRA/ERR/DRR 号取 FASTQ，同时导出 run 级元数据
2. `bulk-rnaseq` 做质控 → 定量 → 差异表达 → 富集分析
3. `scientific-figures` 出火山图、MA 图与富集图

**关键检查点**：`bulk-rnaseq` 的 Validation 章节要求比对 `padj` 是否全为
`NA`——这几乎总是 `design` 公式写错导致的，不是数据问题。

## 2. 单细胞图谱构建

**场景**：拿到 10x 的矩阵，要做质控、聚类和细胞类型注释。

```bash
npx skills add ssinder/bioinformatics-agent-skills \
  --skill single-cell-rnaseq --skill scientific-figures
```

**注意**：`single-cell-rnaseq` 自带 `--backend numpy` 的纯 Python 降级
实现，没有 scanpy 也能产出结果。装了 scanpy 时可用 `--backend scanpy`
走官方实现，两者结果应大体一致。

**关键检查点**：质控过滤率通常在 5%–30%，超过 50% 说明阈值过严。

## 3. 蛋白结构解析与相互作用

**场景**：解析一个复合物结构，找出结合界面的关键残基。

```bash
npx skills add ssinder/bioinformatics-agent-skills \
  --skill protein-structure --skill protein-interactions \
  --skill pymol-structure-visualization
```

**串联方式**：

1. `protein-structure` 解析 PDB、指派二级结构、算 Rg 与形状参数
2. `protein-interactions` 从 PLIP XML 提取残基级相互作用，排序热点残基
3. `pymol-structure-visualization` 计算链间接触几何，输出 PyMOL 脚本

**关键检查点**：`pymol-structure-visualization` 的输出包含每类相互作用的
残基对清单，用于人工复核自动分类是否合理。

## 4. 对接与构象分析

**场景**：预测两个蛋白的复合物结构，并检查构象是否稳定。

```bash
npx skills add ssinder/bioinformatics-agent-skills \
  --skill protein-docking-hdock --skill molecular-dynamics
```

`protein-docking-hdock` 把作业提交到 HDOCK 服务器并记录到 CSV，
可稍后取回结果；`molecular-dynamics` 在本地做轨迹的 RMSD/RMSF/PCA 分析。

**关键检查点**：HDOCK 的 Docking Score 越低越好，但单个分数不能说明
构象可靠——需要用动力学分析补充判断。

## 5. 多肽设计与筛选

**场景**：为某个靶点设计候选多肽并排序。

```bash
npx skills add ssinder/bioinformatics-agent-skills \
  --skill peptide-design --skill ml-peptide-generation
```

`peptide-design` 计算理化性质、按 GO 术语与序列特征判别抗菌性、用 MMseqs2
聚类去冗余；`ml-peptide-generation` 用 PepMLM 生成候选并按伪困惑度排序。

**关键检查点**：生成时固定随机种子，否则两次运行结果不可复现。

## 6. 文献调研闭环

**场景**：为某个课题做系统性的文献检索并精读全文。

```bash
npx skills add ssinder/bioinformatics-agent-skills \
  --skill literature-search --skill paper-fulltext-extract
```

`literature-search` 跨 PubMed / OpenAlex / Semantic Scholar 检索去重；
`paper-fulltext-extract` 把 PDF 转成结构化全文并建语料索引。

**注意**：`paper-fulltext-extract` 需要一个本地 GROBID 服务，
请先阅读该技能的 [Quick Start](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills/paper-fulltext-extract/SKILL.md)。

## 7. 化合物与物种信息查询

```bash
npx skills add ssinder/bioinformatics-agent-skills --skill bio-database-fetch
```

按 PubChem CID 查理化性质与 CAS 号，按 NCBI Taxonomy ID 查物种谱系。

**注意**：NCBI 接口需要项目方注册的 tool/email 身份，脚本已内置；
个人 API key 是可选的，仅用于把速率上限从 3 req/s 提到 10 req/s。

---

## 不确定装哪些？

从 [完整技能清单](https://ssinder.github.io/bioinformatics-agent-skills/skills/)
逐个查看，每个技能的说明里都写了适用与不适用场景。
