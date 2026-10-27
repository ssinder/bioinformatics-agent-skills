---
title: 迁移映射 · 老仓库目录到新技能
description: ssinder/Bioinformatics_pipeline 的 217 个脚本重构为 17 个 Agent Skills 的完整对照表
---

# 迁移映射 · 从老仓库到 Agent Skills

> 老仓库 [`ssinder/Bioinformatics_pipeline`](https://github.com/ssinder/Bioinformatics_pipeline)
> 已被**归档**（read-only）。本表记录它的 217 个脚本、去向与处置理由，
> 供追溯与审计。

## 一句话概览

| | 老仓库 | 新仓库 |
|---|---|---|
| 形态 | 22 个扁平目录，按工具名组织 | 17 个标准 `SKILL.md` 技能包 |
| 脚本 | 217 个 `.py` | 58 个 `.py`（**每个技能都是可独立触发的完整工作流**） |
| 入口 | 人工翻目录找脚本 | `SKILL.md` frontmatter + `biopipe --list` |
| 依赖 | 隐式（`pip install` 后自己摸） | 每个技能独立 `requirements.txt`，核心零第三方依赖 |
| 测试 | 无 | 829 个用例（含反向验证 canary） |

**脚本数从 217 降到 58 不是缩水** —— 减少的是重复与脚手架：
每个技能都补齐了「方法学说明 + 故障排查 + 验证清单 + 测试 + 夹具」，
这些在老仓库里完全不存在。

---

## 完整对照表

### 已改造为技能（18 个老目录 → 17 个技能）

| 老目录 | 脚本 | 行数 | → 新技能 | 关键改造 |
|---|---|---|---|---|
| `rnaseq_python` | 16 | 4026 | `bulk-rnaseq` | DESeq2/edgeR 差异分析；修掉 `padj` 全 NA 的设计矩阵问题 |
| `scrna_python` | 6 | 748 | `single-cell-rnaseq` | 单细胞质控；线粒体比例阈值的合理区间 |
| `rnaseq_download` | 1 | 285 | `ena-download` | ENA 元数据检索与 FASTQ 下载 |
| `crawl_bioproject` | 1 | 338 | `bioproject-search` | NCBI BioProject 元数据抓取 |
| `crawl_database` | 13 | 1974 | `bio-database-fetch` | 通用数据库抓取 |
| `afdb_download` | 2 | 364 | `alphafold-db` | AlphaFold 模型下载；**并入最佳模型提取** |
| `search_pubmed` | 20 | 4519 | `literature-search` | 文献检索 |
| `search_openalex` | 8 | 2063 | `literature-search` | （同上，合并） |
| `paper_extract` | 13 | 2753 | `paper-fulltext-extract` | 全文抽取；补齐 11 个失败测试 |
| `pdb_analyze` | 14 | 2863 | `protein-structure` | 结构解析与分析 |
| `plip_interact` | 10 | 1197 | `protein-interactions` | 相互作用检测 |
| `pymol_show` | 22 | 1486 | `pymol-structure-visualization` | 链间残基接触分析；**几何计算与渲染彻底分离** |
| `dock_lightdock` | 5 | 1413 | `protein-docking-lightdock` | LightDock 对接 |
| `dock_hdock` | 10 | 790 | `protein-docking-hdock` | HADDOCK 对接 |
| `MD_openmm` | 19 | 6881 | `molecular-dynamics` | OpenMM 模拟；平衡协议 + 链分离监测 |
| `peptide_analyze` | 18 | 2595 | `peptide-design` | 肽段性质分析 |
| `design_pepmlm` | 4 | 899 | `ml-peptide-generation` | ML 肽链生成 |
| `data_visualization` | 16 | 3036 | `scientific-figures` | 出版级图表 |

### 未改造（4 个老目录）

| 老目录 | 脚本 | 行数 | 处置 | 理由 |
|---|---|---|---|---|
| `pptx_create` | 9 | 1748 | **归档并导流** | 通用办公工具，不属科研定位；且内容深度绑定单一演示主题，非通用实现。导流至 `workbuddy-pptx` 技能 |
| `feishu_bot` | 10 | 1575 | **安全清理** | 非科研工具。已移除全部明文凭据改用环境变量；README 加安全警告；新增 `.gitignore` |
| `omics_analyze` | 0 | 0 | 空目录 | 无内容 |
| `rnaseq_r` | 0 | 0 | 空目录 | 无内容（Rust/R 版本从未落地） |

---

## 三个"没有一对一迁移"的地方

### 1. `search_pubmed` + `search_openalex` → 一个 `literature-search`

两个目录 28 个脚本做的是同一件事的不同检索后端。拆成两个技能会让
AI 助手在选技能时反复摇摆 —— 它不知道该问「PubMed 检索器」还是
「文献检索器」。

现在 `literature-search` 内部按后端分模块，触发词覆盖两者。

### 2. `afdb_download` 的模型提取被并入 `alphafold-db`

原目录里「下载」和「提取最佳模型」是分开的两步脚本。拆开意味着
用户得记住顺序 —— 而这是个只有一种合理顺序的流程。

### 3. `plip_interact` 的相互作用检测被独立成 `protein-interactions`

它与 `pdb_analyze` 都处理结构，但**触发场景完全不同**：
前者回答「这两个残基之间是什么作用力」，后者回答「这个结构长什么样、
有哪些二级结构」。合在一起会让技能描述失焦。

---

## 改造中修掉的真实缺陷

不是「搬过来」，是**读懂之后重写**。以下是各技能在改造过程中发现并修复的
实质性缺陷（不是风格问题，是会产出错误科学结论的问题）：

| 技能 | 老实现的缺陷 | 后果 |
|---|---|---|
| `pymol-structure-visualization` | `cmd.select(f"{a} within {cutoff} of {b}")` 选出的是**各自靠近**的残基再做笛卡尔积 | 两个关键原子相距 15 Å 也会被报成"接触" |
| `pymol-structure-visualization` | ASP/GLU 被列入氢键供体表 | 羧基在生理 pH 下无氢，系统性虚增氢键数 |
| `pymol-structure-visualization` | 全部 HIS 按正电计 | 咪唑基质子化比例仅约 4%，虚增约 25 倍 HIS 盐桥 |
| `pymol-structure-visualization` | 配体（LIG/HEM）硬套残基分类表 | 苯环配体被报成"疏水残基"、带电配体被报成"负电残基" |
| `bulk-rnaseq` | `design` 设错导致 `padj` 全 NA | 静默产出无意义的差异分析结果 |
| `molecular-dynamics` | 平衡协议不完整 | 未做能量最小化即进入生产模拟 |

每个技能都有 `references/methodology.md` 记录完整的方法学取舍，
`references/troubleshooting.md` 记录常见症状与排查路径。

---

## 老仓库的安全处置

⚠️ **飞书 App Secret 已泄露且无法通过代码修复。**

- 泄露位置：老仓库 `feishu_bot/` 的 10 个文件，提交 `a3e456c`
- 该仓库是**公开仓库**，历史中的凭据**依然可被检出**
- 老仓库现已归档（read-only），代码侧清除**无法推送**
- **唯一有效的补救：在飞书开放平台重置 App Secret**（重置后旧值立即失效）

代码侧的清理已完成并保留在本地分支 `security/credential-cleanup`：
10 个文件全部改为环境变量读取、客户端延迟创建、缺失时明确报错、
新增 `.gitignore` 忽略 `.env` / `*.secret`。

新仓库的 `scripts/ci_security_scan.py` 是这次事故催生的产物 ——
它原先用 `\b(secret)\b` 匹配，而下划线与字母**同属单词字符**，
导致 `APP_SECRET` 完全逃过扫描。现已修复并有 32 个回归测试锁定，
其中包含**反向验证 canary**：把正则改回旧写法，用例立刻失败。
