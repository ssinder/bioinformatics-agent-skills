# 文献挖掘方法学原理与参数选择

本文说明 `literature-search` 技能背后每个算法环节的原理，以及在什么场景下应该选什么参数。

---

## 1. 为什么需要多源检索

单一文献库都有系统性盲区：

| 数据源 | 覆盖优势 | 覆盖劣势 | 引用计数口径 |
| --- | --- | --- | --- |
| PubMed | 生物医学文献最全，含 MeSH 主题词标引 | 覆盖率仅生物医学；会议论文少 | `pmcrefcount`（PubMed 内部引用统计） |
| OpenAlex | 跨学科覆盖最广（预印本、书籍、学位论文、数据集） | 会议论文与中文学位论文覆盖弱 | `cited_by_count`（Crossref + 快照计数） |
| Semantic Scholar | 引文网络最完整，含 influential citation 标记 | 冷门文献收录慢 | `citationCount`（含自引） |

三者合起来才能同时覆盖"经典生物医学证据"和"跨学科最新进展"。

---

## 2. 检索式构造

### 2.1 topic 与 keywords 的关系

- `--topic` 是**主检索式**，原样送入各源搜索引擎（PubMed 走 MeSH + 自由词自动映射，OpenAlex 走全文相关性，Semantic Scholar 走语义检索）。
- `--keywords` 是**收窄条件**，与 topic 做 AND，同时**参与相关性打分**。

因此有两种典型用法：

```
# 用法 A：广撒网 + 事后打分（推荐用于探索性检索）
--topic "antifungal peptide" --keywords "Sclerotinia" "mechanism"
→ 得到大量命中，再靠 relevance 排序收敛

# 用法 B：窄检索（推荐用于系统综述的精确检索阶段）
--topic "antifungal peptide AND mechanism" --keywords "Sclerotinia"
→ 直接在源头收窄，漏检风险低但召回也低
```

### 2.2 检索词应覆盖概念空间

用 `--keywords` 补齐同义词、缩写、旧名与机制词，比加长 `--topic` 更有效：

```
--topic "plant defensin"
--keywords "PR-1" "γ-core" "membrane permeabilization" "fungal"
```

### 2.3 年份过滤

`--date-range` 会被翻译成各源原生语法：

| 数据源 | 翻译结果 |
| --- | --- |
| PubMed | `AND ("2018"[PDAT] : "2024"[PDAT])`（同时把 PDAT 与出版年一起约束） |
| OpenAlex | `filter=from_publication_date:2018-01-01,to_publication_date:2024-12-31` |
| Semantic Scholar | `year=2018-2024` |

注意 PubMed 用的是 `[PDAT]`（Publication Date），比 `[DP]`（Date of Publication）更宽松，能覆盖电子预发表。

---

## 3. 相关性评分公式

```
relevance = Σ_keywords [ 3.0·(命中标题) + 1.0·(命中摘要) ]
          + min(2.0, 0.6·citations^0.25)      # 引用数加成，对数饱和
          + 0.5  if year >= 2021              # 近期加成
```

设计要点：

- **标题权重 3 倍于摘要**：标题中出现关键词是作者主动声明主题相关性的强信号。
- **引用数用四次方根而非线性**：1000 次引用不应等于 10 次引用的 100 倍权重。`0.6·c^0.25` 使 10 次引用得 1.13 分、100 次得 2.0 分（封顶），既奖励影响力又不会让高引老论文垄断榜首。
- **近期加成是固定 0.5 分**：只作为打破平局的次要因素，权重刻意保持很小，避免把领域内被低估的新工作系统性压掉。

**阈值选择建议**：先跑 `--min-relevance 0` 导出全部，人工看 `relevance` 分布直方图，选取 60–70 分位处的值作为阈值。经验起点：

| 目标召回量 | 建议 `--min-relevance` |
| --- | --- |
| 精确探索（约 30 条） | 5.0 以上 |
| 常规综述（约 100–200 条） | 3.0–4.0 |
| 宽口径背景材料 | 1.5–2.5 |

---

## 4. 去重算法

### 4.1 两级策略

```
第一级：DOI 精确匹配
    normalize_doi(): 去 https://doi.org/ 前缀 → 转小写 → 必须以 "10." 开头
    命中即判定重复，并合并信息

第二级：标题模糊匹配
    normalize_title(): 转小写 → 非字母数字转空格 → 压缩空白
    与已保留标题做 difflib.SequenceMatcher 比对
    ratio >= title_similarity (默认 0.92) 判定重复
    ⚠️ 仅当两条记录中至少有一条缺少 DOI 时才生效；
       双方都有 DOI 且不同 → 直接视为不同文献，不做标题比较
```

### 4.2 为什么 DOI 优先

同一篇论文在不同源的 DOI 表现完全一致，是最强标识。标题则可能因为大小写、希腊字母转写（β vs Beta）、副标题增删而不同。

### 4.2.1 DOI 冲突保护（重要）

学术出版中，**不同期刊会刊出完全同名的文章**。实测检索 "antimicrobial peptide" 时，
OpenAlex 返回的前 5 条里有 3 条标题都是 "Antimicrobial Peptides" / "Antimicrobial peptides"，
但 DOI 分别为 `10.3390/ph6121543`、`10.1016/j.cub.2015.11.017`、`10.1016/j.jaad.2004.08.026`
—— 它们是三篇不同论文。

若不加保护，标题相似度会把这些文章全部合并成 1 条，造成严重的文献丢失
（实测 5 条被合并为 1 条，丢失率 80%）。

因此规则为：**双方都有 DOI 且 DOI 不同时，直接判定为不同文献**。
标题模糊匹配只在至少一方缺 DOI 时才作为补充手段。

### 4.3 阈值 0.92 的含义

`difflib.SequenceMatcher` 的 ratio 是"匹配字符数 × 2 / 总字符数"：

- ratio = 1.0 → 完全相同
- ratio = 0.92 → 约 15–20 个字符的差异（对 100 字符标题而言）
- ratio < 0.85 → 通常已经是不同论文

短标题（如 "Analysis of p53" 17 字符）容忍度更高：3 个字符差异就可能达到 0.88。因此对短标题集合，可以把阈值提高到 0.95。

### 4.4 合并而非丢弃

被判定重复的记录不会直接扔掉，其信息会合并回保留记录：

- 补充缺失摘要（例如 PubMed 命中但没摘要，OpenAlex 命中有摘要）
- 引用数取两者最大值
- `source` 列标记为 `PubMed+OpenAlex`，这个信息本身有科研价值——多源命中提示该文献影响力较高

---

## 5. 摘要富集策略

摘要缺失是文献挖掘中最常见的质量缺口。三个源的摘要来源互补：

| 顺序 | 来源 | 手段 | 特点 |
| --- | --- | --- | --- |
| 1 | PubMed | `efetch.fcgi`（`rettype=abstract, retmode=xml`） | 摘要质量最高（作者提交的正式摘要，含分段标签） |
| 2 | OpenAlex | `GET /works/https://doi.org/{doi}` | 摘要以**倒排索引**存储，需重建 |
| 3 | Semantic Scholar | `GET /paper/DOI:{doi}`（`fields=abstract,tldr`） | 有摘要最好，没有则退化为 `tldr` 一句话摘要 |

### 5.1 OpenAlex 倒排索引重建

OpenAlex 不直接返回摘要文本，而是返回 `{"word": [pos0, pos1, ...]}`：

```
输入: {"Helicobacter": [0], "pylori": [1], "infection": [2, 5]}
输出: "Helicobacter pylori infection"      （位置 3、4 缺词，用空串占位）
```

重建算法：建立 `位置 -> 词` 映射，取最大位置 N，遍历 0..N 拼接。

### 5.2 摘要长度阈值

`--min-abstract-for-enrich` 默认 60 字符。低于此长度视为"无效摘要"（可能只是期刊名或一句标题）。设为 0 表示任何非空摘要都算数。

---

## 6. 引用格式规范

### 6.1 GB/T 7714-2015 顺序编码制

```
张三, 李四, 王五, 等. 题名[J]. 刊名, 2020, 12(3): 45-67. DOI: 10.1000/xxx.
```

- 作者 3 人以内全列，超过 3 人取前 3 + "等"
- 期刊论文文献类型标识为 `[J]`
- 顺序编码制要求正文用 `[1] [2]` 角标标注，`export_citations.py --format plain` 直接生成 `[1] ...` 形式

### 6.2 BibTeX

键名规则：首位作者姓（小写）+ 年份 + 标题首个英文词 + 序号，例如 `zhang2020antimicrobial1`。这个规则保证同一篇文献重复导出时键名稳定，不会因排序变化而产生大量重复条目。

### 6.3 RIS

RIS 的 `AU` 字段格式为 `姓, 名`。若数据源只给了全名（`John Smith`），脚本按最后一个空格拆分：`Smith, John`。数据源只给全名而无法区分姓与名时，这是最常见的折中方案。

---

## 7. 速率限制与礼貌池

| 数据源 | 无 Key 限制 | 有 Key 限制 | 本技能默认间隔 |
| --- | --- | --- | --- |
| NCBI E-utilities | 3 req/s | 10 req/s | 0.34 s |
| OpenAlex | 10 req/s（礼貌池建议提供 mailto） | 同 | 0.34 s |
| Semantic Scholar | 1 req/s（未认证） | 提高 | 1.1 s（无 Key 时） |

共享层 `http_session()` 已内置：
- **指数退避重试**：429/500/502/503/504/408 自动重试，退避 `2^attempt + U(0, 0.5)` 秒
- **请求节流**：包装 `session.request`，保证两次请求间至少间隔 `min_interval`

设置环境变量可显著提速：

```bash
export NCBI_API_KEY=xxxx     # NCBI 10 req/s
export S2_API_KEY=xxxx       # Semantic Scholar 大幅提速
```

---

## 8. 推荐的检索工作流

```
阶段 1 探索：--topic 单主题，三源全开，不设过滤 → 看分布
阶段 2 收窄：加入 --keywords 与 --date-range，设 --min-citations 10–20
阶段 3 富化：--enrich-abstracts，让摘要完整度 > 90%
阶段 4 精筛：--min-relevance + --include-terms/--exclude-terms
阶段 5 引文：导出 BibTeX + RIS，人工核查前 50 条
```

每一步都保留中间 CSV，便于回溯"某条文献是在哪一步被淘汰的"。
