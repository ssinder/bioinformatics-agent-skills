# 文献全文抽取方法学

本文档说明 `paper-fulltext-extract` 技能背后的技术原理：TEI 是怎么组织的、GROBID 怎么调用、章节树与首句怎么提取，以及从老仓库 `paper_extract/` 移植时修掉的 5 个缺陷各自的成因与判据。

---

## 1. 为什么需要 GROBID，而不是直接抽 PDF 文本

PDF 存的是**渲染指令**，不是语义结构。直接用 `pdfplumber` / `PyPDF2` 抽出来的是一堆按阅读顺序排列的字符流：章节标题和正文段落混在一起，图表说明与正文交织，参考文献是一段没有结构的连续文本。

科研文献挖掘需要的是**结构化**信息：论文标题、作者、摘要、关键词、章节层级、每张图的图注、每条参考文献的 DOI。这些信息在 PDF 的视觉排版里是"看得见"的（字体变大 = 标题，加粗 = 小标题），但不在字符流里。

GROBID 做的就是这件事：它用版面分析模型重建 PDF 的视觉结构，输出 TEI XML。

**代价**：GROBID 是 Java 服务，需要单独部署，不能像纯 Python 库那样 `pip install`。

---

## 2. GROBID 部署

### 2.1 Docker（推荐）

```bash
docker run -it -p 8070:8070 kermitt/grobid
```

首次启动需要下载模型，耗时较长（数分钟到十几分钟，取决于网络）。等到日志出现：

```
Server started
```

才表示可以接受请求。**在此之前任何请求都会连接失败或超时。**

### 2.2 关键端点

| 端点 | 方法 | 用途 |
| --- | --- | --- |
| `/api/isalive` | GET | 健康检查，返回 `true` |
| `/api/version` | GET | 版本号 |
| `/api/processFulltextDocument` | POST | 全文 PDF → TEI XML |

### 2.3 请求格式

`multipart/form-data`，**文件字段名必须是 `input`**：

```
--boundary
Content-Disposition: form-data; name="consolidateHeader"

true
--boundary
Content-Disposition: form-data; name="flavor"

precision
--boundary
Content-Disposition: form-data; name="input"; filename="paper.pdf"
Content-Type: application/pdf

<PDF 二进制>
--boundary--
```

### 2.4 布尔参数必须是小写字符串

这是 GROBID 的一个经典坑。传 Python 的 `bool` 会被编码成 `"True"` / `"False"`，而 GROBID 的参数解析按**非空即真**处理 —— `"False"` 也是非空字符串，于是被当成 `true`。

结果：想关掉某个选项，参数传下去却**依然是开的**，而且不报任何错。

```python
# 错
params = {"consolidateHeader": str(True)}          # → "True"
# 对
params = {"consolidateHeader": str(True).lower()}   # → "true"
```

本技能的实现（`grobid_client.GrobidClient.build_params`）统一用 `str(bool(x)).lower()`。

### 2.5 consolidateHeader 为什么重要

`consolidateHeader` 控制 GROBID 是否合并参考文献条目里的作者与机构信息。关掉它，参考文献里的作者会被拆成多个 `<persName>`；开着，跨条目的同一作者会被归并。

**标题质量问题多半出在这里没开** —— 很多 PDF 的页眉（期刊名 + 卷期 + 页码）会被 GROBID 误认为标题，`consolidateHeader` 开启时 GROBID 会用更结构化的线索纠正一部分。

### 2.6 flavor

| flavor | 含义 | 权衡 |
| --- | --- | --- |
| `precision` | 高精度模型 | 慢（单篇 10-60 秒），版面还原准 |
| `speed` | 轻量模型 | 快（单篇 1-5 秒），复杂排版会丢结构 |

批量筛读用 `speed`，需要引用与图注准确用 `precision`。

### 2.7 重试策略

GROBID 在批量处理时经常返回 **503（过载）** 或 **500**。这类错误**恰恰应该重试**。

本技能的重试矩阵：

| 错误 | 行为 | 理由 |
| --- | --- | --- |
| 网络异常（URLError） | 重试 | 瞬时问题 |
| 429 / 502 / 503 / 504 | 重试 | 服务端过载或临时不可用 |
| 4xx（除 429） | **不重试** | 客户端错误，重试也没意义 |
| 其他异常 | 不重试 | 未知错误，重试可能放大问题 |

退避是指数的：第 1 次等 `retry_delay`，第 2 次等 `2×retry_delay`，以此类推。

老仓库的实现只在**超时**时重试，5xx 直接当失败 —— 批量处理时遇到一次 503 就整篇丢掉了。

### 2.8 断点续跑

批量处理时，`process_directory` 默认 `skip_existing=True`：已存在 `<name>.tei.xml` 的直接跳过。

这让"处理 500 篇 PDF，中途 Ctrl-C"变成一个安全的操作 —— 重跑会从断点继续，而不是从头再来一遍。

---

## 3. TEI XML 结构

TEI（Text Encoding Initiative）命名空间：`http://www.tei-c.org/ns/1.0`。

### 3.1 整体骨架

```xml
<TEI xmlns="http://www.tei-c.org/ns/1.0" xml:id="doc">
  <teiHeader>              ← 元数据
    <fileDesc>
      <titleStmt>
        <title type="main">论文标题</title>
        <author><persName>…</persName><affiliation>…</affiliation></author>
      </titleStmt>
      <sourceDesc>
        <biblStruct>        ← 本文自己的元数据（期刊名、DOI），**不是**参考文献
          <analytic><idno type="DOI">…</idno></analytic>
          <monogr><title>期刊名</title></monogr>
        </biblStruct>
      </sourceDesc>
    </fileDesc>
    <profileDesc>
      <abstract><div><p>…</p></div></abstract>
      <textClass><keywords><term>…</term></keywords></textClass>
    </profileDesc>
  </teiHeader>

  <text>
    <body>                 ← 正文
      <div><head>Introduction</head><p>…</p></div>
    </body>
    <back>                 ← 参考文献
      <div type="references"><listBibl><biblStruct>…</biblStruct></listBibl></div>
    </back>
</TEI>
```

### 3.2 `<title>` 元素的多义性

TEI 里有**至少四种** `<title>`：

| 位置 | 含义 | 能否作论文标题 |
| --- | --- | --- |
| `<titleStmt><title type="main">` | 论文标题 | ✅ 首选 |
| `<titleStmt><title>`（无 type） | 论文副标题等 | ⚠️ 可作候选 |
| `<monogr><title>` | **期刊名** | ❌ |
| `<biblStruct><analytic><title>` | 文献条目标题 | ❌ |

老仓库的回退逻辑是"遍历所有 `<title>` 元素"。在标题被 GROBID 误填成章节名时，它会回退到**期刊名** —— 比回退到文件名更糟，因为期刊名在检索时毫无区分度。

本技能改为：候选只从 `<titleStmt>` 下取（`tei_parser._candidate_titles`）。

### 3.3 `xml:id` 的真实键名

`xml` 是 XML 规范保留的前缀。ElementTree 会把它**展开**成带命名空间的完整键：

```python
elem.get("xml:id")              # → 永远 None（静默失败）
elem.get("{http://www.w3.org/XML/1998/namespace}id")   # → "fig_1"
```

这是个非常隐蔽的坑 —— 不报错，只是**永远取不到值**。后果是图表 id、参考文献 id 全部变成空字符串。

本技能定义了 `tei_parser.XML_ID_ATTR` 常量统一处理。

---

## 4. 章节树提取

### 4.1 递归解析 `<div>`

GROBID 的正文结构是 `<body><div><head>…</head><p>…</p></div></body>`，`<div>` 可以嵌套：

```xml
<body>
  <div><head>Methods</head><p>…</p>
    <div><head>Data collection</head><p>…</p></div>
    <div><head>Model training</head><p>…</p></div>
  </div>
</body>
```

解析是递归的，`level` 随深度递增。

### 4.2 兜底分支的陷阱

少数 TEI 变体不用 `<div>`，而是把章节平铺在 `<body>` 下（`<head>` 与 `<p>` 互为兄弟）。这种情况需要一个兜底分支。

**但兜底必须只在顶层 `<body>` 上做。**

老仓库的实现里，兜底分支写在递归函数 `_parse_divs` 的末尾 —— 于是它在**每一层**都会触发。后果是每个 `<div>` 自己的 `<head>` 被重复解析一次：

```
Introduction
Methods
Data collection
Methods        ← 重复
Data collection ← 重复
Model training  ← 重复
Results
Conclusion
Results          ← 重复
Conclusion       ← 重复
```

11 个章节 vs 实际 6 个。

修法：把兜底逻辑上移到 `extract_sections`，只处理顶层 `<body>` 的平铺结构。

### 4.3 章节名归一化

PDF 抽取常把章节名弄出异常空格：`"INTRODUC TI ON"`、`"DISCUSS I ON"`、`"MATER I ALS AND METHODS"`。

成因是渲染层的**字间距被当成了空格** —— 排版上字母之间本来就有 kerning，PDF 抽取时如果对字距敏感就会切开。

**修复靠字典而不是猜正则**：把标题里的空格全部去掉，看是否正好等于某个已知章节名。

| 输入 | 去空格 | 是否命中字典 | 结果 |
| --- | --- | --- | --- |
| `"INTRODUC TI ON"` | `introduction` | ✅ | `introduction` |
| `"DISCUSS I ON"` | `discussion` | ✅ | `discussion` |
| `"Materials  and   Methods"` | `materialsandmethods` | ❌ | 保持原样，只压缩空白 |

第三行说明为什么不能用"无脑去空格"：多词章节名去空格后不是任何单词。必须靠字典判断。

---

## 5. 标题有效性判断

### 5.1 GROBID 的标题误填

某些 PDF（尤其是页眉里有大字号的排版）会让 GROBID 把**章节名**填进 `<title type="main">`：

```xml
<title type="main">INTRODUCTION</title>   ← 误填
```

### 5.2 判据：精确相等，不是包含

老仓库的判断：

```python
if any(marker.lower() in title.lower() for marker in ["introduction", "methods", "abstract"]):
    # 认为标题无效
```

**"包含"匹配会误杀正常标题**：

| 正常论文标题 | 含 "abstract"？ | 老代码判定 | 正确判定 |
| --- | --- | --- | --- |
| `"Abstract Methods for Targeting X"` | 是 | ❌ 误杀 | ✅ 有效 |
| `"Introduction to Materials Chemistry"` | 是 | ❌ 误杀 | ✅ 有效 |
| `"Results of a Systematic Review"` | 是 | ❌ 误杀 | ✅ 有效 |
| `"INTRODUCTION"` | 是 | ✅ 判对 | ✅ 判对 |

正确判据是**归一化后与章节名精确相等**。

### 5.3 两层检查

`is_valid_title` 做两层检查：

1. **与已知章节名词典精确相等** → 无效（无论有没有传 `headings`）
2. **与文档中真实出现的章节名精确相等** → 无效（需要传 `headings`）

第 1 层必须在 `if headings:` 块**外面** —— 否则调用方一旦漏传 `headings`，章节名标题就会被判成有效。

### 5.4 其他无效条件

| 条件 | 理由 |
| --- | --- |
| 长度 < 10 | `"Short"` 显然不是论文标题 |
| 含 `\|` 且长度 < 20 | PDF 侧栏文字被拼进来了（如 `"Cell | Nature"`） |
| 无字母/汉字 | 纯数字或标点 |

---

## 6. 关键词抽取

### 6.1 三种形态

TEI 里的关键词有三种形态，优先级从高到低：

```xml
<!-- 形态 1：<term> 子元素（GROBID 标准输出） -->
<keywords>
  <term>guide RNA design</term>
  <term>deep learning</term>
</keywords>

<!-- 形态 2：逗号/分号分隔的纯文本 -->
<keywords>tumor necrosis factor, glycolysis pathway</keywords>

<!-- 形态 3：空格分隔（最后兜底） -->
<keywords>CRISPR Cas9</keywords>
```

### 6.2 老仓库的 bug

```python
keywords = keywords_text.split()   # ← 按空格切
```

`"tumor necrosis factor"` 被拆成 `["tumor", "necrosis", "factor"]` 三个关键词 —— 对下游的关键词共现分析、聚类完全是噪声。

**正确做法**：优先按 `<term>` 子元素取，其次按逗号/分号切分，最后才按空格切。

---

## 7. 引言首句提取

### 7.1 用途

引言的**每个段落的第一句**通常点明了这篇文章要解决什么问题。读 3 个段落的 3 句，比读完整引言快 5-10 倍 —— 这是批量筛读文献最有效的技巧。

### 7.2 缩写感知的断句

不能简单按 `"."` 切：

```python
text.split(".")[0]   # "Studies of e.g"  ← 半句话
```

必须排除的缩写：

| 缩写 | 例 |
| --- | --- |
| `e.g.` / `i.e.` | `"Studies of e.g. off-target effects..."` |
| `et al.` | `"See Smith et al. for..."` |
| `Fig.` / `Fig` / `Eq` / `No` / `ref` / `approx` | `"As shown in Fig. 3..."` |

实现思路：找到句末标点后，若该处匹配缩写模式则跳过，继续找下一个句末标点。

### 7.3 章节名匹配

用 `normalize_heading` 归一化后**精确相等**匹配，所以 `"Introduction"`、`"INTRODUCTION"`、`"INTRODUC TI ON"` 都能匹配上。

---

## 8. 语料索引与检索

### 8.1 索引结构

```json
{
  "files": {
    "normal.paper.pdf.tei.xml": {
      "title": "...",
      "abstract": "...",
      "keywords": ["..."],
      "sections": [{"title": "...", "level": 0, "paragraphs": ["..."]}]
    }
  }
}
```

### 8.2 三种检索范围

| scope | 检索对象 | 用途 |
| --- | --- | --- |
| `paragraph` | 正文段落 | 找具体论述（默认） |
| `section` | 章节标题 | 找有某章节的论文 |
| `title` | 论文标题 | 按主题找论文 |
| `all` | 三者 | 综合 |

### 8.3 单文件失败不中断整批

索引时若某个 TEI 文件 XML 损坏，**跳过并记录到 `errors`**，不中断整批。批量处理里这是必须的 —— 一篇坏文件不该让 200 篇白跑。

---

## 9. 修掉的 5 个老仓库缺陷

| # | 缺陷 | 现象 | 根因 | 修法 |
| --- | --- | --- | --- | --- |
| 1 | 关键词被空格切开 | `"tumor necrosis factor"` → 3 个词 | `keywords_text.split()` | 优先按 `<term>` 子元素取 |
| 2 | 标题误判 | `"Abstract Methods for Targeting X"` 被当成章节名 | 用"包含"匹配章节名 | 归一化后**精确相等** |
| 3 | 章节树重复计数 | 6 个章节变 11 个 | 兜底分支在每一层都触发 | 兜底只放顶层 `<body>` |
| 4 | 参考文献多算 1 条 | 实际 2 条统计成 3 条 | 全文遍历 `<biblStruct>`，把 `<sourceDesc>` 里本文自己的元数据也算进去 | 只取 `<back>` 下 |
| 5 | `xml:id` 取不到值 | 图表/文献 id 全为空 | `elem.get("xml:id")` 键名错 | 用 `{http://www.w3.org/XML/1998/namespace}id` |

另有两个老仓库的**硬编码个人路径**问题（`C:\Users\DELL\Desktop\trae\paperdata\...`）与 **`store_true` 配 `default=True` 关不掉选项**问题，在本技能中一并消除。

---

## 10. 输出文件名的推导

输出文件名必须从**文件名**推导，不能用完整路径。

`ParsedPaper.source_file` 存的是 TEI 文件的**绝对路径**。直接拿它拼：

```python
output_dir / "json" / f"{paper.source_file}.json"
```

在 Windows 上 `Path("out") / "json" / "C:\\Users\\...\\paper.json"` 会**锚定到盘符根目录** —— 文件写到输出目录**之外**去了。而中间目录 `out/json/` 仍会被创建出来，只断言"目录存在"根本发现不了这个问题。

正确做法：

```python
stem = Path(paper.source_file).name        # 只取文件名
stem = re.sub(r'[<>:"/\\|?*\x00-\x1f]', "_", stem)   # 净化非法字符
```

---

## 11. 参考文献抽取

### 11.1 范围

只取 `<back>` / `<listBibl>` 下的 `<biblStruct>`。

`<sourceDesc>` 里也有 `<biblStruct>`，但那是**当前这篇文献自己的元数据**（期刊名、DOI、卷期），不是参考文献。

### 11.2 字段

| 字段 | 来源 |
| --- | --- |
| `id` | `<biblStruct>` 的 `xml:id` |
| `raw` | 整个 `<biblStruct>` 的纯文本 |
| `title` | `<analytic><title>` |
| `doi` | `<idno type="DOI">` |

### 11.3 与引文网络分析对接

`doi` 字段是引文网络分析的接入点：把每篇论文的参考文献 DOI 与语料里的 DOI 对照，就能画出子图谱。只用 `raw` 纯文本则需要正则提取 DOI，容易出错。
