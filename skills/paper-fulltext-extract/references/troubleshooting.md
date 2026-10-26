# 排错手册

按"症状 → 原因 → 修法"组织。每条都对应 `paper-fulltext-extract` 技能里的一个真实坑。

---

## 1. GROBID 相关

### 1.1 `无法连接 GROBID 服务 http://127.0.0.1:8070`

**最常见的原因**：容器还在下模型，**没有启动完成**。

GROBID 首次启动要下载模型，日志里要等这一行出现：

```
Server started
```

在此之前任何请求都会连接失败。**这不是配置问题，是还没就绪。**

排查顺序：

```bash
# 1) 看容器日志，确认是否出现 "Server started"
docker logs -f <container>

# 2) 单独探测健康检查端点
curl http://127.0.0.1:8070/api/isalive     # 应输出 true
```

其他可能：

| 现象 | 原因 | 修法 |
| --- | --- | --- |
| `Connection refused` | 服务未启动或端口映射错 | 确认 `-p 8070:8070` |
| `Connection refused` 但日志正常 | 服务绑了别的网卡 | 改用 `docker run -p 8070:8070` 而非 `-p 127.0.0.1:8070:8070` |
| 超时 | 首次下载模型 | 等模型下完 |
| 404 | URL 拼接错 | 见 1.2 |

### 1.2 请求返回 404，路径看着却又像对的

**根因**：URL 拼接时少了一个斜杠。

```python
# 错
url = f"{server}api/processFulltextDocument"
# server = "http://127.0.0.1:8070"  →  "http://127.0.0.1:8070api/processFulltextDocument"
```

老仓库的默认值带尾斜杠（`http://127.0.0.1:8070/`）所以"能用"，但用户一旦传不带斜杠的地址就 404，而且报错完全指不到真正原因。

**修法**：用 `join_url` 统一处理，本技能已内置：

```python
join_url("http://127.0.0.1:8070", "api/x")   # → "http://127.0.0.1:8070/api/x"
join_url("http://127.0.0.1:8070/", "api/x")  # → "http://127.0.0.1:8070/api/x"
join_url("http://h:8070///", "//api/x")      # → "http://h:8070/api/x"
```

### 1.3 `consolidateHeader` 关不掉

**根因**：布尔参数必须传**小写字符串**。

传 Python `bool` 会被编码成 `"True"` / `"False"`，GROBID 按**非空即真**处理 —— `"False"` 也是非空字符串，于是当成 `true`。**不报任何错。**

```python
params = {"consolidateHeader": str(bool(v)).lower()}   # "true" / "false"
```

本技能的实现已统一处理。

### 1.4 部分论文处理失败，HTTP 500/503

**根因**：GROBID 处理特定版式的 PDF 时会崩或过载。

| 状态码 | 含义 | 应对 |
| --- | --- | --- |
| 500 | 该 PDF 让 GROBID 处理崩溃 | 换 `--flavor speed` 重试；仍失败则跳过 |
| 503 | 服务过载 | 加大 `--pause`（如 `3`），或降低并发 |
| 429 | 限速 | 同上 |

本技能对 429/502/503/504 会自动重试（指数退避），对 4xx 不重试。

**批量跑的建议参数**：

```bash
python skills/paper-fulltext-extract/scripts/grobid_client.py \
    --input ./pdfs --output ./tei --pause 2 --max-retries 5
```

### 1.5 处理的 PDF 一半是空的

**根因**：`grobid_client.py` 的 `--input` 只扫描**当前目录**（非递归）。

子目录里的 PDF 不会被处理。要处理子目录需先展平，或分别指定。

---

## 2. TEI 解析相关

### 2.1 论文标题显示成 `INTRODUCTION`

**根因**：GROBID 误把章节名填进了 `<title type="main">`。常见于页眉有大字号内容的排版。

**这是 GROBID 的已知行为，不是 bug。** 本技能会自动检测并回退到文件名：

```xml
<title type="main">INTRODUCTION</title>   →  标题变成 "paper_name"
```

**回退到文件名是正确的**。如果你的 TEI 里有其他 `<titleStmt>` 下的候选标题，脚本会优先用那个。

**注意**：回退**不会**取到期刊名。`<monogr><title>` 是期刊名，误用它比用文件名更糟 —— 期刊名在检索时毫无区分度。

### 2.2 正常标题被误杀（标题变成了文件名）

**症状**：一篇标题明确正常的论文，输出里标题却是文件名。

**根因**：判断逻辑用了"包含"匹配 —— 标题里含有 `introduction` / `methods` / `abstract` 之一就被当成章节名。

| 标题 | 是否含章节名 | 是否应保留 |
| --- | --- | --- |
| `"Abstract Methods for Targeting X"` | 含 abstract | ✅ 应保留 |
| `"Introduction to Materials Chemistry"` | 含 introduction | ✅ 应保留 |
| `"INTRODUCTION"` | 含 introduction | ❌ 应判无效 |

**正确判据是归一化后精确相等**。本技能已修正。若你看到误杀，检查是否有人把 `is_valid_title` 改回了"包含"匹配。

### 2.3 关键词被拆成一个个单词

**症状**：`"tumor necrosis factor"` 变成了 `tumor` / `necrosis` / `factor` 三条。

**根因**：走了**空格分隔**的兜底分支 —— 说明 TEI 的 `<keywords>` 下**没有** `<term>` 子元素，且内容里**没有逗号/分号**。

**排查**：打开对应的 TEI 看：

```xml
<keywords>tumor necrosis factor, glycolysis pathway</keywords>
```

如果长这样，说明是**纯文本形态**，本技能会按逗号切（正确）。如果确实被拆开，看内容是不是：

```xml
<keywords>tumor necrosis factor glycolysis pathway</keywords>
```

这种**无分隔符的多词关键词在信息论上就是歧义的** —— 无法可靠还原边界。此时应该修数据源，或接受拆分。

### 2.4 参考文献条数比实际多 1

**症状**：论文末尾明明 2 条参考文献，统计出 3 条。

**根因**：把 `<sourceDesc>` 里**本文自己的元数据**（`<biblStruct>`）也当成了参考文献条目 —— 它排在最前面。

**修法**：只从 `<back>` / `<listBibl>` 下取条目。本技能已修正。

### 2.5 章节标题重复出现

**症状**：章节树里每个标题出现两次，6 个章节变 11 个。

**根因**：处理"无 `<div>` 变体"的兜底分支写在了递归函数里，导致**每一层**都触发一次，把每个 `<div>` 自己的 `<head>` 又解析了一遍。

**修法**：兜底只在顶层 `<body>` 做。本技能已修正。

**自查**：

```bash
python skills/paper-fulltext-extract/scripts/search_corpus.py \
    --xml-dir ./tei --show-structure normal.paper.pdf.tei.xml
```

标题不该有重复项。

### 2.6 章节名匹配不上（`INTRODUC TI ON`）

**症状**：`--intro-section introduction` 提不出首句，但 TEI 里明明有引言。

**根因**：PDF 抽取把章节名的字间距当成了空格。

**修法**：归一化时去空格后与已知章节名词典比较。本技能已修正，`INTRODUC TI ON` / `DISCUSS I ON` 都能匹配。

如果你有非常规的章节名（如 `Materials & Experimental Protocols`），用 `--intro-section` 传准确的名字。

### 2.7 图表/文献的 id 全是空

**症状**：JSON 里 `figures[].id` 与 `references[].id` 都是 `""`。

**根因**：`elem.get("xml:id")` —— 键名错。ElementTree 里 `xml:id` 的真实键是 `{http://www.w3.org/XML/1998/namespace}id`。

**不报错，只是永远取不到值**，非常隐蔽。本技能已用 `XML_ID_ATTR` 常量修正。

### 2.8 XML 解析报 `mismatched tag`

**症状**：某个 TEI 文件报语法错误。

**可能原因**：

1. GROBID 输出被截断（处理中断 / 磁盘满）—— 重新用 GROBID 处理该 PDF
2. 文件传输过程中损坏
3. 文件本身不是 TEI（扩展名是 `.xml` 但内容是别的）

**影响范围可控**：单篇失败不会中断整批，损坏文件记入警告并跳过，索引里也会记录在 `errors`。

---

## 3. 输出相关

### 3.1 输出目录里出现奇怪的嵌套层级

**症状**：`out/json/` 下面出现了 `tests/fixtures_tei/normal.paper.json` 这样的嵌套目录。

**根因**：输出文件名用了 TEI 的**完整路径**而不是文件名。

在 Windows 上更糟 —— `Path("out") / "json" / "C:\\Users\\...\\x.json"` 会**锚定到盘符根目录**，文件直接写到输出目录**之外**去，而中间的 `out/json/` 目录仍然会被创建出来。

**修法**：

```python
stem = Path(paper.source_file).name    # 只取文件名
stem = re.sub(r'[<>:"/\\|?*\x00-\x1f]', "_", stem)
```

本技能已修正。**自查**：

```bash
ls out/json/            # 应该是平铺的 .json 文件，没有子目录
```

### 3.2 HTML 里的图片不显示

**症状**：HTML 生成了，但 `<img>` 是裂图。

**根因**：TEI 里的图形路径是相对路径（如 `figures/fig_1.png`），这些文件在 **GROBID 的输出目录**里，不在 HTML 旁边。

**修法**：把 GROBID 输出的图片目录拷到 HTML 旁边，然后用 `--image-base` 指过去：

```bash
cp -r tei/figures out/html/figures
python ... --image-base figures
```

### 3.3 `--formats` 里写了不认识的格式

**症状**：报错 `未知的输出格式: ['pdf']`。

**可用格式**只有 `json` / `markdown` / `html`。要 PDF 请先用 Markdown 再转。

### 3.4 CSV 用 Excel 打开是乱码

**本技能导出的 CSV 已带 UTF-8 BOM**（`utf-8-sig`），Excel 直接双击打开就是正常的中文。若仍乱码，说明是别处导出的文件。

---

## 4. 检索相关

### 4.1 检索命中数为 0

按顺序排查：

1. **范围选错** —— 默认 `--scope paragraph` 只搜正文。搜章节名要用 `--scope section`，搜论文标题用 `--scope title`，全都要用 `--scope all`
2. **大小写** —— 默认不区分大小写，`CRISPR` 能搜到 `crispr`。加 `--case-sensitive` 后区分
3. **词形不同** —— 搜 `"gene"` 搜不到 `"genes"`、`"genetic"`。本技能是**子串匹配**，不做词形还原
4. **文件没进索引** —— 看 `corpus_index.json` 里的 `files` 键，或跑 `--list-files`

### 4.2 章节检索匹配到断字标题

**症状**：`--scope section --keyword introduc` 命中了标题为 `INTRODUC TI ON` 的章节。

**这是正确的**。子串匹配天然能命中断字形式。若要做精确的章节名查找，看 `--show-structure` 的输出。

### 4.3 索引建立很慢

索引要解析每个 TEI 的全部正文。几百篇文献时正常。

优化：用 `--flavor speed` 处理 PDF（TEI 本身更小）、或分批处理。索引结果可以 `--index-file` 保存，之后用 `--load-index` 直接复用，不必重建。

---

## 5. 自检清单

改动本技能代码后，按这个顺序自检：

| # | 检查项 | 命令 / 方法 |
| --- | --- | --- |
| 1 | 单元测试全绿 | `python -m pytest tests/test_paper_fulltext_extract.py -q` |
| 2 | CLI 能打出用法 | `python skills/paper-fulltext-extract/scripts/extract_fulltext.py --help` |
| 3 | 无依赖环境能启动 | `python scripts/ci_smoke_no_deps.py` |
| 4 | 技能规范合规 | `python scripts/validate_skills.py` |
| 5 | 流水线端到端跑通 | `python .../extract_fulltext.py --tei-dir tests/fixtures_tei --output /tmp/o --formats json,markdown,html --build-index` |
| 6 | 产物平铺在格式目录下 | `ls /tmp/o/json/`（不应有子目录） |
| 7 | 关键词未被拆开 | 看 `json/*.json` 里的 `keywords` 数组 |
| 8 | 参考文献条数正确 | 看 `json/*.json` 里的 `references` 长度 vs 论文实际条目数 |
| 9 | 章节无重复 | 看 `corpus_index.json` 里的 `sections[].title` |
| 10 | 引言首句完整 | 看 `intro_sentences.csv`，每句应以句号结尾 |
