# 文献挖掘常见问题与排错

按"症状 → 原因 → 解决"组织。多数问题可先用 `-v`（DEBUG 日志）定位。

---

## 1. 完全检索不到结果（0 条）

### 1.1 检索式太窄

**症状**：三个源都返回 0 条。

**原因**：`--topic` 与 `--keywords` 是 AND 关系，两者同时要求命中；再叠加 `--date-range` 就更容易归零。

**解决**：
```bash
# 先单独试 topic，不带 keywords 和日期
python skills/literature-search/scripts/search_literature.py \
    --topic "antimicrobial peptide" --source pubmed --max-results 10 --output /tmp/t.csv
```
若这一步有结果，逐步加回 `--keywords`、再逐步加回 `--date-range`，定位是哪个条件卡死的。

### 1.2 年份过滤区间反了

**症状**：报 `起始年份 X 大于结束年份 Y`。

**原因**：`--date-range` 格式写成了 `2024:2018`。

**解决**：改成 `2018:2024`。注意 `:` 左边是起始、右边是结束。

### 1.3 PubMed 的 `[PDAT]` 过滤过严

PubMed 的 `[PDAT]` 字段只覆盖正式出版日期，在线先发表但尚未编入卷期的文献会被漏掉。

**解决**：放宽起始年份，或先不设日期范围再在本地用 `--min-relevance` 筛。

---

## 2. HTTP 429 / 请求失败

### 2.1 Semantic Scholar 限速最严

**症状**：日志出现 `请求失败(第 N 次)，... 后重试`，最终 `HTTP 429`。

**原因**：无 API Key 时 S2 严格限制约 1 req/s。本技能已自动把 S2 的间隔提高到 1.1 s，但 `--max-results` 大（>100）时总量仍然很大。

**解决**：
```bash
export S2_API_KEY=xxxx          # 最有效
# 或减少查询量
--source pubmed openalex        # 先不查 S2
--max-results 30
```

### 2.2 NCBI API Key

无 Key 时 NCBI 限 3 req/s。技能默认间隔 0.34 s 恰好卡在上限，批量 `esummary` 分片（200 个 PMID/次）时容易触发 429。

**解决**：`export NCBI_API_KEY=xxxx`（免费申请，10 req/s）。

### 2.3 重试耗尽

共享层默认 `retries=4`、`backoff=2.0`，最坏情况会等待 `2^0+2^1+2^2+2^3 = 15` 秒再失败。网络被墙或代理异常时所有源都会失败。

**诊断**：
```bash
python skills/literature-search/scripts/search_literature.py -v ... 2>&1 | head -50
```
若所有源同时报连接错误而非 429，那是网络层问题，不是配额问题。

---

## 3. OpenAlex 请求返回 403

**原因**：OpenAlex 的"礼貌池"要求提供 `mailto` 参数，缺少时可能对高频请求降级处理。

**解决**：`--email your.address@example.com`

同时 OpenAlex 端点已不再支持无 `select` 的全量字段返回；本脚本已固定 `select` 字段列表，避免 400。

---

## 4. 摘要列大量为空

### 4.1 忘了开 `--enrich-abstracts`

PubMed 的 `esummary` **不返回摘要**（只有题录），OpenAlex 的 `works` 列表端点在未指定时也不返回摘要。本脚本已对 OpenAlex 请求了 `abstract_inverted_index`，但很多文献在源库本身就无摘要。

**解决**：加 `--enrich-abstracts`：
```bash
--enrich-abstracts --min-abstract-for-enrich 60
```

### 4.2 OpenAlex 摘要重建后为空

**症状**：`abstract` 列出现空白，但 `abstract_inverted_index` 在原始 JSON 中确实存在。

**原因**：倒排索引重建时按位置拼接，若源数据的位置数组有空洞，理论上会丢词（但不致全空）。

**更可能的原因**：这些文献在 OpenAlex 中根本没有摘要字段（`abstract_inverted_index` 缺失或为 null）。

**解决**：接受现实，或改用 `--enrich-abstracts` 让 PubMed/S2 兜底。

### 4.3 PubMed efetch 返回的是 XML 不是 JSON

本脚本的 `_parse_pubmed_abstracts` 用 `xml.etree.ElementTree` 解析 `retmode=xml` 的响应，抽取 `.//Abstract/AbstractText` 节点并保留 `Label`（如 `BACKGROUND:`、`METHODS:`）。

若你手工调用 efetch，注意**不能用 `retmode=json`**——efetch 没有 JSON 模式。

---

## 5. 去重不干净 / 过度去重

### 5.1 同一文献仍出现多条

**原因**：跨源重复的文献没有共同 DOI（部分早期文献只有 PMID 没有 DOI），此时只能靠标题模糊匹配。

**解决**：
```bash
--title-similarity 0.88     # 放宽阈值
```
注意阈值越低越容易把"同一课题组发的系列工作"误判为重复。

### 5.2 不同文献被合并成一条

**症状**：`source` 列出现 `PubMed+OpenAlex`，但标题其实不同。

**原因**：标题过短（如 "Analysis of p53"）或副标题被不同源截断。

**解决**：
```bash
--title-similarity 0.96
```
并检查是否有 DOI 字段为空导致走了模糊匹配路径。

**注意**：本技能已内置 DOI 冲突保护 —— 两条记录都带 DOI 且不同时，不会被标题相似度合并。
实测检索 "antimicrobial peptide" 时 OpenAlex 前 5 条里有 3 条标题同为 "Antimicrobial Peptides"
但 DOI 不同，属三篇不同论文，已全部保留。
若看到标题相同但 DOI 不同的记录被合并，说明 DOI 字段在导入时丢失了（需检查 CSV 的 `doi` 列名与取值）。

### 5.3 `source` 列出现 `A+B` 是正常现象

这表示该文献在两个源都被命中，信息合并了，不是 bug。

---

## 6. 相关性得分不合理

### 6.1 全是 0 分

**原因**：没有传 `--keywords`，而相关性打分完全依赖关键词命中 + 引用数 + 年份。若既无关键词又都是零引用零近期文献，得分必然为 0。

**解决**：传 `--keywords`，或改用 `--min-citations` 做筛选。

### 6.2 得分普遍过高

**原因**：`--keywords` 里的词太宽泛（如 "protein" 会命中几乎所有文献）。

**解决**：换成更有区分度的术语（"molecular docking" 而非 "docking"），或收紧 `--min-relevance`。

---

## 7. 中文乱码

### 7.1 Excel 打开 CSV 乱码

本脚本用共享层 `write_csv` 写出，编码是 **UTF-8 with BOM**（`utf-8-sig`），Excel 双击打开不会乱码。

若你用其他工具重新保存成了无 BOM 的 UTF-8，Excel 就会误判编码。可用共享层 `ensure_unicode_csv()` 修复。

### 7.2 Windows 控制台输出乱码

日志走 stderr，若终端代码页不是 UTF-8，Windows 上可能出现乱码字符。日志内容均为 ASCII 关键字（中文日志在 Windows 终端下同样可能出现编码问题）。不影响文件写出。

### 7.3 引用格式文件里的中文

`export_citations.py` 默认 `--encoding utf-8-sig`（带 BOM），若目标工具（LaTeX/BibTeX）不认 BOM，请显式传 `--encoding utf-8`。

---

## 8. 文件与路径问题

### 8.1 输出目录不存在

脚本会自动创建 `--output` 的父目录。若未指定 `--output`，默认写到 `<root>/results/literature_results.csv`，`<root>` 取自环境变量 `BIOPIPE_PROJECT_ROOT` 或当前工作目录。

### 8.2 "不要在脚本里执行网络请求的默认行为"

本技能遵守规范：**只有不传 `--input` 时才会联网**。传 `--input` 走纯离线路径（去重、打分、筛选、导出），不会发出任何请求。

---

## 9. 自检命令

```bash
cd <repo>
export PYTHONPATH="$PWD/packages"

# frontmatter 校验
python -c "import yaml; d=yaml.safe_load(open('skills/literature-search/SKILL.md').read().split('---')[1]); print(d['name'], '|', d['description'][:80])"

# CLI 自检
python skills/literature-search/scripts/search_literature.py --help
python skills/literature-search/scripts/export_citations.py --help

# 离线端到端（不联网）
python skills/literature-search/scripts/search_literature.py \
    --input <你的csv> --keywords "demo" --output /tmp/out.csv
python skills/literature-search/scripts/export_citations.py \
    --input /tmp/out.csv --format bibtex --output /tmp/out.bib
```

---

## 10. 结果可信度自查清单

在把结果写进论文或综述前，逐项确认：

- [ ] 三个源的返回条数都 > 0（某个源为 0 说明检索式对该源不适用）
- [ ] `source` 列的分布合理（不是单一源独占）
- [ ] 摘要完整度 > 80%（低于此值说明需要 `--enrich-abstracts`）
- [ ] 抽查 `relevance` 前 20 与后 20 条，确认排序有效
- [ ] 引用格式导出后导入 Zotero 验证作者与卷期解析
- [ ] 关键结论类文献**手工阅读全文**，不要仅依赖自动化筛选
