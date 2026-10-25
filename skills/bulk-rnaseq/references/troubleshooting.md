# bulk RNA-seq 排错手册

按"症状 → 原因 → 解决"的顺序组织。先看第 1 节定位问题类别。

---

## 0. 快速索引

| 症状 | 跳转 |
|---|---|
| `未找到外部工具 fastqc/fastp/salmon` | [1.1](#11-外部工具缺失) |
| `Salmon 索引不存在，跳过定量` | [1.2](#12-salmon-索引问题) |
| `tx2gene 文件缺少必需列` | [1.3](#13-tx2gene-映射不匹配) |
| `counts 与元数据的样本名完全无法匹配` | [1.4](#14-样本名对不上) |
| `分组样本量不足` | [1.5](#15-分组标签问题) |
| `R 后端 DESeq2 运行失败` | [1.6](#16-r-后端问题) |
| 差异基因数为 0 | [2.1](#21-没有显著差异基因) |
| 几乎所有基因都显著 | [2.2](#22-几乎所有基因都显著) |
| log2FC 全是 0 或缺失 | [2.3](#23-log2fc-异常) |
| 火山图只有一条水平带 | [2.4](#24-所有-padj-相同) |
| 富集分析没有结果 | [2.5](#25-富集分析没有结果) |
| 富集结果全部显著 | [2.6](#26-富集结果全是显著的) |
| size factor 极端 | [3.1](#31-size-factor-异常) |
| counts 全是 0 | [3.2](#32-counts-矩阵全零) |
| 内存溢出 | [3.3](#33-内存溢出) |

---

## 1. 报错类问题

### 1.1 外部工具缺失

**日志**：`未找到外部工具 fastqc，跳过依赖它的步骤`

**原因**：`fastqc` / `fastp` / `salmon` 未安装或不在 `PATH`。

**解决**：

```bash
# conda
conda install -c bioconda fastqc fastp salmon

# 检查
which fastqc fastp salmon
```

**重要**：这三个工具缺失**不会**中断流程，纯计算步骤（counts 矩阵、DESeq2、富集、可视化）依然可跑。若已有上游产物，用 `--skip-qc --skip-trimming --skip-quant` 直接进入统计阶段。

---

### 1.2 Salmon 索引问题

**日志**：`Salmon 索引不存在，跳过定量` 或 `未找到可用的转录组 FASTA`

**原因**：
- 未提供 `--transcriptome`
- FASTA 路径不存在或格式不对
- 索引目录存在但为空（上次构建中断）

**解决**：

```bash
# 手动验证 FASTA
head -2 /ref/transcriptome.fa
# 应该有 > 开头的标题行

# 手动构建索引
salmon index -t /ref/transcriptome.fa -i /tmp/salmon_index -k 31 -p 8
```

**FASTA 标题必须包含转录本 ID**，且与 `tx2gene.tsv` 的 `transcript_id` 完全一致：

```bash
# 对比两边 ID
grep '^>' /ref/transcriptome.fa | sed 's/>//' | cut -f1 -d' ' | sort > /tmp/tx.txt
cut -f1 /ref/tx2gene.tsv | tail -n +2 | sort > /tmp/tx2gene.txt
comm -23 /tmp/tx2gene.txt /tmp/tx.txt | head   # 有输出 = 有转录本对不上
```

---

### 1.3 tx2gene 映射不匹配

**日志**：`tx2gene 文件缺少必需列 transcript_id / gene_id`

**原因**：TSV 的列名不是这两个。常见的是 `tx` / `gene`、`transcript` / `gene_name`、`#transcript_id`（带 `#` 前缀）。

**解决**：检查并重命名列名：

```bash
head -1 /ref/tx2gene.tsv
# 期望: transcript_id<TAB>gene_id
```

另外注意：基因 ID 空间必须与注释文件一致。例如注释文件用 Ensembl ID（`Bca003g000001`）而 tx2gene 用 NCBI（`BCA_G000001`），会导致富集分析完全匹配不上。

---

### 1.4 样本名对不上

**日志**：`counts 与元数据的样本名完全无法匹配`

**原因**：counts 的列名与 metadata 的 `sample` 列值不同。最常见的是：

- counts 列名带后缀（`sample_1` vs `sample1`）
- metadata 用 Excel 打开过，样本名尾部多了空格
- metadata 的样本顺序/名称与实际不一致

**解决**：

```bash
python - <<'EOF'
import pandas as pd
c = pd.read_csv("04_counts/counts.csv", index_col=0)
m = pd.read_csv("metadata.csv")
print("counts cols:", list(c.columns))
print("meta sample:", list(m["sample"]))
print("only in counts:", set(c.columns) - set(m["sample"]))
print("only in meta:", set(m["sample"]) - set(c.columns))
EOF
```

修好 metadata 后重跑。脚本会自动取交集并给出警告，只要不是完全为空就能继续。

---

### 1.5 分组标签问题

**日志**：`分组样本量不足: treated=0, control=3`

**原因**：`--group-a` / `--group-b` 的值与 metadata 里 `condition` 列的实际取值不匹配（大小写、下划线、空格）。

**解决**：

```bash
python -c "import pandas as pd; print(pd.read_csv('metadata.csv')['condition'].unique())"
```

用输出中的**精确字符串**作为 `--group-a` / `--group-b`。

**另一类症状**：分组标签写反（treated 和 control 互换）。表现为所有 log2FC 符号颠倒、上调基因变成下调。用火山图的对称性检查，或对比一个已知上调基因的符号来确认。

---

### 1.6 R 后端问题

**日志**：`未找到 Rscript，无法使用 r 后端`

**解决**：

```bash
conda install -c conda-forge r-base
Rscript -e 'if (!requireNamespace("DESeq2", quietly=TRUE)) BiocManager::install("DESeq2")'
```

**日志**：`R DESeq2 执行失败`

常见原因与排查：

| 报错片段 | 原因 | 解决 |
|---|---|---|
| `there is no package called 'DESeq2'` | Bioconductor 包没装 | `BiocManager::install("DESeq2")` |
| `DESeqDataSetFromMatrix: some variables in design formula are characters` | design 公式里的因子是 character 而非 factor | 脚本已自动转 factor，若仍报错检查 `--design` 是否引用了 metadata 中不存在的列 |
| `contrasts can be applied only to factors with 2 or more levels` | condition 只有一个取值 | 检查分组标签 |
| `libcurl error` / 无法访问 Bioconductor | 网络问题 | 用 `--backend python` |

**验证 R 环境**：

```bash
Rscript -e 'suppressPackageStartupMessages(library(DESeq2)); cat("DESeq2 OK\n")'
```

---

## 2. 结果异常类问题

### 2.1 没有显著差异基因

**症状**：`差异表达完成: N 个基因, 0 个显著`

**按可能性排序**：

1. **分组标签写反或写错**（最常见）。检查 `size_factors.csv` 与 `deseq2_results_full.csv` 中 log2FC 的符号是否符合预期。
2. **`--lfc-threshold` 太大**。默认 1.0 表示要求 2 倍以上变化。若预期是 1.5 倍，用 `--lfc-threshold 0.585`（= log2(1.5)）。
3. **`--alpha` 太严**。用 `--alpha 0.1` 看趋势。
4. **生物学上确实无差异**。查看 `deseq2_results_full.csv` 中 pvalue 最小的那批基因，如果 pvalue 都在 0.05 附近，说明功效不足而非无差异。
5. **样本量不足**。每组 2 个重复时，功效非常有限。

**诊断脚本**：

```bash
python - <<'EOF'
import pandas as pd
r = pd.read_csv("05_deseq2/deseq2_results_full.csv", index_col=0)
print("min pvalue:", r["pvalue"].min())
print("n padj<0.1:", (r["padj"] < 0.1).sum())
print("n padj<0.05:", (r["padj"] < 0.05).sum())
print(r.nsmallest(10, "pvalue")[["baseMean","log2FoldChange","pvalue","padj"]])
EOF
```

---

### 2.2 几乎所有基因都显著

**症状**：100 个基因里 95 个 `padj < 0.05`

**原因**（按可能性排序）：

1. **批次与条件完全重合**。若所有 treated 样本同时来自一个测序批次，批次效应被误判为生物学效应。
2. **对照不是真正对照**。例如对照来自细胞系而实验组来自组织。
3. **元数据条件列填错**，导致所有样本被归为两组但组内实际差异巨大。

**解决**：

- 检查 `07_figures/Figure4_PCA.png`：若两组在 PC1 上完美分离但组内也沿 PC1 展开，说明是连续变化而非二分类。
- 用 PCA 前几个 PC 做相关性检验，评估批次效应强度。
- 排除批次：R 后端可用 `--design "~batch + condition"`。

---

### 2.3 log2FC 异常

**症状 A**：全是 0 或 NaN

**原因**：设计矩阵奇异（所有样本在同一组，或组别与某个协变量完全共线）。

**解决**：确认 `--group-a` 与 `--group-b` 各有样本。

**症状 B**：log2FC 数值大得不合理（比如 15）

**原因**：某样本的 counts 极低（近似 0），导致比值爆炸。

**解决**：

```bash
python -c "import pandas as pd; c=pd.read_csv('04_counts/counts.csv',index_col=0); print(c.sum())"
```

若某样本总量远低于其他样本，应剔除该样本或排查建库失败。

---

### 2.4 所有 padj 相同

**症状**：火山图只有一条水平带，所有点 `-log10(padj)` 一样。

**原因**：BH 校正在极端情况下会给出相同的 q 值（例如所有 p 值都远小于最小 p 值时）。这本身不一定是错误，但如果 padj 完全相同且 padj > 0.05，说明所有 p 值都很接近。

**排查**：看 `pvalue` 列是否有区分度。若 pvalue 也有大量相同值，说明数值精度问题。

---

### 2.5 富集分析没有结果

**日志**：`没有差异基因，跳过富集分析` 或 `未产生富集结果`

**按顺序排查**：

1. **没有显著 DEG** → 回到 [2.1](#21-没有显著差异基因)。
2. **注释文件没匹配上基因** → 检查日志中"注释基因 N 个"，若远小于你的基因总数，说明 ID 空间不匹配。
3. **基因集大小过滤太严** → `--min-gene-set 3 --max-gene-set 1000` 再试。
4. **DEG 数量太少** → 少于 `min_gene_set` 时超几何检验没有功效。
5. **注释文件列名不对** → 脚本自动识别 `EnsemblPlants`/`gene_id`/`GeneID`/`gene`/`Gene` 作为基因列；GO 列必须是 `Gene Ontology (biological process)` 等完整名称。列名不符时该本体会被整体跳过。

**检查注释匹配率**：

```bash
python - <<'EOF'
import pandas as pd
c = pd.read_csv("04_counts/counts.csv", index_col=0)
go = pd.read_csv("go_annotation.csv")
col = "EnsemblPlants" if "EnsemblPlants" in go.columns else go.columns[0]
gids = set(go[col].astype(str).str.rsplit(".", n=1).str[0])
print("注释基因数:", len(gids))
print("与 counts 交集:", len(gids & set(c.index)))
EOF
```

---

### 2.6 富集结果全是显著的

**症状**：几乎每个 GO term 的 `padj` 都 < 0.05

**原因**：**背景基因集选错**。本实现默认用注释文件中的全部基因作背景（N 很大），这会高估富集显著性。

**正确做法**：用"实验中实际检测到（且通过过滤）的基因"作为背景，即 `counts.csv` 的行集合。

**解决**：改造 enrichment 步骤，把背景限定为 counts 中出现的基因：

```bash
python - <<'EOF'
import pandas as pd
go = pd.read_csv("go_annotation.csv")
col = "EnsemblPlants"
c = pd.read_csv("04_counts/counts.csv", index_col=0)
measured = set(c.index)
go = go[go[col].isin(measured)]
go.to_csv("go_annotation_filtered.csv", index=False)
print("过滤后注释基因数:", len(go))
EOF
python skills/bulk-rnaseq/scripts/run_pipeline.py ... --go-annotation go_annotation_filtered.csv
```

---

## 3. 数据质量问题

### 3.1 size factor 异常

**症状**：`size_factors.csv` 中某个值远大于 1 或远小于 1

**判断标准**：正常范围约 0.7–1.4。超出说明：

| 现象 | 原因 |
|---|---|
| 某样本 sf 明显偏大 | 该样本测序深度低（文库上样量不足） |
| 某样本 sf 明显偏小 | 该样本测序深度高，或存在大量重复序列/污染 |
| 组内 sf 差异大组间一致 | 正常深度差异 |
| 组间 sf 系统性差异 | 可能提示两组建库方式不同，需警惕混杂 |

**处置**：
- 深度差异大：可在分析时保留（DESeq2 的 size factor 就是为此设计的）
- 极端异常样本：剔除后重跑，并在文章中说明剔除标准

---

### 3.2 counts 矩阵全零

**症状**：`counts.csv` 全为 0，或基因数极少

**原因**：
1. Salmon 定量失败（`quant.sf` 里 NumReads 全 0）→ 检查 Salmon 的日志
2. tx2gene 的基因 ID 全对不上 → 汇总后所有基因名变成原始转录本名，与注释不匹配
3. `--min-counts` 设得太大

**排查**：

```bash
# 直接看 Salmon 原始输出
head -5 03_salmon/sample1/quant.sf
python -c "import pandas as pd; c=pd.read_csv('04_counts/counts.csv',index_col=0); print(c.sum(axis=1).describe())"
```

---

### 3.3 内存溢出

**症状**：`MemoryError` 或进程被系统杀死

**原因**：差异表达时对每个基因做网格搜索 + IRLS，基因数 × 网格点数会占用较多内存；富集分析会构建基因集字典。

**解决**：

1. **降低 `--min-counts`**，去掉大量低表达基因（同时加快速度）
2. **分块处理**：本实现已按基因逐个处理，内存主要来自 counts 矩阵本身；基因数超过 5 万时建议先做一次方差过滤
3. **减少重复基因**：确认 tx2gene 正确，避免把 3 万个转录本当成 3 万个基因
4. **分块运行富集**：把显著 DEG 分成几批分别跑

---

## 4. 环境与依赖问题

### 4.1 `ModuleNotFoundError: No module named 'numpy'`

```bash
pip install numpy scipy pandas matplotlib
```

### 4.2 `ModuleNotFoundError: No module named 'yaml'`（仅自检脚本需要）

```bash
pip install pyyaml
```

### 4.3 `No module named 'biopipe_core'`

脚本会把仓库的 `packages/` 目录加入 `sys.path`，但需要正确的目录层级。若是脱离仓库单独拷贝脚本使用，需要手工设置：

```bash
export PYTHONPATH=/path/to/repo/packages
```

### 4.4 matplotlib 中文/字体缺失

图例使用英文标签（避免字体问题）。若要改中文，需安装中文字体并设置 `font.sans-serif`。

---

## 5. 结果可重复性

- 每次运行会在输出根目录生成 `run_manifest.json`，记录输入路径、后端、对比方向与时间戳。
- DESeq2 的 python 后端是确定性算法（无随机采样），同样的输入必然得到同样的输出。
- 可视化中的 PCA 使用 SVD 分解（确定性），不含随机种子依赖。

---

## 6. 何时该怀疑是"流程 bug"而非"数据问题"

用下面的对照实验快速区分：

1. **造一份合成数据跑通全流程**（counts 矩阵 + 元数据），确认能稳定产出显著结果 —— 排除代码路径问题。
2. **同一份数据用 `--backend python` 和 `--backend r` 各跑一次**：
   - 两者结论一致 → 流程可信，是数据/设计问题
   - 结论矛盾 → 优先怀疑 python 后端的实现细节
3. **检查 `size_factors.csv`**：若数值明显异常，问题在 counts 矩阵构建阶段而非统计阶段。
