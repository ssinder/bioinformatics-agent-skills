# bulk RNA-seq 方法学原理与参数选择

本文解释 `bulk-rnaseq` 技能中每一步的统计原理、公式，以及参数应该如何选。

---

## 1. 质控与过滤

### FastQC 关注点

| 指标 | 期望 | 异常含义 |
|---|---|---|
| Per base sequence quality | 5' 端略低，R2 端最低（Illumina 双端测序特性），其余 Q30+ | 全序列低质量 → 文库降解或仪器问题 |
| Per tile sequence quality | 单一 tile 偏离 → flow cell 局部问题 | 大量 tile 异常 → flow cell 缺陷 |
| Adapter content | 无接头或含量极低 | R1/R2 接头残留 → 需要 trim |
| GC content | 单一峰，与物种基因组一致 | 双峰 → rRNA 污染或污染 DNA |
| Sequence length | 分布集中 | 长度分布过宽 → 片段化严重 |

### fastp 参数

- `--qualified_quality_phred 20`：一端质量低于 Q20 的比例超过 40% 则丢弃该端。RNA-seq 常用 Q15–Q20。数据质量差时降到 Q15，质量好可提到 Q28。
- `--length_required 36`：短于 36 bp 的 read 难以唯一比对到转录本，直接丢弃。Illumina 常见 R1 50 bp / R2 75 bp 时，36 是安全下限。
- `--detect_adapter_for_pe`：双端自动检测接头，比手工指定更稳健。

**权衡**：过滤越严保留的有效片段越少，测序深度越浅。过滤后保留比例低于 70% 时应回退阈值重新评估。

---

## 2. 定量：为什么用 Salmon 而不是 STAR

Salmon 采用基于 k-mer 索引的 EM 算法做转录本级比对与计数，得到：

- `NumReads`：估计的映射到该转录本的 read 数（已做长度/GC 偏差校正）
- `TPM`：每百万转录本 reads 的转录本表达量（样本内可加和）
- `EffectiveLength`：校正 read 长度后该转录本的有效长度
- `Length`：转录本原始长度

Salmon 比 STAR 快 10–50 倍，且不需要比对到基因组，是差异分析的首选方案。

### 差异表达必须用 counts 而不是 TPM

TPM 是样本内归一化量，**样本间不可直接比较**（每个样本的 TPM 总和都被拉到 1e6）。DESeq2 需要的是原始整数 counts，再由模型自己估计 size factor。因此本技能输出 `counts.csv`（整数）作为差异分析的输入，TPM/FPKM 仅用于样本间可视化与展示。

### 转录本→基因汇总

用 `tx2gene.tsv` 把同一基因的多个转录本合并：

- **counts**：求和（一个 read 落在该基因任一转录本上都算该基因的 read）
- **TPM**：求和（Salmon 的 TPM 满足 sum(转录本 TPM) = 1e6，汇总到基因后同样满足）
- **FPKM**：转录本水平先算 `(count × 1e9) / (length × total_reads)`，再求和到基因

---

## 3. DESeq2 统计模型

DESeq2 建模每个基因的表达量：

```
K_gi ~ NB(mu_gi, alpha_g)
mu_gi = s_j * q_gi
q_gi = 2^(x_i' * beta_g)
```

其中 `K_gi` 是基因 g 在样本 i 的 count，`s_j` 是样本 j 的 size factor（测序深度校正），`beta_g` 是回归系数，`alpha_g` 是离散度（越大表示样本间波动越强）。

### 3.1 Size factor 归一化

采用 median-of-ratios：

```
s_j = median_g ( K_gj / geomean_g(K_.j) )
```

**为什么不用 TPM**：TPM 用的是长度加权的分母，会让高表达长基因主导归一化。median-of-ratios 对每个基因等权，且对少数极端表达基因稳健（取中位数）。

size factor 应集中在 1 附近。偏离说明测序深度不同（正常）或某样本质量异常（需排查）。

### 3.2 基因离散度估计（gene-wise dispersion）

对每个基因，在离散度网格上最大化 **Cox-Reid 调整后的 profile 对数似然**：

```
logL_CR(alpha) = logL(alpha, beta_hat(alpha)) - 0.5 * log|X'WX|
```

其中 `X'WX` 是 NB-GLM 的 Fisher 信息矩阵。Cox-Reid 调整的作用是扣除参数估计自由度带来的乐观偏差，这是 DESeq2 与简单 MLE 的关键区别。

本实现的做法：
1. NB-GLM 用 IRLS 迭代加权最小二乘拟合（log 链接，size factor 作为 offset）
2. 在 `1e-7` 到 `10` 的对数网格上粗搜，profile 似然取最大
3. 在最优点附近按 3 倍范围细化搜索

### 3.3 离散度趋势

真实数据中离散度随均值系统性下降（低表达基因噪声更大）。DESeq2 的 parametric 家族：

```
alpha_trend(mu) = a0 + a1 / mu
```

本实现在 log 空间拟合该曲线，使用 **Huber 加权迭代最小二乘** 抵抗离群基因，并**剔除落在离散度下限的基因**（它们的 log 值是人为截断的，混入会把趋势拉向无穷大）。

### 3.4 MAP 离散度

离散度先验取 log 空间正态分布：`log alpha ~ N(log alpha_trend, sigma^2)`。`sigma^2` 由趋势拟合残差的 MAD 估计（乘以正态一致性因子 1.4826）。

后验为：

```
log post(alpha) = logL_CR(alpha) - 0.5 * (log alpha - log alpha_trend(mu))^2 / sigma^2
```

**为什么需要收缩**：低表达基因的离散度估计极不稳定（可能因为 1 个 read 的波动就跳到极大值）。向趋势收缩能把这些基因的离散度拉到合理区间，从而得到可用的 p 值。DESeq2 的最终 Wald 检验用的就是这个 MAP 离散度。

### 3.5 Wald 检验

对每个基因在 MAP 离散度下重新拟合 NB-GLM，用 Fisher 信息矩阵近似系数协方差：

```
se^2(beta) = (X'WX)^{-1}
z = (x'beta) / se(x'beta)
p = 2 * (1 - Phi(|z|))
```

`log2FoldChange = x'beta / ln(2)`，即 `group_a` 相对 `group_b` 的对数倍数变化。

### 3.6 BH 校正

Benjamini-Hochberg 控制 FDR：

```
q_(i) = min_{j>=i} ( n / j ) * p_(j)
```

BH 比 Bonferroni 更有功效，因为它允许一部分假阳性以换取更高的真阳性检出率。

### 3.7 双后端对比

| | `--backend python` | `--backend r` |
|---|---|---|
| 实现 | 内置 NumPy/SciPy | Bioconductor DESeq2 |
| 依赖 | numpy, scipy | Rscript + DESeq2 |
| 输出列名 | 与 DESeq2 一致 | 与 DESeq2 一致 |
| 特有能力 | 无 LFC shrinkage | apeGLM/ashr shrinkage、独立过滤 |

两个后端在 log2FC 上高度一致（差异通常 < 0.1），p 值可能因离散度估计细节略有不同。**若两者结论冲突，优先信任 R 后端**（它是官方参考实现）。

---

## 4. 富集分析

### 超几何检验

设背景基因集大小 N，某个 GO term 含 K 个基因，从中随机抽取 n 个（DEG 数量），问至少抽到 k 个属于该 term 的概率：

```
P(X >= k) = hypergeom.sf(k-1, N, K, n)
```

用 `sf`（生存函数）而不是 `cdf`，因为我们关心的是"至少 k 个"这个单侧富集事件。

### 背景集的选择

**这是富集分析最容易出错的地方。** 本实现用注释文件中所有基因作为背景（N = `len(gene_go)`）。更严谨的做法是用"本次实验中所有被检测到的基因"作背景——因为如果某个基因在实验里根本没测到，它不该参与竞争。若 counts 矩阵覆盖的基因数少于注释基因数，应改用实测基因集作背景。

### 基因集大小过滤

- `min_gene_set`（默认 5）：太小的 term 没有统计功效，容易产生大量假阳性。
- `max_gene_set`（默认 500）：太大的 term（如"metabolic process"）几乎必然包含 DEG，缺乏区分度。

### FDR 校正

对 BP / CC / MF 各自做 BH 校正（也可合并后统一校正，更保守）。

---

## 5. 可视化

| 图 | 用途 | 判读 |
|---|---|---|
| 火山图 | x=log2FC，y=-log10(padj) | 显著点呈双翼；翼越陡说明效应量越大 |
| MA 图 | x=log10(baseMean)，y=log2FC | 应平坦无趋势；有下弯说明低表达基因效应量被高估 |
| PCA | 样本间距离 | 组间分离 = 生物学信号；组内混杂 = 批次效应 |
| 热图 | Top DEG 的表达模式 | 应呈现按分组的块状结构 |

所有图同时输出 PDF（矢量，用于投稿）与 PNG（300 DPI，用于预览）。

---

## 6. 参数选择速查

| 场景 | 建议调整 |
|---|---|
| 测序深度低（<20M reads/样本） | `--min-counts 5`，增大 `--min-gene-set` 敏感性检查 |
| 数据质量差 | `--quality-phred 15`，`--min-length 30` |
| 生物学效应小 | `--lfc-threshold 0.5`（但需在文章中声明） |
| 样本数多（>10） | 考虑批次校正或在 metadata 中加入 batch 因子 |
| 只有 2 重复/组 | 统计功效极低，谨慎解读；不要用 `--lfc-threshold 0` |
| 需要发表级 LFC | 用 R 后端配合 `lfcShrink(apeGLM)` |

---

## 参考文献

1. Love MI, Huber W, Anders S. *Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.* Genome Biology, 2014, 15:550.
2. Cox DR. *Inference and asymptotic expansions for the maximum likelihood estimator of a regression coefficient.* Biometrika, 1971, 58(3):583-592.
3. Patro R, et al. *Preferential degradation of differentially expressed mRNAs in the transcriptome of human and yeast.* Nature Methods, 2018, 15:747-751.
4. Anders S, McCarthy DJ, Love MI. *Count-based differential expression analysis with edgeR.* F1000Research, 2014, 3:176.
5. Subramanian A, et al. *Gene set enrichment analysis.* PNAS, 2005, 102(43):15545-15550.
6. Chen H, et al. *fastp: an ultra-fast all-in-one FASTQ preprocessor.* Bioinformatics, 2018, 34(17):i884-i890.
7. Patro R, et al. *Salmon provides fast and accurate quantification of transcript expression.* Nature Methods, 2018, 15:201-205.
