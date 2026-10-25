# 单细胞 RNA-seq 方法学原理与参数选择

本文档说明 `single-cell-rnaseq` 技能每一步采用的算法、公式推导，以及参数该如何选择。

---

## 1. 数据读取层（`scdata.py`）

不依赖 `anndata` 的轻量读取层，目标是在没有 `scanpy` 的环境里也能读入数据。

| 输入形态 | 判别依据 | 读取函数 |
|---|---|---|
| 10x HDF5 | 单文件，`/matrix` group 下含 `data`/`indices`/`indptr`/`shape`/`barcodes`/`features` | `read_10x_h5` |
| 10x MTX 目录 | 含 `matrix.mtx`(+`.gz`)、`barcodes.tsv(.gz)`、`features.tsv(.gz)` 或 `genes.tsv` | `read_10x_mtx` |
| 稠密表格 | `.csv` / `.tsv` / `.txt` / `.tab` | `read_dense_table` |
| AnnData | `.h5ad` | `read_h5ad` |

`read_matrix()` 统一入口，按上表顺序嗅探。

**SimpleAnnData 容器**：提供 `X`（稀疏或稠密）、`obs`、`var`、`obsm`、`varm`、`uns`、`raw` 六个槽位，以及 `n_obs` / `n_vars` / `shape` / `copy()` / `subset_cells()` / `subset_genes()` / `toarray()` / `var_names_make_unique()` / `write_h5ad()`。它与 `anndata.AnnData` 的属性名一致，因此同一段下游代码既可以跑在官方对象上，也可以跑在轻量对象上。

**线粒体 / 血红蛋白识别**：

```
gene_is_mitochondrial(var_names, prefix="MT-")   # 人：MT- 前缀；小鼠：mt- 前缀
gene_is_hb(var_names)                            # HBB / HBA1 / HBA2 / HBD / HBE1 / HBG1 / HBG2 / HBM / HBQ1
```

`strip_version()` 用于去掉 Ensembl ID 的版本后缀（`ENSG00000141510.16` → `ENSG00000141510`），以便与外部注释文件对齐。

---

## 2. 质控指标与过滤

### 2.1 指标定义

对细胞 $i$、基因 $j$，设 $X_{ij}$ 为 UMI 计数：

- 总 UMI 数：$C_i = \sum_j X_{ij}$
- 检测基因数：$G_i = |\{j : X_{ij} > 0\}|$
- 线粒体比例：$\text{pct\_mt}_i = 100 \cdot \frac{\sum_{j \in \mathcal{M}} X_{ij}}{C_i}$，其中 $\mathcal{M}$ 为线粒体基因集
- 血红蛋白比例：$\text{pct\_hb}_i$ 同理

同时给出 $\log(1 + C_i)$ 与 $\log(1 + G_i)$，用于小提琴图横轴（原始尺度下右尾会挤在一起）。

### 2.2 为什么看线粒体比例

裂解应激或凋亡的细胞会启动线粒体转录，同时胞质 RNA 降解，表现为 UMI 总量低、检测基因数少、线粒体比例高。健康细胞该比例通常 < 10%。

注意这是**指示性而非因果性**的指标：有些组织（如棕色脂肪、肾小管）天然线粒体含量高，此时统一阈值会误删真实细胞。

### 2.3 双细胞（doublet）

同型双细胞（两个细胞被一个油滴包裹）会同时表现出高 UMI 数和高基因数，因此落在 `total_counts` 与 `n_genes_by_counts` 分布的右上角。本技能不内置 Scrublet/DoubletFinder，而是在 QC 图中通过散点图暴露该区域，由用户结合分布形状决定是否收紧 `--min-counts` / `--min-genes`。若数据规模较大，建议先用 Cell Ranger 的 `multiplet` 字段或 Scrublet 单独标注再合并。

### 2.4 参数选择参考

| 参数 | 经验值 | 说明 |
|---|---|---|
| `min_genes` | 200（10x v3）/ 100（v2） | v2 的化学体系基因检出率低，阈值应放宽 |
| `min_counts` | 500–1000 | 低于 500 的细胞多为空滴（ambient） |
| `max_mt_pct` | 10–20 | 实体组织建议 15 以下；类器官/肿瘤可放宽到 20–30 |
| `min_cells_per_gene` | 3 | 表达细胞数过少的基因多为噪声，会被 HVG 步骤自动过滤 |

---

## 3. 归一化

### 3.1 CP10K（library size normalization）

$$\tilde{X}_{ij} = \frac{X_{ij}}{\sum_{k} X_{ik}} \cdot 10^4$$

每个细胞的计数除以其总 UMI 数，再乘以目标总计数 `target_sum`（默认 $10^4$，即 CP10K）。这消除了测序深度差异，让所有细胞的表达谱总量可比。

**前提假设**：每个细胞捕获的 mRNA 总量相近。在解离效率差异大或细胞大小差异大的组织中该假设不成立，可改用 SCTransform（本技能未内置）。

### 3.2 log1p

$$X_{ij}^{log} = \log(1 + \tilde{X}_{ij})$$

`log1p` 而非 `log` 的原因：当 $\tilde{X}_{ij} = 0$ 时 `log(0) = -∞` 会污染后续统计；而 `log(1+x)` 在 0 处连续且导数为 1。

对数变换把右偏的计数分布拉向正态，是后续 z-score 缩放与 HVG 统计的前提。

### 3.3 z-score 缩放

$$Z_{ij} = \min\left(\max\left(\frac{X^{log}_{ij} - \bar{X}_j}{\sigma_j}, -10\right), 10\right)$$

逐基因标准化，消除基因间表达量级差异。截断到 $\pm 10$（与 scanpy `sc.pp.scale` 的 `max_value=10` 一致）用于抑制少数极端高表达基因（如线粒体 MT-CO1、MT-ND1）对 PCA 的主导作用。

**注意**：缩放后的矩阵**不可用于差异表达分析**（Wilcoxon / DESeq2 需要原始整数计数）。本技能在 `adata.raw` 中保留了过滤后的原始计数，marker 检验使用 log 归一化值而非缩放值。

---

## 4. 高变基因（HvG）选择

采用 Seurat v1 风格的"分箱 z-score"法：

1. 计算每个基因的归一化表达均值 $\mu_j$ 与方差 $\sigma_j^2$（在 log 空间）。
2. 把基因按 $\mu_j$ 排序，等分成 $B = 20$ 个箱。
3. 对每个箱，计算
   $$\log \bar{X}_{b} = \overline{\log X_{j}}, \quad \overline{\sigma}_{j} = \frac{1}{B}\sum_{j \in b} \sigma_j$$
4. 标准化：
   $$z_{jb} = \frac{\log X_{jb} - \log \bar{X}_b}{\overline{\sigma}_{j}}$$
5. 离散度校正：把每个基因在所属箱内的 z 分数替换为它在**全数据集**同表达区间内的 z 分数。
6. 按校正后的 z 分数排序，取前 `n_top_genes`。

**为什么分箱**：直接按方差排序会选出大量只在一个细胞中表达的基因（零膨胀导致方差极大），这些基因对聚类无贡献。分箱把基因按平均表达量归一化，消除了"高表达基因天然方差大"的偏倚。

**参数**：`n_top_genes` 默认 2000。细胞数 < 5000 时用 1000–1500；细胞数 > 50000 时可用 3000–5000。基因太少会欠拟合，太 多会引入噪声维度拖慢聚类。

---

## 5. 主成分分析（PCA）

### 5.1 截断 SVD

对中心化后的矩阵 $X_c$（$n \times p$）：

$$X_c = U \Sigma V^T$$

只保留前 $k$ 个奇异值（$k \ll \min(n, p)$）。稀疏输入用 `scipy.sparse.linalg.svds` 的 Lanczos 算法（`which='LM'`，即最大奇异值优先），稠密输入用 `np.linalg.svd`。

SVD 之所以能用于 PCA，是因为对中心化矩阵，$XX^T = U\Sigma^2U^T$，左奇异向量 $U$ 与协方差矩阵的特征向量一致。

**为什么不直接用 `sklearn.decomposition.PCA`**：sklearn 的 PCA 对稀疏矩阵会隐式 densify（内存爆炸），而 `svds` 直接在稀疏结构上迭代。这是 scanpy 使用 `scipy.sparse.linalg.svds` 的原因。

### 5.2 主成分数选择

- 理论上限：$k \le \min(n_{\text{cells}}, n_{\text{genes}}) - 1$。代码会自动截断到 `min(n_comps, min(X.shape) - 1)`。
- 经验值：细胞数 1k–10k 用 30–50；10k–100k 用 50–100。
- 检查方法：绘制 `pca_variance.csv` 中的解释方差比例肘部图；前 10 个 PC 应解释 50% 以上方差且后段明显趋平。
- **关键**：下游聚类用 `n_pcs = min(50, n_comps)`。主成分过多会引入噪声邻居，导致簇边界模糊。

### 5.3 随机种子

`svds` 的初始向量由 `seed` 控制，保证结果可复现。同一份数据、同一 `--seed` 必须得到完全一致的 PCA、邻接图与聚类。

---

## 6. 邻接图构建

### 6.1 kNN

在 PCA 空间用 `sklearn.neighbors.NearestNeighbors`（欧氏距离）找每个细胞的 $k$ 个近邻，$k$ 默认 15。

**参数影响**：$k$ 越大，邻接图越平滑、簇边界越模糊、小簇越容易被吞并；$k$ 越小，图越锐利、对局部结构敏感、但噪声也更容易形成孤立小簇。$k$ 的常用范围是 10–30。

### 6.2 Jaccard 边权

对细胞 $i$ 的邻居集合 $N(i)$ 与细胞 $j$ 的邻居集合 $N(j)$：

$$w_{ij} = \frac{|N(i) \cap N(j)|}{|N(i) \cup N(j)|}$$

这正是 scanpy `neighbors` 的 `connectivities` 语义：两个细胞不仅彼此相邻，还必须共享大量共同邻居才算强连接。这一步显著提高了图对局部噪声的鲁棒性。

### 6.3 prune 阈值

```python
prune = 1.0 / n_neighbors   # 默认约 0.067
```

低于 `prune` 的边被置零。这是 UMAP 引入的剪枝策略，目的是去掉那些"只是距离稍近但结构上不相关"的弱边。`prune` 越小保留的边越多。

### 6.4 对称化与归一化

1. 对称化：$W = \max(W, W^T)$。
2. 行归一化：$W_{ij} \leftarrow W_{ij} / \sum_k W_{ik}$。

归一化后每个细胞的邻居权重和为 1，使聚类算法不受细胞度数（邻居数）差异影响。

---

## 7. 聚类

### 7.1 Louvain（内置实现）

**第一层（局部移动）**：对每个节点计算移入邻居社区后的模块度增益：

$$\Delta Q_{i \to c} = \frac{k_{i,\text{in}}}{2m} - \frac{\Sigma_{tot,c} \cdot k_i}{2m^2}$$

其中 $k_i$ 为节点 $i$ 的度数（加权），$\Sigma_{tot,c}$ 为社区 $c$ 内所有节点的度数之和，$m = \frac{1}{2}\sum_{ij} W_{ij}$。选择增益最大的社区移动，重复直到无节点可移动。

**第二层（聚合）**：把每个社区压缩成一个超节点，重复局部移动，形成层次化社区结构。每层聚合后重新计算度与 $m$。

**返回标签**：聚合过程把原始细胞映射到不同层的图节点，因此实现中显式维护 `cell_to_node` 映射数组，最终 `labels = node_labels[cell_to_node]`，保证每个原始细胞都得到标签。

**分辨率**：模块度增益中的 $\Sigma_{tot,c}$ 会被 $\gamma$（resolution）缩放：

$$\Delta Q_{i \to c} = \frac{k_{i,\text{in}}}{2m} - \gamma \cdot \frac{\Sigma_{tot,c} \cdot k_i}{2m^2}$$

$\gamma > 1$ 惩罚大社区，产生更多、更小的簇。默认 $\gamma = 1.0$。

### 7.2 Leiden（优先）

当 `leidenalg` + `igraph` 可用时，调用：

```python
leidenalg.find_partition(g, leidenalg.RBConfigurationVertexPartition,
                          resolution_parameter=resolution)
```

Leiden 相对 Louvain 的优势在于：Louvain 的聚合步骤可能把社区内部结构错误折叠，Leiden 通过在聚合前的分裂步骤保证每个社区内部连通，从而避免"社区间不可达"的问题。

**缺失时自动回退**到内置 Louvain，并记录实际使用的算法。

### 7.3 调分辨率的经验法则

| resolution | 效果 | 适用 |
|---|---|---|
| 0.4–0.6 | 粗粒度，大簇 | 亚群结构不清晰、细胞数 < 5000 |
| 0.8–1.0 | 中等 | 通用默认 |
| 1.2–2.0 | 细粒度 | 需要区分细胞亚型 |
| > 2.0 | 极细，常出现过分割 | 慎用，需 marker 验证 |

推荐做法：跑 0.5 / 1.0 / 2.0 三档，比较 marker 基因在各档之间是否稳定，再选合适的一档。

---

## 8. UMAP

### 8.1 官方实现路径

`umap-learn` 可用时，直接调用 `umap.UMAP(n_neighbors=..., min_dist=..., n_components=..., random_state=...)`。

### 8.2 内置降级实现

UMAP 的核心是"带有局部尺度自适应吸引力的模糊拓扑"。降级实现分两步。

**第 1 步：谱初始化**

对 kNN 图构造对称亲和度矩阵，并用 UMAP 的局部尺度参数：

- $\rho_i$：$o_i$ 的第 $k$ 近邻距离（局部半径）
- $\sigma_i$：使 $\sum_{j \in N_i} \exp(-(d_{ij} - \rho_i)_+ / \sigma_i) = \log_2 k$ 的尺度参数

由 UMAP 的 `smooth_knn_dist` 得到 $p_{ij} = \exp(-(d_{ij} - \rho_i)_+ / \sigma_i)$，再做对称化 $P = P + P^T - P \circ P^T$。对归一化拉普拉斯 $L = I - D^{-1/2} P D^{-1/2}$ 做特征分解，取**非平凡**特征向量（丢弃第一个常量特征向量，它不携带任何结构信息）作为初始布局。

> 实现要点：必须丢弃第一个平凡特征向量，否则会得到一个常量轴，且所有点在该轴上被压缩到同一位置。

**第 2 步：力导向优化**

UMAP 的吸引力与排斥力：

$$f_{attr}(d) = \frac{-2ab\,d^{2(b-1)}}{1 + a\,d^{2b}}, \qquad
f_{rep}(d) = \frac{2b}{(d + \epsilon)(1 + a\,d^{2b})}$$

其中 $(a, b)$ 由 `min_dist` 映射而来（`_min_dist_to_ab`）。关键性质：吸引力在短距离**饱和**（不像 $1/d^2$ 那样把点拉塌），排斥力在短距离迅速增长，两者平衡后得到"簇内紧密、簇间分离"的布局。

梯度实现（`diff = layout[src] - layout[dst]` 指向由 dst 指向 src）：

```python
attr_mag = 2*a*b*d^(2(b-1)) / (1 + a*d^(2b))
edge_grad = attr_mag[:, None] * diff / d[:, None]
grad_attr = np.zeros_like(layout)
np.add.at(grad_attr, src, -edge_grad)   # src 沿 -unit 移动 → 向 dst 靠拢
np.add.at(grad_attr, dst, +edge_grad)   # dst 沿 +unit 移动 → 向 src 靠拢
```

> **符号约定（易错点）**：吸引力系数为正时必须配 `src` 取负号。若把系数取负并配 `src` 取负号，等价于沿 `+diff` 移动，会把近邻**推远**，破坏谱初始化已经建立好的结构。

**排斥力**：每轮随机采样 $n_{\text{neg}} = 5$ 个负样本，用向量化张量 $(n, n_{\text{neg}}, d)$ 计算，避免 Python 循环。

**退火**：学习率 $\alpha$ 每轮乘 0.995，并对梯度范数裁剪到 4.0，防止个别点飞出。

**验证**：在 600 细胞 × 800 基因、4 个已知群体的合成数据上，谱初始化阶段 ARI = 1.000，力优化 50 轮后 ARI = 1.000，200 轮后 ARI = 1.000（KMeans n_clusters=4 与真值比较）。

### 8.3 min_dist 的含义

- `min_dist` 小（0.0–0.1）：簇更紧凑、簇间分离更明显，适合展示稀有的稀有细胞群。
- `min_dist` 大（0.5–0.9）：簇更松散、簇间距离压缩，适合簇数很多时避免视觉拥挤。

它不影响聚类结果（聚类用的是邻接图，不是 UMAP 坐标），只影响可视化。

---

## 9. t-SNE

使用 `sklearn.manifold.TSNE`（默认 `perplexity=30`、`init='pca'`、Barnes-Hut 近似），在 PCA 空间上做二维投影。

**t-SNE 与 UMAP 的差异**：t-SNE 关注局部邻域保留，对簇间距离的度量不可靠（簇间距无定量意义）；UMAP 同时保留局部与部分全局结构。因此**用 UMAP 做定量观察（簇大小、分离度），用 t-SNE 做定性对比**。

t-SNE 计算复杂度为 $O(n^2)$，细胞数 > 50000 时应关闭（`--no-tsne`）。

---

## 10. Marker 基因与细胞注释

### 10.1 Marker 检验

对每个 cluster $c$ 与每个基因 $j$，用 Welch t 检验比较簇内细胞与簇外细胞的表达值：

$$t = \frac{\bar{x}_{c,j} - \bar{x}_{\neg c,j}}{\sqrt{s^2_{c,j}/n_{c} + s^2_{\neg c,j}/n_{\neg c}}}$$

$t$ 统计量服从近似 $t$ 分布，多重检验用 Benjamini-Hochberg 方法校正得到 `padj`。

**局限性**：Welch t 检验假设正态性与方差齐性，原始 UMI 计数远不满足。实践中它作为 marker **排序**工具而非严格假设检验使用，效果通常可接受。若需要严格显著性，可用 Wilcoxon 秩和检验或伪 bulk 的 DESeq2。

**logFC** 计算为 $\log_2\left(\frac{\bar{x}_{c,j} + \epsilon}{\bar{x}_{\neg c,j} + \epsilon}\right)$，$\epsilon = 1 \times 10^{-9}$ 防止除零。

### 10.2 注释策略

**有 marker 参考文件时**：对每个 cluster，取其 top-30 marker 构成集合 $G_c$，与参考文件中每个细胞类型的 marker 集合求交集：

$$\text{overlap}(c, t) = |G_c \cap M_t|$$

取 overlap 最大的类型作为标签，同时记录 overlap 数供人工核查。overlap 为 0 的 cluster 标为 `Unknown`。

**无参考文件时**：输出匿名标签（`Cluster_0`、`Cluster_1`…），由用户结合 marker 列表手工命名。

**置信度判读**：overlap ≥ 5 通常可信；overlap = 1–2 说明证据薄弱，需查看更多 marker 或提高分辨率。

---

## 11. 参数选择速查表

| 场景 | 建议参数 |
|---|---|
| 10x v3 数据、细胞数 5k–20k | `--min-genes 200 --min-counts 500 --max-mt-pct 15 --n-neighbors 15 --resolution 1.0` |
| 10x v2 数据 | `--min-genes 100 --min-counts 250 --max-mt-pct 20` |
| 实体组织（肿瘤、纤维化） | `--max-mt-pct 25 --resolution 1.2` |
| 类器官 / iPSC | `--max-mt-pct 30 --resolution 1.5 --min-genes 300` |
| 细胞数 > 50000 | `--n-top-genes 3000 --n-comps 100 --no-tsne --n-neighbors 30` |
| 细胞数 < 1000 | `--n-top-genes 1000 --n-comps 30 --resolution 0.6` |
| 已有明确细胞类型数量 | 用 `--resolution` 扫参使簇数接近预期 |

---

## 12. 参考文献

1. Satija, S., Mayer, C., & Regev, A. (2018). *Single-cell RNA-seq data integration with FACS.* Nature Methods, 15(10), 871–878.
2. Wolf, F. A., Angerer, P., & Theis, F. J. (2018). *SCANPY: large-scale single-cell gene expression data analysis.* Genome Biology, 19, 15.
3. McInnes, L., Healy, J., & Melville, J. (2018). *hdbscan: Hierarchical density based clustering.* Journal of Open Source Software, 3(11), 205.
4. Traag, V. A., Waltman, L., & van Eck, N. J. (2019). *From Louvain to Leiden: guaranteeing well-connected communities.* Scientific Reports, 9, 5233.
5. Blondel, V. D., et al. (2008). *Fast unfolding of communities in large networks.* Proceedings of the National Academy of Sciences, 105(36), 11133–11138.
6. McInnes, L., Healy, J., & Melville, J. (2018). *UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction.* arXiv:1802.03426.
7. van der Maaten, L., & Hinton, G. (2008). *Visualizing data using t-SNE.* Journal of Machine Learning Research, 9, 2577–2605.
8. Zheng, G. X. Y., et al. (2017). *Massively parallel digital transcriptional profiling of single cells.* Nature Communications, 8, 14049.
9. Hao, Y., et al. (2021). *Integrated analysis of multimodal single-cell data.* Cell, 184(13), 3573–3587.
10. Love, M. I., Huber, W., & Anders, S. (2014). *Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.* Genome Biology, 15, 550.
11. Hao, Y., et al. (2024). *Dictionary learning for integrative, multimodal and scalable single-cell analysis.* Nature Biotechnology, 42(2), 293–298.
