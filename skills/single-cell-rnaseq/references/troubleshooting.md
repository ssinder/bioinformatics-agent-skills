# 单细胞 RNA-seq 常见问题与排错手册

按症状组织。**先用症状索引定位，再看对应条目。**

---

## 快速索引

| 症状 | 跳转 |
|---|---|
| 找不到输入文件 / 格式无法识别 | [1. 输入读取问题](#1-输入读取问题) |
| 质控后细胞数几乎为 0 | [2. 质控过滤过度](#2-质控过滤过度) |
| 质控几乎没过滤掉任何细胞 | [3. 质控过滤不足](#3-质控过滤不足) |
| `ValueError: 质控后仅剩 N 个细胞` | [2. 质控过滤过度](#2-质控过滤过度) |
| UMAP 上所有点挤成一团 | [4. 降维塌陷](#4-降维塌陷) |
| UMAP 上出现极端离群细胞 | [5. 离群点](#5-离群点) |
| 簇之间混在一起分不开 | [6. 簇不分离](#6-簇不分离) |
| 簇数远多于/少于预期 | [7. 聚类粒度问题](#7-聚类粒度问题) |
| 出现大量极小簇（< 20 细胞） | [7. 聚类粒度问题](#7-聚类粒度问题) |
| marker 基因全是核糖体蛋白 | [8. marker 无意义](#8-marker-无意义) |
| 注释结果全是 `Unknown` | [9. 注释失败](#9-注释失败) |
| 报错 `No module named 'scanpy'` | [10. 依赖缺失](#10-依赖缺失) |
| 报错 `Leiden 不可用` / 回退 Louvain | [10. 依赖缺失](#10-依赖缺失) |
| 报错 `scales not aligned` / 数组长度不匹配 | [11. 代码缺陷或数据不一致](#11-代码缺陷或数据不一致) |
| 结果不可复现 | [12. 可重复性](#12-可重复性) |
| 图上出现中文乱码 | [13. 其他](#13-其他) |

---

## 1. 输入读取问题

### 1.1 `输入路径不存在`

```
SkillError: 输入路径不存在: /data/xxx
```

检查路径拼写与挂载点。注意容器内路径与宿主机路径不同。

### 1.2 `无法识别的输入格式`

`read_matrix()` 支持四种格式，按以下顺序嗅探：

1. 单个 `.h5` 文件 → 10x HDF5
2. 含 `matrix.mtx` / `matrix.mtx.gz` 的目录 → 10x MTX
3. `.csv` / `.tsv` / `.txt` / `.tab` → 稠密表格
4. `.h5ad` → AnnData

若目录结构不标准，典型情况：

- 10x v3 解压后目录名为 `filtered_feature_bc_matrix/`，内含 `barcodes.tsv.gz`、`features.tsv.gz`、`matrix.mtx.gz`。**注意**：特征文件叫 `features.tsv.gz`（v3）或 `genes.tsv.gz`（v2），两者都支持。
- 若只有 `raw_feature_bc_matrix/`（未过滤），细胞数会明显偏多，属于正常现象。

**排查**：

```bash
ls -R /data/filtered_feature_bc_matrix
```

### 1.3 读进来维度反了

细胞数应该远大于基因数。若发现 `n_cells < n_genes`，说明矩阵转置了。用稠密 CSV 时常见：很多工具导出的是"基因 × 细胞"。

**处理**：把 CSV 转置后再读，或改用 10x 原始格式。

### 1.4 HDF5 读取报 `indices and data should have the same size`

说明输入的 HDF5 不是 CellRanger 标准格式（可能已被 `scipy.io.savemat` 之类工具改写）。

**处理**：从原始 `filtered_feature_bc_matrix` 重新导出，或改用 MTX 目录。

---

## 2. 质控过滤过度

**症状**：过滤后细胞数骤降（< 30% 保留率），或直接触发 `ValueError: 质控后仅剩 N 个细胞`。

**排查顺序**：

1. **看 `01_qc/qc_summary.json` 里的 `n_cells_before`**，确认读进来的细胞总数是否合理。若只有几百个细胞但数据是 10x 完整输出，说明读取路径选错了（可能读到了 filtered 目录里的某个子集）。
2. **看 `qc_metrics.csv` 的 `n_genes_by_counts` 分布**。如果中位数只有 50 而不是 300+，那是数据本身质量问题（解离差或测序浅），不是阈值问题——放宽阈值也救不回来。
3. **如果中位数正常（> 500）但过滤率仍 > 60%**，说明 `--max-mt-pct` 太严。检查 `pct_counts_mt` 分布。

**处理**：

```bash
--min-genes 100 --min-counts 250 --max-mt-pct 25
```

10x v2 化学体系应使用更低的 `min_genes`（100 左右）。

---

## 3. 质控过滤不足

**症状**：过滤率 < 2%，但 UMAP 上有明显的低质量细胞团。

**典型原因**：该组织天然线粒体含量高（棕色脂肪、肾小管、心肌），`max_mt_pct=20` 形同虚设。

**处理**：

1. 看 `qc_scatter.png` 中 `total_counts` vs `n_genes_by_counts` 的右上角区域是否有离群点群（双细胞）。
2. 适度收紧 `--min-counts`（例如从 100 提到 2000）比动 `max_mt_pct` 更安全。
3. 若确实存在一批双细胞，本技能内置的质控不能可靠去除，需要单独跑 Scrublet / DoubletFinder 后合并标注。

---

## 4. 降维塌陷

**症状**：UMAP 图上所有细胞重叠成一个点。

**可能原因与处理**：

| 原因 | 判断方法 | 处理 |
|---|---|---|
| 高变基因太少 | `02_normalized/var_hvg.csv` 中 `highly_variable=True` 的数量很少 | 提高 `--n-top-genes`（如 500 → 2000） |
| PCA 维度不足 | `pca_variance.csv` 只有 1–2 行 | 提高 `--n-comps` |
| 数据本身无结构 | PCA 前 5 PC 解释方差 < 10% | 检查质控与归一化，或数据确实高度同质 |
| `--n-neighbors` 太大 | 邻接图退化成全连接 | 降到 10–15 |

**注意**：使用内置降级后端时，UMAP 的力优化是**在谱初始化之上**做的微调。如果谱初始化本身是好的（簇结构清晰），力优化不会塌陷。如果发现塌陷，先确认 `--backend numpy` 与 `--backend scanpy` 的差异；若 scanpy 正常而 numpy 塌陷，说明是降级实现的问题。

---

## 5. 离群点

**症状**：UMAP 上有几个细胞被甩到很远的位置。

**排查**：

1. 这些细胞的 `total_counts` 与 `n_genes_by_counts` 是否极端（双细胞特征）？是则收紧质控。
2. 是否某个基因在少数细胞中表达量极高（会造成 1 个主成分被主导）？看 `var_hvg.csv` 的 `means` 极值，并考虑降低 `--n-top-genes` 或检查该基因是否为批次相关。
3. 是否真的只有几个细胞？少于 10 个的孤立细胞对聚类影响有限，可忽略。

**处理**：收紧质控过滤掉；或在热图里忽略它们。

---

## 6. 簇不分离

**症状**：UMAP 上群之间界限模糊、多个群连成一片。

**原因与处理**：

- **分辨率太低**：把 `--resolution` 从 1.0 提到 1.5–2.0 重跑。
- **簇内确实没有足够差异**：合成数据显示连续谱系（如单核细胞向单核细胞过渡态）在 UMAP 上本就连成一条带，这是真实生物学。
- **细胞类型差异主要在非高变基因上**：增加 `--n-top-genes` 或减少（避免引入噪声），一般增加有效。
- **数据被过度归一化**：CP10K 把测序深度差异全消掉，但也压掉了真实丰度差异。检查 `--target-sum` 是否合理。

**关键检查**：聚类结果与 UMAP 是两件事。**聚类只看邻接图**，不看 UMAP 坐标。若邻接图（`04_clustered/`）中的簇是干净的，只是 UMAP 显示重叠，那是可视化问题，不是聚类问题。

---

## 7. 聚类粒度问题

### 7.1 簇太多

扫描 `--resolution` 0.5 → 1.0 → 2.0，找到 marker 稳定的分辨率。同时检查：

- 是否把不该合并的不同细胞类型强行聚在一起（分辨率不足）。
- 是否存在批次效应被当成生物学差异（分辨率过高会把每个批次切成独立簇）。

**批次效应排查**：把批次信息写进 `obs`，在 UMAP 上按批次着色，如果簇与批次一一对应，说明聚类捕获的是批次而非细胞类型。可用 Harmony / BBKNN 校正后再聚类（本技能未内置）。

### 7.2 簇太少

- 降低 `--resolution`。
- 提高 `--n-neighbors`（更平滑的图 → 更容易合并）。
- 检查是否 HVG 太少导致信息丢失。

### 7.3 大量极小簇

出现很多 < 20 细胞的小簇，通常是：

- 双细胞；
- 极低质量细胞（空滴）；
- 噪声维度过高（HVG 里混入了环境 RNA）。

处理：收紧质控、降低 `--n-top-genes`、略微提高 `--resolution`（会让小簇更碎）或降低分辨率（会吞并小簇）。

---

## 8. marker 无意义

**症状**：某簇的 top marker 全是 `RPL*` / `RPS*`（核糖体蛋白）或 `MALAT1`、`MIR*`（长非编码）。

**原因**：这类基因要么 ubiquitous、要么与细胞大小相关而非细胞身份相关。当簇内缺乏区分信号时，排序算法只能给出这些噪声基因。

**处理**：

1. 提高 `--resolution`，让该簇被拆细——亚群一旦分开，真正的 marker 就会出现。
2. 检查是否 `n_top_genes` 过大：HVG 里混入了低信息量的高变基因。
3. 在 marker 计算时排除核糖体蛋白与高表达 housekeeping 基因（本技能默认保留，可用 `rank_genes_groups` 的输入自行过滤）。

---

## 9. 注释失败

**症状**：所有 cluster 都标为 `Unknown`，或 `cluster_annotations.csv` 里 `overlap` 全是 0。

**原因**：

1. **marker 参考文件的基因 ID 空间与数据不一致**（最常见）。例如数据用 `ENSG00000141510`，参考文件用 `Gene1` 符号，或一个带版本号一个不带。
2. **marker 列名不匹配**：参考文件的 marker 列必须叫 `subclass.markers.combo`、`markers`、`marker_genes` 或 `markers.combo` 之一。列名不对会打 warning 并跳过注释。
3. **基因符号大小写 / 别名差异**（`CD3D` vs `Cd3d`）。

**处理**：

```bash
head -2 /ref/cell_marker_reference.csv
```

确认列名与分隔符。然后统一 ID 空间——用 `strip_version()` 去版本号，或直接用 gene symbol 做注释。

**另外**：参考文件需要至少两列标签（`class_label` 与 `subclass_label`）。只有一列时，代码会退化取第 1、2 列作为标签，可能不符合预期。

---

## 10. 依赖缺失

### 10.1 scanpy 不可用

```
[INFO] scanpy 不可用，使用 numpy 降级后端
```

这不是错误。降级后端实现了全部核心算法：CP10K、log1p、Seurat 风格 HVG、截断 SVD PCA、Jaccard kNN 图、Louvain 聚类、UMAP 谱初始化 + 力优化、t-SNE、Welch t 检验 marker。

**差异**：
- 降级后端用 Louvain 而非 Leiden（除非装了 `leidenalg`）。
- 降级后端的 UMAP 是近似实现，簇内精细结构可能与官方 UMAP 略有差异，但簇的划分一般一致。

要装 scanpy：

```bash
pip install scanpy anndata
```

### 10.2 leidenalg 不可用

```
[WARNING] leidenalg/igraph 不可用，回退到 Louvain
```

```bash
pip install leidenalg igraph
```

Louvain 与 Leiden 在大多数数据上给出几乎相同的划分。

### 10.3 umap-learn 不可用

自动使用内置降级 UMAP，功能完整。若要官方实现：

```bash
pip install umap-learn
```

### 10.4 matplotlib 不可用

图表步骤被跳过，其他产出正常。要出图：

```bash
pip install matplotlib
```

---

## 11. 代码缺陷或数据不一致

### 11.1 `operands could not be broadcast together with shapes (N,2) (N,)`

这是 UMAP 排斥力计算中的维度错误：向量化张量必须是 `(n, n_neg, d)`，用 `np.where` 把 `(n, d)` 降到 `(n,)` 就会触发。当前实现已使用显式的多负样本张量广播；若你修改了这段代码，请保持三轴结构。

### 11.2 `Length of values (4) does not match length of index (600)`

这是 Louvain 多层聚合时的细胞→节点映射错误：聚合层结束后必须返回**每个原始细胞**的标签，而不是社区数量。当前实现用 `cell_to_node` 数组显式追踪映射。若你改写了 `cluster_louvain`，请确认最终返回 `labels = node_labels[cell_to_node]`。

### 11.3 `Length of values (4) does not match length of index` 出现在 PCA 之后

检查 `X[:, hvg_mask]` 的布尔掩码长度是否等于 `X.shape[1]`。若 `highly_variable_genes` 返回的是"基因名列表"而不是"布尔掩码"，索引会错位。

### 11.4 UMAP 结果与聚类结果不匹配

**这是正常现象，不是 bug**。UMAP 是二维投影，簇在投影中可能被压缩或旋转。判断聚类是否成功要看 `04_clustered/` 里的标签与 `05_annotated/marker_genes.csv`，而不是看 UMAP 坐标。

---

## 12. 可重复性

- 始终显式传 `--seed`（默认 0）。
- `svds`（PCA）、负样本采样（UMAP）、t-SNE 都会受随机性影响。
- **同一份数据 + 同一个 `--seed` + 同一个 `--backend` 必须得到完全一致的聚类标签。** 若不一致，检查是否在中间步骤引入了未固定的随机性（例如自己调用了 `np.random.rand` 而没有用 `np.random.default_rng(seed)`）。
- 后端之间不应要求逐点一致：`scanpy` 与 `numpy` 的簇编号顺序可能不同（编号本身是任意的），应比较簇的大小分布与 marker 集合。

---

## 13. 其他

### 13.1 图上中文乱码

matplotlib 默认字体不含 CJK。设置环境变量或修改脚本里的字体配置：

```python
import matplotlib
matplotlib.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "Noto Sans CJK SC"]
matplotlib.rcParams["axes.unicode_minus"] = False
```

### 13.2 内存不足

10x 5 万细胞的数据，稠密化是最主要的内存杀手。

- `scale()` 会把矩阵转成稠密：`n_cells × n_genes × 8` 字节。5 万 × 2 万基因 = 8 GB。建议先减小 `--n-top-genes`。
- `read_dense_table` 一次性读入整个 CSV，10x 规模请务必用 MTX/HDF5。
- 可加 `--no-tsne` 节省大量内存与时间（t-SNE 是 $O(n^2)$）。

### 13.3 运行很慢

耗时排序（5 万细胞量级）：t-SNE > UMAP > 邻接图 > PCA > 聚类。依次尝试 `--no-tsne`、减小 `--n-comps`、减小 `--n-top-genes`。

### 13.4 何时该怀疑是流程 bug 而不是数据问题

按以下顺序自检：

1. 用 `quality_control.py` 单独跑一遍质控，看 `qc_summary.json` 的数字是否合理。
2. 检查 `02_normalized/var_hvg.csv` 里被选中的基因是否确实高变（`means` 与离散度是否单调）。
3. 检查 `03_dimensionality/pca_variance.csv` 的前 5 个 PC 是否解释了大量方差（若 < 5%，说明输入矩阵有问题）。
4. 检查 `05_annotated/marker_genes.csv` 的 `padj` 分布：应绝大多数 > 0.05，若几乎全显著，说明分组或矩阵有问题。
5. 对照 `qc_metrics.csv` 抽查几个离群细胞的原始计数，确认不是读取错位。

若以上都正常但结论仍不合常理，优先怀疑数据本身（解离质量、批次、样本标签）。

---

## 附：常用排查命令

```bash
# 质控统计
cat 01_qc/qc_summary.json

# 每细胞指标（快速看分布）
python -c "import pandas as pd; d=pd.read_csv('01_qc/obs_qc.csv', index_col=0); print(d.describe())"

# 主成分解释方差
cat 03_dimensionality/pca_variance.csv

# 每簇细胞数与 marker 数
cat 05_annotated/cluster_summary.csv

# 每簇 top 10 marker
python -c "import pandas as pd; d=pd.read_csv('05_annotated/marker_genes.csv'); print(d.groupby('cluster').head(10)[['cluster','gene','log2FC','padj']])"

# marker 的 padj 分布（应绝大多数 > 0.05）
python -c "import pandas as pd; d=pd.read_csv('05_annotated/marker_genes.csv'); print(d['padj'].describe())"
```
