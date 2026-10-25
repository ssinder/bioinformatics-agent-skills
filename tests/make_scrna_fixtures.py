#!/usr/bin/env python3
"""生成 single-cell-rnaseq 流程的自检用合成数据（10x HDF5 + marker 参考表）。

合成 4 个细胞群，每群有独特的 marker 基因高表达，用于验证：
质控过滤 -> 归一化 -> HVG -> PCA -> 聚类 -> marker 识别 是否都能跑通。
"""

import os

import h5py
import numpy as np
import pandas as pd

rng = np.random.default_rng(42)
BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures_scrna")
os.makedirs(BASE, exist_ok=True)

n_cells = 600
n_genes = 800

# 基因名：前 60 个是各群 marker
gene_names = []
for g in range(4):
    for m in range(15):
        gene_names.append(f"MARKER{g}_{m}")
for g in range(n_genes - 60):
    gene_names.append(f"GENE{g:04d}")

clusters = np.repeat([0, 1, 2, 3], n_cells // 4)
rng.shuffle(clusters)

# 每群的 marker 基因索引
marker_idx = {g: [g * 15 + m for m in range(15)] for g in range(4)}

# 基础表达（gamma 分布）
base = rng.gamma(shape=2.0, scale=8.0, size=n_genes)

counts = np.zeros((n_genes, n_cells), dtype=np.int32)
for j in range(n_cells):
    c = clusters[j]
    mu = base.copy()
    for idx in marker_idx[c]:
        mu[idx] *= 25.0
    # 加一些噪声与线粒体污染
    libsize = rng.uniform(0.6, 1.4)
    mu *= libsize
    counts[:, j] = rng.poisson(mu)

barcodes = [f"AAACCTGAGAAACCAT{i:05d}" for i in range(n_cells)]

# 写 10x HDF5（CSR，列为细胞，行为基因 —— 与 10x 官方格式一致）
h5_path = os.path.join(BASE, "test_matrix.h5")
import scipy.sparse as sp

csr = sp.csr_matrix(counts)  # (n_genes, n_cells)
with h5py.File(h5_path, "w") as f:
    g = f.create_group("matrix")
    g.create_dataset("barcodes", data=np.array([b.encode() for b in barcodes]))
    g.create_dataset("data", data=csr.data.astype(np.int32))
    g.create_dataset("indices", data=csr.indices.astype(np.int64))
    g.create_dataset("indptr", data=csr.indptr.astype(np.int64))
    g.create_dataset("shape", data=np.array([n_genes, n_cells], dtype=np.int32))
    fg = g.create_group("features")
    fg.create_dataset("name", data=np.array([n.encode() for n in gene_names]))
    fg.create_dataset("id", data=np.array([n.encode() for n in gene_names]))

# marker 参考表
rows = []
for g in range(4):
    rows.append({
        "class_label": f"Type{g}",
        "subclass_label": f"Sub{g}",
        "subclass.markers.combo": ",".join(gene_names[i] for i in marker_idx[g]),
    })
pd.DataFrame(rows).to_csv(os.path.join(BASE, "markers.csv"), index=False)

# 真值 cluster（供自检脚本比对）
pd.DataFrame({"cell_id": barcodes, "true_cluster": clusters}).to_csv(
    os.path.join(BASE, "truth_clusters.csv"), index=False
)

print("fixtures written to", BASE)
print("cells:", n_cells, "genes:", n_genes, "clusters:", 4)
