#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""生成 bulk-rnaseq 流程的自检用合成数据（Salmon 定量 + 元数据 + GO 注释）。"""

import os
import numpy as np
import pandas as pd

rng = np.random.default_rng(7)
BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures_bulk")
QUANT = os.path.join(BASE, "salmon")
os.makedirs(QUANT, exist_ok=True)

samples = []
for grp in ("ctrl", "trt"):
    for rep in (1, 2, 3):
        samples.append((f"{grp}_rep{rep}", grp))

n_tx = 300
tx_ids = [f"TX{i:04d}" for i in range(n_tx)]
gene_of = {t: f"GENE{i // 3:03d}" for i, t in enumerate(tx_ids)}
genes = sorted({g for g in gene_of.values()})

# 每个基因的 baseline 表达与抽样噪声
base_mu = rng.gamma(shape=3.0, scale=60.0, size=len(genes))
gene_index = {g: i for i, g in enumerate(genes)}

# 前 15 个基因在 trt 组上调 4 倍
up_genes = set(genes[:15])

rows_meta = []
for sname, grp in samples:
    sf = rng.uniform(0.85, 1.2)
    mult = np.ones(len(genes))
    if grp == "trt":
        for g in up_genes:
            mult[gene_index[g]] = 4.0
    gene_mu = base_mu * mult * sf
    # 分配到转录本：同一基因的 3 个转录本共享权重
    weights = rng.dirichlet(np.ones(3), size=len(genes))
    per_tx_mu = (gene_mu[:, None] * weights).ravel()
    counts = rng.negative_binomial(n=10.0, p=10.0 / (10.0 + per_tx_mu))
    length = rng.integers(300, 3000, size=n_tx)
    rate = counts / length
    tpm = rate / rate.sum() * 1e6
    df = pd.DataFrame({"Name": tx_ids, "TPM": tpm, "NumReads": counts.astype(float), "Length": length})
    d = os.path.join(QUANT, sname)
    os.makedirs(d, exist_ok=True)
    df.to_csv(os.path.join(d, "quant.sf"), sep="\t", index=False)
    rows_meta.append({"sample": sname, "condition": grp, "replicate": sname.split("_")[-1].replace("rep", "")})

pd.DataFrame(rows_meta).to_csv(os.path.join(BASE, "metadata.csv"), index=False)

pd.DataFrame({"transcript_id": tx_ids, "gene_id": [gene_of[t] for t in tx_ids]}).to_csv(
    os.path.join(BASE, "tx2gene.tsv"), sep="\t", index=False
)

# GO 注释：让前 30 个基因共享一个 GO term
go_rows = []
for i, g in enumerate(genes):
    terms = "ribosome assembly [GO:0005840]" if i < 30 else "oxidative phosphorylation [GO:0006119]"
    go_rows.append({
        "EnsemblPlants": g,
        "Gene Ontology (biological process)": terms,
        "Gene Ontology (cellular component)": f"cytosol [GO:0005829]" if i < 30 else "",
        "Gene Ontology (molecular function)": "",
    })
pd.DataFrame(go_rows).to_csv(os.path.join(BASE, "go_annotation.csv"), index=False)

print("fixtures written to", BASE)
