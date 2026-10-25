---
name: peptide-design
description: This skill should be used when analyzing or designing antimicrobial peptides — computing physicochemical properties such as molecular weight, isoelectric point, Kyte-Doolittle hydrophobicity and instability index, classifying peptides as antimicrobial or non-antimicrobial from GO terms and sequence features, deduplicating peptide libraries by clustering, and exporting FASTA or CSV. 用于多肽理化性质计算、抗菌肽判别、序列聚类去冗余与导出的场景。
license: MIT
compatibility: Python 3.9+；理化性质计算为自包含纯 Python 实现（无需 Biopython/numpy）；MMseqs2 为可选外部依赖
metadata:
  version: "1.0"
  last-reviewed: "2026-10-25"
  skill-author: sinder
---

## Overview

从 FASTA 或 CSV 读取多肽序列，用**自包含的公式实现**计算全部理化性质（分子量、等电点、疏水性、净电荷、不稳定指数、二级结构倾向、溶解度、Boman 值），再结合 UniProt GO 术语与关键词做抗菌肽判别打分，最后用聚类去冗余提取代表序列并导出 FASTA / CSV。

理化性质不依赖 Biopython：分子量按氨基酸残基质量表累加加水的质量；等电点用 Henderson-Hasselbalch 在 pH 0–14 的电荷曲线上二分求零点；疏水性用 Kyte-Doolittle 标准值表；不稳定指数用 Guruprasad 的 DIHEX 二面角公式。

## When to Use This Skill

- 需要批量计算多肽的理化性质（分子量、pI、疏水性、净电荷、稳定性）作为筛选依据
- 需要判断一批候选肽中哪些是抗菌肽、哪些是非抗菌肽（用于构建负样本集）
- 需要对肽库做去冗余（MMseqs2 不可用时有内置聚类回退）
- 需要为下游模型准备带特征列的训练数据
- 需要把筛选结果导出为 FASTA 或带理化性质的 CSV

## Quick Start

```bash
export PYTHONPATH="$PWD/packages"

# 1) 计算 FASTA 中所有序列的理化性质
python skills/peptide-design/scripts/analyze_peptides.py \
    --input peptides.fasta --output results/peptides_props.csv \
    --report results/summary.txt

# 2) 判别抗菌肽并导出筛选结果 + FASTA
python skills/peptide-design/scripts/analyze_peptides.py \
    --input peptides.csv --sequence-column sequence \
    --classify --antimicrobial-threshold 5 \
    --label antimicrobial \
    --export-fasta results/amp.fasta --output results/amp.csv

# 3) 提取高质量样本（需 UniProt protein_existence 注释）
python skills/peptide-design/scripts/analyze_peptides.py \
    --input uniprot_peptides.csv --classify --min-quality 60 \
    --require-complete-metadata --output results/high_quality.csv

# 4) 聚类去冗余（优先 MMseqs2，缺失时自动回退内置算法）
python skills/peptide-design/scripts/cluster_peptides.py \
    --input peptides.fasta --method auto --min-identity 0.3 \
    --output results/representatives.fasta --output-csv results/cluster_members.csv
```

## Workflow

1. **读取序列**：FASTA（多行序列、`;` 注释、CRLF 均可）或 CSV（自动嗅探列名，缺 `sequence` 列时依次回退到 `seq`/`peptide`/`aa_seq`）
2. **理化性质计算**（`peptide_properties.py`）
   - 分子量 = Σ残基平均质量 + 18.0153
   - 等电点：pH 0–14 步长 0.1 扫描电荷曲线 → 找符号翻转区间 → 二分 60 次求精确零点
   - 净电荷：Henderson-Hasselbalch 分别计算 N 端/Lys/Arg/His 的质子化与 C 端/Asp/Glu/Cys/Tyr 的去质子化比例
   - 疏水性：Kyte-Doolittle 平均值 + Hobbs 疏水矩（单位向量合成）
   - 不稳定指数：Guruprasad DIHEX 权重表，`II = 100/L × Σ DIHEX(di-1, di)`
   - 构象倾向：Chou-Fasman 螺旋/折叠/拐角参数 + 螺旋终止残基修正
   - 溶解度：Zhao (2009) 残基级 log-s 贡献 + 长度项 + 螺旋项
   - Boman 值、芳香性、脂肪族指数、消光系数、估算 kD
3. **判别打分**：注释层（GO 术语 5 组 + 关键词 5 组加权）× 0.6 + 理化层（8 条规则）× 0.4
4. **质量评估**：UniProt 蛋白质存在水平 + 序列完整性 + 注释完整度 + 理化合理性
5. **筛选**：长度区间 / 未知残基 / 元数据完整度 / 质量分 / 判别类别
6. **聚类去冗余**（`cluster_peptides.py`）：MMseqs2 `easy-cluster` 优先，回退到内置 k-mer Jaccard + LCS identity + average-linkage 聚类
7. **导出**：CSV（UTF-8 BOM）/ FASTA / JSON / 文本报告

## Parameters

| 参数 | 说明 | 默认值 |
| --- | --- | --- |
| `--input` | 输入 FASTA 或 CSV | 必填 |
| `--format` | 输入格式（`auto` 按扩展名与内容嗅探） | `auto` |
| `--sequence-column` | CSV 序列列名 | `sequence` |
| `--id-column` / `--name-column` | CSV 的 ID 列 / 名称列 | `id` / `protein_name` |
| `--go-column` / `--keyword-column` | CSV 的 GO 术语列 / 关键词列 | `go_terms` / `keywords` |
| `--classify` | 启用抗菌肽判别 | 关闭 |
| `--antimicrobial-threshold` | 抗菌肽判定阈值 | 5.0 |
| `--non-antimicrobial-threshold` | 非抗菌肽判定阈值 | 2.0 |
| `--annotation-weight` | 注释层在综合分中的权重 | 0.6 |
| `--label` | 保留哪一类：`antimicrobial` / `non_antimicrobial` / `uncertain` / `all` | `antimicrobial` |
| `--min-length` / `--max-length` | 长度过滤区间 | 0（不限） |
| `--drop-unknown` | 剔除含 X 等未知残基的序列 | 关闭 |
| `--require-complete-metadata` | 只保留同时有 GO 与关键词的记录 | 关闭 |
| `--min-quality` | 最低质量分 | 0 |
| `--export-fasta` / `--report` / `--json-output` | 额外导出 FASTA / 报告 / JSON | 无 |
| `--method`（聚类） | `auto` / `mmseqs` / `internal` | `auto` |
| `--min-identity`（聚类） | 簇内最小序列 identity | 0.3 |
| `--representative`（聚类） | `medoid` / `longest` / `first` | `medoid` |

## Outputs

| 文件 | 说明 |
| --- | --- |
| `peptide_analysis.csv` | 主结果，含 24 个理化性质列 + `label` / `combined_score` / `quality_score` / `qc_flags` |
| `selected.fasta` | 筛选后序列（`--export-fasta`） |
| `summary.txt` | 统计报告：长度/MW/pI/疏水性分布、稳定性占比、判别分布、GO 覆盖度、物种分布 |
| `cluster_members.csv` | 聚类明细：`cluster_id` / `representative` / `member_id` / `is_representative` / `cluster_size` |
| `representatives.fasta` | 代表序列（每簇一条） |
| `cluster_report.txt` | 聚类统计：簇数、冗余率、最大/中位簇大小、单序列簇数 |

## Validation

- **分子量合理性**：20 aa 肽的 MW 应在 2000–2500 Da 区间；含 C 残基时每形成一对二硫键减少约 2 Da
- **等电点合理性**：酸性肽（Asp/Glu 多）pI < 7，碱性肽（Lys/Arg 多）pI > 7；等电点应落在 pH 1–14 之间
- **疏水性**：Kyte-Doolittle 平均值 > 1.8 强疏水，< 0 强亲水；典型两性抗菌肽在 0.5–2.5 且疏水矩 > 0.4
- **稳定性**：不稳定指数 < 40 为预测稳定，> 40 提示可能快速降解（此时 `qc_flags` 会自动标记）
- **未知残基**：含 X 的序列所有定量指标可靠性下降，`qc_flags` 与 `has_unknown` 会标出
- **判别结果**：抽查各 label 的样本，确认 `hit_keywords` 与 `classification_reasons` 合理；纯 FASTA 输入无注释时只能靠理化层，此时 `uncertain` 占比升高属正常
- **聚类合理性**：冗余率应与序列多样性相符——短肽库通常冗余率 > 50%；若全部聚成单序列簇，说明 `--min-identity` 过严
- **代表序列**：medoid 代表序列应与簇内其他成员有肉眼可见的相似性；可用 `--representative longest` 交叉验证

## References

- `references/methodology.md` —— 各理化性质的公式推导、参数选择与判别规则设计原理
- `references/troubleshooting.md` —— 常见报错、数值异常与排错步骤
