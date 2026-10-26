---
layout: default
title: 技能总览
---

# 技能总览 / Skills Overview

> 本页由 `scripts/build_docs_index.py` 自动生成，请勿手工编辑。

下表列出本仓库提供的全部 Agent Skills。

| 技能 / Skill | 说明 / Description |
|---|---|
| [`alphafold-db`](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills/alphafold-db/SKILL.md) | Downloads predicted protein structures from the EBI AlphaFold Database API and exports per-protein prediction metadata. Resolves the current PDB/CIF download URL from the API instead of hard-coding a model version, supports batch accession lists, resumable downloads, pLDDT-based quality filtering before spending bandwidth, and validates response content types so HTML error pages are never mistaken for structure files. This skill should be used when fetching AlphaFold models, retrieving predicted structures by UniProt accession, or collecting AlphaFold confidence metrics (pLDDT) is requested. |
| [`bulk-rnaseq`](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills/bulk-rnaseq/SKILL.md) | This skill should be used when running bulk RNA-seq differential expression analysis: FastQC/trimming, Salmon quantification, counts matrix construction, DESeq2 statistics, GO/KEGG enrichment and publication figures. 用于 bulk 转录组（RNA-seq）差异表达分析全流程。 |
| [`literature-search`](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills/literature-search/SKILL.md) | This skill should be used when mining academic literature from PubMed, OpenAlex or Semantic Scholar — searching papers by topic or keywords, deduplicating multi-source results, filtering by relevance or citation count, enriching missing abstracts, and exporting RIS, BibTeX or GB/T 7714 citations. 用于多源学术文献检索、去重、相关性筛选、摘要富集与参考文献格式导出的场景。 |
| [`molecular-dynamics`](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills/molecular-dynamics/SKILL.md) | Analyses molecular-dynamics trajectories with a self-contained pure-numpy DCD and PDB trajectory parser plus built-in Kabsch-Sander alignment, SVD-based PCA, hierarchical clustering of RMSD matrices, geometric hydrogen-bond detection, RMSF, radius of gyration, Ramachandran analysis, inter-chain contacts and Shrake-Rupley SASA. This skill should be used when trajectory RMSD or RMSF, PCA, conformational clustering, hydrogen-bond occupancy, Ramachandran plots or solvent-accessible surface analysis is requested. |
| [`peptide-design`](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills/peptide-design/SKILL.md) | This skill should be used when analyzing or designing antimicrobial peptides — computing physicochemical properties such as molecular weight, isoelectric point, Kyte-Doolittle hydrophobicity and instability index, classifying peptides as antimicrobial or non-antimicrobial from GO terms and sequence features, deduplicating peptide libraries by clustering, and exporting FASTA or CSV. 用于多肽理化性质计算、抗菌肽判别、序列聚类去冗余与导出的场景。 |
| [`protein-interactions`](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills/protein-interactions/SKILL.md) | Analyzes protein-ligand and protein-protein interactions from PLIP XML reports and AlphaFold model packages using a self-contained standard-library XML parser: extracts residue-level salt bridges, hydrogen bonds, hydrophobic, pi-stacking, water bridges, halogen bonds and metal complexes, ranks hotspot residues, and emits Mermaid network graphs plus PyMOL annotation scripts. This skill should be used when interaction detection, hotspot ranking, PLIP report parsing, interaction network visualisation or AlphaFold confidence analysis is requested. |
| [`protein-structure`](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills/protein-structure/SKILL.md) | Analyzes protein structures from PDB files: parses ATOM/HETATM records with a dependency-free pure-Python parser, computes radius of gyration, gyration tensor and dipole moments, assigns secondary structure via DSSP or a built-in Kabsch-Sander hydrogen-bond implementation, and optionally retrieves UniProt metadata. This skill should be used when secondary structure assignment, Rg, shape analysis, DSSP, PDB parsing or protein structure annotation is requested. |
| [`scientific-figures`](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills/scientific-figures/SKILL.md) | This skill should be used when generating publication-ready scientific figures — ROC curves, confusion matrix heatmaps, docking binding-energy heatmaps, pLDDT confidence histograms, SHAP summary beeswarms, class distribution charts and other chart types with a unified journal-grade visual style. 用于生成顶刊级出版图表（ROC 曲线、混淆矩阵、热图、pLDDT 分布、SHAP 摘要等）的场景。 |
| [`single-cell-rnaseq`](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills/single-cell-rnaseq/SKILL.md) | This skill should be used when running single-cell RNA-seq analysis: reading 10x/AnnData matrices, QC filtering, normalization, HVG selection, PCA, clustering, UMAP/t-SNE and cell type annotation. 用于单细胞转录组（scRNA-seq）质控、降维、聚类与细胞注释全流程。 |

## 安装方式

将 `skills/` 目录复制到宿主的技能目录即可被发现：

```bash
git clone https://github.com/ssinder/bioinformatics-agent-skills.git
cp -r bioinformatics-agent-skills/skills/* ~/.claude/skills/
```

详见 [快速上手](/getting-started.html)。
