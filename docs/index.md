---
layout: default
title: BioinfoAgent Skills
---

# BioinfoAgent Skills

面向计算生物学的**标准 Agent Skills 工具库** —— 遵循 [Agent Skills 开放标准](https://agentskills.io/)，
把 RNA-seq、单细胞、结构生物学、分子动力学、多肽设计与文献挖掘等分析流程，
封装成任何 AI Agent 都能即插即用的技能。

<div class="hero">
  <p class="tagline">把可复现的科研分析流程，变成 Agent 可以直接调用的技能</p>
</div>

## 快速开始

```bash
git clone https://github.com/ssinder/bioinformatics-agent-skills.git
cp -r bioinformatics-agent-skills/skills/* ~/.claude/skills/
```

重启宿主后技能即被自动发现。完整说明见 [快速上手]({{ '/getting-started.html' | relative_url }})。

## 为什么是 Agent Skills

| 特性 | 说明 |
|---|---|
| **开放标准** | `SKILL.md` + `scripts/` + `references/`，遵循 [agentskills.io](https://agentskills.io/) 规范，兼容 Claude Code、Cursor、Codex 等宿主 |
| **优雅降级** | `--help` 属于元信息查询，不会因缺少 numpy / pandas 等科学计算栈而崩溃；真正用到该功能时才提示安装 |
| **零重型依赖** | 核心安装仅需 `requests` + `pyyaml`，其余按需通过 extras 安装 |
| **跨版本兼容** | CI 覆盖 Python 3.9 / 3.11 / 3.13 |
| **质量保障** | 6 条流水线：规范校验、语法与 lint、单测、依赖与凭据安全扫描、PR 技能扫描、发布 |

## 技能总览

<!-- SKILLS_TABLE_START -->
| [`alphafold-db`](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills/alphafold-db/SKILL.md) | Downloads predicted protein structures from the EBI AlphaFold Database API and exports per-protein prediction metadata. Resolves the current PDB/CIF download URL from the API instead of hard-coding a model version, supports batch accession lists, resumable downloads, pLDDT-based quality filtering before spending bandwidth, and validates response content types so HTML error pages are never mistaken for structure files. This skill should be used when fetching AlphaFold models, retrieving predicted structures by UniProt accession, or collecting AlphaFold confidence metrics (pLDDT) is requested. |
| [`bio-database-fetch`](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills/bio-database-fetch/SKILL.md) | Fetches compound properties and species taxonomy from the official PubChem PUG-REST and NCBI Taxonomy APIs. Retrieves molecular formula, weight, SMILES, InChI, InChIKey, IUPAC name, XLogP, TPSA and CAS numbers for a batch of PubChem CIDs, and full lineage classification (domain through species) for a batch of NCBI Taxonomy IDs. Supports batched queries (100 IDs per request), CSV/JSON export, and requires a contact email as the NCBI usage policy mandates. This skill should be used when looking up properties of known compounds by CID, retrieving CAS registry numbers, resolving taxonomy IDs to species names, or collecting chemical and taxonomic reference data in bulk is requested. |
| [`bioproject-search`](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills/bioproject-search/SKILL.md) | Retrieves NCBI BioProject metadata — title, description, organism, submitter and release dates — for a batch of BioProject accessions, using the official NCBI E-utilities API rather than scraping the BioProject web pages. Supports keyword filtering over titles and descriptions, CSV/JSON export, and requires a contact email as the NCBI usage policy mandates. This skill should be used when looking up what a BioProject or PRJNA/PRD accession contains, retrieving sequencing project descriptions, or collecting BioProject metadata in bulk is requested. |
| [`bulk-rnaseq`](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills/bulk-rnaseq/SKILL.md) | This skill should be used when running bulk RNA-seq differential expression analysis: FastQC/trimming, Salmon quantification, counts matrix construction, DESeq2 statistics, GO/KEGG enrichment and publication figures. 用于 bulk 转录组（RNA-seq）差异表达分析全流程。 |
| [`ena-download`](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills/ena-download/SKILL.md) | Downloads RNA-seq and other sequencing FASTQ data from the ENA (European Nucleotide Archive) Portal API and exports run-level metadata. Resolves download URLs and file sizes from the structured API rather than scraping archive web pages, supports batch accession lists, HTTP Range resume for interrupted transfers, and MD5 verification against the checksums ENA publishes. This skill should be used when fetching raw sequencing reads by SRA/ERR/DRR accession, retrieving FASTQ files for an RNA-seq project, or collecting run metadata and read counts is requested. |
| [`literature-search`](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills/literature-search/SKILL.md) | This skill should be used when mining academic literature from PubMed, OpenAlex or Semantic Scholar — searching papers by topic or keywords, deduplicating multi-source results, filtering by relevance or citation count, enriching missing abstracts, and exporting RIS, BibTeX or GB/T 7714 citations. 用于多源学术文献检索、去重、相关性筛选、摘要富集与参考文献格式导出的场景。 |
| [`ml-peptide-generation`](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills/ml-peptide-generation/SKILL.md) | Generates candidate binder peptides for a target protein using the PepMLM-650M masked language model — appends <mask> tokens to the target sequence and samples amino acids at each position, with illegal residues (X, B, Z) excluded from the sampling candidates. Reports pseudo-perplexity for each candidate, supports reproducible generation via an explicit random seed, and provides offline filtering and parameter-sweep analysis that needs no GPU. This skill should be used when designing peptides or protein fragments that should bind a target, generating a ranked candidate pool for a target sequence, or comparing how peptide length and top-k sampling settings affect candidate quality. |
| [`molecular-dynamics`](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills/molecular-dynamics/SKILL.md) | Analyses molecular-dynamics trajectories with a self-contained pure-numpy DCD and PDB trajectory parser plus built-in Kabsch-Sander alignment, SVD-based PCA, hierarchical clustering of RMSD matrices, geometric hydrogen-bond detection, RMSF, radius of gyration, Ramachandran analysis, inter-chain contacts and Shrake-Rupley SASA. This skill should be used when trajectory RMSD or RMSF, PCA, conformational clustering, hydrogen-bond occupancy, Ramachandran plots or solvent-accessible surface analysis is requested. |
| [`peptide-design`](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills/peptide-design/SKILL.md) | This skill should be used when analyzing or designing antimicrobial peptides — computing physicochemical properties such as molecular weight, isoelectric point, Kyte-Doolittle hydrophobicity and instability index, classifying peptides as antimicrobial or non-antimicrobial from GO terms and sequence features, deduplicating peptide libraries by clustering, and exporting FASTA or CSV. 用于多肽理化性质计算、抗菌肽判别、序列聚类去冗余与导出的场景。 |
| [`protein-docking-hdock`](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills/protein-docking-hdock/SKILL.md) | Submits protein-protein docking jobs to the HDOCK server and retrieves the resulting Docking Score, Confidence Score, ligand RMSD and ranked complex models. Uses plain HTTP form submission against the HDOCK server rather than browser automation, supports PDB files, PDB ID:ChainID references and FASTA sequences for both receptor and ligand, records submitted jobs to a CSV so results can be collected later, and can poll for completion. This skill should be used when running protein-protein complex structure prediction on the HDOCK server, checking the docking score of a submitted complex, or downloading top-ranked docked models is requested. |
| [`protein-interactions`](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills/protein-interactions/SKILL.md) | Analyzes protein-ligand and protein-protein interactions from PLIP XML reports and AlphaFold model packages using a self-contained standard-library XML parser: extracts residue-level salt bridges, hydrogen bonds, hydrophobic, pi-stacking, water bridges, halogen bonds and metal complexes, ranks hotspot residues, and emits Mermaid network graphs plus PyMOL annotation scripts. This skill should be used when interaction detection, hotspot ranking, PLIP report parsing, interaction network visualisation or AlphaFold confidence analysis is requested. |
| [`protein-structure`](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills/protein-structure/SKILL.md) | Analyzes protein structures from PDB files: parses ATOM/HETATM records with a dependency-free pure-Python parser, computes radius of gyration, gyration tensor and dipole moments, assigns secondary structure via DSSP or a built-in Kabsch-Sander hydrogen-bond implementation, and optionally retrieves UniProt metadata. This skill should be used when secondary structure assignment, Rg, shape analysis, DSSP, PDB parsing or protein structure annotation is requested. |
| [`scientific-figures`](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills/scientific-figures/SKILL.md) | This skill should be used when generating publication-ready scientific figures — ROC curves, confusion matrix heatmaps, docking binding-energy heatmaps, pLDDT confidence histograms, SHAP summary beeswarms, class distribution charts and other chart types with a unified journal-grade visual style. 用于生成顶刊级出版图表（ROC 曲线、混淆矩阵、热图、pLDDT 分布、SHAP 摘要等）的场景。 |
| [`single-cell-rnaseq`](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills/single-cell-rnaseq/SKILL.md) | This skill should be used when running single-cell RNA-seq analysis: reading 10x/AnnData matrices, QC filtering, normalization, HVG selection, PCA, clustering, UMAP/t-SNE and cell type annotation. 用于单细胞转录组（scRNA-seq）质控、降维、聚类与细胞注释全流程。 |
<!-- SKILLS_TABLE_END -->

## 与普通脚本仓库的区别

普通脚本仓库要求使用者自己读文档、拼命令、判断参数。
Agent Skills 把这些知识**写进了技能描述里** —— Agent 读一眼 `SKILL.md`
就知道该用哪个技能、怎么调、要注意什么。```yaml
---
name: bulk-rnaseq
description: >
  This skill should be used when running bulk RNA-seq differential
  expression analysis: FastQC/trimming, Salmon quantification, counts
  matrix construction, DESeq2 statistics, GO/KEGG enrichment and
  publication figures. 用于 bulk 转录组（RNA-seq）差异表达分析全流程。
---
```

这段描述就是 Agent 的路由依据。

## 文档

- [快速上手]({{ '/getting-started.html' | relative_url }}) —— 安装、配置、首次使用
- [技能总览]({{ '/skills/' | relative_url }}) —— 全部技能及其适用场景
- [技能规范]({{ '/SKILL_SPEC.html' | relative_url }}) —— 如何编写一个新技能
- [贡献指南](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/CONTRIBUTING.md)

## 引用

如果你在研究中使用了本项目，请引用：

```bibtex
@misc{bioinfoagentskills2026,
  title        = {BioinfoAgent Skills: An Open Standard Agent Skills Library for Computational Biology},
  author       = {ssinder},
  year         = {2026},
  howpublished = {\url{https://github.com/ssinder/bioinformatics-agent-skills}},
  note         = {Accessed: 2026-10-25}
}
```

## 许可证

MIT，见 [LICENSE](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/LICENSE)。
