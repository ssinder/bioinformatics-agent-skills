---
layout: default
title: BioinfoAgent Skills
---

# BioinfoAgent Skills

<img src="{{ '/assets/logo.svg' | relative_url }}" alt="BioinfoAgent Skills" width="348" style="max-width:100%;">

A **standard library of Agent Skills for computational biology** — built on the open
[Agent Skills standard](https://agentskills.io/). It packages RNA-seq, single-cell, structural
biology, molecular dynamics, peptide design, and literature-mining workflows into drop-in
skills that any AI agent can pick up and run.

<div class="hero">
  <p class="tagline">Turn reproducible research workflows into skills your agent can call directly</p>
</div>

> Refactored from a private archive of **217 real research scripts** into a discoverable,
> triggerable, reproducible skill library that follows an open standard. Each skill carries the
> scientific reasoning behind its parameter choices, the root causes of its common failures,
> and an explicit validation step — so an agent does not just run the code, it runs it correctly.

## Quick start

```bash
npx skills add ssinder/bioinformatics-agent-skills
```

Works with 80+ agent hosts (Claude Code, Codex, Cursor, OpenCode, …) and handles
per-host paths automatically. To install only what a given workflow needs, name the
skills explicitly:

```bash
npx skills add ssinder/bioinformatics-agent-skills \
  --skill ena-download --skill bulk-rnaseq --skill scientific-figures
```

> **Install selectively.** The agent reads each skill's description at startup to
> decide whether to load one. Installing too many at once makes the host drop
> some descriptions to fit its context budget — and a skill without a description
> can never be matched, so it never fires. Two to three per workflow is the
> effective range; the [skill guides]({{ '/skill-guides/' | relative_url }}) list
> recommended combinations.

See [Getting Started]({{ '/getting-started.html' | relative_url }}) for the full walkthrough.

## Why Agent Skills

| Capability | Details |
|---|---|
| **Open standard** | `SKILL.md` + `scripts/` + `references/`, following the [agentskills.io](https://agentskills.io/) spec — works with Claude Code, Cursor, Codex, and other hosts |
| **Graceful degradation** | `--help` is a metadata query: it never crashes because numpy / pandas are missing. You are prompted to install a dependency only when a feature actually needs it |
| **Zero heavy dependencies** | The core install needs only `requests` + `pyyaml`; everything else comes in through optional extras |
| **Cross-version compatible** | CI covers Python 3.9 / 3.11 / 3.13 |
| **Tested** | 900 unit tests covering skill CLIs, validators, and scientific algorithms |
| **Quality guaranteed** | 7 CI pipelines: spec validation, syntax & lint, unit tests, dependency & credential security scans, PR skill scans, docs site deployment, and release |

## Skills overview

<!-- SKILLS_TABLE_START -->
| [`alphafold-db`](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills/alphafold-db/SKILL.md) | Use when fetching AlphaFold predicted structures, retrieving models by UniProt accession, or collecting AlphaFold confidence metrics (pLDDT). Downloads PDB/CIF models from the EBI AlphaFold Database API, resolving the current URL from the API instead of hard-coding a model version; supports batch accession lists, resumable downloads, pLDDT filtering, and content-type validation. 用于下载 AlphaFold 结构与 pLDDT。 |
| [`bio-database-fetch`](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills/bio-database-fetch/SKILL.md) | Use when looking up compound properties by PubChem CID, retrieving CAS registry numbers, resolving NCBI Taxonomy IDs to species names, or collecting chemical and taxonomic reference data in bulk. Fetches molecular formula, weight, SMILES, InChI, InChIKey, IUPAC name, XLogP, TPSA and CAS numbers from PubChem PUG-REST, plus full taxonomy lineage from NCBI Taxonomy; batched with CSV/JSON export. 用于化合物性质、CAS 号与物种谱系查询。 |
| [`bioproject-search`](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills/bioproject-search/SKILL.md) | Use when looking up what a BioProject or PRJNA/PRD accession contains, retrieving sequencing project descriptions, or collecting BioProject metadata in bulk. Retrieves title, description, organism, submitter and release dates from the official NCBI E-utilities API rather than scraping web pages; supports keyword filtering with CSV/JSON export. 用于 BioProject 元数据检索与批量导出。 |
| [`bulk-rnaseq`](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills/bulk-rnaseq/SKILL.md) | Use when running bulk RNA-seq differential expression analysis: FastQC/trimming, Salmon quantification, counts matrix construction, DESeq2 statistics, GO/KEGG enrichment and publication figures. 用于 bulk 转录组（RNA-seq）差异表达分析全流程。 |
| [`ena-download`](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills/ena-download/SKILL.md) | Use when fetching raw sequencing reads by SRA/ERR/DRR accession, retrieving FASTQ files for an RNA-seq project, or collecting run metadata and read counts. Downloads FASTQ from the ENA Portal API, resolving URLs and sizes from the structured API rather than scraping archive pages; supports batch accession lists, HTTP Range resume, and MD5 verification. 用于下载原始 reads 与运行元数据。 |
| [`literature-search`](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills/literature-search/SKILL.md) | Use when mining academic literature from PubMed, OpenAlex or Semantic Scholar — searching by topic or keywords, deduplicating multi-source results, filtering by relevance or citation count, enriching missing abstracts, and exporting RIS, BibTeX or GB/T 7714 citations. 用于多源学术文献检索、去重、相关性筛选、摘要富集与参考文献格式导出的场景。 |
| [`ml-peptide-generation`](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills/ml-peptide-generation/SKILL.md) | Use when designing peptides or protein fragments that should bind a target, generating a ranked candidate pool for a target sequence, or comparing how peptide length and top-k sampling affect candidate quality. Generates candidates with the PepMLM-650M masked language model, excluding illegal residues (X, B, Z); reports pseudo-perplexity per candidate, supports reproducible seeds, and needs no GPU. 用于多肽生成与候选排序。 |
| [`molecular-dynamics`](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills/molecular-dynamics/SKILL.md) | Use when trajectory RMSD or RMSF, PCA, conformational clustering, hydrogen-bond occupancy, Ramachandran plots or solvent-accessible surface analysis is requested. Analyses trajectories with a self-contained pure-numpy DCD and PDB parser plus built-in Kabsch-Sander alignment, SVD-PCA, hierarchical clustering, geometric hydrogen-bond detection, radius of gyration, inter-chain contacts and Shrake-Rupley SASA. 用于轨迹 RMSD/RMSF/PCA 与构象分析。 |
| [`paper-fulltext-extract`](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills/paper-fulltext-extract/SKILL.md) | Use when extracting structured full text from scientific PDF papers, parsing section trees and figures, building a searchable corpus across many papers, or pulling lead sentences for rapid paper triage. Converts PDFs to TEI XML via a local GROBID service and parses titles, authors, abstracts, keywords, sections, figures and references into JSON, Markdown or HTML. 用于 PDF 全文抽取与语料检索。 |
| [`peptide-design`](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills/peptide-design/SKILL.md) | Use when analyzing or designing antimicrobial peptides — computing physicochemical properties, classifying peptides as antimicrobial or non-antimicrobial, deduplicating peptide libraries by clustering, and exporting FASTA or CSV. Computes molecular weight, isoelectric point, hydrophobicity and instability index; classifies peptides from GO terms and sequence features; clusters libraries with MMseqs2. 用于多肽性质计算、抗菌肽判别与去冗余。 |
| [`protein-docking-hdock`](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills/protein-docking-hdock/SKILL.md) | Use when running protein-protein complex structure prediction on the HDOCK server, checking the docking score of a submitted complex, or downloading top-ranked docked models. Submits jobs by plain HTTP rather than browser automation, accepting PDB files, PDB ID:ChainID references and FASTA sequences; retrieves Docking Score, Confidence Score, ligand RMSD and ranked models, records submissions to CSV, and polls for completion. 用于 HDOCK 对接与结果获取。 |
| [`protein-docking-lightdock`](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills/protein-docking-lightdock/SKILL.md) | Use when protein-protein docking must be performed on the local machine, when a swarm-level score ranking is needed, or when the best docked complex model must be generated as a PDB file. Runs lightdock3 locally and extracts the best complex with correct GSO field semantics, numeric step selection, and a single consistent coordinate frame; supports batch grids, resumable execution and DFIRE ranking. 用于本地 LightDock 对接。 |
| [`protein-interactions`](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills/protein-interactions/SKILL.md) | Use when interaction detection, hotspot ranking, PLIP report parsing, interaction network visualisation or AlphaFold confidence analysis is requested. Parses PLIP XML reports and AlphaFold packages with a standard-library XML parser; extracts residue-level salt bridges, hydrogen bonds, hydrophobic, pi-stacking, water bridges, halogen bonds and metal complexes, ranks hotspot residues, and emits Mermaid graphs plus PyMOL scripts. 用于相互作用解析与热点鉴定。 |
| [`protein-structure`](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills/protein-structure/SKILL.md) | Use when secondary structure assignment, radius of gyration, shape analysis, DSSP, PDB parsing or protein structure annotation is requested. Parses ATOM/HETATM records with a dependency-free pure-Python parser, computes radius of gyration, gyration tensor and dipole moments, assigns secondary structure via DSSP or a built-in Kabsch-Sander implementation, and optionally retrieves UniProt metadata. 用于 PDB 解析与二级结构指派。 |
| [`pymol-structure-visualization`](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills/pymol-structure-visualization/SKILL.md) | Use when analyzing inter-chain contacts in protein structures, identifying binding-interface hotspot residues, or building interaction networks and PyMOL visualization scripts. Computes residue-level contact geometry with a uniform-grid spatial index, classifies each contact into hydrogen bonds, salt bridges, hydrophobic, pi-pi stacking and cation-pi interactions, and exports CSV/JSON with Mermaid networks and PyMOL scripts. 用于链间接触分析、热点鉴定与 PyMOL 脚本。 |
| [`scientific-figures`](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills/scientific-figures/SKILL.md) | Use when generating publication-ready scientific figures: ROC curves, confusion matrix heatmaps, docking binding-energy heatmaps, pLDDT confidence histograms, SHAP summary beeswarms, or class distribution charts. Renders each chart type with a unified journal-grade visual style. 用于生成顶刊级出版图表（ROC 曲线、混淆矩阵、热图、SHAP 摘要等）。 |
| [`single-cell-rnaseq`](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills/single-cell-rnaseq/SKILL.md) | Use when running single-cell RNA-seq analysis: reading 10x/AnnData matrices, QC filtering, normalization, HVG selection, PCA, clustering, UMAP/t-SNE and cell type annotation. 用于单细胞转录组（scRNA-seq）质控、降维、聚类与细胞注释全流程。 |
<!-- SKILLS_TABLE_END -->

## How this differs from a script dump

A plain script repo makes you read the docs, assemble the right commands, and judge the
parameters yourself. Agent Skills writes that knowledge **into the skill description** — an agent
reads one line of `SKILL.md` and knows which skill to use, how to invoke it, and what to watch
out for.

```yaml
---
name: bulk-rnaseq
description: >
  This skill should be used when running bulk RNA-seq differential
  expression analysis: FastQC/trimming, Salmon quantification, counts
  matrix construction, DESeq2 statistics, GO/KEGG enrichment and
  publication figures. 用于 bulk 转录组（RNA-seq）差异表达分析全流程。
---
```

That description is the routing signal the agent dispatches on.

## Coverage

### 🧬 Transcriptomics

- `bulk-rnaseq` — RNA-seq quantification and differential expression (FastQC, Salmon, DESeq2, GO/KEGG enrichment)
- `ena-download` — ENA raw read download with resumable transfers and MD5 verification
- `single-cell-rnaseq` — the standard scRNA-seq workflow: QC, normalization, clustering, UMAP, annotation

### 🔬 Structural biology

- `protein-structure` — PDB parsing, secondary structure assignment, shape analysis
- `protein-interactions` — PLIP report parsing, residue-level interactions, PyMOL annotations
- `molecular-dynamics` — RMSD/RMSF/PCA, hydrogen bonds, Ramachandran, trajectory clustering
- `alphafold-db` — AlphaFold DB structure retrieval with pLDDT confidence filtering
- `protein-docking-hdock` — HDOCK server docking jobs, scores, and model download
- `protein-docking-lightdock` — local LightDock swarms with DFIRE-based ranking
- `pymol-structure-visualization` — inter-chain contact geometry and ready-to-run PyMOL scripts

### 🧪 Protein engineering

- `peptide-design` — antimicrobial peptide properties, classification, and library deduplication
- `ml-peptide-generation` — PepMLM masked language model generation with reproducible seeds

### 📚 Literature mining

- `literature-search` — multi-source search across PubMed, OpenAlex, and Semantic Scholar
- `bioproject-search` — NCBI BioProject metadata via the official E-utilities API
- `paper-fulltext-extract` — structured full text from scientific PDFs via GROBID/TEI
- `bio-database-fetch` — PubChem compound properties and NCBI Taxonomy lineages

### 📊 Scientific visualization

- `scientific-figures` — publication-grade ROC curves, confusion matrices, heatmaps, and SHAP summaries

## Documentation

- [Getting Started]({{ '/getting-started.html' | relative_url }}) — install, configure, and run your first skill
- [Skill Guides]({{ '/skill-guides/' | relative_url }}) — recommended skill combinations per research workflow
- [Skills Overview]({{ '/skills/' | relative_url }}) — every skill and the scenarios it fits
- [Benchmark]({{ '/BENCHMARK.html' | relative_url }}) — measured correctness: backends agree, routing is tested
- [Skill Spec]({{ '/SKILL_SPEC.html' | relative_url }}) — how to write a new skill
- [Roadmap](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/ROADMAP.md) — what is planned and what is deliberately not
- [Contributing](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/CONTRIBUTING.md)
- [Support](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/SUPPORT.md)

## Citation

If you use this project in your research, please cite:

```bibtex
@misc{bioinfoagentskills2026,
  title        = {BioinfoAgent Skills: An Open Standard Agent Skills Library for Computational Biology},
  author       = {ssinder},
  year         = {2026},
  howpublished = {\url{https://github.com/ssinder/bioinformatics-agent-skills}},
  note         = {Accessed: 2026-10-25}
}
```

## License

MIT — see [LICENSE](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/LICENSE).
