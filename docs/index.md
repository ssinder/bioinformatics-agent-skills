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
git clone https://github.com/ssinder/bioinformatics-agent-skills.git
cp -r bioinformatics-agent-skills/skills/* ~/.claude/skills/
```

Restart your host and the skills are discovered automatically. See
[Getting Started]({{ '/getting-started.html' | relative_url }}) for the full walkthrough.

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
| [`alphafold-db`](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills/alphafold-db/SKILL.md) | Downloads predicted protein structures from the EBI AlphaFold Database API and exports per-protein prediction metadata. Resolves the current PDB/CIF download URL from the API instead of hard-coding a model version, supports batch accession lists, resumable downloads, pLDDT-based quality filtering before spending bandwidth, and validates response content types so HTML error pages are never mistaken for structure files. This skill should be used when fetching AlphaFold models, retrieving predicted structures by UniProt accession, or collecting AlphaFold confidence metrics (pLDDT) is requested. |
| [`bio-database-fetch`](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills/bio-database-fetch/SKILL.md) | Fetches compound properties and species taxonomy from the official PubChem PUG-REST and NCBI Taxonomy APIs. Retrieves molecular formula, weight, SMILES, InChI, InChIKey, IUPAC name, XLogP, TPSA and CAS numbers for a batch of PubChem CIDs, and full lineage classification (domain through species) for a batch of NCBI Taxonomy IDs. Supports batched queries (100 IDs per request), CSV/JSON export, and requires a contact email as the NCBI usage policy mandates. This skill should be used when looking up properties of known compounds by CID, retrieving CAS registry numbers, resolving taxonomy IDs to species names, or collecting chemical and taxonomic reference data in bulk is requested. |
| [`bioproject-search`](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills/bioproject-search/SKILL.md) | Retrieves NCBI BioProject metadata — title, description, organism, submitter and release dates — for a batch of BioProject accessions, using the official NCBI E-utilities API rather than scraping the BioProject web pages. Supports keyword filtering over titles and descriptions, CSV/JSON export, and requires a contact email as the NCBI usage policy mandates. This skill should be used when looking up what a BioProject or PRJNA/PRD accession contains, retrieving sequencing project descriptions, or collecting BioProject metadata in bulk is requested. |
| [`bulk-rnaseq`](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills/bulk-rnaseq/SKILL.md) | This skill should be used when running bulk RNA-seq differential expression analysis: FastQC/trimming, Salmon quantification, counts matrix construction, DESeq2 statistics, GO/KEGG enrichment and publication figures. 用于 bulk 转录组（RNA-seq）差异表达分析全流程。 |
| [`ena-download`](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills/ena-download/SKILL.md) | Downloads RNA-seq and other sequencing FASTQ data from the ENA (European Nucleotide Archive) Portal API and exports run-level metadata. Resolves download URLs and file sizes from the structured API rather than scraping archive web pages, supports batch accession lists, HTTP Range resume for interrupted transfers, and MD5 verification against the checksums ENA publishes. This skill should be used when fetching raw sequencing reads by SRA/ERR/DRR accession, retrieving FASTQ files for an RNA-seq project, or collecting run metadata and read counts is requested. |
| [`literature-search`](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills/literature-search/SKILL.md) | This skill should be used when mining academic literature from PubMed, OpenAlex or Semantic Scholar — searching papers by topic or keywords, deduplicating multi-source results, filtering by relevance or citation count, enriching missing abstracts, and exporting RIS, BibTeX or GB/T 7714 citations. 用于多源学术文献检索、去重、相关性筛选、摘要富集与参考文献格式导出的场景。 |
| [`ml-peptide-generation`](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills/ml-peptide-generation/SKILL.md) | Generates candidate binder peptides for a target protein using the PepMLM-650M masked language model — appends <mask> tokens to the target sequence and samples amino acids at each position, with illegal residues (X, B, Z) excluded from the sampling candidates. Reports pseudo-perplexity for each candidate, supports reproducible generation via an explicit random seed, and provides offline filtering and parameter-sweep analysis that needs no GPU. This skill should be used when designing peptides or protein fragments that should bind a target, generating a ranked candidate pool for a target sequence, or comparing how peptide length and top-k sampling settings affect candidate quality. |
| [`molecular-dynamics`](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills/molecular-dynamics/SKILL.md) | Analyses molecular-dynamics trajectories with a self-contained pure-numpy DCD and PDB trajectory parser plus built-in Kabsch-Sander alignment, SVD-based PCA, hierarchical clustering of RMSD matrices, geometric hydrogen-bond detection, RMSF, radius of gyration, Ramachandran analysis, inter-chain contacts and Shrake-Rupley SASA. This skill should be used when trajectory RMSD or RMSF, PCA, conformational clustering, hydrogen-bond occupancy, Ramachandran plots or solvent-accessible surface analysis is requested. |
| [`paper-fulltext-extract`](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills/paper-fulltext-extract/SKILL.md) | This skill should be used when extracting structured full text from scientific PDF papers — converting PDFs to TEI XML via a local GROBID service, parsing titles, authors, abstracts, keywords, section trees, figures and references into JSON, Markdown or HTML, building a searchable corpus index across many papers, and pulling the first sentences of each Introduction paragraph for rapid paper triage. 用于科研 PDF 全文结构化抽取、TEI 解析、文献语料全文检索与引言精读的场景。 |
| [`peptide-design`](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills/peptide-design/SKILL.md) | This skill should be used when analyzing or designing antimicrobial peptides — computing physicochemical properties such as molecular weight, isoelectric point, Kyte-Doolittle hydrophobicity and instability index, classifying peptides as antimicrobial or non-antimicrobial from GO terms and sequence features, deduplicating peptide libraries by clustering, and exporting FASTA or CSV. 用于多肽理化性质计算、抗菌肽判别、序列聚类去冗余与导出的场景。 |
| [`protein-docking-hdock`](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills/protein-docking-hdock/SKILL.md) | Submits protein-protein docking jobs to the HDOCK server and retrieves the resulting Docking Score, Confidence Score, ligand RMSD and ranked complex models. Uses plain HTTP form submission against the HDOCK server rather than browser automation, supports PDB files, PDB ID:ChainID references and FASTA sequences for both receptor and ligand, records submitted jobs to a CSV so results can be collected later, and can poll for completion. This skill should be used when running protein-protein complex structure prediction on the HDOCK server, checking the docking score of a submitted complex, or downloading top-ranked docked models is requested. |
| [`protein-docking-lightdock`](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills/protein-docking-lightdock/SKILL.md) | Runs LightDock protein-protein docking locally by driving the lightdock3_setup.py and lightdock3.py executables, then extracts the best docked complex. Parses the GSO swarm output files with the correct field semantics (glowworm id is the line number, not the LigID column), selects the latest gso_<step>.out by numeric step rather than dictionary order, and generates conformations from the origin-shifted structures that the setup stage produced so the receptor and ligand stay in a consistent coordinate frame. Supports batch receptor-ligand grids, resumable execution, incremental score logging, and DFIRE-based ranking. This skill should be used when protein-protein docking must be performed on the local machine, when a swarm-level score ranking is needed, or when the best docked complex model must be generated as a PDB file. |
| [`protein-interactions`](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills/protein-interactions/SKILL.md) | Analyzes protein-ligand and protein-protein interactions from PLIP XML reports and AlphaFold model packages using a self-contained standard-library XML parser: extracts residue-level salt bridges, hydrogen bonds, hydrophobic, pi-stacking, water bridges, halogen bonds and metal complexes, ranks hotspot residues, and emits Mermaid network graphs plus PyMOL annotation scripts. This skill should be used when interaction detection, hotspot ranking, PLIP report parsing, interaction network visualisation or AlphaFold confidence analysis is requested. |
| [`protein-structure`](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills/protein-structure/SKILL.md) | Analyzes protein structures from PDB files: parses ATOM/HETATM records with a dependency-free pure-Python parser, computes radius of gyration, gyration tensor and dipole moments, assigns secondary structure via DSSP or a built-in Kabsch-Sander hydrogen-bond implementation, and optionally retrieves UniProt metadata. This skill should be used when secondary structure assignment, Rg, shape analysis, DSSP, PDB parsing or protein structure annotation is requested. |
| [`pymol-structure-visualization`](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills/pymol-structure-visualization/SKILL.md) | This skill should be used when analyzing inter-chain contacts in protein structures — computing residue-level contact geometry from PDB files with a uniform-grid spatial index, classifying each contact into hydrogen bonds, salt bridges, hydrophobic contacts, pi-pi stacking and cation-pi interactions using residue chemistry, identifying hotspot residues, and exporting CSV/JSON results together with Mermaid interaction networks and ready-to-run PyMOL visualization scripts. 用于 PDB 链间残基接触分析、结合界面热点残基鉴定、相互作用网络构建与 PyMOL 可视化脚本生成的场景。 |
| [`scientific-figures`](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills/scientific-figures/SKILL.md) | This skill should be used when generating publication-ready scientific figures — ROC curves, confusion matrix heatmaps, docking binding-energy heatmaps, pLDDT confidence histograms, SHAP summary beeswarms, class distribution charts and other chart types with a unified journal-grade visual style. 用于生成顶刊级出版图表（ROC 曲线、混淆矩阵、热图、pLDDT 分布、SHAP 摘要等）的场景。 |
| [`single-cell-rnaseq`](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/skills/single-cell-rnaseq/SKILL.md) | This skill should be used when running single-cell RNA-seq analysis: reading 10x/AnnData matrices, QC filtering, normalization, HVG selection, PCA, clustering, UMAP/t-SNE and cell type annotation. 用于单细胞转录组（scRNA-seq）质控、降维、聚类与细胞注释全流程。 |
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
- [Skills Overview]({{ '/skills/' | relative_url }}) — every skill and the scenarios it fits
- [Skill Spec]({{ '/SKILL_SPEC.html' | relative_url }}) — how to write a new skill
- [Contributing](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/CONTRIBUTING.md)

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
