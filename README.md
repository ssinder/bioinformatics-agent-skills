<div align="center">

<img src="docs/assets/logo.svg" alt="BioinfoAgent Skills" width="348">

**Turn any AI coding assistant into a computational biologist's research partner**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python](https://img.shields.io/badge/Python-3.9%2B-blue.svg)](https://www.python.org/)
[![Skills](https://img.shields.io/badge/Skills-17-brightgreen.svg)](#available-skills)
[![Agent Skills](https://img.shields.io/badge/Standard-Agent_Skills-blueviolet.svg)](https://agentskills.io/)
[![CI](https://github.com/ssinder/bioinformatics-agent-skills/actions/workflows/skill-spec-validation.yml/badge.svg)](https://github.com/ssinder/bioinformatics-agent-skills/actions/workflows/skill-spec-validation.yml)
[![Tests](https://github.com/ssinder/bioinformatics-agent-skills/actions/workflows/tests.yml/badge.svg)](https://github.com/ssinder/bioinformatics-agent-skills/actions/workflows/tests.yml)
[![Security](https://github.com/ssinder/bioinformatics-agent-skills/actions/workflows/security-scan.yml/badge.svg)](https://github.com/ssinder/bioinformatics-agent-skills/actions/workflows/security-scan.yml)
[![Pages](https://img.shields.io/badge/Docs-GitHub_Pages-blue.svg)](https://ssinder.github.io/bioinformatics-agent-skills/)

📖 **[Full documentation](https://ssinder.github.io/bioinformatics-agent-skills/)** ——
quick start, skill catalogue, methodology references, and troubleshooting

</div>

---

## What is this

A curated collection of **Agent Skills for computational biology**. Every skill is a standardized
`SKILL.md` package that gives Claude Code, Cursor, Codex, and any other host supporting the
open [Agent Skills standard](https://agentskills.io/) **domain expertise, reproducible workflows, and
validation checks** — instead of another generalist that merely knows how to call an API.

Coverage spans transcriptomics, single-cell analysis, protein structure, molecular dynamics, protein
docking, inter-chain contact analysis with PyMOL visualization, antimicrobial peptide design, ML
peptide generation, literature mining, NCBI BioProject metadata, and publication-grade figures.

> Refactored from a private archive of **217 real research scripts** into a discoverable,
> triggerable, reproducible skill library that follows an open standard. The per-script mapping —
> what was ported, what was rewritten, and why — is documented in the
> **[migration map](docs/migration-map.md)**.

---

## Why

Most AI assistants can write `scanpy` code. They do not know:

- how much mitochondrial gene percentage to drop from your single-cell data before the filter starts
  eating real cells
- why DESeq2's `padj` suddenly came back all `NA` — almost always a mis-specified `design`
- which flexible loop the first spike in an RMSF curve corresponds to, and what experiment to run on it
- which of 3,000 retrieved papers are genuinely *target discovery* rather than *functional validation*
- where to truncate the axes of a reviewer-requested ROC curve so it does not invite criticism

**Those are the things this project carries.** Each skill ships a `references/methodology.md`
recording the scientific reasoning behind parameter choices, a `references/troubleshooting.md`
recording root causes for common failures, and a Validation section in `SKILL.md` describing how to
check whether your results are right.

---

## Quick start

### Option 1 — install as Claude Code skills

```bash
git clone https://github.com/ssinder/bioinformatics-agent-skills.git
cp -r bioinformatics-agent-skills/skills/* ~/.claude/skills/
```

Restart Claude Code and the skills are discovered automatically. Ask it to *"run a bulk RNA-seq
differential expression analysis"* and it will load the `bulk-rnaseq` skill on its own.

### Option 2 — use as a Python package

```bash
pip install -e .
biopipe --list                      # list every skill
biopipe bulk-rnaseq --help         # usage for one skill
```

### Option 3 — call the scripts directly

```bash
export PYTHONPATH="$PWD/packages"

# peptide physicochemical properties
python skills/peptide-design/scripts/analyze_peptides.py \
    --input peptides.fasta --output results/peptide_stats.csv

# ROC curve
python skills/scientific-figures/scripts/make_figures.py roc \
    --input metrics.csv --output figures/roc.png

# literature search
python skills/literature-search/scripts/search_literature.py \
    --topic "antimicrobial peptide delivery" --source pubmed --max-results 50
```

---

## Available skills

### 🧬 Transcriptomics

| Skill | Description | Key capabilities |
|------|------|----------|
| [`bulk-rnaseq`](skills/bulk-rnaseq/) | RNA-seq quantification and differential expression | FastQC QC · Salmon quantification · DESeq2 differential analysis · GO/KEGG enrichment |
| [`ena-download`](skills/ena-download/) | ENA raw data download | Portal API FASTQ retrieval · resumable transfers · MD5 verification · metadata export |
| [`single-cell-rnaseq`](skills/single-cell-rnaseq/) | Standard single-cell RNA-seq workflow | QC filtering · CP10K normalization · HVG · PCA · Leiden clustering · UMAP · cell annotation |

### 🔬 Structural biology

| Skill | Description | Key capabilities |
|------|------|----------|
| [`protein-structure`](skills/protein-structure/) | PDB structure parsing and secondary structure | Built-in PDB parser · DSSP secondary structure · geometry computation · UniProt metadata |
| [`protein-interactions`](skills/protein-interactions/) | Protein-ligand interactions | PLIP XML parsing · residue-level extraction · interaction networks · PyMOL annotations |
| [`molecular-dynamics`](skills/molecular-dynamics/) | Molecular dynamics simulation and trajectory analysis | RMSD/RMSF/PCA · hydrogen bonds · Ramachandran · trajectory clustering · built-in DCD parser |
| [`alphafold-db`](skills/alphafold-db/) | AlphaFold DB structure download | API link resolution · batch download · pLDDT confidence filtering · metadata export |
| [`protein-docking-hdock`](skills/protein-docking-hdock/) | HDOCK protein docking | HTTP job submission · Docking Score · confidence formula · model download · queue-status distinction |
| [`protein-docking-lightdock`](skills/protein-docking-lightdock/) | LightDock local protein docking | GSO swarm optimization · DFIRE scoring parsing · best-complex assembly · resumable batch runs |
| [`pymol-structure-visualization`](skills/pymol-structure-visualization/) | Inter-chain residue contacts and PyMOL visualization | Grid spatial partitioning · hydrogen bond/salt bridge/hydrophobic/π-π/cation-π classification · hotspot residues · Mermaid network diagrams · PyMOL script generation |

### 🧪 Protein engineering

| Skill | Description | Key capabilities |
|------|------|----------|
| [`peptide-design`](skills/peptide-design/) | Antimicrobial peptide design and screening | Physicochemical properties · hydrophobicity/pI · stability · MMseqs2 clustering for redundancy removal |
| [`ml-peptide-generation`](skills/ml-peptide-generation/) | ML peptide generation and screening | PepMLM masked language model · illegal residue filtering · pseudo-perplexity ranking · reproducible seeds |

### 🧬 Cheminformatics

| Skill | Description | Key capabilities |
|------|------|----------|
| [`bio-database-fetch`](skills/bio-database-fetch/) | PubChem compound and taxonomy lookup | PUG-REST batch properties · CAS number extraction · full taxonomy lineage · no browser dependency |

### 📚 Literature mining

| Skill | Description | Key capabilities |
|------|------|----------|
| [`literature-search`](skills/literature-search/) | Multi-source scholarly literature search | PubMed · OpenAlex · Semantic Scholar · DOI deduplication · RIS/BibTeX export |
| [`bioproject-search`](skills/bioproject-search/) | NCBI BioProject metadata retrieval | E-utilities API · title/description/species · keyword filtering · CSV/JSON export |
| [`paper-fulltext-extract`](skills/paper-fulltext-extract/) | Structured full-text extraction from PDFs | GROBID to TEI · section tree/figure/reference parsing · JSON/MD/HTML output · corpus search · lead sentences |

### 📊 Scientific visualization

| Skill | Description | Key capabilities |
|------|------|----------|
| [`scientific-figures`](skills/scientific-figures/) | Publication-grade figure generation | ROC · confusion matrices · heatmaps · volcano plots · SHAP summaries · unified publication styling |

---

## Design principles

**Self-contained** — core algorithms ship as pure Python implementations. No scanpy installed? The
single-cell workflow still runs. No OpenMM? DCD trajectories still parse, and RMSD/PCA still compute.
Dependencies are optional; capability is not.

**Reproducible** — fixed random seeds, recorded parameters, pinned versions, and run metadata written
alongside the results.

**Verifiable** — every skill defines its own quality-control stage, so there is an explicit way to
check whether the output is correct.

**No hardcoded paths** — all paths are parameterized and all output follows a standard directory
convention. Not a single `/Users/xxx/...` is baked into the code.

**Zero setup** — `pip install -e .` is enough. No cloud services, API keys, or databases to configure.

---

## Layout

```
bioinformatics-agent-skills/
├── skills/                    # 17 independently installable skills
│   ├── bulk-rnaseq/
│   │   ├── SKILL.md           # trigger description + workflow + validation
│   │   ├── scripts/           # executable scripts
│   │   └── references/        # methodology + troubleshooting guide
│   └── ...
├── packages/
│   ├── biopipe_core/          # shared layer: config/logging/HTTP retry/CLI/IO
│   └── biopipe/               # unified CLI entry point and skill registry
├── scripts/
│   ├── validate_skills.py     # skill spec validator (CI core)
│   └── init_skill.py          # skill scaffolding generator
├── tests/                     # unit tests and validator self-tests
├── docs/                      # usage guides and skill catalogue
└── .github/workflows/         # 7 CI pipelines
```

---

## Development

```bash
git clone https://github.com/ssinder/bioinformatics-agent-skills.git
cd bioinformatics-agent-skills
pip install -e ".[dev]"

python -m pytest tests -v              # run the test suite
python scripts/validate_skills.py      # validate every skill
python scripts/validate_skills.py --strict   # treat warnings as failures
```

See [CONTRIBUTING.md](CONTRIBUTING.md) before adding a new skill.

### CI pipelines

| Workflow | Purpose |
|--------|------|
| `skill-spec-validation` | Validates that every skill conforms to the Agent Skills spec |
| `script-syntax` | Multi-version Python syntax check + Ruff linting |
| `tests` | Unit tests + `--help` smoke test for every skill |
| `security-scan` | Dependency vulnerability audit + hardcoded credential scan |
| `pr-skill-scan` | Targeted validation and reporting for skills changed in a PR |
| `pages` | Builds and deploys the documentation site |
| `release` | Automated release packaging |

---

## Citation

If you use this project in your research, please cite:

```bibtex
@software{bioinfoagent2026,
  title  = {BioinfoAgent Skills: Agent Skills for Computational Biology},
  author = {sinder},
  year   = {2026},
  url    = {https://github.com/ssinder/bioinformatics-agent-skills}
}
```

See [CITATION.cff](CITATION.cff) for details.

---

## Contributing

Contributions of new skills, bug fixes, and documentation improvements are welcome. See
[CONTRIBUTING.md](CONTRIBUTING.md). Read [docs/SKILL_SPEC.md](docs/SKILL_SPEC.md) before submitting a
new skill.

Participation in this project implies agreement with the [Code of Conduct](CODE_OF_CONDUCT.md).

## License

MIT © sinder. See [LICENSE](LICENSE).

## Disclaimer

This project provides analysis and computation tooling. Its outputs require independent verification
by the researcher and do not constitute medical, clinical, or regulatory advice.

---

<div align="center">

**[Star this project ★](https://github.com/ssinder/bioinformatics-agent-skills)** — help more
researchers discover these reproducible research skills.

</div>
