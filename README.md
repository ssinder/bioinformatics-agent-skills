<div align="center">

<img src="docs/assets/logo.svg" alt="BioinfoAgent Skills" width="348">

**Turn any AI coding assistant into a computational biologist's research partner**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python](https://img.shields.io/badge/Python-3.9%2B-blue.svg)](https://www.python.org/)
[![Skills](https://img.shields.io/badge/Skills-17-brightgreen.svg)](#available-skills)
[![Agent Skills](https://img.shields.io/badge/Standard-Agent_Skills-blueviolet.svg)](https://agentskills.io/)
[![skills.sh](https://skills.sh/b/ssinder/bioinformatics-agent-skills?style=for-the-badge&height=28)](https://skills.sh/ssinder/bioinformatics-agent-skills)
[![PyPI](https://img.shields.io/pypi/v/biopipe-agent-skills.svg)](https://pypi.org/project/biopipe-agent-skills/)
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

### Option 1 — install with the skills CLI (recommended)

```bash
npx skills add ssinder/bioinformatics-agent-skills
```

This works with 80+ agent hosts (Claude Code, Codex, Cursor, OpenCode, …) and
handles per-host paths automatically. It also supports installing a single
skill or a subset:

```bash
# install one skill
npx skills add ssinder/bioinformatics-agent-skills --skill bulk-rnaseq

# install a workflow subset (2-3 skills is the sweet spot -- see below)
npx skills add ssinder/bioinformatics-agent-skills \
  --skill ena-download --skill bulk-rnaseq --skill scientific-figures
```

> **Install selectively.** Skills carry a description that the agent reads at
> startup to decide whether to load one. Install too many at once and the
> host has to drop some descriptions to fit its context budget — a skill
> without a description can never be matched, so it never fires. Two to three
> skills per workflow is the effective range; the
> [skill guides](docs/skill-guides/) list recommended combinations for common
> research tasks.

### Option 2 — install as Claude Code skills manually

```bash
git clone https://github.com/ssinder/bioinformatics-agent-skills.git
cp -r bioinformatics-agent-skills/skills/* ~/.claude/skills/
```

Restart Claude Code and the skills are discovered automatically. Ask it to *"run a bulk RNA-seq
differential expression analysis"* and it will load the `bulk-rnaseq` skill on its own.

Other hosts scan different directories — `~/.claude/skills/` is Claude Code's
personal path only. See the [supporting guide](SUPPORT.md) for the per-host paths.

### Option 3 — use as a Python package

```bash
pip install biopipe-agent-skills
bioagent --list                      # list every skill
bioagent bulk-rnaseq --help         # usage for one skill
bioagent scientific-figures roc --input metrics.csv --output figures/roc
```

The wheel ships the skill packages themselves, so the installed command runs the
real workflows — no checkout required.

> The console script is `bioagent`. The old `biopipe` name still works as an
> alias, because `biopipe` is already taken on PyPI by an unrelated package.

### Option 4 — call the scripts directly

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

**Zero setup for most skills** — `pip install -e .` is enough. No cloud services, no accounts, no
API keys. Two exceptions are documented rather than hidden:

| Skill | What it needs | Why |
|---|---|---|
| `bio-database-fetch`, `bioproject-search`, `literature-search` | Nothing — the NCBI identity (`tool` + developer email) is registered and shipped with the skills | NCBI's usage policy requires a *developer's* email, not the end user's. An optional `NCBI_API_KEY` only raises the rate limit from 3 to 10 req/s. |
| `paper-fulltext-extract` | A local [GROBID](https://grobid.readthedocs.io/) service for PDF → TEI conversion | GROBID is a Java service; the skill talks to it over HTTP. See the skill's Quick Start. |

Every other skill runs with the standard library alone, or degrades gracefully with an explicit
install hint when an optional scientific stack is missing.

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
│   ├── bioagent_core/          # shared layer: config/logging/HTTP retry/CLI/IO
│   └── bioagent/               # unified CLI entry point and skill registry
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
| `backend-agreement` | Cross-checks the built-in Python statistics against Bioconductor DESeq2 |
| `scrna-backend-agreement` | Cross-checks the built-in single-cell pipeline against scanpy |
| `security-scan` | Dependency vulnerability audit + hardcoded credential scan |
| `pr-skill-scan` | Targeted validation and reporting for skills changed in a PR |
| `pages` | Builds and deploys the documentation site |
| `release` | Automated release packaging |

---

## Verification

The correctness claims in this repository are measured, not asserted.

| What | How | Where |
|---|---|---|
| Built-in Python statistics vs Bioconductor DESeq2 | Paired run on the same fixture; Spearman ρ on log2FC ranking, Jaccard on significant genes, direction agreement | [`tests/test_backend_agreement.py`](tests/test_backend_agreement.py) |
| Built-in single-cell pipeline vs scanpy | Paired run plus comparison against 600 hand-labelled truth clusters (ARI) | [`tests/test_scrna_backend_agreement.py`](tests/test_scrna_backend_agreement.py) |
| Agent routing | 87 evaluation cases across all 17 skills, each with at least one near-miss negative case | [`skills/*/evals/evals.json`](skills/) |
| Spec conformance | Agent Skills spec validation plus cross-harness load checks | [`scripts/validate_skills.py`](scripts/validate_skills.py) |

Full numbers, thresholds and how to reproduce them are in
**[BENCHMARK.md](BENCHMARK.md)**.

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
