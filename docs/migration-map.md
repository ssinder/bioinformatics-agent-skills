---
title: Migration map · From the legacy archive to Agent Skills
description: Complete mapping of 217 scripts from a private legacy archive into 17 Agent Skills
---

# Migration map · From the legacy archive to Agent Skills

> The source repository this project was refactored from is now **private and archived**
> (read-only), so its scripts cannot be browsed publicly. This table records all 217 of
> those scripts, their destinations, and the reasoning behind every decision — for
> auditing and traceability.

## At a glance

| | Legacy archive | This repository |
|---|---|---|
| Shape | 22 flat directories, organized by tool name | 17 standard `SKILL.md` skill packages |
| Scripts | 217 `.py` files | 58 `.py` files (**each skill is a complete, independently triggerable workflow**) |
| Entry point | Browse directories to find a script | `SKILL.md` frontmatter + `biopipe --list` |
| Dependencies | Implicit (install and explore on your own) | Per-skill `requirements.txt`; core has zero third-party dependencies |
| Tests | None | 900 test cases (including reverse-verification canaries) |

**Dropping from 217 scripts to 58 is not a reduction in capability** — what was removed
is duplication and scaffolding. Every skill gained what the legacy archive never had:
a methodology write-up, a troubleshooting guide, a validation checklist, tests, and
fixtures.


---

## Complete mapping

### Converted into skills (18 legacy directories → 17 skills)

| Legacy directory | Scripts | Lines | → New skill | Key changes |
|---|---|---|---|---|
| `rnaseq_python` | 16 | 4026 | `bulk-rnaseq` | DESeq2/edgeR differential analysis; fixed the design-matrix bug that made every `padj` come back `NA` |
| `scrna_python` | 6 | 748 | `single-cell-rnaseq` | Single-cell QC; sensible ranges for the mitochondrial percentage threshold |
| `rnaseq_download` | 1 | 285 | `ena-download` | ENA metadata retrieval and FASTQ download |
| `crawl_bioproject` | 1 | 338 | `bioproject-search` | NCBI BioProject metadata scraping |
| `crawl_database` | 13 | 1974 | `bio-database-fetch` | Generic database scraping |
| `afdb_download` | 2 | 364 | `alphafold-db` | AlphaFold model download; **best-model extraction merged in** |
| `search_pubmed` | 20 | 4519 | `literature-search` | Literature search |
| `search_openalex` | 8 | 2063 | `literature-search` | (same skill, merged) |
| `paper_extract` | 13 | 2753 | `paper-fulltext-extract` | Full-text extraction; filled in 11 failing tests |
| `pdb_analyze` | 14 | 2863 | `protein-structure` | Structure parsing and analysis |
| `plip_interact` | 10 | 1197 | `protein-interactions` | Interaction detection |
| `pymol_show` | 22 | 1486 | `pymol-structure-visualization` | Inter-chain residue contacts; **geometry computation fully separated from rendering** |
| `dock_lightdock` | 5 | 1413 | `protein-docking-lightdock` | LightDock docking |
| `dock_hdock` | 10 | 790 | `protein-docking-hdock` | HADDOCK docking |
| `MD_openmm` | 19 | 6881 | `molecular-dynamics` | OpenMM simulation; equilibration protocol + chain-separation monitoring |
| `peptide_analyze` | 18 | 2595 | `peptide-design` | Peptide property analysis |
| `design_pepmlm` | 4 | 899 | `ml-peptide-generation` | ML peptide generation |
| `data_visualization` | 16 | 3036 | `scientific-figures` | Publication-grade figures |

### Not converted (4 legacy directories)

| Legacy directory | Scripts | Lines | Disposition | Rationale |
|---|---|---|---|---|
| `pptx_create` | 9 | 1748 | **Archived and redirected** | General-purpose office tooling, outside the research scope; also deeply tied to a single presentation topic rather than a general implementation. Redirected to the `workbuddy-pptx` skill |
| `feishu_bot` | 10 | 1575 | **Security cleanup** | Not a research tool. All plaintext credentials replaced with environment variables; README gained a security warning; a `.gitignore` was added |
| `omics_analyze` | 0 | 0 | Empty directory | No content |
| `rnaseq_r` | 0 | 0 | Empty directory | No content (the Rust/R rewrite never landed) |

---

## Three places with no one-to-one mapping

### 1. `search_pubmed` + `search_openalex` → a single `literature-search`

The 28 scripts across these two directories did the same thing against different search
backends. Splitting them into two skills makes an AI assistant hesitate at selection
time — it cannot tell whether to reach for a "PubMed searcher" or a "literature searcher".

`literature-search` now organizes by backend internally, and its trigger keywords cover
both.

### 2. `afdb_download`'s model extraction was merged into `alphafold-db`

In the original directory, "download" and "extract the best model" were separate scripts.
Keeping them apart means the user has to remember the order — for a workflow that has
exactly one sensible order.

### 3. `plip_interact`'s interaction detection became its own `protein-interactions` skill

It also handles structures, like `pdb_analyze`, but the **trigger scenarios are entirely
different**: the former answers "what force acts between these two residues", the latter
answers "what does this structure look like, and what secondary structures does it
contain". Merging them would blur the skill descriptions.

---

## Real defects fixed during the refactor

This was not "copying files over" — it was **reading and then rewriting**. Below are the
substantive defects found and fixed in each skill during the refactor (not stylistic
issues — these produced scientifically wrong conclusions):

| Skill | Defect in the original implementation | Consequence |
|---|---|---|
| `pymol-structure-visualization` | `cmd.select(f"{a} within {cutoff} of {b}")` selected residues that were *near* each other, then took the Cartesian product | Two key atoms 15 Å apart would still be reported as a "contact" |
| `pymol-structure-visualization` | ASP/GLU were listed as hydrogen-bond donors | Carboxylates carry no hydrogen at physiological pH, systematically inflating hydrogen-bond counts |
| `pymol-structure-visualization` | All HIS residues were counted as positively charged | Imidazole is only ~4% protonated, inflating HIS salt bridges by roughly 25× |
| `pymol-structure-visualization` | Ligands (LIG/HEM) were forced through the residue classification table | A benzene ligand was reported as a "hydrophobic residue", a charged ligand as a "negative residue" |
| `bulk-rnaseq` | A mis-specified `design` produced all-`NA` `padj` values | Silently yielded meaningless differential-expression results |
| `molecular-dynamics` | Incomplete equilibration protocol | Entered production simulation without energy minimization |

Every skill has a `references/methodology.md` recording the full methodological
trade-offs, and a `references/troubleshooting.md` recording common symptoms and
diagnostic paths.

---

## Security handling of the legacy repository

⚠️ **A Feishu App Secret was leaked in the legacy repository's history.**

- Location: 10 files under the legacy `feishu_bot/` directory, introduced in commit `a3e456c`
- The legacy repository was **public**, so the credential was **extractable from its history**
- The repository has since been made **private and archived** (read-only), so the historical
  credential is no longer publicly reachable
- **The effective remedy was to reset the App Secret in the Feishu developer console** — done,
  which immediately invalidated the old value
- The code-side cleanup was completed and preserved on a `security/credential-cleanup` branch:
  all 10 files now read from environment variables, create their clients lazily, fail with a
  clear message when configuration is missing, and a `.gitignore` excludes `.env` / `*.secret`

This repository's `scripts/ci_security_scan.py` is a direct product of that incident — it
originally matched on `\b(secret)\b`, and because underscores and letters are **both word
characters**, `APP_SECRET` slipped through the scan entirely. The matcher is now fixed, and
32 regression tests lock the behaviour in — including a **reverse-verification canary**:
revert the regex to the old form and the test fails immediately.

