---
name: bio-database-fetch
description: "Use when looking up compound properties by PubChem CID, retrieving CAS registry numbers, resolving NCBI Taxonomy IDs to species names, or collecting chemical and taxonomic reference data in bulk. Fetches molecular formula, weight, SMILES, InChI, InChIKey, IUPAC name, XLogP, TPSA and CAS numbers from PubChem PUG-REST, plus full taxonomy lineage from NCBI Taxonomy; batched with CSV/JSON export. 用于化合物性质、CAS 号与物种谱系查询。"
license: MIT
compatibility: Python 3.9+; needs requests for network access; standard library only otherwise
metadata:
  version: "1.0"
  last-reviewed: "2026-10-26"
  skill-author: sinder
  edam-topic: "topic_0637"
  edam-operation: "operation_2421"
  edam-input-data: "data_1086, data_1179"
  edam-output-data: "data_2087, data_3028"
  edam-input-format: "format_3752"
  edam-output-format: "format_3752, format_3464"

# BioDatabase Reference Retrieval


---

## Overview

Retrieves reference data from two official public APIs:

| Database | API | Script |
|---|---|---|
| PubChem compounds | PUG-REST | `scripts/fetch_compounds.py` |
| NCBI Taxonomy | E-utilities `efetch` | `scripts/fetch_taxonomy.py` |

Both are plain HTTP/JSON-or-XML endpoints — **no browser automation is required**, unlike the pipeline this skill was extracted from, where scraping breaks on template changes and cannot be reproduced.

Six PubChem/NCBI traps are handled here (the `CanonicalSMILES` alias, HTTP 200 for a missing CID, `xrefs/RegistryID` not a CAS list, `LineageEx` holding only ancestors, ranks beyond the eight, `root.iter("Taxon")` over-collecting) — all derived in `references/methodology.md`.

## When to Use This Skill

- Getting molecular properties of compounds you already know by CID
- Looking up CAS registry numbers for a set of CIDs
- Converting NCBI Taxonomy IDs into species names and lineages
- Annotating a hits table with compound descriptors or taxonomic classification
- Verifying that a CID exists in PubChem

## Quick Start

```bash
export PYTHONPATH="$PWD/packages"
# No email setup needed — a registered identity ships with the project

# Batch fetch properties + synonyms + CAS
python scripts/fetch_compounds.py --cids cids.txt --output-dir results/

# Direct CIDs, no synonyms
python scripts/fetch_compounds.py --cids 2244 3672 --output-dir results/ --no-synonyms

# Batch fetch lineages
python scripts/fetch_taxonomy.py --tax-ids 9606 10090 --output-dir results/
```

## Workflow

### Compounds (`fetch_compounds.py`)

1. **Collect CIDs** — `--cids` takes CID strings or text files (one per line; `#` comments and header rows skipped); duplicates removed, order kept.
2. **Fetch properties** — one PUG-REST request per 100 CIDs; default set: formula, weight, exact mass, XLogP, TPSA, charge, both SMILES variants, InChI, InChIKey, IUPAC name, H-bond counts, rotatable bonds, heavy atoms.
3. **Fetch synonyms** — skipped with `--no-synonyms`; keeps the first `--max-synonyms` (default 20), reports true `synonym_count`.
4. **Extract CAS numbers** — filtered from synonyms by CAS format, topped up from `RegistryID`; skipped with `--no-cas`.
5. **Export** — `compounds.csv`, `compounds.json`.

### Taxonomy (`fetch_taxonomy.py`)

1. **Collect IDs** — same as CIDs.
2. **Fetch** — `efetch.fcgi?db=taxonomy&id=...&retmode=xml`; taxonomy IDs are primary keys, so no `esearch` step (unlike BioProject).
3. **Parse** — scientific name, common name, rank, lineage.
4. **Export** — `taxonomy.csv` (one row per taxon, `domain`…`species` as columns plus a `full_lineage` string) and `taxonomy.json`.

## Parameters

### `fetch_compounds.py`

| Parameter | Default | Description |
|---|---|---|
| `--cids` | required | CID list or path(s) to text files of CIDs |
| `--output-dir` | `compound_results` | Output dir for the CSV/JSON results |
| `--no-synonyms` | off | Skip synonyms (may be hundreds per compound) |
| `--no-cas` | off | Skip CAS extraction |
| `--max-synonyms` | `20` | Maximum synonyms kept per record |
| `--timeout` | `60` | Per-request timeout in seconds |

### `fetch_taxonomy.py`

| Parameter | Default | Description |
|---|---|---|
| `--tax-ids` | required | Taxonomy ID list or path(s) to text files of IDs |
| `--output-dir` | `taxonomy_results` | Output dir for the CSV/JSON results |
| `--email` | registered project identity | NCBI-required contact email; a registered identity ships with the project, so overriding is optional |
| `--api-key` | `$NCBI_API_KEY` | NCBI API key; raises the rate limit from 3 to 10 req/s |
| `--timeout` | `60` | Per-request timeout in seconds |

## Outputs

| File | Content |
|---|---|
| `compounds.csv` | Per CID: found, IUPAC name, formula, molecular weight, exact mass, XLogP, TPSA, charge, SMILES, ConnectivitySMILES, InChI, InChIKey, H-bond counts, rotatable bonds, heavy atoms, CAS, synonyms |

| `compounds.json` | Same content as nested JSON |
| `taxonomy.csv` | Per taxon: tax_id, scientific name, common name, rank, `domain`…`species`, lineage depth, full lineage string |
| `taxonomy.json` | Full nested records with every lineage entry's rank and taxid |

## Validation

- **Check `found` first.** `found=False` means PubChem has no such CID — not transient; the ID is wrong or the compound was removed.
- **Several CAS numbers may be legitimate.** Salts, hydrates and stereoisomers each get one; `CAS` is semicolon-separated — match the primary only when joining a single-registry table.
- **CAS candidates are format-validated, not registry-validated.** The pattern check drops non-CAS registry IDs but cannot detect a re-assigned number; cross-check the InChIKey.
- **`superkingdom` is empty in current data.** NCBI retired it; "cellular organisms" is now `cellular root`, the top-level split is `domain` (Eukaryota / Bacteria / Archaea) — use `domain`.
- **Not every taxon has a species rank.** Strain-level entries report `no rank` (e.g. *Escherichia coli* str. K-12 substr. MG1655); read `rank`, don't assume `species`.
- **Species names are not stable identifiers.** *E. coli*'s phylum was renamed `Proteobacteria` → `Pseudomonadota`; `tax_id` is stable — join on it.
- **A missing common name is normal.** `GenbankCommonName` exists only for well-studied taxa.
- **Taxonomy IDs are stable, names are not.** Never key a cross-study join on a scientific name.

## References

- `references/methodology.md` — API protocols, field semantics, the six PubChem/NCBI traps, CAS extraction, external docs
- `references/troubleshooting.md` — common failures and diagnosis
