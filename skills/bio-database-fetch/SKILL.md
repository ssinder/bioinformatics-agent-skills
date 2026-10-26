---
name: bio-database-fetch
description: >-
  Fetches compound properties and species taxonomy from the official PubChem
  PUG-REST and NCBI Taxonomy APIs. Retrieves molecular formula, weight, SMILES,
  InChI, InChIKey, IUPAC name, XLogP, TPSA and CAS numbers for a batch of PubChem
  CIDs, and full lineage classification (domain through species) for a batch of
  NCBI Taxonomy IDs. Supports batched queries (100 IDs per request), CSV/JSON
  export, and requires a contact email as the NCBI usage policy mandates. This
  skill should be used when looking up properties of known compounds by CID,
  retrieving CAS registry numbers, resolving taxonomy IDs to species names, or
  collecting chemical and taxonomic reference data in bulk is requested.
license: MIT
compatibility: Python 3.9+; needs requests for network access; standard library only otherwise
metadata:
  version: "1.0"
  last-reviewed: "2026-10-26"
  skill-author: sinder
---

# BioDatabase Reference Retrieval

## Overview

This skill retrieves reference data from two official, well-documented public
APIs:

| Database | API | Script |
|---|---|---|
| PubChem compounds | PUG-REST | `scripts/fetch_compounds.py` |
| NCBI Taxonomy | E-utilities `efetch` | `scripts/fetch_taxonomy.py` |

Both are plain HTTP/JSON-or-XML endpoints. **No browser automation is
required** — that is the main difference from the pipeline this skill was
extracted from, where every database script drove a Selenium browser session.
Scraping a data portal for structured records is strictly worse than asking the
API: it breaks whenever a template changes, runs an order of magnitude slower,
and produces output that cannot be reproduced from a documented request.

### What this skill gets right that a naive implementation does not

**CAS numbers come from the right endpoint.** The `xrefs/RegistryID` endpoint
does *not* return CAS numbers. It returns several hundred registry IDs from
many different databases mixed together — `CHEBI:15365`, `DB00945`,
`LGCAMP0133.00-11`, and the occasional `50-78-2`. A CAS number has a fixed
format (`NNNNN-NN-N`), so the reliable way to get one is to filter the
`synonyms` list by that pattern. For CID 2244 (aspirin) the synonyms endpoint
yields three CAS candidates; `RegistryID` yields only two, and neither is
guaranteed to be the primary registry number.

**Missing CIDs are detectable.** PubChem answers a request for a non-existent
CID with **HTTP 200** and a body of `{"CID": 999999999}` — just the CID, no
properties. Checking the status code alone will therefore report a
successful lookup of a compound that does not exist. Every record here carries
an explicit `found` flag.

**Lineage includes the taxon itself.** NCBI's `LineageEx` element contains
*only ancestors*. The species level is available solely from the top-level
`<Rank>` and `<ScientificName>`. Reading only `LineageEx` leaves the species
name — the field you usually want — empty.

**All ranks are retained.** `LineageEx` contains dozens of rank values beyond
the textbook eight: `domain`, `subphylum`, `superclass`, `clade`, `suborder`,
`infraorder`, `parvorder`, `superfamily`, `subfamily` and more. Selecting only
the eight standard ranks silently drops `domain` (Eukaryota / Bacteria /
Archaea), which is the single most useful top-level discriminator.

## When to Use This Skill

- Getting molecular properties of compounds you already know by CID
- Looking up CAS registry numbers for a set of CIDs
- Converting a list of NCBI Taxonomy IDs into species names and lineages
- Annotating a hits table with compound descriptors or taxonomic classification
- Verifying that a CID actually exists in PubChem

## Quick Start

```bash
export PYTHONPATH="$PWD/packages"
export BIOPIPE_EMAIL="you@your-institution.edu"   # NCBI 要求的联系邮箱

# 化合物：批量取属性 + 同义词 + CAS
python scripts/fetch_compounds.py --cids cids.txt --output-dir results/

# 化合物：直接给 CID，只要基本属性（不要几百个同义词）
python scripts/fetch_compounds.py --cids 2244 3672 --output-dir results/ --no-synonyms

# 分类学：批量取 lineage
python scripts/fetch_taxonomy.py --tax-ids 9606 10090 --output-dir results/
```

## Workflow

### Compounds (`fetch_compounds.py`)

1. **Collect CIDs** — `--cids` accepts CID strings or paths to text files
   (one per line; `#` comments and a header row are skipped). Duplicates are
   removed while preserving order.
2. **Fetch properties** — one PUG-REST request per 100 CIDs. The default
   property set covers molecular formula, weight, exact mass, XLogP, TPSA,
   formal charge, both SMILES variants, InChI, InChIKey, IUPAC name, H-bond
   donors/acceptors, rotatable bonds and heavy atom count.
3. **Fetch synonyms** — skipped with `--no-synonyms`. A well-known compound
   has several hundred synonyms; the output keeps the first `--max-synonyms`
   (default 20) but reports the true `synonym_count`.
4. **Extract CAS numbers** — filtered out of the synonyms by CAS format, then
   topped up from `RegistryID` for any that the synonyms missed. Skipped
   entirely with `--no-cas`.
5. **Export** — `compounds.csv` and `compounds.json`.

### Taxonomy (`fetch_taxonomy.py`)

1. **Collect IDs** — same file/argument handling as CIDs.
2. **Fetch** — `efetch.fcgi?db=taxonomy&id=...&retmode=xml`. Taxonomy IDs are
   primary keys, so no `esearch` step is needed (unlike BioProject).
3. **Parse** — scientific name, common name, rank, and the full lineage.
4. **Export** — `taxonomy.csv` (one row per taxon, with `domain` … `species`
   flattened into columns plus a `full_lineage` string) and `taxonomy.json`
   (the full nested record, including every lineage entry).

## Parameters

### `fetch_compounds.py`

| Parameter | Default | Description |
|---|---|---|
| `--cids` | required | CID list or path(s) to text files containing CIDs |
| `--output-dir` | `compound_results` | Output directory for the CSV/JSON results |
| `--no-synonyms` | off | Skip synonyms (each compound may have hundreds) |
| `--no-cas` | off | Skip CAS extraction |
| `--max-synonyms` | `20` | Maximum synonyms kept per record |
| `--timeout` | `60` | Per-request timeout in seconds |

### `fetch_taxonomy.py`

| Parameter | Default | Description |
|---|---|---|
| `--tax-ids` | required | Taxonomy ID list or path(s) to text files containing IDs |
| `--output-dir` | `taxonomy_results` | Output directory for the CSV/JSON results |
| `--email` | `$BIOPIPE_EMAIL` | Contact email (required by NCBI's usage policy) |
| `--api-key` | `$NCBI_API_KEY` | NCBI API key; raises the rate limit from 3 to 10 req/s |
| `--timeout` | `60` | Per-request timeout in seconds |

## Outputs

| File | Content |
|---|---|
| `compounds.csv` | One row per CID: CID, found, IUPAC name, formula, molecular weight, exact mass, XLogP, TPSA, charge, SMILES, ConnectivitySMILES, InChI, InChIKey, H-bond counts, rotatable bonds, heavy atoms, CAS, synonym count, synonyms |
| `compounds.json` | Same content as nested JSON |
| `taxonomy.csv` | One row per taxon: tax_id, scientific name, common name, rank, `domain`…`species`, lineage depth, full lineage string |
| `taxonomy.json` | Full nested records including every lineage entry with its rank and taxid |

## Validation

- **Check `found` before trusting a row.** `found=False` means PubChem has no
  compound with that CID. This is not a transient error — the identifier is
  wrong or refers to a compound that was removed.
- **A compound may legitimately have several CAS numbers.** Salts, hydrates
  and stereoisomers each get their own. The `CAS` column holds them
  semicolon-separated. Match on the primary one only if you are joining
  against a single-registry table.
- **CAS candidates are format-validated, not registry-validated.** The
  pattern check removes non-CAS registry IDs, but it cannot tell you whether a
  CAS number has been re-assigned. For critical identity work, cross-check
  the InChIKey.
- **`superkingdom` is empty in current data.** NCBI retired that rank;
  "cellular organisms" is now `cellular root` and the top-level split is
  `domain` (Eukaryota / Bacteria / Archaea). Use the `domain` column.
- **Not every taxon has a species rank.** Strain-level entries report
  `no rank` (e.g. *Escherichia coli* str. K-12 substr. MG1655). The `rank`
  column tells you what the entry actually is — do not assume `species`.
- **Species names are not stable identifiers.** *Escherichia coli*'s phylum
  was renamed from `Proteobacteria` to `Pseudomonadota`; the `tax_id` is
  stable, the name is not. Join on `tax_id`.
- **A missing common name is normal.** `GenbankCommonName` is only present
  for well-studied organisms; most taxa have none.
- **Taxonomy IDs are stable, names are not.** Never key a cross-study join on
  the scientific name.

## References

- `references/methodology.md` — API protocols, field semantics, the four
  PubChem and NCBI traps, CAS extraction logic
- `references/troubleshooting.md` — common failures and how to diagnose them
- PubChem PUG-REST documentation: <https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest>
- NCBI Taxonomy: <https://www.ncbi.nlm.nih.gov/taxonomy>
- NCBI E-utilities documentation: <https://www.ncbi.nlm.nih.gov/books/NBK25501/>
