---
name: bioproject-search
description: "Use when looking up what a BioProject or PRJNA/PRD accession contains, retrieving sequencing project descriptions, or collecting BioProject metadata in bulk. Retrieves title, description, organism, submitter and release dates from the official NCBI E-utilities API rather than scraping web pages; supports keyword filtering with CSV/JSON export. 用于 BioProject 元数据检索与批量导出。"
license: MIT
compatibility: Python 3.9+; needs requests for network access; standard library only otherwise
metadata:
  version: "1.0"
  last-reviewed: "2026-10-26"
  skill-author: sinder
  edam-topic: "topic_3170"
  edam-operation: "operation_2421, operation_2422"
  edam-input-data: "data_0849"
  edam-output-data: "data_2080"
  edam-input-format: "format_3752"
  edam-output-format: "format_3752, format_3464"

# BioProject Metadata Retrieval


---

## Overview

Retrieves NCBI BioProject records through the official **E-utilities** API and exports them as CSV/JSON.

Scraping the BioProject web pages is fragile: selectors like `div:nth-child(5)` encode the *current* number of wrapper elements, so a layout change silently empties the description column — indistinguishable from "no description". E-utilities returns structured XML with stable field names, so the query survives a redesign (`references/methodology.md` §4).

NCBI's usage policy also applies: every caller must identify itself with a `tool` name and a **developer** email — explicitly *not* a third-party end user's address — registered in advance. This project ships that identity, so the scripts send `tool` and `email` automatically and **you do not need to configure anything**. An optional `NCBI_API_KEY` raises the rate limit from 3 to 10 req/s; override the identity with `--email` or `BIOAGENT_EMAIL` if your IT policy needs a local one. Every output directory also gets an `NCBI_DISCLAIMER.txt` (`references/methodology.md` §3).

## When to Use This Skill

- Finding out what a specific BioProject / PRJNA / PRJEB accession contains
- Collecting titles, descriptions and organisms for a list of projects
- Screening projects by keyword before deciding which to download
- Building a project catalogue for a methods or dataset section
- Looking up submitter and release dates for citation purposes

## Quick Start

```bash
export PYTHONPATH="$PWD/packages"
# No email setup needed — the project ships a registered identity

python scripts/search_bioprojects.py --ids prj_ids.txt --output-dir results/
python scripts/search_bioprojects.py --ids PRJNA172563 PRJNA212358 --output-dir results/
python scripts/search_bioprojects.py --ids prj_ids.txt --output-dir results/ --filter CRISPR
python scripts/search_bioprojects.py --ids prj_ids.txt --output-dir results/ --api-key $NCBI_API_KEY
```

## Workflow

1. **Collect accessions** — `--ids` takes accession strings or text files (one per line, `#` for comments). Entries that do not look like BioProject accessions are reported and skipped, not silently dropped.
2. **Resolve to UIDs** — `esearch` translates accessions into NCBI UIDs using the `[Project Accession]` qualifier; without it a string can match unrelated records.
3. **Fetch from two endpoints** — `esummary` and `efetch` are called together and merged. **No single endpoint carries every field**: `esummary` alone gives organism, taxonomy ID, project type and submitter organisation; `efetch` alone gives release date and description. Calling one silently drops half the metadata.
4. **Parse** — records carry accession, title, name, description, organism, taxonomy ID, project type, submitter, registration/release/last-update dates; the three NCBI date formats normalise to `YYYY-MM-DD`.
5. **Filter (optional)** — `--filter` keeps records whose title, name or description contains the keyword (case-insensitive); `--drop-noisy` drops navigation-like descriptions.
6. **Export** — `bioproject_data.csv` and `bioproject_data.json`.

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `--ids` | required | Accession list or path(s) to text files of accessions |
| `--output-dir` | `bioproject_results` | Output directory for the CSV/JSON results |
| `--email` | registered project identity | Developer contact email required by NCBI's usage policy; a registered identity ships with the project, so overriding is optional |
| `--api-key` | `$NCBI_API_KEY` | Optional NCBI API key; raises the rate limit from 3 to 10 req/s |
| `--filter` | none | Keep only records whose title/name/description contains this keyword |
| `--drop-noisy` | off | Drop descriptions that look like page navigation |
| `--timeout` | `60` | Per-request timeout in seconds |

## Outputs

| File | Content |
|---|---|
| `bioproject_data.csv` | One row per project: accession, title, name, description, organism, taxid, study type, submitter, registration date, release date, last-update date, UID |
| `bioproject_data.json` | Same content as nested JSON |

## Validation

- **A description should be a paragraph, not a fragment.** Empty or one-word descriptions usually mean the record genuinely lacks one — some submitters give only a title. Cross-check the linked BioSample or SRA.
- **The accession prefix indicates the archive** — `PRJNA` (NCBI), `PRJEB` (EBI/ENA), `PRJDB` (DDBJ); an ENA-submitted project may be reachable through ENA's portal.
- **Organism names come from the submitter** and are not normalised. A project may list several organisms, so the field holds a comma-separated list, not one taxon — do not read it as one taxid. `organism` holds a *name*; the numeric ID is in `taxid`.
- **An empty date means "not recorded", not "zero".** NCBI uses the placeholder `1/01/01 00:00` for "date unknown" and deliberately does not convert it, so the field stays empty. Do not read an empty `registration_date` as a bug — check whether the project is dbGaP-backed (those strip submitter details).
- **Release date may be in the future or absent** for embargoed projects, whose metadata returns to authorised submitters only; the public record stays sparse.
- **Search term semantics** — `[Project Accession]` matches the accession field exactly. To search descriptions by keyword across *all* projects rather than a known list, use NCBI Entrez with `bioproject[Filter]` plus a term.

## References

- `references/methodology.md` — E-utilities protocol, field semantics, rate limits
- `references/troubleshooting.md` — common failures and their diagnosis
