---
name: bioproject-search
description: >-
  Retrieves NCBI BioProject metadata — title, description, organism, submitter and
  release dates — for a batch of BioProject accessions, using the official NCBI
  E-utilities API rather than scraping the BioProject web pages. Supports
  keyword filtering over titles and descriptions, CSV/JSON export, and requires
  a contact email as the NCBI usage policy mandates. This skill should be used when
  looking up what a BioProject or PRJNA/PRD accession contains, retrieving
  sequencing project descriptions, or collecting BioProject metadata in bulk is
  requested.
license: MIT
compatibility: Python 3.9+; needs requests for network access; standard library only otherwise
metadata:
  version: "1.0"
  last-reviewed: "2026-10-26"
  skill-author: sinder
---

# BioProject Metadata Retrieval

## Overview

This skill retrieves NCBI BioProject records through the official
**E-utilities** API and exports them as CSV/JSON.

The alternative — scraping the BioProject web pages — is fragile in a way that
is easy to underestimate. Page-scraping code ends up with selectors like
`div:nth-child(5)` that encode the *current* number of wrapper elements in
the page. When NCBI adjusts its layout, the selector silently stops matching
and the description column comes back empty, which is indistinguishable from
"this project has no description".

E-utilities returns structured XML with stable field names. The same query
that works today keeps working after a website redesign.

The skill also follows NCBI's usage policy properly: a **contact email is
required** (via `--email` or the `BIOPIPE_EMAIL` environment variable).
NCBI uses the email to identify callers and to contact them when problems
occur; fabricating placeholder addresses defeats that purpose and can get
your IP rate-limited.

## When to Use This Skill

- Finding out what a specific BioProject / PRJNA / PRJEB accession contains
- Collecting titles, descriptions and organisms for a list of projects
- Screening projects by keyword before deciding which ones to download
- Building a project catalogue for a methods or dataset section
- Looking up submitter and release dates for citation purposes

## Quick Start

```bash
export PYTHONPATH="$PWD/packages"
export BIOPIPE_EMAIL="you@your-institution.edu"

# 从文件读取 accession
python scripts/search_bioprojects.py --ids prj_ids.txt --output-dir results/

# 直接给 accession
python scripts/search_bioprojects.py --ids PRJNA172563 PRJNA212358 --output-dir results/

# 只保留描述里提到 CRISPR 的项目
python scripts/search_bioprojects.py --ids prj_ids.txt --output-dir results/ --filter CRISPR

# 有 NCBI API key 时可提高限速（10 req/s）
python scripts/search_bioprojects.py --ids prj_ids.txt --output-dir results/ --api-key $NCBI_API_KEY
```

## Workflow

1. **Collect accessions** — `--ids` accepts accession strings or paths to
   text files (one per line, `#` for comments). Entries that do not look like
   BioProject accessions are reported and skipped rather than silently
   dropped.
2. **Resolve to UIDs** — `esearch` translates accessions into NCBI's internal
   UIDs. The query uses the `[Project Accession]` field qualifier, without
   which an accession string can match unrelated records.
3. **Fetch details from two endpoints** — `esummary` and `efetch` are called
   together and their results merged. This is not redundancy: **no single
   endpoint carries every field**. `esummary` is the only source for the
   structured organism name, taxonomy ID, project type and submitter
   organisation; `efetch` is the only source for the release date and the
   full description text. Calling only one of them silently drops about
   half the metadata without raising an error.
4. **Parse** — the merged records carry accession, title, name, description,
   organism, taxonomy ID, project type, submitter, registration date, release
   date and last-update date. All three NCBI date formats are normalised to
   `YYYY-MM-DD`.
5. **Filter (optional)** — `--filter` keeps records whose title, name or
   description contains the keyword (case-insensitive). `--drop-noisy`
   removes descriptions that look like page navigation rather than content.
6. **Export** — `bioproject_data.csv` and `bioproject_data.json` are written
   to the output directory.

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `--ids` | required | Accession list or path(s) to text files containing accessions |
| `--output-dir` | `bioproject_results` | Output directory for the CSV/JSON results |
| `--email` | `$BIOPIPE_EMAIL` | Contact email (required by NCBI's usage policy) |
| `--api-key` | `$BIOPIPE_API_KEY` | NCBI API key; raises the rate limit from 3 to 10 req/s |
| `--filter` | none | Keep only records whose title/name/description contains this keyword |
| `--drop-noisy` | off | Drop descriptions that look like page navigation |
| `--timeout` | `60` | Per-request timeout in seconds |

## Outputs

| File | Content |
|---|---|
| `bioproject_data.csv` | One row per project: accession, title, name, description, organism, taxid, study type, submitter, registration date, release date, last-update date, UID |
| `bioproject_data.json` | Same content as nested JSON |

## Validation

- **A description should be a paragraph, not a sentence fragment.** Empty or
  one-word descriptions usually mean the project record genuinely lacks one —
  some submitters provide only a title. Cross-check with the linked BioSample
  or SRA records.
- **The accession prefix indicates the archive** — `PRJNA` (NCBI), `PRJEB`
  (EBI/ENA), `PRJDB` (DDBJ). A project submitted to ENA may be reachable
  through ENA's own portal rather than NCBI's.
- **Organism names come from the submitter** and are not normalised. A
  project may list several organisms, in which case the field contains a
  comma-separated list rather than a single taxon — do not treat it as a
  single taxid. The `organism` field always holds a *name*; the numeric
  identifier belongs in `taxid`.
- **An empty date means "not recorded", not "zero".** NCBI uses the
  placeholder `1/01/01 00:00` for "date unknown"; it is deliberately not
  converted into a real date, so the field stays empty. Do not interpret an
  empty `registration_date` as a bug — check whether the project is
  dbGaP-backed (those records have their submitter details stripped).
- **Release date may be in the future or absent** for embargoed projects.
  An embargoed project returns its metadata to authorised submitters only;
  the public record will be sparse.
- **Search term semantics** — `[Project Accession]` matches the accession
  field exactly. To search descriptions by keyword across *all* projects
  rather than a known list, use NCBI Entrez directly with
  `bioproject[Filter]` and a keyword term.

## Validation

- **A description should be a paragraph, not a sentence fragment.** Empty or
  one-word descriptions usually mean the project record genuinely lacks one —
  some submitters provide only a title. Cross-check with the linked BioSample
  or SRA records.
- **The accession prefix indicates the archive** — `PRJNA` (NCBI), `PRJEB`
  (EBI/ENA), `PRJDB` (DDBJ). A project submitted to ENA may be reachable
  through ENA's own portal rather than NCBI's.
- **Organism names come from the submitter** and are not normalised. A
  project may list several organisms, in which case the field contains a
  comma-separated list rather than a single taxon — do not treat it as a
  single taxid.
- **Release date may be in the future or absent** for embargoed projects.
  An embargoed project returns its metadata to authorised submitters only;
  the public record will be sparse.
- **Search term semantics** — `[Project Accession]` matches the accession
  field exactly. To search descriptions by keyword across *all* projects
  rather than a known list, use NCBI Entrez directly with
  `bioproject[Filter]` and a keyword term.

## References

- `references/methodology.md` — E-utilities protocol, field semantics, rate limits
- `references/troubleshooting.md` — common failures and how to diagnose them
- NCBI E-utilities documentation: <https://www.ncbi.nlm.nih.gov/books/NBK25501/>
- BioProject: <https://www.ncbi.nlm.nih.gov/bioproject/>
