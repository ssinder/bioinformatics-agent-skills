---
name: ena-download
description: >-
  Downloads RNA-seq and other sequencing FASTQ data from the ENA (European
  Nucleotide Archive) Portal API and exports run-level metadata. Resolves
  download URLs and file sizes from the structured API rather than scraping
  archive web pages, supports batch accession lists, HTTP Range resume for
  interrupted transfers, and MD5 verification against the checksums ENA
  publishes. This skill should be used when fetching raw sequencing reads by
  SRA/ERR/DRR accession, retrieving FASTQ files for an RNA-seq project, or
  collecting run metadata and read counts is requested.
license: MIT
compatibility: Python 3.9+; needs requests for network access; standard library only otherwise
metadata:
  version: "1.0"
  last-reviewed: "2026-10-26"
  skill-author: sinder
---

# ENA Sequencing Data Retrieval

## Overview

This skill retrieves raw sequencing data (FASTQ) from the ENA Portal API and
exports run-level metadata. It replaces a common but fragile pattern: scraping
archive web pages to extract download links and file sizes.

The ENA Portal API returns **structured** fields — `fastq_ftp`, `fastq_bytes`,
`fastq_md5`, `library_strategy`, `read_count`, `base_count` — for every run.
That removes three failure modes at once:

- **Fragile HTML parsing.** Page-scraping approaches depend on absolute XPath
  positions (`/html/body/div[2]/div[2]/div[2]/...`) that break whenever the
  archive changes its layout. The API's schema is stable.
- **Only getting `.sra`.** The archive's own download links point at the
  original submitted format (usually `.sra`), which then needs a separate
  `fastq-dump` conversion step. ENA already ships the converted FASTQ, with
  per-file sizes and MD5 checksums.
- **No integrity check.** A truncated download looks like a successful one
  until a downstream tool fails on malformed reads. ENA publishes the MD5 of
  every FASTQ; this skill verifies each file after transfer and deletes it if
  the checksum does not match.

Downloads are **resumable**: an interrupted transfer continues from the byte
it stopped at using an HTTP `Range` request, instead of starting over.

## When to Use This Skill

- Fetching FASTQ for a list of SRA/ERR/DRR accessions for an RNA-seq project
- Collecting run metadata (library strategy, read counts, organism) before
  deciding what to download
- Resuming a large download that was interrupted
- Verifying that downloaded FASTQ files are complete and uncorrupted
- Sizing a dataset before committing bandwidth (metadata-only mode)

## Quick Start

```bash
export PYTHONPATH="$PWD/packages"

# 批量下载
python scripts/download_runs.py --ids SRR390728 SRR1975568 --output-dir fastq/

# 从文件读取 accession 列表
python scripts/download_runs.py --ids runs.txt --output-dir fastq/

# 只查元数据不下载（先摸清数据规模）
python scripts/download_runs.py --ids runs.txt --metadata-only --output-dir meta/

# 只处理 RNA-Seq 文库，跳过 WGS 等其它类型
python scripts/download_runs.py --ids runs.txt --strategy RNA-Seq --output-dir fastq/

# 单端数据只取第一个 FASTQ
python scripts/download_runs.py --ids runs.txt --single-end --output-dir fastq/
```

## Workflow

1. **Collect accessions** — `--ids` accepts accession strings directly or
   paths to text files (one accession per line, `#` for comments). Unlike the
   original script, there is no separate "index" configuration file mixing
   paths and URLs together; the accession list is the only input.
2. **Query the API** — `ena_client.fetch_run_metadata()` builds an
   `run_accession=X OR run_accession=Y` query and posts it to the ENA
   `search` endpoint. (The `filereport` endpoint rejects multi-accession
   queries, so `search` is used for batching.)
3. **Export metadata** — every run is written to `ena_run_metadata.csv` and
   `.json`, including `library_strategy`, `read_count`, `base_count` and the
   FASTQ paths, sizes and MD5s.
4. **Filter (optional)** — `--strategy` keeps only runs with a given library
   strategy, e.g. `RNA-Seq`.
5. **Download** — each FASTQ is fetched as a stream. Existing files are
   skipped unless `--force` is given; partial files are continued with a
   `Range` header.
6. **Verify** — the file's MD5 is compared against the ENA-published value.
   A mismatch deletes the file (so it is not mistaken for complete on the
   next run) and counts the run as failed.
7. **Report** — counts of downloaded / skipped / failed files, plus the list
   of failed runs.

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `--ids` | required | Accession list or path(s) to text files containing accessions |
| `--output-dir` | `ena_fastq` | Download directory for FASTQ files and metadata |
| `--metadata-only` | off | Query and export metadata only; skip downloads |
| `--single-end` | off | Download only the first FASTQ of each run |
| `--strategy` | none | Keep only runs with this `library_strategy` (e.g. `RNA-Seq`) |
| `--force` | off | Re-download files that already exist |
| `--no-resume` | off | Disable resume; restart interrupted transfers from zero |
| `--no-verify` | off | Skip MD5 verification (not recommended) |
| `--timeout` | `600` | Per-request timeout in seconds |

## Outputs

| File | Content |
|---|---|
| `{SRRxxxxxxx}_1.fastq.gz` | Downloaded FASTQ (ENA's original filenames are kept) |
| `ena_run_metadata.csv` | One row per run: accessions, organism, instrument, library strategy/layout, read and base counts, FASTQ paths, sizes and MD5s |
| `ena_run_metadata.json` | Same content as nested JSON |

## Validation

- **The FASTQ is valid gzip** — `gzip -t file.fastq.gz` should pass. A file
  that fails this check was likely truncated; re-download with `--force`.
- **MD5 matches** — the metadata CSV records ENA's expected MD5 per file.
  The skill already verifies this, and deletes mismatching files.
- **File size matches `fastq_bytes`** — the metadata records ENA's byte
  count for each file. A size mismatch means the transfer was incomplete
  even if the MD5 was skipped with `--no-verify`.
- **Read count is plausible** — `read_count` in the metadata is the number of
  reads (not read *pairs*). For paired-end data the number of FASTQ files
  should be 2, and the total base count should be roughly
  `read_count × read_length`.
- **Library strategy is what you expect** — WGS runs also live in SRA/ENA.
  If you expected RNA-Seq and see `WGS`, the accession may belong to a
  different project; check `--strategy` filtering output.
- **Spot-check the first reads** — `zcat file.fastq.gz | head -4` should show
  a header line starting with `@` and a quality line starting with `+`.

## References

- `references/methodology.md` — API semantics, format choice, integrity verification strategy
- `references/troubleshooting.md` — common failures and how to diagnose them
- ENA Portal API: <https://www.ebi.ac.uk/ena/portal/api>
