---
name: ena-download
description: "Use when fetching raw sequencing reads by SRA/ERR/DRR accession, retrieving FASTQ files for an RNA-seq project, or collecting run metadata and read counts. Downloads FASTQ from the ENA Portal API, resolving URLs and sizes from the structured API rather than scraping archive pages; supports batch accession lists, HTTP Range resume, and MD5 verification. 用于下载原始 reads 与运行元数据。"
license: MIT
compatibility: Python 3.9+; needs requests for network access; standard library only otherwise
metadata:
  version: "1.0"
  last-reviewed: "2026-10-26"
  skill-author: sinder
  edam-topic: "topic_3170"
  edam-operation: "operation_2422, operation_3921"
  edam-input-data: "data_0849"
  edam-output-data: "data_0850"
  edam-input-format: "format_3752"
  edam-output-format: "format_1930, format_3752"

# ENA Sequencing Data Retrieval


---

## Overview

Retrieves raw sequencing data (FASTQ) from the ENA Portal API and exports run-level metadata, replacing the fragile pattern of scraping archive pages for download links and file sizes.

The ENA Portal API returns **structured** fields — `fastq_ftp`, `fastq_bytes`, `fastq_md5`, `library_strategy`, `read_count`, `base_count` — for every run. That removes three failure modes at once:

- **Fragile HTML parsing.** Page-scraping depends on absolute XPath positions (`/html/body/div[2]/div[2]/div[2]/...`) that break when the archive changes its layout. The API's schema is stable.
- **Only getting `.sra`.** The archive's own links point at the original submitted format (usually `.sra`), which then needs a separate `fastq-dump` conversion. ENA already ships the converted FASTQ, with per-file sizes and MD5 checksums.
- **No integrity check.** A truncated download looks successful until a downstream tool fails on malformed reads. ENA publishes the MD5 of every FASTQ; this skill verifies each file and deletes it if the checksum does not match.

Downloads are **resumable**: an interrupted transfer continues from the byte it stopped at via an HTTP `Range` request.

## When to Use This Skill

- Fetching FASTQ for a list of SRA/ERR/DRR accessions for an RNA-seq project
- Collecting run metadata (library strategy, read counts, organism) before deciding what to download
- Resuming a large download that was interrupted
- Verifying that downloaded FASTQ files are complete and uncorrupted
- Sizing a dataset before committing bandwidth (metadata-only mode)

## Quick Start

```bash
export PYTHONPATH="$PWD/packages"

python scripts/download_runs.py --ids SRR390728 SRR1975568 --output-dir fastq/
python scripts/download_runs.py --ids runs.txt --output-dir fastq/
python scripts/download_runs.py --ids runs.txt --metadata-only --output-dir meta/
python scripts/download_runs.py --ids runs.txt --strategy RNA-Seq --output-dir fastq/
python scripts/download_runs.py --ids runs.txt --single-end --output-dir fastq/
```

## Workflow

1. **Collect accessions** — `--ids` accepts accession strings or text files (one accession per line, `#` for comments). Unlike the original script, there is no separate "index" configuration file mixing paths and URLs; the accession list is the only input.
2. **Query the API** — `ena_client.fetch_run_metadata()` builds an `run_accession=X OR run_accession=Y` query and posts it to the ENA `search` endpoint. (The `filereport` endpoint rejects multi-accession queries, so `search` does the batching.)
3. **Export metadata** — every run goes to `ena_run_metadata.csv` and `.json`, with `library_strategy`, `read_count`, `base_count` and the FASTQ paths, sizes and MD5s.
4. **Filter (optional)** — `--strategy` keeps only runs with a given library strategy, e.g. `RNA-Seq`.
5. **Download** — each FASTQ is fetched as a stream. Existing files are skipped unless `--force`; partial files continue with a `Range` header.
6. **Verify** — the file's MD5 is compared against the ENA-published value. A mismatch deletes the file (so it is not mistaken for complete next run) and counts the run as failed.
7. **Report** — downloaded / skipped / failed counts, plus the failed runs.

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
| `ena_run_metadata.csv` | One row per run: accessions, organism, instrument, library strategy/layout, read and base counts, FASTQ paths, sizes, MD5s |
| `ena_run_metadata.json` | Same content as nested JSON |

## Validation

- **The FASTQ is valid gzip** — `gzip -t file.fastq.gz` should pass. A file that fails was likely truncated; re-download with `--force`.
- **MD5 matches** — the metadata CSV records ENA's expected MD5 per file; the skill verifies this and deletes mismatching files.
- **File size matches `fastq_bytes`** — the metadata records ENA's byte count per file. A size mismatch means an incomplete transfer even if the MD5 was skipped via `--no-verify`.
- **Read count is plausible** — `read_count` is the number of reads (not read *pairs*). For paired-end data there should be 2 FASTQ files, and the total base count should be roughly `read_count × read_length`.
- **Library strategy is what you expect** — WGS runs also live in SRA/ENA. If you expected RNA-Seq and see `WGS`, the accession may belong to another project; check `--strategy` filtering output.
- **Spot-check the first reads** — `zcat file.fastq.gz | head -4` should show a header line starting with `@` and a quality line starting with `+`.

## References

- `references/methodology.md` — API semantics, format choice, integrity verification strategy
- `references/troubleshooting.md` — common failures and how to diagnose them
