# Troubleshooting: ENA Data Retrieval

## No records returned for an accession

**Symptom**: the log says `ENA 未返回记录（N 个 accession）` and the run is
reported as failed.

**Cause**: the accession is not in ENA's mirror. ENA covers SRA (NCBI) and
DRA (DDBJ); some submissions exist only in the originating archive, and
records can be suppressed by the submitter.

**Fix**: verify the accession against NCBI SRA:

```bash
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRR390728&result=read_run&fields=run_accession&format=tsv"
```

An empty body (header only) means the accession is unknown to ENA. Check
whether the accession is correct — a common cause is a **run accession being
used where an experiment (SRX) or study (PRJNA) accession is meant**. ENA
accepts all of them, but the returned record set differs.

## The API returns an error message instead of TSV

**Symptom**: the script logs a warning and no records come back, but a manual
`curl` shows a plain-text explanation of supported accession types.

**Cause**: a query was malformed. ENA's `search` endpoint uses a
field-qualified expression syntax:

```
run_accession=SRR390728 OR run_accession=SRR1975568
```

not `run_accession:SRR390728` (colon) and not a comma-separated list.

**Fix**: this skill builds the expression for you. If you query manually,
use the `=` form with `OR` between clauses.

## A run has no FASTQ files

**Symptom**: `没有可用的 FASTQ（可能未转换或被删除）`.

**Cause**: ENA has not converted this run's submitted files to FASTQ. This
happens for some large WGS and long-read submissions where the original
format is kept as-is.

**Fix**: use the `submitted_ftp` field instead. This skill prefers FASTQ
because it avoids the `fastq-dump` step, but the submitted files are always
present. Query them directly:

```bash
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRR390728&result=read_run&fields=run_accession,submitted_ftp&format=tsv"
```

## MD5 verification fails

**Symptom**: `MD5 不匹配（期望 …，实际 …）` and the file is deleted.

**Causes and fixes**:

1. **Truncated transfer.** The connection dropped. Re-run — the file was
   deleted precisely so the next attempt starts clean.
2. **The file was modified in place.** Something else wrote to it. Re-run
   with `--force`.
3. **The archive updated the file.** Rare, but if the run was re-submitted
   the published MD5 changes. Re-fetch the metadata with `--force` to pick
   up the new checksum.

**If MD5 checks keep failing while the size matches**, the transfer path
(mirror, CDN) is likely serving a modified copy. Report it to ENA.

## Downloads are extremely slow or stall

**Cause**: the file is large and the single connection is the bottleneck.

**Fixes**:
- Use `--metadata-only` first to check `fastq_bytes` — a 10 GB file is a
  different commitment from a 200 MB one.
- Prefer FASTQ over `.sra`: the submitted `.sra` container is usually
  *larger* than the converted FASTQ, and requires an extra conversion step
  afterwards.
- If the network supports it, download several runs in parallel with
  separate processes, each handling a disjoint slice of the accession list.

## The download completes but the FASTQ is corrupt

**Symptom**: `gzip -t` fails, or the aligner reports malformed records.

**Fix**: the file was not fully written. Delete it and re-download:

```bash
rm fastq/SRR390728_1.fastq.gz
python scripts/download_runs.py --ids SRR390728 --output-dir fastq/ --force
```

If MD5 verification is enabled this should not happen — it means
`--no-verify` was used.

## A previous run left a partial file and it is being skipped

**Symptom**: the file exists but is much smaller than the metadata's
`fastq_bytes`, and the run is reported as skipped.

**Cause**: the "already exists" check only looks at file presence. If an
earlier run was interrupted and left a partial file, it looks complete.

**Fix**: re-run with `--force` to overwrite, or delete the partial files
manually. A more robust check would compare the on-disk size against
`fastq_bytes` — worth adding if partial files are common in your workflow.

## `requests` is not installed

**Symptom**: `requests 不可用，无法查询 ENA`.

**Fix**: `pip install requests`. The skill degrades gracefully rather than
raising a traceback, so a missing dependency produces an empty result
instead of a crash.

## Accessions are silently ignored

**Symptom**: you pass 10 accessions but only 7 are queried.

**Causes**:
- **Duplicates** — the list is deduplicated case-insensitively. `SRR390728`
  and `srr390728` are the same accession and count as one.
- **Comment lines** — lines starting with `#` are treated as comments.
- **Trailing text** — only the first whitespace-separated token of a line is
  used, so `SRR390728 human liver` parses as `SRR390728`.

**Fix**: check the log line `从文件读取 accession: …` followed by the query
count, and verify the list has no duplicates.

## Behind a corporate proxy the request fails

**Symptom**: connection errors or certificate verification failures.

**Fix**: this skill trusts environment variables by default. Set:

```bash
export HTTPS_PROXY=http://proxy.internal:3128
export REQUESTS_CA_BUNDLE=/path/to/internal-ca.pem
```

Do **not** set `trust_env = False` — that would disable both the proxy and
the custom CA bundle, which is exactly the configuration that fails in a
corporate network.
