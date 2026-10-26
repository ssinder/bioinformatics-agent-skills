# Methodology: ENA Data Retrieval

## 1. Why the API instead of web scraping

The ENA (European Nucleotide Archive) mirrors NCBI's SRA and DDBJ's DRA, and
exposes a queryable portal API at `https://www.ebi.ac.uk/ena/portal/api`.

A scraping-based workflow (the common pattern) has three structural problems:

**Layout coupling.** Archive pages are rendered for humans. Any selector that
works today — e.g. an absolute path like
`/html/body/div[2]/div[2]/div[2]/div/div/table/tbody/tr[5]/td/span[1]/span`
— encodes the current number of nested wrapper divs, the current table row
order and the current span structure. Any layout change breaks it, and the
failure mode is "file size column silently empty", not an exception.

**Wrong artifact.** The archive's own download buttons serve the *submitted*
format, which for most modern submissions is a `.sra` container holding
SRA-compressed reads. Getting usable FASTQ then requires running
`fastq-dump` (from the SRA Toolkit) on every file, which needs the toolkit
installed, is slow, and produces a second copy on disk.

**No checksum.** A page-scraped size is advisory text. Nothing verifies that
the bytes actually received match what was offered.

The API removes all three: it returns the converted FASTQ directly, with exact
byte counts and MD5 checksums, in a stable schema.

## 2. The two endpoints, and why `search` is used

ENA offers two query endpoints:

| Endpoint | Multi-accession | Use for |
|---|---|---|
| `filereport` | **No** — rejects `A;B` with an explanatory text body | single accession lookups |
| `search` | **Yes** — accepts a query expression | batch lookups (this skill) |

Requesting `filereport?accession=SRR390728;SRR1975568` returns:

```
Accession(s) SRR390728;SRR1975568 not valid for search requests on read_run data.
Supported accession types are: sample [...], study [...], submission [...], ...
```

Note the error text itself says "not valid for **search** requests" — the
endpoint routes `filereport` through the same search machinery but does not
accept a multi-value `accession` parameter.

The `search` endpoint instead takes a `query` expression:

```
GET /ena/portal/api/search
    ?result=read_run
    &query=run_accession=SRR390728 OR run_accession=SRR1975568
    &fields=run_accession,scientific_name,...
    &format=tsv
    &limit=0
```

`limit=0` means "no limit" — the default is 25, which silently truncates
large requests.

This skill batches at 50 accessions per query. ENA limits query-string
length, and larger batches mean fewer round-trips.

## 3. Response format

The API returns TSV. Two details matter:

**Empty results still include a header.** A query matching nothing returns a
single line containing only the column names, terminated by a newline. Naive
`split("\n")` + `split("\t")` parsing would produce one bogus row of all-empty
fields. This skill's `parse_tsv()` treats a header-only response as an empty
result.

**Multi-value fields are semicolon-separated.** `fastq_ftp`, `fastq_bytes` and
`fastq_md5` are *parallel arrays*: the n-th element of each refers to the
n-th file of the same run. This correspondence is what makes size and
checksum verification possible, and it is why the three fields must be split
together and zipped by index — never matched by name.

## 4. Choosing which files to download

ENA exposes two file families per run:

| Field | Content | When to use |
|---|---|---|
| `fastq_ftp` | FASTQ, already converted from the submitted format | **default** — directly usable by aligners and quantifiers |
| `submitted_ftp` | The original submitted files (often `.sra`) | only when the original format is needed (e.g. preserving the submitter's exact reads, or a format ENA did not convert) |

The `fastq_*` path is preferred because it removes the `fastq-dump` step.

### Single-end vs paired-end

`library_layout` reports `SINGLE` or `PAIRED`. For a paired run, ENA lists
two FASTQ files (`_1.fastq.gz`, `_2.fastq.gz`). `--single-end` downloads only
the first, which is useful for quick inspection but **not** for quantification —
most tools expect both mates.

A useful check: the number of `fastq_ftp` entries should match the layout
(1 for SINGLE, 2 for PAIRED). A run with 3 files is a multi-sample or
multi-lane run.

## 5. Integrity verification

ENA publishes an MD5 for every FASTQ file. The verification procedure is:

1. Stream the file to disk in 1 MiB chunks.
2. Compute the MD5 over the complete file.
3. Compare (case-insensitively) against the value from `fastq_md5`.
4. On mismatch, **delete the file** and report the run as failed.

Step 4 is the important one. If a corrupt file were left on disk, the next
run would see it as "already downloaded" and skip it — turning a transient
error into a permanent one. Deleting it makes the failure recoverable.

MD5 is used here purely as a transfer checksum (detecting truncation and
bit-rot), not as a security primitive; it is what the archive publishes.

## 6. Resumable transfers

Large FASTQ files frequently exceed what a single connection completes
reliably. The transfer is therefore resumable:

1. If the destination exists, its current size is the resume offset.
2. A `Range: bytes=<offset>-` header is sent.
3. The response body is **appended** to the file.
4. If the server replies `200 OK` instead of `206 Partial Content`, the server
   ignored the Range header — in that case the file must be rewritten from
   scratch, not appended to. Appending would corrupt it by duplicating the
   already-present prefix.

The final MD5 check catches any residual inconsistency, so a resume that
went wrong cannot pass silently.

## 7. Proxy and environment handling

The HTTP client is created through the shared `biopipe_core.http.http_session()`,
which by default **trusts environment variables** — `HTTPS_PROXY`,
`REQUESTS_CA_BUNDLE` and friends. This matters in two ways:

- Setting `trust_env = False` (a common "disable system proxy" pattern) also
  disables custom CA bundle configuration, which breaks TLS verification
  behind corporate proxies that use internal certificate authorities.
- In networks that require a proxy, disabling environment proxy settings
  makes the client simply unable to reach the archive.

The shared session also provides exponential-backoff retries on transient
status codes (429, 500, 502, 503, 504) and a minimum inter-request interval
of ~0.34 s, which keeps the request rate inside ENA's politeness envelope.

## 8. Rate and etiquette

ENA asks that tools stay within ~15 requests per second and identify
themselves. This skill:

- limits the request rate via the shared session's throttle;
- batches up to 50 accessions per query, so the request count is proportional
  to the *number of batches*, not the number of accessions;
- downloads files sequentially rather than in parallel. Parallel downloads of
  multi-GB files would compete for bandwidth and risk being throttled; the
  bottleneck for a large dataset is usually the transfer itself, not the
  request count.

## References

- ENA Portal API documentation: <https://www.ebi.ac.uk/ena/portal/api>
- ENA data model: <https://www.ebi.ac.uk/ena/browser/home>
- SRA Toolkit (`fastq-dump`, only needed for `.sra` files): <https://ncbi.github.io/sra-tools/>
- MD5 checksum practice for data integrity: <https://en.wikipedia.org/wiki/MD5>
