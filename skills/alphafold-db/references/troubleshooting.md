# Troubleshooting: AlphaFold DB Retrieval

## Every download fails with HTTP 400

**Symptom**: every accession reports `Invalid identifier format. Please use a
UniProt accession or a supported AlphaFold DB ID.`

**Cause**: the identifier in the file is not what you think it is. Common
culprits:

- The file still contains the header row (`id`, `uniprot_id`, `Accession`).
  This skill auto-detects and skips a header, but only if the first line does
  **not** look like an accession.
- The file uses **Ensembl gene IDs** (`ENSG00000139618`) instead of UniProt
  accessions. AlphaFold DB indexes UniProt accessions; map via UniProt first.
- The file contains **isoform** identifiers (`P69905-1`, `P69905-2`). The API
  wants the canonical accession without the isoform suffix.
- The file contains trailing whitespace, tabs, or a BOM. Check with
  `head -c 32 ids.txt | xxd`.
- The file has Windows `CRLF` line endings mixed with the content — usually
  harmless, but combined with a BOM it breaks the first line.

**Fix**: verify the identifiers independently.

```bash
curl -s "https://alphafold.ebi.ac.uk/api/prediction/P69905" | head -c 200
```

If that returns a JSON array starting with `[`, the accession form is correct.

## Some accessions are missing from the results

**Symptom**: `共 N 条记录` is smaller than the number of IDs requested, and the
log says "无预测模型".

**Cause**: the protein genuinely has no AlphaFold model. The AlphaFold DB covers
UniProt reference proteomes; it does not cover every isoform, every predicted
gene, or unreviewed TrEMBL entries added after the last release.

**Fix**: cross-check against UniProt to distinguish "no model" from "wrong ID".

```bash
curl -s "https://rest.uniprot.org/uniprotkb/P69905?format=json" | head -c 200
```

If UniProt knows the accession but AlphaFold does not, the protein is simply
not covered yet. Do not treat it as an error.

## The downloaded file is HTML, not a structure

**Symptom**: the file exists but a PDB parser fails immediately; `head` shows
`<!DOCTYPE html>`.

**Cause**: an interstitial page — corporate proxy, WAF, or a cached error
response.

**Fix**:
1. Inspect the file directly to identify the blocker.
2. The script already validates `Content-Type` against PDB types, so this
   usually means the proxy returned a `text/plain` error page — check the
   first line of the file for the actual message.
3. For a manual retry, use `curl -L -o file.pdb <pdbUrl>` to see the raw
   response and headers.

## The URL contains `model_v6` but files 404

**Symptom**: manually constructed URLs like
`https://alphafold.ebi.ac.uk/files/AF-P69905-F1-model_v6.pdb` return 404.

**Cause**: the model version has moved on (e.g. to `model_v7`).

**Fix**: never construct the URL yourself. Query the API and take `pdbUrl`:

```bash
curl -s "https://alphafold.ebi.ac.uk/api/prediction/P69905" \
  | python -c "import json,sys; print(json.load(sys.stdin)[0]['pdbUrl'])"
```

If you must work offline with a stored URL, re-query the API after any
AlphaFold DB release announcement rather than assuming the version is stable.

## `KeyError: 'uniprotAccession'`

**Symptom**: a KeyError when parsing the API response.

**Cause**: the response was the error object `{"error": "..."}` rather than a
record — the identifier was malformed (see the 400 case above).

**Fix**: use this skill's `alphafold_client.fetch_predictions()`, which
normalises list / dict / error shapes. If you wrote your own parsing, always
check for the `error` key **before** treating a dict as a record.

## The download is extremely slow

**Cause**: AlphaFold DB is rate-limited, and there is no batch endpoint, so
retrieval is inherently serial.

**Fix**:
- Use `--metadata-only` first to confirm which accessions actually have
  models, then download only those.
- Apply `--min-plddt` to skip low-confidence models.
- Split large ID lists across several machines; the per-accession cost is
  dominated by the API round-trip, not by transfer size.

## A file has zero bytes or is truncated

**Symptom**: the file exists but is empty or cuts off mid-record.

**Cause**: the connection dropped mid-transfer and the file was still created.

**Fix**: delete the file and re-run with `--force`. The script writes the file
only after the full response body has been received, so an existing file should
normally be complete — but a manual `curl` with a truncated transfer can still
produce one.

To check integrity quickly:

```bash
tail -n 1 structures/P69905.pdb     # should be "END"
```

## The metadata CSV has empty pLDDT columns

**Cause**: the record was fetched from an older cached API response, or the
accession maps to a fragment record without a global metric.

**Fix**: check `model_created_date` in the CSV. If it predates the current
AlphaFold release, delete `alphafold_metadata.csv` and re-run to refresh.

## The first accession in my file is skipped

**Symptom**: one fewer structure than expected.

**Cause**: older scripts unconditionally dropped the first line of the ID file
assuming it was a header. This skill auto-detects headers by pattern — a first
line matching the accession format (letter followed by digits) is kept.

**Fix**: if your file really does have a header, the detection handles it. To
confirm what was parsed, run with logging enabled:

```bash
python scripts/download_structures.py --ids ids.txt --output-dir structures/ 2>&1 | head -5
```

## `requests` is not installed

**Symptom**: `requests 不可用，无法查询 AlphaFold DB`.

**Fix**: `pip install requests`. The script degrades gracefully (returns empty
rather than raising a traceback), so a missing dependency never crashes the
run — it just produces no results.
