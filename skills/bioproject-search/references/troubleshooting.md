# Troubleshooting: BioProject Retrieval

## "必须提供 --email"

**Symptom**: the script exits with code 2 and this message.

**Cause**: NCBI's usage policy requires callers to identify themselves with a
valid email address. There is no default, by design.

**Fix**:

```bash
export BIOPIPE_EMAIL="you@your-institution.edu"
python scripts/search_bioprojects.py --ids prj_ids.txt --output-dir results/
```

or pass `--email` directly. If you have an NCBI API key, set
`BIOPIPE_API_KEY` too — it raises the rate limit from 3 to 10 requests per
second.

## HTTP 429 (Too Many Requests)

**Symptom**: requests start failing with 429, especially on a large list.

**Cause**: the rate limit was exceeded. Without an API key the limit is 3
requests/second.

**Fix**:
- Get an NCBI API key (free, from the account settings page) and pass
  `--api-key`.
- The shared HTTP client already throttles to ~3 requests/second and retries
  with exponential backoff, so a single run should not normally hit this.
  Sustained 429s usually mean several runs are executing in parallel.

## No records returned

**Symptom**: `没有查到任何 BioProject 记录`.

**Causes and fixes**:

1. **The accession does not exist.** Verify directly:

   ```bash
   curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=bioproject&term=PRJNA172563%5BProject%20Accession%5D&retmode=json"
   ```

   An empty `idlist` confirms it.

2. **Wrong accession type.** A `PRJNA` is a *project*. Passing a BioSample
   (`SAMN…`) or SRA run (`SRR…`) returns nothing. This skill's
   `--ids` validation reports entries that do not look like project
   accessions and skips them with a warning.

3. **The project is in another archive.** `PRJEB` projects are submitted to
   EBI and `PRJDB` to DDBJ. NCBI mirrors most of them, but very recent
   submissions may not be mirrored yet — query ENA's portal instead:

   ```bash
   curl -s "https://www.ebi.ac.uk/ena/portal/api/search?result=study&query=study_accession=PRJEB123&fields=study_accession,study_title&format=tsv"
   ```

4. **The project is embargoed.** Metadata for embargoed projects is
   restricted until release.

## The description column is empty

**Cause**: the submitter supplied only a title, or the description is stored
in a secondary `<Description>` element.

**Fix**: the parser joins multiple description elements. If it is still empty,
the record genuinely has no description — check the linked BioSample or SRA
records, which often carry more context.

## The organism column contains a number

**Symptom**: `organism` shows `28901` instead of `Salmonella enterica`.

**Cause**: the parser is reading the `species` **attribute** of the
`<Organism>` element. In NCBI's XML that attribute holds the *species-level
taxid*, not a name:

```xml
<Organism species="28901" taxID="90370">
  <OrganismName>Salmonella enterica subsp. enterica serovar Typhi</OrganismName>
</Organism>
```

**Fix**: the name is in the `<OrganismName>` child element. The current
parser reads it from there; if you see this symptom, you are running an
older build of the script. A quick check: `organism` should never be a bare
integer — the numeric identifier belongs in the `taxid` column.

## Every date column is empty

**Symptom**: `registration_date`, `release_date` and `last_update_date` are
all blank for every record.

**Cause**: the parser is looking for date *elements* when NCBI puts them in
XML *attributes*:

```xml
<Submission last_update="2013-10-23" submitted="2013-07-16" submission_id="SUB365098">
```

A parser that only walks child elements finds no `<submitted>` node and
reports an empty date.

**Fix**: read `Submission/@submitted` and `Submission/@last_update` as
attributes. The release date is a real element
(`Project/ProjectDescr/ProjectReleaseDate`) but only exists for records that
have been released — an empty release date on some rows is normal.

## The organism, project type or submitter column is empty

**Symptom**: `organism`, `taxid`, `study_type` or `submitter` is blank while
title and description are fine.

**Cause**: the parser only called `efetch`. Those four fields live in the
`esummary` response, not in `efetch` — the two endpoints carry disjoint
subsets of the metadata.

**Fix**: call both endpoints and merge them (`merge_records()` in
`bioproject_client.py` does this). The rule is efetch-first with esummary
filling the gaps.

## A date shows as 1/01/01 or the field is blank for a dbGaP project

**Symptom**: `registration_date` is empty for a project you know has a
submission date.

**Cause**: NCBI uses the placeholder `1/01/01 00:00` for "date unknown".
Projects hosted in dbGaP (like `PRJNA172563`) have their submitter
information stripped, and the public record frequently carries only the
placeholder date.

**Fix**: nothing to fix — the date is genuinely not public. The parser
deliberately does not convert the placeholder into `0001-01-01`, because a
fake date is worse than an empty cell.

## A description contains HTML tags

**Symptom**: the description shows `<i>Homo sapiens</i>` or similar.

**Cause**: the submitter wrote markup in the description field. E-utilities
returns it verbatim.

**Fix**: strip it downstream if plain text is required:

```python
import re
text = re.sub(r"<[^>]+>", "", description)
```

## Some accessions are skipped with a warning

**Symptom**: `以下条目不像 BioProject accession，已跳过: …`.

**Cause**: the list contains something that is not a `PRJNA`/`PRJEB`/`PRJDB`
accession — often a header line like "ID" or a note.

**Fix**: correct the list. The validation pattern is
`^PRJ(?:NA|EB|DB)\d+$` (case-insensitive), so `PRJNA172563` passes and
`PRJNA` or `172563` does not.

## The output CSV opens with garbled text in Excel

**Symptom**: non-ASCII characters (species names with diacritics, Chinese
project titles) appear as mojibake.

**Cause**: Excel's default CSV encoding on Windows is not UTF-8.

**Fix**: import the CSV with UTF-8 explicitly (Excel: *Data → From Text/CSV*
→ set the origin to UTF-8), or read it with the JSON version.

## An error mentions "UID" or "efetch"

**Symptom**: the request returns an error about invalid identifiers.

**Cause**: the accession was passed to `efetch` without going through
`esearch` first, or the UID list is malformed.

**Fix**: use `search_bioprojects()` rather than calling `efetch` directly.
The function performs the whole flow — `esearch` for the UID conversion,
then `esummary` + `efetch` for the details.

## A query returns thousands of records

**Symptom**: a keyword query returns far more projects than expected.

**Cause**: `esearch` defaults to 25 results, but the skill's request sets
`retmax` based on the input size, so a *known accession list* is exact. A
free-text keyword query would need pagination.

**Fix**: for keyword sweeps across all projects, use `retstart`/`retmax`
pagination on `esearch`, and add field qualifiers to narrow the scope
(`[Organism]`, `[Description]`, `[Project Data Type]`).

## Connection errors behind a proxy

**Symptom**: timeouts or TLS errors in a corporate network.

**Fix**: the shared HTTP client trusts environment variables by default. Set:

```bash
export HTTPS_PROXY=http://proxy.internal:3128
export REQUESTS_CA_BUNDLE=/path/to/internal-ca.pem
```

Do not disable environment handling — that would disable the proxy *and* the
custom CA bundle, which is exactly the configuration that fails in a
corporate network.
