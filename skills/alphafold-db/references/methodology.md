# Methodology: AlphaFold DB Retrieval

## 1. The API, precisely

### Endpoint

```
GET https://alphafold.ebi.ac.uk/api/prediction/{accession}
Accept: application/json
```

`{accession}` accepts a UniProt accession (e.g. `P69905`) or an AlphaFold DB ID
(e.g. `AF-P69905-F1`).

### Response shape — the part that breaks naive scripts

For a **single** accession the endpoint returns a **JSON array**, not an object:

```json
[
  {
    "entryId": "AF-P69905-F1",
    "uniprotAccession": "P69905",
    "uniprotId": "HBA_HUMAN",
    "uniprotDescription": "Hemoglobin subunit alpha",
    "gene": "HBA1",
    "organismScientificName": "Homo sapiens",
    "taxId": "9606",
    "chainId": "A",
    "modelCreatedDate": "2025-08-01T00:00:00Z",
    "latestVersion": 6,
    "globalMetricValue": 98.06,
    "fractionPlddtVeryLow": 0.0,
    "fractionPlddtLow": 0.0,
    "fractionPlddtConfident": 0.0146,
    "fractionPlddtVeryHigh": 0.9854,
    "sequence": "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR",
    "pdbUrl": "https://alphafold.ebi.ac.uk/files/AF-P69905-F1-model_v6.pdb",
    "cifUrl": "https://alphafold.ebi.ac.uk/files/AF-P69905-F1-model_v6.cif",
    "bcifUrl": "https://alphafold.ebi.ac.uk/files/AF-P69905-F1-model_v6.bcif",
    "paeImageUrl": "https://alphafold.ebi.ac.uk/files/AF-P69905-F1-predicted_aligned_error.json",
    "isReviewed": true,
    "isComplex": false
  }
]
```

Three distinct shapes must all be handled:

| Shape | When | Correct handling |
|---|---|---|
| `[ {...} ]` | normal single-accession query | iterate the list |
| `{ ... }` | occasional direct-object response | treat as one record |
| `{"error": "Invalid identifier format..."}` | malformed / unknown identifier | **not** a record — must be distinguished before falling into the dict branch |

The error object is the trap: it *is* a dict, so a naive `isinstance(payload, dict)`
check would treat the error message as metadata. This skill's `_iter_entries()`
requires `uniprotAccession` to be present before accepting a bare dict as a record.

### Batch queries are not supported

`GET /api/prediction/P69905,P04637` returns:

```json
{"error": "Invalid identifier format. Please use a UniProt accession or a supported AlphaFold DB ID."}
```

There is no comma-separated batch endpoint, so retrieval is necessarily one request
per accession. Rate limiting is therefore the binding constraint, not bandwidth.

### Not found

An accession with no prediction returns **HTTP 400** with the same
`Invalid identifier format` error body — note that this is a 400, *not* a 404, even
though the identifier may be syntactically valid. Accessions absent from the
AlphaFold DB simply have no model.

## 2. Why model versions must not be hard-coded

The AlphaFold DB naming scheme is:

```
AF-{ACCESSION}-F{FRAGMENT}-model_{VERSION}.{ext}
```

`{VERSION}` increments when the database retrains and republishes its models.
A script that constructs `AF-{ACCESSION}-F1-model_v6.pdb` therefore:

- works today,
- returns HTTP 404 the day the next model version is published,
- fails in a way that looks like "the protein has no structure" rather than
  "the script has a stale version constant".

The robust strategy is:

1. Query `/api/prediction/{accession}`.
2. Read `pdbUrl` / `cifUrl` / `bcifUrl` from the response.
3. Download from whichever of those the user asked for.

Because the URL comes from the server, it tracks version bumps automatically. This
skill additionally implements a *fallback* URL builder that probes a short list of
candidate versions — but it is only used when the API itself is unreachable, and
never overrides the server-provided URL.

### Fragment suffix `-F1`

`{FRAGMENT}` is `F1` for single-fragment (single-chain) predictions. Multi-fragment
models use `F2`, `F3`, … The API returns one entry per fragment, with `chainId`
identifying the chain. This skill takes the first entry as the primary record and
logs a debug message when more are present.

## 3. Choosing a coordinate format

| Format | Endpoint field | When to use |
|---|---|---|
| PDB | `pdbUrl` | Maximum compatibility: every downstream tool reads it |
| mmCIF | `cifUrl` | Preferred for large models — PDB's fixed columns cannot represent >9999 residues or >62 chains |
| BinaryCIF | `bcifUrl` | Compact; needs a BinaryCIF-capable reader |

For proteins whose sequence exceeds ~2700 residues, PDB format becomes unreliable and
mmCIF is the safer download target. The AlphaFold DB itself serves all three from the
same prediction.

## 4. pLDDT interpretation

pLDDT (predicted local distance difference test) is AlphaFold's per-residue
confidence, on a 0–100 scale. The API exposes:

- `globalMetricValue` — the mean pLDDT over the whole model.
- `fractionPlddtVeryLow` (< 50) — residues considered unreliable, typically
  intrinsically disordered regions.
- `fractionPlddtLow` (50–70) — backbone generally fine, side chains unreliable.
- `fractionPlddtConfident` (70–90) — typical well-ordered backbone.
- `fractionPlddtVeryHigh` (> 90) — very high confidence, usually folded cores.

Practical thresholds:

| Global pLDDT | Interpretation |
|---|---|
| > 90 | Folded, high confidence — safe to use for docking/structural work |
| 70–90 | Generally reliable; inspect low-confidence fractions |
| 50–70 | Caution — large flexible or disordered regions |
| < 50 | Low confidence overall; the model may be mostly disorder |

The four confidence fractions sum to 1.0. That invariant is a cheap self-check that
the metadata record was parsed completely.

**Filtering before download** is the practical win: metadata queries are small JSON
responses, while a PDB file for a 500-residue protein is roughly 500 KB. Screening a
few thousand proteomes on metadata alone is orders of magnitude cheaper than
downloading first and discarding afterwards. This is what `--min-plddt` does.

## 5. Download correctness

Three details separate a working downloader from a subtly broken one:

**Binary write.** The response body must be written with `write_bytes()` /
`response.content`. Decoding to `response.text` and writing in text mode applies
platform newline translation, which rewrites `\n` inside the PDB's fixed-column
records. On Windows this silently changes every line ending of the file.

**Content-Type validation.** An error page returned with `text/html` is not a
structure. Checking the response `Content-Type` before writing catches proxy
interstitials, WAF blocks and captive portals, which otherwise produce a file that
looks plausible to a directory listing but fails at the first `ATOM` parse.

**Resumability.** Skipping files that already exist makes a batch run restartable
after a network interruption. `--force` overrides this when a refresh is intended
(for example after AlphaFold DB publishes a new model version).

## 6. Naming convention

The API's original filename embeds the model version
(`AF-P69905-F1-model_v6.pdb`), which is informative but noisy when a project
accumulates structures from several accessions. The default output name is
`{accession}.pdb`, which is stable across model versions and sorts naturally.
Use `--keep-original-name` when the version in the filename needs to be visible —
for instance to keep two vintages of the same protein side by side.

## References

- AlphaFold DB API: <https://alphafold.ebi.ac.uk/api>
- AlphaFold DB documentation: <https://alphafold.ebi.ac.uk/>
- AlphaFold paper (DeepMind, 2021): <https://doi.org/10.1038/s41586-021-03819-2>
- pLDDT and PAE interpretation: <https://alphafold.ebi.ac.uk/faq>
