# Troubleshooting: BioDatabase Reference Retrieval

## Compound retrieval

### The SMILES column is empty

PubChem renamed `CanonicalSMILES` to `SMILES`. If you request
`CanonicalSMILES` and read `CanonicalSMILES` back, the field is missing.

```bash
# wrong — response key is "SMILES", not "CanonicalSMILES"
curl "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/property/CanonicalSMILES/JSON"

# right
curl "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/property/SMILES,ConnectivitySMILES/JSON"
```

This skill requests both `SMILES` and `ConnectivitySMILES`, so the columns
should be populated. If they are empty for a *specific* compound but populated
for others, that compound genuinely has no standard structure — check the
`found` column.

### A compound that does not exist still "succeeded"

```
GET /compound/cid/999999999/property/MolecularWeight/JSON
→ HTTP 200, {"PropertyTable":{"Properties":[{"CID":999999999}]}}
```

PubChem returns HTTP 200 with a body containing only the CID. The status code
tells you nothing. Check the `found` column in `compounds.csv`; `found=False`
means the CID is unknown to PubChem.

Common reasons: a typo in the CID, a CID from a different registry (CID and
CAS and InChIKey are different identifier spaces), or a compound that was
removed/merged.

### The CAS column is empty but the compound exists

Possible causes, in order of likelihood:

1. **`--no-cas` was passed.** Check the invocation.
2. **The synonyms request failed** while the properties request succeeded.
   The client logs a warning on failure and continues, so the run still
   exits 0. Look for `PubChem 同义词查询失败` in the log.
3. **The compound genuinely has no CAS number** in PubChem. Rare, but it
   happens for obscure entries.
4. **The compound's synonyms list is empty** even though the compound
   exists — check `synonym_count` in the output.

### The CAS column contains something that is not a CAS number

`RegistryID` values are mixed registry identifiers. Only entries matching
`^\d{2,7}-\d{2}-\d$` are treated as CAS numbers. If you see entries like
`CHEBI:15365` or `DB00945` in the CAS column, the filter was bypassed — most
likely by reading the `RegistryID` endpoint directly rather than the
extracted `CAS` column.

Do not use `xrefs/RegistryID` as a CAS source. It is a cross-reference dump
across hundreds of registries.

### A compound has several CAS numbers

Expected. Salts, hydrates and stereoisomers carry distinct registry numbers.
For CID 2244 the skill finds three: `50-78-2` (aspirin), `98201-60-6` and
one more for the related forms.

When joining against a table with a single CAS column, decide explicitly
which one is the primary — the first in the list is the one PubChem lists
first, which is usually the anhydrous form, but that is a convention, not a
guarantee.

### Requests are slow

PubChem's limit is 5 requests/second and 400/minute. The skill throttles to
0.25 s between requests, so 1,000 CIDs take roughly 30 s (10 property
requests plus synonym and registry requests).

If it is much slower, the likely causes are:

- `--no-synonyms` not passed — the synonyms endpoint returns large payloads
  (aspirin alone is ~700 synonyms).
- A very large CID list with a low timeout causing retries.
- Network issues causing the retry logic to kick in.

### The request fails with a connection error

The client retries with exponential backoff (4 attempts). Persistent failure
usually means the network is unreachable rather than the API being down —
PubChem itself is highly available. Check whether you need a proxy:

```bash
export https_proxy=http://127.0.0.1:7897
export http_proxy=http://127.0.0.1:7897
```

### `found` is `False` for every row

Every CID came back with no properties. Check the request URL: a malformed
property list (a trailing comma, an unknown property name) makes PubChem
return an error payload that the client parses as empty. Verify with curl:

```bash
curl -s "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/property/MolecularWeight/JSON"
```

## Taxonomy retrieval

### "必须提供 email" / email missing

NCBI requires a contact email on every request. Supply it with `--email` or:

```bash
export BIOPIPE_EMAIL="you@your-institution.edu"
```

Do not substitute a placeholder address. The field exists so that NCBI can
reach a real person when a request misbehaves; a fake address defeats that
purpose and breaches NCBI's terms of use.

### The output contains far more rows than IDs requested

This is the `root.iter("Taxon")` trap. `iter()` is recursive and descends into
each record's `LineageEx`, collecting every ancestor as if it were a separate
taxon. A request for 3 species yields 75 rows: the 3 species plus 72
ancestors named `cellular organisms`, `Eukaryota`, `Opisthokonta` and so on.

The parser must iterate the root's **direct children** only:

```python
# wrong
for taxon in root.iter("Taxon"):
    ...

# right
for taxon in root.findall("Taxon"):
    ...
```

If you are writing your own parser, this is the single most likely bug.

### The `species` column is empty

`LineageEx` contains only ancestors — it stops at `genus` for most taxa. The
species level comes from the top-level `<Rank>` and `<ScientificName>`
elements. A parser that reads only `LineageEx` produces an empty species
column and a lineage that stops one level too high.

Note this also applies to the record's own rank: a taxon at `no rank` (strain
level) has no species entry to fall back on, and `species` stays empty for a
legitimate reason.

### The `superkingdom` column is empty

Expected on current data. NCBI retired the `superkingdom` rank — "cellular
organisms" is now `cellular root`, and the top-level split is expressed as
`domain` (`Eukaryota` / `Bacteria` / `Archaea`). Use the `domain` column.

### The `domain` column is empty

Check that the taxon actually has a domain-level ancestor. Most do, but
records queried at a very high level (e.g. taxid 131567, "cellular
organisms") legitimately have no domain ancestor.

### The `common_name` column is empty

`GenbankCommonName` is only populated for well-studied organisms. Most taxa
have no common name — this is normal, not a failure.

### `rank` says `no rank`

Strain-level entries are reported as `no rank`. For example *Escherichia coli*
str. K-12 substr. MG1655 (taxid 511145). This is correct behaviour: NCBI
does not assign species rank to a strain.

### The lineage depth looks too shallow

`lineage_depth` counts the entries in `LineageEx` plus the taxon itself. A
deep vertebrate has ~30; a yeast has ~12. If it looks truncated, check
whether you are reading `LineageEx` from the right record — remember that
`iter()` over the root mixes records and lineages together.

### The scientific name changed between runs

NCBI renames clades. *E. coli*'s phylum was `Proteobacteria` and is now
`Pseudomonadota`. The taxid is stable; the name is not. Join on `tax_id`,
never on the name.

### Requests are rate-limited (HTTP 429)

The client throttles to 3 req/s, which is NCBI's limit without an API key.
If you are hitting 429:

1. Add an API key: `--api-key $NCBI_API_KEY` (raises the limit to 10 req/s).
2. Make sure the email is set — requests without an email are throttled
   more aggressively.
3. Check whether another process is hitting NCBI from the same IP.

## Output problems

### The CSV opens with a strange first column name

The CSV is written with `utf-8-sig` (a BOM) so that Excel on Windows
detects UTF-8 correctly. If you are parsing it with a tool that does not strip
the BOM, the first column appears as `\ufeffCID`. Python's `csv` module and
pandas handle this correctly; `pd.read_csv(..., encoding="utf-8-sig")` is the
explicit form.

### The output directory was created somewhere unexpected

`--output-dir` is resolved relative to the current working directory. Pass an
absolute path if you need a specific location; the skill logs the resolved
absolute path on completion.

## Still stuck

1. Run with `--timeout 120` — the default 60 s can be tight behind a slow
   proxy.
2. Check the log for `PubChem 查询失败` / `NCBI Taxonomy 查询失败`, which
   carry the underlying exception message.
3. Reproduce a single identifier with `curl` — it separates "the API
   returned something unexpected" from "the client has a bug".
