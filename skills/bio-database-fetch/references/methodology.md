# Methodology: BioDatabase Reference Retrieval

## 1. PubChem PUG-REST

### The protocol

PUG-REST is a path-style REST API. The general shape is:

```
/rest/pug/<domain>/<namespace>/<identifiers>/<operation>/<output>
```

Two operations matter here:

| Operation | Path | Returns |
|---|---|---|
| properties | `/compound/cid/<cids>/property/<props>/JSON` | `PropertyTable.Properties[]` |
| synonyms | `/compound/cid/<cids>/synonyms/JSON` | `InformationList.Information[].Synonym[]` |
| registry IDs | `/compound/cid/<cids>/xrefs/RegistryID/JSON` | `InformationList.Information[].RegistryID[]` |

The GET interface accepts **up to 100 identifiers per request**. A request for
1,000 CIDs is therefore 10 HTTP calls rather than 1,000 — this is the single
biggest practical difference from a naive implementation, which usually
loops over identifiers one at a time.

```
cids.txt (1000 lines)
        │
        ├── chunk 1: 2244,3672,1983,... (100)  → GET /compound/cid/2244,3672,.../property/.../JSON
        ├── chunk 2: ...                        → GET ...
        └── chunk 10: ...                       → GET ...
                ↓
        PropertyTable.Properties[] → merge by CID → rows with found flag
```

The identifier list is a comma-separated path segment, which is why the CID
values are joined directly into the URL rather than passed as a query
parameter.

### Trap 1 — `CanonicalSMILES` comes back as `SMILES`

PubChem renamed the property. Requesting `CanonicalSMILES` and reading
`CanonicalSMILES` out of the response yields an empty value. The current
mapping is:

| Requested property | Key in the response |
|---|---|
| `CanonicalSMILES` | `SMILES` |
| `IsomericSMILES` | `SMILES` (with `ConnectivitySMILES` also present) |
| `ConnectivitySMILES` | `ConnectivitySMILES` |

This skill requests **both** `SMILES` and `ConnectivitySMILES` so the
connectivity-only (stereo-free) form is always available alongside the
stereo-aware one. A wrapper library that predates the rename will
silently produce empty SMILES columns.

### Trap 2 — a non-existent CID returns HTTP 200

```
GET /compound/cid/999999999/property/MolecularWeight/JSON

HTTP/1.1 200 OK
{"PropertyTable": {"Properties": [{"CID": 999999999}]}}
```

The request **succeeds**. The body contains the CID and nothing else. A
status-code check reports a successful lookup of a compound that does not
exist; a naive loop that appends whatever the API returned will produce a row
full of empty properties with no indication that anything went wrong.

The reliable test is the **shape of the response object**: a real compound
has more than one key. The skill sets `found = len(item) > 1` for each
record and exposes it as a column.

```
{"CID": 2244, "MolecularFormula": "C9H8O4", ...}   → found = True
{"CID": 999999999}                                   → found = False
```

### Trap 3 — `RegistryID` is not a list of CAS numbers

`xrefs/RegistryID` returns cross-references to *every* registry PubChem knows
about. For CID 2244 (aspirin) that is **461 entries**, of which only a couple
happen to look like CAS numbers:

```
CHEBI:15365
DB00945
LGCAMP0133.00-11
50-78-2              ← this one is a CAS number
...
```

The CAS registry number is *among* them but not *identified* among them. The
same compound's `synonyms` endpoint returns **699 entries**, among which
**three** match the CAS pattern — and those three are the meaningful ones
(aspirin's own CAS plus the CAS of related salt/hydrate forms).

### CAS extraction

CAS registry numbers have a fixed format:

```
NNNNN-NN-N
└─2-7─┘ └2┘ └1┘
```

- 2–7 leading digits
- exactly 2 digits
- exactly 1 digit

The last digit is a check digit computed from the rest, so the pattern alone
catches every CAS number and rejects essentially everything else — `CHEBI:15365`
has colons, `DB00945` has no dashes, `LGCAMP0133.00-11` has letters and two
dashes but the wrong shape.

The extraction is a regex filter over the synonym list:

```python
CAS_PATTERN = re.compile(r"^\d{2,7}-\d{2}-\d$")
```

**It is a format check, not a registry lookup.** The check digit can be
verified for extra confidence; the skill does not do so because a small
fraction of CAS numbers in the wild carry transcription errors in PubChem
itself, and rejecting them would silently drop legitimate entries. For
identity-critical work, cross-check the InChIKey instead.

### Property set

The default set was chosen to cover the descriptors that show up in compound
tables and QSAR panels:

| Group | Properties |
|---|---|
| Composition | `MolecularFormula`, `MolecularWeight`, `ExactMass` |
| Structure | `SMILES`, `ConnectivitySMILES`, `InChI`, `InChIKey` |
| Name | `IUPACName` |
| Physicochemical | `XLogP`, `TPSA`, `Charge` |
| Counts | `HBondDonorCount`, `HBondAcceptorCount`, `RotatableBondCount`, `HeavyAtomCount` |

Note that `MolecularWeight` (average, monoisotopic-averaged) and `ExactMass`
(monoisotopic) are different numbers and are frequently confused. They agree
to the first two or three decimals for light elements and diverge as
molecular weight increases.

## 2. NCBI Taxonomy

### The protocol

Taxonomy records are fetched directly with `efetch`:

```
GET /entrez/eutils/efetch.fcgi?db=taxonomy&id=9606,10090,4932&retmode=xml&email=...
```

**No `esearch` step is needed.** Unlike BioProject — where `PRJNA` accessions
must be translated into internal UIDs — a Taxonomy ID *is* the primary key.
This makes the flow one request per chunk instead of two.

### Trap 4 — `LineageEx` holds only ancestors

```xml
<Taxon>
  <TaxId>9606</TaxId>
  <ScientificName>Homo sapiens</ScientificName>
  <Rank>species</Rank>
  <LineageEx>
    <Taxon><TaxId>131567</TaxId><Rank>cellular root</Rank>...</Taxon>
    <Taxon><TaxId>2759</TaxId><Rank>domain</Rank>...</Taxon>
    ...
    <Taxon><TaxId>9605</TaxId><Rank>genus</Rank><ScientificName>Homo</ScientificName></Taxon>
    <!-- note: no species entry here -->
  </LineageEx>
</Taxon>
```

The `<LineageEx>` block stops at `genus`. The species level — usually the
level you actually want — exists **only** at the top of the record. A parser
that reads only `LineageEx` will produce a lineage that ends at the genus and
a species column that is empty.

The fix is to add the top-level `<Rank>`/`<ScientificName>` into the
rank-keyed map, which is what `_parse_taxon` does.

### Trap 5 — the ranks go well beyond the textbook eight

`LineageEx` contains far more than the standard
superkingdom/kingdom/phylum/class/order/family/genus/species. A real
human lineage contains 30 entries with ranks including:

```
cellular root · domain · kingdom · phylum · class · order · family · genus
· clade · subphylum · superclass · subclass · suborder · infraorder
· parvorder · superfamily · subfamily · tribe · subgenus
```

Selecting only the eight textbook ranks discards `domain` — which for
*Eukaryota* vs *Bacteria* vs *Archaea* is exactly the distinction most
queries want.

`superkingdom` deserves a specific note: **NCBI retired that rank**. In
current data "cellular organisms" has rank `cellular root` and the top-level
split uses `domain`. A `superkingdom` column is therefore legitimately empty
against the live API — it is not a parsing bug. The column is kept only for
compatibility with older data snapshots.

### Trap 6 — `root.iter("Taxon")` over-collects

This is the subtlest of the set, because it produces a *plausible-looking*
result rather than an obvious error.

```python
root.iter("Taxon")        # 75 nodes for a 3-taxon request
root.findall("Taxon")     # 3 nodes — the actual requested taxa
```

`iter()` is recursive, so it descends into each record's `LineageEx` and
collects the ancestors as if they were separate taxa. The result is 3
requested species plus 72 "records" named `cellular organisms`,
`Eukaryota`, `Opisthokonta`, `Metazoa` and so on. The output looks like a
complete taxonomic dump; it is simply 24× too large and full of duplicates
of lineage entries.

**Always iterate the root's direct children**: `root.findall("Taxon")`.

### Lineage representation

Each lineage entry keeps three fields — `tax_id`, `rank`, `name`:

```json
{
  "tax_id": "9605",
  "rank": "genus",
  "name": "Homo"
}
```

The taxid matters more than the name. NCBI renames clades over time — the
bacterial phylum formerly called `Proteobacteria` is now `Pseudomonadota`
— while the taxid stays fixed. Any join keyed on scientific name will break
when NCBI renames a clade.

The CSV flattens the common ranks into columns and keeps the entire chain in
`full_lineage` as a ` > `-separated string; the JSON keeps every entry
structurally, so ranks beyond the flattened ones (`clade`, `superfamily`,
`tribe`, …) are not lost.

## 3. Rate limits and etiquette

| Service | Limit | With API key |
|---|---|---|
| PubChem PUG-REST | 5 req/s, 400 req/min | n/a |
| NCBI E-utilities | 3 req/s | 10 req/s |

The clients here throttle to 0.25 s between PubChem requests and 0.34 s
between NCBI requests (3 req/s), with exponential backoff on failure.

NCBI's usage policy asks for a contact `email` on every request so that
problems can be reported to a real person. The skill makes the email
**mandatory** for the Taxonomy path and reads it from `--email` or the
`BIOPIPE_EMAIL` environment variable.

Fabricating placeholder addresses is not acceptable here. It defeats the
purpose of the field — nobody can be contacted — and it is a violation of
NCBI's terms of use. The pipeline this skill came from prepared five
`userN@example.com` strings and picked one at random per request; that
pattern was removed.

## 4. Why not a browser?

The original scripts drove a Selenium browser to reach these same pages.
For a data portal that publishes a documented REST API, that is strictly
worse:

| | Selenium scraping | REST API |
|---|---|---|
| Throughput | ~1 record / 2–5 s | 100 records / request |
| Fragility | breaks on any template change | versioned, documented contract |
| Reproducibility | depends on browser state | depends only on the request |
| Dependencies | driver + browser binary | `requests` |
| Rate limiting | looks like a user, hits limits fast | uses the published quota |

The exception in the original collection was the docking servers, which have
no public API — but that is a different class of tool and belongs in a
different skill.

## References

- PubChem PUG-REST: <https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest>
- PubChem property table: <https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest/property>
- NCBI Taxonomy: <https://www.ncbi.nlm.nih.gov/taxonomy>
- NCBI E-utilities: <https://www.ncbi.nlm.nih.gov/books/NBK25501/>
- NCBI Taxonomy XML DTD: <https://www.ncbi.nlm.nih.gov/entrez/query/DTD/taxon.dtd>
