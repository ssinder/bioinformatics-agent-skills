# Methodology: BioProject Retrieval

## 1. The E-utilities protocol

NCBI exposes a single API surface, E-utilities, for all its databases. The
BioProject-relevant endpoints are three:

| Endpoint | Purpose | Input | Output |
|---|---|---|---|
| `esearch.fcgi` | text query → UIDs | `term` | JSON with `idlist` |
| `efetch.fcgi` | UIDs → full records | `id` | XML `RecordSet` |
| `esummary.fcgi` | UIDs → brief summaries | `id` | XML `DocumentSummarySet` |

`efetch` and `esummary` accept UIDs, not accessions, so `esearch` is always
required first. (Skipping esearch and passing accessions to efetch does not
work.)

### The three-step flow

```
accessions → esearch.fcgi (term=PRJNA172563[Project Accession] OR ...)
           → uid list
           ├→ esummary.fcgi (id=172563,212358, retmode=xml) → structured fields
           └→ efetch.fcgi   (id=172563,212358, rettype=xml) → full description
           → merged records
```

### Why both esummary and efetch

**This is the single most important design decision in the skill, and getting
it wrong is silent.** The two endpoints carry *disjoint* subsets of the
metadata:

| Field | `esummary` | `efetch` |
|---|---|---|
| accession, title, name | ✅ | ✅ |
| description | ✅ (may be truncated) | ✅ (full text) |
| **organism name** | ✅ | ⚠️ only via `<OrganismName>` |
| **taxonomy ID** | ✅ | ✅ (as `taxID` attribute) |
| **project type** | ✅ | ⚠️ only the method |
| **submitter organisation** | ✅ | ✅ (nested deeper) |
| registration date | ✅ (slash format) | ✅ (as `submitted` attribute) |
| **release date** | ❌ | ✅ |

Calling only `efetch` loses the organism, project type and submitter for the
records where they are absent from the descriptive XML — and nothing raises
an error; the CSV just has empty cells. Calling only `esummary` loses the
release date and truncates some descriptions.

The merge rule is **efetch first, esummary fills the gaps**: `efetch` has the
fuller description and ISO dates, so it wins whenever both have a value.

### The `[Project Accession]` field qualifier

Without a field qualifier, `PRJNA172563` is searched as a free-text term and
can match records where the string merely appears somewhere in the text. The
qualifier restricts the match to the accession field:

```
PRJNA172563[Project Accession]
```

Multiple accessions are combined with `OR`. The resulting term is what the
old scraping code was trying to approximate by regex-matching `ID:\s*(\d+)`
in page text — but structured, and it handles `PRJNA`/`PRJEB`/`PRJDB`
prefixes uniformly.

### Common query fields for BioProject

| Field tag | Matches |
|---|---|
| `[Project Accession]` | the PRJ accession itself |
| `[Organism]` | the organism name |
| `[Title]` | project title |
| `[Description]` | project description text |
| `[Submitter]` | submitting institution |
| `[Project Data Type]` | e.g. "BioProject" |

For a keyword sweep across *all* projects rather than a known list:

```
"CRISPR"[Description] AND "Homo sapiens"[Organism]
```

## 2. Response structure

### `efetch` (real shape)

```xml
<RecordSet><DocumentSummary uid="355343">
    <Project>
        <ProjectID>
            <ArchiveID accession="PRJNA355343" archive="NCBI" id="355343"/>
        </ProjectID>
        <ProjectDescr>
            <Name>Salmonella enterica subsp. enterica serovar Typhi strain:429038</Name>
            <Title>Salmonella enterica genome sequencing</Title>
            <Description>Cephalosporin resistant Salmonella Typhi</Description>
            <ProjectReleaseDate>2016-12-15T00:00:00Z</ProjectReleaseDate>
        </ProjectDescr>
        <ProjectType>
            <ProjectTypeSubmission>
                <Target capture="eOther" material="eOther" sample_scope="eMultiisolate">
                    <Organism species="28901" taxID="90370">
                        <OrganismName>Salmonella enterica subsp. enterica serovar Typhi</OrganismName>
                        <Supergroup>eEukaryotes</Supergroup>
                    </Organism>
                </Target>
                <Method method_type="eOther"/>
            </ProjectTypeSubmission>
        </ProjectType>
    </Project>
    <Submission last_update="2016-11-29" submitted="2016-11-29" submission_id="SUB365098">
        <Description>
            <Organization role="owner" type="institute">
                <Name>Christian Medical College</Name>
            </Organization>
            <Access>public</Access>
        </Description>
    </Submission>
</DocumentSummary></RecordSet>
```

Field locations used by the parser:

| Field | Path |
|---|---|
| accession | `Project/ProjectID/ArchiveID/@accession` |
| title | `Project/ProjectDescr/Title` |
| name | `Project/ProjectDescr/Name` |
| description | `Project/ProjectDescr/Description` |
| release date | `Project/ProjectDescr/ProjectReleaseDate` |
| organism | `Project/ProjectType/.../Target/Organism/OrganismName` |
| taxonomy ID | `.../Organism/@taxID` (or `@TaxonomyID`) |
| submitter | `Submission/Description/Organization/Name` |
| registration date | `Submission/@submitted` |
| last update | `Submission/@last_update` |

### Three traps in that XML

**1. `species` is a taxid, not a name.** Look again at
`<Organism species="28901" taxID="90370">`. The `species` attribute holds
the *species-level* taxid (28901 = *Salmonella enterica*), while `taxID`
holds the strain-level one (90370). The human-readable name is in the
`<OrganismName>` child element. Reading `species` as the organism name
produces `organism=28901` — a number where a name belongs, with no error
raised.

**2. The dates are attributes, not elements.** `<Submission submitted="2013-07-16">`
and `<Submission last_update="2013-10-23">` carry the dates as XML
*attributes*. A parser that only looks at child elements finds nothing and
reports an empty date for every record.

**3. `ProjectReleaseDate` is optional.** Only some records carry it — the
release date is set when the project was actually released, so embargoed and
old records frequently lack it. Its absence is normal, not an error.

### `esummary` (real shape)

```xml
<eSummaryResult>
  <DocumentSummarySet status="OK">
    <DocumentSummary uid="355343">
      <TaxId>90370</TaxId>
      <Project_Acc>PRJNA355343</Project_Acc>
      <Project_Type>Primary submission</Project_Type>
      <Registration_Date>2016/11/29 00:00</Registration_Date>
      <Project_Name>Salmonella enterica subsp. enterica serovar Typhi</Project_Name>
      <Project_Title>Salmonella enterica genome sequencing</Project_Title>
      <Project_Description>Cephalosporin resistant Salmonella Typhi</Project_Description>
      <Organism_Name>Salmonella enterica subsp. enterica serovar Typhi</Organism_Name>
      <Submitter_Organization>Christian Medical College</Submitter_Organization>
      <Submitter_Organization_List><string>Christian Medical College</string></Submitter_Organization_List>
    </DocumentSummary>
  </DocumentSummarySet>
</eSummaryResult>
```

Field names here use underscores (`Project_Title`, not `Title`) and the
organism name lives in `Organism_Name` — **not** in a `species` attribute.
`Submitter_Organization` holds a single organisation;
`Submitter_Organization_List` holds a `<string>` list when there are
several.

### Date normalisation

Three endpoints, three date formats — all normalised to `YYYY-MM-DD`:

| Source | Raw format |
|---|---|
| `efetch` `ProjectReleaseDate` | `2016-12-15T00:00:00Z` |
| `efetch` `Submission/@submitted` | `2013-07-16` |
| `esummary` `Registration_Date` | `2013/07/16 00:00` |

The `esummary` date is occasionally the placeholder **`1/01/01 00:00`**,
which NCBI uses for "date unknown" (common for dbGaP-backed projects whose
submitter details are stripped). It is deliberately *not* converted into a
real date — the field stays empty rather than showing a bogus `0001-01-01`.

**The description may be split across multiple `<Description>` elements.**
The parser joins them with a space when the first one is empty, so a
description that is only present in a secondary element is not lost.

**A `<Description>` may contain embedded HTML** (submitters often write
`<i>`, `<b>` or links). The parser keeps the raw text; strip markup
downstream if the output must be plain text.

## 3. The email requirement

NCBI's E-utilities policy states that all users should identify themselves
with a valid email address, and NCBI reserves the right to block callers
that do not. With an email (and optionally an API key), the rate limit is:

| Credentials | Rate limit |
|---|---|
| none | 3 requests/second |
| email only | 3 requests/second |
| email + API key | 10 requests/second |

Consequences for a tool that fakes the address:

- The limit applies per IP regardless. A random `userN@example.com` gains
  nothing over a real address, and it makes it impossible for NCBI to reach
  you when something goes wrong.
- Fabricated addresses in a public repository are a privacy and professionalism
  problem: they invite other users to reuse them, and they get blocked.

This skill therefore makes `--email` mandatory and reads `BIOPIPE_EMAIL` from
the environment as a convenience. It never invents an address.

## 4. Why not scrape the web page?

For comparison, the scraping approach needs, for a single project:

1. An HTTP request with a browser-like User-Agent (NCBI blocks obvious bots).
2. A selector for the description container — the position of which changes
   when the page layout changes.
3. A separate selector for the organism, the title and the dates.
4. HTML-entity decoding and whitespace normalisation.
5. A retry policy (NCBI throttles aggressively).

The API needs no selectors at all. The XML field names are part of NCBI's
documented schema; the CSS structure is not.

A second, subtler cost of scraping: it cannot distinguish "the field is
empty" from "the selector broke". With the API, a missing element means the
submitter did not supply the field.

A third cost applies to the API as well, and is worth understanding: the
*structured* fields and the *descriptive* fields live in different places
(see §2). A hand-written parser that assumes one endpoint carries everything
will produce a CSV with empty organism and submitter columns — a failure
that produces no error message, only quietly missing data.

## 5. Handling multi-project queries

E-utilities accepts up to ~500 UIDs per request. This skill issues one
`esearch` for the whole accession list, then one `esummary` and one `efetch`
for all resulting UIDs — **3 requests** for any realistic batch, well inside
the rate limit (3 req/s without an API key, 10 req/s with one).

For very large lists (thousands of accessions), split the input and merge
the results, or use the `retstart`/`retmax` pagination parameters on
`esearch` for keyword queries.

## 6. Output format

Results are written as both CSV and JSON:

- **CSV** for spreadsheets and quick inspection. The description column may
  contain newlines — the file is written with `newline=""` and the standard
  `csv` module, so embedded newlines are quoted correctly rather than
  corrupting the row structure.
- **JSON** for programmatic consumption, where the description keeps its
  exact text without CSV quoting.

The old implementation rewrote the entire CSV after every single record,
which is O(n²) in disk writes. Collecting first and writing once is both
faster and produces a consistent file (no partially-written state if the
process is interrupted mid-run).

## References

- NCBI E-utilities usage policy: <https://www.ncbi.nlm.nih.gov/books/NBK25497/>
- E-utilities documentation: <https://www.ncbi.nlm.nih.gov/books/NBK25501/>
- BioProject field list: <https://www.ncbi.nlm.nih.gov/bioproject/label/BioProject%20Data%20Fields>
- ENA's equivalent portal (for PRJEB projects): <https://www.ebi.ac.uk/ena/portal/api>
