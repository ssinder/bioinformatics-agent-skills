# Methodology: HDOCK Docking via HTTP

## 1. Why plain HTTP works

HDOCK has no public REST API. The interface is a PHP web form. But a web
form *is* an HTTP interface, and this particular one is unusually simple:

- **No CSRF token.** No `__VIEWSTATE`, no `__EVENTVALIDATION` — this is a
  plain PHP script, not ASP.NET.
- **No hidden fields.** Every input is a visible form control.
- **No JavaScript signature.** The `OnSubmitForm()` handler only validates
  inputs client-side and rewrites the form action.
- **No session cookie.** The server does not set one.

The consequence: a `requests` POST reproduces exactly what the browser does,
minus the browser.

### The form

```html
<form id="form1" name="form1" enctype="multipart/form-data"
      onSubmit="OnSubmitForm();" action="" method="post">
```

`enctype="multipart/form-data"` means the file inputs and the plain fields
travel in one multipart body. That is what `requests` does when you pass
both `data=` and `files=`.

### The action rewrite — the one real gotcha

The `action` attribute is empty. The actual target is chosen by JavaScript:

```javascript
function OnSubmitForm(){
    let sym = document.getElementById("symmetry");
    ...
    if( ! sym.value ){
        ...
        document.form1.action = "runhdock.php";     // ordinary docking
    }else{
        ...
        document.form1.action = "runhsymdock.php";  // symmetric multimer
    }
    return true;
}
```

**Posting to `http://hdock.phys.hust.edu.cn/` returns the homepage.**
Byte for byte. HTTP 200, the same 16,787-byte document, no error, no
redirect, no job. A script that checks the status code will happily report
"submitted" while nothing happened at all.

The correct endpoints are:

| Case | Endpoint |
|---|---|
| Ordinary docking | `http://hdock.phys.hust.edu.cn/rundock.php` → **`runhdock.php`** |
| Symmetric multimer | `http://hdock.phys.hust.edu.cn/runhsymdock.php` |

(`runhdock.php` — the `h` is part of the name. `rundock.php` is a 404.)

## 2. The form fields

All fields travel in one multipart body. Empty values are harmless; the
server only checks that exactly one input per molecule is non-empty.

### Receptor — exactly one of

| Field | Type | Meaning |
|---|---|---|
| `pdbfile1` | file | PDB file upload |
| `pdbid1` | text | `PDBID:ChainID`, e.g. `1CGI:E` |
| `fastaseq1` | text | FASTA sequence pasted inline |
| `fastafile1` | file | FASTA file upload |

### Ligand — exactly one of

| Field | Type | Meaning |
|---|---|---|
| `pdbfile2` | file | PDB file upload |
| `pdbid2` | text | `PDBID:ChainID` |
| `fastaseq2` | text | FASTA sequence |
| `fastafile2` | file | FASTA file upload |

### Shared options

| Field | Type | Meaning |
|---|---|---|
| `email` | text | Notification address; the server emails you when the job finishes |
| `jobname` | text | Job label (the server replaces whitespace with underscores) |
| `ligtyp` | select | `protein`, `rna` (ssRNA), `dna` (dsDNA), `rna.d` (dsRNA), `ssdna` (ssDNA) |
| `symmetry` | text | `C2`/`C3` (cyclic) or `D2`/`D3` (dihedral); **changes the endpoint** |
| `sitenum1` / `sitenum2` | text | Binding-site residue numbers, comma separated |
| `pdbsite` | text | Binding-site specification |
| `restrnum` | text | Residue distance restraints |
| `sitefile1` / `sitefile2` | file | Binding-site files |
| `restrfile` | file | Restraint file |
| `saxsfile` | file | SAXS experimental data |
| `upload` | submit | Submit button |

### Client-side validation, reproduced server-side

`OnSubmitForm()` refuses to submit when the receptor or ligand input is
empty, and sanitises the job name:

```javascript
if ( inputjobname.value ) {
    inputjobname.value = inputjobname.value.trim().replace(/\s+/g, '_');
}
```

The skill applies the same job-name normalisation locally, so the name in
`jobs.csv` matches what the server will show.

## 3. The result page state machine

The job URL has the shape `http://hdock.phys.hust.edu.cn/data/<job_id>/`,
where `job_id` is a 12-character hex-ish token.

The same URL serves three different pages:

```
GET /data/<job_id>/
   │
   ├── contains "is QUEUED" ────────────► QUEUED
   │      "Your HDOCK job <name> is QUEUED."
   │      "This page is automatically refreshed every 10 seconds."
   │
   ├── contains "ERROR:" ───────────────► FAILED
   │      "ERROR: No structure found for the Receptor molecule."
   │
   └── contains a "Docking Score" table ► DONE
```

The queued page says the server refreshes every 10 seconds, so polling
faster than that is pointless. The skill polls every 20 seconds by
default.

**Distinguishing QUEUED from FAILED is not cosmetic.** The old Selenium
script extracted the scores with
`//th[text()='Docking Score']/following-sibling::td` and, when the element
was absent, wrote `N/A` and **aborted the entire batch**. One unfinished
job at position 1 of 50 meant positions 2–50 were never checked.

### The finished page is a ranking matrix

```html
<tr><th>Rank</th><th>1</th><th>2</th>...<th>10</th></tr>
<tr><th>Docking Score</th><td>-444.95</td><td>-244.21</td>...</tr>
<tr><th>Ligand rmsd (&Aring;)</th><td>0.42</td><td>18.23</td>...</tr>
```

Rows are metrics, columns are ranks. Parsing this as "label → first value"
loses information: the full row tells you how sharply the score drops off
between rank 1 and rank 10, which is often more informative than the best
score alone. The skill keeps the whole row in `rankings`.

Note the `&Aring;` entity in "Ligand rmsd (&Aring;)" — the label is
*Ligand rmsd (Å)*, HTML-escaped.

## 4. The scoring functions

### Docking Score

Computed by a knowledge-based iterative scoring function (ITScorePP for
protein-protein, ITScorePR for protein-RNA). **More negative is better.**
For reference, protein complexes in the PDB normally score around −200 or
better.

The documentation is careful about what this means:

> the score should not be treated as the true binding affinity of two
> molecules because it has not been calibrated to the experimental data.

It is a **ranking** score, not an affinity prediction.

### Confidence Score

HDOCK derives it from the Docking Score with an explicitly published
closed form:

```
Confidence_score = 1.0 / [1.0 + e^(0.02 * (Docking_Score + 150))]
```

For the example above (Docking Score = −444.95):

```
= 1.0 / (1.0 + e^(0.02 * (-444.95 + 150)))
= 1.0 / (1.0 + e^(-5.899))
= 1.0 / (1.0 + 0.00274)
= 0.9973
```

The documentation states that "when the confidence score is above 0.7, the
two molecules would most likely interact".

Two consequences:

1. **Confidence Score is a monotone function of Docking Score.** It carries
   no information the score does not already contain; it only rescales it
   onto an interpretable 0–1 scale.
2. **It can always be computed.** When the results page does not carry the
   column, computing it from the Docking Score is exact — strictly better
   than waiting for an element that may never appear, which is what the old
   Selenium script did.

The skill uses the page value when present and falls back to the formula.

## 5. Output files

The full docking output package contains a file named like
`hdock_5c984053e4b83.out` — all 4,392 sampled binding modes. Its format:

```
Grid spacing: 1.200
Angle step: 15.000
Initial rotation: 0.00000 0.00000 0.00000
1CGI_r_b.pdb 23.562 26.523 22.675      <- receptor file + centre of geometry
1CGI_l_b.pdb 47.776 34.961 33.826      <- ligand file + centre of geometry
1.27246 0.01055 5.02167 -0.328 -0.164 0.264 -445.20 0.45 1.00 2.80075 ...
```

From the 6th line onward each row is one predicted binding mode: three
translations, three rotations, the binding score, the RMSD from the initial
ligand orientation, and the translational ID for the rotation.

The server pre-generates the **top 100** binding models. Users wanting more
than 100 can download `createpl_linux` and run it locally on the `.out`
file.

By default HDOCK samples 4,392 rotations (15° Euler angle step) and keeps
the best 100 models after RMSD-based clustering with a 5 Å cutoff.

## 6. Rate limiting and etiquette

HDOCK is a free academic service with no published rate limit. That makes
courtesy the only constraint, and the old scripts were not courteous: they
launched one browser instance per file with `threading.Thread` and no
delay between submissions.

This skill:

- submits **serially**, one job per command invocation;
- sleeps at least 1 second between requests within a command;
- keeps the number of concurrent jobs under the user's control.

For a large campaign, submit in batches with a pause, and collect results
later with `--jobs-file` rather than holding a long poll open.

## 7. Why not Selenium

| | Selenium | HTTP |
|---|---|---|
| Submission latency | ~3–8 s (browser startup) | ~1–2 s |
| Dependencies | browser binary + matching driver | `requests` |
| Portability | driver version must match the browser exactly | none |
| Fragility | XPath selectors break silently on a redesign | explicit field names |
| Resource use | ~200 MB per browser instance | negligible |
| Batch throughput | serialised by browser startup | trivially parallel |

The decisive point is that the old scripts hard-coded
`C:\Users\DELL\Desktop\protein\msedgedriver.exe`. That path exists on one
machine. The HTTP approach works everywhere.

## References

- HDOCK server: <http://hdock.phys.hust.edu.cn/>
- HDOCK help: <http://hdock.phys.hust.edu.cn/help.php>
- Yan Y. et al. Nature Protocols 15, 2020 (server protocol).
  <https://doi.org/10.1038/s41596-020-0312-x>
- Yan Y. et al. Nucleic Acids Res 45(W1):W365, 2017 (algorithm).
  <https://doi.org/10.1093/nar/gkx357>
- Huang S-Y, Zou X. Nucleic Acids Res 42:e55, 2014 (ITScorePR).
