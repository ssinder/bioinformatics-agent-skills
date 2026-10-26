# Troubleshooting: HDOCK Docking

## Submission

### The command returns "submitted" but nothing ever appears

Almost certainly the form was posted to the homepage instead of
`runhdock.php`. The homepage's `<form>` has an **empty** `action`
attribute; the real target is set by JavaScript at submit time.

The symptom is distinctive: the response is HTTP 200 with a body that is
**byte-identical to the homepage** (16,787 bytes). No error, no redirect,
no job.

```bash
# wrong — silently returns the homepage
curl -s -X POST -d "fastaseq1=..." http://hdock.phys.hust.edu.cn/

# right
curl -s -X POST -d "fastaseq1=..." http://hdock.phys.hust.edu.cn/runhdock.php
```

The skill posts to `runhdock.php` (`runhsymdock.php` when `--symmetry` is
given), so this should not occur — but if you adapt the code, note that
the spelling is `runhdock.php`. **`rundock.php` is a 404.**

### HTTP 404 on submit

Same cause: a misspelled endpoint. The two valid paths are `runhdock.php`
and `runhsymdock.php`.

### The job never appears on the server

Check that the response URL contains `/data/<job_id>/`. If the URL is still
the homepage, the submission was not accepted — see the first item above.

### "Receptor input is required" / "Ligand input is required"

The server found no usable input for one of the molecules. Remember:

- each molecule needs **exactly one** of the four input forms;
- the field names are `pdbfile1`/`pdbid1`/`fastaseq1`/`fastafile1` for the
  receptor and the `2` variants for the ligand;
- an **empty string does not count** as input.

If you pass a FASTA *file path* to `--receptor-seq`, the skill reads the
file and sends its content as `fastaseq1` — that is correct. But if you
pass a path that does not exist, the skill treats the path string itself as
the sequence, and the server rejects it.

### The job fails with "No structure found for the Receptor molecule"

The server could not build a structure from the receptor input. Usual
causes:

1. **Wrong chain ID.** The format is `PDBID:ChainID` — a four-character
   PDB ID, a colon, then the chain. `1CGI:E` is valid; `1CGI` alone is
   not; `1CGI:E:A` is not.
2. **Non-standard residues in the sequence.** Only the 20 standard amino
   acids are accepted. Remove `X`, `B`, `Z`, `J` and any non-letter
   characters.
3. **A PDB file the server cannot parse.** Some tools emit PDB files with
   unusual line endings or missing headers. Re-export or use `PDBID:ChainID`.
4. **A sequence too long for the text field.** The `<textarea>` declares
   `maxlength="5500"`. Use `fastafile1` for long sequences.

The error is asymmetric — the receptor error appears first, and the ligand
error appears only if the receptor was fine. Fix them one at a time.

### DNA/RNA sequence input fails

HDOCK's own notice (July 2026) records a period when sequence input for
DNA/RNA was broken; it was fixed, but nucleotide input remains the least
exercised path. For nucleic acids, upload PDB structures instead.

## Retrieval

### The result page says "is QUEUED"

The job is still being computed server-side. This is **not** an error.
HDOCK's documentation says a job "can normally complete within 30 min";
small peptide pairs finish faster, large complexes take much longer.

- Re-run `fetch_docking_results.py` later.
- Do not resubmit — a resubmission starts a fresh job and wastes server
  time.
- `--wait` will poll for you, but note that the command then holds a
  connection open for up to `--max-wait` seconds.

### The Docking Score column is empty

Three different situations produce this, and they need different responses:

| What the page shows | Meaning | Action |
|---|---|---|
| `is QUEUED` | Job not finished | Wait and retry |
| `ERROR: ...` | Job failed | Fix the input and resubmit |
| Table present but no `Docking Score` row | Unrecognised page layout | Check `references/methodology.md` for the expected structure |

The skill reports these as three distinct statuses (`queued`, `failed`,
`unknown`) rather than collapsing them, precisely because the old script
collapsed "not ready" into "failed" and aborted the batch.

### The Confidence Score is 0 or missing

If the page provides the column, it is used directly. If not, the skill
computes it from the Docking Score using HDOCK's published formula:

```
Confidence_score = 1.0 / [1.0 + e^(0.02 * (Docking_Score + 150))]
```

If it is still missing, the Docking Score itself is missing — check the
status first.

### The model download finds nothing

`--download-models` only works on a **finished** job. On a queued or failed
job it is skipped with a warning.

Also note that the results page exposes several kinds of resources: the
ranked complex models (`top_N.pdb`), the full docking output
(`hdock_<job_id>.out`), and assorted text files. The download step filters
for `.pdb` files and takes the top `--top-n` by the rank number in the
filename.

### The job URL returns 404

Results are stored on the server for **two weeks**. Beyond that the URL is
dead and the job must be resubmitted. The `job_id` in the URL is the
stable identifier — record it in `jobs.csv` at submission time, because
the job *name* is normalised by the server and is not guaranteed unique.

## Score interpretation

### A score of −444 is excellent, −150 is mediocre

More negative is better. Protein complexes in the PDB normally score
around −200 or better, so:

| Docking Score | Rough reading |
|---|---|
| < −400 | Strong; confidence > 0.99 |
| −300 to −400 | Good; confidence > 0.9 |
| −200 to −300 | Typical for a real complex |
| −100 to −200 | Weak; confidence < 0.7 |
| > −100 | Very weak |

Confidence > 0.7 is the documented threshold for "the two molecules would
most likely interact".

### The score is not a binding affinity

HDOCK's documentation is explicit that the score "should not be treated as
the true binding affinity of two molecules because it has not been
calibrated to the experimental data". Use it to rank poses and candidates;
do not convert it to a Kd or report it as one.

### Rank 1 is much better than rank 2 — what does that mean?

A steep drop-off between rank 1 and rank 2 usually means the first pose is
a clearly preferred binding mode. A flat curve means the algorithm found
no strongly preferred mode, and the ranking should be treated with more
caution. The full top-10 rows are available in the `rankings` field.

## Network problems

### Connection timeouts

The default timeout is 180 s for submission and 60 s for retrieval.
Behind a slow proxy, raise them:

```bash
python scripts/submit_docking.py ... --timeout 300
```

If the machine needs a proxy:

```bash
export http_proxy=http://127.0.0.1:7897
export https_proxy=http://127.0.0.1:7897
```

### The server is unreachable

`hdock.phys.hust.edu.cn` is a single academic host and does occasionally
go down for maintenance. There is no mirror. The skill logs the
underlying exception and returns a failed result rather than raising — the
job record is still written, so the submission can be retried later.

### A 429 or the server refuses connections

HDOCK does not publish a rate limit, but it is a shared academic
resource. Submit serially, not in parallel. The old scripts launched one
browser per file with no delay; that pattern gets a user IP throttled and
slows down everyone else's jobs.

## Still stuck

1. Check the job state with the retrieval command — it distinguishes
   queued, failed and unknown.
2. Compare against the server's own example page:
   <http://hdock.phys.hust.edu.cn/data/example/> — it shows the finished
   page layout with real numbers.
3. Re-read `references/methodology.md` §3 for the state machine and §4
   for the scoring functions.
4. If you are adapting the code, check the form field names against the
   live page — `curl -s http://hdock.phys.hust.edu.cn/ | grep -oE
   'name="[^"]+"' | sort -u`.
