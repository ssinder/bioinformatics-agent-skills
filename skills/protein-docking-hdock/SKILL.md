---
name: protein-docking-hdock
description: >-
  Submits protein-protein docking jobs to the HDOCK server and retrieves the
  resulting Docking Score, Confidence Score, ligand RMSD and ranked complex
  models. Uses plain HTTP form submission against the HDOCK server rather
  than browser automation, supports PDB files, PDB ID:ChainID references and
  FASTA sequences for both receptor and ligand, records submitted jobs to a
  CSV so results can be collected later, and can poll for completion. This
  skill should be used when running protein-protein complex structure
  prediction on the HDOCK server, checking the docking score of a submitted
  complex, or downloading top-ranked docked models is requested.
license: MIT
compatibility: Python 3.9+; needs requests for network access; standard library only otherwise
metadata:
  version: "1.0"
  last-reviewed: "2026-10-26"
  skill-author: sinder
---

# HDOCK Protein Docking

## Overview

This skill submits protein-protein docking jobs to the
**HDOCK server** (<http://hdock.phys.hust.edu.cn/>) and collects the results:
Docking Score, Confidence Score, ligand RMSD, and the ranked complex models.

### It uses HTTP, not a browser

HDOCK publishes no REST API — only a web interface. The natural assumption
is that this therefore requires Selenium. It does not.

The submission page is a **standard HTML form POST**. There is no CSRF
token, no hidden field, no JavaScript-computed signature, and no session
cookie. The form posts to `runhdock.php`, and the results page is an
ordinary HTML table. Everything the old Selenium scripts did is reachable
with `requests`.

The one thing worth knowing is that the form's `action` attribute is
**empty**, and the real target is set by JavaScript at submit time:

```javascript
if( ! sym.value ){
    document.form1.action = "runhdock.php";     // ordinary docking
}else{
    document.form1.action = "runhsymdock.php";  // symmetric multimer
}
```

Posting to the homepage instead of `runhdock.php` returns a **byte-identical
copy of the homepage** — HTTP 200, no error, no result. This is the single
most likely way to waste an hour of debugging.

### Submission and retrieval are separate steps

HDOCK's own documentation says a job "can normally complete within 30 min".
Baking the wait into the submission command means the user stares at a blank
terminal for half an hour. This skill therefore splits the workflow:

```
submit_docking.py     →  returns a job URL in ~2 seconds
        ↓ (minutes to tens of minutes, server-side)
fetch_docking_results.py →  reads the score
```

Submitted jobs are appended to `jobs.csv`, so a batch can be collected later
in one pass.

## When to Use This Skill

- Predicting a protein-protein complex structure with HDOCK
- Getting the Docking Score / Confidence Score for a candidate complex
- Screening many receptor-ligand pairs through the HDOCK server
- Downloading the top-ranked docked models for further analysis

## Quick Start

```bash
export PYTHONPATH="$PWD/packages"

# 用 PDB ID:ChainID 提交（最快，不用准备文件）
python scripts/submit_docking.py \
    --receptor-pdb-id 1CGI:E --ligand-pdb-id 1CGI:I \
    --job-name my_docking --output-dir results/

# 用 FASTA 序列提交
python scripts/submit_docking.py \
    --receptor-seq receptor.fasta --ligand-seq ligand.fasta \
    --output-dir results/

# 稍后取结果（查一次）
python scripts/fetch_docking_results.py \
    --url http://hdock.phys.hust.edu.cn/data/xxxxxxxxxxxx/ \
    --output-dir results/

# 批量收集所有已提交作业
python scripts/fetch_docking_results.py \
    --jobs-file results/jobs.csv --output-dir results/ --download-models
```

## Workflow

### Submission

1. **Provide one input per molecule.** Each of the receptor and the ligand
   accepts exactly one of four input forms: a PDB file, a
   `PDBID:ChainID` string, a FASTA sequence, or a FASTA file. The
   `--*-seq` arguments accept either inline sequence text or a path to a
   FASTA file — the skill detects which.
2. **Submit.** The skill GETs the homepage first (matching browser
   behaviour), then POSTs the multipart form to `runhdock.php` — or
   `runhsymdock.php` when `--symmetry` is given.
3. **Record.** The job name, ID and URL are appended to `jobs.csv` and
   `jobs.json` in the output directory. The command returns immediately.

### Retrieval

1. **Read the page.** The result page has three possible states, which the
   skill distinguishes explicitly:

   | State | Marker on the page | Reported status |
   |---|---|---|
   | Queued | `is QUEUED` | `queued` |
   | Finished | `Docking Score` table present | `done` |
   | Failed | `ERROR:` message | `failed` |

   This distinction matters. The old script treated "not ready yet" as
   "failed" and aborted the whole batch on the first unfinished job.
2. **Parse the ranking matrix.** The results page is a matrix: rows are
   metrics (Docking Score, Ligand rmsd), columns are ranks 1–10. The skill
   extracts the full row, so the top-10 scores are all available, not just
   the best one.
3. **Compute Confidence Score if needed.** HDOCK's documentation gives the
   closed form:

   ```
   Confidence_score = 1.0 / [1.0 + e^(0.02 * (Docking_Score + 150))]
   ```

   When the page does not include the column, the skill computes it from the
   Docking Score. The formula is public and exact — there is no reason to
   wait for an element that may never appear.
4. **Optionally download models.** `--download-models` fetches the
   top-ranked complex models.

## Parameters

### `submit_docking.py`

| Parameter | Default | Description |
|---|---|---|
| `--receptor-pdb-id` | — | Receptor as `PDBID:ChainID` (e.g. `1CGI:E`) |
| `--receptor-seq` | — | Receptor FASTA sequence, inline text or file path |
| `--receptor-pdb` | — | Receptor PDB file path |
| `--ligand-pdb-id` | — | Ligand as `PDBID:ChainID` |
| `--ligand-seq` | — | Ligand FASTA sequence, inline text or file path |
| `--ligand-pdb` | — | Ligand PDB file path |
| `--job-name` | timestamp | Job label |
| `--email` | none | Notification address; the server emails you when the job finishes |
| `--ligand-type` | `protein` | `protein`, `rna`, `dna`, `rna.d`, `ssdna` |
| `--symmetry` | none | Symmetry such as `C2` or `D3`; switches to the symmetric endpoint |
| `--receptor-site` | none | Receptor binding-site residues, comma separated |
| `--ligand-site` | none | Ligand binding-site residues, comma separated |
| `--restraint` | none | Residue distance restraints |
| `--output-dir` | `hdock_results` | Output directory for the job records |
| `--timeout` | `180` | Per-request timeout in seconds |

### `fetch_docking_results.py`

| Parameter | Default | Description |
|---|---|---|
| `--url` | — | Single job result URL |
| `--jobs-file` | — | Job record CSV to process in batch |
| `--output-dir` | `hdock_results` | Output directory for the results |
| `--wait` | off | Poll until the job finishes instead of checking once |
| `--poll-interval` | `20` | Polling interval in seconds |
| `--max-wait` | `3600` | Maximum wait in seconds |
| `--download-models` | off | Download the top-ranked complex models |
| `--top-n` | `5` | How many models to download |
| `--timeout` | `60` | Per-request timeout in seconds |

## Outputs

| File | Content |
|---|---|
| `jobs.csv` / `jobs.json` | Submitted jobs: name, ID, URL, status |
| `results.csv` / `results.json` | Results: status, Docking Score, Confidence Score, error |

## Validation

- **A Docking Score is not a binding affinity.** HDOCK's documentation is
  explicit: the score comes from a knowledge-based iterative scoring
  function and "should not be treated as the true binding affinity of two
  molecules because it has not been calibrated to the experimental data."
  It is useful for *ranking* candidate poses, not for predicting a Kd.
- **A more negative score is better.** For comparison, complexes in the PDB
  normally score around −200 or better.
- **Confidence Score above 0.7 means "likely to interact".** The score is a
  monotone function of the Docking Score, so it adds no information beyond
  the score itself — it just puts the result on an interpretable scale.
- **A queued job is not a failed job.** `status=queued` means the server is
  still computing. HDOCK jobs commonly take 5–30 minutes; re-run the
  retrieval command later rather than resubmitting.
- **"No structure found" means the input was unusable**, not that the
  interaction is impossible. Common causes: a chain ID that does not exist
  in that PDB entry, a sequence containing characters other than the 20
  standard amino acids, or a PDB file that the server cannot parse. Check
  the chain ID format (`1CGI:E` — four-character PDB ID, colon, chain).
- **Results are stored for two weeks** by the server. Retrieve them within
  that window or resubmit.
- **Sequence input is modelled, not docked as-is.** HDOCK builds a structure
  from the sequence first (via homology modelling), then docks it. Upload
  PDB structures when you have them — the accuracy is better and the job is
  faster.
- **Bind the score to the job, not to the URL alone.** Job names are
  normalised by the server (spaces become underscores); the `job_id` in the
  URL is the stable identifier.

## References

- `references/methodology.md` — the submission protocol, the form fields,
  the page state machine, the scoring functions
- `references/troubleshooting.md` — common failures and how to diagnose them
- HDOCK server: <http://hdock.phys.hust.edu.cn/>
- HDOCK help page: <http://hdock.phys.hust.edu.cn/help.php>
- Yan Y. et al. *The HDOCK server for integrated protein-protein docking.*
  Nature Protocols (2020). <https://doi.org/10.1038/s41596-020-0312-x>
- Yan Y. et al. *HDOCK: a web server for protein-protein and protein-DNA/RNA
  docking.* Nucleic Acids Research (2017).
  <https://doi.org/10.1093/nar/gkx357>
