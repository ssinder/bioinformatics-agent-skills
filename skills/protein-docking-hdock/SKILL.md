---
name: protein-docking-hdock
description: "Use when running protein-protein complex structure prediction on the HDOCK server, checking the docking score of a submitted complex, or downloading top-ranked docked models. Submits jobs by plain HTTP rather than browser automation, accepting PDB files, PDB ID:ChainID references and FASTA sequences; retrieves Docking Score, Confidence Score, ligand RMSD and ranked models, records submissions to CSV, and polls for completion. 用于 HDOCK 对接与结果获取。"
license: MIT
compatibility: Python 3.9+; needs requests for network access; standard library only otherwise
metadata:
  version: "1.0"
  last-reviewed: "2026-10-26"
  skill-author: sinder
  edam-topic: "topic_0176"
  edam-operation: "operation_3899"
  edam-input-data: "data_1460"
  edam-output-data: "data_1460, data_1537"
  edam-input-format: "format_1475, format_1929"
  edam-output-format: "format_1475"

# HDOCK Protein Docking


---

## Overview

Submits protein-protein docking jobs to the **HDOCK server** (<http://hdock.phys.hust.edu.cn/>) and collects Docking Score, Confidence Score, ligand RMSD and the ranked models.

The submission page is a **standard HTML form POST** (no REST API, no Selenium). The form's `action` attribute is **empty** — JavaScript sets the target at submit time (`runhdock.php` for ordinary docking, `runhsymdock.php` for symmetric multimer). Posting to the homepage returns a **byte-identical copy of it**: HTTP 200, no error, no result.

Submission and retrieval are separate: a job "can normally complete within 30 min", so submission returns a job URL in ~2 s and retrieval reads the score later.

## When to Use This Skill

- Predicting a protein-protein complex structure with HDOCK
- Getting the Docking Score / Confidence Score for a candidate complex
- Screening many receptor-ligand pairs through the server
- Downloading the top-ranked docked models

## Quick Start

```bash
export PYTHONPATH="$PWD/packages"
python scripts/submit_docking.py \
    --receptor-pdb-id 1CGI:E --ligand-pdb-id 1CGI:I \
    --job-name my_docking --output-dir results/
python scripts/submit_docking.py \
    --receptor-seq receptor.fasta --ligand-seq ligand.fasta \
    --output-dir results/
python scripts/fetch_docking_results.py \
    --url http://hdock.phys.hust.edu.cn/data/xxxxxxxxxxxx/ \
    --output-dir results/
python scripts/fetch_docking_results.py \
    --jobs-file results/jobs.csv --output-dir results/ --download-models
```

## Workflow

1. **Submission — one input per molecule.** Each of the receptor and the ligand takes exactly one input form (see Parameters); `--*-seq` takes inline text or a path.
2. **Submit.** The skill GETs the homepage first (matching browser behaviour), then POSTs the multipart form to `runhdock.php` — or `runhsymdock.php` with `--symmetry`.
3. **Record.** Job name, ID and URL go to `jobs.csv` and `jobs.json`; the command returns at once.
4. **Retrieval — read the page.** Three states: Queued (`is QUEUED`) → `queued`; Finished (`Docking Score` table present) → `done`; Failed (`ERROR:` message) → `failed`. The old script treated "not ready yet" as "failed" and aborted the batch.
5. **Parse the ranking matrix.** Rows are metrics (Docking Score, Ligand rmsd), columns are ranks 1–10; the full row is extracted, so all top-10 scores are available.
6. **Compute Confidence Score if needed.** HDOCK gives `Confidence_score = 1.0 / [1.0 + e^(0.02 * (Docking_Score + 150))]`; when the page omits it the skill computes it from the Docking Score.
7. **Optionally download models.** `--download-models` fetches the top-ranked models.

## Parameters

`submit_docking.py`

- `--receptor-pdb-id` (—) — as `PDBID:ChainID` (e.g. `1CGI:E`); `--receptor-seq` (—) — FASTA; `--receptor-pdb` (—) — PDB path
- `--ligand-pdb-id` (—) — as `PDBID:ChainID`; `--ligand-seq` (—) — FASTA; `--ligand-pdb` (—) — PDB path
- `--job-name` (timestamp) — job label; `--email` (none) — server emails on completion
- `--ligand-type` (`protein`) — `protein`, `rna`, `dna`, `rna.d`, `ssdna`; `--symmetry` (none) — e.g. `C2`, `D3`; symmetric endpoint
- `--receptor-site` / `--ligand-site` (none) — binding-site residues, comma separated; `--restraint` (none) — distance restraints
- `--output-dir` (`hdock_results`) — job records; `--timeout` (`180`) — timeout, s

`fetch_docking_results.py`

- `--url` (—) — job result URL; `--jobs-file` (—) — job CSV, in batch; `--output-dir` (`hdock_results`) — output dir
- `--wait` (off) — poll until finished; `--poll-interval` (`20`) — polling interval, s; `--max-wait` (`3600`) — max wait, s
- `--download-models` (off) — download the top-ranked models; `--top-n` (`5`) — how many; `--timeout` (`60`) — timeout, s

## Outputs

- `jobs.csv` / `jobs.json` — submitted jobs: name, ID, URL, status
- `results.csv` / `results.json` — results: status, Docking Score, Confidence Score, error

## Validation

- **A Docking Score is not a binding affinity** — HDOCK is explicit: the score "should not be treated as the true binding affinity of two molecules because it has not been calibrated to the experimental data." Use it to *rank* poses, not to predict a Kd.
- **A more negative score is better** — PDB complexes normally score around −200 or better.
- **Confidence Score above 0.7 means "likely to interact"** — a monotone function of the Docking Score, so it adds no information, just an interpretable scale.
- **A queued job is not a failed job** — `status=queued` means the server is still computing. Jobs take 5–30 min; re-run retrieval later, don't resubmit.
- **"No structure found" means the input was unusable**, not that the interaction is impossible — usually a bad chain ID or non-standard sequence characters. Check the chain ID format (`1CGI:E` — 4-char PDB ID, colon, chain); other causes in `references/troubleshooting.md`.
- **Results are stored for two weeks** by the server. Retrieve within that window or resubmit.
- **Sequence input is modelled, not docked as-is** — HDOCK builds a structure from the sequence first (homology modelling), then docks it. Upload PDB structures when you have them: better accuracy, faster job.
- **Bind the score to the job, not the URL alone** — job names are normalised by the server (spaces → underscores); `job_id` is stable.

## References

- `references/methodology.md` — submission protocol, form fields, page state machine, scoring functions, HDOCK papers
- `references/troubleshooting.md` — common failures and diagnosis
