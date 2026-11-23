---
name: alphafold-db
description: "Use when fetching AlphaFold predicted structures, retrieving models by UniProt accession, or collecting AlphaFold confidence metrics (pLDDT). Downloads PDB/CIF models from the EBI AlphaFold Database API, resolving the current URL from the API instead of hard-coding a model version; supports batch accession lists, resumable downloads, pLDDT filtering, and content-type validation. 用于下载 AlphaFold 结构与 pLDDT。"
license: MIT
compatibility: Python 3.9+; needs requests for network access; no third-party package required for offline metadata processing; biopython optional for CIF→PDB conversion (a built-in fallback converter is used when it is absent)
metadata:
  version: "1.1"
  last-reviewed: "2026-10-25"
  skill-author: sinder
  edam-topic: "topic_0176, topic_0203"
  edam-operation: "operation_0474"
  edam-input-data: "data_1460"
  edam-output-data: "data_2537"
  edam-input-format: "format_1475"
  edam-output-format: "format_1475"

# AlphaFold DB Structure Retrieval


---
## Overview

The model version suffix of AlphaFold DB files is **not stable** — the DB retrains and increments it (`model_v6` → `model_v7`), so hard-coded `AF-{ACCESSION}-F1-model_v6.pdb` URLs silently start 404ing. It calls `GET /api/prediction/{accession}` and takes `pdbUrl` from the response, so the server always points at the current model; the response also carries free pLDDT metrics (global plus per-band fractions) for filtering before download. It also selects the best model of a multi-sample prediction and converts it to PDB.

## When to Use This Skill

- Fetching predicted structures for a list of UniProt accessions
- Collecting pLDDT confidence metrics without downloading structures
- Screening a proteome for low-confidence regions before downloads
- Refreshing a collection after AlphaFold DB publishes new models
- Populating a docking or interaction-analysis project with coordinates
- Selecting the best-confidence model from an AlphaFold Server / ColabFold / AlphaFold DB package, as PDB

## Quick Start

```bash
export PYTHONPATH="$PWD/packages"

python scripts/download_structures.py --ids ids.txt --output-dir structures/
python scripts/download_structures.py --ids P69905 P04637 --output-dir structures/
python scripts/download_structures.py --ids ids.txt --metadata-only --output-dir meta/
python scripts/download_structures.py --ids ids.txt --min-plddt 70 --output-dir structures/
python scripts/download_structures.py --ids ids.txt --output-dir structures/ --keep-original-name
python scripts/select_best_model.py --input-dir folds/ --output-dir best/
python scripts/select_best_model.py --input-dir folds/ --rank-only
```

## Workflow

1. **Collect accessions** — `--ids` file (one per line, `#` comments) or accession strings; a header line is auto-detected by pattern, so a pure-ID file keeps its first entry.
2. **Query the API** — `alphafold_client.fetch_predictions()`, per accession via `http_session` (backoff retry + rate limiting); response is a **JSON array**; list / dict / error shapes normalise to one record type.
3. **Export metadata** — each record to `alphafold_metadata.csv` and `.json`, with `global_plddt` and the four band fractions.
4. **Filter (optional)** — `--min-plddt`: models below the threshold are counted `filtered`, never downloaded.
5. **Download** — each `pdbUrl` fetched with `Content-Type` checked against expected PDB types, payload written as **bytes**; existing files skipped unless `--force`.
6. **Report** — downloaded / skipped / filtered / failed counts plus failures.

## Parameters

### `download_structures.py`

| Parameter | Default | Description |
|---|---|---|
| `--ids` | required | Accession list, or path(s) to accessions files |
| `--output-dir` | `alphafold_structures` | Output dir for structures and metadata |
| `--metadata-only` | off | Metadata only; no downloads |
| `--min-plddt` | none | Skip models whose global pLDDT is below this value (0–100) |
| `--force` | off | Re-download existing files |
| `--keep-original-name` | off | Keep `AF-{ID}-F1-model_vN.pdb` instead of `{ID}.pdb` |
| `--timeout` | `120` | Per-request timeout in seconds |

### `select_best_model.py`

| Parameter | Default | Description |
|---|---|---|
| `--input-dir` | required | Prediction package dir (structures + `summary_confidences_*.json`) |
| `--output-dir` | `<input>_selected` | Where the selected model and ranking tables go |
| `--rank-only` | off | Ranking table only; convert nothing |
| `--rank-by` | `ranking_score` | Ordering metric: `ranking_score` / `ptm` / `iptm` / `mean_plddt` |
| `--keep-all` | off | Convert every model, not only the best |
| `--no-convert` | off | Copy originals without CIF → PDB conversion |

## Outputs

| File | Content |
|---|---|
| `{accession}.pdb` | Downloaded prediction (or `AF-{ID}-F1-model_vN.pdb` with `--keep-original-name`) |
| `alphafold_metadata.csv` | Accession, gene, description, organism, model version, global pLDDT, bands, URLs |
| `alphafold_metadata.json` | Same as nested JSON |
| `model_ranking.csv` | Position, rank, file names, ranking_score, pTM, ipTM, mean pLDDT, disorder, clash |
| `model_ranking.json` | Same as nested JSON (plus `rank_by`) |
| `{structure_stem}.pdb` | Selected model as PDB |

## Validation

- **File is really a PDB** — first lines contain `HEADER` or `ATOM`, last line `END`; HTML (`<!DOCTYPE html>`) = error page — check `Content-Type` and small files.
- **pLDDT range** — `global_plddt` (mean per-residue pLDDT): globular 70–95, <50 = low-confidence regions (disorder), >95 unusual.
- **Bands sum to 1** — `fraction_plddt_very_low + fraction_plddt_low + fraction_plddt_confident + fraction_plddt_very_high` = 1.0 (within rounding); mismatch = truncated record.
- **UniProt cross-check** — metadata `sequence_length` must equal the UniProt entry length; mismatch = wrong isoform.
- **Coverage** — ~230 million entries cover the UniProt reference proteomes; many "no prediction" hits ⇒ accessions must be canonical, not isoform IDs.
- **ranking_score ranks, not accuracy** — `0.8 × ipTM + 0.2 × pTM + 0.5 × fraction_disordered − 100 × has_clash`; it orders predictions, so `ranking_score = 1.0` = best of its batch, not 100%-accurate.
- **PDB column alignment** — fixed-width records: chain ID col 22, residue 23–26, coordinates 31–38 / 39–46 / 47–54; misaligned records are silently misread downstream — check `awk 'substr($0,22,1)'`.

## References

- `references/methodology.md` — API semantics, version resolution, pLDDT interpretation
- `references/troubleshooting.md` — common failures and diagnosis
