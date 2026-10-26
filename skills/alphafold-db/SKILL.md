---
name: alphafold-db
description: >-
  Downloads predicted protein structures from the EBI AlphaFold Database API and
  exports per-protein prediction metadata. Resolves the current PDB/CIF download
  URL from the API instead of hard-coding a model version, supports batch
  accession lists, resumable downloads, pLDDT-based quality filtering before
  spending bandwidth, and validates response content types so HTML error pages
  are never mistaken for structure files. This skill should be used when fetching
  AlphaFold models, retrieving predicted structures by UniProt accession, or
  collecting AlphaFold confidence metrics (pLDDT) is requested.
license: MIT
compatibility: Python 3.9+; needs requests for network access; no third-party package required for offline metadata processing; biopython optional for CIF→PDB conversion (a built-in fallback converter is used when it is absent)
metadata:
  version: "1.1"
  last-reviewed: "2026-10-25"
  skill-author: sinder
---

# AlphaFold DB Structure Retrieval

## Overview

This skill retrieves AlphaFold protein structure predictions from the EBI AlphaFold
Database and exports the accompanying confidence metadata. It addresses the single
biggest fragility in naive AlphaFold download scripts: **the model version number is
not stable**. The AlphaFold DB retrains models periodically and increments the
version suffix (`model_v6` → `model_v7`, …), so any script that builds the URL as
`AF-{ACCESSION}-F1-model_v6.pdb` will silently start returning 404s.

Instead, this skill calls `GET /api/prediction/{accession}` and takes the `pdbUrl`
field straight from the response — the server always points at the current valid
model. The same response also carries per-protein confidence metrics (global pLDDT
and the residue fractions in each confidence band) which are **free** and can be used
to filter low-quality models *before* downloading anything.

The skill also handles the second half of the AlphaFold workflow: **selecting the
best model from a multi-sample prediction** and converting it to PDB format.

## When to Use This Skill

- Fetching predicted structures for a list of UniProt accessions
- Collecting AlphaFold confidence metrics (pLDDT) without downloading structures
- Screening a proteome for low-confidence regions before committing to downloads
- Refreshing an existing structure collection after AlphaFold DB publishes new models
- Populating a docking or interaction-analysis project with predicted coordinates
- Picking the highest-confidence model out of an AlphaFold Server / ColabFold /
  AlphaFold DB prediction package and converting it to PDB

## Quick Start

```bash
export PYTHONPATH="$PWD/packages"

# 批量下载（ID 文件，每行一个）
python scripts/download_structures.py --ids ids.txt --output-dir structures/

# 直接给 ID
python scripts/download_structures.py --ids P69905 P04637 --output-dir structures/

# 只要元数据，不要 PDB（最快，用于先做质量筛选）
python scripts/download_structures.py --ids ids.txt --metadata-only --output-dir meta/

# 先按 pLDDT 过滤再下载，跳过低质量模型
python scripts/download_structures.py --ids ids.txt --min-plddt 70 --output-dir structures/

# 保留 AlphaFold 原始文件名
python scripts/download_structures.py --ids ids.txt --output-dir structures/ --keep-original-name

# 从多模型预测包里选最佳模型并转成 PDB
python scripts/select_best_model.py --input-dir folds/ --output-dir best/

# 只要排名表，不做转换
python scripts/select_best_model.py --input-dir folds/ --rank-only
```

## Workflow

1. **Collect accessions** — `--ids` accepts either a text file (one accession per line,
   `#` for comments) or accession strings directly. A leading header line is detected
   by pattern rather than skipped unconditionally, so a pure-ID file never loses its
   first entry.
2. **Query the API** — `alphafold_client.fetch_predictions()` issues one request per
   accession through the shared `http_session` (exponential-backoff retry + rate
   limiting). The API returns a **JSON array** for a single accession, not an object;
   the client normalises list / dict / error-object shapes into one record type.
3. **Export metadata** — every record is written to `alphafold_metadata.csv` and
   `.json`, including `global_plddt` and the four confidence-band fractions.
4. **Filter (optional)** — with `--min-plddt`, models below the threshold are counted
   as `filtered` and never downloaded.
5. **Download** — each `pdbUrl` is fetched with the response `Content-Type` checked
   against the expected PDB types and the payload written as **bytes**. Files that
   already exist are skipped unless `--force` is given.
6. **Report** — a summary line reports downloaded / skipped / filtered / failed counts
   and lists the failed accessions.

## Parameters

### `download_structures.py`

| Parameter | Default | Description |
|---|---|---|
| `--ids` | required | Accession list or path(s) to text files containing accessions |
| `--output-dir` | `alphafold_structures` | Output directory for structures and metadata |
| `--metadata-only` | off | Query and export metadata only; skip file downloads |
| `--min-plddt` | none | Skip models whose global pLDDT is below this value (0–100) |
| `--force` | off | Re-download files that already exist |
| `--keep-original-name` | off | Keep `AF-{ID}-F1-model_vN.pdb` instead of `{ID}.pdb` |
| `--timeout` | `120` | Per-request timeout in seconds |

### `select_best_model.py`

| Parameter | Default | Description |
|---|---|---|
| `--input-dir` | required | Directory holding the prediction package (structures + `summary_confidences_*.json`) |
| `--output-dir` | `<input>_selected` | Where the selected model and ranking tables are written |
| `--rank-only` | off | Only write the ranking table; convert nothing |
| `--rank-by` | `ranking_score` | Metric used for ordering: `ranking_score` / `ptm` / `iptm` / `mean_plddt` |
| `--keep-all` | off | Convert every model instead of only the best one |
| `--no-convert` | off | Copy the original files without converting CIF → PDB |

## Outputs

| File | Content |
|---|---|
| `{accession}.pdb` | Downloaded AlphaFold prediction (or `AF-{ID}-F1-model_vN.pdb` with `--keep-original-name`) |
| `alphafold_metadata.csv` | One row per protein: accession, gene, description, organism, model version, global pLDDT, confidence-band fractions, structure URLs |
| `alphafold_metadata.json` | Same content as nested JSON |
| `model_ranking.csv` | One row per model: position, rank, file names, ranking_score, pTM, ipTM, mean pLDDT, disorder, clash |
| `model_ranking.json` | Same content as nested JSON (plus the `rank_by` key used) |
| `{structure_stem}.pdb` | The selected model converted to PDB format |

## Validation

- **A downloaded file is really a PDB** — the first lines must contain `HEADER`
  or `ATOM` records, and the last line should be `END`. A file containing HTML
  (`<!DOCTYPE html>`) means the server returned an error page; the script's
  Content-Type check rejects most such cases, but verify if a file is unusually small.
- **pLDDT sanity** — AlphaFold's mean per-residue pLDDT is `global_plddt`. Globular
  proteins typically sit at 70–95; values below 50 indicate extensive low-confidence
  regions (often disorder), and values above 95 are unusual.
- **Confidence bands sum to 1** — `fraction_plddt_very_low + fraction_plddt_low +
  fraction_plddt_confident + fraction_plddt_very_high` should equal 1.0 (within
  rounding). A mismatch means the record was truncated.
- **Cross-check against UniProt** — the `sequence_length` in the metadata must equal
  the UniProt entry length. A mismatch indicates the accession was mapped to a
  different isoform.
- **Coverage** — roughly 230 million entries cover the UniProt reference proteomes;
  if many accessions report "no prediction", confirm the accessions are canonical
  rather than isoform or isoform-variant identifiers.
- **ranking_score is a ranking metric, not an accuracy metric** — the official
  definition is `0.8 × ipTM + 0.2 × pTM + 0.5 × fraction_disordered − 100 × has_clash`.
  It is designed to *order* predictions, not to be interpreted as "how good" the
  structure is. A structure with `ranking_score = 1.0` is the best of its batch,
  not a 100%-accurate model.
- **PDB column alignment** — PDB records are fixed-width. A converted file must have
  the chain ID at column 22, the residue number at 23–26, and the coordinates at
  31–38 / 39–46 / 47–54. Misaligned records are silently misread by every
  downstream tool, so check with `awk 'substr($0,22,1)'` if a structure looks off.

## References

- `references/methodology.md` — API semantics, version-resolution strategy, pLDDT interpretation
- `references/troubleshooting.md` — common failures and how to diagnose them
