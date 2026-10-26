---
name: protein-docking-lightdock
description: >-
  Runs LightDock protein-protein docking locally by driving the lightdock3_setup.py
  and lightdock3.py executables, then extracts the best docked complex. Parses the
  GSO swarm output files with the correct field semantics (glowworm id is the line
  number, not the LigID column), selects the latest gso_<step>.out by numeric step
  rather than dictionary order, and generates conformations from the
  origin-shifted structures that the setup stage produced so the receptor and
  ligand stay in a consistent coordinate frame. Supports batch receptor-ligand
  grids, resumable execution, incremental score logging, and DFIRE-based ranking.
  This skill should be used when protein-protein docking must be performed on the
  local machine, when a swarm-level score ranking is needed, or when the best
  docked complex model must be generated as a PDB file.
license: MIT
compatibility: Python 3.9+; requires the lightdock package (not available on Windows); no third-party Python package required for parsing and scoring
metadata:
  version: "1.0"
  last-reviewed: "2026-10-25"
  skill-author: sinder
---

# LightDock Local Protein Docking

## Overview

This skill runs [LightDock](https://lightdock.org/), a glowworm-swarm-optimization
protein-protein docking engine, **on the local machine** and turns its raw swarm
output into a ranked score table and a set of docked complex structures.

LightDock works differently from the HDOCK skill in this same repository. HDOCK
is a **web service**: you submit a job and the server computes. LightDock is a
**local program**: `pip install lightdock`, and the search runs on your own CPUs.
A single job with 50 swarms × 50 steps takes minutes to hours depending on core
count and protein size. A batch of 100 pairs takes days. This skill is therefore
built around **long-running, resumable, interrupt-safe** batch execution rather
than around request/response handling.

The hard part of driving LightDock from a script is not the subprocess calls —
it is the **output semantics**. The `gso_<step>.out` files and the
`lightdock_<i>.pdb` conformations have three separate traps that all produce
plausible-looking but wrong results, and this skill fixes all three.

## When to Use This Skill

- Running protein-protein docking locally with full control over swarm/step counts
- Ranking a batch of receptor-ligand pairs by DFIRE score
- Recovering the best complex model from an existing LightDock run
- Producing docked complex structures for downstream MD, interaction analysis, or figure preparation
- Comparing against the `protein-docking-hdock` skill when a **local** (no data
  leaves the machine) prediction is required

## Quick Start

```bash
export PYTHONPATH="$PWD/packages"

# 一个受体 vs 一组配体
python scripts/run_docking.py --receptor rec.pdb --ligands ligands/*.pdb \
    --output-dir results/

# 多个受体 vs 多个配体（笛卡尔积）
python scripts/run_docking.py --receptors recs/*.pdb --ligands ligands/*.pdb \
    --output-dir results/

# 快速试跑（先验证流程再上规模）
python scripts/run_docking.py --receptor rec.pdb --ligands a.pdb b.pdb \
    --swarms 5 --glowworms 10 --steps 20 --cores 4 --output-dir test/

# 从已有结果生成最优复合物
python scripts/generate_conformations.py --job-dir results/rec__lig/

# 汇总分数排名
python scripts/extract_scores.py --jobs-root results/ --output-dir scores/
```

## Workflow

1. **Locate the executables** — `resolve_executables()` finds `lightdock3_setup.py`
   and `lightdock3.py` by checking, in order: an explicit `--setup-exe` /
   `--run-exe` path, the current interpreter's `<prefix>/bin` (where conda and
   venv installs put them — **these are not on `PATH`**, which is why the old
   scripts hard-coded an absolute path), and finally the system `PATH`.
2. **Setup stage** — `lightdock3_setup.py` copies both structures into the job
   directory, strips waters/oxygens/hydrogens per the flags, computes the
   minimum-volume ellipsoid of the receptor, places the swarms on its surface,
   and writes `setup.json`. Critically, it **translates both structures to the
   origin** and saves them as `lightdock_receptor.pdb` / `lightdock_ligand.pdb`.
3. **Simulation** — `lightdock3.py setup.json <steps> -c <cores>` runs the glowworm
   swarm optimisation. Each swarm writes `gso_<step>.out` every 10 steps.
4. **Score extraction** — `collect_swarm_results()` picks the **highest-step**
   `gso_*.out` per swarm by numeric step, parses every glowworm line, and finds
   the minimum-scoring entry. **Lower DFIRE score is better.**
5. **Conformation generation** — `lgd_generate_conformations.py` writes
   `lightdock_<i>.pdb` files from the GSO coordinates, using the
   origin-shifted structures so receptor and ligand stay aligned.
6. **Reporting** — scores go to `docking_scores.csv` / `.json`, incrementally
   after every job.

## Parameters

### `run_docking.py`

| Parameter | Default | Description |
|---|---|---|
| `--receptor` | — | Single receptor PDB path |
| `--receptors` | — | Multiple receptor PDB paths (forms a grid with `--ligands`) |
| `--ligands` | required | Ligand PDB paths, glob patterns, or directories |
| `--output-dir` | `lightdock_results` | Output root |
| `--swarms` | `50` | Number of swarms (search space breadth) |
| `--glowworms` | `50` | Glowworms per swarm (sampling density) |
| `--steps` | `50` | Simulation steps |
| `--cores` | `4` | Parallel CPU cores |
| `--resume` | off | Skip jobs that already have a valid result |
| `--keep-hydrogens` / `--keep-oxt` / `--keep-water` | off | Keep stripped atoms |
| `--setup-timeout` | `600` | Setup stage timeout (seconds) |
| `--run-timeout` | `86400` | Simulation timeout (seconds) |
| `--setup-exe` / `--run-exe` | auto | Explicit paths to the LightDock scripts |

### `generate_conformations.py`

| Parameter | Default | Description |
|---|---|---|
| `--job-dir` | — | A single docking job directory |
| `--jobs-root` | — | Root whose subdirectories are treated as jobs |
| `--output-dir` | in-place | Where to copy the generated conformations |
| `--top` | `1` | How many conformations to generate (cost scales linearly) |
| `--ligand-chain` | — | Rewrite the ligand chain ID in a **copy** (original untouched) |
| `--conformations-exe` | auto | Explicit path to `lgd_generate_conformations.py` |

### `extract_scores.py`

| Parameter | Default | Description |
|---|---|---|
| `--job-dir` / `--jobs-root` | — | Which jobs to summarise |
| `--output-dir` | next to jobs | Where the score tables are written |
| `--per-swarm` | off | One row per swarm instead of one row per job |

## Outputs

| File | Content |
|---|---|
| `docking_scores.csv` | One row per job: job id, receptor, ligand, status, best score, best swarm, best glowworm |
| `docking_scores.json` | Same content as nested JSON |
| `conformations.csv` | Per job: chosen swarm, glowworm, score, generated PDB file names |
| `{job}/swarm_<n>/gso_<step>.out` | Raw glowworm swarm output |
| `{job}/swarm_<n>/lightdock_<i>.pdb` | Generated docked complex (receptor + ligand merged) |

## Validation

- **Score sign** — DFIRE scores are **negative and lower is better**. A positive
  "best score" usually means the parser picked the wrong column.
- **Score magnitude** — a converged small-complex DFIRE score typically lands
  around `-30` to `-60`. A score near `0` means the search barely moved from
  the random start.
- **glowworm index vs LigID** — the glowworm id is the **line number** in the
  GSO file. In single-structure docking, `RecID` and `LigID` are both `0`. If
  your report says "best glowworm = 0" for every job, the parser is reading the
  wrong field.
- **GSO step selection** — the file to use is the one with the **largest numeric
  step** (`gso_100.out` over `gso_20.out`), not the last one alphabetically.
- **Conformation count** — `lgd_generate_conformations.py` clips the requested
  count to the number of coordinate lines actually found, and says so in its log.
- **Complex sanity** — open the generated PDB and confirm the two parts touch
  but do not interpenetrate. A ligand flying off to the origin means the
  origin-shifted structures were not used.
- **RMSD** — LightDock does not report RMSD by default; run `lgd_rank.py` if you
  need it, or compare against a known holo structure.

## References

- `references/methodology.md` — GSO format semantics, coordinate-frame handling, score extraction, official command reference
- `references/troubleshooting.md` — common failures and how to diagnose them
