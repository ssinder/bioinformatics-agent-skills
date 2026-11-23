---
name: protein-docking-lightdock
description: "Use when protein-protein docking must be performed on the local machine, when a swarm-level score ranking is needed, or when the best docked complex model must be generated as a PDB file. Runs lightdock3 locally and extracts the best complex with correct GSO field semantics, numeric step selection, and a single consistent coordinate frame; supports batch grids, resumable execution and DFIRE ranking. 用于本地 LightDock 对接。"
license: MIT
compatibility: Python 3.9+; requires the lightdock package (not available on Windows); no third-party Python package required for parsing and scoring
metadata:
  version: "1.0"
  last-reviewed: "2026-10-25"
  skill-author: sinder
  edam-topic: "topic_0176"
  edam-operation: "operation_3899"
  edam-input-data: "data_1460"
  edam-output-data: "data_1460, data_1537"
  edam-input-format: "format_1475, format_1929"
  edam-output-format: "format_1475"

# LightDock Local Protein Docking


---

## Overview

Runs [LightDock](https://lightdock.org/), a glowworm-swarm-optimization protein-protein docking engine, **on the local machine**, turning its raw swarm output into a ranked score table and docked complexes.

LightDock differs from the HDOCK skill here: HDOCK is a **web service** (submit a job, the server computes), LightDock is a **local program** (`pip install lightdock`, the search runs on your CPUs). One job with 50 swarms × 50 steps takes minutes to hours; 100 pairs take days — hence the **long-running, resumable, interrupt-safe** design.

The hard part of driving LightDock from a script is not the subprocess calls — it is the **output semantics**. The `gso_<step>.out` files and the `lightdock_<i>.pdb` conformations carry three traps that all produce plausible but wrong results; this skill fixes them.

## When to Use This Skill

- Running protein-protein docking locally with full control over swarm/step counts
- Ranking a batch of receptor-ligand pairs by DFIRE score
- Recovering the best complex model from an existing LightDock run
- Producing docked complexes for downstream MD, interaction analysis, or figures
- Comparing against `protein-docking-hdock` when a **local** (no data leaves the machine) run is needed

## Quick Start

```bash
export PYTHONPATH="$PWD/packages"

python scripts/run_docking.py --receptor rec.pdb --ligands ligands/*.pdb \
    --output-dir results/
python scripts/run_docking.py --receptors recs/*.pdb --ligands ligands/*.pdb \
    --output-dir results/
python scripts/run_docking.py --receptor rec.pdb --ligands a.pdb b.pdb \
    --swarms 5 --glowworms 10 --steps 20 --cores 4 --output-dir test/
python scripts/generate_conformations.py --job-dir results/rec__lig/
python scripts/extract_scores.py --jobs-root results/ --output-dir scores/
```

## Workflow

1. **Locate the executables** — `resolve_executables()` finds `lightdock3_setup.py` and `lightdock3.py` by checking, in order: an explicit `--setup-exe` / `--run-exe` path, the interpreter's `<prefix>/bin` (where conda and venv installs put them — **not on `PATH`**, which is why the old scripts hard-coded a path), then `PATH`.
2. **Setup stage** — `lightdock3_setup.py` copies both structures into the job directory, strips waters/oxygens/hydrogens per the flags, computes the receptor's minimum-volume ellipsoid, places the swarms on its surface, and writes `setup.json`. Critically, it **translates both structures to the origin**, saving them as `lightdock_receptor.pdb` / `lightdock_ligand.pdb`.
3. **Simulation** — `lightdock3.py setup.json <steps> -c <cores>` runs the glowworm swarm optimisation; each swarm writes `gso_<step>.out` every 10 steps.
4. **Score extraction** — `collect_swarm_results()` picks the **highest-step** `gso_*.out` per swarm by numeric step, parses every glowworm line, and finds the minimum-scoring entry. **Lower DFIRE score is better.**
5. **Conformation generation** — `lgd_generate_conformations.py` writes `lightdock_<i>.pdb` files from the GSO coordinates, using the origin-shifted structures so receptor and ligand stay aligned.
6. **Reporting** — scores go to `docking_scores.csv` / `.json`, written incrementally after each job.

## Parameters

`run_docking.py`

- `--receptor` (—) — single receptor PDB path; `--receptors` (—) — multiple receptor paths (a grid with `--ligands`); `--ligands` (required) — ligand paths, globs, or directories; `--output-dir` (`lightdock_results`) — output root
- `--swarms` (`50`) — swarms (search breadth); `--glowworms` (`50`) — per swarm (sampling density); `--steps` (`50`) — simulation steps; `--cores` (`4`) — CPU cores
- `--resume` (off) — skip jobs with a valid result; `--keep-hydrogens` / `--keep-oxt` / `--keep-water` (off) — keep stripped atoms
- `--setup-timeout` (`600`) — setup timeout (s); `--run-timeout` (`86400`) — simulation timeout (s); `--setup-exe` / `--run-exe` (auto) — paths to the LightDock scripts

`generate_conformations.py`

- `--job-dir` (—) — a single job directory; `--jobs-root` (—) — root whose subdirectories are jobs; `--output-dir` (in-place) — where to copy conformations
- `--top` (`1`) — how many conformations (cost scales linearly); `--ligand-chain` (—) — rewrite the ligand chain ID in a **copy** (original untouched); `--conformations-exe` (auto) — path to `lgd_generate_conformations.py`

`extract_scores.py`

- `--job-dir` / `--jobs-root` (—) — which jobs to summarise; `--output-dir` (next to jobs) — where score tables go; `--per-swarm` (off) — one row per swarm, not per job

## Outputs

- `docking_scores.csv` — one row per job: job id, receptor, ligand, status, best score, best swarm, glowworm
- `docking_scores.json` — same content as nested JSON
- `conformations.csv` — per job: chosen swarm, glowworm, score, generated PDB names
- `{job}/swarm_<n>/gso_<step>.out` — raw glowworm swarm output
- `{job}/swarm_<n>/lightdock_<i>.pdb` — docked complex (receptor + ligand merged)

## Validation

- **Score sign** — DFIRE scores are **negative and lower is better**. A positive "best score" means the parser picked the wrong column.
- **Score magnitude** — a converged small-complex DFIRE score lands around `-30` to `-60`; near `0` means the search barely moved from the random start.
- **glowworm index vs LigID** — the glowworm id is the **line number** in the GSO file. In single-structure docking `RecID` and `LigID` are both `0`. "best glowworm = 0" everywhere means the parser reads the wrong field.
- **GSO step selection** — use the file with the **largest numeric step** (`gso_100.out` over `gso_20.out`), not the last alphabetically.
- **Conformation count** — `lgd_generate_conformations.py` clips the requested count to the coordinate lines found, and logs it.
- **Complex sanity** — open the generated PDB and confirm the two parts touch but do not interpenetrate. A ligand flying to the origin means the origin-shifted structures were not used.
- **RMSD** — LightDock does not report RMSD by default; run `lgd_rank.py` if needed, or compare against a known holo structure.

## References

- `references/methodology.md` — GSO format semantics, coordinate-frame handling, score extraction, command reference
- `references/troubleshooting.md` — common failures and how to diagnose them
