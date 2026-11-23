---
name: molecular-dynamics
description: "Use when trajectory RMSD or RMSF, PCA, conformational clustering, hydrogen-bond occupancy, Ramachandran plots or solvent-accessible surface analysis is requested. Analyses trajectories with a self-contained pure-numpy DCD and PDB parser plus built-in Kabsch-Sander alignment, SVD-PCA, hierarchical clustering, geometric hydrogen-bond detection, radius of gyration, inter-chain contacts and Shrake-Rupley SASA. 用于轨迹 RMSD/RMSF/PCA 与构象分析。"
license: MIT
compatibility: Python 3.9+; requires numpy; optional scipy, MDAnalysis, OpenMM
metadata:
  version: "1.1"
  last-reviewed: "2026-10-25"
  skill-author: sinder
  edam-topic: "topic_0176, topic_3892"
  edam-operation: "operation_2476, operation_0244, operation_3890"
  edam-input-data: "data_3870, data_1460"
  edam-output-data: "data_3842, data_1544, data_1710"
  edam-input-format: "format_3867, format_1475"
  edam-output-format: "format_3752, format_3464"

# Molecular Dynamics Trajectory Analysis


---
## Overview

Turns MD trajectories into conformational descriptors. Pure-numpy DCD/PDB/XYZ parser (DCD = Fortran record stream: 84-byte header then per-frame coordinates) reads NAMD/CHARMM/GROMACS/LAMMPS trajectories without MDAnalysis (used when present); `run_simulation.py` drives OpenMM.

## When to Use This Skill

- RMSD / RMSF / Rg from a DCD trajectory
- PCA of conformational motion; clustering into states
- H-bond counts and occupancy; Ramachandran / secondary structure
- Inter-chain stability; SASA, buried-surface estimation
- Reading DCD without MDAnalysis
- Short MD simulation with OpenMM

## Quick Start

```bash
python scripts/run_analysis.py --input traj.dcd --topology system.pdb --output results/
python scripts/run_analysis.py --input traj.dcd --topology system.pdb --analysis rmsd,rmsf --output results/
python scripts/run_analysis.py --input traj.dcd --topology system.pdb --start 500 --stride 2 --output results/
python scripts/run_analysis.py --input traj.dcd --topology system.pdb --analysis cluster --n-clusters 5 --output results/
python scripts/run_simulation.py --input system.pdb --output results/ --steps 100000
python scripts/run_simulation.py --input complex.pdb --output results/ --nvt-steps 50000 --npt-steps 50000 --restraints backbone --monitor-chain-separation --steps 500000
python scripts/run_simulation.py --input complex.pdb --output results/ --nvt-steps 50000 --restraints backbone --production-restraint-force-constant 10
```

Cmd 3 = production only; 6 = NVT→NPT→production; 7 = k=10 kJ/mol/Å² restraints.

## Workflow

1. **Load** — `trajectory_io.py`: MDAnalysis else built-in parsers.
2. **Select** — `protein`, `name CA`, `resname LIG`, `resid 1-50`, `chain A` → indices, masses.
3. **Superpose** — `kabsch_align` per frame → rotation, translation, RMSD.
4. **Analyse** — `--analysis` set; RMSD per-residue.
5. **Report** — `md_analysis.json`, `md_analysis_frames.csv`.

Simulation path (`references/methodology.md` §12) — avoids chain drift, never-relaxed packing:

6. **Minimise** — `simulation.minimizeEnergy()`.
7. **NVT** (`--nvt-steps`) — constant volume, velocities redrawn at T; `--restraints` holds backbone.
8. **NPT** (`--npt-steps`) — `MonteCarloBarostat` added, simulation **rebuilt** so coupling hits density.
9. **Production** — restraints off (`--production-restraint-force-constant` keeps weak); `--monitor-chain-separation` tracks distance.
10. **Stability** — NaN / Inf / blow-up screened; if found, rollback to input (logged).

> **Implementation note.** Position restraints must be added to the `System` **before** the `Simulation`/`Context` is created — modifying `system` afterwards does not rebuild the context, so the force would silently never take effect (same for the barostat and for removing restraints; `references/methodology.md` §12.4).

## Parameters

`run_analysis.py`:

|Parameter|Default|Description|
|---|---|---|
|`--input`|required|DCD/PDB/XYZ|
|`--topology`|—|Topology|
|`--output`|`results/`|Output dir|
|`--start` / `--stop` / `--stride`|`0` / — / `1`|Frame range|
|`--selection`|`protein`|Selection|
|`--analysis`|`all`|`rmsd,rmsf,rg,hbond,pca,dihedral,ramachandran,cluster,interchain,sasa`|
|`--fit-reference`|`0`|RMSD ref|
|`--hbond-angle`|`120.0`|Min D-H···A (deg)|
|`--hbond-h-dist`|`2.5`|Max H···A (Å)|
|`--hbond-da-dist`|`3.5`|Max D···A (Å)|
|`--n-clusters`|auto|Clusters|
|`--cluster-threshold`|auto|RMSD thr (Å); beats `--n-clusters`|
|`--contact-dist`|`4.5`|Contact (Å)|
|`--sasa-points`|`100`|SASA pts/atom|
|`--max-cluster-frames`|`200`|Max frames|
|`--format`|`both`|`json` / `csv` / `both`|

`run_simulation.py`:

|Parameter|Default|Description|
|---|---|---|
|`--force-field`|`amber14`|`amber14` / `amber99` / `charmm36`|
|`--solvent`|`implicit`|`none` / `implicit` / `vacuum`|
|`--steps`|`1000`|Prod steps|
|`--equilibration-steps`|`100`|After NVT/NPT|
|`--minimize-steps`|`200`|Min steps|
|`--temperature`|`300.0`|(K)|
|`--friction`|`1.0`|(/ps)|
|`--timestep`|`0.002`|(ps)|
|`--platform`|—|CUDA / CPU / OpenCL|
|`--nvt-steps` / `--npt-steps`|`0` / `0`|0 = skip|
|`--pressure`|`1.0`|(bar)|
|`--restraints`|`none`|`none` / `backbone` / `ca` / `heavy` / `all`|
|`--restraint-force-constant`|`100.0`|kJ/mol/Å²|
|`--production-restraint-force-constant`|`0.0`|>0 keeps|
|`--monitor-chain-separation`|off|Chain dist|
|`--chain-separation-threshold`|`30.0`|(Å)|
|`--chain-monitor-interval`|`10000`|(steps)|
|`--check-only`|off|OpenMM check|

## Outputs

Analysis: `md_analysis.json`, `md_analysis_frames.csv` (frame, time, RMSD, Rg, H-bonds, SASA). Simulation: `trajectory.dcd`, `final_state.pdb`, `state.csv` (energy, temperature, volume), `simulation_report.json` (`stages`/`restraints`/`stability`/`chain_separation`).

## Validation

- **RMSD** — `--fit-reference` frame RMSD ≈ 0.
- **RMSD/drift** — 0.02 Å/frame ⇒ linear growth; flat ⇒ superposition absorbs it.
- **RMSF** — Cα RMSF 0.5-2 Å; >5 Å ⇒ unstable/badly parameterised, <0.01 Å ⇒ frozen frames.
- **Rg** — `Rg ≈ 2.2 × N^0.38` Å globular; jump ⇒ unfolding.
- **H-bonds** — 50-200 backbone-sidechain + backbone-backbone H-bonds/frame; zero with donors/acceptors ⇒ topology lacks H (X-ray).
- **PCA** — first component 30-60% of variance; >95% ⇒ translation/rotation not removed.
- **Clustering** — count matches the RMSD matrix; one wide-RMSD cluster ⇒ threshold too large.
- **Ramachandran** — <0.5% outliers in allowed regions; more ⇒ force-field or backbone problems.
- **SASA** — 5000-7000 Å² for a 100-residue globular protein; much lower ⇒ collapsed/over-packed.
- **DCD** — `n_atoms` must match topology count; mismatch ⇒ DCD and topology differ.
- **Stages** — `stages` lists every requested stage; absent `--nvt-steps` ⇒ silently skipped.
- **Chain separation** — block exists, `final_check.separated` must be `false`; `true` ⇒ chains drifted apart (never-relaxed pose).
- **Restraint effect** — NVT backbone RMSD smaller than unrestrained run; identical ⇒ restraints never took effect (force added to `System` after `Context`).
- **Stability** — `stability.should_rollback` must be `false`; if `true`, `stability.reason` records NaN / Inf / blow-up; NaN in first frames ⇒ bad force field or missing H.

## References

- `references/methodology.md` — algorithms, formulas
- `references/troubleshooting.md` — failures, diagnosis
