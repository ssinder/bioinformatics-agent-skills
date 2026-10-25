---
name: molecular-dynamics
description: >-
  Analyses molecular-dynamics trajectories with a self-contained pure-numpy DCD
  and PDB trajectory parser plus built-in Kabsch-Sander alignment, SVD-based
  PCA, hierarchical clustering of RMSD matrices, geometric hydrogen-bond
  detection, RMSF, radius of gyration, Ramachandran analysis, inter-chain
  contacts and Shrake-Rupley SASA. This skill should be used when trajectory
  RMSD or RMSF, PCA, conformational clustering, hydrogen-bond occupancy,
  Ramachandran plots or solvent-accessible surface analysis is requested.
license: MIT
compatibility: Python 3.9+; requires numpy; optional scipy, MDAnalysis, OpenMM
metadata:
  version: "1.0"
  last-reviewed: "2026-10-25"
  skill-author: sinder
---

# Molecular Dynamics Trajectory Analysis

## Overview

This skill turns raw MD trajectories into quantitative conformational descriptors.
It ships a **self-contained pure-numpy DCD parser** (DCD is a Fortran-style
fixed-length record stream: an 84-byte header record followed by per-frame
coordinate records), so a trajectory produced by NAMD, CHARMM, GROMACS or LAMMPS
can be analysed without installing MDAnalysis. A multi-model PDB parser and an
XYZ parser are provided as well, and MDAnalysis is used automatically when it
is available.

All core algorithms are implemented directly and are numerically well-defined:

* **Kabsch-Sander** optimal rigid-body superposition via SVD of the covariance
  matrix (with the determinant correction that prevents reflection).
* **PCA** by eigendecomposition of the coordinate covariance matrix.
* **Hierarchical clustering** (average/complete/ward linkage) of the pairwise
  RMSD matrix.
* **Hydrogen bonds** by a distance + angle geometric criterion.
* **Shrake-Rupley** solvent-accessible surface area.
* **Ramachandran** region assignment from backbone φ/ψ dihedrals.

`run_simulation.py` drives a real MD run through OpenMM when OpenMM is
installed (force-field setup, Langevin thermostat, minimisation, equilibration,
production run, DCD output).

## When to Use This Skill

- RMSD / RMSF / radius-of-gyration time series from a DCD trajectory
- Principal component analysis of conformational motion
- Clustering a trajectory into conformational states (RMSD matrix + linkage)
- Hydrogen-bond counts and occupancy across frames
- Ramachandran distribution and secondary-structure validation of a trajectory
- Inter-chain contact stability in a protein–protein or protein–ligand complex
- Solvent-accessible surface area and buried-surface estimation
- Reading DCD files without installing MDAnalysis
- Running a short MD simulation with OpenMM

## Quick Start

```bash
# 基础分析（DCD + 拓扑 PDB）
python scripts/run_analysis.py --input traj.dcd --topology system.pdb \
    --output results/

# 只算 RMSD 与 RMSF
python scripts/run_analysis.py --input traj.dcd --topology system.pdb \
    --analysis rmsd,rmsf --output results/

# 跳过平衡阶段（只看生产段）
python scripts/run_analysis.py --input traj.dcd --topology system.pdb \
    --start 500 --stride 2 --output results/

# 按簇数做构象聚类
python scripts/run_analysis.py --input traj.dcd --topology system.pdb \
    --analysis cluster --n-clusters 5 --output results/

# 用 OpenMM 跑一段模拟
python scripts/run_simulation.py --input system.pdb --output results/ --steps 100000
```

## Workflow

1. **Load** — `trajectory_io.py` reads the trajectory. MDAnalysis is tried
   first; otherwise the built-in parsers handle `.dcd` (Fortran record stream),
   `.pdb` (multi-model) and `.xyz`. Topology (atom names, elements, chains,
   residues) comes from the `--topology` PDB when supplied.
2. **Select** — the selection string (`protein`, `name CA`, `resname LIG`,
   `resid 1-50`, `chain A`) is resolved to an atom index array; masses are
   attached for mass-weighted quantities.
3. **Superpose** — every frame is aligned to the reference frame by
   `kabsch_align`, giving the rotation matrix, translation and residual RMSD.
4. **Analyse** — RMSD (global + per-residue Cα), RMSF, radius of gyration,
   hydrogen bonds, PCA, dihedral series, Ramachandran, conformational
   clustering, inter-chain contacts and SASA.
5. **Report** — writes `md_analysis.json` and `md_analysis_frames.csv`.

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `--input` | required | 轨迹文件路径（DCD / PDB / XYZ） |
| `--topology` | — | 拓扑 PDB，提供原子名 / 元素 / 链 / 残基信息 |
| `--output` | `results/` | 输出目录 |
| `--start` / `--stop` / `--stride` | `0` / — / `1` | 帧范围与步长 |
| `--selection` | `protein` | 分析选择串 |
| `--analysis` | `all` | 分析项逗号分隔：`rmsd,rmsf,rg,hbond,pca,dihedral,ramachandran,cluster,interchain,sasa` |
| `--fit-reference` | `0` | RMSD 叠合参考帧 |
| `--hbond-angle` | `120.0` | 氢键最小 D-H···A 角度（度） |
| `--hbond-h-dist` | `2.5` | H···A 最大距离（Å） |
| `--hbond-da-dist` | `3.5` | D···A 最大距离（Å） |
| `--n-clusters` | 自动 | 聚类簇数 |
| `--cluster-threshold` | 自动 | 聚类 RMSD 阈值（Å），给出时优先 |
| `--contact-dist` | `4.5` | 链间接触距离阈值（Å） |
| `--sasa-points` | `100` | SASA 每原子球面采样点数 |
| `--max-cluster-frames` | `200` | 参与聚类的最大帧数 |
| `--format` | `both` | `json` / `csv` / `both` |

`run_simulation.py` 额外参数：`--force-field`（`amber14` / `amber99` / `charmm36`）、
`--solvent`（`none` / `implicit` / `vacuum`）、`--steps`、`--equilibration-steps`、
`--minimize-steps`、`--temperature`、`--friction`、`--timestep`、`--platform`、
`--check-only`。

## Outputs

| File | Content |
|---|---|
| `md_analysis.json` | 完整结果：RMSD/RMSF/Rg/氢键/PCA/二面角/Ramachandran/聚类/链间接触/SASA |
| `md_analysis_frames.csv` | 逐帧汇总表（帧号、时间、RMSD、Rg、氢键数、SASA） |
| `trajectory.dcd` | OpenMM 模拟输出轨迹（`run_simulation.py`） |
| `final_state.pdb` | 模拟末态结构（`run_simulation.py`） |
| `state.csv` | 模拟过程状态报告（能量、温度、体积） |
| `simulation_report.json` | 模拟参数与结果摘要 |

## Validation

- **RMSD must start near zero**: the frame used as `--fit-reference` must have
  RMSD ≈ 0 by construction. If not, the alignment code path is broken.
- **RMSD is bounded by the drift**: a trajectory with a global translation of
  0.02 Å/frame must show a roughly linear RMSD growth after superposition;
  a flat curve means the superposition is absorbing the motion.
- **RMSF sanity**: for a well-folded globular domain, mean Cα RMSF is
  typically 0.5-2 Å; values above 5 Å indicate an unstable or badly
  parameterised simulation. Values below 0.01 Å mean the trajectory is frozen
  or the frames are identical.
- **Rg vs structure**: `Rg ≈ 2.2 × N^0.38` Å is a common empirical rule for
  globular proteins. A sudden jump in Rg usually indicates unfolding.
- **Hydrogen bonds**: a stable protein in water shows 50-200 backbone-sidechain
  and backbone-backbone hydrogen bonds per frame. Zero hydrogen bonds with
  non-zero donors and acceptors means the topology lacks hydrogen atoms
  (X-ray structures) — add them before analysis.
- **PCA**: the first principal component of a globular protein usually
  explains 30-60% of the variance. A single component explaining >95% usually
  means global translation/rotation was not removed before the analysis.
- **Clustering**: the number of clusters must be consistent with the RMSD
  matrix structure. If all frames land in one cluster while the RMSD range is
  wide, the linkage threshold is too large.
- **Ramachandran**: a well-behaved simulation has fewer than 0.5% outliers in
  the allowed regions. More outliers indicates either force-field problems or
  an incorrect topology (e.g. broken backbone connectivity).
- **SASA**: for a 100-residue globular protein, total SASA is typically
  5000-7000 Å². Much lower values indicate collapsed/over-packed structures.
- **DCD parsing**: the parsed `n_atoms` must match the topology atom count;
  a mismatch means the DCD and topology do not correspond.

## References

- `references/methodology.md` — algorithms, formulas and parameter rationale
- `references/troubleshooting.md` — common failures and how to diagnose them
