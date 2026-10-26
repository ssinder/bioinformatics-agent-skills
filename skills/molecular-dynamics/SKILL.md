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
  version: "1.1"
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

# 用 OpenMM 跑一段模拟（生产期）
python scripts/run_simulation.py --input system.pdb --output results/ --steps 100000

# 复合物的分阶段平衡协议（NVT → NPT → 生产期，带位置限制）
python scripts/run_simulation.py --input complex.pdb --output results/ \
    --nvt-steps 50000 --npt-steps 50000 --restraints backbone \
    --monitor-chain-separation --steps 500000

# 生产期也保持弱位置限制（k=10 kJ/mol/Å²）
python scripts/run_simulation.py --input complex.pdb --output results/ \
    --nvt-steps 50000 --restraints backbone --production-restraint-force-constant 10
```

## Workflow

### Analysis path

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

### Simulation path (`run_simulation.py`)

The simulation follows a **staged equilibration protocol**, which is the
standard way to avoid the two classic failure modes of complex simulations:
chains drifting apart during equilibration, and sampling a structure that
never relaxed out of its crystal-packing geometry.

1. **Minimise** — `simulation.minimizeEnergy()` removes the steric clashes
   inherited from the input structure.
2. **NVT equilibration** (`--nvt-steps`) — constant volume, velocities
   re-drawn at the target temperature. Optional position restraints
   (`--restraints`) hold the backbone in place so only side chains and
   solvent relax.
3. **NPT equilibration** (`--npt-steps`) — a `MonteCarloBarostat` is added
   and the simulation is **rebuilt** so the pressure coupling takes effect
   on the density.
4. **Production** — position restraints are removed by default (use
   `--production-restraint-force-constant` to keep a weak version). Chain
   separation is monitored throughout if `--monitor-chain-separation` is set.
5. **Stability check** — the final coordinates are screened for NaN / Inf /
   numerical blow-up; if found, the simulation is rolled back to the input
   structure and the reason is recorded in the report.

> **Implementation note.** Position restraints must be added to the
> `System` **before** the `Simulation` (and therefore the `Context`) is
> created — modifying `system` afterwards does not rebuild the context, so
> the force would silently never take effect. The same applies when adding
> the barostat or removing the restraints: the simulation is rebuilt and
> positions plus velocities are carried over.

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

`run_simulation.py` 额外参数：

| Parameter | Default | Description |
|---|---|---|
| `--force-field` | `amber14` | 力场预设：`amber14` / `amber99` / `charmm36` |
| `--solvent` | `implicit` | 溶剂处理：`none` / `implicit` / `vacuum` |
| `--steps` | `1000` | 生产期步数 |
| `--equilibration-steps` | `100` | NVT/NPT 之后的通用平衡步数 |
| `--minimize-steps` | `200` | 能量最小化步数 |
| `--temperature` | `300.0` | 温度（K） |
| `--friction` | `1.0` | Langevin 摩擦系数 /ps |
| `--timestep` | `0.002` | 积分步长（ps） |
| `--platform` | — | OpenMM 平台名（CUDA / CPU / OpenCL） |
| `--nvt-steps` | `0` | NVT 平衡步数（0 = 跳过） |
| `--npt-steps` | `0` | NPT 平衡步数（0 = 跳过） |
| `--pressure` | `1.0` | NPT 目标压力（bar） |
| `--restraints` | `none` | 平衡期位置限制：`none` / `backbone` / `ca` / `heavy` / `all` |
| `--restraint-force-constant` | `100.0` | 位置限制力常数（kJ/mol/Å²） |
| `--production-restraint-force-constant` | `0.0` | 生产期位置限制力常数（>0 表示保留限制） |
| `--monitor-chain-separation` | off | 监测链间距离，检测复合物解离 |
| `--chain-separation-threshold` | `30.0` | 链间距告警阈值（Å） |
| `--chain-monitor-interval` | `10000` | 链分离检查间隔（步） |
| `--check-only` | off | 仅检查 OpenMM 是否可用并退出 |

## Outputs

| File | Content |
|---|---|
| `md_analysis.json` | 完整结果：RMSD/RMSF/Rg/氢键/PCA/二面角/Ramachandran/聚类/链间接触/SASA |
| `md_analysis_frames.csv` | 逐帧汇总表（帧号、时间、RMSD、Rg、氢键数、SASA） |
| `trajectory.dcd` | OpenMM 模拟输出轨迹（`run_simulation.py`） |
| `final_state.pdb` | 模拟末态结构（`run_simulation.py`） |
| `state.csv` | 模拟过程状态报告（能量、温度、体积） |
| `simulation_report.json` | 模拟参数与结果摘要（含 `stages` / `restraints` / `stability` / `chain_separation`） |

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
- **Equilibration stages**: `simulation_report.json` → `stages` must list every
  stage that was requested. A requested `--nvt-steps` that does not appear
  means the stage was silently skipped.
- **Chain separation**: when `--monitor-chain-separation` is on, the
  `chain_separation` block must exist and `final_check.separated` must be
  `false` for a stable complex. A `true` value means the chains drifted
  apart — check the starting structure (usually a docking pose that was
  never relaxed) before trusting any downstream analysis.
- **Restraint effect**: with `--restraints backbone` on, the NVT trajectory
  must show a noticeably smaller backbone RMSD than an unrestrained run. If
  the numbers are identical, the restraints never took effect — this happens
  when the force is added to the `System` after the `Context` was created.
- **Stability**: `stability.should_rollback` must be `false` for a healthy
  run. If it is `true`, the report's `stability.reason` explains whether the
  coordinates went NaN, hit Inf, or blew up numerically — a NaN in the very
  first frames usually points at a bad force field or missing hydrogens
  rather than at the dynamics themselves.

## References

- `references/methodology.md` — algorithms, formulas and parameter rationale
- `references/troubleshooting.md` — common failures and how to diagnose them
