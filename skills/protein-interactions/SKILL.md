---
name: protein-interactions
description: >-
  Analyzes protein-ligand and protein-protein interactions from PLIP XML reports
  and AlphaFold model packages using a self-contained standard-library XML
  parser: extracts residue-level salt bridges, hydrogen bonds, hydrophobic,
  pi-stacking, water bridges, halogen bonds and metal complexes, ranks hotspot
  residues, and emits Mermaid network graphs plus PyMOL annotation scripts. This
  skill should be used when interaction detection, hotspot ranking, PLIP report
  parsing, interaction network visualisation or AlphaFold confidence analysis
  is requested.
license: MIT
compatibility: Python 3.9+; optional PLIP, MDAnalysis, requests
metadata:
  version: "1.0"
  last-reviewed: "2026-10-25"
  skill-author: sinder
---

# Protein Interaction Analysis

## Overview

This skill turns PLIP (Protein-Ligand Interaction Profiler) reports into structured,
residue-level interaction records and visualisations. The XML parsing is implemented
with the Python standard library (`xml.etree.ElementTree`) only — no `lxml`, no
Biopython, no network access — so a report produced on a workstation can be analysed
anywhere. On top of the parsed records the skill computes per-type and per-residue
statistics, renders an interaction network as a Mermaid flowchart with
per-interaction-type colouring, and generates a ready-to-run PyMOL script that
selects and labels the interacting residues in the 3D structure.

AlphaFold model packages (`*.zip` from the AlphaFold Protein Structure Database) can
be fed in directly: the best model is picked by predicted TM-score confidence
(`pTM + 0.5 × ipTM`) and the per-model confidence JSON is summarised.

## When to Use This Skill

- Parsing a PLIP `report.xml` into machine-readable residue-level interactions
- Ranking hotspot residues by interaction count on receptor and ligand sides
- Detecting salt bridges / hydrogen bonds / π-stacking / water bridges without PLIP
- Selecting the best AlphaFold model from a multi-model package and reporting pLDDT
- Generating a Mermaid interaction graph for reports, slides or documentation
- Generating PyMOL commands that visualise and label the interacting residues

## Quick Start

```bash
# 解析已有 PLIP 报告
python scripts/analyze_interactions.py --input report.xml --output results/

# 同时生成 Mermaid 网络图与 PyMOL 标注脚本
python scripts/analyze_interactions.py --input report.xml \
    --output results/ --make-mermaid --make-pymol

# 直接运行 PLIP 分析结构（需已安装 plip）
python scripts/analyze_interactions.py --input 1crn.pdb --run-plip \
    --plip-structure 1crn.pdb --output results/

# 分析 AlphaFold 模型包（自动挑选置信度最高的模型）
python scripts/analyze_interactions.py --input AF-P69905-F1-model_v4.zip \
    --output results/ --make-mermaid
```

## Workflow

1. **Collect** — `--input` accepts a PLIP `report.xml`, a directory containing
   `report.xml` files, or an AlphaFold `*.zip` package. Directories are scanned
   recursively for `report.xml`.
2. **Parse** — `plip_parser.py` walks the XML with `xml.etree.ElementTree`,
   handling both the modern attribute-style schema (`<site ... onf="A:ARG+102:NE">`
   with child `<salt_bridge>` elements) and the legacy element-style schema
   (`<resnr>`, `<restype>`, `<reschain>` with `_lig` suffixed variants).
   Residue labels are decoded into `chain / resname / resnr / atom` fields,
   distinguishing protein residues (`A:ARG+102:NE`) from ligand atoms
   (`LIG:1:CL1`).
3. **Summarise** — `summarize_interactions` counts interactions per type; the
   entry script additionally ranks receptor and ligand residues by participation
   count to identify hotspots.
4. **Visualise** — `interaction_network.py` renders the Mermaid flowchart
   (nodes = residues grouped into receptor / ligand subgraphs, edges = individual
   interactions coloured by type via `linkStyle`) and the PyMOL script
   (`select` / `show sticks` / `color` / `label` / `distance` commands).
5. **Report** — writes `interaction_analysis.json` and
   `interactions_residue_level.csv`, one row per interaction.

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `--input` | required | `report.xml`、包含报告的目录，或 AlphaFold `*.zip` |
| `--output` | `results/` | 输出目录 |
| `--run-plip` | off | 对输入结构运行 PLIP |
| `--plip-structure` | — | 传给 PLIP 的结构文件（PDB） |
| `--plip-executable` | 自动查找 `plip` | PLIP 可执行文件路径 |
| `--plip-timeout` | `600` | PLIP 超时秒数 |
| `--make-mermaid` | off | 生成 Mermaid 相互作用网络图 |
| `--make-pymol` | off | 生成 PyMOL 标注脚本 |
| `--format` | `both` | `json` / `csv` / `both` |
| `--top-residues` | `10` | 统计中展示的热门残基数量 |

## Outputs

| File | Content |
|---|---|
| `interaction_analysis.json` | 每个报告的完整结果：相互作用记录、按类型统计、热门残基 |
| `interactions_residue_level.csv` | 残基级扁平表（受体残基、配体残基、类型、原子、距离） |
| `*_interaction_network.md` | Mermaid `flowchart` 网络图（启用 `--make-mermaid`） |
| `*_annotation.pml` | PyMOL 标注脚本（启用 `--make-pymol`） |

## Validation

- **Record count sanity**: a PLIP report on a ligand with a buried binding site
  normally yields tens to a few hundred interactions; `0` records for a
  non-empty `report.xml` means the schema was not recognised — inspect the file.
- **Residue label integrity**: decoded labels must keep `resnr` as an integer and
  `resname` as a three-letter code; protein labels contain a `+` separator
  (`A:ARG+102`), ligand labels do not.
- **Hotspot plausibility**: a binding hotspot usually shows ≥ 3 interactions
  (e.g. `ARG102` engaging the ligand in 3 distinct contacts). A flat histogram
  of one interaction per residue usually means the receptor/ligand roles were
  swapped or hydrogens were stripped from the input structure.
- **Type cross-check**: salt bridges require oppositely charged side chains
  (ARG/LYS versus ASP/GLU); a salt bridge between two ARG residues is a parsing
  error, not a result.
- **AlphaFold confidence**: `pLDDT > 90` is high confidence, `70-90` is
  acceptable, `< 70` marks a disordered or unreliable region; a mean pLDDT below
  70 means the model should not be used for interaction analysis at all.
- **Network graph integrity**: the Mermaid file must contain one `linkStyle`
  line per edge; a mismatch means edges were dropped during rendering.

## References

- `references/methodology.md` — PLIP XML schema, interaction definitions, statistics and rendering
- `references/troubleshooting.md` — common failures and how to diagnose them
