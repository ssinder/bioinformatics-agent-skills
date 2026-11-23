---
name: protein-interactions
description: "Use when interaction detection, hotspot ranking, PLIP report parsing, interaction network visualisation or AlphaFold confidence analysis is requested. Parses PLIP XML reports and AlphaFold packages with a standard-library XML parser; extracts residue-level salt bridges, hydrogen bonds, hydrophobic, pi-stacking, water bridges, halogen bonds and metal complexes, ranks hotspot residues, and emits Mermaid graphs plus PyMOL scripts. 用于相互作用解析与热点鉴定。"
license: MIT
compatibility: Python 3.9+; optional PLIP, MDAnalysis, requests
metadata:
  version: "1.0"
  last-reviewed: "2026-10-25"
  skill-author: sinder
  edam-topic: "topic_0176"
  edam-operation: "operation_0248"
  edam-input-data: "data_0906, data_2537"
  edam-output-data: "data_1566, data_1547"
  edam-input-format: "format_2054, format_1475"
  edam-output-format: "format_3752, format_3464"

# Protein Interaction Analysis


---

## Overview

Turns PLIP (Protein-Ligand Interaction Profiler) reports into structured, residue-level interaction records and visualisations. The XML parsing uses the Python standard library (`xml.etree.ElementTree`) only — no `lxml`, no Biopython, no network access — so a report produced elsewhere can be analysed anywhere. On top of the parsed records it computes per-type and per-residue statistics, renders an interaction network as a Mermaid flowchart with per-type colouring, and generates a PyMOL script that selects and labels interacting residues in 3D.

AlphaFold model packages (`*.zip` from the AlphaFold Protein Structure Database) can be fed in directly: the best model is picked by predicted TM-score confidence (`pTM + 0.5 × ipTM`) and the per-model confidence JSON summarised.

## When to Use This Skill

- Parsing a PLIP `report.xml` into machine-readable residue-level interactions
- Ranking hotspot residues by interaction count on receptor and ligand sides
- Detecting salt bridges / hydrogen bonds / π-stacking / water bridges without PLIP
- Selecting the best AlphaFold model from a multi-model package and reporting pLDDT
- Generating a Mermaid interaction graph for reports, slides or documentation
- Generating PyMOL commands that visualise and label interacting residues

## Quick Start

```bash
python scripts/analyze_interactions.py --input report.xml --output results/
python scripts/analyze_interactions.py --input report.xml \
    --output results/ --make-mermaid --make-pymol

# Run PLIP on the structure directly (requires plip to be installed)
python scripts/analyze_interactions.py --input 1crn.pdb --run-plip \
    --plip-structure 1crn.pdb --output results/

# AlphaFold model package (the most confident model is picked automatically)
python scripts/analyze_interactions.py --input AF-P69905-F1-model_v4.zip \
    --output results/ --make-mermaid
```

## Workflow

1. **Collect** — `--input` accepts a PLIP `report.xml`, a directory containing `report.xml` files, or an AlphaFold `*.zip` package; directories are scanned recursively.
2. **Parse** — `plip_parser.py` walks the XML with `xml.etree.ElementTree`, handling both the modern attribute-style schema (`<site ... onf="A:ARG+102:NE">` with child `<salt_bridge>` elements) and the legacy element-style schema (`<resnr>`, `<restype>`, `<reschain>` with `_lig` suffixed variants). Labels decode into `chain / resname / resnr / atom` fields, separating protein residues (`A:ARG+102:NE`) from ligand atoms (`LIG:1:CL1`).
3. **Summarise** — `summarize_interactions` counts interactions per type; the entry script also ranks receptor and ligand residues by participation count for hotspots.
4. **Visualise** — `interaction_network.py` renders the Mermaid flowchart (nodes = residues in receptor / ligand subgraphs, edges = interactions coloured by type via `linkStyle`) and the PyMOL script (`select` / `show sticks` / `color` / `label` / `distance`).
5. **Report** — writes `interaction_analysis.json` and `interactions_residue_level.csv` (one row per interaction).

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `--input` | required | `report.xml`, a directory of reports, or an AlphaFold `*.zip` |
| `--output` | `results/` | Output directory |
| `--run-plip` | off | Run PLIP on the input structure |
| `--plip-structure` | — | Structure file handed to PLIP (PDB) |
| `--plip-executable` | `plip` found automatically | Path to the PLIP executable |
| `--plip-timeout` | `600` | PLIP timeout in seconds |
| `--make-mermaid` | off | Generate the Mermaid interaction network graph |
| `--make-pymol` | off | Generate the PyMOL annotation script |
| `--format` | `both` | `json` / `csv` / `both` |
| `--top-residues` | `10` | Number of hotspot residues to show in the statistics |

## Outputs

| File | Content |
|---|---|
| `interaction_analysis.json` | Complete results per report: interaction records, per-type statistics, hotspot residues |
| `interactions_residue_level.csv` | Residue-level flat table (receptor residue, ligand residue, type, atom, distance) |
| `*_interaction_network.md` | Mermaid `flowchart` network graph (with `--make-mermaid`) |
| `*_annotation.pml` | PyMOL annotation script (with `--make-pymol`) |

## Validation

- **Record count sanity**: a PLIP report on a ligand with a buried binding site normally yields tens to a few hundred interactions; `0` records for a non-empty `report.xml` means the schema was unrecognised.
- **Residue label integrity**: decoded labels must keep `resnr` as an integer and `resname` as a three-letter code; protein labels carry a `+` separator (`A:ARG+102`), ligand labels do not.
- **Hotspot plausibility**: a binding hotspot usually shows ≥ 3 interactions (e.g. `ARG102` engaging the ligand in 3 distinct contacts). A flat histogram of one interaction per residue usually means the receptor/ligand roles were swapped or hydrogens stripped.
- **Type cross-check**: salt bridges need oppositely charged side chains (ARG/LYS versus ASP/GLU); a salt bridge between two ARG residues is a parsing error.
- **AlphaFold confidence**: `pLDDT > 90` is high confidence, `70-90` is acceptable, `< 70` marks a disordered or unreliable region; a mean pLDDT below 70 means the model should not be used for interaction analysis.
- **Network graph integrity**: the Mermaid file must contain one `linkStyle` line per edge; a mismatch means edges were dropped.

## References

- `references/methodology.md` — PLIP XML schema, interaction definitions, statistics, rendering
- `references/troubleshooting.md` — common failures and how to diagnose them
