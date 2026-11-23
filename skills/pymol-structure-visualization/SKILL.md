---
name: pymol-structure-visualization
description: "Use when analyzing inter-chain contacts in protein structures, identifying binding-interface hotspot residues, or building interaction networks and PyMOL visualization scripts. Computes residue-level contact geometry with a uniform-grid spatial index, classifies each contact into hydrogen bonds, salt bridges, hydrophobic, pi-pi stacking and cation-pi interactions, and exports CSV/JSON with Mermaid networks and PyMOL scripts. 用于链间接触分析、热点鉴定与 PyMOL 脚本。"
license: MIT
compatibility: Python 3.9+；核心分析仅需标准库（零第三方依赖）；PyMOL 为可选的渲染端
metadata:
  version: "1.0"
  last-reviewed: "2026-10-26"
  skill-author: sinder
  edam-topic: "topic_0176"
  edam-operation: "operation_0392, operation_0248"
  edam-input-data: "data_1460, data_1743"
  edam-output-data: "data_1547, data_1710"
  edam-input-format: "format_1475"
  edam-output-format: "format_3752, format_3464"


---
## Overview

Computes **inter-chain residue contacts** from PDB files, classified by chemistry not distance alone: a uniform-grid index finds residue pairs by minimum atom distance, then residue chemistry (donor/acceptor, charge, aromaticity, hydrophobicity) sets the type — distance is only a cutoff. Exports CSV / JSON / Mermaid network / PyMOL script — a generated artifact, never a runtime dependency (zero third-party deps, no PyMOL GUI). Ported from legacy `pymol_show/`, 6 defects fixed (§6): the Cartesian product of two `within`-selected residue sets ("each residue came close", not "these two residues touch") and the `pymol.finish_launching(['pymol'])` GUI dependency, leaving it unusable without PyMOL, untestable in CI.

## When to Use This Skill

- What holds two chains together (H-bonds? salt bridges? hydrophobic?)
- Binding-interface hotspot residues (mutagenesis candidates)
- Interaction network graph for figures/slides
- Reproducible script, not GUI selection
- Docking visualisation scripts when PyMOL is **not** installed

**Not applicable**: H-bond geometry → `protein-interactions` (PLIP); docking scoring → `protein-docking-hdock`/`protein-docking-lightdock`; structure quality → `protein-structure`.

## Quick Start

```bash
export PYTHONPATH="$PWD/packages"
python skills/pymol-structure-visualization/scripts/analyze_contacts.py --input complex.pdb --chains A B --output results/
python skills/pymol-structure-visualization/scripts/analyze_contacts.py --input complex.pdb --output results/ --make-mermaid --make-pymol --make-pymol-contacts
python skills/pymol-structure-visualization/scripts/analyze_contacts.py --input complex.pdb --output results/ --max-distance 4.0
python skills/pymol-structure-visualization/scripts/analyze_contacts.py --input complex.pdb --output results/ --exclude-hetero
```

## Workflow

1. **Prepare input** — readable PDB, target chains; multi-MODEL files (NMR, trajectories) use only the **first MODEL**; chain ids from file, not `A`/`B`.
2. **Quality control** — check atom/chain/residue counts; no hydrogens is normal (H-bonds use **residue chemistry**, not H coordinates); waters/ligands by default (`--exclude-hetero` drops them).
3. **Core analysis** — uniform-grid index (cell `max_distance × 1.0001`, 3×3×3 scan) + chemistry classification; each type has its own threshold.
4. **Summarise** — CSV / JSON with per-type counts and hotspot residues.
5. **Validate and plot** — Mermaid network graph and PyMOL script.
6. **Report** — combine with domain knowledge to judge which contacts matter.

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `--input` | required | Input PDB file |
| `--output` | `<input-stem>_contacts/` | Output directory |
| `--chains` | all pairs of chains | Exactly 2 pairs, e.g. `A B` |
| `--max-distance` | `6.0` | Max contact distance (Å) |
| `--exclude-hetero` | off | Ignore HETATM (water, ions, ligands) |
| `--include-hydrogen` | off | Count hydrogens too (off by default) |
| `--include-histidine` | off | Treat HIS as positively charged |
| `--include-ambiguous-donors` | off | Treat ASP/GLU as weak donors |
| `--min-distance` | `0.0` | Only contacts ≥ this apart |
| `--make-mermaid` | off | Generate Mermaid network graph |
| `--make-pymol` | off | Generate PyMOL label script |
| `--make-pymol-contacts` | off | Generate `cmd.distance` script |
| `--json` | `<output-dir>/contacts.json` | Full-result JSON path |
| `-v/--verbose` | off | Emit DEBUG-level logs |
| `-q/--quiet` | off | Only emit warnings and errors |

## Contact Types and Criteria

| Type | Criteria | Max distance |
|---|---|---|
| `Hydrogen_bond` | Donor one side + acceptor the other | 3.5 Å |
| `Salt_bridge` | Charged (LYS/ARG) + negative (ASP/GLU) | 4.0 Å |
| `Hydrophobic` | Both sides hydrophobic | 4.5 Å |
| `Pi_Pi` | Both sides aromatic | 5.0 Å |
| `Cation_Pi` | Positive residue + aromatic residue | 6.0 Å |
| `Ligand_contact` | Either side a non-standard residue (ligand and so on) | `--max-distance` |

A pair can carry several types (ARG–ASP: hydrogen bond *and* salt bridge): `interaction_types` lists all, `interaction_type` the highest-priority (§2).

Three departures from legacy (§3); each **inflated** a count:

1. **ASP/GLU are not hydrogen-bond donors** (carboxyl pKa ≈ 4, fully deprotonated at physiological pH) — `--include-ambiguous-donors` restores them.
2. **HIS is not positively charged by default** (imidazole pKa ≈ 6.0, only ~10% protonated) — `--include-histidine` restores it.
3. **Ligands do not follow the residue-level table** (aromaticity/polarity/charge unknown) — tagged separately as `Ligand_contact`, distance reported as usual; use PLIP for the type.

## Outputs

Under `<output-dir>/`:

- `contacts.csv` — contact detail (UTF-8 BOM, Excel-ready); columns: `Chain_Pair`, `ChainA`, `ResnA`, `ResiA`, `ChainB`, `ResnB`, `ResiB`, `Min_Distance`, `N_Atom_Pairs`, `Interaction_Type` (primary), `All_Types` (semicolon-separated)
- `contact_summary.csv` — per-type counts
- `contacts.json` — full results (contacts + summary + rules)
- `contact_network.mmd` — Mermaid network graph (`--make-mermaid`)
- `visualize_contacts.pml` — PyMOL residue-label script (`--make-pymol`)
- `draw_contacts.pml` — PyMOL dashed contact script (`--make-pymol-contacts`)

## Validation

- **Contact count** — typical protein-protein interface carries 5-30 hydrogen bonds; hundreds ⇒ `--max-distance` too large.
- **Hotspots** — interface hotspots should fall on known active-site or conserved positions.
- **Distances** — hydrogen bonds 2.6-3.2 Å, salt bridges 2.6-3.5 Å; many contacts < 2.0 Å ⇒ steric clashes or precision problem.
- **Atom pair** — re-derive the CSV distance with `cmd.distance` inside PyMOL.
- **Ligand contacts** — a structure containing a ligand ⇒ `by_type` contains `Ligand_contact`.
- **Multi-type** — a pair reported as both salt bridge and hydrogen bond is normal (ARG–ASP guanidinium–carboxylate).

**Note**: supporting analysis only. The thresholds and chemical classification are **general approximations**, not exact force-field criteria; for publication-grade assignments use PLIP.

## References

- [Methodology](references/methodology.md) — criteria (§2), grid algorithm, chemistry tables (§3), the 6 fixed defects
- [FAQ and troubleshooting](references/troubleshooting.md) — symptoms, root causes, solutions
