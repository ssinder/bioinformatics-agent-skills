---
name: protein-structure
description: "Use when secondary structure assignment, radius of gyration, shape analysis, DSSP, PDB parsing or protein structure annotation is requested. Parses ATOM/HETATM records with a dependency-free pure-Python parser, computes radius of gyration, gyration tensor and dipole moments, assigns secondary structure via DSSP or a built-in Kabsch-Sander implementation, and optionally retrieves UniProt metadata. 用于 PDB 解析与二级结构指派。"
license: MIT
compatibility: Python 3.9+; optional Biopython, DSSP (mkdssp), requests
metadata:
  version: "1.0"
  last-reviewed: "2026-10-25"
  skill-author: sinder
  edam-topic: "topic_0176"
  edam-operation: "operation_2406, operation_0320"
  edam-input-data: "data_1460, data_2537"
  edam-output-data: "data_2956, data_0897"
  edam-input-format: "format_1475"
  edam-output-format: "format_3752, format_3464"

# Protein Structure Analysis


---

## Overview

Turns raw PDB coordinates into quantitative structural descriptors without heavyweight dependencies. It ships a self-contained PDB parser that handles `ATOM` / `HETATM` records, multi-model files, alternate locations and insertion codes, plus a built-in secondary-structure assigner based on the Kabsch-Sander hydrogen-bond energy model. Geometry (radius of gyration, gyration tensor principal axes, shape descriptors, dipole moment) is computed in pure Python, and UniProt metadata can be attached on request.

## When to Use This Skill

- Secondary structure assignment of a protein (`DSSP`, `H`/`E`/`C` content, helix percentage)
- Radius of gyration, gyration tensor, shape anisotropy of a folded protein
- Parsing PDB files without installing Biopython
- Batch statistics over a directory of PDB files (length distribution, average Rg)
- Attaching UniProt protein names, organism and function to structural results
- Preliminary characterisation of AlphaFold models or predicted structures

## Quick Start

```bash
# Single structure: parsing + geometry + secondary structure
python scripts/analyze_structure.py --input 1crn.pdb --output results/

# Force the built-in DSSP implementation (when mkdssp is not installed)
python scripts/analyze_structure.py --input 1crn.pdb --output results/ --dssp-method builtin

# Process a whole directory in batch (4 threads)
python scripts/batch_analyze.py --input ./pdb_data --output results/ --workers 4

# Attach UniProt metadata (requires network access)
python scripts/analyze_structure.py --input P69905.pdb --output results/ --query-uniprot
```

## Workflow

1. **Parse** — `pdb_parser.py` reads fixed-column PDB records into a `Structure → Model → Chain → Residue → Atom` tree. Alternate locations keep only the primary conformer (`' '`, `'A'`, `'1'`); `MODEL`/`ENDMDL` splits multi-model files.
2. **Summarise** — chain/residue/atom counts, element composition, heteroatom inventory, per-chain one-letter sequences.
3. **Geometry** — `geometry.py` computes mass-weighted `Rg` (and `Rg` over Cα only), the gyration tensor with a Jacobi eigen-decomposition, shape descriptors (asphericity, acylindricity, relative shape anisotropy) and a partial-charge dipole estimate.
4. **Secondary structure** — `secondary_structure.py` either shells out to `mkdssp` or runs the built-in assigner: rebuild missing amide hydrogens, evaluate Kabsch-Sander hydrogen-bond energies, and classify `H`/`G`/`I`/`E`/`B`/`T`/`S`/`-`.
5. **Annotate** — optionally query `rest.uniprot.org` for protein name, organism, sequence length and function via the shared `http_session` (retry + rate limit).
6. **Report** — write `structure_analysis.json` and `structure_analysis.csv`.

## Parameters

| Parameter | Script | Default | Description |
|---|---|---|---|
| `--input` | both | required | PDB file or directory of PDB files |
| `--output` | both | `results/` | Output directory for reports |
| `--dssp-method` | both | `auto` | `auto` / `external` / `builtin` secondary-structure backend |
| `--dssp-executable` | both | auto-detected | Path to `mkdssp` or `dssp` |
| `--hbond-energy` | both | `-0.5` | Hydrogen-bond energy threshold (kcal/mol) for the built-in assigner |
| `--model` | analyze | `1` | Which model to analyse (1-based) |
| `--query-uniprot` | both | off | Enable UniProt metadata lookup (network access) |
| `--format` | analyze | `both` | `json` / `csv` / `both` |
| `--workers` | batch | `4` | Worker threads for batch processing |
| `--pattern` | batch | `*.pdb` | Glob pattern for batch discovery |
| `--recursive` | batch | off | Recurse into subdirectories |
| `--min-residues` | batch | `0` | Skip structures shorter than this |
| `--max-files` | both | `0` | Cap the number of files (0 = unlimited) |

## Outputs

| File | Content |
|---|---|
| `structure_analysis.json` | Full nested report: structure summary, chains, geometry, per-residue secondary structure, UniProt annotation |
| `structure_analysis.csv` | Flat one-row-per-structure table (Rg, shape class, helix/sheet percentages) |
| `batch_summary.json` | Aggregate statistics over a directory plus per-structure detail |
| `batch_results.csv` | Flat table of the batch run |

## Validation

- **Hydrogen-bond geometry**: in a genuine alpha helix the `O(i)···N(i+4)` distance must be ~2.9-3.1 Å. The built-in assigner only reports `H` when that contact exists.
- **Backbone sanity**: `Cα(i)-Cα(i+1)` should be ~3.8 Å and `C(i)-N(i+1)` ~1.33 Å. Distances far from these indicate a malformed or non-standard PDB.
- **Helix content**: globular proteins typically show 20-50% helix. A value near 0% or 100% for a folded protein is a red flag for the hydrogen-bond threshold.
- **Cross-check**: compare `builtin` against `external` DSSP on the same file (`--dssp-method builtin` vs `--dssp-method external`). Agreement within a few percent on helix content is expected; large divergence usually means the structure lacks resolved hydrogens or is highly non-standard.
- **Rg plausibility**: for a 100-residue globular domain `Rg` is roughly 11-14 Å (`Rg ≈ 2.2 × N^0.38` Å is a common empirical rule).
- **Element counts**: hydrogen atoms are usually absent from X-ray structures, so `element_counts` containing only C/N/O is expected, not a parse failure.

## References

- `references/methodology.md` — algorithms, formulas and parameter selection rationale
- `references/troubleshooting.md` — common failures and how to diagnose them
