---
name: peptide-design
description: "Use when analyzing or designing antimicrobial peptides — computing physicochemical properties, classifying peptides as antimicrobial or non-antimicrobial, deduplicating peptide libraries by clustering, and exporting FASTA or CSV. Computes molecular weight, isoelectric point, hydrophobicity and instability index; classifies peptides from GO terms and sequence features; clusters libraries with MMseqs2. 用于多肽性质计算、抗菌肽判别与去冗余。"
license: MIT
compatibility: Python 3.9+；理化性质计算为自包含纯 Python 实现（无需 Biopython/numpy）；MMseqs2 为可选外部依赖
metadata:
  version: "1.0"
  last-reviewed: "2026-10-25"
  skill-author: sinder
  edam-topic: "topic_0203"
  edam-operation: "operation_3646, operation_0291"
  edam-input-data: "data_1254"
  edam-output-data: "data_2016, data_1519, data_1520"
  edam-input-format: "format_1929"
  edam-output-format: "format_1929, format_3752"


---

## Overview

Computes physicochemical properties of peptide sequences from FASTA or CSV with a **self-contained formula implementation** (molecular weight, isoelectric point, hydrophobicity, net charge, instability index, secondary structure propensity, solubility, Boman value), scores candidates as antimicrobial or non-antimicrobial from UniProt GO terms and keywords, and deduplicates the library by clustering. No Biopython needed; formulas in `references/methodology.md`.

## When to Use This Skill

- Physicochemical properties (MW, pI, hydrophobicity, net charge, stability) for a peptide batch as a screening criterion
- Deciding which candidates are antimicrobial and which are not (for a negative sample set)
- Deduplicating a peptide library (MMseqs2 optional; built-in fallback)
- Training data with feature columns for a downstream model
- Exporting results as FASTA or CSV with physicochemical properties

## Quick Start

```bash
export PYTHONPATH="$PWD/packages"

python skills/peptide-design/scripts/analyze_peptides.py --input peptides.fasta --output results/peptides_props.csv --report results/summary.txt
python skills/peptide-design/scripts/analyze_peptides.py --input peptides.csv --sequence-column sequence --classify --antimicrobial-threshold 5 --label antimicrobial --export-fasta results/amp.fasta --output results/amp.csv
python skills/peptide-design/scripts/analyze_peptides.py --input uniprot_peptides.csv --classify --min-quality 60 --require-complete-metadata --output results/high_quality.csv
python skills/peptide-design/scripts/cluster_peptides.py --input peptides.fasta --method auto --min-identity 0.3 --output results/representatives.fasta --output-csv results/cluster_members.csv
```

## Workflow

1. **Read sequences** — FASTA (multi-line, `;` comments, CRLF accepted) or CSV (columns sniffed; without `sequence` → `seq`/`peptide`/`aa_seq`)
2. **Property computation** (`peptide_properties.py`) — every property in the Overview, each self-contained: MW = Σ residue masses + 18.0153, pI by pH 0–14 sweep + 60 bisections, instability index `II = 100/L × Σ DIHEX(di-1, di)` — derivations in `references/methodology.md`
3. **Classification scoring** — annotation layer (5 weighted GO term groups + 5 keyword groups) × 0.6 + physicochemical layer (8 rules) × 0.4
4. **Quality assessment** — UniProt protein existence level + sequence completeness + annotation completeness + physicochemical plausibility
5. **Filtering** — length, unknown residues, metadata completeness, quality score, class label
6. **Clustering** (`cluster_peptides.py`) — MMseqs2 `easy-cluster` preferred, else built-in k-mer Jaccard + LCS identity + average-linkage
7. **Export** — CSV (UTF-8 BOM) / FASTA / JSON / text report

## Parameters

- `--input` (required) — FASTA or CSV input; `--format` (`auto`) sniffs extension/content
- `--sequence-column` (`sequence`), `--id-column`/`--name-column` (`id`/`protein_name`), `--go-column`/`--keyword-column` (`go_terms`/`keywords`) — CSV columns
- `--classify` (off) — enable antimicrobial classification; `--antimicrobial-threshold` (`5.0`), `--non-antimicrobial-threshold` (`2.0`), `--annotation-weight` (`0.6`)
- `--label` (`antimicrobial`) — keep: `antimicrobial` / `non_antimicrobial` / `uncertain` / `all`
- `--min-length`/`--max-length` (`0`, unlimited), `--drop-unknown` (off), `--require-complete-metadata` (off), `--min-quality` (`0`) — filters
- `--export-fasta`/`--report`/`--json-output` (none) — extra FASTA / report / JSON
- Clustering: `--method` (`auto`/`mmseqs`/`internal`), `--min-identity` (`0.3`), `--representative` (`medoid`/`longest`/`first`)

## Outputs

- `peptide_analysis.csv` — 24 physicochemical columns plus `label` / `combined_score` / `quality_score` / `qc_flags`
- `selected.fasta` — sequences that passed the filter (`--export-fasta`)
- `summary.txt` — length/MW/pI/hydrophobicity distributions, stability share, classification distribution, GO coverage, species distribution
- `cluster_members.csv` — `cluster_id` / `representative` / `member_id` / `is_representative` / `cluster_size`
- `representatives.fasta` — one representative per cluster
- `cluster_report.txt` — cluster count, redundancy rate, largest/median size, singletons

## Validation

- **Molecular weight** — a 20 aa peptide lands in the 2000–2500 Da range; each disulfide bridge from C residues removes ~2 Da
- **Isoelectric point** — many Asp/Glu ⇒ pI < 7, many Lys/Arg ⇒ pI > 7; pI falls between pH 1 and 14
- **Hydrophobicity** — Kyte-Doolittle mean > 1.8 strongly hydrophobic, < 0 hydrophilic; amphipathic AMPs sit at 0.5–2.5, hydrophobic moment > 0.4
- **Stability** — instability index < 40 predicts a stable peptide, > 40 rapid degradation (`qc_flags` flags it)
- **Unknown residues** — every metric is less reliable for sequences containing X; `qc_flags` and `has_unknown` flag them
- **Classification** — spot-check each label's samples; `hit_keywords` and `classification_reasons` should be sensible. Pure FASTA input without annotation ⇒ physicochemical layer only, expect more `uncertain`
- **Clustering** — redundancy rate should match sequence diversity; a short library usually exceeds 50%. All-singleton clusters ⇒ `--min-identity` too strict
- **Representatives** — a medoid should look like its cluster's other members; cross-check with `--representative longest`

## References

- `references/methodology.md` — derivation of each property, parameter choices, classification-rule rationale
- `references/troubleshooting.md` — common errors, out-of-range values, debugging
