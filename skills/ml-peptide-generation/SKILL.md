---
name: ml-peptide-generation
description: >-
  Generates candidate binder peptides for a target protein using the PepMLM-650M
  masked language model — appends <mask> tokens to the target sequence and
  samples amino acids at each position, with illegal residues (X, B, Z)
  excluded from the sampling candidates. Reports pseudo-perplexity for each
  candidate, supports reproducible generation via an explicit random seed,
  and provides offline filtering and parameter-sweep analysis that needs no GPU.
  This skill should be used when designing peptides or protein fragments that
  should bind a target, generating a ranked candidate pool for a target
  sequence, or comparing how peptide length and top-k sampling settings affect
  candidate quality.
license: MIT
compatibility: Python 3.9+; generation needs torch + transformers; analysis needs only the standard library
metadata:
  version: "1.0"
  last-reviewed: "2026-10-26"
  skill-author: sinder
---

# ML Peptide Generation

## Overview

This skill generates candidate binder peptides for a target protein with
**PepMLM-650M**, a masked language model trained on protein sequences.

The mechanism is straightforward: append N `<mask>` tokens to the target
sequence, ask the model for its amino-acid distribution at each masked
position, and sample from it. Repeating that for N positions produces an
N-residue peptide.

Two things this skill does that a naive implementation does not:

- **Reproducible by construction.** Every run takes an explicit `--seed`.
  The same seed and the same input always produce the same peptides, so a
  parameter sweep is a controlled experiment rather than a lottery.
- **Illegal residues are excluded at sampling time, not afterwards.** The
  model's vocabulary contains `X` (unknown), `B` and `Z` (ambiguous amino
  acids). Sampling from the raw vocabulary produces peptides that are not
  chemically meaningful. Here those candidates are removed *before* the
  draw, so the sampling budget is not wasted on rejects.

## When to Use This Skill

- Designing a peptide or protein fragment intended to bind a target protein
- Producing a ranked candidate pool for a set of targets
- Comparing peptide lengths (10 / 12 / 15 aa) or top-k settings
- Screening candidates before committing to synthesis or wet-lab validation

## Quick Start

```bash
export PYTHONPATH="$PWD/packages"

# 生成（需要 torch + transformers）
python scripts/generate_peptides.py \
    --input-fasta targets.fasta \
    --output-dir results/ \
    --peptide-lengths 10,12,15 \
    --num-binders 20 \
    --seed 42

# 分析（只用标准库，无需 GPU）
python scripts/analyze_peptides.py \
    --input results/peptides.csv \
    --output-dir analysis/ \
    --max-ppl 5.0
```

## Workflow

1. **Collect targets** — `--input-fasta` (FASTA), `--input-csv` (needs `id`
   and `seq` columns) or `--input-seq` (a single sequence).
2. **Set the seed** — `--seed` fixes `random`, `numpy` and `torch`. The seed
   is written into every output row so a run can be reproduced later.
3. **Generate** — for each target × peptide length × top-k combination,
   sample `--num-binders` peptides and score each by pseudo-perplexity.
4. **Analyse** — `analyze_peptides.py` filters invalid peptides, applies
   optional PPL/length thresholds, and reports per-configuration statistics.
5. **Select** — rank by pseudo-perplexity and carry the top candidates into
   downstream validation.

## Parameters

### `generate_peptides.py`

| Parameter | Default | Description |
|---|---|---|
| `--input-fasta` | — | Target protein FASTA file |
| `--input-csv` | — | Target CSV file (needs `id` and `seq` columns) |
| `--input-seq` | — | A single target sequence |
| `--output-dir` | `peptide_results` | Output directory |
| `--per-target-dir` | off | Also write one CSV per target |
| `--peptide-lengths` | `10,12,15` | Peptide lengths to generate |
| `--top-ks` | `3` | Top-k values to sweep |
| `--num-binders` | `10` | Peptides per (target × length × top-k) |
| `--seed` | `42` | Random seed — same seed, same results |
| `--model` | `TianlaiChen/PepMLM-650M` | Pretrained model |
| `--device` | auto | Force `cuda` or `cpu` |
| `--allow-invalid-amino-acids` | off | Allow sampling `X`/`B`/`Z` |

### `analyze_peptides.py`

| Parameter | Default | Description |
|---|---|---|
| `--input` | required | One or more result CSV files |
| `--output-dir` | `peptide_analysis` | Output directory |
| `--max-ppl` | none | Keep only peptides with PPL ≤ this value |
| `--min-length` / `--max-length` | none | Peptide length bounds |
| `--keep-invalid` | off | Keep peptides containing non-standard residues |
| `--top-n` | `20` | How many to print |

## Outputs

| File | Content |
|---|---|
| `peptides.csv` | All candidates: rank, target, length, top-k, binder, pseudo-perplexity, validity flag, seed |
| `peptides.json` | Same content as JSON |
| `filtered_peptides.csv` | Candidates passing the filters, sorted by PPL |
| `parameter_analysis.csv` / `.json` | Per-configuration totals, mean/median/min/max PPL |

## Validation

- **Pseudo-perplexity measures the model's own belief, not binding affinity.**
  A low PPL means "this looks like a peptide the model would have generated",
  not "this binds the target". It is a triage filter for wet-lab follow-up,
  nothing more. Any claim about affinity needs an experimental or
  orthogonal computational check.
- **PPL is only comparable within one (target, length) pair.** The scale
  depends on the target sequence; do not compare PPL across different
  targets and conclude one target is "easier".
- **Check the validity flag.** `has_unknown_aa=True` means the peptide
  contains a residue the model could not decide on — such a peptide is not
  a real design. The skill excludes these at sampling time, but the flag
  stays in the output so the filter can be audited.
- **Verify reproducibility.** Re-run with the same `--seed` and diff
  `peptides.csv`; the files must be byte-identical apart from the
  `generated_at` timestamp.
- **Peptides shorter than about 8 residues rarely bind meaningfully.** Use
  `--min-length` to exclude them from a synthesis list.
- **The model truncates at 552 tokens.** Target sequences longer than that
  are rejected by the skill rather than silently truncated — a truncated
  target produces peptides for the wrong protein.

## References

- `references/methodology.md` — PepMLM mechanism, pseudo-perplexity derivation, sampling strategy
- `references/troubleshooting.md` — common failures and how to diagnose them
- Model card: <https://huggingface.co/TianlaiChen/PepMLM-650M>
- E-utilities context for target discovery: the `bioproject-search` skill
