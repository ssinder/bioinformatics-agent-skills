# Methodology: ML Peptide Generation

## 1. What PepMLM actually does

PepMLM is a **masked language model** for proteins, the same idea as BERT
applied to sequences. Training works like this: take a protein sequence,
hide some residues behind a `<mask>` token, and ask the model to predict
what was there. The prediction is a probability distribution over the
amino-acid vocabulary at each masked position.

Generation inverts that:

```
target sequence:   M V L S P A D K T N V K A A W G K V G A
masked input:      M V L S P A D K T N V K A A W G K V G A <mask> <mask> <mask>
                                                                    ↓
                                            model predicts 3 distributions,
                                            one per masked position
                                                                    ↓
                                            sample one residue from each
                                                                    ↓
candidate peptide:   W D K
```

Each masked position is predicted **independently** — the model does not
condition later positions on earlier samples. This is a deliberate trade-off:
it makes generation fast and every position a draw from the model's marginal
distribution, at the cost of not enforcing peptide-level coherence. The
pseudo-perplexity score below is what catches the incoherence.

## 2. Top-k sampling

At each masked position the model produces a distribution over its whole
vocabulary (tens of thousands of tokens, including `X`, `B`, `Z`,
`<cls>`, `<pad>` and so on). Sampling from the raw distribution would
occasionally pick a special token, which is meaningless as a residue.

**Top-k sampling** restricts the draw to the k most likely tokens, then
samples proportionally within them:

```
logits → top-k (k=3) → softmax over those 3 → sample
```

`k` trades diversity against fidelity:

| k | Effect |
|---|---|
| 1 | Greedy — always the most likely residue. Highly reproducible, low diversity, easily falls into repetition. |
| 3 | The default here. A reasonable balance. |
| 8–10 | More diverse; more low-confidence residues per candidate. |
| > 20 | Approximates unrestricted sampling; the model starts producing dubious residues. |

### Why illegal residues are excluded before sampling

The old implementation sampled from the raw top-k and only *afterwards*
flagged peptides containing non-standard residues (`Has_Unknown_AA`). Two
problems with that:

1. **Wasted budget.** A draw that lands on `X` is a wasted draw. Over 10
   peptides per configuration, a non-trivial fraction is lost.
2. **The residue still influenced the peptide.** The other positions in the
   same peptide were drawn independently, so removing the `X` afterwards
   leaves a peptide that the model never actually proposed as a unit.

This skill masks out the illegal tokens from the candidate set *before* the
softmax, renormalising over what remains. When every candidate at a
position turns out to be illegal (an extreme case), the position degrades
to a uniform distribution and a warning is logged rather than crashing.

### The valid residue set

```
A C D E F G H I K L M N P Q R S T V W Y
```

`X` (unknown), `B` (asparagine/aspartate), `Z` (glutamine/glutamate) and
`J` (leucine/isoleucine) are **ambiguity codes**, not residues. They appear
in the vocabulary but a designed peptide containing them cannot be
synthesised or interpreted.

## 3. Pseudo-perplexity

Perplexity is the standard language-model metric: `exp(mean loss)`, where
the loss is the negative log-likelihood the model assigned to the observed
tokens. It is normally computed over a corpus.

**Pseudo-perplexity** adapts it to the design setting: instead of asking
"how surprised is the model by this sequence?", it asks "how surprised is
the model by this peptide **given the target**?"

The procedure:

1. Concatenate target + candidate peptide, encode to token ids.
2. For each position of the peptide, mask **that position alone** and run
   the model.
3. Collect the cross-entropy loss at each masked position against the true
   token.
4. Average the losses, exponentiate.

```
pseudo_perplexity = exp( (1/L) Σ_i CE(target + peptide, position i) )
```

Interpretation: a PPL of 1.0 means the model assigned the peptide exactly
the probability mass it would have assigned to its own most likely
continuation. Values above 1 mean the peptide is less likely than the
model's best guess; larger values mean more surprising.

### What it does and does not tell you

**Does:** rank candidates by how plausible they are *to this model*. It
filters out the obviously broken peptides — sequences the model considers
implausible in the target's context.

**Does not:** measure binding affinity, stability, solubility, or
immunogenicity. The model was trained on sequence statistics, not on
binding data. A peptide with excellent PPL may still be a poor binder.

**Comparison caveat:** PPL depends on the target sequence. The same
peptide scores differently against different targets. Only compare PPL
values within the same target (and preferably the same length).

### Implementation note

The old code used the model's own `labels=` argument, which internally
computes a mean cross-entropy over non-ignored positions. Different
`transformers` versions differ in how they handle that reduction, so the
result can shift between library versions. This implementation computes
the loss explicitly with `torch.nn.functional.cross_entropy` at each
masked position, which pins the semantics.

## 4. Reproducibility

Sampling is stochastic by nature. Without a seed, two runs on the same
input produce different peptides, which makes every parameter sweep
confounded by sampling noise.

This skill seeds all three relevant sources before any model work:

| Source | What it controls |
|---|---|
| `random` | Python-level randomness |
| `numpy` | array shuffling, dropout masks in some code paths |
| `torch` | the sampling itself (`torch.multinomial`) and model dropout |

The seed is set **before** the model is loaded, because dropout is active
during model construction in some configurations and would otherwise consume
randomness in an order-dependent way.

Every output row carries the seed, so a result file is self-describing: six
months later you can still tell how it was produced.

## 5. Why the target is placed first

The peptide is appended to the **C-terminus** of the target, not the
N-terminus. PepMLM was trained with peptides in this orientation, so
appending at the C-terminus puts the peptide in the position the model
expects. Reversing the order (peptide first) produces markedly worse
scores — the model treats the leading peptide as the "context" and the
trailing target as the thing to predict.

## 6. Context length

PepMLM supports up to **552 tokens**. The target and the peptide must fit
together:

```
len(target) + len(peptide) ≤ 552
```

Because tokenisation is not one-character-per-token (some amino acids map
to multiple tokens), the check is conservative — it uses residue counts as
an upper bound, which can only reject borderline cases, never accept
over-long ones.

Older implementations silently truncated the sequence instead, which is
worse than failing: the model then generates peptides for a *different
protein* than intended, and nothing indicates a problem. This skill raises
an explicit error.

## 7. Parameter sweeps

The natural use of this skill is to find which settings produce the best
candidates. Sweep two parameters:

- **peptide length** (10 / 12 / 15 aa) — longer peptides have more
  interaction surface but lower hit rates and worse synthesability.
- **top-k** (3 / 5 / 8) — higher k explores more but produces more noise.

Sweep them **together**, not one at a time: length and top-k interact.
A configuration that looks good at k=3 may be the worst at k=8.

`analyze_peptides.py` groups results by the `config` column
(`{length}aa_topk{k}`) and reports the mean, median, min and max PPL per
group. **Use the median, not the mean**, when a group contains outliers —
a single catastrophic peptide can drag the mean noticeably, and the median
is what describes the typical candidate.

## 8. Filtering

The default filter removes peptides containing non-standard residues. On top
of that:

- `--max-ppl` sets a ceiling. Choose it by looking at the distribution
  first (the parameter-analysis table gives min/median/max per group), not
  by guessing a round number.
- `--min-length` excludes peptides too short to bind meaningfully. Below
  about 8 residues there is usually not enough surface for a specific
  interaction.

## References

- PepMLM model card: <https://huggingface.co/TianlaiChen/PepMLM-650M>
- Devlin et al., *BERT* (the masked LM formulation): <https://arxiv.org/abs/1810.04805>
- Pseudo-perplexity as a sequence-scoring measure: Dekker et al., *ESP*
  <https://doi.org/10.1093/nar/gkab446>
