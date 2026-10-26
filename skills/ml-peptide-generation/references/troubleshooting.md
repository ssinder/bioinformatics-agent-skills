# Troubleshooting: ML Peptide Generation

## "缺少依赖: torch"

**Symptom**: `MissingDependencyError: 缺少依赖: torch`.

**Cause**: the generation script needs `torch` and `transformers`. They are
heavy (the CUDA build alone is several gigabytes), so they are declared as
an *optional* dependency group rather than a hard requirement of the
package.

**Fix**:

```bash
pip install torch transformers --index-url https://download.pytorch.org/whl/cpu
```

CPU-only is fine for small peptide sets. If you only need to **analyse**
existing results, use `analyze_peptides.py` — it uses nothing but the
standard library.

## The script hangs while loading the model

**Symptom**: no output after "正在加载模型 …" for several minutes.

**Cause**: the model is being downloaded from HuggingFace. PepMLM-650M is
roughly 2.5 GB; the first run downloads it.

**Fix**: wait for the download, or pre-fetch it:

```bash
huggingface-cli download TianlaiChen/PepMLM-650M
```

After the first download the model is cached under
`~/.cache/huggingface` and loads in seconds.

**In China**: set the mirror endpoint before running:

```bash
export HF_ENDPOINT=https://hf-mirror.com
```

## Two runs produce different peptides

**Symptom**: the same command twice gives different `peptides.csv`.

**Cause**: either `--seed` was not passed (the default is 42, so this is
unlikely), or the two runs used different `--num-binders` / `--top-ks` /
`--peptide-lengths`. The seed fixes the *random stream*, not the number of
draws — changing how many draws happen changes everything downstream.

**Fix**: keep the full parameter set identical, and check that the `seed`
column in the output is the same in both files.

## A peptide contains `X`

**Symptom**: `binder` shows `X`, or `has_unknown_aa` is `True`.

**Cause**: `X` is in the model's vocabulary as an "unknown residue" token.

**Fix**: by default the skill excludes such tokens from the sampling
candidates, so this should not happen. If it does:

- Check that `--allow-invalid-amino-acids` was **not** passed.
- If you passed it deliberately, the peptides are usable only as a rough
  ranking signal — filter them out with `analyze_peptides.py` (the default
  behaviour) before acting on them.

## The script exits with "超过模型上下文上限 552"

**Symptom**: `蛋白长度 (600) + 肽链长度 (15) = 615 超过模型上下文上限 552`.

**Cause**: PepMLM's context window is 552 tokens; the target plus the
peptide must fit.

**Fix**: use a shorter target fragment. If the binding site is known, take
the relevant structural segment rather than the whole protein. Note that
the check is conservative — it counts residues, and tokenisation can use
more than one token per residue, so a target of 530 residues may still
fail at load time.

## "模型输出里找不到 <mask> 位置"

**Symptom**: `RuntimeError: 模型输出里找不到 <mask> 位置`.

**Cause**: the tokenizer does not recognise `<mask>` — usually a model
whose vocabulary differs from PepMLM's, or an input that was truncated
away entirely.

**Fix**: pass the correct `--model`. If you substituted a different masked
LM, verify that its tokenizer has a mask token:

```python
from transformers import AutoTokenizer
tok = AutoTokenizer.from_pretrained("<你的模型>")
print(tok.mask_token, tok.mask_token_id)
```

## CUDA out of memory

**Symptom**: `torch.OutOfMemoryError`.

**Cause**: the model plus the activation memory for long targets exceeds
the GPU.

**Fix**:
- Use `--device cpu` (slow but no memory limit).
- Shorten the target.
- Close other GPU processes; the model occupies about 2.5 GB in fp32.

## The analysis script reports "筛选后没有任何肽链"

**Symptom**: exit code 1 and `筛选后没有任何肽链（输入 N 条）`.

**Cause**: the filters removed everything — usually `--max-ppl` was set
below the minimum observed PPL, or every peptide contained an invalid
residue.

**Fix**: drop `--max-ppl` and look at the raw distribution first:

```bash
python scripts/analyze_peptides.py --input results/peptides.csv --output-dir analysis/
```

The `parameter_analysis.csv` it writes shows the actual min/median/max PPL
per configuration — pick the threshold from that table, not by guessing.

## The `config` column says `unknown`

**Symptom**: all results end up in one `unknown` group.

**Cause**: the input CSV has neither `peptide_length`/`top_k` columns nor a
filename matching the `{length}_{topk}` pattern.

**Fix**: regenerate with `generate_peptides.py`, which always writes both
columns. For hand-built CSVs, include `peptide_length` and `top_k`
columns, or name the files `results_10_3.csv` style.

## `--help` works but the script fails immediately

**Symptom**: `--help` prints fine, then a real run fails with an import
error.

**Cause**: this is expected behaviour and not a bug — the model and
`transformers` are imported lazily, so argument parsing and validation work
without them. Only the actual generation needs them.

**Fix**: install the dependencies (see the first section).

## The results look identical across seeds

**Symptom**: `--seed 1` and `--seed 42` give the same peptides.

**Cause**: with `top_k=1` the sampling is deterministic — there is only one
candidate, so no randomness is involved. (Note that top-k=1 is not exposed
as a separate flag; a top-k of 1 in `--top-ks` produces this.)

**Fix**: use `--top-ks 3` or higher if you want the seed to matter.

## Peptides are valid but do not bind

**Symptom**: good PPL scores, no experimental activity.

**Cause**: pseudo-perplexity measures plausibility to the language model,
not affinity. The model has no knowledge of the specific interaction.

**Fix**: this is a limitation of the method, not a bug. Use the output as a
**triage list** for experimental follow-up, and validate binding with a
binding assay or an orthogonal computational method (docking, MD, or a
specific predictor). Do not report PPL as evidence of binding.

## Connection errors when downloading the model

**Symptom**: `ConnectionError` or a timeout while fetching from
HuggingFace.

**Fix**: use the mirror (`HF_ENDPOINT=https://hf-mirror.com`) or set a
proxy:

```bash
export HTTPS_PROXY=http://proxy.internal:3128
```

The shared HTTP client honours proxy environment variables — do not
disable that handling, since it would also disable the custom CA bundle
that many corporate networks require.
