---
name: literature-search
description: "Use when mining academic literature from PubMed, OpenAlex or Semantic Scholar — searching by topic or keywords, deduplicating multi-source results, filtering by relevance or citation count, enriching missing abstracts, and exporting RIS, BibTeX or GB/T 7714 citations. 用于多源学术文献检索、去重、相关性筛选、摘要富集与参考文献格式导出的场景。"
license: MIT
compatibility: Python 3.9+；需要 requests（联网检索时）；无第三方库也可对本地 CSV 做去重打分
metadata:
  version: "1.0"
  last-reviewed: "2026-10-25"
  skill-author: sinder
  edam-topic: "topic_3068"
  edam-operation: "operation_0305, operation_0306"
  edam-input-data: "data_0857"
  edam-output-data: "data_0972"
  edam-input-format: "format_3752"
  edam-output-format: "format_3752, format_3464"


---

## Overview

Searches papers by topic or keyword across PubMed, OpenAlex and Semantic Scholar, normalizes the results into one schema, deduplicates them at two levels (exact DOI plus fuzzy title), filters by relevance, citation count and year, fills in missing abstracts from the other sources in turn, and exports citations. Nothing is bound to a specific disease or project, so it serves systematic-review searches in any direction.

## When to Use This Skill

- Mapping the literature landscape on PubMed for a topic (e.g. `antimicrobial peptide`, `protein structure prediction`)
- Broader coverage across several databases, merged and deduplicated
- An existing literature CSV needing rescoring, refiltering and cleaning
- Records without abstracts to complete from OpenAlex / Semantic Scholar
- A reference list to export as RIS / BibTeX / GB/T 7714

## Quick Start

```bash
export PYTHONPATH="$PWD/packages"
python skills/literature-search/scripts/search_literature.py --topic "antimicrobial peptide" --source pubmed --date-range 2020: --max-results 50 --output results/amp_pubmed.csv
python skills/literature-search/scripts/search_literature.py --topic "protein structure prediction" --keywords "AlphaFold" "molecular dynamics" --source pubmed openalex semantic-scholar --date-range 2018:2024 --min-citations 20 --enrich-abstracts --sort-by relevance --output results/prediction.csv
python skills/literature-search/scripts/search_literature.py --input results/raw.csv --keywords "antifungal" --exclude-terms "review" "protocol" --min-relevance 3 --output results/filtered.csv
python skills/literature-search/scripts/export_citations.py --input results/prediction.csv --format bibtex ris gbt7714 --output-dir results/citations
```

## Workflow

1. **Search** — per `--source`: PubMed `esearch.fcgi` → `esummary.fcgi`; OpenAlex `works`; Semantic Scholar `graph/v1/paper/search`.
2. **Normalization** — every record unified into the schema of `literature_results.csv` (see Outputs).
3. **Scoring** — keyword hits (title weight 3.0, abstract weight 1.0) + log-saturated citation count + recency bonus → `relevance`.
4. **Deduplication** — exact DOI first, then `difflib` fuzzy matching on the normalized title (threshold 0.92); dropped records' data is merged into the kept one (abstract, citations, multi-source).
5. **Filtering** — relevance and citation-count thresholds + include/exclude terms + year interval.
6. **Abstract enrichment** (optional) — PubMed (efetch XML) → OpenAlex (rebuilt from the inverted index) → Semantic Scholar (abstract/tldr), in order.
7. **Export** — normalized-schema CSV (UTF-8 BOM) plus optional JSON.

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `--topic` | none | Search topic (primary query) |
| `--keywords` | empty | Extra keywords; AND-narrow search, feed the score |
| `--source` | `pubmed` | Multi-select: `pubmed` / `openalex` / `semantic-scholar` |
| `--date-range` | unlimited | Year interval `START:END`, e.g. `2018:2024` |
| `--max-results` | `50` | Max records per source per query |
| `--input` | none | Existing CSV (offline) |
| `--output` | `<root>/results/literature_results.csv` | Output CSV path |
| `--json-output` | none | Path for the JSON dump |
| `--min-citations` | `0` | Minimum citation count |
| `--min-relevance` | `0` | Minimum relevance score |
| `--min-abstract-length` | `0` | Minimum abstract length |
| `--include-terms` / `--exclude-terms` | empty | Required / excluded keywords |
| `--title-similarity` | `0.92` | Fuzzy title dedup threshold |
| `--enrich-abstracts` | off | Enable abstract enrichment |
| `--sort-by` | `relevance` | `relevance` / `citations` / `year` |
| `--email` | none | Contact email (OpenAlex polite pool) |
| `--ncbi-api-key` / `--s2-api-key` | none | API keys (`NCBI_API_KEY`/`S2_API_KEY` env vars win) |

## Outputs

| File | Description |
|---|---|
| `literature_results.csv` | Results, UTF-8 BOM, columns: `title, authors, year, doi, journal, abstract, citations, source, pmid, url, type, oa_status, relevance, query` |
| `<name>.json` | Full records (incl. parsed fields) |
| `references.bib` | BibTeX, key `@article{surname+year+firstword+index}` |
| `references.ris` | RIS, imports into EndNote / Zotero / Mendeley |
| `references.txt` | GB/T 7714-2015 sequential numbering, or plain-text refs |

## Validation

- **Search coverage**: confirm each source returns the expected record count; a source returning 0 means the year filter is too narrow or a missing API key rate-limited it
- **Deduplication effectiveness**: the count difference before and after dedup should roughly match the true overlap; an `A+B` in `source` means several sources hit it
- **Abstract completeness**: after enrichment, check the null rate of the `abstract` column; OpenAlex abstracts are rebuilt from the inverted index, so records with no upstream abstract stay empty
- **Relevance sanity**: spot-check the top 20 rows by `relevance` against the topic; set `--min-relevance` to 0 first to see the distribution, then tighten
- **Citation format**: import into Zotero/EndNote and check that author names, volume/issue and pages parse correctly; GB/T 7714 output needs manual review
- **Citation count semantics**: `citations` is defined differently per source (PubMed: `pmcrefcount`, OpenAlex: `cited_by_count`), so across-source comparisons are relative only

## References

- `references/methodology.md` — search strategy, per-source API parameters, scoring and dedup
- `references/troubleshooting.md` — common errors (rate limiting, empty results, encoding, duplicates) and debugging
