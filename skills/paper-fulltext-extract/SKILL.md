---
name: paper-fulltext-extract
description: "Use when extracting structured full text from scientific PDF papers, parsing section trees and figures, building a searchable corpus across many papers, or pulling lead sentences for rapid paper triage. Converts PDFs to TEI XML via a local GROBID service and parses titles, authors, abstracts, keywords, sections, figures and references into JSON, Markdown or HTML. 用于 PDF 全文抽取与语料检索。"
license: MIT
compatibility: Python 3.9+；仅需标准库（GROBID 服务本身需用 Docker 单独部署）
metadata:
  version: "1.0"
  last-reviewed: "2026-10-26"
  skill-author: sinder
  edam-topic: "topic_3068"
  edam-operation: "operation_0306"
  edam-input-data: "data_0972"
  edam-output-data: "data_0972"
  edam-input-format: "format_2021"
  edam-output-format: "format_2021, format_3464"


---
## Overview

Turns PDF batches into searchable, citable full text. GROBID is a locally deployed Java service (not a cloud API) restoring PDF layout into TEI XML (hierarchy, figures, references); parsing, conversion, indexing, search follow. Ported from legacy `paper_extract/`; fixes 5 defects: keyword whitespace splitting, misassigned titles, duplicate section trees, over-counted references, unreadable `xml:id`.

## When to Use This Skill

- Structured metadata (title/authors/abstract/keywords/DOI/journal)
- Convert to Markdown/HTML for reading, translation or LLM input
- Cross-paper full-text search ("which papers mention off-target")
- Skim **lead sentences of the introduction**; extract **reference lists** for citations

## Quick Start

```bash
export PYTHONPATH="$PWD/packages"

# Deploy GROBID (one-time; downloads models, wait for "Server started")
docker run -it -p 8070:8070 kermitt/grobid

python skills/paper-fulltext-extract/scripts/extract_fulltext.py --pdf-dir ./pdfs --output ./out --formats json,markdown,html --build-index
python skills/paper-fulltext-extract/scripts/extract_fulltext.py --tei-dir ./tei --output ./out --intro-only
python skills/paper-fulltext-extract/scripts/tei_parser.py --input ./tei/paper.tei.xml --json out/paper.json --markdown out/paper.md --html out/paper.html
python skills/paper-fulltext-extract/scripts/grobid_client.py --input ./pdfs --output ./tei --check-only
python skills/paper-fulltext-extract/scripts/search_corpus.py --xml-dir ./tei --keyword "CRISPR" --scope all --export-csv out/corpus_metadata.csv
```

Cmd 2 = intro lead sentences (triage); 3 = parse TEI; 4 = GROBID batch (resumable); 5 = search; chain: `grobid_client` → `tei_parser` → `search_corpus` → `extract_first_sentences`.

## Workflow

1. **PDF → TEI** (needs GROBID) — `POST /api/processFulltextDocument`, field `input`; `flavor` `precision`/`speed` (accurate/fast); batch **resumes** (skips existing `<name>.tei.xml`).
2. **TEI → data** — metadata (title, authors, affiliations, abstract, keywords, DOI, year, journal); section tree (`<div>` nests, hierarchy kept); figures/tables (number, caption, path); references (`<back>` only).
3. **Output formats** — JSON (machine) / Markdown (reading, LLM) / HTML (captions).
4. **Indexing + search** — scopes: section title / body paragraph / paper title; metadata as CSV.
5. **Lead sentences** — first sentence of each introduction paragraph, as CSV.

## Parameters

`extract_fulltext.py`:

|Parameter|Default|Description|
|---|---|---|
|`--pdf-dir`|none|PDF dir (GROBID)|
|`--tei-dir`|none|TEI dir (no GROBID)|
|`--output`|none|Out dir (required)|
|`--grobid-server`|`http://127.0.0.1:8070`|GROBID URL|
|`--flavor`|`precision`|`precision`/`speed`|
|`--no-consolidate-header`|off|Disable `consolidateHeader` (default on)|
|`--no-consolidate-citations`|off|Disable `consolidateCitations` (default on)|
|`--pause`|`0`|Delay (s)|
|`--formats`|`json,markdown`|`json`/`markdown`/`html`|
|`--image-base`|`images`|HTML img prefix|
|`--intro-only`|off|Intro lead sents|
|`--intro-section`|`introduction`|Intro section|
|`--max-intro-sentences`|`3`|Max sentences|
|`--build-index`|off|Build index|
|`--search`|none|Search term (needs index)|
|`--search-scope`|`paragraph`|`paragraph`/`section`/`title`/`all`|
|`--skip-grobid-check`|off|Skip GROBID check|

`tei_parser.py`:

|Parameter|Default|Description|
|---|---|---|
|`--input`|none|TEI XML (required)|
|`--json` / `--markdown` / `--html`|none|Output paths|
|`--image-base`|empty|HTML img prefix|

`grobid_client.py`:

|Parameter|Default|Description|
|---|---|---|
|`--input` / `--output`|none|PDF dir / TEI dir|
|`--server`|`http://127.0.0.1:8070`|GROBID URL|
|`--timeout`|`300`|Timeout (s)|
|`--max-retries`|`3`|Attempts incl 1st|
|`--retry-delay`|`5`|Delay before retry|
|`--flavor`|`precision`|`precision` / `speed`|
|`--pause`|`0`|Delay (s)|
|`--check-only`|off|Check reach|

`search_corpus.py`:

|Parameter|Default|Description|
|---|---|---|
|`--xml-dir`|none|TEI dir (indexed)|
|`--index-file` / `--load-index`|none|Index save/load|
|`--keyword`|none|Search term|
|`--scope`|`paragraph`|`paragraph`/`section`/`title`/`all`|
|`--case-sensitive`|off|Case-sensitive|
|`--context`|`100`|Context chars|
|`--export-csv`|none|Metadata CSV|
|`--show-structure`|none|Show structure|

## Outputs

Paper: `out/json/<name>.json` (full data), `out/markdown/<name>.md`, `out/html/<name>.html`, `out/tei/<name>.tei.xml` (raw GROBID, `--pdf-dir` only). Corpus: `out/corpus_index.json` (`--load-index`), `out/corpus_metadata.csv` (metadata, UTF-8 BOM), `out/search_hits.json` (hits), `out/intro_sentences.csv`.

`intro_sentences.csv`: `file, title, sentence_index, first_sentence`
`corpus_metadata.csv`: `file, title, keywords, n_sections, n_paragraphs, abstract_chars`

## Validation

- **GROBID** — probes `/api/isalive`; prints the deployment command when unreachable; ready once `Server started` appears in the log.
- **Parse rate** — summary line reports papers parsed successfully and failed (still in Chinese); one corrupt XML doesn't abort the batch — corrupt files become warnings.
- **Title** — `INTRODUCTION` or `Methods` means GROBID took a section name and the script falls back to the filename. **Falling back to the filename is correct behaviour**, and better than the journal name.
- **Keywords** — multi-word keywords (`tumor necrosis factor`) must stay intact; split into three words ⇒ whitespace-splitting fallback ran — check for `<term>` child elements.
- **Reference count** — must equal the entries at the paper's end; if 1 higher, the extra is the paper's own metadata (`<biblStruct>` inside `<sourceDesc>`).
- **Hierarchy** — `sections[].level` in `corpus_index.json` should reflect nesting; repeated headings ⇒ duplicate section trees.
- **Output dir** — artifacts flat under `out/json/`, `out/markdown/`, `out/html/`; nested subdirs ⇒ filename derived from an absolute path.
- **Lead sentence** — must be complete (ending with a period); mid-sentence truncation ⇒ an abbreviation (such as `e.g.`) was mistaken for the end.

## References

- `references/methodology.md` — TEI structure, GROBID deploy/params, section tree + lead-sentence algorithms, 5 defects, inventory
- `references/troubleshooting.md` — errors (unreachable, 404, misassigned title, split keywords, over-counted refs, escaping) + debugging
