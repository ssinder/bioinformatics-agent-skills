#!/usr/bin/env python3
"""生成 TEI XML 测试夹具。

夹具刻意覆盖了几种真实场景：

1. ``normal`` —— 标准 GROBID 输出：正常标题、多个作者、<term> 形式的
   关键词、多段摘要、嵌套章节、图表、参考文献。
2. ``title_is_section`` —— GROBID 把章节名误填进 title 元素（老仓库
   用来检测这种情况，但用的是"包含"匹配，会误杀正常标题）。
3. ``broken_words`` —— 标题含 PDF 断字（``INTRODUC TI ON``）。
4. ``space_keywords`` —— 关键词是纯文本形式，用逗号分隔。
5. ``no_abstract`` —— 没有摘要（老代码会直接崩溃或返回 None）。
6. ``malformed`` —— XML 语法错误（测试单文件失败不中断整批）。

用法::

    python tests/make_tei_fixtures.py
"""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
FIXTURES_DIR = REPO_ROOT / "tests" / "fixtures_tei"

TEI_NS = "http://www.tei-c.org/ns/1.0"


def _wrap(body: str, front: str = "") -> str:
    return f"""<?xml version="1.0" encoding="UTF-8"?>
<TEI xmlns="{TEI_NS}" xml:id="doc">
{front}
{body}
</TEI>
"""


NORMAL = _wrap(
    front="""
  <teiHeader>
    <fileDesc>
      <titleStmt>
        <title type="main">Deep learning frameworks for CRISPR guide RNA design</title>
        <author>
          <persName><forename type="first">Jane</forename><surname>Smith</surname></persName>
          <affiliation key="aff0">Department of Biology, Example University</affiliation>
        </author>
        <author>
          <persName><forename type="first">Bob</forename><surname>Jones</surname></persName>
          <affiliation key="aff1">Institute of Computation</affiliation>
        </author>
      </titleStmt>
      <publicationStmt>
        <publisher>Example Press</publisher>
        <availability><licence/></availability>
      </publicationStmt>
      <sourceDesc>
        <biblStruct>
          <analytic>
            <idno type="DOI">10.1234/example.2024.001</idno>
          </analytic>
          <monogr>
            <title>Journal of Example Biology</title>
            <imprint><date when="2024" /></imprint>
          </monogr>
        </biblStruct>
      </sourceDesc>
    </fileDesc>
    <profileDesc>
      <abstract>
        <div xmlns="http://www.tei-c.org/ns/1.0">
          <p>Guide RNA design pipelines have become a bottleneck in genome editing
          workflows, requiring both accuracy and interpretability.</p>
          <p>We present a benchmark of six deep learning architectures across
          eleven organisms, evaluating both predictive power and calibration.</p>
        </div>
      </abstract>
      <textClass>
        <keywords>
          <term>CRISPR</term>
          <term>guide RNA design</term>
          <term>deep learning</term>
          <term>benchmark</term>
        </keywords>
      </textClass>
    </profileDesc>
  </teiHeader>
""",
    body="""
  <text>
    <body>
      <div xmlns="http://www.tei-c.org/ns/1.0">
        <head n="1">Introduction</head>
        <p>CRISPR-Cas9 genome editing has transformed biological research.</p>
        <p>Despite progress, guide RNA design remains time-consuming for
        non-specialists. e.g. off-target effects complicate downstream assays.</p>
        <figure xml:id="fig_1">
          <head>Figure 1</head>
          <label>1</label>
          <figDesc>Comparison of architectures by AUROC across organisms.</figDesc>
          <graphic coords="1,1,100,100" url="figures/fig_1.png"/>
        </figure>
      </div>
      <div xmlns="http://www.tei-c.org/ns/1.0">
        <head n="2">Methods</head>
        <p>We evaluated six architectures on eleven organisms with three splits.</p>
        <div xmlns="http://www.tei-c.org/ns/1.0">
          <head n="2.1">Data collection</head>
          <p>Genomic sequences were downloaded from RefSeq and aligned with MAFFT.</p>
        </div>
        <div xmlns="http://www.tei-c.org/ns/1.0">
          <head n="2.2">Model training</head>
          <p>All models were trained with the same data split and seed.</p>
        </div>
      </div>
      <div xmlns="http://www.tei-c.org/ns/1.0">
        <head n="3">Results</head>
        <p>The best model achieved an AUROC of 0.89 across all organisms.</p>
        <figure xml:id="fig_2">
          <head>Figure 2</head>
          <label>2</label>
          <figDesc>Calibration curves for the top three models.</figDesc>
          <graphic coords="1,1,100,100" url="figures/fig_2.png"/>
        </figure>
      </div>
      <div xmlns="http://www.tei-c.org/ns/1.0">
        <head n="4">Conclusion</head>
        <p>Interpretable architectures matched black-box models on this task.</p>
      </div>
    </body>
    <back>
      <div type="references">
        <listBibl>
          <biblStruct xml:id="b0">
            <analytic>
              <title>Foundation models for genomics</title>
              <idno type="DOI">10.1038/s41586-023-00001-1</idno>
            </analytic>
          </biblStruct>
          <biblStruct xml:id="b1">
            <analytic>
              <title>Efficient transformer variants</title>
              <idno type="DOI">10.1038/s41586-023-00002-2</idno>
            </analytic>
          </biblStruct>
        </listBibl>
      </div>
    </back>
  </text>
""",
)

TITLE_IS_SECTION = _wrap(
    front="""
  <teiHeader>
    <fileDesc>
      <titleStmt>
        <title type="main">INTRODUCTION</title>
        <author>
          <persName><forename type="first">Carol</forename><surname>White</surname></persName>
        </author>
      </titleStmt>
      <publicationStmt><publisher>Example Press</publisher></publicationStmt>
      <sourceDesc>
        <biblStruct>
          <monogr><title>Another Journal</title></monogr>
        </biblStruct>
      </sourceDesc>
    </fileDesc>
    <profileDesc>
      <textClass>
        <keywords>
          <term>single cell</term>
          <term>trajectory inference</term>
        </keywords>
      </textClass>
    </profileDesc>
  </teiHeader>
""",
    body="""
  <text>
    <body>
      <div xmlns="http://www.tei-c.org/ns/1.0">
        <head n="1">Introduction</head>
        <p>Single-cell RNA sequencing reveals cellular heterogeneity at scale.</p>
        <p>Trajectory inference methods differ substantially in their assumptions.</p>
      </div>
      <div xmlns="http://www.tei-c.org/ns/1.0">
        <head n="2">Results</head>
        <p>Our benchmark spans fourteen datasets.</p>
      </div>
    </body>
  </text>
""",
)

BROKEN_WORDS = _wrap(
    front="""
  <teiHeader>
    <fileDesc>
      <titleStmt>
        <title type="main">INTRODUC TI ON</title>
        <author>
          <persName><forename type="first">Dave</forename><surname>Brown</surname></persName>
        </author>
      </titleStmt>
      <publicationStmt><publisher>Example Press</publisher></publicationStmt>
      <sourceDesc>
        <biblStruct><monogr><title>Third Journal</title></monogr></biblStruct>
      </sourceDesc>
    </fileDesc>
    <profileDesc>
      <textClass>
        <keywords>
          <term>protein folding</term>
        </keywords>
      </textClass>
    </profileDesc>
  </teiHeader>
""",
    body="""
  <text>
    <body>
      <div xmlns="http://www.tei-c.org/ns/1.0">
        <head n="1">INTRODUC TI ON</head>
        <p>Protein folding prediction has advanced rapidly in recent years.</p>
      </div>
    </body>
  </text>
""",
)

SPACE_KEYWORDS = _wrap(
    front="""
  <teiHeader>
    <fileDesc>
      <titleStmt>
        <title type="main">Multi-omic integration for cancer subtyping</title>
        <author>
          <persName><forename type="first">Eve</forename><surname>Green</surname></persName>
        </author>
      </titleStmt>
      <publicationStmt><publisher>Example Press</publisher></publicationStmt>
      <sourceDesc>
        <biblStruct><monogr><title>Fourth Journal</title></monogr></biblStruct>
      </sourceDesc>
    </fileDesc>
    <profileDesc>
      <textClass>
        <keywords>tumor necrosis factor, glycolysis pathway, immune infiltration</keywords>
      </textClass>
    </profileDesc>
  </teiHeader>
""",
    body="""
  <text>
    <body>
      <div xmlns="http://www.tei-c.org/ns/1.0">
        <head n="1">Introduction</head>
        <p>Multi-omic integration reveals tumour cell states.</p>
      </div>
    </body>
  </text>
""",
)

NO_ABSTRACT = _wrap(
    front="""
  <teiHeader>
    <fileDesc>
      <titleStmt>
        <title type="main">A minimal record without abstract or keywords</title>
        <author>
          <persName><forename type="first">Frank</forename><surname>Miller</surname></persName>
        </author>
      </titleStmt>
      <publicationStmt><publisher>Example Press</publisher></publicationStmt>
      <sourceDesc>
        <biblStruct><monogr><title>Fifth Journal</title></monogr></biblStruct>
      </sourceDesc>
    </fileDesc>
  </teiHeader>
""",
    body="""
  <text>
    <body>
      <div xmlns="http://www.tei-c.org/ns/1.0">
        <head n="1">Introduction</head>
        <p>This document has neither abstract nor keyword elements.</p>
      </div>
    </body>
  </text>
""",
)

MALFORMED = """<?xml version="1.0" encoding="UTF-8"?>
<TEI xmlns="http://www.tei-c.org/ns/1.0">
  <teiHeader>
    <fileDesc>
      <titleStmt>
        <title type="main">Broken document
  </teiHeader>
</TEI>
"""


FIXTURE_FILES = {
    "normal.paper.pdf.tei.xml": NORMAL,
    "title_is_section.paper.pdf.tei.xml": TITLE_IS_SECTION,
    "broken_words.paper.pdf.tei.xml": BROKEN_WORDS,
    "space_keywords.paper.pdf.tei.xml": SPACE_KEYWORDS,
    "no_abstract.paper.pdf.tei.xml": NO_ABSTRACT,
    "malformed.paper.pdf.tei.xml": MALFORMED,
}


def main() -> int:
    """写出全部夹具。"""
    FIXTURES_DIR.mkdir(parents=True, exist_ok=True)
    for name, content in FIXTURE_FILES.items():
        path = FIXTURES_DIR / name
        path.write_text(content, encoding="utf-8")
        print(f"已写出 {path.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
