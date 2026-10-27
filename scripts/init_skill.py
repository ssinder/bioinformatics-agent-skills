#!/usr/bin/env python3
"""Generate a skill skeleton that conforms to the Agent Skills specification.

Usage:
    python scripts/init_skill.py my-new-skill
    python scripts/init_skill.py my-new-skill --category Transcriptomics --force

Generated structure:
    skills/my-new-skill/
    ├── SKILL.md
    ├── scripts/run.py
    ├── references/methodology.md
    ├── references/troubleshooting.md
    ├── requirements.txt
    └── README.md
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Optional

REPO_ROOT = Path(__file__).resolve().parents[1]
SKILLS_DIR = REPO_ROOT / "skills"

SKILL_MD_TEMPLATE = """---
name: {name}
description: This skill should be used when a user needs to {summary}, including analysis of {data_desc} with reproducible and validated scientific methods. 用于{description_cn}的场景。
license: MIT
compatibility: Requires Python 3.9+ with numpy and pandas.
metadata:
  version: "1.0"
  last-reviewed: "{today}"
  skill-author: {author}
---

# {title}

## Overview

{summary}.

This skill encapsulates the analysis workflow for {data_desc}. It provides reproducible
computational steps, the scientific rationale behind parameter choices, and methods to
verify that the results are correct.

## When to Use This Skill

Use this skill when the user wants to:
- Perform {summary_cn} on {data_desc}
- Process many samples in batch while keeping the parameters consistent
- Understand the scientific rationale behind the parameter choices, not just obtain the results
- Record the analysis as a reproducible workflow

For **other** kinds of data, use the dedicated skill that matches them. The authoritative
list of skill names is the output of `biopipe --list`.

## Quick Start

```bash
export PYTHONPATH="$PWD/packages"

# Show the full usage
python skills/{name}/scripts/{entry} --help

# Basic usage
python skills/{name}/scripts/{entry} \\
    --input {example_input} \\
    --output {example_output}
```

## Workflow

```mermaid
flowchart TD
    A["Prepare input data"] --> B["Quality control and checks"]
    B --> C["Core analysis step"]
    C --> D["Summarize results"]
    D --> E["Validate and plot"]
    E --> F["Write the final report"]
```

1. **Prepare input data** -- Confirm the data format and quality, record the source and version
2. **Quality control and checks** -- Complete the basic checks before the core step, record the key metrics
3. **Core analysis step** -- Follow the standard workflow and fix the random seed
4. **Summarize results** -- Emit structured intermediate results and statistics
5. **Validate and plot** -- Confirm the results are reliable using the methods of the Validation section
6. **Write the final report** -- Collect every artifact and record the run parameters

## Parameters

| Parameter | Type | Default | Description |
|------|------|--------|------|
| `--input` | Path | required | Path to the input file |
| `--output` | Path | `results/` | Output directory |
| `--config` | Path | none | Project configuration file (JSON or key=value) |
| `--log-file` | Path | none | Log output file |
| `-v/--verbose` | flag | off | Emit DEBUG-level logs |
| `-q/--quiet` | flag | off | Emit warnings and errors only |

## Outputs

| File | Description |
|------|------|
| `results/{name}_report.txt` | Run summary and key statistics |
| `results/{name}_results.csv` | Structured results table |
| `figures/` | Generated figures |
| `logs/` | Full run logs |

## Validation

Result reliability checklist:

- [ ] The number of input samples matches the expectation, with nothing missing
- [ ] Key metrics fall inside reasonable ranges (record the values measured in this run)
- [ ] The number of result rows matches the number of input records
- [ ] Every source of randomness is fixed (record the random seed)
- [ ] Outliers are flagged instead of being silently dropped

**Note**: this skill produces supporting analysis results, so the final judgement must
combine them with domain knowledge and the experimental design.

## References

- [Methodology](references/methodology.md) -- The rationale and scientific background behind the parameter choices
- [Troubleshooting](references/troubleshooting.md) -- Symptoms, root causes, and solutions
"""

SCRIPT_TEMPLATE = '''#!/usr/bin/env python3
"""{title} -- entry script.

Usage examples:
    python skills/{name}/scripts/{entry} --help
    python skills/{name}/scripts/{entry} --input example.txt --output results/
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "packages"))

from biopipe_core import bootstrap_cli, get_config, get_logger


def main(args) -> int:
    """Main entry point of the skill.

    Args:
        args: The parsed command-line arguments.

    Returns:
        The exit code; 0 means success.
    """
    logger = get_logger("biopipe.{name}")
    config = get_config(args.config)

    logger.info("Skill started: {name}")
    logger.info("Project root: %s", config.root)
    logger.info("Output directory: %s", config.output_dir)

    # TODO: implement the core logic
    logger.info("The core logic is not implemented yet")
    return 0


if __name__ == "__main__":
    sys.exit(bootstrap_cli("{title}", "{name}", main))
'''

METHODOLOGY_TEMPLATE = """# {title} -- Methodology

This document explains the **scientific rationale** and the **parameter choices** behind
every step of this skill, so that researchers can judge whether the results are reasonable.

## Overview

{summary}.

The analysis workflow follows the standard practice of the field. The key decision
points are described below.

## Key Parameters and Their Rationale

### Parameter One

- **Valid range**: the reasonable interval and what it means
- **Why it was chosen**: the scientific reason
- **Sensitivity and how to check it**: how to verify whether this choice affects the conclusion

### Parameter Two

- **Valid range**: the reasonable interval and what it means
- **Why it was chosen**: the scientific reason
- **Sensitivity and how to check it**: how to verify whether this choice affects the conclusion

## Limitations of the Method

Be honest about the boundaries of this method:

- The situations in which the results are unreliable
- The confounders that were not taken into account
- The extra validation needed to support the conclusion

## References

- Standard methodology literature of the field
"""

TROUBLESHOOTING_TEMPLATE = """# {title} -- Troubleshooting

## Input and Data Problems

### Symptom: the input file cannot be read

**Possible causes**

1. The file path does not exist or is misspelled
2. The file encoding is not UTF-8
3. The file format differs from the expectation (changed header column names)

**Solutions**

1. Run with `--verbose` to see the path that was actually resolved
2. Convert the file encoding: `iconv -f gbk -t utf-8 input.csv > input_utf8.csv`
3. Check whether the header column names still match the documentation

## Computation and Results Problems

### Symptom: the results are empty or the number of rows is unexpected

**Debugging steps**

1. Confirm the number of input records
2. Check whether any record was filtered out by the quality control step
3. Inspect the filtering ratio in the log; a high filtering rate usually means the parameters are unsuitable

### Symptom: the numerical results are clearly abnormal

**Possible causes**

1. The parameters fall outside the reasonable range
2. The input data format is wrong
3. The behaviour changed because of a version difference

**Solutions**

1. Compare against the parameter table in `references/methodology.md`
2. Validate the pipeline with a minimal sample dataset
3. Record the environment versions and compare against the reference results

## Performance Problems

### Symptom: the run takes too long or uses too much memory

**Optimization suggestions**

1. Process the data in batches
2. Lower the default sampling density
3. Check whether unnecessary intermediate results are written to disk

## How to Ask Questions

For anything not covered here, please open an Issue and provide:

- The complete command line and the `--verbose` log
- A sample of the input data (de-identified)
- The difference between the expected and the actual results
- The Python and key dependency versions
"""

REQUIREMENTS_TEMPLATE = """# {title} dependencies
requests>=2.28
numpy>=1.23
pandas>=1.5
"""

README_TEMPLATE = """# {title}

{summary}.

## Quick Start

```bash
export PYTHONPATH="$PWD/packages"
python skills/{name}/scripts/{entry} --help
```

See [SKILL.md](SKILL.md) for the full description.

- Methodology: [references/methodology.md](references/methodology.md)
- Troubleshooting guide: [references/troubleshooting.md](references/troubleshooting.md)
"""


def main(argv: Optional[list[str]] = None) -> int:
    """Generate a skill skeleton."""
    parser = argparse.ArgumentParser(
        description="Generate a skill skeleton that conforms to the Agent Skills specification"
    )
    parser.add_argument("name", help="Skill name (kebab-case)")
    parser.add_argument(
        "--summary", default="perform domain-specific analysis", help="Short summary of what the skill does"
    )
    parser.add_argument(
        "--summary-cn",
        default="domain-specific analysis",
        help="Short Chinese summary of what the skill does (the Chinese trigger wording used in the skill body)",
    )
    parser.add_argument(
        "--description-cn",
        default="domain-specific analysis",
        help="Chinese wording for the description field (the Chinese trigger phrase appended after the English text)",
    )
    parser.add_argument("--data-desc", default="domain data", help="Description of the data type being processed")
    parser.add_argument("--title", default=None, help="Skill title")
    parser.add_argument("--category", default="Other", help="Skill category")
    parser.add_argument("--author", default="sinder", help="Skill author")
    parser.add_argument("--entry", default="run.py", help="Entry script file name")
    parser.add_argument("--force", action="store_true", help="Overwrite an existing directory")
    args = parser.parse_args(argv)

    if not re.match(r"^[a-z0-9]+(?:-[a-z0-9]+)*$", args.name):
        print(
            f"Invalid skill name: {args.name} (only lowercase letters, digits, and hyphens are allowed)",
            file=sys.stderr,
        )
        return 2

    skill_dir = SKILLS_DIR / args.name
    if skill_dir.exists():
        if not args.force:
            print(f"Directory already exists: {skill_dir} (use --force to overwrite it)", file=sys.stderr)
            return 1
        import shutil

        shutil.rmtree(skill_dir)

    title = args.title or args.name.replace("-", " ").title()
    from datetime import date

    today = date.today().isoformat()

    (skill_dir / "scripts").mkdir(parents=True)
    (skill_dir / "references").mkdir(parents=True)

    ctx = {
        "name": args.name,
        "title": title,
        "summary": args.summary,
        "summary_cn": args.summary_cn,
        "description_cn": args.description_cn,
        "data_desc": args.data_desc,
        "author": args.author,
        "entry": args.entry,
        "today": today,
        "example_input": "input/data.csv",
        "example_output": "results/",
    }

    (skill_dir / "SKILL.md").write_text(SKILL_MD_TEMPLATE.format(**ctx), encoding="utf-8")
    (skill_dir / "scripts" / args.entry).write_text(SCRIPT_TEMPLATE.format(**ctx), encoding="utf-8")
    (skill_dir / "references" / "methodology.md").write_text(
        METHODOLOGY_TEMPLATE.format(**ctx), encoding="utf-8"
    )
    (skill_dir / "references" / "troubleshooting.md").write_text(
        TROUBLESHOOTING_TEMPLATE.format(**ctx), encoding="utf-8"
    )
    (skill_dir / "requirements.txt").write_text(REQUIREMENTS_TEMPLATE.format(**ctx), encoding="utf-8")
    (skill_dir / "README.md").write_text(README_TEMPLATE.format(**ctx), encoding="utf-8")

    print(f"Skill skeleton generated: {skill_dir}")
    print(f"Category: {args.category}")
    print("\nNext steps:")
    print(f"  1. Edit the description of {skill_dir / 'SKILL.md'} (it decides when the skill is triggered)")
    print(f"  2. Implement the core logic of {skill_dir / 'scripts' / args.entry}")
    print("  3. Write the two reference documents under references/")
    print(f"  4. Validate: python scripts/validate_skills.py --skill {args.name} --strict")
    return 0


if __name__ == "__main__":
    sys.exit(main())
