# Skill Refactoring Specification Contract

All skill sub-agents MUST strictly comply with this document. The refactored skill directory structure is fixed as follows:

```
skills/<skill-name>/
├── SKILL.md              # Required, Agent Skills standard format
├── scripts/              # Executable scripts, path-parameterized + logging-based
│   └── <entry>.py
├── references/           # Reference documents (optional but recommended)
│   └── *.md
├── requirements.txt      # Dependencies (optional)
└── README.md             # Short human-readable description (optional)
```

## SKILL.md Format (Strictly Enforced)

The YAML frontmatter MUST contain:

```yaml
---
name: <skill-name>          # Lowercase with hyphens, matching the directory name
description: <third-person description; must contain the sentence "This skill should be used when..." explaining when to trigger>
license: MIT
compatibility: <runtime requirements, e.g. Python 3.9+ plus the key dependencies>
metadata:
  version: "1.0"
  last-reviewed: "2026-10-25"
  skill-author: sinder
---
```

Requirements for `description`:

- Write it in the third person, as "This skill should be used when..."
- It MUST state both **what the skill does** and **when to trigger** (including trigger keywords; English keywords must appear, because users may ask in English)
- Keep the length within 100-400 characters
- Never use second-person phrasing such as "Use this skill when"

The body MUST contain the following sections (as `##` level-2 headings):

1. `## Overview` -- what problem this skill solves, 2-4 sentences
2. `## When to Use This Skill` -- a list of concrete usage scenarios (use `-`)
3. `## Quick Start` -- command examples that can be copied and run directly
4. `## Workflow` -- step-by-step instructions (numbered list or mermaid flowchart)
5. `## Parameters` -- a table of the key parameters
6. `## Outputs` -- a table of the produced files
7. `## Validation` -- how to verify that the results are correct (quality-control checkpoints)
8. `## References` -- pointers to the documents under references/

## Script Rules (Strictly Enforced)

1. **Import the shared layer**: add this at the top of the script
   ```python
   sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "packages"))
   from biopipe_core import get_logger, get_config, SkillError
   ```
   (`parents[3]` is used because the script lives at `skills/<name>/scripts/x.py`, three levels up reaches the repository root)

2. **Never hardcode absolute paths**: every input and output path comes from a command-line argument; defaults are derived from `get_config()`

3. **Never use bare `print()`**: use `logger = get_logger("biopipe.<skill-name>")` and emit output with `logger.info(...)`

4. **Unified CLI**: use `biopipe_core.cli.bootstrap_cli(description, skill_name, main)` as the entry point; `main(args)` returns the exit code. It MUST expose business parameters such as `--input` and `--output`.

5. **File header**: every script MUST carry a shebang + encoding declaration + module docstring (including a usage example)

6. **MUST have an `if __name__ == "__main__":` entry point**

7. **The script MUST exit cleanly with `--help` when no input file is present** (it must not depend on real data)

8. **Modules that do not rely on heavy third-party libraries MUST degrade gracefully at import time**: for example, when scanpy/matplotlib are missing, report a clear error instead of an ImportError stack trace

## Content Migration Sources

The original repository `Bioinformatics_pipeline/<corresponding directory>/scripts/*.py` is the migration source.
**Preserve the original scripts' scientific computation logic**, but:
- Remove the hardcoded paths and turn them into parameters
- Remove the bare `print()` calls and switch them to the logger
- Add argparse
- Keep the original core algorithms and comments

Original directory location: `C:/Users/DELL/WorkBuddy/2026-10-25-20-38-30/Bioinformatics_pipeline/`

## Prohibited Practices

- Do not modify the original repository `Bioinformatics_pipeline/` (read-only reference)
- Do not create temporary files other than `README.md`
- Do not write real names, local paths, QQ numbers, email addresses, or other private information in SKILL.md
- Do not use emoji (none inside SKILL.md; the README top level excepted)
- Do not make network requests the default behavior of a script (networking requires an explicit parameter)

## Pre-delivery Self-check

Before finishing, you MUST run the following self-check and report the results:
```bash
cd <repo>
python -c "import yaml,sys; d=yaml.safe_load(open('skills/<name>/SKILL.md').read().split('---')[1]); print(d['name'], '|', d['description'][:80])"
python skills/<name>/scripts/<entry>.py --help
```
