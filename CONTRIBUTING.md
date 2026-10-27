# Contributing Guide

Thank you for considering contributing a skill to BioinfoAgent Skills. This project follows the
[Agent Skills open standard](https://agentskills.io/), and we welcome skill submissions from
researchers in any field.

---

## What to Contribute

We welcome three kinds of contributions:

| Type | Description | Who it suits |
|------|------|--------|
| **New Skill** | Submit a brand-new domain skill | Researchers who have a concrete analysis pipeline |
| **Improve Skill** | Fix bugs, add parameters, improve the documentation | All users |
| **Infrastructure** | Improve the shared tooling layer, CI, and documentation | Developers who want to help with the engineering |

## Quick Start: Contributing a New Skill

### 1. Fork the repository and create a branch

```bash
git clone https://github.com/ssinder/bioinformatics-agent-skills.git
cd bioinformatics-agent-skills
git checkout -b feature/my-new-skill
```

### 2. Generate the skill skeleton

```bash
python scripts/init_skill.py my-new-skill
```

This creates the `skills/my-new-skill/` directory containing the `SKILL.md` template, `scripts/`, `references/`, and `assets/`.

### 3. Write the SKILL.md

**This is the most important step.** The SKILL.md decides when the AI triggers your skill.

```yaml
---
name: my-new-skill
description: This skill should be used when a user needs to perform X analysis on Y data, including A and B steps.
license: MIT
compatibility: Requires Python 3.9+ with numpy and pandas.
metadata:
  version: "1.0"
  last-reviewed: "2026-10-25"
  skill-author: your-name
---
```

The three iron rules for `description`:

1. **Third person** -- write "This skill should be used when...", never "You should use this..."
2. **State clearly when to trigger** -- include the keywords users would actually type (include both English and Chinese, since users may ask in either language)
3. **Be specific** -- name the data formats, tools, and analysis steps; never write empty filler

### 4. Write the scripts

All scripts must follow the project conventions:

```python
#!/usr/bin/env python3
"""Skill description and usage example."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "packages"))

from biopipe_core import get_logger, get_config, bootstrap_cli


def main(args) -> int:
    logger = get_logger("biopipe.my-new-skill")
    config = get_config(args.config)
    logger.info("Processing input: %s", args.input)
    return 0


if __name__ == "__main__":
    sys.exit(bootstrap_cli("description", "my-new-skill", main))
```

**Hard rules**:

- ❌ Never hardcode absolute paths (`/Users/xxx/...`, `C:/Users/xxx/...`)
- ❌ Never use bare `print()`; always use `logger`
- ✅ Always provide a command-line interface with `argparse`
- ✅ Always include an `if __name__ == "__main__":` entry point
- ✅ Heavy dependencies (scanpy, OpenMM, ...) must degrade gracefully with a clear error message
- ✅ Import the shared helpers from `biopipe_core`

### 5. Validate locally

```bash
python scripts/validate_skills.py --skill my-new-skill
python scripts/validate_skills.py --skill my-new-skill --strict   # warnings count as failures
```

Zero errors and zero warnings are required before a skill can be merged.

### 6. Open a Pull Request

```bash
git add .
git commit -m "feat: add my-new-skill for X analysis"
git push origin feature/my-new-skill
```

Then open the Pull Request. CI automatically runs the specification validation, syntax checks, dependency audit, and smoke tests.

## Quality Bar

| Requirement | Description |
|------|------|
| **Scientific correctness** | The algorithms must be genuinely usable, not pseudocode placeholders |
| **Reproducibility** | Pin the random seeds and record the parameters and versions |
| **A fallback path** | The core algorithm should embed a pure-Python implementation, so it does not hard-depend on a heavy library |
| **Complete documentation** | You MUST provide `references/methodology.md` (the methodological rationale) and `references/troubleshooting.md` (debugging) |
| **A validation stage** | The Validation section of SKILL.md explains how to check that the results are correct |

Scientific code must be genuinely runnable, not a stub.

## Need Help?

- Questions about the skill specification -- open an Issue and label it `question`
- Not sure how to write the `description` -- open an Issue describing your analysis scenario and we will help you write it
- Want more examples -- see `docs/skill-guides/`
