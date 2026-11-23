# Getting Started

This guide walks you through installing, configuring, and using BioinfoAgent Skills for the first time.

## Installation

### With the skills CLI (recommended)

```bash
npx skills add ssinder/bioinformatics-agent-skills
```

This works with 80+ agent hosts (Claude Code, Codex, Cursor, OpenCode, …) and
handles the per-host skill directories automatically.

To install only what a given workflow needs, name the skills explicitly:

```bash
npx skills add ssinder/bioinformatics-agent-skills \
  --skill ena-download --skill bulk-rnaseq --skill scientific-figures
```

> **Install selectively.** The agent reads every installed skill's description at
> startup to decide whether to load one. When too many are installed, the host
> drops some descriptions to fit its context budget — and a skill without a
> description can never be matched, so it never fires. Two to three per workflow
> is the effective range. The [skill guides]({{ '/skill-guides/' | relative_url }})
> list recommended combinations per research task.

Restart your host and the skills are discovered automatically.

> Other hosts that follow the [Agent Skills standard](https://agentskills.io/)
> are supported as well — Cursor, Codex, Google Antigravity, and others. The CLI
> handles their directory conventions; if you prefer to install manually, see
> [Support](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/SUPPORT.md)
> for the per-host paths.

### As a Python package

```bash
pip install biopipe-agent-skills          # from PyPI
# or from a checkout
git clone https://github.com/ssinder/bioinformatics-agent-skills.git
cd bioinformatics-agent-skills
pip install -e ".[all]"     # install all scientific computing dependencies
# or install only what you need
pip install -e ".[rna-seq]"
pip install -e ".[single-cell]"
pip install -e ".[structure]"
pip install -e ".[visualization]"
```

The wheel ships the skill packages themselves, so the installed command runs the
real workflows — no checkout required.

Verify the installation:

```bash
bioagent --list
```

> The console script is `bioagent`. The old `biopipe` name still works as an
> alias, because `biopipe` is already taken on PyPI by an unrelated package.

## Configuration

### Project layout

All skills follow the same directory convention:

```
my_project/
├── data/         # raw inputs
├── work/         # intermediate artifacts
├── results/      # final results
├── figures/      # plots and charts
└── logs/         # run logs
```

Run your scripts from the project root — the directories are created automatically.

### Configuration file

Place a `bioagent.yaml` in the project root to customise the directory layout:

```yaml
root: /path/to/project
work_dir: intermediate
output_dir: final_results
figures_dir: plots
```

JSON is supported as well, and you can point at a different file with the `--config` command-line option.

### Environment variables

| Variable | Purpose |
|------|------|
| `BIOPIPE_PROJECT_ROOT` | Override the project root directory |
| `NCBI_API_KEY` | NCBI E-utilities API key (optional, raises the rate limit) |
| `OPENALEX_MAILTO` | Contact email for the OpenAlex polite pool (optional) |

## First run

```bash
export PYTHONPATH="$PWD/packages"

# 1. List the available skills
bioagent --list

# 2. Show the full usage of a specific skill
bioagent bulk-rnaseq --help

# 3. Try a small step first to review the parameters
python skills/peptide-design/scripts/analyze_peptides.py --help
```

## Logging and reproducibility

Every script accepts the same logging options:

```bash
# verbose logging + write to a file
python skills/<name>/scripts/<entry>.py --input ... --log-file logs/run.log -v

# quiet mode (suitable for automated pipelines)
python skills/<name>/scripts/<entry>.py --input ... -q
```

Logs are written to stderr and to the given file at the same time, using this format:

```
2026-10-25 20:46:22 | INFO    | bioagent.bulk-rnaseq | starting QC
```

## Common installation issues

### Q: A warning about missing scanpy / OpenMM / MDAnalysis

**A**: This is by design. Every core algorithm ships with a built-in pure-Python implementation, so the skill prints a notice and falls back automatically:

```
WARNING | bioagent.single-cell-rnaseq | scanpy is not available, switching to the numpy fallback backend
```

Install the matching optional dependency whenever you need the full functionality.

### Q: matplotlib renders CJK characters as boxes

**A**: The default matplotlib fonts contain no CJK glyphs, so Chinese, Japanese, and Korean text is drawn as empty boxes. Install a CJK font on your system (for example Noto Sans CJK), then list it first in the sans-serif family. The quickest way is to edit the `font.sans-serif` entry of `figure_style.py`:

```python
import matplotlib
matplotlib.rcParams["font.sans-serif"] = ["Noto Sans CJK JP", "Microsoft YaHei", "SimHei", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False
```

Note that matplotlib only picks up fonts it can find on the system, so install the font first and restart your session before rebuilding the font cache.

### Q: Network requests fail

**A**: The shared layer already implements exponential-backoff retries and rate limiting. If failures persist:

1. Check whether a proxy is required: `export https_proxy=http://127.0.0.1:7897`
2. For services such as NCBI, set `NCBI_API_KEY`
3. Run with `-v` to see the exact HTTP status codes

### Q: Port or service name not found (external services such as HDOCK)

**A**: Some docking services have to run locally. These scripts never reach the network by default — you must specify the server endpoint explicitly.

## Next steps

- See what each skill can do → [Skills catalog]({{ '/skills/' | relative_url }})
- Read the skill specification and contribute a new skill → [SKILL_SPEC.md]({{ '/SKILL_SPEC.html' | relative_url }})
- Explore the shared utility layer API → [AGENTS.md](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/AGENTS.md)
