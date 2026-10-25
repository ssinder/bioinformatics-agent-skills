# 快速上手 / Getting Started

本指南介绍如何安装、配置并首次使用 BioinfoAgent Skills。

## 安装

### 作为 Claude Code 技能

```bash
git clone https://github.com/ssinder/bioinformatics-agent-skills.git
cp -r bioinformatics-agent-skills/skills/* ~/.claude/skills/
```

重启 Claude Code 后技能即可被自动发现。

> 也支持其他遵循 [Agent Skills 标准](https://agentskills.io/) 的宿主:
> Cursor、Codex、Google Antigravity 等。把 `skills/` 目录放到对应宿主读取技能的位置即可。

### 作为 Python 包

```bash
git clone https://github.com/ssinder/bioinformatics-agent-skills.git
cd bioinformatics-agent-skills
pip install -e ".[all]"     # 安装全部科学计算依赖
# 或按需安装
pip install -e ".[rna-seq]"
pip install -e ".[single-cell]"
pip install -e ".[structure]"
pip install -e ".[visualization]"
```

验证安装:

```bash
biopipe --list
```

## 配置

### 项目目录约定

所有技能遵循统一的目录约定:

```
my_project/
├── data/         # 原始输入
├── work/         # 中间产物
├── results/      # 最终结果
├── figures/      # 图表
└── logs/         # 运行日志
```

在项目根目录运行脚本,目录会被自动创建。

### 配置文件

在项目根目录放 `biopipe.yaml`,可自定义目录:

```yaml
root: /path/to/project
work_dir: intermediate
output_dir: final_results
figures_dir: plots
```

也支持 JSON 格式,或通过命令行 `--config` 指定其他路径。

### 环境变量

| 变量 | 用途 |
|------|------|
| `BIOPIPE_PROJECT_ROOT` | 覆盖项目根目录 |
| `NCBI_API_KEY` | NCBI E-utilities API Key(可选,提升速率限制) |
| `OPENALEX_MAILTO` | OpenAlex polite pool 邮箱(可选) |

## 首次运行

```bash
export PYTHONPATH="$PWD/packages"

# 1. 查看技能清单
biopipe --list

# 2. 查看某个技能的完整用法
biopipe bulk-rnaseq --help

# 3. 小步试跑,先看参数说明
python skills/peptide-design/scripts/analyze_peptides.py --help
```

## 日志与可复现

所有脚本支持统一日志参数:

```bash
# 详细日志 + 写入文件
python skills/<name>/scripts/<entry>.py --input ... --log-file logs/run.log -v

# 静默模式(适合自动化流水线)
python skills/<name>/scripts/<entry>.py --input ... -q
```

日志同时写到 stderr 与指定文件,格式为:

```
2026-10-25 20:46:22 | INFO    | biopipe.bulk-rnaseq | 开始质控
```

## 常见安装问题

### Q: 提示缺少 scanpy / OpenMM / MDAnalysis

**A**: 这是设计如此。核心算法都有内置的纯 Python 实现,会打印一条提示后自动降级:

```
WARNING | biopipe.single-cell-rnaseq | 未检测到 scanpy，使用内置降级后端
```

想要完整功能,安装对应可选依赖即可。

### Q: matplotlib 中文显示为方块

**A**: matplotlib 默认字体不含中文字形。设置环境变量或修改 `figure_style.py`:

```python
import matplotlib
matplotlib.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False
```

### Q: 网络请求失败

**A**: 共享层已内置指数退避重试与限速。若仍失败:

1. 检查是否需要代理:`export https_proxy=http://127.0.0.1:7897`
2. NCBI 等服务建议设置 `NCBI_API_KEY`
3. 用 `-v` 查看具体的 HTTP 状态码

### Q: 端口或服务名找不到（HDOCK 等外部服务）

**A**: 部分对接服务需要本地运行。这类脚本默认不会联网,必须显式指定服务端点。

## 下一步

- 了解各技能能做什么 → [技能目录]({{ '/skills.html' | relative_url }})
- 了解技能规范并贡献新技能 → [SKILL_SPEC.md]({{ '/SKILL_SPEC.html' | relative_url }})
- 了解共享工具层 API → [AGENTS.md](https://github.com/ssinder/bioinformatics-agent-skills/blob/main/AGENTS.md)
