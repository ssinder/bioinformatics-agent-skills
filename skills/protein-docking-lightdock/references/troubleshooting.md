# LightDock 故障排查

按"症状 → 原因 → 处理"的顺序整理。全部问题都来自实际跑 LightDock 时的
典型故障模式。

---

## 找不到 lightdock3_setup.py / lightdock3.py

**症状**

```
[ERROR] 找不到 LightDock 可执行脚本: lightdock3_setup.py, lightdock3.py。
请先 `pip install lightdock==0.9.4`（Windows 官方不支持，请用 Linux/macOS），
或用 --setup-exe/--run-exe 显式指定路径。
```

**原因**

LightDock 通过 `setup.py` 的 `scripts=[...]` 声明安装脚本，它们是**裸文件**，
落在 `<env>/bin/` 下，但**不在 `PATH` 里**（只有 `pip install --user` 或系统级
安装才可能进 PATH）。

**处理**

```bash
# 看脚本到底在哪
python -c "import sys; print(sys.prefix)"
ls "$(python -c 'import sys; print(sys.prefix)')/bin" | grep lightdock

# 显式指定
python scripts/run_docking.py ... --setup-exe /path/to/bin/lightdock3_setup.py \
                                 --run-exe /path/to/bin/lightdock3.py
```

本技能已经会自动检查 `<sys.prefix>/bin`，所以在**激活了 lightdock 环境的
Python** 里直接跑通常就能找到。关键是**要用那个环境的 Python**：

```bash
conda activate lightdock
python scripts/run_docking.py ...
```

---

## Windows 上装不上

**症状**：`pip install lightdock` 编译 Cython 扩展失败。

**原因**：LightDock 官方**不支持 Windows**。它的 Cython 扩展用了
`lightdock/gso/searchspace/cython/` 下的 `.c` 文件，理论上能编译，但 ProDy、
freesasa 等依赖在 Windows 上的组合问题很多。

**处理**：用 WSL2、Linux 或 macOS。WSL2 里记得 `apt install build-essential`，
并且**不要**在 Windows 侧与 WSL 侧混用路径（PDB 文件用 WSL 的路径）。

---

## 准备阶段报 `LightdockError`

**症状**

```
LightDock 准备阶段失败（退出码 1）:
[ERROR] Exception while reading structure ...
```

**常见原因与处理**

| 报错关键词 | 原因 | 处理 |
|---|---|---|
| `not recognized` / `parse` | PDB 文件格式有问题 | 用 `pdb-tools` 或 Biopython 校验；检查是否有 CRLF 换行 |
| `no atoms` | PDB 里没有 `ATOM`/`HETATM` 记录 | 确认文件不是 `.cif` 改名 |
| `already exists, please delete previous setup generated files` | 上次的 `lightdock_*.pdb` 还在 | 删掉作业目录重跑，或换 `--output-dir` |
| `NoneType has no attribute` | 结构里有非标准残基 | 加 `--keep-water` 试试；或用 `pdb2pqr` 预处理 |

**注意**：`already exists` 这条错误来自 LightDock 源码：

```python
if moved_file_name.exists():
    raise LightDockError(
        f"{moved_file_name} already exists, please delete previous setup generated files"
    )
```

所以**同一个作业目录不能重复跑 setup**，要么先清理，要么换目录。

---

## 模拟阶段跑很久 / 看起来卡住

**症状**：`lightdock3.py` 长时间没有输出。

**原因**：这是正常的。LightDock 的日志在**每个 swarm 结束时**才打印一行
（`[lightdock] INFO: GSO cluster 0 - Minion 5`）。50 个 swarm 要等
**全部**跑完才有完整输出。

**处理**

- 看进程：`top` / `htop`，确认 CPU 真的在跑（LightDock 是纯 CPU 计算）。
- 看输出目录：`ls swarm_*/gso_*.out | wc -l`，**文件数量就是进度**。
- 中途想停：直接 Ctrl-C。已写入的 GSO 文件不丢，下次用 `--resume` 续跑。

**加速建议**

```bash
# 用更快的评分函数
python scripts/run_docking.py ... --run-exe /path/to/lightdock3.py
# 在 setup.json 里改 scoring_function 为 fastdfire

# 减少 swarm 数（cores 不超过 swarms）
--swarms 20 --cores 8
```

---

## 模拟阶段报 `ModuleNotFoundError: No module named 'lightdock.scoring.xxx'`

**症状**

```
ModuleNotFoundError: No module named 'lightdock.scoring.xyz'
```

**原因**：setup.json 里指定的评分函数名不存在。

**处理**：合法的评分函数名是 `dfire`、`fastdfire`、`vdw`、`pisa`、`cpydock`。
检查 `setup.json` 里的 `scoring_function` 字段。

---

## 分数全是 0 或提取不到

**症状**：`docking_scores.csv` 里 `best_score` 为空，或者所有作业的
`best_glowworm` 都是 0。

**诊断步骤**

```bash
# 1. 看 GSO 文件是否存在
ls results/rec__lig/swarm_0/

# 2. 看文件内容
head -3 results/rec__lig/swarm_0/gso_100.out
```

正常应该看到：

```
(31.4171143, 1.8570079, ...) 0 0 11.25395618 0 4.200 7.52800101
```

**如果 `best_glowworm` 恒为 0**：说明用了错误版本的解析器（把 LigID 当成了
glowworm 编号）。本技能的 `parse_gso_line()` 已修正——glowworm 编号来自
`enumerate()` 的行号。

**如果文件是空的**：说明模拟没跑完就被中断了，重跑即可（GSO 文件是增量写入的）。

**如果打分是正数**：DFIRE 打分应为负。正数说明取错了列——检查是否误取了
`Luciferin`（那个是正的）。

---

## 生成的复合物看起来是错的

**症状**：打开 `lightdock_0.pdb`，配体飞到很远的地方，或者和受体完全分离。

**原因**：几乎肯定是**用了原始结构而不是 setup 平移后的结构**。

**处理**

```bash
# 检查这两个文件存在吗
ls results/rec__lig/lightdock_receptor.pdb results/rec__lig/lightdock_ligand.pdb
```

不存在说明准备阶段没成功。存在就要确认生成命令用的是它们而不是
`receptor.pdb` / `ligand.pdb`。

本技能的 `generate_conformations()` 里：

```python
receptor = moved_structure_path(output_dir / "receptor.pdb")
ligand = moved_structure_path(output_dir / "ligand.pdb")
```

并且会检查文件存在，缺失时明确报错而不是静默用错文件。

---

## 想给配体换链 ID

**注意**：**不要**原地修改输入文件。老仓库的做法是"备份 → 原地改写 → 还原"，
靠 `finally` 兜底；一旦中途异常退出，配体文件就永久处于被改坏的状态，
后续所有对接都会用错误的输入。

本技能的做法是**生成一份副本**：

```bash
python scripts/generate_conformations.py --job-dir results/rec__lig/ --ligand-chain B
```

原文件不动，副本用完即删。

---

## 磁盘空间不够

**症状**：`No space left on device`。

**原因**：生成的构象 PDB 很大。每个 `lightdock_<i>.pdb` 含受体 + 配体的全部原子，
一个 500 残基的复合物约 2–5 MB。100 个 glowworm × 50 个 swarm 就是几 GB。

**处理**

- 只生成需要的构象：`--top 1` 而不是 `--top 100`。
- 批量生成后清理：`extract_scores.py` 之后删掉 `swarm_*/lightdock_*.pdb`，
  保留 GSO 输出（那才是原始数据）。
- 先跑 `--swarms 5 --glowworms 10 --steps 20` 估算规模再上生产参数。

---

## 作业跑了一半想中断

**直接 Ctrl-C 即可。** 本技能的结果是**每个作业完成就增量写**到 CSV 的
（`handle.flush()`），中断不会丢失已完成作业的结果。

续跑：

```bash
python scripts/run_docking.py --receptor rec.pdb --ligands ligands/*.pdb \
    --output-dir results/ --resume
```

`--resume` 会跳过已有有效结果的作业。

**注意**：中断的作业目录里可能留有半截的 `lightdock_*.pdb`（准备阶段产物）。
重跑前需要清掉，否则 setup 阶段会报 `already exists`。

---

## 内存不足

**症状**：`MemoryError` 或进程被 OOM killer 杀掉。

**处理**

- 减少 `--glowworms`（每个 glowworm 都在内存里有一份坐标）。
- 减少 `--cores`（每个进程是一份独立的 swarm 副本）。
- 用 `fastdfire` 或 `vdw` 评分函数（CPyDock 需要 freesasa，内存开销最大）。
- 加 `--keep-hydrogens` 是反方向，会更占内存，通常应该保留默认的移除选项。

---

## 怎么确认结果是不是真的最优

1. **看分数分布**。把所有作业的分数排开序，如果最高分和最低分几乎一样
   （比如都在 -30 附近），说明搜索没充分展开——加 swarms/steps。
2. **看 swarm 之间的差异**。`extract_scores.py --per-swarm` 给出每个 swarm
   的最优分。如果不同 swarm 之间差异巨大，说明某个 swarm 塌缩了。
3. **看收敛曲线**。`gso_0.out` 到 `gso_100.out` 的最优分应该单调下降。
   解析器之外，也可以直接看文件最后一行。
4. **结构合理性检查**。用 `protein-interactions` 技能（PLIP）看界面残基是否
   合理，用 `protein-structure` 技能做二级结构标注。
