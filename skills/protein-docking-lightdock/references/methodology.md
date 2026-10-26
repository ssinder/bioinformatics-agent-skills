# LightDock 方法学

本文记录 LightDock 的输出语义、坐标系处理与评分提取的原理。所有内容依据
LightDock 0.9.4 源码（`lightdock/util/analysis.py`、`lightdock/prep/simulation.py`、
`lightdock/util/parser.py`、`bin/` 下的脚本）与官方教程核对。

---

## 1. LightDock 的计算流程

```
receptor.pdb + ligand.pdb
        │
        ▼
lightdock3_setup.py            ← 预处理：去水/氧/氢、算椭球、放 swarm、写 setup.json
        │
        ├── lightdock_receptor.pdb   ← 平移到原点后的受体
        ├── lightdock_ligand.pdb     ← 平移到原点后的配体
        ├── setup.json
        ├── swarm_0/ ... swarm_N/
        │
        ▼
lightdock3.py setup.json N -c C ← 模拟：GSO 群体优化
        │
        ├── swarm_0/gso_0.out, gso_10.out, gso_20.out, ...
        │
        ▼
lgd_generate_conformations.py   ← 生成复合物 PDB
        │
        └── swarm_0/lightdock_0.pdb, lightdock_1.pdb, ...
```

**LightDock 与 HDOCK 的根本区别**：HDOCK 是网络服务，计算在服务器上；LightDock
是本地程序，计算在用户机器上，且**必须自己准备环境**（conda / venv + 编译好的
Cython 扩展）。所以脚本设计的重点是"长时间、可中断、可恢复"的批处理，而不是
"提交完就拿结果"。

---

## 2. GSO 输出文件的真实格式

`swarm_<n>/gso_<step>.out` 里**每行一个 glowworm**。官方
`lightdock/util/analysis.py` 的 `read_lightdock_output()` 给出的字段顺序是：

```
(坐标) RecID LigID Luciferin Neighbors VisionRange Scoring
```

官方示例（来自 lightdock.org/tutorials/2UUY）：

```
(31.4171143, 1.8570079, -6.3956223, -0.1058407, -0.4849369, 0.5997430, -0.6276482) 0 0 11.25395618 0 4.200 7.52800101
```

各字段含义：

| 字段 | 含义 | 说明 |
|---|---|---|
| 坐标括号 | `(tx, ty, tz, qw, qx, qy, qz)` | 平移 + 四元数。开启 ANM 时会追加法向模态系数 |
| `RecID` | 受体结构索引 | 单体对接恒为 0；多结构对接时指向输入列表中的第几个 |
| `LigID` | 配体结构索引 | 同上 |
| `Luciferin` | 萤光素值 | GSO 算法里与打分成正比，是"吸引力"标量 |
| `Neighbors` | 视野内邻居数 | 反映构象聚集程度 |
| `VisionRange` | 当前视野半径 | GSO 的自适应搜索半径 |
| `Scoring` | 打分函数值 | 默认 DFIRE。**越负越好** |

### 2.1 陷阱一：glowworm 编号是行号

官方代码里写得很直白：

```python
# lightdock/util/analysis.py
id_glowworm=id_line,              # ← 行号
receptor_id=int(rest[0]),         # ← RecID
ligand_id=int(rest[1]),           # ← LigID
```

**glowworm 编号既不是 `RecID` 也不是 `LigID`，而是该行在文件中的行号（从 0 起）。**

老仓库 `dock_lightdock/scripts/` 里的正则：

```python
pattern = r'\([^\)]+\)\s+\d+\s+(\d+)\s+[\d.-]+\s+\d+\s+[\d.]+\s+([\d.-]+)'
glowworm_idx = int(match.group(1))  # 注释写的是 "LigID"
```

捕获组 1 是 **LigID**——单体对接时恒为 `0`。于是"最佳构象"永远是
glowworm #0，与真实得分最高的 glowworm 毫无关系。

**这个 bug 的恶劣之处在于它不产生任何报错**：文件正常生成、分数正常输出、
构象正常写出，只是**挑出来的结构从来不是最优的那个**。

### 2.2 陷阱二：GSO 文件要按步数数值挑

LightDock **每 10 步存一次**结果，文件名不补零：

```
gso_0.out  gso_10.out  gso_100.out  gso_110.out  gso_20.out  gso_30.out ...
```

字典序下 `gso_100.out` 排在 `gso_20.out` **前面**。所以 `sorted(glob(...))[-1]`
取到的"最后一个"文件，在 100 步的运行里可能是 `gso_90.out` 甚至更早——
**用了一个还没收敛的中间结果**。

正确做法是解析文件名里的数字取最大：

```python
best = max(candidates, key=lambda p: int(re.search(r"gso_(\d+)\.out$", p.name).group(1)))
```

### 2.3 陷阱三：`Scoring` 越负越好

DFIRE 打分是能量相关的量，**越负表示能量越低、构象越好**。所以"最优"是
`min()` 而不是 `max()`。

老代码写的是 `if score < best_score`，方向是对的——但因为捕获组取的是
LigID 而非 scoring，实际比较的是 LigID 的数值，这个"正确的方向"根本没作用在
正确的字段上。

---

## 3. 坐标系：为什么必须用平移后的结构

`lightdock3_setup.py` 在准备阶段做了这件事：

```python
# lightdock/structure/complex.py
def move_to_origin(self):
    """Moves the structure to the origin of coordinates"""
    translation = [-1 * c for c in self.center_of_coordinates()]
    self.translate(translation)
    return translation
```

然后把平移后的结构另存（`lightdock/prep/simulation.py`）：

```python
moved_file_name = Path(file_name).parent / Path(DEFAULT_LIGHTDOCK_PREFIX % Path(file_name).name)
# DEFAULT_LIGHTDOCK_PREFIX = "lightdock_%s"
```

所以 `receptor.pdb` → `lightdock_receptor.pdb`，`ligand.pdb` → `lightdock_ligand.pdb`。

**GSO 文件里的平移量 `tx, ty, tz` 正是相对这个原点记录的。**

生成构象时官方用法是：

```bash
cd swarm_0
lgd_generate_conformations.py ../lightdock_receptor.pdb ../lightdock_ligand.pdb gso_100.out 10
```

老仓库 `generate_conformations.py` 传的是：

```python
cmd = [
    str(gen_conf_exe),
    "../receptor.pdb",      # ← 原始文件，坐标系不对
    "../ligand.pdb",        # ← 原始文件，坐标系不对
    Path(best_gso).name,
    "1",
]
```

配体因此被平移了两次（一次到原点、一次由 GSO 记录），复合物的两个部分
**根本对不上**。生成的 PDB 里配体在离谱的位置，或者与受体完全脱离接触。

这类错误同样**不会报错**，只会让结果"看起来像那么回事"。

---

## 4. 构象文件命名

`lgd_generate_conformations.py` 的输出：

```python
write_pdb_to_file(receptor, os.path.join(destination_path, "lightdock_%s.pdb" % i), receptor_pose)
write_pdb_to_file(ligand,  os.path.join(destination_path, "lightdock_%s.pdb" % i), ligand_pose)
```

**i 从 0 开始，不补零。** 所以字典序：

```
lightdock_0.pdb  lightdock_1.pdb  lightdock_10.pdb  lightdock_11.pdb  lightdock_2.pdb  ...
```

`i` 同时就是**该构象对应的 glowworm 编号**（= GSO 文件里的行号）。所以：

- 要拿"最优构象"，必须先知道最优 glowworm 的编号，再去找 `lightdock_<编号>.pdb`；
- 不能靠 `sorted(...)[0]`（碰巧能拿到 0 号，但不是"最优"）；
- 不能靠 `sorted(...)[-1]`（拿到字典序最后的，等于编号最大的，不是最优）。

本技能的 `sort_conformation_files()` 按**数值**排序，`_pick_best_file()` 直接按
编号匹配。

---

## 5. 官方命令速查（0.9.4）

这些脚本通过 `setup.py` 的 `scripts=[...]` 声明安装，是**裸文件**而非
console_scripts 入口点——所以 `pip show lightdock` 的 entry point 列表里看不到
它们，但它们确实落在 `<env>/bin/` 下。

| 命令 | 用途 | 关键参数 |
|---|---|---|
| `lightdock3_setup.py` | 准备 | `receptor.pdb ligand.pdb [-s N] [-g N] [--noh] [--noxt] [--now] [--seed_points N] [-anm]` |
| `lightdock3.py` | 模拟 | `setup.json <steps> [-c cores] [-s scoring_fn] [-sg seed] [-f config]` |
| `lgd_generate_conformations.py` | 生成构象 | `receptor.pdb ligand.pdb gso_<step>.out <n> [--setup setup.json]` |
| `lgd_rank.py` | 排名 | `<num_swarms> <steps> [-c clashes_cutoff] [--ignore_clusters]` |
| `lgd_top.py` | 提取 top 模型 | `receptor.pdb ligand.pdb <ranking_file> <top> [--setup setup.json]` |
| `lgd_generate_trajectory.py` | 生成轨迹 | 用于动画 |
| `lgd_map_contacts.py` | 接触图 | — |
| `lgd_gso_to_csv.py` | GSO 转 CSV | — |

### 参数默认值（`lightdock/constants.py`）

```python
DEFAULT_NUM_SWARMS = 400
DEFAULT_NUM_GLOWWORMS = 200
DEFAULT_LIGHTDOCK_PREFIX = "lightdock_%s"
DEFAULT_REC_NM_FILE = "lightdock_rec.nm"
DEFAULT_LIG_NM_FILE = "lightdock_lig.nm"
```

注意官方默认是 **400 swarms × 200 glowworms**——那是生产级配置，单个作业在
普通笔记本上要跑很久。老仓库用的 50×50 是"快速配置"，方向是对的。

### 评分函数可选值

`dfire`（默认）、`fastdfire`、`vdw`、`pisa`、`cpydock`。速度排序：
`vdw` > `fastdfire` > `dfire` > `pisa` > `cpydock`。**调试阶段用 `fastdfire`
可以显著提速**，正式结果用 `dfire`。

---

## 6. 参数选择建议

| 场景 | swarms | glowworms | steps | 说明 |
|---|---|---|---|---|
| 流程验证 | 5 | 10 | 20 | 几分钟内跑完，确认环境没问题 |
| 常规筛选 | 20 | 30 | 50 | 老仓库的"优化版"配置，平衡速度与覆盖 |
| 标准对接 | 50 | 50 | 50 | 老仓库 `dock_simple.py` 的配置 |
| 精细优化 | 100 | 100 | 100 | 接近官方默认值，耗时成倍增长 |

**核数选择**：LightDock 用多进程跑不同 swarm，所以 `cores` 不应超过 swarms
数（多余的核空转）。官方文档建议 `optimal_cores = min(available_cores, num_swarms)`。

---

## 7. 限速与礼貌

LightDock 是本地计算，没有服务端限速的问题。但要注意：

- **不要一口气开几十个并行作业**。每个作业会开 `cores` 个进程，
  50 个作业 × 8 核 = 400 个进程，会把机器拖垮。
- **磁盘占用可观**。每个 swarm 目录有 GSO 输出（不大），但生成构象后
  `lightdock_*.pdb` 每个几十 MB。100 个 glowworm 的构象就是几 GB。
- 老仓库的 `dock_completess.py` 会把 `design×completess` 的**笛卡尔积**全部
  跑一遍。如果有 10 个 design 和 500 个配体，那是 **5000 个作业**。跑之前
  先用 `--swarms 5 --glowworms 10 --steps 20` 在两三个组合上验证。

---

## 8. 从老仓库继承并修掉的问题清单

| # | 问题 | 后果 | 修法 |
|---|---|---|---|
| 1 | 把 LigID 当 glowworm 索引 | **挑出的结构永远不是最优** | 按行号编号 |
| 2 | `sorted(glob)[-1]` 取最新 GSO | 用未收敛的中间结果 | 按步数数值取最大 |
| 3 | 生成构象用原始 PDB | 复合物坐标系错配 | 用 `lightdock_*.pdb` |
| 4 | `sorted(glob)[0]` 取"最佳"构象 | 多构象时取错 | 按编号数值排序 + 精确匹配 |
| 5 | 硬编码 `/Users/jung/...` 路径 | 换机即失效 | PATH / 前缀自动探测 + CLI 参数 |
| 6 | 三份逐字重复的脚本 | 改一处漏两处 | 合并为单一实现 |
| 7 | 配体链 ID 原地改写 + 备份还原 | 中途崩溃留下损坏文件 | 改副本，原文件不动 |
| 8 | 每 500 作业才写一次 CSV | 中断丢失全部结果 | 每个作业增量写 |
| 9 | 单作业失败中止整个批次 | 一次失败毁掉全部 | 单作业失败不影响其余 |
| 10 | `design×ligand` 笛卡尔积无预演 | 直接开 5000 个作业 | 支持小规模快速试跑参数 |
