---
name: pymol-structure-visualization
description: This skill should be used when analyzing inter-chain contacts in protein structures — computing residue-level contact geometry from PDB files with a uniform-grid spatial index, classifying each contact into hydrogen bonds, salt bridges, hydrophobic contacts, pi-pi stacking and cation-pi interactions using residue chemistry, identifying hotspot residues, and exporting CSV/JSON results together with Mermaid interaction networks and ready-to-run PyMOL visualization scripts. 用于 PDB 链间残基接触分析、结合界面热点残基鉴定、相互作用网络构建与 PyMOL 可视化脚本生成的场景。
license: MIT
compatibility: Python 3.9+；核心分析仅需标准库（零第三方依赖）；PyMOL 为可选的渲染端
metadata:
  version: "1.0"
  last-reviewed: "2026-10-26"
  skill-author: sinder
---

## Overview

从 PDB 结构文件里算出**链间残基接触**，并按化学类型分类：

```
PDB ──解析──► 原子列表 ──网格空间划分──► 残基对接触（最小距离）
                                              │
                                     残基化学性质分类
                                              │
              ┌───────────────┬───────────────┼───────────────┐
           CSV 明细      JSON 完整结果    Mermaid 网络图    PyMOL 脚本
```

核心技术主张：**接触类型不是靠几何距离猜的，而是靠参与接触的两侧残基的化学性质组合判定的**。距离只起"截断"作用，真正决定类型的是供体/受体、电荷、芳香性、疏水性这些化学属性。

技术上从老仓库 `pymol_show/` 移植而来（22 个脚本），修掉了 6 个真实缺陷 —— 最核心的两个是：

- **老实现报告的是"两个残基各自靠近了对方"，不是"这两个残基互相接触了"**。它用 `cmd.select(f"{sel1} within {cutoff} of {sel2}")` 选出残基再做笛卡尔积，一个供体残基和一个受体残基只要各自内部有原子在截断距离内就算一对，与两残基之间是否真有原子对接触无关。
- **必须启动 PyMOL GUI**（`pymol.finish_launching(['pymol'])`），导致无 PyMOL 的机器上完全跑不了、CI 里无法测试。

本技能把几何计算与 PyMOL 渲染彻底分离：**PyMOL 脚本是生成物，不是运行依赖**。

## When to Use This Skill

- 有一个蛋白复合物的 PDB，想知道**两条链之间靠什么结合**（氢键？盐桥？疏水？）
- 需要鉴定**结合界面的热点残基** —— 参与相互作用最多的那几个残基往往就是突变实验的候选位点
- 想画**相互作用网络图**，放进论文的补充材料或汇报 slides
- 需要**可复现**的接触分析脚本，而不是在 PyMOL GUI 里手动点选
- 想为已有的对接结果准备可视化脚本，但**手头没有装 PyMOL**（本技能照样能出脚本）

**不适用**：

- 需要精确的**氢键几何判据**（供受体原子 + 氢原子角度）→ 用 `protein-interactions`（PLIP）
- 需要**对接打分** → 用 `protein-docking-hdock` / `protein-docking-lightdock`
- 需要**结构质量评估**（Ramachandran、Twiggle、steric clash）→ 用 `protein-structure`

## Quick Start

```bash
export PYTHONPATH="$PWD/packages"

# 1) 基本用法：分析两条链的接触
python skills/pymol-structure-visualization/scripts/analyze_contacts.py \
    --input complex.pdb --chains A B --output results/

# 2) 所有链两两分析 + 生成 PyMOL 可视化脚本 + Mermaid 网络图
python skills/pymol-structure-visualization/scripts/analyze_contacts.py \
    --input complex.pdb --output results/ \
    --make-mermaid --make-pymol --make-pymol-contacts

# 3) 收紧距离阈值（只看强接触）
python skills/pymol-structure-visualization/scripts/analyze_contacts.py \
    --input complex.pdb --output results/ --max-distance 4.0

# 4) 排除配体与水
python skills/pymol-structure-visualization/scripts/analyze_contacts.py \
    --input complex.pdb --output results/ --exclude-hetero
```

## Workflow

1. **准备输入** —— 确认 PDB 文件可读、包含目标链
   - 多 MODEL 文件（NMR、轨迹）只取**第一个 MODEL**
   - 链标识从文件里读，不要假设是 `A`/`B`
2. **质控** —— 确认原子数、链数、残基数合理
   - 结构里没有氢是正常的：氢键判定基于**残基化学性质**，不需要氢坐标
   - 水与配体默认参与分析（`--exclude-hetero` 可排除）
3. **核心分析** —— 网格空间划分 + 残基化学分类
   - 网格边长取 `max_distance × 1.0001`，扫 3×3×3 邻域
   - 每个接触类型有独立的距离阈值（见下表）
4. **结果汇总** —— 输出 CSV / JSON，含按类型统计与热点残基
5. **验证与出图** —— 生成 Mermaid 网络图与 PyMOL 脚本
6. **产出报告** —— 结合领域知识判断哪些接触是关键的

## Parameters

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `--input` | Path | 必填 | 输入 PDB 结构文件 |
| `--output` | Path | `<输入名>_contacts/` | 输出目录 |
| `--chains` | str+ | 全部链两两 | 要分析的链对（恰好 2 个），如 `A B` |
| `--max-distance` | float | `6.0` | 最大接触距离（Å） |
| `--exclude-hetero` | flag | off | 忽略 HETATM（水、离子、配体） |
| `--include-hydrogen` | flag | off | 把氢原子也算进去（默认忽略） |
| `--include-histidine` | flag | off | 把 HIS 当作正电残基 |
| `--include-ambiguous-donors` | flag | off | 把 ASP/GLU 也当作弱氢键供体 |
| `--min-distance` | float | `0.0` | 只报告距离不小于此值的接触 |
| `--make-mermaid` | flag | off | 生成 Mermaid 交互网络图 |
| `--make-pymol` | flag | off | 生成 PyMOL 残基标注脚本 |
| `--make-pymol-contacts` | flag | off | 生成基于 `cmd.distance` 的接触绘制脚本 |
| `--json` | Path | `<输出目录>/contacts.json` | 完整结果 JSON 路径 |
| `-v/--verbose` | flag | off | 输出 DEBUG 级日志 |
| `-q/--quiet` | flag | off | 仅输出警告与错误 |

## 接触类型与判据

| 类型 | 中文 | 判据 | 距离上限 |
|------|------|------|----------|
| `Hydrogen_bond` | 氢键 | 一侧供体残基 + 另一侧受体残基 | 3.5 Å |
| `Salt_bridge` | 盐桥 | 正电残基（LYS/ARG）+ 负电残基（ASP/GLU） | 4.0 Å |
| `Hydrophobic` | 疏水接触 | 两侧均为疏水残基 | 4.5 Å |
| `Pi_Pi` | π-π 堆积 | 两侧均为芳香族残基 | 5.0 Å |
| `Cation_Pi` | 阳离子-π | 正电残基 + 芳香族残基 | 6.0 Å |
| `Ligand_contact` | 配体接触 | 任一侧为非标准残基（配体等） | 由 `--max-distance` 决定 |

**一个接触对可以有多个类型**（例如 ARG–ASP 同时是氢键和盐桥）。`interaction_types` 字段给出全部类型，`interaction_type` 字段给出优先级最高的那个。

### 三处与老实现不同、需要注意的地方

1. **ASP/GLU 不是氢键供体**。羧基 pKa ≈ 4，生理 pH 下几乎完全去质子化。老实现把 ASP/GLU 列进供体表，会**虚增氢键数**。需要时可开 `--include-ambiguous-donors`。
2. **HIS 默认不算正电残基**。咪唑基 pKa ≈ 6.0，生理 pH 下仅约 10% 质子化。老实现直接把它算进正电列表，会**虚增盐桥数**。需要时可开 `--include-histidine`。
3. **配体不套用残基级分类表**。我们不知道配体是芳香还是极性、是否带电，硬套会得出无意义的结论。所以配体接触单独标为 `Ligand_contact`，距离照常报告 —— 具体类型请用 `protein-interactions`（PLIP）。

## Outputs

| 文件 | 说明 |
|------|------|
| `<输出目录>/contacts.csv` | 接触明细（UTF-8 BOM，Excel 可直接打开） |
| `<输出目录>/contact_summary.csv` | 按类型的计数统计 |
| `<输出目录>/contacts.json` | 完整结果（接触 + 汇总 + 规则描述） |
| `<输出目录>/contact_network.mmd` | Mermaid 交互网络图（`--make-mermaid`） |
| `<输出目录>/visualize_contacts.pml` | PyMOL 残基标注脚本（`--make-pymol`） |
| `<输出目录>/draw_contacts.pml` | PyMOL 接触虚线绘制脚本（`--make-pymol-contacts`） |

### contacts.csv 列说明

| 列 | 说明 |
|----|------|
| `Chain_Pair` | 链对（如 `A-B`） |
| `ChainA` / `ResnA` / `ResiA` | 第一侧链 / 残基名 / 残基号 |
| `ChainB` / `ResnB` / `ResiB` | 第二侧链 / 残基名 / 残基号 |
| `Min_Distance` | 两侧最近原子间距（Å） |
| `N_Atom_Pairs` | 构成该接触的原子对数量 |
| `Interaction_Type` | 主要类型（优先级最高的） |
| `All_Types` | 全部类型，分号分隔 |

## Validation

结果可靠性检查清单：

- [ ] **接触对数量合理** —— 一个典型蛋白-蛋白界面的氢键数在 5-30 量级；突然出现几百个通常意味着 `--max-distance` 设得过大
- [ ] **热点残基符合结构生物学直觉** —— 结合界面的热点残基往往落在已知的活性位点或保守位点上
- [ ] **距离分布合理** —— 氢键应该在 2.6-3.2 Å，盐桥 2.6-3.5 Å；出现大量 < 2.0 Å 的接触说明结构有立体冲突，或坐标精度有问题
- [ ] **抽查一个接触的原子对** —— 在 PyMOL 里用 `cmd.distance` 复核 CSV 里记录的距离
- [ ] **配体接触已识别** —— 如果结构里有配体，`by_type` 里应出现 `Ligand_contact`
- [ ] **多类型接触合理** —— 同一个残基对同时报出盐桥 + 氢键是正常的（ARG–ASP 胍基–羧基的典型情况）

**注意**：本技能输出的是辅助分析结果。距离阈值与化学分类是**通用近似**，不是精确的力场判据；发表级别的相互作用判定请用 PLIP。

## References

- [方法学原理](references/methodology.md) —— 判据来源、网格算法、化学分类表的取舍理由、修掉的 6 个缺陷
- [常见问题与排错](references/troubleshooting.md) —— 故障现象、根因与解决方案
