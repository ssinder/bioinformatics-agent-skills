# PyMOL 结构可视化与链间接触分析

从 PDB 结构文件分析链间残基接触，按化学类型分类（氢键 / 盐桥 / 疏水 /
π-π 堆积 / 阳离子-π），输出 CSV、JSON、Mermaid 交互网络图与可直接运行的
PyMOL 可视化脚本。

**核心分析零第三方依赖** —— PyMOL 只是可选的渲染端，脚本是生成物而非运行依赖。

## 快速开始

```bash
export PYTHONPATH="$PWD/packages"

# 基本用法
python skills/pymol-structure-visualization/scripts/analyze_contacts.py \
    --input complex.pdb --chains A B --output results/

# 全链对分析 + Mermaid 网络图 + PyMOL 脚本
python skills/pymol-structure-visualization/scripts/analyze_contacts.py \
    --input complex.pdb --output results/ \
    --make-mermaid --make-pymol --make-pymol-contacts
```

## 脚本

| 脚本 | 作用 |
|------|------|
| `scripts/analyze_contacts.py` | 主入口：解析 → 接触检测 → 分类 → 输出 |
| `scripts/contact_geometry.py` | 几何引擎：均匀网格空间划分 + 残基对接触检测 |
| `scripts/contact_chemistry.py` | 知识核心：残基化学分类表与接触类型判定 |

## 接触类型

| 类型 | 中文 | 判据 | 距离上限 |
|------|------|------|----------|
| `Hydrogen_bond` | 氢键 | 供体残基 + 受体残基 | 3.5 Å |
| `Salt_bridge` | 盐桥 | 正电残基 + 负电残基 | 4.0 Å |
| `Hydrophobic` | 疏水接触 | 两侧均为疏水残基 | 4.5 Å |
| `Pi_Pi` | π-π 堆积 | 两侧均为芳香族残基 | 5.0 Å |
| `Cation_Pi` | 阳离子-π | 正电残基 + 芳香族残基 | 6.0 Å |
| `Ligand_contact` | 配体接触 | 任一侧为非标准残基 | 随 `--max-distance` |

一个接触对可有多个类型（如 ARG–ASP 同时是氢键与盐桥）。

完整说明见 [SKILL.md](SKILL.md)。

- 方法学原理：[references/methodology.md](references/methodology.md)
- 排错手册：[references/troubleshooting.md](references/troubleshooting.md)
