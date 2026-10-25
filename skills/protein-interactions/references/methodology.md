# 方法学：PLIP 报告解析与相互作用分析

## 1. 数据来源与前提

| 输入 | 说明 |
|---|---|
| PLIP `report.xml` | PLIP 输出目录（`<job>/report.xml`）中的 XML 报告 |
| AlphaFold `*.zip` | AlphaFold Protein Structure Database 下载包，内含 `*_model_N.cif` 与 `*_summary_confidences_N.json` |
| 结构文件（PDB） | 仅在 `--run-plip` 时使用，需先由 PLIP 处理 |

输入结构必须**带氢原子**。PLIP 的氢键、盐桥、疏水接触判据全部依赖氢坐标；
若结构由 X 射线实验得到且未加氢，PLIP 会跳过大量接触，输出会显著偏少。

## 2. PLIP XML 两种 schema

### 2.1 新版（属性式，PLIP ≥ 2020）

```xml
<pocket name="pocket1" residues_count="12">
  <site chain="A" resnr="102" restype="ARG" onf="A:ARG+102:NE">
    <salt_bridge dist="2.9" onp="LIG:1:CL1" acc="A:ARG+102:NE" don="LIG:1:CL1"/>
  </site>
  <site chain="A" resnr="93" restype="TYR" onf="A:TYR+93:OH">
    <pi_stack dist="3.6" onp="LIG:1:CL1" group1="..." group2="..."/>
  </site>
</pocket>
```

* 残基身份写在 `<site>` 的属性上：`chain` / `resnr` / `restype` / `onf`。
* 单条相互作用是 `<site>` 的子元素，标签名即相互作用类型（单数形式）。
* 参与原子通过 `onp`（配体侧）与 `acc` / `don`（受体侧）给出。
* `dist` 为 PLIP 计算的接触距离（Å）。

### 2.2 旧版（子元素式）

```xml
<pocket name="pocket1">
  <site>
    <salt_bridge>
      <resnr>102</resnr>
      <restype>ARG</restype>
      <reschain>A</reschain>
      <restype_lig>LIG</restype_lig>
      <resnr_lig>1</resnr_lig>
      <dist>2.9</dist>
    </salt_bridge>
  </site>
</pocket>
```

* 残基字段写在子元素里；配体侧字段带 `_lig` 后缀。
* 解析器对两种 schema 同时兼容：先看 `<site>` 是否有属性，有则走新版分支；
  否则读取子元素走旧版分支。

## 3. 残基标签解析

PLIP 用 `onf` / `onp` 给出带原子名的残基标签，格式为
`chain:resname+resnr:atom`（蛋白）或 `resname:resnr:atom`（配体 / 杂原子）。

区分依据是**中间段是否包含 `+`**：

| 标签 | 中间段 | 判定 |
|---|---|---|
| `A:ARG+102:NE` | `ARG+102` | 蛋白残基：chain=`A`, resname=`ARG`, resnr=`102` |
| `LIG:1:CL1` | `1` | 配体原子：chain=`LIG`, resname=`LIG`, resnr=`1`, is_ligand=True |

若把配体标签误判为蛋白残基，网络图节点会显示成 `1?`（残基号缺失）——
这是本技能早期版本修掉的解析 bug，现在以 `is_ligand` 字段显式区分。

## 4. 相互作用类型与物理判据

技能实现 8 类 PLIP 相互作用的抽取：

| 类型 | XML 标签 | PLIP 判据摘要 |
|---|---|---|
| 盐桥 | `salt_bridges` | 带正电侧链（ARG/LYS/HIS）与带负电侧链（ASP/GLU）距离 < 5 Å，且有可成键的 N/O 原子 |
| 疏水接触 | `hydrophobic_interactions` | 两个疏水残基（ALA/VAL/LEU/ILE/PHE/MET/PRO/TRP/TYR 等）侧链碳原子距离 < 4.5 Å |
| 氢键 | `hydrogen_bonds` | 供体 D–H 与受体 A 距离 < 3.5 Å，且 D–H···A 角度 > 120° |
| π-阳离子 | `pi_cation_interactions` | 芳环质心与带电基团距离 < 6 Å，质心到最近芳香原子距离在环平面法向投影合理范围内 |
| π-堆积 | `pi_stacks` | 两个芳环环心距离 < 5.5 Å，环平面夹角 < 30°（平行）或 > 150°（交叉 T 型） |
| 水桥 | `water_bridges` | 两个水分子各自与供体/受体成键，且水–水距离 < 3.5 Å |
| 卤键 | `halogen_bonds` | X···受体距离 < 3.5 Å 且 C–X···受体 角度 > 120° |
| 金属复合 | `metal_complexes` | 金属离子与配位残基（ASP/GLU/HIS/CYS）距离 < 3.0 Å |

技能本身**不重新计算**这些几何判据——它忠实抽取 PLIP 已完成的判定结果，
保证与 PLIP 官方结果一致；判据细节列在上表仅供解释与排错。

## 5. 统计汇总

对抽取出的记录集合计算：

1. **类型计数** — 每种相互作用的出现次数。
2. **受体侧热门残基** — 按残基聚合参与次数，取 top-N。
3. **配体侧热门残基** — 同上，标注为 ligand。
4. **链级分布** — 按残基所属链统计交互数，识别主要结合面。

网络图节点以残基为单位（不是原子单位），因此同一残基参与多种相互作用时，
图中是多条边指向同一个节点——这正好保留了"一个热点残基接触多个配体原子"的信息。

## 6. Mermaid 网络图渲染

* 语法：`flowchart LR`，按 receptor / ligand 分成两个 `subgraph`。
* 节点 id 做归一化（`L_` + `链_残基名_残基号`），避免 Mermaid 保留字与特殊字符问题。
* 边标签为人类可读的相互作用类型名。
* 着色：先输出全部边，再用 `linkStyle <index> stroke:<color>,stroke-width:1.5px;`
  按边序号逐条指定颜色。Mermaid 的 `linkStyle default` 会把所有边设成同一颜色，
  无法表达"按类型着色"，因此必须逐条写。
* 每种类型一个固定色（`INTERACTION_COLORS`），保证多次运行结果配色一致。

## 7. PyMOL 标注脚本

生成的 `.pml` 包含：

* `select <name>, resi <resnr> and chain <chain>` —— 选出参与相互作用的残基；
* `show sticks, <name>` / `color <type-color>, <name>` —— 棒状显示并按类型着色；
* `label <name>, "<resname><resnr> (<type>)"` —— 标注残基与相互作用类型；
* `distance ...` —— 在配体与受体关键残基之间画距离标尺。

脚本是纯文本命令序列，可直接 `pymol -r annotation.pml` 批处理运行。

## 8. AlphaFold 模型选择

AlphaFold DB 下载包通常包含多个模型（recycled / BFE 变体）。选择规则：

```
score = pTM + 0.5 × ipTM
```

`pTM` 反映整体折叠可信度，`ipTM` 反映多链界面的可信度；对单体蛋白两者接近，
对复合物 ipTM 是决定界面是否可信的关键。`pTM + 0.5·ipTM` 是常见的加权折中。

置信度摘要输出：`ptm`、`iptm`、`mean_plddt`、`max_pae`、
`high_confidence_fraction`（pLDDT > 90 的残基比例）。

pLDDT 解读：> 90 高可信（骨架位置明确）；70-90 可接受（骨架大致正确）；
50-70 低（可能错误折叠）；< 70 很高（基本不可用）。

## 9. 参数选择

| 场景 | 建议 |
|---|---|
| 只想看统计 | `--format csv` |
| 用于报告 / 幻灯片 | `--make-mermaid`（Markdown 渲染器可直接内嵌） |
| 需要结构可视化 | `--make-pymol`，在 PyMOL 中 `run annotation.pml` |
| AlphaFold 复合物 | 先看 `iptm`，ipTM < 0.6 说明界面预测不可靠 |
| 大批量报告 | 按目录传入 `--input`，脚本递归收集所有 `report.xml` |
