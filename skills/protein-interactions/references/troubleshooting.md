# 排错指南：protein-interactions

## 1. 解析类问题

### 1.1 `0` 条相互作用 / 记录数为零

**症状**：JSON 里 `interactions: []`，但 `report.xml` 明显有内容。

**原因与对策**：

1. **选错了文件**。PLIP 输出目录里除了 `report.xml` 还有 `report.html`、`report.txt`。
   确认传入的是 `report.xml`。
2. **schema 不匹配**。先用文本编辑器打开确认根标签结构：
   * 属性式 → `<site chain="A" resnr="102" ...>`；
   * 子元素式 → `<site><salt_bridge><resnr>102</resnr>...`。
   两者本技能都支持；若都不匹配，说明文件被截断或经过了自定义改写。
3. **报告是"空报告"**。若结构中没有可识别的配体（无 HETATM 或配体残基名不在
   PLIP 的常见配体库中），PLIP 会正常返回但报告里没有任何 pocket。

**验证**：把 `report.xml` 里的前几十行贴出来数一下 `<pocket>` 元素个数。

### 1.2 残基号显示为 `?` 或节点名异常

**症状**：Mermaid 图里出现 `1?`、`LIG?` 之类的节点。

**原因**：配体标签（`LIG:1:CL1`）被当成蛋白标签（`chain:resname+resnr:atom`）解析。

**对策**：本技能用"中间段是否含 `+`"判定，`parse_residue_label` 返回的
`is_ligand` 字段应为 `True`。若发现不成立，检查 `plip_parser.py` 是否被改动。

### 1.3 盐桥出现 ARG–ARG 这类不合理配对

**症状**：`Salt Bridge` 类型里出现同号电荷残基组合。

**原因**：PLIP 报告本身没错，但你把 **protein–protein** 报告当成了 protein–ligand
报告解读，或 receptor/ligand 角色被混淆。

**对策**：查看记录的 `receptor_residue` / `ligand_residue` 字段确认角色；
在 protein–protein 场景下，两侧都是蛋白残基属正常。

## 2. PLIP 运行类问题

### 2.1 找不到 `plip` 可执行文件

**症状**：报 `plip` 不可用。

**对策**：

```bash
conda install -c conda-forge plip     # 或 pip install plip
```

或用 `--plip-executable /path/to/plip` 显式指定。

### 2.2 PLIP 运行超时或崩溃

**症状**：超时退出（默认 600 s）。

**对策**：

* `--plip-timeout 1800` 放宽限制（大复合物首次运行较慢）；
* 确认结构文件已加氢（PLIP 需自行补氢时会显著变慢甚至不收敛）；
* PLIP 对超大结构（> 10⁵ 原子）会内存爆炸，先截取结合位点附近的残基。

### 2.3 报告里氢键数量明显偏少

**原因**：结构缺氢。X 射线结构通常没有 H 原子。

**对策**：用 `pdb2pqr` 或 `reduce` 加氢：

```bash
pdb2pqr --ff=AMBER --with-ph=7.4 1crn.pdb 1crn_pH.pdb
```

## 3. 统计类问题

### 3.1 热门残基全部只有 1 次相互作用

**症状**：每个残基都只参与 1 个相互作用。

**原因**：受体 / 配体角色被互换，或只统计了配体侧。

**对策**：JSON 中 `hot_residues.receptor` 与 `hot_residues.ligand` 分别给出；
protein–ligand 体系中结合热点通常在受体侧出现 ≥ 3 次。

### 3.2 相互作用总数比 PLIP HTML 报告少

**原因**：解析器只处理 8 类已知的相互作用标签。PLIP 报告中还有
`hbond_lig` / `hbond_prot` / `halogen` 等细分类型，未覆盖的会被跳过。

**对策**：查看 `summarize_interactions` 的 `skipped_types` 字段（若有），
或在 `plip_parser.py` 的 `PLIP_INTERACTION_TYPES` 中补充新标签。

## 4. AlphaFold 相关

### 4.1 zip 包里找不到 `*_model_N.cif`

**症状**：选不到模型。

**原因**：包结构与预期不同（某些版本把模型放在 `*.cif` 而置信度 JSON 名字不同）。

**对策**：

```bash
unzip -l AF-P69905-F1-model_v4.zip
```

确认实际文件名，再调整 `alphafold_models.py` 的匹配模式。

### 4.2 `mean_plddt` 很低（< 70）

**原因**：模型本身质量差（无序区、错误折叠），或选错了模型。

**对策**：用 `pTM + 0.5·ipTM` 排序确认是否选到了最佳模型；
若所有模型 pLDDT 都低，说明该靶点不适合 AlphaFold 建模，不要用它做相互作用分析。

### 4.3 ipTM 很低但 pTM 很高

**原因**：单链折叠可靠，多链界面不可靠。

**对策**：不要把该模型的界面相互作用当结论，只用单体内部数据。

## 5. 可视化类问题

### 5.1 Mermaid 图不渲染 / 报错

**症状**：Markdown 预览中图形不显示。

**原因**：节点 id 或标签含非法字符（`:`、`+`、`(`）。

**对策**：节点 id 已做归一化（`L_A_ARG_102`）。若自行扩展了解析逻辑，
务必保持"节点 id 只用字母数字下划线，标签用引号包裹"。

### 5.2 所有边都是同一种颜色

**症状**：`linkStyle default` 生效导致统一配色。

**对策**：删除 `linkStyle default` 那行，改为逐条 `linkStyle <index> stroke:<color>`。

### 5.3 PyMOL 脚本执行报语法错误

**症状**：`syntax error` 在某一行。

**原因**：残基号带插入码（如 `102A`）时 `resi 102A` 需要写成 `resi 102A`
（PyMOL 支持），但若字符串拼接产生了空残基号则会失败。

**对策**：确认选中集合非空；可用 `print count(sel)` 检查。

## 6. 自检清单

```bash
# 1) 帮助可用
python skills/protein-interactions/scripts/analyze_interactions.py --help

# 2) 解析器自检（内置样例报告）
python skills/protein-interactions/scripts/plip_parser.py

# 3) 全流程 + 图生成
python skills/protein-interactions/scripts/analyze_interactions.py \
    --input work/plip_test/report.xml --output work/out_inter \
    --make-mermaid --make-pymol

# 4) 规范校验
python scripts/validate_skills.py
```
