# 常见问题与排错

本文按"症状 → 原因 → 处理"的格式整理 `protein-structure` 技能的实际故障模式。

## 1. 解析类问题

### 症状：`未解析到任何原子记录`

**原因**：输入文件不是 PDB 格式，或使用了列偏移不同的格式。

**排查**：

```bash
head -5 your_file.pdb        # 前几行应含 ATOM/HETATM
awk 'substr($0,1,4)=="ATOM"' your_file.pdb | wc -l
```

若 `ATOM` 记录存在但仍解析失败，检查文件是否为 **mmCIF**（`.cif` / `.pdbx`）格式。
mmCIF 使用 `ATOM` 关键字而非 PDB 的固定列，本技能不解析 mmCIF。

**处理**：先转换成 PDB 格式（`gemmi convert`、`pdb2pqr` 或
`Bio.PDB.MMCIFParser` → `PDBIO`）。

### 症状：残基数远少于预期

**原因**：交替位置过滤过于激进，或文件只含 `MODEL` 记录而无 `ATOM`。

**排查**：统计不同 altLoc 的原子数：

```bash
awk '/^ATOM/ {print substr($0,17,1)}' your_file.pdb | sort | uniq -c
```

**处理**：本技能默认保留主构象（altLoc 为空/`A`/`1`）。若结构只有 `B` 构象，
需修改 `pdb_parser.py` 中 `_parse_atom_line` 之后的 altLoc 过滤条件。

### 症状：链被拆成多段，或链 ID 混乱

**原因**：文件使用 `TER` 记录分段但复用同一 chainID，或采用 `segid` 而非 `chainID`。

**排查**：

```bash
grep -c "^TER" your_file.pdb
awk '/^ATOM/ {print substr($0,22,1)}' your_file.pdb | sort -u
```

**处理**：本技能按 `(chain_id, res_seq, icode)` 唯一标识残基，因此同一 chainID
下重复的残基号会被合并。遇到 `TER` 分段且 chainID 复用时，建议先用
`pdb_selaltloc` / `pdb_split` 之类工具清理文件。

### 症状：元素统计里缺少某些原子 / 元素识别错误

**原因**：`element` 列（76-78）为空，且原子名推断失败。

**处理**：本技能的推断规则（首字符为数字时取第二位）覆盖绝大多数情况。
遇到特殊命名（如金属配位的 `HG` 表示汞）时，用 `awk` 检查：

```bash
awk '/^HETATM/ {print substr($0,77,2)}' your_file.pdb | sort -u
```

### 症状：序列里出现 `X`

**原因**：非标准残基（HETATM 修饰残基、MSI、SEC 等）。

**处理**：属于预期行为。本技能把 MSE 映射为 M，其余非标准残基为 X。
若需完整序列，应输出 `HETATM` 修饰残基。

## 2. 二级结构类问题

### 症状：内置 DSSP 报告的螺旋含量为 0%，但结构明显是螺旋

**原因**（按可能性排序）：

1. **结构缺少氢原子且重建失败**——若残基缺少 `N` 或 `CA`，无法重建氢。
2. **链被 `TER` 打断**——`pos` 序列在片段边界不连续，`i → i+4` 的判定错位。
3. **坐标是 nm 而非 Å**——Cα-Cα 距离变成 0.38 Å，氢键距离判据全部失效。
4. **阈值过严**——`--hbond-energy -0.5` 对某些低分辨率结构太严格。

**排查**：

```bash
# 检查 Cα-Cα 距离（应为 ~3.8 Å）
python - <<'EOF'
import sys; sys.path.insert(0, "skills/protein-structure/scripts")
from pdb_parser import parse_pdb_file
import math
s = parse_pdb_file("your_file.pdb")
for chain in s.model.chains:
    ca = [r.get("CA") for r in chain if r.get("CA")]
    ds = [math.dist(a.coord, b.coord) for a, b in zip(ca, ca[1:])]
    print(chain.chain_id, "CA-CA mean", round(sum(ds)/len(ds), 2) if ds else "N/A")
EOF
```

若均值在 0.38 左右，说明单位是 nm（OpenMM 默认 nm 输出）。

**处理**：

- 单位是 nm → 统一乘 10 转为 Å
- 阈值太严 → `--hbond-energy -0.2`
- 链被打断 → 检查 `TER` 记录

### 症状：内置 DSSP 把 β-折叠股判成 `T`（转角）或 `-`

**原因**：折叠股内部的判定依赖 β 区的 φ/ψ。若 φ/ψ 计算有偏差，或结构经过
低分辨率精修导致二面角偏离，就会落到转角/卷曲分支。

**排查**：查看输出的 `phi` / `psi` 字段。理想 β-折叠股应为 φ ≈ -139°、ψ ≈ +135°。

**处理**：使用外部 DSSP（`--dssp-method external`）得到权威结果。

### 症状：外部 DSSP 与内置实现结果差异很大

**常见原因**：

1. 外部 DSSP 用了不同的能量模型与更完整的氢键网络计算（含水桥）
2. 外部 DSSP 使用文件中的显式氢原子（若有），而内置实现重建了虚拟氢
3. 内置实现把 β-桥简化为需要相邻双氢键，DSSP 对孤立接触也计为桥

**判断是否可接受**：同一结构上两者的**总二级结构分布**（螺旋 + 折叠 + 卷曲）
通常应相差 5-10% 以内。超过 20% 说明结构异常或阈值需要调整。

### 症状：`mkdssp` 找不到

**原因**：未安装，或不在 PATH 中。

**排查**：

```bash
which mkdssp || which dssp
mkdssp --version
```

**处理**：

```bash
conda install -c bioconda dssp
# 或指定路径
python scripts/analyze_structure.py --input x.pdb --output out/ \
    --dssp-method external --dssp-executable /opt/conda/bin/mkdssp
```

注意新版 `mkdssp` 的参数是位置参数 + `-o`；旧版 `dssp` 用 `-o` 输出文件。
本技能两种命令形式都尝试。

### 症状：外部 DSSP 超时

**原因**：大结构（>5000 残基）+ 慢速存储。

**处理**：增大超时不可配置时，建议拆分为单结构处理，或直接用内置实现
（纯 Python 但无需子进程启动开销，小结构反而更快）。

## 3. 几何计算类问题

### 症状：Rg 与文献值对不上

**排查清单**：

1. **定义不同**：全原子 Rg vs Cα-only Rg，两者可差 5-10%。输出中 `rg` 与
   `rg_ca` 都给出了，请确认引用的是哪个。
2. **氢原子处理**：全原子包含 H 的结构其 Rg 略大于只含重原子的结构。
3. **是否已对齐**：Rg 是**平移与旋转不变**的量，但前提是只统计同一结构的
   同一组原子。若把不同模型或不同选择（`protein` vs `name CA`）混在一起，
   结果没有意义。
4. **质量加权 vs 等权**：本技能是质量加权的。

### 症状：`shape_class` 总是 `sphere` 或总是 `coil`

**原因**：判据阈值（`b/Rg² < 0.01`、`c/Rg² < 0.05`）偏保守。

**处理**：查看原始描述符 `asphericity`、`acylindricity`、
`relative_shape_anisotropy`（κ）。κ > 0.4 已是明显各向异性，
即使分类落在 `coil` 也说明形状并不球形。

### 症状：偶极矩数值明显偏大

**原因**：元素级简化电荷对金属离子（Zn²⁺ 取 1.30、Fe 取 1.80 等）
与阴离子残基（Asp/Glu 取 -0.42）会累积出很大的贡献。

**处理**：该数值只用于**同一体系不同构象之间**的相对比较。
需要绝对值时请用带电荷的力场拓扑。

## 4. UniProt 查询类问题

### 症状：`未安装 requests，HTTP 功能不可用`

**处理**：`pip install requests`。注意不加 `--query-uniprot` 时**完全不联网**，
脚本不会因为缺少 requests 而失败。

### 症状：查询全部失败 / 返回 0 条

**排查**：

```bash
curl -s "https://rest.uniprot.org/uniprotkb/P69905.json?fields=accession" | head -c 200
```

- 网络不通 → 检查代理 / 防火墙
- 返回 404 → 文件名不是 UniProt accession（PDB 的四位代码如 `1CRN`
  不是 accession）。此时应改用 `--dssp-method` 无关的元数据来源，
  或手工提供 accession 列表
- 返回 429（限流）→ 共享层已内置重试与 0.34 s 节流；若持续出现，
  减少单批数量

### 症状：部分 accession 有结果、部分没有

**原因**：UniProt 中这些条目可能是次要条目（TrEMBL 的冗余条目已被合并）、
未被注释、或已被删除。

**处理**：查询结果会少于提交数量，这是正常现象。日志中会给出实际返回条数。

## 5. 输出与批量问题

### 症状：批量处理报"没有成功分析任何结构"

**排查**：逐个文件手动测试，找出共同失败原因：

```bash
for f in dir/*.pdb; do
  python scripts/analyze_structure.py --input "$f" --output /tmp/out -v 2>&1 | tail -2
done
```

### 症状：输出目录被写到了意料之外的位置

**原因**：未传 `--output` 时，默认取 `biopipe_core.get_config()` 的
`output_dir`，其根目录来自 `BIOPIPE_PROJECT_ROOT` 环境变量或当前工作目录。

**处理**：显式传 `--output`，或设置环境变量：

```bash
export BIOPIPE_PROJECT_ROOT=/path/to/project
```

### 症状：批量运行很慢

**优化建议**：

1. 减小 `--max-files` 先做小规模验证
2. 用 `--dssp-method builtin` 避免每帧启动外部进程
3. 提高 `--workers`（线程池在纯 Python 计算下加速有限，但 IO 密集时有效）

## 6. 自检清单

提交结果前确认：

- [ ] `python scripts/analyze_structure.py --help` 正常退出
- [ ] 对目标文件 `CA(i)-CA(i+1)` 距离在 3.6-4.0 Å
- [ ] `Cα(i)-CA(i+4)` 在 α-螺旋结构上约 6.0-6.4 Å
- [ ] 螺旋结构上 `O(i)···N(i+4)` 约 2.9-3.1 Å
- [ ] 内置与外部 DSSP（若都可用）的螺旋含量相差 < 15%
- [ ] 残基数与 PDB 文件头部（`REMARK 3` 或 SEQRES）大致相符
- [ ] 输出 JSON 中的 `n_atoms` 与 `grep -c "^ATOM"` 的结果一致（若含 HETATM 则相应更大）
