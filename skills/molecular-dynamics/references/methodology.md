# 方法学：分子动力学轨迹分析

## 1. 轨迹格式与内置解析器

### 1.1 DCD：Fortran 风格定长记录流

DCD（CHARMM/NAMD trajectory）不是"固定 84 字节一帧"的简单格式，而是
**Fortran unformatted record 流**：每条记录都由 4 字节长度标记 + 数据 + 4 字节
长度标记组成。文件结构：

```
[header 记录]           84 字节 = 21 个 int32
[标题记录]              可选，CHARMM 版，长度是 80 的整数倍
[natoms 记录]           1 个 int32
[自由原子列表记录]      可选，长度 = nfixed * 4
重复 n_frames 次:
    [X 记录][Y 记录][Z 记录]     每个 natoms * 4 字节（float32）
    [unit cell 记录]             可选，48 字节（6 个 float64）
```

**header 记录（21 个 int32）的语义**：

| 索引 | 含义 | 索引 | 含义 |
|---|---|---|---|
| 0 | 记录长度（84 或 164） | 9-10 | DELTA（float32，存放于两个 int32 槽） |
| 1 | NSET（帧数） | 11 | 有约束的原子数 |
| 2 | NSAVC（写帧间隔） | 12-20 | 备用 / 单元胞标志 |
| 3 | NSTEP（总步数） | | |
| 8 | NAMNF（固定原子数） | | |

**关键设计决策**：不同发行版（NAMD / CHARMM / X-PLOR / VMD）对 header 内字段
的偏移并不统一。本实现**不依赖 header 字段定位 natoms**，而是利用记录的
**自描述性**做嗅探：

* 长度为 4 的记录 → natoms 记录；
* 长度为 `natoms * 4` 的记录 → X / Y / Z 单轴坐标块（每帧 3 条）；
* 长度为 `natoms * 3 * 4` 的记录 → 单记录帧（内部按 `[x…, y…, z…]` 排列）；
* 长度为 48 的记录 → unit cell；
* 其余 4 的倍数长度的记录 → 辅助记录（标题、自由原子列表），直接跳过。

两条坐标布局都被支持，因此 GROMACS 经 `convert` 转换出的 DCD 与
CHARMM 直接写出的 DCD 都能正确读取。

### 1.2 PDB 多模型轨迹

`MODEL` / `ENDMDL` 之间是一帧；各帧原子数必须一致才能堆叠成轨迹数组。
单模型 PDB 也可读取，此时轨迹长度为 1 帧。

### 1.3 XYZ 轨迹

标准 XYZ 格式：第 1 行原子数，第 2 行注释，之后每行 `元素 x y z`。
逐帧堆叠，用于快速测试与简单可视化。

### 1.4 读取优先级

`load_trajectory` 先尝试 MDAnalysis（若已安装且能解析），失败则回退到内置
解析器。返回的 `Trajectory` 对象统一提供 `coords` `(n_frames, n_atoms, 3)`、
`topology`、`times` 与 `reader` 字段。

## 2. Kabsch-Sander 叠合

### 2.1 推导

给定待叠合坐标 `P`（mobile）与参考坐标 `Q`（reference），求旋转 `R` 与
平移 `t` 使 `Σ wᵢ |R pᵢ + t − qᵢ|²` 最小。

1. 加权质心：`c_P = Σ wᵢ pᵢ / Σ wᵢ`，`c_Q = Σ wᵢ qᵢ / Σ wᵢ`；
2. 去质心：`p̃ᵢ = pᵢ − c_P`，`q̃ᵢ = qᵢ − c_Q`；
3. 协方差矩阵：`C = Σ wᵢ p̃ᵢ q̃ᵢᵀ`（3×3）；
4. SVD：`C = U S Vᵀ`；
5. 旋转：`R = V D Uᵀ`，其中 `D = diag(1, 1, det(V Uᵀ))`。

`D` 的第三项保证 `det(R) = +1`，即结果是**纯旋转**而非"旋转 + 镜像"。
这是 Kabsch 算法中经常被遗漏的一步——漏掉它会让 RMSD 在某些构型下偏大。

6. 平移：`t = c_Q − R c_P`；
7. RMSD：`√(Σ wᵢ |R pᵢ + t − qᵢ|² / Σ wᵢ)`。

本实现用 `numpy.linalg.svd` 完成步骤 4，`numpy.linalg.det` 求行列式。

### 2.2 RMSD 的两种用法

| 用法 | 是否叠合 | 含义 |
|---|---|---|
| 叠合后 RMSD | 是 | 内部涨落（internal RMSD / RMSF 的前提） |
| 原始 RMSD | 否 | 含整体平移与旋转，用于检测整体漂移 |

逐残基 RMSD 的做法：先把所有帧叠合到参考帧，再对每个 Cα 计算其相对参考
位置的位移范数。这与"先叠合再逐原子算 RMSF"等价，但能保留逐残基的空间分布。

## 3. RMSF

均方根涨落：

```
RMSF_i = sqrt( (1/T) Σ_t |x_i(t) − x̄_i|² )
```

关键在于 `x̄_i` 是**叠合后**轨迹的平均位置。若不对齐，全局漂移会污染
RMSF（表现为所有原子涨落巨大）。本实现在 `rmsf(coords, align_reference=...)`
中先做逐帧 Kabsch 叠合再统计。

## 4. 回旋半径

质量加权回旋半径：

```
Rg = sqrt( Σᵢ wᵢ |rᵢ − r_cm|² / Σᵢ wᵢ ),   r_cm = Σᵢ wᵢ rᵢ / Σᵢ wᵢ
```

元素质量取自标准原子量表（H=1.008, C=12.011, N=14.007, O=15.999, S=32.06 …）。
未知元素默认按碳处理。

经验标度：`Rg ≈ 2.2 × N^0.38` Å（N 为残基数），可用于快速判断折叠状态。

## 5. PCA

1. 取所有帧作为样本，展平为 `X` `(n_frames, 3·n_atoms)`；
2. 去均值：`X_c = X − mean(X)`；
3. 协方差矩阵 `C = X_cᵀ X_c / (n_frames − 1)`（3n×3n）；
4. `eigh(C)` 得到特征值（升序）与特征向量；
5. 按特征值降序重排，取前 k 个主成分；
6. 投影：`P = X_c @ V_k`（每个帧在主成分空间的坐标）；
7. 解释方差比例：`λᵢ / Σλ`。

**前置要求**：必须先叠合。否则整体平移/旋转会占据 PC1，
使结果失去构象意义（表现为 PC1 解释 >95% 方差）。入口脚本在调用 PCA 前
自动完成逐帧叠合。

## 6. 层次聚类

1. 对所有帧对计算 Kabsch 叠合后的 RMSD，得到 `n×n` 对称矩阵
   （帧数超过 `--max-cluster-frames` 时均匀抽样以控制 O(N²) 计算量）；
2. `scipy.cluster.hierarchy.linkage(squareform(matrix), method='average')`；
3. `fcluster(tree, t=k, criterion='maxclust')` 得到簇标签
   （给定 `--cluster-threshold` 时用 `criterion='distance'`）。

linkage 方法选择：

| 方法 | 特点 | 适用 |
|---|---|---|
| `average` | 簇间平均距离 | 默认，稳健 |
| `complete` | 簇间最大距离 | 对紧凑簇更敏感 |
| `single` | 链状效应 | 容易产生"单点大簇" |
| `ward` | 类方差最小化 | 只对欧氏距离有意义 |

scipy 不可用时退化为**阈值连通分量**划分（按中位数距离建图取连通分量），
结果不如 linkage 稳健但仍可用。

## 7. 氢键：距离 + 角度几何判据

本实现采用经典的几何判据（不计算部分电荷）：

对供体重原子 D（元素为 N/O/S）、其上的氢 H、受体 A（元素为 N/O/S/F）：

1. `d(H···A) ≤ 2.5 Å`（默认）；
2. `d(D···A) ≤ 3.5 Å`（默认）；
3. `∠D-H···A ≥ 120°`（默认）。

角度定义在 H 上：`v1 = D − H`，`v2 = A − H`，`angle = arccos(v1·v2 / (|v1||v2|))`。

**拓扑筛选**：用第一帧确定"哪个 H 属于哪个 D"（距离 < 1.35 Å 且索引相邻），
再枚举候选三元组，避免每帧重新搜索导致的性能问题与误配。

**排除规则**：同一残基内的 D···A 对被排除（否则每个残基内部必然出现
"自氢键"）；可通过 `exclude_pairs` 进一步排除 1-4 键等。

**存在率（occupancy）**：某氢键在多少比例的帧中存在。占据率 > 0.5 的氢键
通常认为是"稳定氢键"。

**重要前提**：轨迹必须含氢原子。X 射线结构一般没有 H，此时氢键数必然为 0。
解决：用 `pdb2pqr` 或 `reduce` 加氢。

## 8. 二面角与 Ramachandran

### 8.1 二面角计算

采用 IUPAC 约定的投影法（与 `praxeolitic` 的公式逐位一致）：

```
b0 = p0 − p1;  b1 = p2 − p1;  b2 = p3 − p2
axis = b1 / |b1|
v = b0 − (b0·axis) axis       # b0 在垂直于 b1 平面上的投影
w = b2 − (b2·axis) axis       # b2 的投影
x = v·w;  cr = axis × v;  y = cr·w
dihedral = atan2(y, x)
```

符号约定必须与蛋白化学主链定义一致：α-螺旋 φ ≈ −57°、ψ ≈ −47°；
β-折叠 φ ≈ −139°、ψ ≈ +135°。符号写反会让整张 Ramachandran 图镜像。

主链二面角定义：

* `φ(i) = C(i−1) − N(i) − CA(i) − C(i)`
* `ψ(i) = N(i) − CA(i) − C(i) − N(i+1)`

### 8.2 区域划分

按 (φ, ψ) 的粗略区域分为四类：

| 区域 | φ 范围 | ψ 范围 |
|---|---|---|
| α 区 | −160 … −20 | −90 … 0 |
| α 区（右螺旋） | −120 … 50 | 60 … 180 |
| β 区 | −180 … −50 | 60 … 180 |
| αL 区 | −90 … 30 | −60 … 100 |
| 其余 | — | coil |

区域判定优先 α → β → αL，最后归入 coil。**离群点**定义为既不在 α 也不在
β 也不在 αL 区域的采样点。

## 9. 链间相互作用

对每一对不同链，逐帧计算所有跨链原子对的距离；距离 < `--contact-dist`
（默认 4.5 Å）记为一次接触。统计每对原子接触的出现率，只报告
`occupancy ≥ 0.5` 的稳定接触。

典型解读：受体–配体界面的稳定接触数在轨迹中保持不变，说明界面未解离；
若接触数在后半段下降，说明复合物正在解离。

## 10. SASA（Shrake-Rupley）

1. 每个原子的范德华半径取自元素表（H 1.20, C 1.70, N 1.55, O 1.52,
   S 1.80, P 1.80 …），探针半径默认 1.4 Å（水）。
2. 在原子球面上按**黄金螺旋**均匀撒 `n_points` 个点
   （默认 100）：`φᵢ = arccos(1 − 2(i+0.5)/n)`，`θᵢ = π(1+√5)(i+0.5)`。
3. 对每个采样点，检查它是否落在其它原子的"扩展球"（半径 + 探针半径）内；
   被覆盖的点不可及。
4. 该原子的 SASA = `4π r² × (不可及点数 / 总点数)`。

结果再分解为极性（N/O）与非极性贡献。典型数值：100 残基球状蛋白的总
SASA 约 5000-7000 Å²。

## 11. OpenMM 模拟

`run_simulation.py` 在 OpenMM 可用时真实执行：

1. `app.PDBFile` 读入结构；
2. 力场：默认 `amber14-all.xml` + `amber14/tip3pfb.xml`，隐式溶剂追加
   `app.GBSAOBCForceGenerator(water_model="tip3p")`；
3. `createSystem(nonbondedMethod=NoCutoff, constraints=HBonds, rigidWater=True)`；
4. `LangevinMiddleIntegrator(300 K, 1/ps, 0.002 ps)`；
5. 最小化 → 平衡 → 生产性 MD；
6. 输出 `trajectory.dcd`（DCDReporter）、`state.csv`（StateDataReporter）、
   `final_state.pdb`。

时间步长选择：刚性氢键约束下 2 fs 是标准做法；4 fs 需要约束 C-H 键，
本技能默认使用 2 fs。

## 12. 分阶段平衡协议

### 12.1 为什么要分阶段

复合物模拟最常见的失败不是"算不动"，而是**平衡阶段就散架**：

- **界面未松弛** —— 对接得到的姿态里往往有原子重叠、侧链穿插。
  直接进入生产期，能量会在前几百皮秒里剧烈释放，轨迹前段全是瞬态。
- **晶格记忆** —— 晶体结构的堆积方式在真空/隐式溶剂下没有物理意义，
  蛋白会朝晶格方向继续"结晶"。NVT 阶段给侧链松弛的时间，
  再进 NPT 让密度回到合理值，能显著改善这一问题。
- **链分离** —— 界面主要由疏水作用维持时，300 K 的热涨落足以在
  平衡期把两条链拉开。位置限制是标准对策。

### 12.2 四个阶段

| 阶段 | 参数 | 物理意义 |
|---|---|---|
| 最小化 | `--minimize-steps` | 消除初始结构中的原子重叠 |
| NVT | `--nvt-steps` | 定容，逐步升温到目标温度 |
| NPT | `--npt-steps` | 定压，密度松弛到真实值 |
| 生产期 | `--steps` | 采样真实动力学 |

温度不用"逐步升温"也能到位：先用 `setVelocitiesToTemperature` 一次性
给速度，然后靠 Langevin 恒温器维持。NVT 阶段的意义不在于升温曲线，
而在于给体系一段**无压力耦合的松弛时间**。

### 12.3 位置限制的力学

势能形式（各向同性谐振子）：

```
U = 0.5 * k * ((x - x0)^2 + (y - y0)^2 + (z - z0)^2)
```

- `k` 是力常数（kJ/mol/Å²）；
- `(x0, y0, z0)` 是参考位置，默认取输入结构中的坐标。

在 `k = 100 kJ/mol/Å²` 时，把原子拉离参考位置 1 Å 需要约 50 kJ/mol，
与一个典型氢键强度相当；拉离 2 Å 则需要 200 kJ/mol，限制足够强。
生产期若要保留限制，通常降到 `k = 10 kJ/mol/Å²` —— 只压制大幅摆动，
不干扰正常涨落。

选择模式：

| 模式 | 目标原子 | 典型用途 |
|---|---|---|
| `backbone` | 主链 N / CA / C | **默认推荐**，让侧链与溶剂松弛 |
| `ca` | 仅 Cα | 更弱，只固定整体构象 |
| `heavy` | 所有非氢原子 | 强限制，用于极不稳定的初始姿态 |
| `all` | 所有原子 | 等价于"冻结" |

### 12.4 一个必须知道的 OpenMM 陷阱

**修改 `System` 不会重建已有的 `Context`。**

OpenMM 的 `Context` 在创建时把 `System` 的力场"编译"进计算核。之后
再调 `system.addForce(...)` 或 `system.removeForce(...)`，对正在跑的
`Context` **没有任何影响** —— 力不会生效，也不会报错，只是静默失效。

所以：

- 加位置限制 → 必须在 `app.Simulation(...)` **之前**加到 `system`；
- 加气压计（`MonteCarloBarostat`）→ 加完后**重建** Simulation；
- 撤除位置限制 → `removeForce` 后**重建** Simulation。

重建时必须显式搬运状态，否则体系会跳回初始构象：

```python
state = simulation.context.getState(getPositions=True, getVelocities=True)
simulation = app.Simulation(topology, system, integrator, platform)
simulation.context.setPositions(state.getPositions())
simulation.context.setVelocities(state.getVelocities())
```

如果不搬速度，重建后的体系速度为零，等于在生产期开头做了一次
隐式的"零动能重启"，前几百步的行为完全失真。

### 12.5 链分离监测

用各链**质心间距**作为解离指标：

```
d(A,B) = |centroid_A - centroid_B|
```

- 判定阈值默认 30 Å。这个数字的来历是：两个各约 300 残基的球状
  蛋白，其质心间距在结合态通常在 20-30 Å；超过 30 Å 基本可以认为
  界面已经断开。
- 同时监测**相对初始值的变化量**，因为某些体系（如链间有共价交联
  或巨大配体）的绝对距离本来就大，只看绝对值会误报。
- 链 ID 取自 OpenMM `Chain.id`（即 PDB 文件里的真实链标识），
  **不是** `chr(65 + chain.index)` —— PDB 里链标识可以是任意字符，
  按 index 猜出来的 A/B/C 往往和文件里的实际标识对不上。

### 12.6 数值稳定性回滚

生产期结束后检查末帧坐标：

- 出现 `NaN` → 通常是力场参数缺失、氢原子没加、或初始结构有重叠；
- 出现 `Inf` → 能量爆炸的终态；
- 坐标量级超过 `1e6 Å` → 数值发散。

三者任一成立就回滚到输入结构，并把原因写进报告的 `stability.reason`。
回滚目标本身也要检查是否健全 —— 回滚到一个同样坏的结构没有意义。

## 13. 参数选择速查

| 场景 | 建议 |
|---|---|
| 长时间轨迹（>10⁵ 帧） | `--stride 10` 抽稀，`--max-cluster-frames 200` |
| 只关心柔性 | `--selection "name CA"` |
| 关注结合界面 | `--selection "resid 1-40"` 配合 `--analysis interchain` |
| 界面稳定性 | `--analysis interchain --contact-dist 4.0` |
| 水合氢键 | `--selection "resname HOH"` 配合 `--analysis hbond` |
| 采样不足 | 加大 `--steps`；RMSD 平台期未出现说明轨迹未收敛 |
