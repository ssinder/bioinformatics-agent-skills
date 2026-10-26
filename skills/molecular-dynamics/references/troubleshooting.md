# 排错指南：molecular-dynamics

## 1. 轨迹读取问题

### 1.1 `不是有效的 DCD 文件（header 长度 X）`

**症状**：header 记录的首个 int32 不是 84 或 164。

**原因**：

* 文件被截断或损坏；
* 文件其实是其它格式（如 `.xtc`、`.trr`、`.nc`）；
* 文件是 GROMACS 的 `.dcd` 但经过了文本转换。

**对策**：

```bash
xxd traj.dcd | head -3      # 确认首 4 字节是 0x54 0x00 0x00 0x00
file traj.dcd
```

GROMACS 原生格式需先转换：

```bash
gmx trjconv -s md.tpr -f md.trr -o traj.dcd
```

### 1.2 `DCD 文件中未找到合法的 natoms 记录`

**症状**：嗅探不到 natoms。

**原因**：文件头被破坏，或使用了带扩展头的私有格式。

**对策**：用 MDAnalysis 交叉验证：

```python
import MDAnalysis as mda
u = mda.Universe("traj.dcd")
print(u.atoms.n_atoms, u.trajectory.n_frames)
```

若 MDAnalysis 能读，说明本实现的嗅探需要调整（记录长度判别条件）。

### 1.3 `DCD 记录长度不匹配: 头 X != 尾 Y`

**症状**：文件在读取中途报长度不一致。

**原因**：文件损坏，或存在本实现未覆盖的记录布局。

**对策**：检查文件大小是否为 4 的倍数；用 `mdcheck`（MDAnalysis 自带）
验证：

```bash
python -m MDAnalysis.tests.datafiles.mdcheck traj.dcd
```

### 1.4 原子数与拓扑不匹配

**症状**：`拓扑原子数 (X) 与轨迹原子数 (Y) 不一致`。

**原因**：DCD 与 PDB 不是同一次模拟的产物，或 PDB 经过了后处理
（加氢、去水、去配体）而 DCD 没有。

**对策**：使用与 DCD 同时输出的拓扑文件；确认 PDB 的原子顺序与 DCD 一致
（同一份 topology 用 `gmx trjconv` 抽取，不要手工改）。

### 1.5 PDB 轨迹解析出 1 帧

**症状**：`n_frames = 1`。

**原因**：PDB 文件没有 `MODEL`/`ENDMDL` 标记（单模型结构）。

**对策**：这是正常行为——单模型 PDB 就是一帧。若需要多帧，请提供 DCD。

### 1.6 轨迹帧数比预期少

**症状**：读到 40 帧，但预期 1000 帧。

**原因**：DCD 的 NSAVC（写帧间隔）大于 1，或文件被截断。

**对策**：`dcd_header_info()` 返回 `nsavc` 与 `n_frames`；真实时间轴为
`time = frame_index × nsavc × delta`。

## 2. 分析结果异常

### 2.1 RMSD 全为 0 或恒定不变

**症状**：所有帧 RMSD 相同。

**原因**：轨迹各帧坐标完全相同（模拟未推进），或叠合算法退化为恒等变换。

**对策**：

* 用 `--analysis rmsf` 检查是否有涨落；
* 用 `trajectory_io.py` 打印坐标范围，确认不同帧坐标有差异；
* 检查 DCD 的 NSAVC 是否导致读到的都是同一帧。

### 2.2 RMSD 持续线性增大到很大

**症状**：RMSD 单调上升到 10 Å 以上。

**原因**：蛋白在解折叠（这是真实信号，不是 bug）。

**对策**：结合 Rg（是否增大）、SASA（是否增大）、氢键（是否减少）判断。
若 RMSD 大而 Rg 与 SASA 稳定，多半是柔性环区在动——看逐残基 RMSD 定位。

### 2.3 RMSD 数值大但结构明显没变

**症状**：结构叠合良好，RMSD 却很高。

**原因**：选择了含大量溶剂的原子，或参考帧选的是异常帧。

**对策**：

* `--selection "name CA"` 或 `"protein"`；
* 换 `--fit-reference` 到 RMSD 中位数附近的那一帧；
* 检查是否有整体平移未叠合（用 `--analysis rmsd` 的原始模式对比）。

### 2.4 RMSF 全部接近 0

**症状**：RMSF < 0.01 Å。

**原因**：轨迹未采样到真实涨落，或读到的帧重复。

**对策**：确认 `--stride` 没有把有效帧全过滤掉；确认轨迹时长足够
（至少 ns 量级）。

### 2.5 RMSF 全部很大（>5 Å）

**症状**：所有原子 RMSF 都很高。

**原因**：

1. 采样不足（未平衡）；
2. 力场问题；
3. 温度过高。

**对策**：延长平衡时间；检查 `state.csv` 中的温度是否稳定在目标值 ±5 K；
检查能量是否有漂移。

### 2.6 氢键数为 0

**症状**：`mean_count = 0`，但拓扑里明明有 N/O。

**原因**：

1. **结构没有氢原子**（最常见）——X 射线结构通常只有重原子；
2. 选区把所有 N/O 排除了；
3. 判据过严。

**对策**：

```bash
pdb2pqr --ff=AMBER --with-ph=7.4 top.pdb top_pH.pdb
```

加氢后重新分析。也可放宽判据：

```bash
--hbond-h-dist 2.8 --hbond-da-dist 4.0 --hbond-angle 100
```

### 2.7 氢键数过多（>500）

**症状**：每帧氢键数异常高。

**原因**：判据过松，或溶剂水被选入（水–水氢键极多）。

**对策**：

* `--selection "protein"` 排除水；
* 恢复默认判据（2.5 / 3.5 / 120°）；
* 检查是否有多个重叠 altLoc 结构被当作不同原子。

### 2.8 PCA 第一主成分解释 >95% 方差

**症状**：PC1 几乎解释全部运动。

**原因**：整体平移/旋转未被去除，或体系整体在漂移。

**对策**：入口脚本已自动叠合；若仍出现，检查参考帧选择是否合理
（选一个"中间"帧而不是第 0 帧）。

### 2.9 聚类只得到 1 个簇

**症状**：`n_clusters = 1`，但 RMSD 范围很宽。

**原因**：`--cluster-threshold` 太大。

**对策**：不给阈值（让算法按簇数走），或显式 `--n-clusters 5`；
也可查看 RMSD 矩阵的分布确定合理阈值。

### 2.10 聚类得到过多簇（>20）

**症状**：几乎每帧一个簇。

**原因**：阈值太小，或轨迹未收敛（帧间噪声大于构象差异）。

**对策**：增大 `--cluster-threshold`；先确认 RMSD 已进入平台期。

### 2.11 Ramachandran 离群点 >5%

**症状**：`outliers` 比例很高。

**原因**：

1. 拓扑的主链连接错误（残基顺序错乱会导致 φ/ψ 完全乱）；
2. 力场问题；
3. 采样不足。

**对策**：

* 用 PyMOL 检查 backbone 是否连续（`show cartoon`）；
* 确认 `--topology` 与轨迹来自同一次模拟；
* 若拓扑缺失 N 原子，本技能的二面角计算会跳过该残基——检查
  `md_analysis.json` 中 `dihedral.phi_per_frame` 的非零条目数。

### 2.12 SASA 接近 0

**症状**：`mean_total` 几乎为 0 Å²。

**原因**：原子密度极高（坐标未归一化），或范德华半径查表失败（元素名异常）。

**对策**：检查坐标量级（正常蛋白在几十 Å 范围内）；
确认 PDB 的元素列（第 77-78 列）正确，空的话会退化成按原子名推断。

### 2.13 链间接触为空

**症状**：`interchain.pairs` 为空。

**原因**：拓扑只有一条链；或链 ID 判定错误。

**对策**：确认 PDB 的第 22 列（链 ID）非空——用编辑器逐字符检查很容易
看漏，`ATOM` 行的链 ID 位于第 22 列。

## 3. 模拟相关问题

### 3.1 `运行分子动力学模拟需要 OpenMM`

**症状**：报错退出。

**对策**：

```bash
conda install -c conda-forge openmm
# 或
pip install openmm
```

仅做轨迹分析不需要 OpenMM。

### 3.2 `ValueError: No force field for residue` / 参数无法赋原子

**症状**：OpenMM 在 `createSystem` 时报错。

**原因**：结构里含非标准残基（水、离子、辅因子、共价修饰），力场没有对应模板。

**对策**：

* 用 `pdb2pqr` 或 PDBFixer 标准化结构（补末端、加氢、去掉不需要的配体）；
* 或者只保留蛋白部分做模拟；
* 为特定配体准备相应的 XML 参数文件。

### 3.3 模拟几步后就发散

**症状**：能量爆炸，温度飙升。

**原因**：初始结构能量高（未最小化）、存在原子重叠、或时间步长过大。

**对策**：

* 增加 `--minimize-steps`；
* 用 `pdb2pqr` 或 OpenMM 的 `PDBFile` + `Modeller` 重新加氢；
* 减小 `--timestep` 到 0.001 ps。

### 3.4 模拟极慢

**症状**：`run_simulation.py` 耗时超预期。

**对策**：

* `--platform CUDA`（有 GPU 时）；
* 减小体系（截取残基范围）；
* `--solvent vacuum`（仅用于快速测试，不建议做正式分析）。

### 3.5 位置限制看起来没生效

**症状**：开了 `--restraints backbone`，但 NVT 阶段的骨架 RMSD 和
不加限制时几乎一样。

**原因**：位置限制加到了 `System` 上，但 `Simulation`（以及 `Context`）
**已经创建**。OpenMM 的 Context 在创建时就把力场编译进计算核，
之后改 `System` 不会重建 Context —— 力**静默失效**，不报错也不生效。

**对策**：位置限制必须在 `app.Simulation(...)` 之前加到 `System` 上。
本技能的实现已经按这个顺序执行；如果你自己写了模拟脚本，
请检查这个顺序。同样的道理适用于加气压计和撤除限制，
两者都需要**重建 Simulation** 并显式搬运 positions + velocities。

验证方法：对比开/关限制时 NVT 轨迹的骨架 RMSD，应有明显差异。

### 3.6 链在模拟中分离了

**症状**：`chain_separation.final_check.separated` 为 `true`，
或模拟途中出现链间距告警。

**对策**：

1. **先检查输入结构** —— 大多数情况下解离的根源是起始姿态不物理
   （对接姿态没松弛、界面有原子穿插）。先用 `--minimize-steps`
   加大最小化，再跑 NVT。
2. **加位置限制** —— `--nvt-steps 50000 --restraints backbone`，
   让侧链和溶剂先松弛，把界面稳住。
3. **降低初始温度冲击** —— 先用较短、温度较低的 NVT 段
   （例如 100 K → 300 K 分两步），再进生产期。
4. **检查是不是误报** —— 阈值 30 Å 对超大体系可能偏小，
   可用 `--chain-separation-threshold` 调整。若两个体系本来就离得远
   （例如有共价交联），绝对距离会一直很大，这时应关注
   `distance_history` 里的**相对变化**而非绝对值。

### 3.7 重建 Simulation 后轨迹跳变

**症状**：NPT 阶段或撤除限制后，轨迹出现一个明显的坐标跳变。

**原因**：重建 `Simulation` 时没有搬运当前的位置和速度。新 Context
的初始位置是 0，跑到第一次 `step` 之前都停在原点。

**对策**：重建后必须显式搬运状态：

```python
state = simulation.context.getState(getPositions=True, getVelocities=True)
simulation = app.Simulation(topology, system, integrator, platform)
simulation.context.setPositions(state.getPositions())
simulation.context.setVelocities(state.getVelocities())
```

**注意速度也要搬** —— 只搬位置的话，体系速度为零，等于在生产期
开头做了一次隐式的"零动能重启"，前几百步的行为完全失真。

## 4. 自检清单

```bash
# 1) 四个入口的帮助
python skills/molecular-dynamics/scripts/run_analysis.py --help
python skills/molecular-dynamics/scripts/run_simulation.py --help
python skills/molecular-dynamics/scripts/trajectory_io.py --help
python skills/molecular-dynamics/scripts/md_equilibration.py --help

# 2) 轨迹解析器自检
python skills/molecular-dynamics/scripts/trajectory_io.py --input work/md_test_trajectory.dcd

# 3) 核心算法自检（合成轨迹）
python skills/molecular-dynamics/scripts/md_analysis.py

# 4) 全流程
python skills/molecular-dynamics/scripts/run_analysis.py \
    --input work/md_test_trajectory.dcd \
    --topology work/md_test_topology.pdb \
    --output work/out_md

# 5) 平衡协议自检（仅检查 OpenMM 是否可用）
python skills/molecular-dynamics/scripts/run_simulation.py \
    --input work/md_test_topology.pdb --check-only

# 6) 规范校验
python scripts/validate_skills.py
```
