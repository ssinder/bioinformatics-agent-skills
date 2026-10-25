# 多肽分析常见问题与排错

---

## 1. 输入与解析问题

### 1.1 "输入文件没有可分析的记录"

**原因**：CSV 中没有找到序列列。

**解决**：脚本按 `--sequence-column`（默认 `sequence`）找列，找不到时会依次
回退到 `seq` / `Sequence` / `peptide` / `aa_seq`。若你的列名不同，显式指定：

```bash
--sequence-column aa_sequence
```

也可用 `-v` 查看 DEBUG 日志确认读了多少行。

### 1.2 FASTA 被误判成 CSV

**原因**：文件扩展名不是 `.fasta`，且首行不是 `>`。

**解决**：显式指定 `--format fasta`。

### 1.3 一条序列被拆成多行 / 被合并

本脚本用共享层 `read_fasta` 读取，它支持标准多行 FASTA（续行直接拼接）、
`;` 开头的注释行和 CRLF 换行。如果你的文件用**其它格式**记录续行
（如某些数据库导出的 `&` 分隔），需要先转成标准 FASTA。

### 1.4 "跳过 N 条无法处理的记录"

**原因**：序列为空，或 `--drop-unknown` 打开后仍含非标准残基（X/B/Z/U）。

**解决**：序列为空通常是 CSV 解析时列名不匹配；含 X 的序列如果确实要保留，
去掉 `--drop-unknown`（此时这些序列仍会被计算，只是 `qc_flags` 会标记）。

---

## 2. 数值异常

### 2.1 分子量偏大或偏小

| 现象 | 原因 | 处理 |
| --- | --- | --- |
| 偏大约 18 Da | 序列里含 `*` 终止符或 `-` gap 字符 | 脚本已自动去 `*`，但 `-` 不会被去除，需手工清理 |
| 偏小 | 序列被截断（FASTA 续行缺失） | 检查原文件 |
| 差一个整数量级 | 序列列混入注释文字 | 确认 `--sequence-column` 指向正确列 |

**校验方法**：20 aa 的典型肽 MW 应在 2000–2500 Da。10 aa 肽约 1100–1200 Da。

### 2.2 等电点超出 0–14 或明显异常

**原因**：序列极短（< 3 aa）时电荷曲线可能不穿过零点；或含大量 His（pKa 对
环境极敏感）。

**解决**：脚本在曲线不穿过零时返回"电荷绝对值最小的点"，这在极端情况下
不是真正的 pI。短肽（< 5 aa）请结合具体组成判断，不要依赖该数值。

### 2.3 不稳定指数与预期差很多

**原因**：Guruprasad 公式来自蛋白数据，对 20 aa 短肽指示性弱。含 Pro/Gly
（DIHEX 权重最低）的序列会显著偏低。

**处理**：`qc_flags` 只在 > 40 时提示"预测不稳定"，不会自动剔除。短肽的
稳定性判断请结合实验数据或专用的短肽稳定性预测工具。

### 2.4 Boman 值为负

**正常**。负值表示序列的带电残基多于疏水残基，多肽本身不是抗菌肽——
这正是构建负样本集时需要的结果。

---

## 3. 判别结果不符合预期

### 3.1 所有序列都被判为 `uncertain`

**原因 1**：阈值太高。默认 `--antimicrobial-threshold 5.0`，若注释层得分低
（只有一两个关键词），综合分很容易落在 2–5 之间。

**原因 2**：纯 FASTA 输入没有 GO/keyword，注释层为 0，综合分上限只有 4.0。

**解决**：
```bash
# 纯序列输入
--antimicrobial-threshold 2.5 --non-antimicrobial-threshold 1.0
# 有注释但关键词少
--antimicrobial-threshold 4.0
```

### 3.2 所有序列都被判为 `antimicrobial`

**原因**：阈值太低，或注释词表太宽（如 "defense" 会命中大量泛防御蛋白）。

**解决**：提高阈值到 5–6；检查 `hit_keywords` 列，看是哪条规则在贡献分数。
若命中的是 "Immunity"/"Defense" 这类泛化词，说明该条目不是直接的抗菌肽，
应提高 `--antimicrobial-threshold` 或用 `--exclude-terms` 做二次排除。

### 3.3 判别结果与实验结论相反

**可能原因**：
1. **注释不可靠**：UniProt 的 keyword 是人工+自动标注，可能存在过度注释
2. **理化外推的局限**：Boman 值与 KD 无法反映空间排列
3. **物种差异**：抗真菌肽与抗细菌肽的理化要求不同（真菌膜含麦角甾醇，
   需要更强的疏水性；细菌膜带负电，需要正电荷）

**处理**：把本技能当作**预筛工具**而非最终判定。所有关键结论必须用实验
（MIC/MFC）验证。可用 `--label all` 导出全部并附分数，人工划定边界。

### 3.4 `--label` 过滤后结果为空

检查 `label` 列的实际取值。`uncertain` 是独立的一个类，
若想保留所有非明确判定的样本，用 `--label uncertain` 或 `--label all`。

---

## 4. 聚类问题

### 4.1 "未找到 MMseqs2 可执行文件"

**原因**：`mmseqs` 不在 PATH 中。

**解决**：
```bash
# 确认安装
which mmseqs
# 或显式指定用内置算法
--method internal
```

若用 `conda install -c bioconda mmseqs2` 安装，PATH 里通常是 `mmseqs`。

### 4.2 MMseqs2 聚类失败并自动回退

脚本会记录 WARNING 并回退到内置算法。若你**必须**用 MMseqs2（不允许回退），
用 `--method mmseqs`，此时失败会直接报错。

常见失败原因：
- 磁盘空间不足（MMseqs2 会建临时索引）
- FASTA 中有非法字符（`--`、数字开头等）
- 线程数超出机器核心数（`--threads` 调小）

### 4.3 内置算法太慢

内置算法最坏情况是 O(n²) 的 LCS 计算。1000 条 × 20 aa 大约几秒；
5000 条 × 60 aa 会到分钟级。

**解决**：
1. 提高 `--min-identity` —— k-mer 预筛会跳过更多对
2. 增大 `--similarity-weight`（如 0.8）—— 预筛下界提高，更多对被跳过
3. 减少 `--kmer-k` 增大 k —— k-mer 更稀疏，Jaccard 更严格
4. 直接安装 MMseqs2

### 4.4 所有序列都聚成单个大簇

**原因**：`--min-identity` 太低（如 0.1），或 k-mer 过短导致 Jaccard 普遍偏高。

**解决**：提高 `--min-identity` 到 0.3–0.5，并检查簇大小分布报告。

### 4.5 所有序列都是单序列簇（零冗余）

**原因**：`--min-identity` 太严（如 0.8 用在短肽上）。

**解决**：短肽库用 0.3–0.4。报告里的"单序列簇数"可以直接告诉你这个情况。

### 4.6 代表序列看起来不像簇内成员

**原因**：默认 `medoid` 用的是**综合相似度**（identity 与 Jaccard 的加权），
当 identity 权重高时 medoid 可能选了个"长度居中但序列不像"的序列。

**解决**：用 `--representative longest` 或 `first` 交叉对比。

---

## 5. 环境与依赖

### 5.1 依赖是可选的

`analyze_peptides.py` 与 `cluster_peptides.py` **完全不需要**第三方库——
理化性质计算是纯 Python 实现（仅用 `math` 与 `statistics`）。

只有 `cluster_peptides.py --method mmseqs` 需要外部 MMseqs2 二进制。

这意味着在受限环境中也能跑通核心流程。

### 5.2 Windows 控制台

脚本日志走 stderr，Windows 终端下非 ASCII 内容可能显示异常，但不影响
文件写出。CSV 用 UTF-8 BOM，Excel 双击打开不乱码。

### 5.3 Python 版本

需要 **Python 3.9+**（`from __future__ import annotations` + 内置泛型标注）。
3.9 以下不支持 `list[str]` 这类写法。

---

## 6. 自检命令

```bash
cd <repo>
export PYTHONPATH="$PWD/packages"

# frontmatter 校验
python -c "import yaml; d=yaml.safe_load(open('skills/peptide-design/SKILL.md').read().split('---')[1]); print(d['name'], '|', d['description'][:80])"

# CLI 自检
python skills/peptide-design/scripts/analyze_peptides.py --help
python skills/peptide-design/scripts/cluster_peptides.py --help

# 端到端：造一个小样例
python - <<'PY'
from pathlib import Path
Path("/tmp/demo.fasta").write_text(
    ">pep1 antimicrobial\nKKRLKAFLKRL\n"
    ">pep2 defensin\nCRLVRGTTVRR\n"
    ">pep3 repeat of pep1\nKKRLKAFLKRL\n"
)
PY
python skills/peptide-design/scripts/analyze_peptides.py \
    --input /tmp/demo.fasta --output /tmp/demo_out.csv --report /tmp/demo_report.txt
python skills/peptide-design/scripts/cluster_peptides.py \
    --input /tmp/demo.fasta --method internal --min-identity 0.8 \
    --output /tmp/demo_rep.fasta
```

---

## 7. 结果可信度自查清单

- [ ] 分子量与 pI 抽检 5 条，与 BioPython `ProtParam` 或 Expasy 对照
- [ ] 检查 `qc_flags` 列，有未知残基的记录单独处理
- [ ] 判别结果抽查各 label 的 `hit_keywords` 与 `classification_reasons`
- [ ] 负样本集必须确认 `protein_existence` 包含 "Evidence at protein level"
- [ ] 聚类的簇大小分布合理（不是全 1 也不是全 N）
- [ ] 关键肽的实验活性验证完成后再下结论
