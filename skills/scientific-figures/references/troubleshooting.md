# 科学图表常见问题与排错

---

## 1. matplotlib 缺失

### 症状

```
FigureStyleError: matplotlib 未安装，无法绘图。请先安装：
    pip install matplotlib
```

这是**有意的降级设计**——脚本在导入时捕获 ImportError，保存到模块级变量
`MATPLOTLIB_AVAILABLE`，在真正需要绘图时才抛出带安装指引的异常，
而不是让用户看到一串 `ImportError` 堆栈。

### 解决

```bash
pip install matplotlib
# 或
pip install -r skills/scientific-figures/requirements.txt
```

---

## 2. 列名不匹配

### 症状

```
SkillError: ROC 输入需要 Model / FPR / TPR 列
```

### 原因与解决

脚本的列名匹配是**大小写不敏感**的，并且对常见变体做了候选列表：

| 图类型 | 首选列名 | 可接受的变体 |
| --- | --- | --- |
| roc | `Model` | `model`, `method`, `classifier` |
| roc | `FPR` | `fpr`, `false_positive_rate` |
| roc | `TPR` | `tpr`, `true_positive_rate`, `sensitivity` |
| roc | `AUC` | `auc`（可选，缺失时自行计算） |
| confusion | `TN` `FP` `FN` `TP` | `true_negative` 等 |
| plddt | `Global pLDDT` | `plddt`, `pLDDT`, `confidence`, `score` |
| shap | 任意数值列 | 自动识别（排除 id 类列） |
| class-distribution | `class` | `Class`, `纲`；`count` / `数量` / `n` |
| heatmap | 首列为行标签，其余为数值 | 用 `--row-labels none` 关闭行标签 |

如果你的列名不在候选列表里，改列名是最快的办法：

```bash
python - <<'PY'
import csv
with open("data/roc.csv", newline="", encoding="utf-8-sig") as f:
    rows = list(csv.DictReader(f))
for r in rows:
    r["Model"] = r.pop("algo")
    r["FPR"] = r.pop("fp")
    r["TPR"] = r.pop("tp")
with open("data/roc_std.csv", "w", newline="", encoding="utf-8-sig") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)
PY
```

---

## 3. 热图相关

### 3.1 "热图输入没有数值"

**原因**：首列被当成行标签（默认 `--row-labels auto`），但数据实际是纯矩阵
（没有标签列）。

**解决**：`--row-labels none`。此时首列也会被当作数值列读取。

### 3.2 颜色几乎全部一样

**原因 1**：色阶被极端离群值拉宽。
**解决**：`--percentile-clip`（2%/98% 分位裁剪）。

**原因 2**：数据的动态范围本来就小（例如全是 -8 到 -7 kcal/mol）。
**解决**：这属于真实情况，不要人为拉伸。在图注中写明色阶范围即可。

### 3.3 X 轴标签重叠

**原因**：列数太多。脚本按 `n//20` 抽样，但列名本身很长时仍会挤。

**解决**：
- 缩短列名（用受体的简称而非全名）
- 增大 `--width`
- 或改用不带标签的矩阵图（`--row-labels none` + 自己在图注中给出顺序）

### 3.4 矩阵不是方阵

脚本支持任意形状的矩阵（行数与列数可不同），但热图用于方阵时最易读。
非方阵（受体 × 配体的结合能矩阵）也很常见，属正常用法。

---

## 4. ROC 相关

### 4.1 AUC 值不对

**原因 1**：输入的 (0,0) 与 (1,1) 端点缺失，梯形积分面积偏小。
**解决**：确认 CSV 的每条曲线都包含 FPR=0 与 FPR=1 两行。

**原因 2**：CSV 里有 `AUC` 列且值与曲线不自洽——脚本优先用 CSV 的值。
**解决**：删除 `AUC` 列让脚本自己算，或检查上游生成 CSV 的脚本。

### 4.2 Z 检验的 p 值不可信

**原因**：`--n-samples` 默认 500 是假设值。Hanley-McNeil 公式中的 n 必须是
**真实验证集样本数**。

**解决**：
```bash
--n-samples 1200    # 换成你的实际验证集大小
```

样本量小于 30 时这个近似本身就不可靠，不建议报告 p 值。

### 4.3 `--bootstrap-ci` 很慢

每次 bootstrap 都要重算一次 AUC（对 141 个点排序 + 积分）。
默认 200 次，对大曲线尚可；若曲线有上万点会明显变慢。
减小点数或去掉该选项。

---

## 5. 混淆矩阵相关

### 5.1 只有部分子图有颜色

**原因**：输入 CSV 中同一模型出现了多行（脚本为每行建一个子图）。

**解决**：确认 CSV 中每个模型只有一行。若上游输出有重复，先去重：

```bash
python - <<'PY'
import csv
seen, rows = set(), []
with open("data/cm.csv", newline="", encoding="utf-8-sig") as f:
    for r in csv.DictReader(f):
        if r["Model"] in seen: continue
        seen.add(r["Model"]); rows.append(r)
with open("data/cm_clean.csv", "w", newline="", encoding="utf-8-sig") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)
PY
```

### 5.2 颜色深浅在不同模型间不可比

脚本用全局 vmin/vmax，正常情况下不会出现这个问题。若出现，检查是否有
子图的值超出全局范围（正常会被截断到边界色，不会造成"深浅不一"）。
真正的"深浅不一"通常意味着每个子图各自归一化了——那不是本脚本的行为。

---

## 6. SHAP 相关

### 6.1 "SHAP 输入需要至少一个数值特征列"

**原因**：CSV 里的列不是数字（如已经是类别编码的字符串，或 ID 列是数字被误判）。

**解决**：检查列内容。ID 列如果全是数字（如 1,2,3...），脚本无法自动识别，
会把它当成一个特征参与计算。

**缓解**：把 ID 列改名成含 `id`/`sample`/`index` 的名字，脚本会优先排除它。

### 6.2 蜂群图是"均匀抖动"而不是真正的蜂群

脚本用随机抖动模拟蜂群（与原仓库实现一致），点密度不会像 shap 库的
核密度估计那样"挤在一起"。

若需要与 shap 库对比，改用 shap 库；若只是做汇报图，本实现足够。

### 6.3 特征名太长导致 Y 轴标签被裁切

**原因**：`tight bbox` 会裁切，但超长的 y 标签可能超出画布。

**解决**：缩短特征名，或增大画布高度（脚本按特征数自动调整：
`figsize=(10, max(6, 0.45*n+2))`）。

---

## 7. 字体问题

### 7.1 找不到 Arial 的 warning

```
findfont: Font family 'Arial' not found. Falling back to DejaVu Sans.
```

**影响**：图能正常生成，但字形与投稿要求有细微差异。

**解决**：
- Windows：一般自带 Arial，无需处理
- Linux：安装 msttcorefonts 或 `fonts-liberation`（Liberation Sans 是 Arial 的
   metrics 兼容替代）
- macOS：Arial 随 Office 安装；没有则用 Helvetica（已在回退列表中）

### 7.2 中文标签显示为方块

**原因**：所选字体不含中文字形。

**解决**：在 `figure_style.py` 的 `font.sans-serif` 中加入中文字体，
如 `['Arial', 'Microsoft YaHei', 'SimHei', 'Noto Sans CJK SC', 'DejaVu Sans']`。
建议图表标签统一用英文，避免投稿排版问题。

---

## 8. 输出问题

### 8.1 只生成了一个文件

`--formats png pdf` 应该生成两个。若只生成一个，说明 `--output` 里已经
带了扩展名（脚本此时以该扩展名为准，只出一个文件）。

### 8.2 文件名带双扩展名

`--output figures/fig8_roc` + `--formats png pdf` → `fig8_roc.png`、`fig8_roc.pdf`。
若你写成 `--output figures/fig8_roc.png`，则只生成 `fig8_roc.png`。

### 8.3 保存失败

```
FigureStyleError: 保存图片失败 ...
```

**原因**：目录不存在或无写权限。脚本会自动创建 `--output` 的父目录，
若仍失败则是权限问题。

---

## 9. 自检命令

```bash
cd <repo>
export PYTHONPATH="$PWD/packages"

# frontmatter 校验
python -c "import yaml; d=yaml.safe_load(open('skills/scientific-figures/SKILL.md').read().split('---')[1]); print(d['name'], '|', d['description'][:80])"

# CLI 自检
python skills/scientific-figures/scripts/make_figures.py --help

# 端到端：造样例数据并生成 ROC 与混淆矩阵
python - <<'PY'
import csv, math
rows = []
for model, shift in (("ModelA", 0.9), ("ModelB", 0.6)):
    auc = 0
    for i in range(21):
        fpr = i / 20
        tpr = min(1.0, shift + (1 - shift) * (1 - (1 - fpr) ** 2))
        rows.append({"Model": model, "FPR": fpr, "TPR": tpr, "AUC": 0})
        auc += (1 / 20) * tpr if i else 0
    for r in rows:
        if r["Model"] == model:
            r["AUC"] = round(auc, 4)
with open("/tmp/roc.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=["Model", "FPR", "TPR", "AUC"]); w.writeheader(); w.writerows(rows)
with open("/tmp/cm.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=["Model", "TN", "FP", "FN", "TP"])
    w.writeheader()
    w.writerow({"Model": "ModelA", "TN": 85, "FP": 15, "FN": 20, "TP": 80})
    w.writerow({"Model": "ModelB", "TN": 70, "FP": 30, "FN": 35, "TP": 65})
PY
python skills/scientific-figures/scripts/make_figures.py roc --input /tmp/roc.csv --output /tmp/roc
python skills/scientific-figures/scripts/make_figures.py confusion --input /tmp/cm.csv --output /tmp/cm
ls -la /tmp/roc.png /tmp/cm.png
```

---

## 10. 交付前自查清单

- [ ] 所有图的字体一致（同一 `--style` 生成）
- [ ] 图例不遮挡数据（ROC 的图例在右下，混淆矩阵无图例）
- [ ] 坐标轴有明确的单位与含义
- [ ] 色阶范围已在图注说明（尤其用了 `--percentile-clip` 时）
- [ ] 统计方法的 n 值与真实样本量一致
- [ ] PNG 已导出用于预览，PDF 已导出用于投稿
- [ ] 灰度转换后各系列仍可区分（用 Okabe-Ito 调色板的原因）
- [ ] 图注写清了数据来源、样本量、统计方法与缩写
