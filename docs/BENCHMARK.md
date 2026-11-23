# Benchmark

本文件记录本仓库的**实测**结果。所有数字都来自仓库内的自动化测试，
可通过 `make test` 或 `pytest tests` 复现。

本文件采用 [NVIDIA Verified Skills](https://docs.nvidia.com/skills) 的
`BENCHMARK.md` 约定：给出结论、测量到的变化，以及每项测量的方法。

---

## 1. 后端一致性（Backend agreement）

### 结论

内置的纯 Python / NumPy 统计实现，与官方参考实现
（Bioconductor DESeq2、scanpy）在**统计上是等价的**——差异在数值噪声
范围内，不存在系统性偏差。

这是本仓库最重要的正确性保证。43,507 行自研实现若与参考实现有偏离，
所有下游结论都不可信；因此这项对照**在 CI 里强制执行**，不是可选的。

### 方法

- **数据**：仓库内的固定夹具 `tests/fixtures_bulk/`（3 个处理组 + 3 个对照
  组重复，含 counts、tx2gene、样本元数据、GO 注释）
- **对照**：同一份数据分别用 `--backend python` 与 `--backend r` 跑，
  比较输出；单细胞同理比较 `--backend numpy` 与 `--backend scanpy`
- **指标**：log2FC 排序的 Spearman 秩相关；baseMean 的秩相关；
  显著基因集合的 Jaccard 相似度；显著基因调控方向是否一致
- **验证器**：确定性断言，不使用 LLM 判定

### 阈值与实测

| 测量项 | 阈值 | 依据 |
|---|---|---|
| log2FC 排序 Spearman ρ | > 0.95 | 低于此值说明存在系统性偏差 |
| baseMean Spearman ρ | > 0.99 | 同量纲的量，应当几乎完全一致 |
| 显著基因集合 Jaccard | > 0.80 | 显著基因结论不一致等于结论不可用 |
| 显著基因调控方向 | 0 个相反 | 方向相反是最严重的错误 |
| 单细胞簇数差异 | ≤ 1 个簇 | SKILL.md 的明确承诺 |
| 单细胞聚类结构 ARI | > 0.80 | 划分结构应当高度一致 |
| 单细胞对真值聚类 ARI | > 0.60 | 夹具内含 600 个细胞的人工真值标签 |

### 复现

```bash
# 需要 Rscript + Bioconductor DESeq2
pytest tests/test_backend_agreement.py -v

# 需要 scanpy + anndata
pytest tests/test_scrna_backend_agreement.py -v
```

R 或 scanpy 缺失时相关用例会**跳过**而非失败；但 CI 的
`backend-agreement` 与 `scrna-backend-agreement` 两个 job 在依赖缺失时
会**直接失败**——避免出现"永远绿灯却什么都没验证"的假徽章。

---

## 2. 评测用例覆盖（Evaluation coverage）

### 结论

17 个技能全部配备评测用例，共 **87 个**用例，每个都含**至少一个近似负例**
（`should_not_trigger`）。

### 为什么近似负例是必需的

SkillsBench（[arXiv:2602.12670](https://arxiv.org/abs/2602.12670)）在 87
个任务上的配对评测发现，**13 个任务在挂上技能后表现反而更差**；而
[Vercel 的一手评测](https://vercel.com/blog/agents-md-outperforms-skills-in-our-agent-evals)
测到 **56% 的场景里技能从未被调用**——因为两个技能描述重叠，代理选了
另一个，或两个都没选。

近似负例正是防住这两种失败的机制：它把"这个请求该归谁"写成了可测的
断言，而不是留给模型临场判断。

### 方法

- 每个技能目录下 `evals/evals.json`，遵循 Agent Skills 规范对评测的约定
  （[evaluating-skills.mdx](https://github.com/agentskills/agentskills)）
- 每个用例包含 `prompt`（真实用户口吻）、`expected_output`（评分依据）、
  `assertions`（可核查的断言）
- 结构与覆盖度由 `scripts/ci_check_evals.py` 在 CI 中校验

### 用例构成

| 类型 | 数量 | 作用 |
|---|---|---|
| 正常请求 | 每个技能 1+ | 验证基本能力 |
| 科学痛点问题 | 每个技能 1+ | 验证领域知识是否真的被用上 |
| 边界与配置 | 每个技能 1+ | 验证参数与阈值判断 |
| 近似负例 | 每个技能 ≥1 | 验证**不**该抢别人的请求 |
| 降级场景 | 多个技能 | 验证缺依赖时的行为 |

### 复现

```bash
python scripts/ci_check_evals.py --require-all
```

---

## 3. 技能规范符合性（Spec conformance）

### 结论

17 个技能全部通过 Agent Skills 规范的校验，**0 错误 0 警告**。

校验覆盖：frontmatter 必填字段、`name` 与目录名一致且为合法 kebab-case、
`description` 的长度与触发语句、metadata 的 string→string 映射、
八个必需章节、脚本的语法与硬编码路径、references 非空、
`requirements.txt` 存在。

### 额外的仓库级门禁

| 检查 | 脚本 | 拦截的问题 |
|---|---|---|
| 技能列表预算 | `validate_skills.py` | name+description 合计超 8,000 字符时，最少被用的技能会被丢弃描述、从此无法触发 |
| 跨宿主兼容 | `ci_harness_check.py` | description 超 1,024 字符时部分宿主**整个技能不加载**；frontmatter 含规范外字段时打包失败 |
| 版本同步 | `ci_version_sync.py` | plugin.json / pyproject.toml / CITATION.cff 三处版本漂移 |
| 仓库一致性 | `ci_consistency_check.py` | README 徽章、技能清单、REGISTRY 与磁盘不一致 |
| 评测覆盖 | `ci_check_evals.py` | 技能缺少评测文件或缺少近似负例 |

### 复现

```bash
python scripts/validate_skills.py --strict
python scripts/ci_harness_check.py
python scripts/ci_version_sync.py
python scripts/ci_consistency_check.py
python scripts/ci_check_evals.py --require-all
```

---

## 4. 无依赖降级（Graceful degradation）

### 结论

在屏蔽 numpy / pandas / scipy / matplotlib / seaborn 的子进程中，
**全部 17 个技能的入口脚本仍能响应 `--help`**。

这直接支撑了仓库"能力不依赖依赖"的设计主张：没有科学计算栈的环境里，
技能依然可用，只是精度或速度受限。

### 复现

```bash
python scripts/ci_smoke_no_deps.py
```

---

## 5. 尚未覆盖的方面

诚实列出本文件**没有**测量的部分：

- **代理端到端任务成功率**。上述评测是"用例设计层面"的覆盖，还需要在真实
  agent 里跑 with-skill / without-skill 配对，才能得到任务通过率的差值
  （Skill Lift）。这需要模型调用成本，已列入路线图。
- **与 BixBench 等外部基准的对齐**。本仓库的夹具是自己构造的，
  不等同于公共基准的难度分布。
- **大样本性能**。当前夹具规模小（100 基因 / 600 细胞），
  足以验证统计一致性，不足以反映大数据的运行时表现。

---

## 6. 引用方式

引用本仓库时，请同时引用 [CITATION.cff](CITATION.cff) 中的软件条目。
具体技能的方法学依据见各技能的 `references/methodology.md`。
