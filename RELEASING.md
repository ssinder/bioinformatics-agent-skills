# Releasing

本仓库的发布流程，以及每一步**需要你在平台侧做的操作**。

自动化部分已经就绪；少数几步必须在对应平台的网页后台完成。

---

## 一次性准备（首次发布前）

### ✅ 1. GitHub Environment `pypi`（已创建）

已在仓库 Settings → Environments 下创建名为 **`pypi`** 的 environment。

> **名字必须与 workflow 完全一致。** `.github/workflows/release.yml` 里
> `publish-pypi` job 声明了 `environment: pypi`；Trusted Publishing 的
> 三元组里填的 Environment 必须一字不差，否则 OIDC 令牌签发会失败。

### 🔲 2. PyPI 账号与项目（需要你操作）

本仓库**还没有 PyPI 账号**，也没有 `biopipe-agent-skills` 这个项目。

> `biopipe` 这个包名已被第三方占用（dohlee/python-biopipe），所以
> 用了 `biopipe-agent-skills`。这个名字在 PyPI 上目前空闲。

请到 <https://pypi.org/register/> 注册账号（若已有则直接登录），然后：

1. 创建项目 `biopipe-agent-skills`：

   <https://pypi.org/manage/project/biopipe-agent-skills/>

   项目名 **必须**是 `biopipe-agent-skills`——一旦创建就不能改名。

2. 在项目设置 → Publishing → Add a GitHub publisher，填入：

   | 字段 | 值 |
   |---|---|
   | Repository | `ssinder/bioinformatics-agent-skills` |
   | Workflow | `release.yml` |
   | Environment | `pypi` |

   > Environment 这一项填 `pypi`——就是上面已创建好的那个。

### 🔲 3. TestPyPI 演练（强烈建议）

Trusted Publishing 的首次配置**几乎总是**在某处有个笔误（workflow 名、
environment 名、或仓库名），而 **PyPI 不允许重传同名文件**——一旦发错版本，
只能 yank（撤回警告）而不能彻底删除。

演练步骤：

1. 在 <https://test.pypi.org/manage/project/biopipe-agent-skills/> 重复第 2 步
   （TestPyPI 的包名通常需要带版本后缀，如 `biopipe-agent-skills` 即可）
2. 临时把 `release.yml` 里 `publish-pypi` job 的 `environment.name`
   改成 `pypi-test`，在 TestPyPI 侧登记同样的三元组
3. 打一个测试 tag 触发，确认成功后再改回来

### 🔲 4. Zenodo（概念 DOI）

1. 在 <https://zenodo.org/account/settings/applications/> 创建 GitHub 应用，
   授权仓库 `ssinder/bioinformatics-agent-skills`。
2. 在仓库 Settings → General → Repository metadata → 勾选
   **Enable Zenodo GitHub Integration**。

之后每次打 tag 触发 Release，Zenodo 会自动收录并分配 DOI：

- **version DOI** —— 指向该具体版本
- **concept DOI** —— 指向最新版本，引用时应当使用它

元数据来自 `CITATION.cff` 与 `.zenodo.yaml`。

### ✅ 5. skills.sh 徽章（已就绪）

已被自动收录（收录靠安装遥测，无需投稿）。徽章已接好，安装量随安装次数更新。

---

## 发布流程

### 步骤 1 · 确认工作区干净且检查全绿

```bash
make check          # 规范校验 + 静态检查 + 测试
git status          # 应为空
```

### 步骤 2 · 同步版本号

版本号要同时改三个地方：

```bash
python scripts/ci_version_sync.py --set 1.0.1
```

它会统一改写 `pyproject.toml`、`plugin.json`、`CITATION.cff`，
并把 `CITATION.cff` 的 `date-released` 更新为今天。

> 技能自身的 `metadata.version` 是**独立**的——某个技能改了内容，
> 可以只升它的版本，不必升整个包。所以 `--set` 不会碰技能。

### 步骤 3 · 更新 CHANGELOG

在 `CHANGELOG.md` 顶部加新版本段落，写清**用户可见的变化**。

### 步骤 4 · 提交并打 tag

```bash
git add -A
git commit -m "release: v1.0.1"
git tag -a v1.0.1 -m "v1.0.1"
git push origin main --follow-tags
```

打 tag 会同时触发：

| 工作流 | 作用 |
|---|---|
| `release` | 创建 GitHub Release（附 wheel + sdist） |
| `publish-pypi` | 上传 PyPI + 生成 attestations |
| Zenodo 集成 | 归档并分配 DOI |

### 步骤 5 · 验证

```bash
# GitHub Release
gh release view v1.0.1      # 或直接打开仓库的 Releases 页

# PyPI
pip install biopipe-agent-skills==1.0.1
bioagent --list              # 应列出 17 个技能

# 安装后真的能跑
bioagent scientific-figures roc --input metrics.csv --output out/roc
```

### 步骤 6 · 同步各处引用

- `README.md` 的 PyPI 徽章会自动更新（shields.io）
- skills.sh 的安装量在下次安装后更新
- 若引用了 concept DOI，无需改动（它自动指向最新版本）

---

## 回滚

PyPI 不允许删除已发布的版本，只能 yank：

```bash
# 在 https://pypi.org/manage/project/biopipe-agent-skills/releases/ 里操作
```

yank 后 `pip install biopipe-agent-skills` 仍会选到该版本，
但会给出警告。完全撤回需要删除版本——**不可逆**。

GitHub Release 可以直接删除：

```bash
gh release delete v1.0.1
```

---

## 版本号策略

遵循语义化版本：

| 变更 | 版本位 |
|---|---|
| 修 bug、改进文档、改依赖 | patch（1.0.0 → 1.0.1） |
| 新增技能、新增参数、改变默认行为 | minor（1.0.0 → 1.1.0） |
| 移除技能、移除参数、改变既有输出的格式 | major（1.0.0 → 2.0.0） |

**新增技能是 minor，不是 patch**——用户会按版本号判断兼容性，
而新增技能会改变 `bioagent --list` 的输出。

---

## 当前状态

| 项目 | 状态 |
|---|---|
| GitHub Environment `pypi` | ✅ 已创建 |
| `release.yml` 的 Trusted Publishing 配置 | ✅ 已就绪（`id-token: write` + `environment: pypi`） |
| `pypa/gh-action-pypi-publish@release/v1` | ✅ 已接入（自动生成 PEP 740 attestations） |
| `.zenodo.yaml` | ✅ 已就绪 |
| skills.sh 收录 | ✅ 已收录 |
| PyPI 账号 | 🔲 需要你注册 |
| PyPI 项目 `biopipe-agent-skills` | 🔲 需要你创建 |
| PyPI GitHub publisher 三元组 | 🔲 需要你填写 |
| Zenodo GitHub 集成 | 🔲 需要你授权 |

**剩下的三步都必须在 PyPI / Zenodo 的网页后台完成**——PyPI 不支持
第三方登录，且 `biopipe-agent-skills` 这个项目名一旦占用就不能更改。
