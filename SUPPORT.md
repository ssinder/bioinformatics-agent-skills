# 支持与求助

先查这里，再开 issue——大部分问题能就地解决。

## 我该怎么用某个技能？

每个技能的 `SKILL.md` 都有完整的用法说明：

- **When to Use This Skill** —— 什么时候该用它、什么时候不该用
- **Quick Start** —— 可直接复制运行的命令
- **Parameters** —— 全部参数、默认值与含义
- **Outputs** —— 产出文件及其含义
- **Validation** —— 怎么检查结果是否正确（**最重要的一节**）

`references/methodology.md` 讲参数背后的科学依据，
`references/troubleshooting.md` 讲常见报错与排错。

## 报告问题

请开 [issue](https://github.com/ssinder/bioinformatics-agent-skills/issues)，
并使用对应的模板。报告 bug 时请附上：

1. 技能名与执行的完整命令
2. 报错信息或异常结果（完整输出，不要截断关键行）
3. 运行环境：`python --version`、`pip list | grep -E "numpy|pandas|scipy|scanpy"`
4. 输入数据的最小复现样本（如果可以共享）

## 提问与讨论

使用方法类问题、不确定该用哪个技能、想讨论方法学取舍，请到
[Discussions](https://github.com/ssinder/bioinformatics-agent-skills/discussions)。
这里比 issue 更适合开放式讨论。

## 功能建议

想增加新技能或改进现有技能，请开
[feature request](https://github.com/ssinder/bioinformatics-agent-skills/issues/new?template=feature_request.md)。
请说明你实际想解决的科研场景——这比说明技术实现更有助于判断优先级。

## 报告安全问题

**不要**为此开公开 issue。请按 [SECURITY.md](SECURITY.md) 的说明私下报告。

## 常见问题

<details>
<summary>安装后 Claude Code 认不出我的技能</summary>

不同宿主扫描的目录不同。请确认安装到了对应位置：

| 宿主 | 技能目录 |
|---|---|
| Claude Code（个人） | `~/.claude/skills/` |
| Claude Code（项目） | `.claude/skills/` |
| OpenCode / 部分宿主 | `.agents/skills/` |

用 `npx skills add ssinder/bioinformatics-agent-skills` 安装可以自动处理路径。
</details>

<details>
<summary>技能安装了但从不自动触发</summary>

两个常见原因：

1. **装了太多技能。** 建议按工作流每次装 2–3 个。装载的技能过多时，
   宿主会为了控制上下文预算而丢弃部分技能的描述，那些技能便无法被匹配。
2. **描述与你的用词不匹配。** 描述里应包含你实际会用的词（如
   `RNA-seq`、`h5ad`、`BAM`、`DESeq2`），而不是抽象的科学名词。
</details>

<details>
<summary>缺少 numpy / pandas 时脚本报错</summary>

多数技能有纯 Python 降级路径，缺少科学计算栈也能跑完主要流程。
若某个技能确实需要特定依赖，错误信息会给出确切的安装命令。
</details>

<details>
<summary>结果和 R / 原版工具对不上</summary>

先看该技能的 **Validation** 章节，它写明了正常的结果范围与常见偏差来源。
`bulk-rnaseq` 与 `single-cell-rnaseq` 提供双后端，可直接对照官方实现。
</details>
