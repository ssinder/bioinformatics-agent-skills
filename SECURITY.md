# 安全策略 / Security Policy

## 支持的版本 / Supported Versions

| 版本 | 状态 |
|------|------|
| latest | ✅ 支持 |
| < latest | ⚠️ 仅安全修复 |

## 报告漏洞 / Reporting a Vulnerability

**请不要通过公开 Issue 报告安全漏洞。**

请通过 GitHub 的 [Private vulnerability reporting](https://docs.github.com/en/code-security/security-advisories)
功能私下提交,或联系维护者 `ssinder`。

请在报告中说明:

- 受影响的技能与版本
- 漏洞类型与影响范围
- 复现步骤或最小复现样例
- 你认为的修复建议

我们会在收到报告后 72 小时内确认,并在 7 天内给出修复计划或缓解方案。

## 安全设计原则 / Security Design Principles

本项目遵循以下安全设计原则:

### 1. 技能不执行未声明的网络行为

所有联网脚本必须显式传入参数才发起请求。技能在默认状态下不产生外部网络流量。

### 2. 不收集用户数据

技能仅在用户指定的本地目录读写文件,不向任何远程服务上传数据,除非用户明确配置了相应 API 端点。

### 3. 凭据仅从环境变量读取

所有 API 密钥、令牌通过环境变量传入,绝不硬编码在代码中。CI 会扫描硬编码凭据。

```bash
export NCBI_API_KEY="your-key"
export OPENALEX_MAILTO="you@example.com"
```

### 4. 路径隔离

技能只写入用户指定的工作目录,路径通过参数传递,不越界写入系统目录。

### 5. 依赖审计

每次 PR 都运行 `pip-audit` 依赖漏洞扫描,发现高危漏洞将阻止合并。

## 适用范围 / Scope

以下内容在安全审查范围内:

- `skills/` 下所有技能脚本
- `packages/` 共享工具层
- `scripts/` 校验与构建脚本
- CI 工作流配置

## 免责声明 / Disclaimer

本项目提供的是分析与计算工具,其输出结果需要研究者自行验证。工具的设计不构成任何医学、
临床或监管建议。
