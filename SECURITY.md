# Security Policy

## Supported Versions

| Version | Status |
|---------|--------|
| latest | ✅ Supported |
| < latest | ⚠️ Security fixes only |

## Reporting a Vulnerability

**Please do not report security vulnerabilities through public Issues.**

Use GitHub's [Private vulnerability reporting](https://docs.github.com/en/code-security/security-advisories)
feature to submit a report privately, or contact the maintainer `ssinder`.

Please include the following in your report:

- The affected skill(s) and version(s)
- The type of the vulnerability and its scope of impact
- Steps to reproduce, or a minimal reproduction case
- Your suggested remediation, if you have one

We will acknowledge the report within 72 hours, and provide a fix plan or mitigation within
7 days.

## Security Design Principles

This project follows the following security design principles:

### 1. Skills perform no undeclared network activity

Every networked script requires its arguments to be passed explicitly before any request is made.
In their default state, skills generate no outbound network traffic.

### 2. No user data is collected

Skills read and write files only in the local directory you specify, and never upload data to any
remote service unless you have explicitly configured the corresponding API endpoint.

### 3. Credentials are read from environment variables only

All API keys and tokens are passed in through environment variables, and are never hardcoded in the
source. CI scans for hardcoded credentials.

```bash
export NCBI_API_KEY="your-key"
export OPENALEX_MAILTO="you@example.com"
```

### 4. Path isolation

Skills write only to the working directory you specify. Paths are passed in as arguments, and the
skills never write outside that boundary into system directories.

### 5. Dependency auditing

Every PR runs a `pip-audit` dependency vulnerability scan, and a high-severity finding blocks the
merge.

## Scope

The following are in scope for security review:

- All skill scripts under `skills/`
- The shared utility layer under `packages/`
- Validation and build scripts under `scripts/`
- CI workflow configuration

## Disclaimer

This project provides analysis and computation tools. The results they produce must be verified
independently by the researcher. The tools are not designed to provide medical, clinical, or
regulatory advice.
