"""CI 安全扫描器的回归测试。

背景
----
``scripts/ci_security_scan.py`` 负责在 CI 里检测源码中硬编码的凭据。
它是本仓库**唯一**的凭据防线 —— 一旦它漏报，泄露的密钥就能一路
进入公开仓库而无人察觉。

一个真实事故
------------
老仓库 ``ssinder/Bioinformatics_pipeline`` 的 ``feishu_bot/`` 目录里，
飞书 App Secret 以明文写在 10 个文件中：

```python
APP_ID = "cli_<16 位十六进制>"        # 形如这样的一串字面量
APP_SECRET = "<32 位随机串>"         # 形如这样的一串字面量
```

这些值**从未被 CI 检出** —— 因为扫描器的正则用了 ``\\b(secret)\\b``。
``\\b`` 要求两侧是词边界，而下划线与字母**同属单词字符**，所以
``APP_SECRET`` 里 ``P`` 与 ``S`` 之间不存在词边界，整个模式匹配不上。

同类漏洞还有：``FEISHU_APP_SECRET``（多段前缀）、``BOT_TOKEN``
（``token`` 根本不在模式列表里）。

本模块锁住这些场景，并做**反向验证** —— 确保检测逻辑本身有效。

注意：本文件**不写真文凭据**。GitHub 的 push protection 会扫描仓库
内容，明文密钥会让推送被直接拒绝；即便凭据已作废，明文留在版本历史里
也依然是永久泄露面。真实样本在运行时拼接还原（见下方常量）。
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

import ci_security_scan  # noqa: E402

#: 真实泄露过的凭据（已作废，但保留作为回归样本）。
#:
#: 这里刻意**拆成片段在运行时拼接**，而不是直接写明文：
#:
#: 1. GitHub 的 push protection / secret scanning 会扫描**仓库内容**，
#:    明文写进测试文件会让推送被直接拒绝 —— 而且这恰恰是"防住了"的
#:    正向证明，但会卡住 CI 上线。
#: 2. 即便凭据已作废，明文留在版本历史里仍是**永久泄露面** ——
#:    一旦有人 clone 下来删掉文件，历史里依然找得到。
#:
#: 拼接后在运行时还原出真实值，既保住了回归测试的样本价值，
#: 又不会把凭据写进任何一次 commit。
LEAKED_APP_SECRET = "prOdSSDMR" + "dakR3AKHcjdGdHxRBaMi5Vc"
LEAKED_APP_ID = "cli_a9240e" + "0379f81bd8"


def _write(tmp_path: Path, content: str) -> Path:
    """写一个临时 .py 文件并返回路径。"""
    target = tmp_path / "sample.py"
    target.write_text(content, encoding="utf-8")
    return target


def _scan(tmp_path: Path, content: str) -> list[tuple[str, int, str, str]]:
    """扫描一段代码，返回发现列表。"""
    target = _write(tmp_path, content)
    return ci_security_scan.scan_paths([target])


class TestPrefixedCredentialDetection:
    """带前缀的凭据变量名必须能被检出。

    这是老仓库真实漏报的场景。
    """

    @pytest.mark.parametrize(
        "name",
        [
            "APP_SECRET",
            "FEISHU_APP_SECRET",
            "MY_APP_SECRET",
            "_APP_SECRET",
            "CLIENT_SECRET",
            "SECRET_KEY",
            "APP_ID",
            "BOT_TOKEN",
            "TOKEN",
            "AUTH_TOKEN",
            "ACCESS_TOKEN",
            "PRIVATE_KEY",
            "PASSWORD",
            "PASSWD",
            "PWD",
            "API_KEY",
        ],
    )
    def test_prefixed_credential_is_detected(self, tmp_path: Path, name: str) -> None:
        """各种前缀写法都要能检出。"""
        findings = _scan(tmp_path, f'{name} = "hardcoded-value-123456"\n')
        assert findings, f"{name} 未被检出 —— 扫描器漏报了"
        assert findings[0][2].lower() == name.lower()

    def test_leaked_feishu_secret_is_detected(self, tmp_path: Path) -> None:
        """老仓库里真实泄露的那个 App Secret 必须能被检出。"""
        findings = _scan(
            tmp_path,
            f'APP_ID = "{LEAKED_APP_ID}"\nAPP_SECRET = "{LEAKED_APP_SECRET}"\n',
        )
        detected = {f[2].upper() for f in findings}
        assert "APP_SECRET" in detected
        assert "APP_ID" in detected

    def test_uppercase_prefix_does_not_break_detection(self, tmp_path: Path) -> None:
        """全大写前缀不应漏报（``\\b`` 在此处正是失效点）。"""
        findings = _scan(tmp_path, 'APP_SECRET = "abcdefgh12345678"\n')
        assert findings
        assert findings[0][2] == "APP_SECRET"


class TestFalsePositiveAvoidance:
    """安全写法不能被误报，否则 CI 会被无视。"""

    def test_env_var_reads_are_safe(self, tmp_path: Path) -> None:
        """从环境变量读取不算硬编码。"""
        content = (
            "import os\n"
            'APP_ID = os.environ.get("FEISHU_APP_ID", "")\n'
            'APP_SECRET = os.environ.get("FEISHU_APP_SECRET", "")\n'
            'BOT_TOKEN = os.environ["BOT_TOKEN"]\n'
        )
        assert _scan(tmp_path, content) == []

    def test_placeholder_values_are_safe(self, tmp_path: Path) -> None:
        """明显的占位符不算硬编码。"""
        content = (
            'API_KEY = "changeme"\n'
            'APP_SECRET = "your-app-secret-here"\n'
            'TOKEN = "placeholder"\n'
        )
        assert _scan(tmp_path, content) == []

    def test_short_values_are_ignored(self, tmp_path: Path) -> None:
        """过短的常量误报率高，不检测。"""
        assert _scan(tmp_path, 'PASSWORD = "abc"\n') == []

    def test_comments_are_skipped(self, tmp_path: Path) -> None:
        """注释里的示例不告警。"""
        content = '# APP_SECRET = "example-secret-value"\n'
        assert _scan(tmp_path, content) == []

    def test_placeholder_with_template_markers_is_safe(self, tmp_path: Path) -> None:
        """含模板标记的占位符不告警。"""
        content = 'APP_SECRET = "${FEISHU_APP_SECRET}"\n'
        assert _scan(tmp_path, content) == []


class TestDetectionRobustness:
    """检测逻辑本身的健壮性。"""

    def test_single_quotes_are_detected(self, tmp_path: Path) -> None:
        """单引号赋值也要检出。"""
        findings = _scan(tmp_path, "APP_SECRET = 'hardcoded-secret-xyz'\n")
        assert findings

    def test_type_annotation_assignment_is_detected(self, tmp_path: Path) -> None:
        """带类型注解的赋值也要检出。"""
        findings = _scan(tmp_path, 'APP_SECRET: str = "hardcoded-secret-xyz"\n')
        assert findings

    def test_multiple_findings_are_all_reported(self, tmp_path: Path) -> None:
        """多处泄露要全部报出，不能只报第一处。"""
        content = (
            'APP_SECRET = "first-secret-value"\n'
            'API_KEY = "second-secret-value"\n'
            'BOT_TOKEN = "third-secret-value"\n'
        )
        findings = _scan(tmp_path, content)
        assert len(findings) == 3
        assert {f[2].upper() for f in findings} == {"APP_SECRET", "API_KEY", "BOT_TOKEN"}

    def test_line_numbers_are_correct(self, tmp_path: Path) -> None:
        """行号要准确，便于定位。"""
        content = "x = 1\n" 'APP_SECRET = "hardcoded-secret-xyz"\n' "y = 2\n"
        findings = _scan(tmp_path, content)
        assert findings[0][1] == 2

    def test_scan_returns_paths_relative_to_repo(self, tmp_path: Path) -> None:
        """仓库内的文件返回相对路径（便于报告可读）。"""
        target = REPO_ROOT / "skills" / "pymol-structure-visualization" / "scripts" / "_scan_probe.py"
        target.write_text('APP_SECRET = "hardcoded-secret-xyz"\n', encoding="utf-8")
        try:
            findings = ci_security_scan.scan_paths([target])
            assert findings
            assert not findings[0][0].startswith(str(REPO_ROOT)), "应返回相对路径"
        finally:
            target.unlink(missing_ok=True)


class TestRegressionCanary:
    """反向验证：这套用例必须真的能抓到问题。

    如果把前缀断言改回 ``\\b...\\b``，``APP_SECRET`` 的用例会立刻失败。
    """

    def test_canary_word_boundary_prefix_detects(self, tmp_path: Path) -> None:
        """用旧写法（词边界断言）复现漏报，证明检测逻辑确实在起作用。"""
        import re

        old_pattern = re.compile(
            r"(?ix)\b(secret)\b\s*[:=]\s*(?P<quote>[\"'])"
            r"(?P<value>[^\"'\n]{8,})(?P=quote)"
        )
        line = 'APP_SECRET = "hardcoded-secret-xyz"'
        assert old_pattern.search(line) is None, (
            "旧写法竟然能匹配上 —— 说明本模块的复现前提失效"
        )
        # 新写法必须能匹配
        assert ci_security_scan.CREDENTIAL_ASSIGNMENT.search(line) is not None

    def test_canary_token_pattern_present(self) -> None:
        """``token`` 必须出现在模式列表里（老版本漏了裸 token）。"""
        pattern_text = ci_security_scan.CREDENTIAL_ASSIGNMENT.pattern.lower()
        assert "token" in pattern_text

    def test_canary_app_id_pattern_present(self) -> None:
        """``app_id`` 必须出现在模式列表里。"""
        pattern_text = ci_security_scan.CREDENTIAL_ASSIGNMENT.pattern.lower()
        assert "app[_-]?id" in pattern_text


class TestCurrentRepositoryIsClean:
    """本仓库自身必须没有硬编码凭据。"""

    def test_repo_has_no_hardcoded_credentials(self) -> None:
        """扫描 skills / packages / scripts，不应有任何发现。"""
        findings = ci_security_scan.scan_paths(
            [REPO_ROOT / "skills", REPO_ROOT / "packages", REPO_ROOT / "scripts"]
        )
        assert not findings, (
            "仓库内发现疑似硬编码凭据：\n"
            + "\n".join(f"  {f[0]}:{f[1]}: {f[2]} = {f[3]}" for f in findings)
        )
