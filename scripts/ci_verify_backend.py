#!/usr/bin/env python3
"""在 CI 中执行双后端一致性验证，并确保它**真的跑了**。

为什么需要这个脚本
------------------
对照测试的价值完全取决于"它是否真的执行了"。如果环境里缺少参考实现，
pytest 会**跳过**这些用例——CI 显示绿色，但什么都没验证。那是一个
永不失败的摆设徽章，比没有 CI 更危险：它让人误以为后端一致性已被保证。

因此这个脚本做两件事：

1. 先探测参考实现是否可用（DESeq2 / scanpy）；
2. 不可用则**直接失败**（而不是跳过）；
3. 可用则跑 pytest，并把跳过数也视为失败——因为对照用例一旦被跳过，
   就说明这次验证没有覆盖到目标。

用法：
    python scripts/ci_verify_backend.py --reference deseq2 --test tests/test_backend_agreement.py
    python scripts/ci_verify_backend.py --reference scanpy --test tests/test_scrna_backend_agreement.py
"""

from __future__ import annotations

import argparse
import subprocess
import sys

# 参考实现的探测命令。
#
# 每个参考实现对应一条 R / Python 探测语句，返回 "TRUE" 视为可用。
PROBES = {
    "deseq2": {
        "cmd": ["Rscript", "-e", "cat(requireNamespace('DESeq2', quietly=TRUE))"],
        "label": "Bioconductor DESeq2",
    },
    "scanpy": {
        "cmd": [sys.executable, "-c", "import scanpy, anndata; print('TRUE')"],
        "label": "scanpy + anndata",
    },
}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--reference",
        required=True,
        choices=sorted(PROBES),
        help="参考实现（deseq2 / scanpy）",
    )
    parser.add_argument("--test", required=True, help="要执行的 pytest 文件")
    args = parser.parse_args(argv)

    probe = PROBES[args.reference]

    # 1. 探测参考实现
    try:
        result = subprocess.run(
            probe["cmd"],
            capture_output=True,
            text=True,
            timeout=300,
        )
    except (OSError, subprocess.SubprocessError) as exc:
        print(f"::error::无法探测 {probe['label']}: {exc}")
        return 1

    if "TRUE" not in result.stdout:
        print(f"::error::{probe['label']} 不可用，后端一致性验证无法进行")
        if result.stderr.strip():
            print(result.stderr.strip()[-1500:])
        return 1

    print(f"{probe['label']} 可用，开始执行对照测试")

    # 2. 跑测试。跳过数也算失败——对照用例被跳过意味着没验证到。
    proc = subprocess.run(
        [
            sys.executable,
            "-m",
            "pytest",
            args.test,
            "-v",
            "--no-header",
            "-rs",  # 打印跳过原因
        ],
        capture_output=True,
        text=True,
        timeout=3600,
    )

    output = proc.stdout + proc.stderr
    print(output[-6000:])

    if proc.returncode != 0:
        return proc.returncode

    # 3. 检查有没有被跳过的用例
    skip_marker = "skipped"
    summary = [line for line in output.splitlines() if " passed" in line or " skipped" in line]
    if summary:
        last = summary[-1]
        # 形如 "5 passed, 2 skipped in 12.3s"
        if skip_marker in last and "0 skipped" not in last:
            print(f"::error::对照测试有跳过用例，说明未被真正执行: {last}")
            return 1

    print("对照测试全部执行并通过")
    return 0


if __name__ == "__main__":
    sys.exit(main())
