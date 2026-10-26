#!/usr/bin/env python3
"""批量下载 ENA 上的 FASTQ 数据（用于 RNA-seq 原始数据获取）。

核心改进（相对老仓库 ``rnaseq_download/crawl_srr_links.py``）：

1. **用 ENA Portal API 取代 HTML 爬虫**。老脚本用 lxml 解析 NCBI SRA
   页面，依赖脆弱的 XPath 定位、只能拿到 ``.sra`` 原始文件、文件大小靠
   抓页面文本。ENA API 直接给出 ``fastq_ftp`` / ``fastq_bytes`` /
   ``fastq_md5``，结构化且稳定。

2. **不在 import 期执行**。老脚本的所有逻辑都在模块顶层，import 就会
   建目录、读文件、爬网页、下载数据。本脚本所有动作都在 ``main()`` 里。

3. **代理可用**。老脚本设 ``session.trust_env = False`` 显式禁用系统
   代理，在需要代理的网络里直接不可用。本脚本沿用 ``http_session``，
   默认信任环境变量。

4. **断点续传**。老脚本用 ``curl`` 子进程且不带 ``-C -``，中断后只能
   从零重来；本脚本用 HTTP Range 头续传，并按 MD5 校验完整性。

5. **MD5 完整性校验**。ENA 提供每个 FASTQ 的 MD5，下载后核对，
   避免把截断的文件当成完整数据用于后续定量。

Examples:
    # 查询元数据并下载 FASTQ
    python download_runs.py --ids SRR390728 SRR1975568 --output-dir fastq/

    # 从文件读取 accession 列表
    python download_runs.py --ids runs.txt --output-dir fastq/

    # 只查元数据不下载
    python download_runs.py --ids runs.txt --metadata-only --output-dir meta/

    # 只下载单端（每个 run 的第一个 FASTQ）
    python download_runs.py --ids runs.txt --output-dir fastq/ --single-end
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Optional

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "packages"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from biopipe_core import ensure_dir, get_logger  # noqa: E402
from ena_client import build_fastq_url, fetch_run_metadata  # noqa: E402

LOGGER = get_logger("biopipe.ena-download")

# 每次读取的块大小（流式下载）
CHUNK_SIZE = 1024 * 1024


def read_id_file(path: Path) -> list[str]:
    """从文本文件读取 run accession 列表。

    老脚本依赖一个「index」文件里"两行一组"的约定（文件路径 + base URL），
    配置与数据混在一起且路径写死。本技能改为直接读取纯 accession 列表。

    Args:
        path: 文件路径，每行一个 accession。

    Returns:
        去重后的 accession 列表（保持原顺序）。
    """
    if not path.is_file():
        raise FileNotFoundError(f"accession 列表文件不存在: {path}")

    result: list[str] = []
    seen: set[str] = set()
    for line in path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        candidate = stripped.split()[0] if " " in stripped else stripped
        candidate = candidate.rstrip(",;")
        if not candidate or candidate.upper() in seen:
            continue
        seen.add(candidate.upper())
        result.append(candidate)
    return result


def _split_multi(value: str) -> list[str]:
    """拆分 ENA 的多值字段（分号分隔）。

    ENA 的 ``fastq_ftp`` / ``fastq_bytes`` / ``fastq_md5`` 是分号分隔的
    并行数组——第 i 个文件的大小与 MD5 分别对应第 i 个元素。
    """
    if not value:
        return []
    return [item.strip() for item in value.split(";") if item.strip()]


def compute_md5(path: Path, chunk_size: int = CHUNK_SIZE) -> str:
    """计算文件的 MD5。"""
    digest = hashlib.md5()  # noqa: S324 - ENA 用 MD5 做数据校验，非安全用途
    with path.open("rb") as handle:
        while True:
            block = handle.read(chunk_size)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def download_file(
    url: str,
    destination: Path,
    session=None,
    timeout: int = 600,
    expected_md5: Optional[str] = None,
    resume: bool = True,
) -> tuple[bool, str]:
    """下载单个文件，支持断点续传与 MD5 校验。

    Args:
        url: 下载地址（HTTPS）。
        destination: 本地目标路径。
        session: 可复用的 ``requests.Session``。
        timeout: 单次请求超时秒数。
        expected_md5: ENA 提供的 MD5；给出时会校验并删除损坏文件。
        resume: 是否允许断点续传。

    Returns:
        ``(是否成功, 说明信息)``。
    """
    from biopipe_core import http_session

    own_session = session is None
    if own_session:
        session = http_session(min_interval=0.34, retries=3)
    if session is None:
        return False, "requests 不可用"

    destination.parent.mkdir(parents=True, exist_ok=True)

    # 已有部分文件时尝试续传
    offset = destination.stat().st_size if (resume and destination.exists()) else 0
    headers: dict[str, str] = {}
    if offset > 0:
        headers["Range"] = f"bytes={offset}-"

    try:
        response = session.get(url, timeout=timeout, headers=headers, stream=True)

        # 服务端不支持 Range 时会返回 200 而不是 206，此时必须从头写
        if offset > 0 and response.status_code == 200:
            LOGGER.info("服务端不支持续传，重新下载: %s", destination.name)
            offset = 0

        response.raise_for_status()

        mode = "ab" if offset > 0 else "wb"
        written = offset
        with destination.open(mode) as handle:
            for block in response.iter_content(chunk_size=CHUNK_SIZE):
                if not block:
                    continue
                handle.write(block)
                written += len(block)
    except Exception as exc:  # noqa: BLE001
        return False, str(exc)
    finally:
        if own_session:
            session.close()

    actual_size = destination.stat().st_size

    if expected_md5:
        actual_md5 = compute_md5(destination)
        if actual_md5.lower() != expected_md5.lower():
            LOGGER.error("MD5 不匹配: %s（期望 %s，实际 %s）", destination.name, expected_md5, actual_md5)
            # 文件已损坏，删掉避免下次误判为完整
            try:
                destination.unlink()
            except OSError:  # pragma: no cover
                pass
            return False, f"MD5 不匹配（期望 {expected_md5}，实际 {actual_md5}）"

    return True, f"{actual_size} 字节"


def write_metadata(records: Sequence[dict], destination: Path) -> None:
    """把 run 元数据写成 CSV 和 JSON。"""
    destination.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "run_accession",
        "experiment_accession",
        "sample_accession",
        "study_accession",
        "scientific_name",
        "instrument_model",
        "library_strategy",
        "library_source",
        "library_layout",
        "read_count",
        "base_count",
        "fastq_ftp",
        "fastq_bytes",
        "fastq_md5",
    ]

    csv_path = destination.with_suffix(".csv")
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for record in records:
            writer.writerow(record)

    json_path = destination.with_suffix(".json")
    json_path.write_text(
        json.dumps(list(records), ensure_ascii=False, indent=2), encoding="utf-8"
    )

    LOGGER.info("元数据已写入: %s / %s", csv_path, json_path)


def build_parser() -> argparse.ArgumentParser:
    """构造命令行参数解析器。"""
    parser = argparse.ArgumentParser(
        prog="download_runs",
        description="从 ENA 批量下载 FASTQ 数据并导出 run 元数据",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "示例:\n"
            "  python download_runs.py --ids SRR390728 --output-dir fastq/\n"
            "  python download_runs.py --ids runs.txt --output-dir fastq/\n"
            "  python download_runs.py --ids runs.txt --metadata-only --output-dir meta/\n"
            "\n"
            "说明:\n"
            "  元数据来自 ENA Portal API（结构化），不解析 HTML 页面，\n"
            "  因此不受网页改版影响；下载后按 ENA 提供的 MD5 校验完整性。\n"
        ),
    )
    parser.add_argument(
        "--ids",
        nargs="+",
        required=True,
        help="run accession 列表（SRR/ERR/DRR），或包含 accession 的文本文件路径",
    )
    parser.add_argument(
        "--output-dir",
        default="ena_fastq",
        help="下载目录（默认: ena_fastq）",
    )
    parser.add_argument(
        "--metadata-only",
        action="store_true",
        help="只查询并导出元数据，不下载文件",
    )
    parser.add_argument(
        "--single-end",
        action="store_true",
        help="每个 run 只下载第一个 FASTQ（单端数据场景）",
    )
    parser.add_argument(
        "--strategy",
        default=None,
        help="只处理指定文库策略的 run（如 RNA-Seq、WGS）",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="已存在的文件也重新下载",
    )
    parser.add_argument(
        "--no-resume",
        action="store_true",
        help="禁用断点续传（中断后从头下载）",
    )
    parser.add_argument(
        "--no-verify",
        action="store_true",
        help="跳过 MD5 校验（不推荐）",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=600,
        help="单次请求超时秒数（默认: 600）",
    )
    return parser


def _collect_accessions(items: Sequence[str]) -> list[str]:
    """把 --ids 参数解析成 accession 列表。"""
    accessions: list[str] = []
    for item in items:
        try:
            is_file = Path(item).is_file()
        except OSError:  # pragma: no cover
            is_file = False
        if is_file:
            LOGGER.info("从文件读取 accession: %s", item)
            accessions.extend(read_id_file(Path(item)))
        else:
            accessions.append(item.strip())
    return accessions


def main(argv: Optional[list[str]] = None) -> int:
    """命令行主入口。"""
    args = build_parser().parse_args(argv)

    accessions = _collect_accessions(args.ids)
    if not accessions:
        LOGGER.error("没有可用的 run accession")
        return 2

    LOGGER.info("开始查询 %d 个 run 的 ENA 元数据", len(accessions))
    records = fetch_run_metadata(accessions, timeout=args.timeout)

    if not records:
        LOGGER.error("没有查到任何 ENA 记录")
        return 1

    if args.strategy:
        wanted = args.strategy.lower()
        filtered = [
            r for r in records if (r.get("library_strategy") or "").lower() == wanted
        ]
        LOGGER.info(
            "按 library_strategy=%s 过滤: %d/%d 条保留",
            args.strategy,
            len(filtered),
            len(records),
        )
        records = filtered
        if not records:
            LOGGER.error("过滤后没有剩余记录")
            return 1

    output_dir = ensure_dir(args.output_dir)
    write_metadata(records, output_dir / "ena_run_metadata.csv")

    if args.metadata_only:
        LOGGER.info("--metadata-only 已指定，跳过文件下载")
        LOGGER.info("共 %d 条记录，元数据在 %s", len(records), output_dir)
        return 0

    from biopipe_core import http_session

    session = http_session(min_interval=0.34, retries=3)
    if session is None:
        LOGGER.error("requests 不可用，无法下载文件")
        return 1

    stats = {"downloaded": 0, "skipped": 0, "failed": 0, "files": 0}
    failed_runs: list[str] = []

    try:
        for index, record in enumerate(records, 1):
            run = record.get("run_accession", "")
            LOGGER.info("[%d/%d] %s", index, len(records), run)

            ftp_paths = _split_multi(record.get("fastq_ftp", ""))
            sizes = _split_multi(record.get("fastq_bytes", ""))
            md5s = _split_multi(record.get("fastq_md5", ""))

            if not ftp_paths:
                LOGGER.warning("%s 没有可用的 FASTQ（可能未转换或被删除）", run)
                stats["failed"] += 1
                failed_runs.append(run)
                continue

            if args.single_end:
                ftp_paths = ftp_paths[:1]
                sizes = sizes[:1]
                md5s = md5s[:1]

            run_failed = False
            for file_index, ftp_path in enumerate(ftp_paths):
                url = build_fastq_url(ftp_path)
                filename = Path(ftp_path).name
                destination = output_dir / filename

                expected_md5 = md5s[file_index] if file_index < len(md5s) else None
                if args.no_verify:
                    expected_md5 = None

                if destination.exists() and not args.force:
                    LOGGER.info("已存在，跳过: %s", filename)
                    stats["skipped"] += 1
                    continue

                expected_size = ""
                if file_index < len(sizes):
                    try:
                        expected_size = f"{int(sizes[file_index]) / (1024**3):.2f} GB"
                    except ValueError:
                        expected_size = f"{sizes[file_index]} 字节"

                ok, detail = download_file(
                    url,
                    destination,
                    session=session,
                    timeout=args.timeout,
                    expected_md5=expected_md5,
                    resume=not args.no_resume,
                )
                if ok:
                    if expected_size:
                        LOGGER.info("下载成功: %s (%s，预期 %s)", filename, detail, expected_size)
                    else:
                        LOGGER.info("下载成功: %s (%s)", filename, detail)
                    stats["downloaded"] += 1
                    stats["files"] += 1
                else:
                    LOGGER.error("下载失败: %s —— %s", filename, detail)
                    stats["failed"] += 1
                    run_failed = True

            if run_failed:
                failed_runs.append(run)
    finally:
        session.close()

    LOGGER.info("-" * 56)
    LOGGER.info(
        "下载完成：新下载 %d 个文件，跳过 %d，失败 %d",
        stats["downloaded"],
        stats["skipped"],
        stats["failed"],
    )
    LOGGER.info("输出目录: %s", output_dir)

    if failed_runs:
        LOGGER.warning("失败的 run: %s", ", ".join(failed_runs))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
