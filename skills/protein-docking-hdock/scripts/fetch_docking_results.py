#!/usr/bin/env python3
"""获取 HDOCK 对接结果。

HDOCK 服务端计算通常需要 5–30 分钟。这个命令可以：

- **单次查询**（默认）：立刻看一次当前状态，不等待
- **等待完成**（``--wait``）：轮询到作业完成或超时

对应老仓库 ``dock_hdock/scripts/hdock score get can exist.py``——那个脚本
用 Selenium 逐个访问作业页，用 XPath ``//th[text()='Docking Score']``
取分数。这里改成纯 HTTP + 正则解析，并且**区分"排队中"和"已完成"**：
老脚本在作业还没算完时会取到空值并直接中止整个循环。

Examples:
    # 查一次当前状态
    python fetch_docking_results.py --url http://hdock.phys.hust.edu.cn/data/abc123/

    # 等到算完
    python fetch_docking_results.py --url http://hdock.phys.hust.edu.cn/data/abc123/ --wait

    # 批量查 jobs.csv 里记录的所有作业
    python fetch_docking_results.py --jobs-file results/jobs.csv --output-dir results/
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Optional

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "packages"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from biopipe_core import ensure_dir, get_logger  # noqa: E402
from hdock_client import (  # noqa: E402
    JobStatus,
    download_top_models,
    fetch_result_page,
    parse_result_page,
    wait_for_results,
)

LOGGER = get_logger("biopipe.protein-docking-hdock")

FIELDNAMES = [
    "job_name",
    "job_id",
    "url",
    "status",
    "docking_score",
    "confidence_score",
    "error",
]


def build_parser() -> argparse.ArgumentParser:
    """构造命令行参数解析器。"""
    parser = argparse.ArgumentParser(
        prog="fetch_docking_results",
        description="获取 HDOCK 对接结果（分数、状态、可下载模型）",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "示例:\n"
            "  python fetch_docking_results.py --url http://hdock.phys.hust.edu.cn/data/abc/\n"
            "  python fetch_docking_results.py --url http://hdock.phys.hust.edu.cn/data/abc/ --wait\n"
            "  python fetch_docking_results.py --jobs-file results/jobs.csv --output-dir results/\n"
            "\n"
            "说明:\n"
            "  · 默认只查一次，立刻返回当前状态。\n"
            "  · --wait 会轮询到完成（默认最多等 60 分钟）。\n"
            "  · 排队中的作业不会被当成失败，会在结果里标记 status=queued。\n"
        ),
    )
    parser.add_argument("--url", help="单个作业的结果页 URL")
    parser.add_argument(
        "--jobs-file",
        help="作业记录 CSV（submit_docking.py 生成），批量查询其中的作业",
    )
    parser.add_argument(
        "--output-dir", default="hdock_results", help="输出目录（默认: hdock_results）"
    )
    parser.add_argument("--wait", action="store_true", help="轮询等待作业完成")
    parser.add_argument(
        "--poll-interval", type=int, default=20, help="轮询间隔秒数（默认: 20）"
    )
    parser.add_argument(
        "--max-wait", type=int, default=3600, help="最长等待秒数（默认: 3600）"
    )
    parser.add_argument(
        "--download-models",
        action="store_true",
        help="下载排名靠前的复合物模型文件",
    )
    parser.add_argument(
        "--top-n", type=int, default=5, help="下载前 N 个模型（默认: 5）"
    )
    parser.add_argument(
        "--timeout", type=int, default=60, help="单次请求超时秒数（默认: 60）"
    )
    return parser


def _read_jobs_file(path: Path) -> list[dict]:
    """从 jobs.csv 读取作业记录。"""
    with open(path, encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _report(result) -> None:
    """打印单个作业的状态。"""
    if result.status == JobStatus.DONE:
        LOGGER.info(
            "✔ %s —— Docking Score=%s, Confidence Score=%s",
            result.job_name or result.job_id,
            result.docking_score,
            result.confidence_score,
        )
    elif result.status == JobStatus.QUEUED:
        LOGGER.info("⏳ %s —— 仍在排队（HDOCK 服务端计算中）", result.job_name or result.job_id)
    elif result.status == JobStatus.FAILED:
        LOGGER.warning(
            "✘ %s —— 失败: %s", result.job_name or result.job_id, "; ".join(result.errors)
        )
    else:
        LOGGER.warning("? %s —— 状态未知", result.job_name or result.job_id)


def main(argv: Optional[list[str]] = None) -> int:
    """命令行主入口。"""
    args = build_parser().parse_args(argv)

    if not args.url and not args.jobs_file:
        LOGGER.error("必须提供 --url 或 --jobs-file")
        return 2

    targets: list[tuple[str, str]] = []
    if args.url:
        targets.append((Path(args.url).name.strip("/") or "job", args.url))
    if args.jobs_file:
        jobs_path = Path(args.jobs_file)
        if not jobs_path.is_file():
            LOGGER.error("作业文件不存在: %s", jobs_path)
            return 2
        for row in _read_jobs_file(jobs_path):
            url = (row.get("url") or "").strip()
            if url:
                targets.append((row.get("job_name") or "job", url))

    if not targets:
        LOGGER.error("没有可查询的作业")
        return 2

    output_dir = ensure_dir(Path(args.output_dir))
    rows: list[dict] = []

    for job_name, url in targets:
        LOGGER.info("查询作业 %s: %s", job_name, url)

        if args.wait:
            def _on_poll(result) -> None:
                LOGGER.info(
                    "  状态: %s", result.status
                )

            result = wait_for_results(
                url,
                job_name=job_name,
                poll_interval=args.poll_interval,
                max_wait=args.max_wait,
                timeout=args.timeout,
                on_poll=_on_poll,
            )
        else:
            html = fetch_result_page(url, timeout=args.timeout)
            result = parse_result_page(html, job_name=job_name, url=url)

        _report(result)
        rows.append(result.to_row())

        if args.download_models and result.status == JobStatus.DONE:
            model_dir = output_dir / "models" / (result.job_id or job_name)
            downloaded = download_top_models(
                url, model_dir, top_n=args.top_n, timeout=args.timeout
            )
            LOGGER.info("下载了 %d 个模型文件到 %s", len(downloaded), model_dir)

    csv_path = output_dir / "results.csv"
    with open(csv_path, "w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)
    LOGGER.info("结果已写入: %s", csv_path)

    json_path = output_dir / "results.json"
    json_path.write_text(
        json.dumps(rows, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    LOGGER.info("结果已写入: %s", json_path)

    done = sum(1 for row in rows if row["status"] == JobStatus.DONE)
    queued = sum(1 for row in rows if row["status"] == JobStatus.QUEUED)
    failed = sum(1 for row in rows if row["status"] == JobStatus.FAILED)
    LOGGER.info("共 %d 个作业：完成 %d，排队中 %d，失败 %d", len(rows), done, queued, failed)
    return 0


if __name__ == "__main__":
    sys.exit(main())
