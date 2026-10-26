#!/usr/bin/env python3
"""ENA（European Nucleotide Archive）数据检索客户端。

老仓库的 ``rnaseq_download/crawl_srr_links.py`` 用 HTML 爬虫从 NCBI
SRA 页面提取下载链接与文件大小。这套做法的问题：

- 依赖脆弱的 XPath（``/html/body/div[2]/div[2]/div[2]/...``），页面
  一改版就失效；
- 只能拿到 ``.sra`` 原始文件，拿不到已经处理好的 FASTQ；
- 文件大小靠抓页面文本，格式随语言/版本变化；
- 关闭了代理（``trust_env = False``），在需要代理的网络里直接不可用。

ENA Portal API 直接提供结构化的 ``fastq_ftp`` / ``fastq_bytes`` /
``fastq_md5`` 字段——不必解析 HTML，也不必二次转换 SRA。

接口
----
- :func:`fetch_run_metadata` —— 批量查询 run 元数据（search 接口，支持 OR 批量）
- :func:`parse_tsv` —— 解析 ENA 的 TSV 响应
- :func:`build_fastq_url` —— 把 ENA 返回的 ftp 路径转成 HTTPS 下载地址

Usage:
    from ena_client import fetch_run_metadata
    rows = fetch_run_metadata(["SRR390728", "SRR1975568"])
"""

from __future__ import annotations

import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Optional

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "packages"))

from biopipe_core import get_logger, http_session  # noqa: E402

__all__ = [
    "fetch_run_metadata",
    "parse_tsv",
    "build_fastq_url",
    "ENA_SEARCH_URL",
    "RUN_FIELDS",
]

ENA_SEARCH_URL = "https://www.ebi.ac.uk/ena/portal/api/search"

LOGGER = get_logger("biopipe.ena-download")

# 查询 run 时请求的字段。fastq_* 是 ENA 已转换好的 FASTQ，
# 比 submitted_*（原始提交文件，通常是 .sra）更常用。
RUN_FIELDS = (
    "run_accession",
    "experiment_accession",
    "sample_accession",
    "study_accession",
    "secondary_study_accession",
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
    "submitted_ftp",
    "submitted_bytes",
)


def parse_tsv(payload: str) -> list[dict]:
    """解析 ENA 的 TSV 响应。

    ENA 在没有任何记录时返回**只有表头的一行**（以 ``\\n`` 结尾），
    此时应返回空列表而不是一行全空的数据。

    Args:
        payload: TSV 文本。

    Returns:
        每行一个字典的列表。
    """
    lines = [line for line in payload.splitlines() if line.strip()]
    if not lines:
        return []

    header = lines[0].split("\t")
    rows: list[dict] = []
    for line in lines[1:]:
        values = line.split("\t")
        # 列数不足时补空串，避免后续取值 IndexError
        values += [""] * (len(header) - len(values))
        rows.append(dict(zip(header, values)))
    return rows


def build_fastq_url(ftp_path: str) -> str:
    """把 ENA 返回的 FTP 路径转成 HTTPS 下载地址。

    ENA 返回形如 ``ftp.sra.ebi.ac.uk/vol1/fastq/SRR390/SRR390728/...``，
    HTTPS 端点与之同路径，直接替换协议即可——这样不用 FTP 客户端，
    也能复用 HTTP 的重试与代理设置。

    Args:
        ftp_path: ENA 返回的 FTP 路径。

    Returns:
        HTTPS URL；若输入为空则返回空串。
    """
    if not ftp_path:
        return ""
    if ftp_path.startswith(("http://", "https://")):
        return ftp_path
    return f"https://{ftp_path}"


def fetch_run_metadata(
    accessions: Sequence[str],
    session=None,
    timeout: int = 120,
    batch_size: int = 50,
) -> list[dict]:
    """批量查询 run 元数据。

    ENA 的 ``filereport`` 接口**不支持**逗号/分号分隔的批量 accession
    （实测会返回一段说明支持的 accession 类型的纯文本），所以这里用
    ``search`` 接口，把多个 accession 用 ``OR`` 拼成一个查询串。

    Args:
        accessions: SRR/ERR/DRR 等 run accession 列表。
        session: 可复用的 ``requests.Session``。
        timeout: 单次请求超时秒数。
        batch_size: 单个查询串里最多放多少个 accession。ENA 对超长
            查询串有限制，分批更稳妥。

    Returns:
        每个 run 一条记录的列表。字段缺失时为空串。
    """
    unique: list[str] = []
    seen: set[str] = set()
    for acc in accessions:
        cleaned = (acc or "").strip()
        if not cleaned or cleaned.upper() in seen:
            continue
        seen.add(cleaned.upper())
        unique.append(cleaned)

    if not unique:
        return []

    own_session = session is None
    if own_session:
        session = http_session(min_interval=0.34, retries=3)
    if session is None:
        LOGGER.warning("requests 不可用，无法查询 ENA")
        return []

    results: list[dict] = []
    try:
        for start in range(0, len(unique), batch_size):
            chunk = unique[start : start + batch_size]
            query = " OR ".join(f"run_accession={acc}" for acc in chunk)
            params = {
                "result": "read_run",
                "query": query,
                "fields": ",".join(RUN_FIELDS),
                "format": "tsv",
                # limit=0 表示不限制返回条数
                "limit": "0",
            }
            LOGGER.info("查询 ENA: %d 个 run accession", len(chunk))
            try:
                response = session.get(ENA_SEARCH_URL, params=params, timeout=timeout)
                response.raise_for_status()
                payload = response.text
            except Exception as exc:  # noqa: BLE001
                LOGGER.warning("ENA 查询失败（%d 个 accession）: %s", len(chunk), exc)
                continue

            rows = parse_tsv(payload)
            if not rows:
                LOGGER.warning("ENA 未返回记录（%d 个 accession）", len(chunk))
            results.extend(rows)
    finally:
        if own_session:
            session.close()

    LOGGER.info("ENA 返回 %d/%d 条记录", len(results), len(unique))
    return results


def _main(argv: Optional[list[str]] = None) -> int:
    """命令行入口：查询 ENA 并打印结果（自检工具）。"""
    import argparse
    import json

    parser = argparse.ArgumentParser(
        prog="ena_client", description="查询 ENA run 元数据（命令行自检工具）"
    )
    parser.add_argument("accessions", nargs="+", help="run accession（SRR/ERR/DRR）")
    parser.add_argument("--json", action="store_true", help="以 JSON 输出完整记录")
    args = parser.parse_args(argv)

    rows = fetch_run_metadata(args.accessions)
    if not rows:
        LOGGER.warning("没有查到任何记录")
        return 1

    if args.json:
        print(json.dumps(rows, ensure_ascii=False, indent=2))
    else:
        for row in rows:
            LOGGER.info(
                "%s | %s | %s | %s",
                row.get("run_accession", ""),
                row.get("scientific_name", ""),
                row.get("library_strategy", ""),
                row.get("fastq_ftp", "")[:100],
            )
    return 0


if __name__ == "__main__":
    sys.exit(_main())
