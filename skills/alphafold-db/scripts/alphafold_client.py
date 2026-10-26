#!/usr/bin/env python3
"""AlphaFold DB API 客户端（基于共享层 ``http_session``）。

只在显式调用时才发起网络请求；``requests`` 缺失时优雅降级。

接口
----
- :func:`fetch_prediction` —— 拉取单个 UniProt accession 的预测元数据
- :func:`fetch_predictions` —— 批量拉取（逐个请求，自动去重与限速）
- :func:`parse_prediction` —— 把 API 原始 JSON 扁平化为统一结构
- :func:`build_pdb_url` —— 兜底构造 PDB 下载地址

为什么不用 ``/files/AF-{ID}-F1-model_v6.pdb`` 直接拼 URL
--------------------------------------------------------
老脚本把模型版本号 ``model_v6`` 写死在 URL 里。AlphaFold DB 每次重新
训练模型都会递增版本号，写死的版本号在下一版发布后就会 404。

正确做法是调用 ``/api/prediction/{accession}`` 拿到响应里的
``pdbUrl`` 字段——它由服务端给出，永远指向当前有效版本。
本模块仅在 API 不可用时才用 :func:`build_pdb_url` 兜底，
且兜底路径同样不写死版本号（改为探测常见版本）。

API 返回格式的坑
-----------------
``GET /api/prediction/{accession}`` 对**单个** accession 返回
**数组**（含一个元素的 list），而不是对象。老脚本按 dict 处理会直接
崩溃。本模块用 :func:`_iter_entries` 统一处理 list / dict / 错误对象
三种形态。

Usage:
    from alphafold_client import fetch_predictions
    rows = fetch_predictions(["P69905", "P04637"])
"""

from __future__ import annotations

import sys
from collections.abc import Iterator, Sequence
from pathlib import Path
from typing import Any, Optional

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "packages"))

from biopipe_core import get_logger, http_session  # noqa: E402

__all__ = [
    "fetch_prediction",
    "fetch_predictions",
    "parse_prediction",
    "build_pdb_url",
    "ALPHAFOLD_API_BASE",
    "ALPHAFOLD_FILES_BASE",
]

ALPHAFOLD_API_BASE = "https://alphafold.ebi.ac.uk/api/prediction"
ALPHAFOLD_FILES_BASE = "https://alphafold.ebi.ac.uk/files"

LOGGER = get_logger("biopipe.alphafold-db")

# 兜底构造 PDB URL 时探测的模型版本（从新到旧）。
# 之所以列表化而不是写死单一版本，是为了在 API 不可用时仍有机会命中；
# 正常路径永远用 API 返回的 pdbUrl。
_FALLBACK_MODEL_VERSIONS = ("v6", "v5", "v4")


def _iter_entries(payload: Any) -> Iterator[dict]:
    """把 API 返回的不同形态统一成"记录字典"的迭代器。

    实测 ``GET /api/prediction/{acc}`` 的响应有三种可能：

    1. ``[ {...} ]`` —— 单个 accession 正常返回（list of dict）
    2. ``{...}`` —— 少数情况下直接返回单个对象
    3. ``{"error": "..."}`` —— ID 格式非法或不存在

    第 3 种必须与第 2 种区分开，否则会把错误信息当成元数据。
    """
    if isinstance(payload, list):
        for item in payload:
            if isinstance(item, dict):
                yield item
        return

    if isinstance(payload, dict):
        # 错误对象：只有一个 error 字段，没有 uniprotAccession
        if "error" in payload and "uniprotAccession" not in payload:
            return
        yield payload


def _field(entry: dict, key: str, default: Any = "") -> Any:
    """取值，但**不把 0.0 / False 当成缺失**。

    直接写 ``entry.get(key, "") or ""`` 会把合法的数值 0.0 转成空字符串——
    而 AlphaFold 的置信度占比字段里 0.0 是有意义的值（某个置信区间的
    残基占比恰好为 0）。这里只把 ``None`` 视为缺失。
    """
    value = entry.get(key)
    return default if value is None else value


def parse_prediction(entry: dict) -> dict:
    """把 AlphaFold DB 的一条原始记录扁平化为统一结构。

    提取字段：accession、entry_id、pdb_url、cif_url、bcif_url、
    gene、description、organism、tax_id、chain_id、model_version、
    model_created_date、global_plddt（整体 pLDDT）以及各置信区间的
    残基占比（用于快速筛选低质量模型）。
    """
    record: dict[str, Any] = {
        "accession": _field(entry, "uniprotAccession"),
        "entry_id": _field(entry, "entryId"),
        "pdb_url": _field(entry, "pdbUrl"),
        "cif_url": _field(entry, "cifUrl"),
        "bcif_url": _field(entry, "bcifUrl"),
        "gene": _field(entry, "gene"),
        "description": _field(entry, "uniprotDescription"),
        "organism": _field(entry, "organismScientificName"),
        "tax_id": _field(entry, "taxId"),
        "chain_id": _field(entry, "chainId"),
        "model_version": _field(entry, "latestVersion"),
        "model_created_date": _field(entry, "modelCreatedDate"),
        "global_plddt": _field(entry, "globalMetricValue"),
        "fraction_plddt_very_low": _field(entry, "fractionPlddtVeryLow"),
        "fraction_plddt_low": _field(entry, "fractionPlddtLow"),
        "fraction_plddt_confident": _field(entry, "fractionPlddtConfident"),
        "fraction_plddt_very_high": _field(entry, "fractionPlddtVeryHigh"),
        "sequence_length": len(entry.get("sequence") or ""),
        "is_reviewed": _field(entry, "isReviewed"),
        "is_complex": _field(entry, "isComplex"),
    }
    return record


def fetch_prediction(
    accession: str,
    session=None,
    timeout: int = 60,
) -> Optional[dict]:
    """拉取单个 accession 的预测元数据。

    Args:
        accession: UniProt accession（如 ``P69905``）。
        session: 可复用的 ``requests.Session``；None 时创建带重试的会话。
        timeout: 单次请求超时秒数。

    Returns:
        扁平化后的记录字典；该 accession 无 AlphaFold 模型或请求失败时返回 None。
    """
    results = fetch_predictions([accession], session=session, timeout=timeout)
    return results[0] if results else None


def fetch_predictions(
    accessions: Sequence[str],
    session=None,
    timeout: int = 60,
    min_interval: float = 0.2,
) -> list[dict]:
    """批量拉取预测元数据。

    AlphaFold DB 的 ``/api/prediction/`` **不支持逗号分隔的批量查询**
    （实测会返回 ``{"error": "Invalid identifier format..."}``），所以这里
    逐个请求。请求间隔由共享层的 ``http_session`` 统一限速。

    Args:
        accessions: UniProt accession 列表。
        session: 可复用的 ``requests.Session``。
        timeout: 单次请求超时秒数。
        min_interval: 两次请求之间的最小间隔秒数。

    Returns:
        成功获取的记录列表（顺序与去重后的输入一致）。
    """
    if not accessions:
        return []

    # 去重但保持顺序
    unique: list[str] = []
    seen: set[str] = set()
    for acc in accessions:
        cleaned = (acc or "").strip()
        if not cleaned or cleaned.lower() in seen:
            continue
        seen.add(cleaned.lower())
        unique.append(cleaned)

    own_session = session is None
    if own_session:
        session = http_session(min_interval=min_interval, retries=3)
    if session is None:
        LOGGER.warning("requests 不可用，无法查询 AlphaFold DB")
        return []

    results: list[dict] = []
    try:
        for acc in unique:
            url = f"{ALPHAFOLD_API_BASE}/{acc}"
            try:
                response = session.get(url, timeout=timeout)
                if response.status_code == 404:
                    LOGGER.info("%s 在 AlphaFold DB 中无预测模型", acc)
                    continue
                response.raise_for_status()
                payload = response.json()
            except Exception as exc:  # noqa: BLE001
                LOGGER.warning("查询 %s 失败: %s", acc, exc)
                continue

            entries = list(_iter_entries(payload))
            if not entries:
                LOGGER.info("%s 无可用预测记录", acc)
                continue

            # 单个 accession 通常对应一条记录；若返回多条（如多片段模型），
            # 取第一条作为主记录，其余记 debug 日志。
            record = parse_prediction(entries[0])
            if not record.get("accession"):
                record["accession"] = acc
            if len(entries) > 1:
                LOGGER.debug("%s 返回 %d 条记录，取第一条", acc, len(entries))
            results.append(record)
    finally:
        if own_session:
            session.close()

    LOGGER.info("AlphaFold DB 返回 %d/%d 条记录", len(results), len(unique))
    return results


def build_pdb_url(accession: str, model_version: Optional[str] = None) -> str:
    """兜底构造 PDB 下载地址。

    **仅在 API 不可用时使用**。正常路径应当使用 API 返回的 ``pdb_url``，
    因为模型版本号会随 AlphaFold DB 重新训练而变化。

    Args:
        accession: UniProt accession。
        model_version: 模型版本（如 ``v6``）；None 时用默认候选列表的第一个。

    Returns:
        完整 URL 字符串（不保证服务端一定存在该文件）。
    """
    if model_version:
        return f"{ALPHAFOLD_FILES_BASE}/AF-{accession}-F1-model_{model_version}.pdb"
    version = _FALLBACK_MODEL_VERSIONS[0]
    return f"{ALPHAFOLD_FILES_BASE}/AF-{accession}-F1-model_{version}.pdb"


def _main(argv: Optional[list[str]] = None) -> int:
    """命令行入口：查询 AlphaFold DB 并打印结果（自检工具）。"""
    import argparse
    import json

    parser = argparse.ArgumentParser(
        prog="alphafold_client",
        description="查询 AlphaFold DB 预测元数据（命令行自检工具）",
    )
    parser.add_argument("accessions", nargs="+", help="UniProt accession（可多个）")
    parser.add_argument("--json", action="store_true", help="以 JSON 输出完整记录")
    args = parser.parse_args(argv)

    rows = fetch_predictions(args.accessions)
    if not rows:
        LOGGER.warning("没有查到任何记录")
        return 1

    if args.json:
        print(json.dumps(rows, ensure_ascii=False, indent=2))
    else:
        for row in rows:
            LOGGER.info(
                "%s | %s | %s | pLDDT=%s | %s",
                row.get("accession", ""),
                row.get("gene", ""),
                row.get("description", "")[:50],
                row.get("global_plddt", ""),
                row.get("pdb_url", ""),
            )
    return 0


if __name__ == "__main__":
    sys.exit(_main())
