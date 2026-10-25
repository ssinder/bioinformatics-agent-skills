#!/usr/bin/env python3
"""AlphaFold 模型选择与置信度处理。

从 AlphaFold DB 下载的 zip 包中，每个条目包含多个模型（``*_model_1.cif`` …
``*_model_5.cif``）以及对应的置信度 JSON（``*_summary_confidences_0.json`` …
``*_summary_confidences_4.json``）。每个置信度文件含一个模型的逐残基 pLDDT /
PAE 数据，模型级别的汇总分数写在 ``summary_confidences`` 数组中。

本模块提供：

- :func:`select_best_model` —— 按置信度排序选出最佳模型
- :func:`scan_directory` —— 在目录中扫描 (模型文件, 置信度文件) 配对
- :func:`summarize_confidence` —— 提取 pTM / ipTM / 平均 pLDDT 等汇总指标

Usage:
    from alphafold_models import select_best_model
    best = select_best_model(zip_path)
"""

from __future__ import annotations

import json
import re
import zipfile
from pathlib import Path
from typing import Any, Optional, Union

__all__ = [
    "select_best_model",
    "scan_directory",
    "summarize_confidence",
    "load_confidence_json",
]

_MODEL_RE = re.compile(r"_model_(\d+)", re.IGNORECASE)
_CONF_RE = re.compile(r"_summary_confidences_(\d+)", re.IGNORECASE)


def load_confidence_json(source) -> dict:
    """读取置信度 JSON（支持文件路径、zip 内路径或已解析的字典）。

    Args:
        source: ``Path``、zip 内成员名，或已解析的 dict。

    Returns:
        置信度字典；解析失败返回空字典。
    """
    if isinstance(source, dict):
        return source
    try:
        if isinstance(source, str) and "::" in source:
            zip_path, member = source.split("::", 1)
            with zipfile.ZipFile(zip_path) as zf:
                return json.loads(zf.read(member).decode("utf-8"))
        return json.loads(Path(source).read_text(encoding="utf-8"))
    except (OSError, ValueError, KeyError):
        return {}


def summarize_confidence(confidence: dict) -> dict:
    """从置信度 JSON 中提取模型级汇总指标。

    兼容 AlphaFold DB 的两种 JSON 布局：

    - 单模型布局：``{"residue1": [...], "plddt": [...], "max_pae": ...,
      "ptm": 0.8, "iptm": 0.6}``
    - 集合布局：``{"summary_confidences": [ {...}, {...} ], "pae": [...]}``

    Returns:
        含 ``ptm`` / ``iptm`` / ``mean_plddt`` / ``max_pae`` / ``n_residues`` /
        ``high_confidence_fraction`` 的字典，缺失项为 None。
    """
    if not confidence:
        return {
            "ptm": None,
            "iptm": None,
            "mean_plddt": None,
            "max_pae": None,
            "n_residues": 0,
            "high_confidence_fraction": None,
        }

    plddt = confidence.get("plddt") or confidence.get("confidence") or []
    if not isinstance(plddt, list):
        plddt = []

    mean_plddt = sum(plddt) / len(plddt) if plddt else None
    high_fraction = (
        sum(1 for v in plddt if v >= 70) / len(plddt) if plddt else None
    )

    return {
        "ptm": confidence.get("ptm"),
        "iptm": confidence.get("iptm"),
        "mean_plddt": mean_plddt,
        "max_pae": confidence.get("max_pae"),
        "n_residues": len(plddt),
        "high_confidence_fraction": high_fraction,
    }


def _score_of(confidence: dict) -> float:
    """模型排序分数。

    优先使用 ``ptm``；缺失时退回平均 pLDDT / 100；
    再缺失则退回 0。复合物场景下可用 ``ptm + iptm`` 加权，
    这里采用保守的 ``ptm + 0.5 * iptm``（若两者都有）。
    """
    summary = summarize_confidence(confidence)
    ptm = summary.get("ptm")
    iptm = summary.get("iptm")
    if ptm is not None and iptm is not None:
        return float(ptm) + 0.5 * float(iptm)
    if ptm is not None:
        return float(ptm)
    mean_plddt = summary.get("mean_plddt")
    if mean_plddt is not None:
        return float(mean_plddt) / 100.0
    return 0.0


def select_best_model(
    zip_path: Union[str, Path],
    member_prefix: Optional[str] = None,
) -> Optional[dict]:
    """在 AlphaFold DB 的 zip 包中选出置信度最佳的模型。

    Args:
        zip_path: zip 文件路径。
        member_prefix: 只在该前缀（通常是条目名，如 ``AF-P69905-F1``）下查找。

    Returns:
        含 ``model_member`` / ``confidence_member`` / ``model_index`` /
        ``score`` / ``confidence``（原始 JSON）/ ``summary``（汇总指标）的
        字典；找不到模型时返回 None。
    """
    zpath = Path(zip_path)
    if not zpath.exists():
        raise FileNotFoundError(f"zip 不存在: {zpath}")

    models: dict[int, str] = {}
    confidences: dict[int, str] = {}
    folders: dict[str, dict[str, Any]] = {}

    with zipfile.ZipFile(zpath) as zf:
        for name in zf.namelist():
            if name.endswith("/"):
                continue
            folder = name.split("/", 1)[0]
            base = name.rsplit("/", 1)[-1]

            model_match = _MODEL_RE.search(base)
            conf_match = _CONF_RE.search(base)

            if base.endswith((".cif", ".pdb")) and model_match:
                index = int(model_match.group(1))
                models.setdefault(index, name)
                folders.setdefault(folder, {"models": {}, "confidences": {}})
                folders[folder]["models"][index] = name

            if base.endswith(".json") and conf_match:
                index = int(conf_match.group(1))
                confidences.setdefault(index, name)
                folders.setdefault(folder, {"models": {}, "confidences": {}})
                folders[folder]["confidences"][index] = name

        results: list[dict] = []
        for folder, content in folders.items():
            if member_prefix and folder != member_prefix:
                continue
            if not content["models"] or not content["confidences"]:
                continue

            best_index = None
            best_score = -1.0
            for index, member in sorted(content["confidences"].items()):
                try:
                    data = json.loads(zf.read(member).decode("utf-8"))
                except (ValueError, OSError):
                    continue
                score = _score_of(data)
                if score > best_score:
                    best_score = score
                    best_index = index

            if best_index is None:
                # 没有置信度文件时退回第一个模型
                best_index = min(content["models"]) if content["models"] else None
                best_score = 0.0

            if best_index is None or best_index not in content["models"]:
                continue

            conf_member = content["confidences"].get(best_index)
            conf_data = (
                load_confidence_json(f"{zpath}::{conf_member}") if conf_member else {}
            )
            results.append(
                {
                    "folder": folder,
                    "model_member": content["models"][best_index],
                    "confidence_member": conf_member,
                    "model_index": best_index,
                    "score": best_score,
                    "confidence": conf_data,
                    "summary": summarize_confidence(conf_data),
                }
            )

    if not results:
        return None
    # 多个条目时取分数最高者
    return max(results, key=lambda r: r["score"])


def scan_directory(directory: Union[str, Path], pattern: str = "*.zip") -> list[dict]:
    """扫描目录中所有 AlphaFold zip 包并选出每个包的最佳模型。

    Args:
        directory: 目录路径。
        pattern: zip 文件匹配模式。

    Returns:
        每个 zip 一个结果字典，含 ``zip_path`` 与最佳模型信息。
    """
    d = Path(directory)
    if not d.is_dir():
        raise FileNotFoundError(f"目录不存在: {d}")

    out: list[dict] = []
    for zip_path in sorted(d.glob(pattern)):
        try:
            best = select_best_model(zip_path)
        except (FileNotFoundError, zipfile.BadZipFile):
            continue
        if best is None:
            continue
        out.append({"zip_path": str(zip_path), **best})
    return out


def _main(argv: Optional[list[str]] = None) -> int:
    """命令行入口：列出 zip 中各条目的最佳模型。"""
    import argparse

    parser = argparse.ArgumentParser(
        prog="alphafold_models", description="AlphaFold 最佳模型选择（命令行自检工具）"
    )
    parser.add_argument("input", help="AlphaFold zip 文件或包含 zip 的目录")
    args = parser.parse_args(argv)

    target = Path(args.input)
    if target.is_dir():
        results = scan_directory(target)
    else:
        best = select_best_model(target)
        results = [{"zip_path": str(target), **best}] if best else []

    if not results:
        logger.warning("没有找到可用的模型")
        return 1
    for item in results:
        summary = item["summary"]
        logger.info(
            "%s | model_%s | 分数=%.3f | pTM=%s ipTM=%s 平均 pLDDT=%s",
            item.get("folder", "?"),
            item.get("model_index"),
            item.get("score", 0.0),
            summary.get("ptm"),
            summary.get("iptm"),
            (
                round(summary["mean_plddt"], 2)
                if summary.get("mean_plddt") is not None
                else None
            ),
        )
    return 0


if __name__ == "__main__":
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "packages"))
    from biopipe_core import get_logger

    logger = get_logger("biopipe.protein-interactions")
    sys.exit(_main())
