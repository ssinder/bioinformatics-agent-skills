#!/usr/bin/env python3
"""批量下载 AlphaFold DB 预测结构。

核心改进（相对老仓库 ``afdb_download/``）：

1. **不再写死模型版本号**。老脚本拼 ``AF-{ID}-F1-model_v6.pdb``，
   AlphaFold DB 每次重训模型版本号就会递增，写死必然在下一版 404。
   本脚本改为调 API 拿 ``pdbUrl``。
2. **ID 列表不再丢首个**。老脚本无条件跳过第一行（假设是表头），
   用户若提供纯 ID 文件，第一个 ID 会被静默丢弃。
3. **二进制写文件**。老脚本用 ``response.text`` + 文本模式写入，
   会按平台换行转换破坏 PDB 文件内容。
4. **校验 Content-Type**。拿到 HTML 错误页也当 PDB 存下来，
   会让下游解析器报出莫名其妙的错。
5. **所有输出路径来自命令行参数**，不再有写死的相对目录。

Examples:
    # 从 ID 文件批量下载
    python download_structures.py --ids ids.txt --output-dir structures/

    # 直接给 ID
    python download_structures.py --ids P69905 P04637 --output-dir structures/

    # 只查元数据不下载
    python download_structures.py --ids P69905 --metadata-only --output-dir meta/
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Any, Optional

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "packages"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from alphafold_client import fetch_predictions  # noqa: E402
from biopipe_core import ensure_dir, get_logger  # noqa: E402

LOGGER = get_logger("biopipe.alphafold-db")

# 允许的 PDB 内容类型。AlphaFold DB 返回 text/plain；
# 某些代理/CDN 可能返回 octet-stream。
_ALLOWED_CONTENT_TYPES = ("text/plain", "application/octet-stream", "chemical/x-pdb")


def read_id_file(path: Path, skip_header: bool = True) -> list[str]:
    """从文本文件读取 UniProt accession 列表。

    Args:
        path: 文件路径，每行一个 ID（``#`` 开头为注释）。
        skip_header: 是否跳过首个非注释行。老仓库的做法是**无条件**跳过
            第一行，这会在用户提供纯 ID 文件时静默丢掉一个 ID。
            这里改为自动判断——首个非注释行若能匹配 accession 格式就保留。

    Returns:
        去重后的 accession 列表（保持原顺序）。
    """
    if not path.is_file():
        raise FileNotFoundError(f"ID 列表文件不存在: {path}")

    raw_lines: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        raw_lines.append(stripped)

    if not raw_lines:
        return []

    if skip_header and not _looks_like_accession(raw_lines[0]):
        # 首行不像 accession，视为表头丢弃
        LOGGER.info("跳过表头行: %s", raw_lines[0])
        raw_lines = raw_lines[1:]
    elif skip_header:
        LOGGER.debug("首行 %s 看起来就是 accession，不作为表头丢弃", raw_lines[0])

    result: list[str] = []
    seen: set[str] = set()
    for item in raw_lines:
        # 兼容 "P69905, description" 这种带逗号注释的写法
        candidate = item.split(",")[0].strip()
        if not candidate or candidate.lower() in seen:
            continue
        seen.add(candidate.lower())
        result.append(candidate)
    return result


def _looks_like_accession(text: str) -> bool:
    """判断字符串是否像 UniProt accession。

    UniProt accession 的典型形态：1 个字母 + 5 位数字（如 P69905），
    或 3~10 字母 + 数字的组合（如 A0A024RBG1、Q9Y6K9）。
    """
    candidate = text.split(",")[0].strip()
    if not candidate or len(candidate) < 6 or len(candidate) > 12:
        return False
    if not candidate[0].isalpha() or not candidate[0].isupper():
        return False
    return any(ch.isdigit() for ch in candidate)


def download_structure(
    url: str,
    destination: Path,
    session=None,
    timeout: int = 120,
) -> tuple[bool, str]:
    """下载单个结构文件。

    Args:
        url: 结构文件的完整 URL（来自 API 的 ``pdbUrl``）。
        destination: 本地目标路径。
        session: 可复用的 ``requests.Session``。
        timeout: 单次请求超时秒数。

    Returns:
        ``(是否成功, 说明信息)``。
    """
    from biopipe_core import http_session

    own_session = session is None
    if own_session:
        session = http_session(min_interval=0.2, retries=3)
    if session is None:
        return False, "requests 不可用"

    try:
        response = session.get(url, timeout=timeout)
        response.raise_for_status()

        # 校验内容类型——防止把 HTML 错误页当成 PDB 存下来
        content_type = (response.headers.get("Content-Type") or "").split(";")[0].strip()
        if content_type and not content_type.startswith(_ALLOWED_CONTENT_TYPES):
            return False, f"非预期内容类型 {content_type}（可能不是 PDB 文件）"

        # 二进制写文件，避免文本模式的换行转换破坏 PDB 内容
        destination.parent.mkdir(parents=True, exist_ok=True)
        payload = response.content
        if not payload:
            return False, "响应体为空"
        destination.write_bytes(payload)

        return True, f"{len(payload)} 字节"
    except Exception as exc:  # noqa: BLE001
        return False, str(exc)
    finally:
        if own_session:
            session.close()


def write_metadata(records: Sequence[dict], destination: Path) -> None:
    """把预测元数据写成 CSV 和 JSON 两种格式。

    元数据里包含整体 pLDDT 与各置信区间残基占比——这些信息**免费**，
    但在下载 PDB 之前就能用来筛掉低质量模型，省下大量下载时间。
    """
    destination.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "accession",
        "entry_id",
        "gene",
        "description",
        "organism",
        "tax_id",
        "chain_id",
        "model_version",
        "model_created_date",
        "global_plddt",
        "fraction_plddt_very_low",
        "fraction_plddt_low",
        "fraction_plddt_confident",
        "fraction_plddt_very_high",
        "sequence_length",
        "is_reviewed",
        "is_complex",
        "pdb_url",
        "cif_url",
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
        prog="download_structures",
        description="从 AlphaFold DB 批量下载预测结构（PDB）并导出预测元数据",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "示例:\n"
            "  python download_structures.py --ids ids.txt --output-dir structures/\n"
            "  python download_structures.py --ids P69905 P04637 --output-dir structures/\n"
            "  python download_structures.py --ids ids.txt --metadata-only --output-dir meta/\n"
            "\n"
            "说明:\n"
            "  模型版本号不写死——下载地址来自 AlphaFold DB API 的 pdbUrl 字段，\n"
            "  因此模型重新训练后无需修改脚本。\n"
        ),
    )
    parser.add_argument(
        "--ids",
        nargs="+",
        required=True,
        help="UniProt accession 列表，或包含 ID 的文本文件路径（每行一个）",
    )
    parser.add_argument(
        "--output-dir",
        default="alphafold_structures",
        help="结构文件输出目录（默认: alphafold_structures）",
    )
    parser.add_argument(
        "--metadata-only",
        action="store_true",
        help="只查询并导出元数据，不下载 PDB 文件",
    )
    parser.add_argument(
        "--min-plddt",
        type=float,
        default=None,
        help="整体 pLDDT 低于该值的模型跳过下载（0-100）。缺省则全部下载",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="已存在的文件也重新下载",
    )
    parser.add_argument(
        "--keep-original-name",
        action="store_true",
        help="保留 AlphaFold 原始文件名（AF-{ID}-F1-model_v6.pdb）而非简化为 {ID}.pdb",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=120,
        help="单次请求超时秒数（默认: 120）",
    )
    return parser


def _collect_accessions(items: Sequence[str]) -> list[str]:
    """把 --ids 参数解析成 accession 列表（支持文件路径与直接给 ID）。"""
    accessions: list[str] = []
    for item in items:
        candidate = Path(item)
        try:
            is_file = candidate.is_file()
        except OSError:  # pragma: no cover - 路径过长等极端情况
            is_file = False
        if is_file:
            LOGGER.info("从文件读取 ID: %s", candidate)
            accessions.extend(read_id_file(candidate))
        else:
            accessions.append(item.strip())
    return accessions


def _target_name(accession: str, keep_original: bool, pdb_url: str) -> str:
    """决定本地文件名。

    默认简化为 ``{accession}.pdb``；``--keep-original-name`` 时保留
    AlphaFold 的原始文件名（从 URL 里取，这样版本号也跟着正确）。
    """
    if keep_original:
        name = Path(pdb_url).name
        if name:
            return name
    return f"{accession}.pdb"


def main(argv: Optional[list[str]] = None) -> int:
    """命令行主入口。"""
    args = build_parser().parse_args(argv)

    accessions = _collect_accessions(args.ids)
    if not accessions:
        LOGGER.error("没有可用的 UniProt accession")
        return 2

    LOGGER.info("开始查询 %d 个 accession 的 AlphaFold 预测", len(accessions))

    records = fetch_predictions(accessions, timeout=args.timeout)
    if not records:
        LOGGER.error("没有查到任何 AlphaFold 预测记录")
        return 1

    output_dir = ensure_dir(args.output_dir)
    write_metadata(records, output_dir / "alphafold_metadata.csv")

    if args.metadata_only:
        LOGGER.info("--metadata-only 已指定，跳过文件下载")
        LOGGER.info("共 %d 条记录，元数据在 %s", len(records), output_dir)
        return 0

    from biopipe_core import http_session

    session = http_session(min_interval=0.2, retries=3)
    if session is None:
        LOGGER.error("requests 不可用，无法下载结构文件")
        return 1

    stats = {"downloaded": 0, "skipped": 0, "failed": 0, "filtered": 0}
    failed_ids: list[str] = []

    try:
        for index, record in enumerate(records, 1):
            accession = record.get("accession", "")
            LOGGER.info("[%d/%d] %s", index, len(records), accession)

            pdb_url = record.get("pdb_url", "")
            if not pdb_url:
                LOGGER.warning("%s 没有可用的 PDB 下载地址，跳过", accession)
                stats["failed"] += 1
                failed_ids.append(accession)
                continue

            if args.min_plddt is not None:
                plddt = record.get("global_plddt")
                try:
                    if plddt is not None and float(plddt) < args.min_plddt:
                        LOGGER.info(
                            "%s 整体 pLDDT=%s 低于阈值 %.1f，跳过下载",
                            accession,
                            plddt,
                            args.min_plddt,
                        )
                        stats["filtered"] += 1
                        continue
                except (TypeError, ValueError):
                    LOGGER.warning("%s 的 pLDDT 值 %r 无法解析", accession, plddt)

            filename = _target_name(accession, args.keep_original_name, pdb_url)
            destination = output_dir / filename

            if destination.exists() and not args.force:
                LOGGER.info("已存在，跳过: %s", filename)
                stats["skipped"] += 1
                continue

            ok, detail = download_structure(
                pdb_url, destination, session=session, timeout=args.timeout
            )
            if ok:
                LOGGER.info("下载成功: %s (%s)", filename, detail)
                stats["downloaded"] += 1
            else:
                LOGGER.error("下载失败: %s —— %s", accession, detail)
                stats["failed"] += 1
                failed_ids.append(accession)
    finally:
        session.close()

    LOGGER.info("-" * 56)
    LOGGER.info(
        "下载完成：新下载 %d，跳过 %d，阈值过滤 %d，失败 %d",
        stats["downloaded"],
        stats["skipped"],
        stats["filtered"],
        stats["failed"],
    )
    LOGGER.info("输出目录: %s", output_dir)

    if failed_ids:
        LOGGER.warning("失败的 accession: %s", ", ".join(failed_ids))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
