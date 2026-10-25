"""FASTA / CSV 读写工具。

统一处理：

- FASTA 的多行序列、``;`` 注释行、CRLF 换行
- CSV 的 UTF-8 BOM（Excel 打开中文乱码的常见原因）
- 大文件的分块读写
"""

from __future__ import annotations

import csv
from collections.abc import Iterable, Iterator, Sequence
from pathlib import Path
from typing import Any

__all__ = [
    "read_fasta",
    "write_fasta",
    "read_csv_safe",
    "write_csv",
    "ensure_unicode_csv",
]


def read_fasta(path: str | Path) -> Iterator[tuple[str, str]]:
    """读取 FASTA 文件，逐条产出 ``(header, sequence)``。

    Args:
        path: FASTA 文件路径。

    Yields:
        ``(标题行去掉 '>' 与首段描述后的标识, 序列字符串)``。
    """
    header: str | None = None
    chunks: list[str] = []

    with open(path, encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(";"):
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:].split()[0] if line[1:].split() else ""
                chunks = []
            elif header is not None:
                chunks.append(line.replace(" ", ""))

    if header is not None:
        yield header, "".join(chunks)


def write_fasta(records: Iterable[tuple[str, str]], path: str | Path, width: int = 60) -> Path:
    """写出 FASTA 文件。

    Args:
        records: ``(header, sequence)`` 序列。
        path: 输出路径。
        width: 每行序列长度，0 表示不换行。

    Returns:
        输出文件路径。
    """
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", encoding="utf-8", newline="\n") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n")
            if width and width > 0:
                for i in range(0, len(sequence), width):
                    handle.write(sequence[i : i + width] + "\n")
            else:
                handle.write(sequence + "\n")
    return out


def read_csv_safe(path: str | Path, delimiter: str | None = None, encoding: str | None = None) -> list[dict[str, str]]:
    """读取 CSV 为字典列表，自动处理 BOM 与编码回退。

    Args:
        path: CSV 路径。
        delimiter: 分隔符，None 时自动嗅探。
        encoding: 编码，None 时按 utf-8-sig → gbk → latin-1 依次尝试。

    Returns:
        行字典列表。
    """
    p = Path(path)
    encodings = [encoding] if encoding else ["utf-8-sig", "utf-8", "gbk", "latin-1"]
    last_error: Exception | None = None

    for enc in encodings:
        try:
            with open(p, encoding=enc, newline="") as handle:
                sample = handle.read(4096)
                handle.seek(0)
                if delimiter is None:
                    try:
                        delimiter = csv.Sniffer().sniff(sample, delimiters=",;\t|").delimiter
                    except csv.Error:
                        delimiter = ","
                reader = csv.DictReader(handle, delimiter=delimiter)
                return [dict(row) for row in reader]
        except UnicodeDecodeError as exc:
            last_error = exc
            continue

    raise UnicodeDecodeError(  # pragma: no cover
        "unknown", b"", 0, 1, f"无法用候选编码读取 {p}: {last_error}"
    )


def write_csv(
    rows: Sequence[dict[str, Any]],
    path: str | Path,
    fieldnames: Sequence[str] | None = None,
    delimiter: str = ",",
) -> Path:
    """写出 CSV 文件（UTF-8 带 BOM，Excel 直接打开不乱码）。

    Args:
        rows: 行字典列表。
        path: 输出路径。
        fieldnames: 列名，默认取首行的键。
        delimiter: 分隔符。

    Returns:
        输出文件路径。
    """
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    if not fieldnames and rows:
        fieldnames = list(rows[0].keys())

    with open(out, "w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(fieldnames or []), delimiter=delimiter, extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)
    return out


def ensure_unicode_csv(path: str | Path) -> Path:
    """把无 BOM 的 UTF-8 CSV 转换为带 BOM 版本，便于 Excel 正确识别中文。"""
    import tempfile

    p = Path(path)
    rows = read_csv_safe(p)
    fieldnames = list(rows[0].keys()) if rows else []
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8-sig", newline="", delete=False, suffix=".csv"
    ) as tmp:
        tmp_path = Path(tmp.name)
    write_csv(rows, tmp_path, fieldnames=fieldnames)
    tmp_path.replace(p)
    return p
