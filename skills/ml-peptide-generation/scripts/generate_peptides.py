#!/usr/bin/env python3
"""用 PepMLM 生成候选肽链。

对应老仓库 ``design_pepmlm/scripts/batch_generate_peptides.py`` 与
``test_pepmlm.py``（二者功能高度重叠），并修掉它们的共同缺陷：

- 模型不再在 import 时加载——``--help`` 与参数校验无需 torch 即可执行。
- 新增 ``--seed``，同样的种子必然产出同样的肽链。
- 采样时排除 ``X``/``B``/``Z`` 等非法氨基酸，而不是事后才标记。
- 全部路径参数化，不再有 ``/Users/jung/Desktop/...``。
- 修掉老脚本「显示前 3 条结果」引用越界循环变量的 bug。

Examples:
    # 从 FASTA 批量生成
    python generate_peptides.py --input-fasta targets.fasta --output-dir results/

    # 单条序列
    python generate_peptides.py --input-seq MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMF \\
        --peptide-length 15 --num-binders 20 --seed 42

    # 扫参数组合（用于确定 top_k 与长度）
    python generate_peptides.py --input-fasta targets.fasta --output-dir results/ \\
        --peptide-lengths 10,12,15 --top-ks 3,5,8
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from collections.abc import Sequence
from datetime import datetime
from pathlib import Path
from typing import Optional

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "packages"))

from biopipe_core import ensure_dir, get_logger  # noqa: E402
from pepmlm_client import (  # noqa: E402
    PepMLMClient,
    PepMLMConfig,
    PepMLMUnavailableError,
    parse_target_records,
    set_seed,
)

LOGGER = get_logger("biopipe.ml-peptide-generation")

#: 输出 CSV 的列顺序
FIELDNAMES = [
    "rank",
    "target_id",
    "peptide_length",
    "top_k",
    "binder",
    "binder_length",
    "pseudo_perplexity",
    "has_unknown_aa",
    "target_length",
    "config",
    "seed",
    "generated_at",
]


def build_parser() -> argparse.ArgumentParser:
    """构造命令行参数解析器。"""
    parser = argparse.ArgumentParser(
        prog="generate_peptides",
        description="用 PepMLM 掩码语言模型生成候选肽链（可复现）",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "示例:\n"
            "  python generate_peptides.py --input-fasta targets.fasta --output-dir results/\n"
            "  python generate_peptides.py --input-seq MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMF \\\n"
            "      --peptide-length 15 --num-binders 20 --seed 42\n"
            "\n"
            "说明:\n"
            "  · 同样的 --seed 与同样的输入必然产出同样的结果。\n"
            "  · 采样时已排除 X/B/Z 等非标准氨基酸。\n"
            "  · 需要 torch 与 transformers；没有时可用 analyze_peptides.py 做离线分析。\n"
        ),
    )
    source = parser.add_argument_group("输入")
    source.add_argument("--input-fasta", default=None, help="靶标蛋白 FASTA 文件")
    source.add_argument("--input-csv", default=None, help="靶标蛋白 CSV 文件（需含 id 与 seq 列）")
    source.add_argument("--input-seq", default=None, help="直接给出单条蛋白序列")

    output = parser.add_argument_group("输出")
    output.add_argument(
        "--output-dir",
        default="peptide_results",
        help="输出目录（默认: peptide_results）",
    )
    output.add_argument(
        "--per-target-dir",
        action="store_true",
        help="为每个靶标单独建子目录存放 CSV（对应老脚本的按蛋白保存）",
    )

    generation = parser.add_argument_group("生成参数")
    generation.add_argument(
        "--peptide-lengths",
        default="10,12,15",
        help="要生成的肽链长度，逗号分隔（默认: 10,12,15）",
    )
    generation.add_argument(
        "--top-ks",
        default="3",
        help="要尝试的 top-k 值，逗号分隔（默认: 3）",
    )
    generation.add_argument(
        "--num-binders",
        type=int,
        default=10,
        help="每个（靶标 × 长度 × top-k）组合生成多少条肽链（默认: 10）",
    )
    generation.add_argument(
        "--seed",
        type=int,
        default=42,
        help="随机种子——相同种子必然复现相同结果（默认: 42）",
    )
    generation.add_argument(
        "--model",
        default=None,
        help="预训练模型名（默认: TianlaiChen/PepMLM-650M）",
    )
    generation.add_argument(
        "--device",
        default=None,
        help="强制指定设备（cuda / cpu）；默认自动选择",
    )
    generation.add_argument(
        "--allow-invalid-amino-acids",
        action="store_true",
        help="允许采样出 X/B/Z 等非标准氨基酸（默认禁止）",
    )
    generation.add_argument(
        "--timeout",
        type=int,
        default=120,
        help="单次请求超时秒数（下载模型时用；默认: 120）",
    )
    return parser


def parse_int_list(raw: str, label: str) -> list[int]:
    """把 ``"10,12,15"`` 解析成整数列表。"""
    values: list[int] = []
    for piece in raw.split(","):
        piece = piece.strip()
        if not piece:
            continue
        try:
            values.append(int(piece))
        except ValueError as exc:
            raise ValueError(f"{label} 含非法值: {piece!r}（应为整数）") from exc
    if not values:
        raise ValueError(f"{label} 不能为空")
    return values


def collect_targets(args: argparse.Namespace) -> list[dict]:
    """根据参数收集靶标序列。"""
    records: list[dict] = []
    if args.input_fasta:
        records.extend(parse_target_records(fasta_path=Path(args.input_fasta)))
    if args.input_csv:
        records.extend(parse_target_records(csv_path=Path(args.input_csv)))
    if args.input_seq:
        records.append({"id": "target", "seq": args.input_seq.strip().upper()})
    return records


def write_outputs(rows: Sequence[dict], output_dir: Path, per_target_dir: bool) -> None:
    """把结果写成 CSV 与 JSON。"""
    ensure_dir(output_dir)

    csv_path = output_dir / "peptides.csv"
    with open(csv_path, "w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    LOGGER.info("结果已写入: %s", csv_path)

    json_path = output_dir / "peptides.json"
    json_path.write_text(json.dumps(list(rows), ensure_ascii=False, indent=2), encoding="utf-8")
    LOGGER.info("结果已写入: %s", json_path)

    if per_target_dir:
        by_target: dict[str, list[dict]] = {}
        for row in rows:
            by_target.setdefault(row["target_id"], []).append(row)
        for target_id, target_rows in by_target.items():
            target_dir = ensure_dir(output_dir / target_id)
            target_path = target_dir / f"{target_id}_peptides.csv"
            with open(target_path, "w", encoding="utf-8-sig", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, extrasaction="ignore")
                writer.writeheader()
                writer.writerows(target_rows)
            LOGGER.info("靶标 %s 的结果已写入: %s", target_id, target_path)


def main(argv: Optional[list[str]] = None) -> int:
    """命令行主入口。"""
    args = build_parser().parse_args(argv)

    if not (args.input_fasta or args.input_csv or args.input_seq):
        LOGGER.error("必须提供 --input-fasta / --input-csv / --input-seq 之一")
        return 2

    try:
        peptide_lengths = parse_int_list(args.peptide_lengths, "--peptide-lengths")
        top_ks = parse_int_list(args.top_ks, "--top-ks")
    except ValueError as exc:
        LOGGER.error("%s", exc)
        return 2

    if args.num_binders <= 0:
        LOGGER.error("--num-binders 必须为正整数")
        return 2

    try:
        targets = collect_targets(args)
    except (OSError, ValueError) as exc:
        LOGGER.error("读取靶标失败: %s", exc)
        return 2

    if not targets:
        LOGGER.error("没有读到任何靶标序列")
        return 2

    LOGGER.info("共 %d 个靶标，肽长 %s，top-k %s，每组合 %d 条",
                len(targets), peptide_lengths, top_ks, args.num_binders)

    # 种子在**任何**随机操作之前设置——包括模型加载过程中的 dropout。
    set_seed(args.seed)

    config = PepMLMConfig(
        model_name=args.model or PepMLMConfig().model_name,
        peptide_lengths=tuple(peptide_lengths),
        top_k=top_ks[0],
        num_binders=args.num_binders,
        seed=args.seed,
        device=args.device,
        exclude_invalid_amino_acids=not args.allow_invalid_amino_acids,
    )

    try:
        client = PepMLMClient(config)
    except Exception as exc:  # noqa: BLE001
        LOGGER.error("初始化失败: %s", exc)
        return 1

    rows: list[dict] = []
    generated_at = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    for target_index, target in enumerate(targets, start=1):
        target_id = target["id"] or f"target_{target_index}"
        protein_seq = target["seq"]
        LOGGER.info("[%d/%d] %s（蛋白长度 %d）",
                    target_index, len(targets), target_id, len(protein_seq))

        for peptide_length in peptide_lengths:
            for top_k in top_ks:
                try:
                    peptides = client.generate_peptides(
                        protein_seq,
                        peptide_length=peptide_length,
                        top_k=top_k,
                        num_binders=args.num_binders,
                    )
                except ValueError as exc:
                    LOGGER.warning("跳过 %s / %daa / top-k %d: %s",
                                   target_id, peptide_length, top_k, exc)
                    continue
                except PepMLMUnavailableError as exc:
                    LOGGER.error("%s", exc)
                    return 1

                for peptide in peptides:
                    rows.append(
                        {
                            "rank": 0,
                            "target_id": target_id,
                            "peptide_length": peptide_length,
                            "top_k": top_k,
                            "binder": peptide["binder"],
                            "binder_length": peptide["binder_length"],
                            "pseudo_perplexity": peptide["pseudo_perplexity"],
                            "has_unknown_aa": peptide["has_unknown_aa"],
                            "target_length": len(protein_seq),
                            "config": f"{peptide_length}aa_topk{top_k}",
                            "seed": args.seed,
                            "generated_at": generated_at,
                        }
                    )

                LOGGER.info("  %daa / top-k %d: 生成 %d 条", peptide_length, top_k, len(peptides))

    if not rows:
        LOGGER.error("没有生成任何结果")
        return 1

    # 全局按伪困惑度升序排序，并编号。老脚本对每个靶标单独排序，
    # 这里额外给出全局名次，方便跨靶标比较。
    rows.sort(key=lambda row: float(row["pseudo_perplexity"]))
    for index, row in enumerate(rows, start=1):
        row["rank"] = index

    output_dir = Path(args.output_dir)
    write_outputs(rows, output_dir, args.per_target_dir)

    valid = sum(1 for row in rows if not row["has_unknown_aa"])
    LOGGER.info("共 %d 条结果（合法氨基酸 %d/%d），输出目录: %s",
                len(rows), valid, len(rows), output_dir.resolve())
    return 0


if __name__ == "__main__":
    sys.exit(main())
