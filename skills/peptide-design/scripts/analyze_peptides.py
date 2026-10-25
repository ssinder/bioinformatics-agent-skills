#!/usr/bin/env python3
"""多肽序列分析入口：理化性质计算、抗菌判别、筛选与导出。

功能链路：
    序列读取（FASTA / CSV）-> 理化性质计算（自包含公式实现）
    -> 抗菌肽判别打分 -> 高质量样本标记 -> FASTA / CSV 导出 -> 统计报告

理化性质计算完全自包含（``peptide_properties.py``），不依赖 Biopython：
分子量按残基质量表累加 + 水；等电点在 pH 0–14 的电荷曲线上二分求零点；
疏水性用 Kyte-Doolittle 标准值表；不稳定指数用 Guruprasad 公式。

用法示例::

    # 分析 FASTA 并导出带理化性质的 CSV
    python analyze_peptides.py --input peptides.fasta --output peptides_props.csv

    # 判别抗菌肽并导出筛选结果
    python analyze_peptides.py --input peptides.csv --sequence-column sequence \\
        --classify --antimicrobial-threshold 5 --output amp_selected.csv

    # 只输出 FASTA 与统计报告
    python analyze_peptides.py --input peptides.fasta --output out/ \\
        --export-fasta out/selected.fasta --report out/summary.txt
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "packages"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from biopipe_core import (  # noqa: E402
    SkillError,
    cli_parser,
    get_config,
    get_logger,
    read_csv_safe,
    read_fasta,
    resolve_path,
    run_cli,
    write_csv,
    write_fasta,
)
from peptide_filters import (  # noqa: E402
    annotate_go_distribution,
    classify,
    cys_rich,
    find_repeated_motifs,
    gc_content,
    quality_score,
    sequence_complexity,
    split_go_terms,
)
from peptide_properties import SequenceError, compute_properties  # noqa: E402

LOGGER = get_logger("biopipe.peptide-design")

# 输出的理化性质列
PROPERTY_FIELDS = [
    "length",
    "molecular_weight",
    "monoisotopic_weight",
    "isoelectric_point",
    "net_charge_ph7",
    "kyte_doolittle",
    "hydrophobic_moment",
    "hydrophobic_ratio",
    "aromaticity",
    "aliphatic_index",
    "instability_index",
    "predicted_solubility",
    "boman_index",
    "helix_fraction",
    "sheet_fraction",
    "turn_fraction",
    "helix_propensity",
    "gc_content",
    "sequence_complexity",
    "charged_ratio",
    "cysteine_count",
    "disulfide_bonds",
    "extinction_coefficient",
    "estimated_kd",
]


def build_parser() -> argparse.ArgumentParser:
    """构建命令行解析器。"""
    parser = cli_parser(
        prog="peptide-design",
        description="多肽理化性质计算、抗菌判别与序列导出",
    )
    parser.add_argument("--input", required=True, help="输入 FASTA 或 CSV 文件")
    parser.add_argument("--output", help="输出 CSV 路径（或目录，取决于 --output-dir）")
    parser.add_argument("--output-dir", help="输出目录（与 --export-fasta/--report 配合）")
    parser.add_argument("--format", default="auto", choices=["auto", "fasta", "csv"],
                        help="输入文件格式（auto 会按扩展名判断）")
    parser.add_argument("--id-column", default="id", help="CSV 中的序列 ID 列名")
    parser.add_argument("--sequence-column", default="sequence", help="CSV 中的序列列名")
    parser.add_argument(
        "--name-column",
        default="protein_name",
        help="CSV 中的名称列（用于 GO/关键词注释匹配，可缺省）",
    )
    parser.add_argument("--go-column", default="go_terms", help="CSV 中的 GO 术语列")
    parser.add_argument("--keyword-column", default="keywords", help="CSV 中的关键词列")
    parser.add_argument("--protein-existence-column", default="protein_existence",
                        help="CSV 中的蛋白质存在水平列")
    parser.add_argument("--classify", action="store_true", help="启用抗菌肽判别")
    parser.add_argument("--antimicrobial-threshold", type=float, default=5.0,
                        help="抗菌肽判定阈值（综合得分）")
    parser.add_argument("--non-antimicrobial-threshold", type=float, default=2.0,
                        help="非抗菌肽判定阈值（综合得分不高于此值）")
    parser.add_argument("--annotation-weight", type=float, default=0.6,
                        help="注释层得分在综合分中的权重（0–1）")
    parser.add_argument("--label", default="antimicrobial",
                        help="筛选时保留的类别: antimicrobial / non_antimicrobial / uncertain / all")
    parser.add_argument("--min-length", type=int, default=0, help="最短长度过滤")
    parser.add_argument("--max-length", type=int, default=0, help="最长长度过滤（0 表示不限）")
    parser.add_argument("--drop-unknown", action="store_true",
                        help="剔除含 X 等未知残基的序列")
    parser.add_argument("--require-complete-metadata", action="store_true",
                        help="只保留同时具备 GO 术语与关键词的记录")
    parser.add_argument("--min-quality", type=float, default=0.0,
                        help="高质量样本最低质量分")
    parser.add_argument("--export-fasta", help="额外导出 FASTA 的路径")
    parser.add_argument("--report", help="额外导出统计报告（文本）的路径")
    parser.add_argument("--json-output", help="额外导出完整 JSON 的路径")
    parser.add_argument("--id-prefix", default="pep", help="FASTA 记录 ID 前缀")
    return parser


def _detect_format(path: Path, declared: str) -> str:
    """判断输入格式。"""
    if declared != "auto":
        return declared
    suffix = path.suffix.lower()
    if suffix in (".fasta", ".fa", ".fna", ".fas", ".txt"):
        return "fasta"
    if suffix in (".csv", ".tsv"):
        return "csv"
    # 内容嗅探：以 '>' 开头视为 FASTA
    try:
        with open(path, encoding="utf-8", errors="replace") as handle:
            head = handle.read(4096).lstrip()
        return "fasta" if head.startswith(">") else "csv"
    except OSError:
        return "csv"


def load_sequences(args: argparse.Namespace, path: Path) -> list[dict[str, str]]:
    """载入序列记录。

    Returns:
        记录列表，每条含 ``id`` ``sequence`` 与原始列。
    """
    fmt = _detect_format(path, args.format)

    if fmt == "fasta":
        records: list[dict[str, str]] = []
        for index, (header, sequence) in enumerate(read_fasta(path), 1):
            records.append(
                {
                    "id": header or f"{args.id_prefix}{index}",
                    "sequence": sequence,
                    "name": header,
                    "go_terms": "",
                    "keywords": "",
                    "protein_existence": "",
                }
            )
        return records

    rows = read_csv_safe(path)
    records = []
    for index, row in enumerate(rows, 1):
        seq = row.get(args.sequence_column, "") or ""
        if not seq:
            # 退而求其次：找第一个看起来像序列的列
            for candidate in ("seq", "Sequence", "peptide", "aa_seq"):
                if row.get(candidate):
                    seq = row[candidate]
                    break
        records.append(
            {
                "id": row.get(args.id_column, "") or row.get("accession", "") or f"{args.id_prefix}{index}",
                "sequence": seq,
                "name": row.get(args.name_column, "") or row.get("protein_name", "") or "",
                "go_terms": row.get(args.go_column, "") or row.get("go_terms", "") or "",
                "keywords": row.get(args.keyword_column, "") or row.get("keywords", "") or "",
                "protein_existence": row.get(args.protein_existence_column, "") or "",
                "_raw": row,
            }
        )
    return records


def analyze_records(
    records: list[dict[str, str]],
    args: argparse.Namespace,
) -> list[dict]:
    """对每条序列计算性质并做判别。

    Returns:
        分析完成的记录列表。
    """
    results: list[dict] = []
    skipped = 0

    for record in records:
        sequence = str(record.get("sequence", "") or "").strip()
        if not sequence:
            skipped += 1
            continue
        try:
            props = compute_properties(sequence, allow_unknown=not args.drop_unknown)
        except SequenceError as exc:
            LOGGER.warning("跳过序列 %s: %s", record.get("id"), exc)
            skipped += 1
            continue

        row = dict(record.get("_raw", {}) or {})
        row["id"] = record["id"]
        row["sequence"] = sequence
        row["name"] = record.get("name", "")
        row["go_terms"] = record.get("go_terms", "")
        row["keywords"] = record.get("keywords", "")
        row["protein_existence"] = record.get("protein_existence", "")
        row["has_unknown"] = props.get("has_unknown", False)
        row["unknown_residues"] = props.get("unknown_residues", "")

        # GO 三大类别：已有分列直接沿用，否则按 UniProt 的 F:/P:/C: 前缀还原
        buckets = split_go_terms(row["go_terms"])
        row["go_molecular_function"] = row.get("go_molecular_function") or "; ".join(
            buckets["molecular_function"]
        )
        row["go_biological_process"] = row.get("go_biological_process") or "; ".join(
            buckets["biological_process"]
        )
        row["go_cellular_component"] = row.get("go_cellular_component") or "; ".join(
            buckets["cellular_component"]
        )

        for field in PROPERTY_FIELDS:
            row[field] = props.get(field, "")
        row["gc_content"] = gc_content(sequence)
        row["sequence_complexity"] = sequence_complexity(sequence)
        row["cys_rich"] = "Yes" if cys_rich(sequence) else "No"
        repeats = find_repeated_motifs(sequence)
        row["repeated_motifs"] = ";".join(repeats[:5])
        row["qc_flags"] = props.get("qc_flags", "")
        row["quality_score"] = quality_score(record, props)

        if args.classify:
            verdict = classify(
                record,
                props,
                antimicrobial_threshold=args.antimicrobial_threshold,
                non_antimicrobial_threshold=args.non_antimicrobial_threshold,
                annotation_weight=args.annotation_weight,
            )
            row["label"] = verdict["label"]
            row["combined_score"] = verdict["total"]
            row["annotation_score"] = verdict["annotation_score"]
            row["physicochemical_score"] = verdict["physicochemical_score"]
            row["hit_keywords"] = "; ".join(verdict["hit_keywords"])
            row["classification_reasons"] = " | ".join(verdict["reasons"])

        results.append(row)

    if skipped:
        LOGGER.warning("跳过 %d 条无法处理的记录", skipped)
    return results


def apply_filters(rows: list[dict], args: argparse.Namespace) -> list[dict]:
    """按长度、未知残基、元数据完整度、质量分与类别筛选。"""
    selected: list[dict] = []
    for row in rows:
        length = int(row.get("length", 0) or 0)
        if args.min_length and length < args.min_length:
            continue
        if args.max_length and length > args.max_length:
            continue
        if args.drop_unknown and row.get("has_unknown"):
            continue
        if args.require_complete_metadata:
            if not str(row.get("go_terms", "") or "").strip():
                continue
            if not str(row.get("keywords", "") or "").strip():
                continue
        if float(row.get("quality_score", 0.0) or 0.0) < args.min_quality:
            continue
        if args.classify and args.label != "all":
            if row.get("label") != args.label:
                continue
        selected.append(row)
    return selected


def build_report(rows: list[dict], all_rows: list[dict], args: argparse.Namespace) -> str:
    """生成文本统计报告。"""
    lines: list[str] = []
    lines.append("=" * 68)
    lines.append("多肽分析报告")
    lines.append("=" * 68)
    lines.append(f"输入记录总数: {len(all_rows)}")
    lines.append(f"通过筛选记录数: {len(rows)}")
    lines.append("")

    if not all_rows:
        return "\n".join(lines) + "\n"

    lengths = [int(r.get("length", 0) or 0) for r in all_rows]
    mws = [float(r.get("molecular_weight", 0) or 0) for r in all_rows]
    pis = [float(r.get("isoelectric_point", 0) or 0) for r in all_rows]
    kds = [float(r.get("kyte_doolittle", 0) or 0) for r in all_rows]
    instab = [float(r.get("instability_index", 0) or 0) for r in all_rows]

    def _stats(name: str, values: list[float], unit: str = "") -> None:
        values_sorted = sorted(values)
        n = len(values_sorted)
        median = values_sorted[n // 2] if n % 2 else (values_sorted[n // 2 - 1] + values_sorted[n // 2]) / 2
        mean = sum(values) / n
        lines.append(
            f"{name}: 最小 {min(values):.2f}{unit} | 最大 {max(values):.2f}{unit} | "
            f"均值 {mean:.2f}{unit} | 中位数 {median:.2f}{unit}"
        )

    lines.append("--- 理化性质统计 ---")
    _stats("长度", [float(x) for x in lengths], " aa")
    _stats("分子量", mws, " Da")
    _stats("等电点 pI", pis)
    _stats("Kyte-Doolittle 疏水性", kds)
    _stats("不稳定指数", instab)
    lines.append("")

    lines.append("--- 稳定性与组成 ---")
    stable = sum(1 for x in instab if x < 40)
    lines.append(f"不稳定指数 < 40 (预测稳定): {stable} ({stable / len(all_rows) * 100:.1f}%)")
    positive = sum(1 for r in all_rows if float(r.get("net_charge_ph7", 0) or 0) > 0)
    lines.append(f"净正电荷 (pH 7): {positive} ({positive / len(all_rows) * 100:.1f}%)")
    amphipathic = sum(
        1 for r in all_rows
        if float(r.get("hydrophobic_moment", 0) or 0) >= 0.4
        and abs(float(r.get("net_charge_ph7", 0) or 0)) >= 1.0
    )
    lines.append(f"疑似两性 (疏水矩>=0.4 且带电): {amphipathic} ({amphipathic / len(all_rows) * 100:.1f}%)")
    unknown = sum(1 for r in all_rows if r.get("has_unknown"))
    lines.append(f"含未知残基: {unknown} ({unknown / len(all_rows) * 100:.1f}%)")
    lines.append("")

    if args.classify:
        lines.append("--- 抗菌肽判别 ---")
        counts: dict[str, int] = {}
        for row in all_rows:
            label = str(row.get("label", "uncertain"))
            counts[label] = counts.get(label, 0) + 1
        for label in ("antimicrobial", "uncertain", "non_antimicrobial"):
            n = counts.get(label, 0)
            lines.append(f"{label}: {n} ({n / len(all_rows) * 100:.1f}%)")
        lines.append(f"判定阈值: 抗菌 >= {args.antimicrobial_threshold}, "
                     f"非抗菌 <= {args.non_antimicrobial_threshold}")
        lines.append(f"注释层权重: {args.annotation_weight}")
        lines.append("")

    go_stats = annotate_go_distribution([r for r in all_rows if "go_terms" in r])
    if go_stats["total"]:
        lines.append("--- GO 注释覆盖度 ---")
        lines.append(f"分子功能: {go_stats['molecular_function']} ({go_stats['molecular_function_pct']}%)")
        lines.append(f"生物过程: {go_stats['biological_process']} ({go_stats['biological_process_pct']}%)")
        lines.append(f"细胞组分: {go_stats['cellular_component']} ({go_stats['cellular_component_pct']}%)")
        lines.append("")

    organisms: dict[str, int] = {}
    for row in all_rows:
        org = str(row.get("organism", "") or "").strip()
        if org:
            organisms[org] = organisms.get(org, 0) + 1
    if organisms:
        lines.append("--- 物种分布 (Top 10) ---")
        for org, n in sorted(organisms.items(), key=lambda kv: kv[1], reverse=True)[:10]:
            lines.append(f"{org}: {n} ({n / len(all_rows) * 100:.1f}%)")
        lines.append("")

    return "\n".join(lines) + "\n"


def main(args: argparse.Namespace) -> int:
    """主流程。"""
    config = get_config(args.config)
    input_path = resolve_path(args.input)
    if not input_path.exists():
        raise SkillError(f"输入文件不存在: {input_path}")

    LOGGER.info("读取输入: %s", input_path)
    records = load_sequences(args, input_path)
    if not records:
        LOGGER.warning("输入文件没有可分析的记录")
        return 0
    LOGGER.info("载入 %d 条序列", len(records))

    all_rows = analyze_records(records, args)
    if not all_rows:
        LOGGER.warning("没有成功分析的序列")
        return 0

    selected = apply_filters(all_rows, args)
    LOGGER.info("筛选后保留 %d 条", len(selected))

    # ---- 输出 CSV ------------------------------------------------------
    if args.output:
        output_path = resolve_path(args.output)
    elif args.output_dir:
        output_path = resolve_path(args.output_dir) / "peptide_analysis.csv"
    else:
        output_path = config.output_dir / "peptide_analysis.csv"

    if selected:
        fieldnames = list(selected[0].keys())
        write_csv(selected, output_path, fieldnames=fieldnames)
        LOGGER.info("已写出 CSV: %s (%d 行)", output_path, len(selected))
    else:
        LOGGER.warning("筛选后无记录，未写出 CSV")

    # ---- 输出 FASTA ----------------------------------------------------
    if args.export_fasta:
        fasta_path = resolve_path(args.export_fasta)
        write_fasta(
            [(f"{args.id_prefix}_{row['id']}", row["sequence"]) for row in selected],
            fasta_path,
        )
        LOGGER.info("已写出 FASTA: %s (%d 条)", fasta_path, len(selected))

    # ---- 输出报告 ------------------------------------------------------
    if args.report:
        report_path = resolve_path(args.report)
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text(build_report(selected, all_rows, args), encoding="utf-8")
        LOGGER.info("已写出报告: %s", report_path)

    # ---- 输出 JSON ------------------------------------------------------
    if args.json_output:
        json_path = resolve_path(args.json_output)
        json_path.parent.mkdir(parents=True, exist_ok=True)
        with open(json_path, "w", encoding="utf-8") as handle:
            json.dump(selected, handle, ensure_ascii=False, indent=2)
        LOGGER.info("已写出 JSON: %s", json_path)

    return 0


if __name__ == "__main__":
    sys.exit(
        run_cli(
            main,
            parser=build_parser(),
            skill="peptide-design",
        )
    )
