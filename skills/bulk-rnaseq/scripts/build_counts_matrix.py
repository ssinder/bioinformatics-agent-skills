#!/usr/bin/env python3
"""从 Salmon 定量结果构建基因水平 counts / TPM / FPKM 矩阵。

移植自原仓库 ``rnaseq_python/scripts/04_build_counts_matrix.py``，去掉了硬编码路径、
裸 print，并补齐 argparse 与日志。

核心算法：
- 读取每个样本子目录下的 ``quant.sf``（Salmon 转录本定量输出）
- 使用 ``NumReads`` 作为 counts（已经过长度偏差校正）
- 通过 ``tx2gene.tsv`` 将转录本水平汇总到基因水平（counts/TPM/FPKM 分别求和）
- counts 四舍五入为整数并剔除全零基因（下游差异表达分析要求整数 counts）
- 由样本名推断 condition / replicate，生成 ``metadata_template.csv``

Usage:
    python build_counts_matrix.py \\
        --quant-dir /path/to/salmon \\
        --output /path/to/counts \\
        --tx2gene /path/to/tx2gene.tsv
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "packages"))

import pandas as pd  # noqa: E402
from biopipe_core import SkillError, get_config, get_logger  # noqa: E402

LOGGER = get_logger("biopipe.bulk-rnaseq")


def build_counts_matrix(
    quant_dir: Path,
    output_dir: Path,
    tx2gene_file: Path | None,
    min_counts: int = 0,
) -> dict[str, Path]:
    """从 Salmon 定量目录构建基因水平矩阵。

    Args:
        quant_dir: Salmon 定量结果目录，每个样本一个子目录且含 ``quant.sf``。
        output_dir: 输出目录。
        tx2gene_file: 转录本到基因的 TSV 映射（列名 ``transcript_id`` / ``gene_id``）。
        min_counts: 基因在所有样本中的总 counts 下限，低于则剔除。

    Returns:
        产出文件路径字典。
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    tx2gene: dict[str, str] = {}
    if tx2gene_file is not None:
        tx2gene_file = Path(tx2gene_file)
        if not tx2gene_file.exists():
            raise SkillError(f"tx2gene 文件不存在: {tx2gene_file}")
        tx_df = pd.read_csv(tx2gene_file, sep="\t")
        for col in ("transcript_id", "gene_id"):
            if col not in tx_df.columns:
                raise SkillError(f"tx2gene 文件缺少必需列 {col}: {tx2gene_file}")
        tx2gene = dict(zip(tx_df["transcript_id"].astype(str), tx_df["gene_id"].astype(str)))
        LOGGER.info("加载 tx2gene 映射: %d 条", len(tx2gene))
    else:
        LOGGER.warning("未提供 tx2gene，转录本 ID 将直接作为基因 ID")

    quant_dir = Path(quant_dir)
    if not quant_dir.is_dir():
        raise SkillError(f"Salmon 定量目录不存在: {quant_dir}")

    sample_dirs = sorted(d for d in quant_dir.iterdir() if d.is_dir())
    if not sample_dirs:
        raise SkillError(f"在 {quant_dir} 中未找到样本子目录")

    counts_dict: dict[str, pd.Series] = {}
    tpm_dict: dict[str, pd.Series] = {}
    lengths: dict[str, float] = {}

    for sample_dir in sample_dirs:
        quant_file = sample_dir / "quant.sf"
        if not quant_file.exists():
            LOGGER.warning("未找到 quant.sf，跳过样本: %s", sample_dir.name)
            continue
        try:
            df = pd.read_csv(quant_file, sep="\t")
            missing = {"Name", "TPM", "NumReads", "Length"} - set(df.columns)
            if missing:
                raise SkillError(f"{quant_file} 缺少列: {sorted(missing)}")
            sub = df[["Name", "TPM", "NumReads", "Length"]].copy()
            sub.columns = ["transcript_id", "TPM", "NumReads", "Length"]
            sub = sub.set_index("transcript_id")
            sample = sample_dir.name
            counts_dict[sample] = sub["NumReads"].astype(float)
            tpm_dict[sample] = sub["TPM"].astype(float)
            if not lengths:
                lengths = sub["Length"].astype(float).to_dict()
            LOGGER.info("样本 %s: %d 个转录本", sample, len(sub))
        except SkillError:
            raise
        except Exception as exc:  # noqa: BLE001
            LOGGER.error("处理 %s 失败: %s", quant_file, exc)

    if not counts_dict:
        raise SkillError(f"在 {quant_dir} 中没有找到任何可用的 quant.sf")

    counts_matrix = pd.DataFrame(counts_dict)
    tpm_matrix = pd.DataFrame(tpm_dict)

    # 转录本水平 FPKM = (count * 1e9) / (length * total_reads)
    LOGGER.info("计算转录本水平 FPKM ...")
    fpkm_transcript = pd.DataFrame(index=counts_matrix.index)
    for sample in counts_matrix.columns:
        total = float(counts_matrix[sample].sum())
        if total <= 0:
            fpkm_transcript[sample] = 0.0
            continue
        lens = pd.Series({k: lengths.get(k, 0.0) for k in counts_matrix.index}, dtype=float)
        lens = lens.replace(0.0, np_fallback_length())
        fpkm_transcript[sample] = (counts_matrix[sample] * 1e9) / (lens * total)

    # 汇总到基因水平
    def _to_gene(mat: pd.DataFrame) -> pd.DataFrame:
        if tx2gene:
            mat = mat.copy()
            mat.index = mat.index.map(lambda x: tx2gene.get(str(x), str(x)))
        return mat.groupby(mat.index).sum()

    counts_matrix = _to_gene(counts_matrix)
    tpm_matrix = _to_gene(tpm_matrix)
    fpkm_matrix = _to_gene(fpkm_transcript)

    LOGGER.info("汇总到 %d 个基因", len(counts_matrix))

    # counts 取整（差异表达分析要求整数）
    counts_matrix = counts_matrix.round().astype("int64")

    initial = len(counts_matrix)
    keep = counts_matrix.sum(axis=1) > 0
    if min_counts > 0:
        keep &= counts_matrix.sum(axis=1) >= min_counts
    counts_matrix = counts_matrix[keep]
    LOGGER.info("基因过滤: %d -> %d (阈值 min_counts=%d)", initial, len(counts_matrix), min_counts)

    for sample in counts_matrix.columns:
        LOGGER.info("  样本 %s 总 counts: %s", sample, f"{int(counts_matrix[sample].sum()):,}")

    outputs: dict[str, Path] = {}

    counts_path = output_dir / "counts.csv"
    counts_matrix.index.name = "gene_id"
    counts_matrix.to_csv(counts_path)
    outputs["counts"] = counts_path
    LOGGER.info("counts 矩阵已保存: %s", counts_path)

    common = [g for g in counts_matrix.index if g in tpm_matrix.index]
    if common:
        tpm_out = tpm_matrix.loc[common][counts_matrix.columns]
        tpm_out = tpm_out[tpm_out.sum(axis=1) > 0].round(2)
        tpm_out.index.name = "gene_id"
        tpm_path = output_dir / "tpm.csv"
        tpm_out.to_csv(tpm_path)
        outputs["tpm"] = tpm_path
        LOGGER.info("TPM 矩阵已保存: %s", tpm_path)

    common_f = [g for g in counts_matrix.index if g in fpkm_matrix.index]
    if common_f:
        fpkm_out = fpkm_matrix.loc[common_f][counts_matrix.columns]
        fpkm_out = fpkm_out[fpkm_out.sum(axis=1) > 0].round(2)
        fpkm_out.index.name = "gene_id"
        fpkm_path = output_dir / "fpkm.csv"
        fpkm_out.to_csv(fpkm_path)
        outputs["fpkm"] = fpkm_path
        LOGGER.info("FPKM 矩阵已保存: %s", fpkm_path)

    metadata = build_metadata_template(list(counts_matrix.columns))
    metadata_path = output_dir / "metadata_template.csv"
    metadata.to_csv(metadata_path, index=False)
    outputs["metadata"] = metadata_path
    LOGGER.info("元数据模板已保存: %s", metadata_path)

    return outputs


def np_fallback_length() -> float:
    """转录本长度缺失时的回退值（1000 bp）。"""
    return 1000.0


def build_metadata_template(samples: list[str]) -> pd.DataFrame:
    """由样本名推断 condition 与 replicate，生成元数据模板。

    推断规则：取样本名去掉末尾 ``repN`` 段的部分作为 condition，``repN`` 中的
    数字作为 replicate。原脚本同款逻辑，但补上了"样本名不含 rep 段"时的兜底。

    Args:
        samples: 样本名列表。

    Returns:
        含 ``sample`` / ``condition`` / ``replicate`` 三列的 DataFrame。
    """
    rows = []
    for sample in samples:
        parts = str(sample).split("_")
        condition = ""
        replicate = ""
        if len(parts) >= 2:
            last = parts[-1]
            if "rep" in last.lower():
                condition = "_".join(parts[:-1])
                replicate = last.lower().replace("rep", "").strip("_-")
            else:
                condition = "_".join(parts[:-1])
                replicate = str(len([r for r in rows if r["condition"] == condition]) + 1)
        else:
            condition = str(sample)
            replicate = "1"
        rows.append({"sample": sample, "condition": condition, "replicate": replicate})
    return pd.DataFrame(rows)


def main(args) -> int:
    """CLI 入口：构建 counts 矩阵。"""
    quant_dir = Path(args.quant_dir).expanduser()
    if not quant_dir.exists():
        raise SkillError(f"Salmon 定量目录不存在: {quant_dir}")

    output = Path(args.output).expanduser() if args.output else get_config().output_dir / "counts"

    outputs = build_counts_matrix(
        quant_dir=quant_dir,
        output_dir=output,
        tx2gene_file=Path(args.tx2gene).expanduser() if args.tx2gene else None,
        min_counts=args.min_counts,
    )
    for name, path in outputs.items():
        LOGGER.info("%s -> %s", name, path)
    return 0


if __name__ == "__main__":
    import argparse

    from biopipe_core.cli import cli_parser, run_cli

    DESCRIPTION = "从 Salmon 定量结果构建基因水平 counts/TPM/FPKM 矩阵"

    def build_parser() -> argparse.ArgumentParser:
        """构造带业务参数的解析器（沿用 biopipe_core 的公共参数）。"""
        parser = cli_parser(prog="bulk-rnaseq", description=DESCRIPTION, skill="bulk-rnaseq")
        parser.add_argument("--quant-dir", required=True, help="Salmon 定量结果目录（每样本一个子目录）")
        parser.add_argument("--output", help="输出目录（默认 <results>/counts）")
        parser.add_argument("--tx2gene", help="转录本-基因映射 TSV（列 transcript_id/gene_id）")
        parser.add_argument("--min-counts", type=int, default=0, help="基因总 counts 下限（默认 0）")
        return parser

    sys.exit(run_cli(main, parser=build_parser(), skill="bulk-rnaseq"))
