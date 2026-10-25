#!/usr/bin/env python3
"""bulk RNA-seq 端到端分析流程。

整合原仓库 ``rnaseq_python/scripts/`` 下的 8 个步骤：

1. **FastQC** —— 原始 reads 质量评估（外部工具 fastqc + multiqc）
2. **fastp** —— 接头去除与质量过滤（外部工具 fastp）
3. **Salmon** —— 转录本定量（外部工具 salmon）
4. **counts 矩阵构建** —— tx2gene 汇总，得到基因水平整数 counts / TPM / FPKM
5. **DESeq2** —— 差异表达分析（``--backend python`` 走内置纯 Python 实现，
   ``--backend r`` 走 Rscript + Bioconductor DESeq2）
6. **GO/KEGG 富集** —— 超几何检验 + BH FDR 校正
7. **可视化** —— 火山图 / MA 图 / 富集条形图 / PCA 图 / 表达热图
8. **报告** —— 汇总所有结果生成文本报告

所有外部工具（fastqc/fastp/salmon）缺失时会给出清晰的跳过提示，不会中断纯计算步骤。

Usage:
    python run_pipeline.py \\
        --input /path/to/fastq \\
        --output /path/to/results \\
        --transcriptome /path/to/transcriptome.fa \\
        --tx2gene /path/to/tx2gene.tsv \\
        --metadata /path/to/metadata.csv \\
        --condition condition --group-a treated --group-b control
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "packages"))

import numpy as np  # noqa: E402
from biopipe_core import SkillError, get_config, get_logger  # noqa: E402
from biopipe_core.cli import cli_parser, run_cli  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))

LOGGER = get_logger("biopipe.bulk-rnaseq")

# 默认各步骤的子目录名
DIR_QC = "01_qc"
DIR_TRIMMED = "02_trimmed"
DIR_SALMON = "03_salmon"
DIR_COUNTS = "04_counts"
DIR_DE = "05_deseq2"
DIR_ENRICH = "06_enrichment"
DIR_FIG = "07_figures"
DIR_REPORT = "08_report"


# --------------------------------------------------------------------------- #
# 外部工具封装
# --------------------------------------------------------------------------- #
def _require_tool(name: str, logger=LOGGER) -> bool:
    """检查外部命令是否可用，不可用时告警并返回 False。"""
    if shutil.which(name) is None:
        logger.warning("未找到外部工具 %s，跳过依赖它的步骤", name)
        return False
    return True


def _run(cmd: list[str], logger=LOGGER, timeout: int | None = None) -> bool:
    """执行外部命令，返回是否成功。

    Args:
        cmd: 命令参数列表（不使用 shell，避免注入问题）。
        logger: 日志器。
        timeout: 超时秒数。

    Returns:
        命令是否成功返回 0。
    """
    logger.info("执行: %s", " ".join(cmd))
    try:
        proc = subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=timeout)
    except FileNotFoundError:
        logger.warning("命令不存在: %s", cmd[0])
        return False
    except subprocess.CalledProcessError as exc:
        logger.error("命令失败 (%s): %s", cmd[0], exc)
        if exc.stderr:
            logger.error("stderr: %s", exc.stderr.strip()[:2000])
        return False
    except subprocess.TimeoutExpired:
        logger.error("命令超时: %s", cmd[0])
        return False
    if proc.stdout and proc.stdout.strip():
        for line in proc.stdout.strip().splitlines()[-5:]:
            logger.info("  %s", line)
    return True


def step_fastqc(fastq_dir: Path, qc_dir: Path, threads: int) -> None:
    """步骤 1：对原始 fastq 运行 FastQC（可选 multiqc 汇总）。"""
    LOGGER.info("[1/8] FastQC 质量评估")
    qc_dir.mkdir(parents=True, exist_ok=True)
    if not _require_tool("fastqc"):
        return
    fastqs = sorted(list(fastq_dir.glob("*.fastq.gz")) + list(fastq_dir.glob("*.fq.gz")))
    if not fastqs:
        LOGGER.warning("在 %s 中未找到 fastq 文件，跳过 FastQC", fastq_dir)
        return
    cmd = ["fastqc", "-t", str(threads), "-o", str(qc_dir)] + [str(f) for f in fastqs]
    _run(cmd)
    if _require_tool("multiqc"):
        _run(["multiqc", str(qc_dir), "-o", str(qc_dir), "--filename", "multiqc_raw_report"])


def step_trimming(fastq_dir: Path, trimmed_dir: Path, qc_dir: Path, threads: int,
                  qual: int = 20, min_len: int = 36) -> None:
    """步骤 2：fastp 质量过滤与接头去除。"""
    LOGGER.info("[2/8] fastp 过滤与接头去除")
    trimmed_dir.mkdir(parents=True, exist_ok=True)
    qc_dir.mkdir(parents=True, exist_ok=True)
    if not _require_tool("fastp"):
        return

    r1_files = sorted(fastq_dir.glob("*_1.fastq.gz")) + sorted(fastq_dir.glob("*_1.fq.gz"))
    if not r1_files:
        LOGGER.warning("在 %s 中未找到 *_1.fastq.gz 配对文件，跳过 trimming", fastq_dir)
        return

    ok = 0
    for r1 in r1_files:
        sample = r1.name.replace("_1.fastq.gz", "").replace("_1.fq.gz", "")
        r2 = r1.parent / f"{sample}_2.fastq.gz"
        if not r2.exists():
            r2 = r1.parent / f"{sample}_2.fq.gz"
        if not r2.exists():
            LOGGER.warning("样本 %s 缺少 R2 文件，跳过", sample)
            continue
        out1 = trimmed_dir / f"{sample}_1.fq.gz"
        out2 = trimmed_dir / f"{sample}_2.fq.gz"
        cmd = [
            "fastp",
            "-i", str(r1), "-I", str(r2),
            "-o", str(out1), "-O", str(out2),
            "--detect_adapter_for_pe",
            "--qualified_quality_phred", str(qual),
            "--length_required", str(min_len),
            "--thread", str(threads),
            "--json", str(qc_dir / f"{sample}.fastp.json"),
            "--html", str(qc_dir / f"{sample}.fastp.html"),
        ]
        if _run(cmd, timeout=3600):
            ok += 1
    LOGGER.info("fastp 完成: %d/%d 个样本", ok, len(r1_files))


def step_salmon_index(transcriptome: Path, index_dir: Path, threads: int) -> bool:
    """构建 Salmon 转录组索引。"""
    LOGGER.info("[3/8] Salmon 索引与转录本定量")
    if not _require_tool("salmon"):
        return False
    index_dir.mkdir(parents=True, exist_ok=True)
    target = index_dir / "salmon_index"
    if target.exists() and any(target.iterdir()):
        LOGGER.info("Salmon 索引已存在，跳过构建: %s", target)
        return True
    ok = _run(["salmon", "index", "-t", str(transcriptome), "-i", str(target),
               "-k", "31", "-p", str(threads)])
    return ok


def step_salmon_quant(trimmed_dir: Path, index_dir: Path, out_dir: Path, threads: int) -> int:
    """对每个样本运行 Salmon 转录本定量。"""
    out_dir.mkdir(parents=True, exist_ok=True)
    index = index_dir / "salmon_index"
    if not index.exists():
        LOGGER.warning("Salmon 索引不存在，跳过定量")
        return 0

    r1_files = sorted(trimmed_dir.glob("*_1.fq.gz"))
    if not r1_files:
        LOGGER.warning("在 %s 中未找到过滤后的 fastq，跳过定量", trimmed_dir)
        return 0

    ok = 0
    for r1 in r1_files:
        sample = r1.name.replace("_1.fq.gz", "")
        r2 = trimmed_dir / f"{sample}_2.fq.gz"
        if not r2.exists():
            LOGGER.warning("样本 %s 缺少 R2，跳过", sample)
            continue
        cmd = [
            "salmon", "quant",
            "-i", str(index), "-l", "A",
            "-1", str(r1), "-2", str(r2),
            "--gcBias", "--seqBias", "--validateMappings",
            "-p", str(threads),
            "-o", str(out_dir / sample),
        ]
        if _run(cmd, timeout=7200):
            ok += 1
    LOGGER.info("Salmon 定量完成: %d/%d", ok, len(r1_files))
    return ok


def step_counts(quant_dir: Path, counts_dir: Path, tx2gene: Path | None) -> dict[str, Path]:
    """步骤 4：构建基因水平 counts / TPM / FPKM 矩阵。"""
    LOGGER.info("[4/8] 构建 counts 矩阵")
    from build_counts_matrix import build_counts_matrix  # 局部导入，避免模块级耦合

    return build_counts_matrix(quant_dir, counts_dir, tx2gene)


# --------------------------------------------------------------------------- #
# 差异表达
# --------------------------------------------------------------------------- #
def _load_counts_and_metadata(
    counts_file: Path,
    metadata_file: Path,
    min_counts: int,
) -> tuple[object, object, list[str]]:
    """读取 counts 矩阵与样本元数据，返回样本 x 基因矩阵、分组标签与样本名。"""
    import pandas as pd

    counts_df = pd.read_csv(counts_file, index_col=0)
    LOGGER.info("原始 counts 矩阵: %d 基因 x %d 样本", counts_df.shape[0], counts_df.shape[1])
    counts_t = counts_df.T
    counts_t.index = counts_t.index.astype(str)

    meta = pd.read_csv(metadata_file)
    if "sample" not in meta.columns:
        raise SkillError(f"元数据缺少 sample 列: {metadata_file}")
    meta = meta.set_index("sample")
    meta.index = meta.index.astype(str)

    common = sorted(set(counts_t.index) & set(meta.index))
    if not common:
        raise SkillError(f"counts 与元数据的样本名完全无法匹配: {counts_file} / {metadata_file}")
    if len(common) < len(counts_t.index):
        LOGGER.warning("仅使用 %d/%d 个共同样本", len(common), len(counts_t.index))
    counts_t = counts_t.loc[common]
    meta = meta.loc[common]

    if "condition" not in meta.columns:
        raise SkillError(f"元数据缺少 condition 列: {metadata_file}")

    keep = counts_t.sum(axis=0) >= min_counts
    LOGGER.info("低表达基因过滤 (<%d counts): 保留 %d / %d", min_counts, int(keep.sum()), len(keep))
    counts_t = counts_t.loc[:, keep]

    groups = meta["condition"].astype(str).tolist()
    return counts_t, meta, groups


def step_deseq2(
    counts_file: Path,
    metadata_file: Path,
    out_dir: Path,
    group_a: str,
    group_b: str,
    alpha: float,
    lfc_threshold: float,
    min_counts: int,
    backend: str,
    rscript: str,
    design: str,
) -> dict[str, Path]:
    """步骤 5：差异表达分析（python 或 r 后端）。"""
    LOGGER.info("[5/8] 差异表达分析 (backend=%s)", backend)
    out_dir.mkdir(parents=True, exist_ok=True)

    counts_t, meta, groups = _load_counts_and_metadata(counts_file, metadata_file, min_counts)
    LOGGER.info("样本分组: %s", dict(zip(counts_t.index.tolist(), groups)))

    if backend == "r":
        return _deseq2_r_backend(
            counts_t, meta, out_dir, group_a, group_b, alpha, lfc_threshold, rscript, design
        )
    return _deseq2_python_backend(
        counts_t, meta, out_dir, group_a, group_b, alpha, lfc_threshold
    )


def _deseq2_python_backend(
    counts_t, meta, out_dir: Path, group_a: str, group_b: str, alpha: float, lfc_threshold: float
) -> dict[str, Path]:
    """内置纯 Python DESeq2 风格实现。"""
    from deseq2_python import run_deseq2

    counts = counts_t.to_numpy(dtype=float)
    res = run_deseq2(
        counts=counts,
        groups=meta["condition"].astype(str).tolist(),
        group_a=group_a,
        group_b=group_b,
        min_counts=1,
        gene_names=[str(g) for g in counts_t.columns],
    )
    table = res.table
    table["significant"] = (table["padj"] < alpha) & (table["log2FoldChange"].abs() > lfc_threshold)
    table["direction"] = "NS"
    table.loc[(table["padj"] < alpha) & (table["log2FoldChange"] > lfc_threshold), "direction"] = "Up"
    table.loc[(table["padj"] < alpha) & (table["log2FoldChange"] < -lfc_threshold), "direction"] = "Down"

    outputs: dict[str, Path] = {}
    full = out_dir / "deseq2_results_full.csv"
    table.to_csv(full)
    outputs["full"] = full

    sig = table[table["significant"]]
    sig_path = out_dir / "deseq2_results_significant.csv"
    sig.to_csv(sig_path)
    outputs["significant"] = sig_path

    import pandas as pd

    sf = pd.DataFrame({"sample": counts_t.index, "size_factor": res.size_factors})
    sf_path = out_dir / "size_factors.csv"
    sf.to_csv(sf_path, index=False)
    outputs["size_factors"] = sf_path

    disp = pd.DataFrame({"gene": table.index, "dispersion": res.alpha_map})
    disp_path = out_dir / "dispersions.csv"
    disp.to_csv(disp_path, index=False)
    outputs["dispersions"] = disp_path

    summary = pd.DataFrame([{
        "total_genes": len(table),
        "significant_genes": int(table["significant"].sum()),
        "upregulated": int((table["direction"] == "Up").sum()),
        "downregulated": int((table["direction"] == "Down").sum()),
        "alpha": alpha,
        "lfc_threshold": lfc_threshold,
        "design": "~condition",
        "contrast": f"{group_a} vs {group_b}",
        "backend": "python",
    }])
    summary_path = out_dir / "deseq2_summary.csv"
    summary.to_csv(summary_path, index=False)
    outputs["summary"] = summary_path

    (out_dir / "dispersion_trend.json").write_text(
        json.dumps(res.trend_params, indent=2), encoding="utf-8"
    )

    LOGGER.info("差异表达完成: %d 个基因, %d 个显著", len(table), int(table["significant"].sum()))
    return outputs


def _deseq2_r_backend(
    counts_t, meta, out_dir: Path, group_a: str, group_b: str,
    alpha: float, lfc_threshold: float, rscript: str, design: str,
) -> dict[str, Path]:
    """Rscript + Bioconductor DESeq2 后端。"""
    import pandas as pd

    if shutil.which(rscript) is None:
        raise SkillError(
            f"未找到 {rscript}，无法使用 r 后端。可改用 --backend python（内置纯 Python 实现）"
        )

    r_counts = out_dir / "_r_counts.tsv"
    counts_out = counts_t.copy()
    counts_out.index.name = "sample"
    counts_out.to_csv(r_counts, sep="\t")

    r_meta = out_dir / "_r_meta.tsv"
    meta_out = meta.copy()
    meta_out.index.name = "sample"
    meta_out.to_csv(r_meta, sep="\t")

    r_script = out_dir / "_r_deseq2.R"
    result_path = out_dir / "deseq2_results_full.csv"
    r_script.write_text(
        f"""
suppressPackageStartupMessages({{library(DESeq2)}})
counts <- as.matrix(read.delim("{r_counts.as_posix()}", row.names = 1, check.names = FALSE))
meta <- read.delim("{r_meta.as_posix()}", row.names = 1, check.names = FALSE)
meta$condition <- factor(meta$condition)
rownames(meta) <- rownames(counts)
dds <- DESeqDataSetFromMatrix(countData = round(counts), colData = meta, design = {design})
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "{group_a}", "{group_b}"), alpha = {alpha})
resdf <- as.data.frame(res)
resdf <- cbind(gene_id = rownames(resdf), resdf)
write.csv(resdf, "{result_path.as_posix()}", row.names = FALSE)
cat("DEDESEQ2_OK\\n")
""".strip(),
        encoding="utf-8",
    )

    LOGGER.info("调用 Rscript 执行 DESeq2 ...")
    proc = subprocess.run([rscript, str(r_script)], capture_output=True, text=True)
    if proc.returncode != 0 or "DEDESEQ2_OK" not in proc.stdout:
        LOGGER.error("R DESeq2 执行失败: %s", (proc.stderr or proc.stdout)[-2000:])
        raise SkillError("R 后端 DESeq2 运行失败")

    table = pd.read_csv(result_path, index_col=0)
    table["significant"] = (table["padj"] < alpha) & (table["log2FoldChange"].abs() > lfc_threshold)
    table["direction"] = "NS"
    table.loc[(table["padj"] < alpha) & (table["log2FoldChange"] > lfc_threshold), "direction"] = "Up"
    table.loc[(table["padj"] < alpha) & (table["log2FoldChange"] < -lfc_threshold), "direction"] = "Down"
    table.to_csv(result_path)

    outputs = {"full": result_path}
    sig_path = out_dir / "deseq2_results_significant.csv"
    table[table["significant"]].to_csv(sig_path)
    outputs["significant"] = sig_path

    summary = pd.DataFrame([{
        "total_genes": len(table),
        "significant_genes": int(table["significant"].sum()),
        "upregulated": int((table["direction"] == "Up").sum()),
        "downregulated": int((table["direction"] == "Down").sum()),
        "alpha": alpha,
        "lfc_threshold": lfc_threshold,
        "design": design,
        "contrast": f"{group_a} vs {group_b}",
        "backend": "r",
    }])
    summary_path = out_dir / "deseq2_summary.csv"
    summary.to_csv(summary_path, index=False)
    outputs["summary"] = summary_path

    LOGGER.info("R DESeq2 完成: %d 个基因, %d 个显著", len(table), int(table["significant"].sum()))
    return outputs


# --------------------------------------------------------------------------- #
# 富集分析
# --------------------------------------------------------------------------- #
def parse_go_string(go_string: str) -> list[tuple[str, str]]:
    """解析 GO 注释字符串。

    支持两种常见格式：

    - ``"process name [GO:0000001]; other [GO:0000002]"``
    - ``"GO:0000001: process name"``

    Args:
        go_string: GO 注释原始字符串。

    Returns:
        ``[(go_id, go_name), ...]`` 列表。
    """
    terms: list[tuple[str, str]] = []
    for raw in go_string.split(";"):
        term = raw.strip()
        if not term:
            continue
        if "[" in term and "]" in term:
            name_part = term.split("[")[0].strip()
            go_id = term.split("[")[1].replace("]", "").strip()
            if go_id.startswith("GO:"):
                terms.append((go_id, name_part))
        else:
            parts = term.split(":")
            if len(parts) >= 2 and parts[0].startswith("GO:"):
                terms.append((parts[0].strip(), ":".join(parts[1:]).strip()))
            else:
                for word in term.split():
                    if word.startswith("GO:") and len(word) == 10:
                        terms.append((word, term.replace(word, "").strip()))
                        break
    return terms


def load_go_annotation(go_file: Path) -> tuple[dict, dict]:
    """加载本地 GO 注释文件。

    自动识别列名（``EnsemblPlants`` / ``gene_id`` / ``GeneID`` 等作为基因 ID 列，
    ``Gene Ontology (biological process)`` 等作为各本体列）。

    Args:
        go_file: GO 注释 CSV 文件。

    Returns:
        ``(gene_go, go_info)``：基因到 ``(ont, go_id, go_name)`` 集合的字典，
        以及 GO ID 到各本体名称的字典。
    """
    from collections import defaultdict

    import pandas as pd

    df = pd.read_csv(go_file)
    LOGGER.info("GO 注释记录: %d 条", len(df))

    gene_col = None
    for cand in ("EnsemblPlants", "gene_id", "GeneID", "gene", "Gene"):
        if cand in df.columns:
            gene_col = cand
            break
    if gene_col is None:
        gene_col = df.columns[0]

    ont_cols = {
        "BP": "Gene Ontology (biological process)",
        "CC": "Gene Ontology (cellular component)",
        "MF": "Gene Ontology (molecular function)",
    }

    gene_go: dict[str, set] = defaultdict(set)
    go_info: dict[str, dict] = defaultdict(dict)

    for _, row in df.iterrows():
        gene_id = row.get(gene_col)
        if gene_id is None or (isinstance(gene_id, float) and pd.isna(gene_id)):
            continue
        gene_id = str(gene_id).strip()
        if not gene_id:
            continue
        # 去掉 Ensembl 风格的版本后缀（.1）以匹配 GTF 中的基因 ID
        if "." in gene_id and gene_id.rsplit(".", 1)[-1].isdigit():
            gene_id = gene_id.rsplit(".", 1)[0]
        for ont, col in ont_cols.items():
            if col not in df.columns:
                continue
            value = row.get(col)
            if value is None or (isinstance(value, float) and pd.isna(value)):
                continue
            for go_id, go_name in parse_go_string(str(value)):
                gene_go[gene_id].add((ont, go_id, go_name))
                go_info[go_id][ont] = go_name

    LOGGER.info("注释基因 %d 个, GO terms %d 个", len(gene_go), len(go_info))
    return gene_go, go_info


def step_enrichment(
    de_file: Path,
    go_file: Path,
    out_dir: Path,
    min_gene_set: int,
    max_gene_set: int,
    padj_threshold: float,
    direction: str = "all",
) -> dict[str, Path]:
    """步骤 6：GO/KEGG 超几何富集分析 + BH FDR 校正。"""
    LOGGER.info("[6/8] GO/KEGG 富集分析")
    import pandas as pd
    from scipy import stats

    out_dir.mkdir(parents=True, exist_ok=True)

    gene_go, go_info = load_go_annotation(go_file)

    # 读取差异基因列表
    de_file = Path(de_file)
    if de_file.suffix.lower() == ".csv":
        de_df = pd.read_csv(de_file, index_col=0)
        if "significant" in de_df.columns:
            if direction != "all" and "direction" in de_df.columns:
                de_df = de_df[de_df["direction"] == direction]
            de_genes = de_df.index.astype(str).tolist()
        else:
            de_genes = de_df.iloc[:, 0].astype(str).tolist()
    else:
        de_genes = [ln.strip() for ln in de_file.read_text(encoding="utf-8").splitlines() if ln.strip()]

    LOGGER.info("差异表达基因: %d 个", len(de_genes))
    if not de_genes:
        LOGGER.warning("没有差异基因，跳过富集分析")
        return {}

    de_set = set(de_genes)
    background = len(gene_go)
    results_by_ont: dict[str, pd.DataFrame] = {}

    for ontology in ("BP", "CC", "MF"):
        term_genes: dict[str, set] = defaultdict(set)
        for gene in de_set:
            for ont, go_id, _name in gene_go.get(gene, ()):
                if ont == ontology:
                    term_genes[go_id].add(gene)

        rows = []
        for go_id, genes_in_term in term_genes.items():
            k_term = len(genes_in_term)
            if k_term < min_gene_set or k_term > max_gene_set:
                continue
            overlap = genes_in_term & de_set
            k = len(overlap)
            if k == 0:
                continue
            # 单侧超几何检验：P(X >= k)
            p_value = float(stats.hypergeom.sf(k - 1, background, k_term, len(de_set)))
            rows.append({
                "GO_ID": go_id,
                "GO_Name": go_info.get(go_id, {}).get(ontology, go_id),
                "Ontology": ontology,
                "Term_Genes": k_term,
                "Overlap_Genes": k,
                "DE_Genes": len(de_set),
                "Background_Genes": background,
                "P_Value": p_value,
                "Overlap_Gene_List": ",".join(sorted(overlap)),
            })

        if not rows:
            LOGGER.warning("%s 未产生富集结果", ontology)
            continue
        res_df = pd.DataFrame(rows)
        # BH FDR 校正
        p = res_df["P_Value"].to_numpy(dtype=float)
        order = np.argsort(p, kind="mergesort")
        ranked = p[order]
        n = len(ranked)
        q = np.minimum.accumulate((ranked * n / np.arange(1, n + 1))[::-1])[::-1]
        padj = np.empty(n)
        padj[order] = np.clip(q, 0, 1)
        res_df["P_Adj"] = padj
        res_df = res_df.sort_values("P_Value").reset_index(drop=True)

        path = out_dir / f"enrichment_{ontology}.csv"
        res_df.to_csv(path, index=False)
        results_by_ont[ontology] = res_df
        LOGGER.info("%s 富集: %d 个 term, %d 个显著 (padj<%.3f)",
                    ontology, len(res_df), int((res_df["P_Adj"] < padj_threshold).sum()), padj_threshold)

    outputs: dict[str, Path] = {}
    if results_by_ont:
        all_df = pd.concat(list(results_by_ont.values()), ignore_index=True).sort_values("P_Value")
        all_path = out_dir / "enrichment_summary_all.csv"
        all_df.to_csv(all_path, index=False)
        outputs["summary"] = all_path
        for ont, _df in results_by_ont.items():
            outputs[ont] = out_dir / f"enrichment_{ont}.csv"
    return outputs


# --------------------------------------------------------------------------- #
# 可视化
# --------------------------------------------------------------------------- #
def step_visualization(
    counts_file: Path,
    metadata_file: Path,
    results_file: Path,
    enrichment_dir: Path,
    fig_dir: Path,
    padj_threshold: float,
    lfc_threshold: float,
    top_n_genes: int = 30,
) -> list[Path]:
    """步骤 7：生成火山图 / MA 图 / 富集条形图 / PCA 图 / 表达热图。"""
    LOGGER.info("[7/8] 可视化")
    import pandas as pd

    fig_dir.mkdir(parents=True, exist_ok=True)
    produced: list[Path] = []

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        LOGGER.warning("matplotlib 未安装，跳过可视化")
        return produced

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans"],
        "font.size": 10,
        "axes.linewidth": 0.8,
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "figure.facecolor": "white",
    })

    colors = {"NS": "#808080", "Up": "#E64B35", "Down": "#3C8DBC"}

    # --- 火山图 ---
    if Path(results_file).exists():
        df = pd.read_csv(results_file, index_col=0)
        df = df.dropna(subset=["log2FoldChange", "padj"])
        if "padj" in df.columns:
            df["color"] = "NS"
            df.loc[(df["padj"] < padj_threshold) & (df["log2FoldChange"] > lfc_threshold), "color"] = "Up"
            df.loc[(df["padj"] < padj_threshold) & (df["log2FoldChange"] < -lfc_threshold), "color"] = "Down"
            fig, ax = plt.subplots(figsize=(5, 4))
            for cname, cval in colors.items():
                sub = df[df["color"] == cname]
                ax.scatter(sub["log2FoldChange"], -np.log10(sub["padj"].clip(lower=1e-300)),
                           c=cval, alpha=0.6, s=18, label=f"{cname} (n={len(sub)})", edgecolors="none")
            ax.axhline(-np.log10(padj_threshold), color="black", linestyle="--", linewidth=0.5)
            ax.axvline(lfc_threshold, color="black", linestyle="--", linewidth=0.5)
            ax.axvline(-lfc_threshold, color="black", linestyle="--", linewidth=0.5)
            ax.set_xlabel("log$_2$FoldChange")
            ax.set_ylabel("-log$_{10}$(padj)")
            ax.set_title("Volcano Plot")
            ax.legend(loc="upper right", frameon=False)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            for suffix in ("pdf", "png"):
                p = fig_dir / f"Figure1_Volcano_Plot.{suffix}"
                fig.savefig(p, format=suffix)
                produced.append(p)
            plt.close(fig)
            LOGGER.info("火山图已生成")

            # --- MA 图 ---
            fig, ax = plt.subplots(figsize=(5, 4))
            for cname, cval in colors.items():
                sub = df[df["color"] == cname]
                ax.scatter(np.log10(sub["baseMean"] + 1), sub["log2FoldChange"],
                           c=cval, alpha=0.6, s=18, label=f"{cname} (n={len(sub)})", edgecolors="none")
            ax.axhline(0, color="black", linewidth=0.5)
            ax.set_xlabel("mean normalized counts (log$_{10}$)")
            ax.set_ylabel("log$_2$FoldChange")
            ax.set_title("MA Plot")
            ax.legend(loc="upper right", frameon=False)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            for suffix in ("pdf", "png"):
                p = fig_dir / f"Figure2_MA_Plot.{suffix}"
                fig.savefig(p, format=suffix)
                produced.append(p)
            plt.close(fig)
            LOGGER.info("MA 图已生成")

    # --- 富集条形图 ---
    enr_file = Path(enrichment_dir) / "enrichment_summary_all.csv"
    if enr_file.exists():
        e = pd.read_csv(enr_file)
        sig = e[e["P_Adj"] < padj_threshold].head(15)
        if len(sig):
            fig, ax = plt.subplots(figsize=(7, max(3, 0.35 * len(sig))))
            y = np.arange(len(sig))
            ax.barh(y, -np.log10(sig["P_Adj"].clip(lower=1e-300)), color="#E64B35")
            ax.set_yticks(y)
            ax.set_yticklabels(sig["GO_Name"], fontsize=8)
            ax.set_xlabel("-log$_{10}$(adjusted P)")
            ax.set_title("GO Enrichment")
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            for suffix in ("pdf", "png"):
                p = fig_dir / f"Figure3_GO_Enrichment.{suffix}"
                fig.savefig(p, format=suffix)
                produced.append(p)
            plt.close(fig)
            LOGGER.info("富集条形图已生成")

    # --- PCA + 表达热图 ---
    if Path(counts_file).exists() and Path(metadata_file).exists():
        try:
            cdf = pd.read_csv(counts_file, index_col=0)
            meta = pd.read_csv(metadata_file)
            if "sample" in meta.columns:
                meta = meta.set_index("sample")
            common = [s for s in cdf.columns if s in meta.index]
            if len(common) >= 3:
                sub = cdf[common]
                # CPM + log1p
                cpm = sub / sub.sum(axis=0).replace(0, np.nan) * 1e6
                log_cpm = np.log1p(cpm.fillna(0))
                mat = log_cpm.T.to_numpy(dtype=float)
                mat = mat - mat.mean(axis=0, keepdims=True)
                u, s, vt = np.linalg.svd(mat, full_matrices=False)
                pcs = u[:, :2] * s[:2]
                fig, ax = plt.subplots(figsize=(5, 4))
                groups = meta.loc[common, "condition"].astype(str) if "condition" in meta.columns else None
                if groups is not None:
                    for cond in sorted(groups.unique()):
                        m = (groups == cond).to_numpy()
                        ax.scatter(pcs[m, 0], pcs[m, 1], label=cond, alpha=0.8, s=40)
                    ax.legend(frameon=False)
                else:
                    ax.scatter(pcs[:, 0], pcs[:, 1], s=40, color="#3C8DBC")
                ax.set_xlabel("PC1")
                ax.set_ylabel("PC2")
                ax.set_title("PCA of samples (CPM, log1p)")
                ax.spines["top"].set_visible(False)
                ax.spines["right"].set_visible(False)
                for suffix in ("pdf", "png"):
                    p = fig_dir / f"Figure4_PCA.{suffix}"
                    fig.savefig(p, format=suffix)
                    produced.append(p)
                plt.close(fig)
                LOGGER.info("PCA 图已生成")

            # 差异基因表达热图
            if Path(results_file).exists():
                rdf = pd.read_csv(results_file, index_col=0).dropna(subset=["padj", "log2FoldChange"])
                rdf = rdf[rdf["padj"] < padj_threshold]
                if len(rdf):
                    top = rdf.reindex(rdf["log2FoldChange"].abs().sort_values(ascending=False).index).head(top_n_genes)
                    genes = [g for g in top.index if g in cdf.index]
                    if genes:
                        hmat = np.log1p(
                            cdf.loc[genes, common] / cdf.loc[genes, common].sum(axis=0).replace(0, np.nan) * 1e6
                        ).fillna(0)
                        fig, ax = plt.subplots(figsize=(6, max(4, 0.22 * len(genes))))
                        im = ax.imshow(hmat.to_numpy(dtype=float), aspect="auto", cmap="RdBu_r")
                        ax.set_yticks(np.arange(len(genes)))
                        ax.set_yticklabels(genes, fontsize=6)
                        ax.set_xticks(np.arange(len(common)))
                        ax.set_xticklabels(common, rotation=90, fontsize=7)
                        ax.set_title(f"Top {len(genes)} DE genes (log1p CPM)")
                        fig.colorbar(im, ax=ax, shrink=0.6, label="log1p CPM")
                        for suffix in ("pdf", "png"):
                            p = fig_dir / f"Figure5_Expression_Heatmap.{suffix}"
                            fig.savefig(p, format=suffix)
                            produced.append(p)
                        plt.close(fig)
                        LOGGER.info("表达热图已生成")
        except Exception as exc:  # noqa: BLE001
            LOGGER.warning("PCA/热图生成失败: %s", exc)

    return produced


# --------------------------------------------------------------------------- #
# 报告
# --------------------------------------------------------------------------- #
def step_report(
    output_root: Path,
    report_path: Path,
    metadata_file: Path | None,
    group_a: str,
    group_b: str,
    backend: str,
) -> Path:
    """步骤 8：汇总所有结果生成文本报告。"""
    LOGGER.info("[8/8] 生成分析报告")
    import pandas as pd

    lines: list[str] = []
    lines.append("=" * 78)
    lines.append("TRANSCRIPTOME ANALYSIS REPORT")
    lines.append("=" * 78)
    lines.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"Output root: {output_root}")
    lines.append(f"Contrast: {group_a} vs {group_b}   (backend: {backend})")

    if metadata_file and Path(metadata_file).exists():
        meta = pd.read_csv(metadata_file)
        lines.append("")
        lines.append("-" * 78)
        lines.append("1. EXPERIMENTAL DESIGN")
        lines.append("-" * 78)
        for _, row in meta.iterrows():
            lines.append(
                f"  {row.get('sample', '?')}: condition={row.get('condition', '?')}, "
                f"replicate={row.get('replicate', '?')}"
            )

    de_dir = output_root / DIR_DE
    sig_file = de_dir / "deseq2_results_significant.csv"
    full_file = de_dir / "deseq2_results_full.csv"
    if full_file.exists():
        lines.append("")
        lines.append("-" * 78)
        lines.append("2. DIFFERENTIAL EXPRESSION")
        lines.append("-" * 78)
        full = pd.read_csv(full_file, index_col=0)
        lines.append(f"  Total genes tested: {len(full)}")
        if sig_file.exists():
            sig = pd.read_csv(sig_file, index_col=0)
            lines.append(f"  Significant DEGs: {len(sig)}")
            if len(sig):
                lines.append(f"    Up-regulated:   {int((sig['log2FoldChange'] > 0).sum())}")
                lines.append(f"    Down-regulated: {int((sig['log2FoldChange'] < 0).sum())}")
                lines.append("")
                lines.append("  Top 10 up-regulated genes:")
                for gene, row in sig.nlargest(10, "log2FoldChange").iterrows():
                    lines.append(f"    {gene}: log2FC={row['log2FoldChange']:.3f}, padj={row['padj']:.3e}")
                lines.append("")
                lines.append("  Top 10 down-regulated genes:")
                for gene, row in sig.nsmallest(10, "log2FoldChange").iterrows():
                    lines.append(f"    {gene}: log2FC={row['log2FoldChange']:.3f}, padj={row['padj']:.3e}")
    else:
        lines.append("")
        lines.append("2. DIFFERENTIAL EXPRESSION: 结果不可用")

    enr_dir = output_root / DIR_ENRICH
    all_enr = enr_dir / "enrichment_summary_all.csv"
    if all_enr.exists():
        lines.append("")
        lines.append("-" * 78)
        lines.append("3. GO ENRICHMENT (hypergeometric test + BH FDR)")
        lines.append("-" * 78)
        e = pd.read_csv(all_enr)
        sig = e[e["P_Adj"] < 0.05].nsmallest(20, "P_Adj")
        lines.append(f"  Tested terms: {len(e)}, significant: {int((e['P_Adj'] < 0.05).sum())}")
        for _, row in sig.iterrows():
            lines.append(
                f"    [{row['Ontology']}] {row['GO_Name']} ({row['GO_ID']}): "
                f"padj={row['P_Adj']:.3e}, overlap={row['Overlap_Genes']}/{row['Term_Genes']}"
            )

    lines.append("")
    lines.append("-" * 78)
    lines.append("4. METHODS")
    lines.append("-" * 78)
    lines.append("""
  Quality control   : FastQC (+ MultiQC) on raw FASTQ
  Read preprocessing : fastp, Phred Q20, min length 36 bp, PE adapter auto-detection
  Quantification    : Salmon (gcBias, seqBias, validateMappings), transcript level
  Count matrix      : tx2gene rollup, NumReads -> integer gene counts; TPM/FPKM also emitted
  Normalization     : DESeq2 median-of-ratios size factors
  Dispersion        : Cox-Reid adjusted NB profile likelihood; parametric trend a0 + a1/mu;
                      MAP shrinkage with log-normal prior
  Testing           : Wald test on the NB-GLM Fisher information; BH FDR adjustment
  Enrichment        : one-sided hypergeometric test with BH FDR correction
""".rstrip())

    lines.append("")
    lines.append("-" * 78)
    lines.append("5. OUTPUT FILES")
    lines.append("-" * 78)
    for sub in (DIR_QC, DIR_TRIMMED, DIR_SALMON, DIR_COUNTS, DIR_DE, DIR_ENRICH, DIR_FIG):
        p = output_root / sub
        if p.exists():
            n = len([x for x in p.rglob("*") if x.is_file()])
            lines.append(f"  {sub}/  ({n} files)")
    lines.append("")
    lines.append("=" * 78)
    lines.append("END OF REPORT")
    lines.append("=" * 78)

    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(lines), encoding="utf-8")
    LOGGER.info("报告已生成: %s", report_path)
    return report_path


# --------------------------------------------------------------------------- #
# 主流程
# --------------------------------------------------------------------------- #
def main(args) -> int:
    """端到端流程主函数。"""
    config = get_config(args.config)
    output_root = Path(args.output).expanduser() if args.output else config.output_dir
    output_root = output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    qc_dir = output_root / DIR_QC
    trimmed_dir = output_root / DIR_TRIMMED
    salmon_dir = output_root / DIR_SALMON
    index_dir = salmon_dir / "index"
    counts_dir = output_root / DIR_COUNTS
    de_dir = output_root / DIR_DE
    enr_dir = output_root / DIR_ENRICH
    fig_dir = output_root / DIR_FIG
    rep_dir = output_root / DIR_REPORT

    fastq_dir = Path(args.input).expanduser().resolve()
    if not fastq_dir.exists():
        raise SkillError(f"输入目录不存在: {fastq_dir}")

    LOGGER.info("输出根目录: %s", output_root)

    # 1-3: 质控 / 过滤 / 定量
    if not args.skip_qc:
        step_fastqc(fastq_dir, qc_dir, args.threads)
    if not args.skip_trimming:
        step_trimming(fastq_dir, trimmed_dir, qc_dir, args.threads,
                      qual=args.quality_phred, min_len=args.min_length)
    else:
        trimmed_dir = Path(args.trimmed_dir).expanduser().resolve() if args.trimmed_dir else trimmed_dir

    if not args.skip_quant:
        transcriptome = Path(args.transcriptome).expanduser() if args.transcriptome else None
        if transcriptome is None or not transcriptome.exists():
            LOGGER.warning("未提供可用的转录组 FASTA，跳过 Salmon 定量")
        else:
            if step_salmon_index(transcriptome, index_dir, args.threads):
                step_salmon_quant(trimmed_dir, index_dir, salmon_dir, args.threads)
    else:
        salmon_dir = Path(args.quant_dir).expanduser().resolve() if args.quant_dir else salmon_dir

    # 4: counts 矩阵
    counts_file = counts_dir / "counts.csv"
    if salmon_dir.exists() and any(salmon_dir.glob("*/quant.sf")):
        tx2gene = Path(args.tx2gene).expanduser() if args.tx2gene else None
        outputs = step_counts(salmon_dir, counts_dir, tx2gene)
        counts_file = outputs.get("counts", counts_file)
    else:
        LOGGER.warning("没有 Salmon 定量结果，跳过 counts 矩阵构建")

    # 5: 差异表达
    metadata_file = Path(args.metadata).expanduser() if args.metadata else counts_dir / "metadata_template.csv"
    de_outputs: dict[str, Path] = {}
    if counts_file.exists() and Path(metadata_file).exists():
        de_outputs = step_deseq2(
            counts_file=counts_file,
            metadata_file=Path(metadata_file),
            out_dir=de_dir,
            group_a=args.group_a,
            group_b=args.group_b,
            alpha=args.alpha,
            lfc_threshold=args.lfc_threshold,
            min_counts=args.min_counts,
            backend=args.backend,
            rscript=args.rscript,
            design=args.design,
        )
    else:
        LOGGER.warning("缺少 counts 矩阵或元数据，跳过差异表达分析")

    # 6: 富集分析
    go_file = Path(args.go_annotation).expanduser() if args.go_annotation else None
    de_for_enrich = de_outputs.get("significant")
    if de_for_enrich and go_file and go_file.exists():
        step_enrichment(
            de_file=de_for_enrich,
            go_file=go_file,
            out_dir=enr_dir,
            min_gene_set=args.min_gene_set,
            max_gene_set=args.max_gene_set,
            padj_threshold=args.alpha,
            direction=args.enrichment_direction,
        )
    else:
        LOGGER.warning("缺少显著 DEG 文件或 GO 注释文件，跳过富集分析")

    # 7: 可视化
    fig_outputs: list[Path] = []
    if counts_file.exists() and Path(metadata_file).exists():
        fig_outputs = step_visualization(
            counts_file=counts_file,
            metadata_file=Path(metadata_file),
            results_file=de_outputs.get("full", de_dir / "deseq2_results_full.csv"),
            enrichment_dir=enr_dir,
            fig_dir=fig_dir,
            padj_threshold=args.alpha,
            lfc_threshold=args.lfc_threshold,
            top_n_genes=args.top_heatmap_genes,
        )

    # 8: 报告
    step_report(
        output_root=output_root,
        report_path=rep_dir / "analysis_report.txt",
        metadata_file=Path(metadata_file),
        group_a=args.group_a,
        group_b=args.group_b,
        backend=args.backend,
    )

    manifest = {
        "generated": datetime.now().isoformat(timespec="seconds"),
        "input": str(fastq_dir),
        "output_root": str(output_root),
        "backend": args.backend,
        "contrast": {"group_a": args.group_a, "group_b": args.group_b},
        "counts_matrix": str(counts_file) if counts_file.exists() else None,
        "figures": [str(p) for p in fig_outputs],
    }
    manifest_path = output_root / "run_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")

    LOGGER.info("=" * 60)
    LOGGER.info("流程完成。输出根目录: %s", output_root)
    LOGGER.info("=" * 60)
    return 0


if __name__ == "__main__":
    DESCRIPTION = (
        "bulk RNA-seq 端到端流程：FastQC -> fastp -> Salmon -> counts 矩阵 -> "
        "DESeq2 差异表达 -> GO 富集 -> 可视化 -> 报告"
    )

    def build_parser() -> argparse.ArgumentParser:
        """构造带业务参数的解析器（沿用 biopipe_core 的公共参数）。"""
        parser = cli_parser(prog="bulk-rnaseq", description=DESCRIPTION, skill="bulk-rnaseq")

        g_in = parser.add_argument_group("输入")
        g_in.add_argument("--input", required=True, help="原始 FASTQ 目录（*.fastq.gz / *.fq.gz）")
        g_in.add_argument("--transcriptome", help="转录组 FASTA（用于构建 Salmon 索引）")
        g_in.add_argument("--tx2gene", help="转录本-基因映射 TSV")
        g_in.add_argument("--metadata", help="样本元数据 CSV（需含 sample 与 condition 列）")
        g_in.add_argument("--go-annotation", help="本地 GO/KEGG 注释 CSV")
        g_in.add_argument("--trimmed-dir", help="已过滤好的 FASTQ 目录（配合 --skip-trimming）")
        g_in.add_argument("--quant-dir", help="已有的 Salmon 定量目录（配合 --skip-quant）")

        g_out = parser.add_argument_group("输出")
        g_out.add_argument("--output", help="输出根目录（默认 <results>）")

        g_run = parser.add_argument_group("流程控制")
        g_run.add_argument("--threads", type=int, default=8, help="并行线程数（默认 8）")
        g_run.add_argument("--skip-qc", action="store_true", help="跳过 FastQC")
        g_run.add_argument("--skip-trimming", action="store_true", help="跳过 fastp")
        g_run.add_argument("--skip-quant", action="store_true", help="跳过 Salmon 定量")
        g_run.add_argument("--quality-phred", type=int, default=20, help="fastp 质量阈值 Phred（默认 20）")
        g_run.add_argument("--min-length", type=int, default=36, help="fastp 最短保留长度（默认 36）")

        g_de = parser.add_argument_group("差异表达")
        g_de.add_argument("--backend", choices=["python", "r"], default="python",
                          help="DESeq2 后端：python=内置实现，r=Rscript+Bioconductor")
        g_de.add_argument("--group-a", default="treated", help="分子组标签（默认 treated）")
        g_de.add_argument("--group-b", default="control", help="分母组标签（默认 control）")
        g_de.add_argument("--design", default="~condition", help="R 后端的设计公式")
        g_de.add_argument("--rscript", default="Rscript", help="Rscript 可执行文件（默认 Rscript）")
        g_de.add_argument("--alpha", type=float, default=0.05, help="FDR 阈值（默认 0.05）")
        g_de.add_argument("--lfc-threshold", type=float, default=1.0,
                          help="|log2FC| 阈值（默认 1.0）")
        g_de.add_argument("--min-counts", type=int, default=10, help="低表达基因过滤阈值（默认 10）")

        g_enr = parser.add_argument_group("富集分析")
        g_enr.add_argument("--min-gene-set", type=int, default=5, help="最小基因集大小（默认 5）")
        g_enr.add_argument("--max-gene-set", type=int, default=500, help="最大基因集大小（默认 500）")
        g_enr.add_argument("--enrichment-direction", choices=["all", "Up", "Down"], default="all",
                           help="富集分析使用的 DEG 方向（默认 all）")

        g_viz = parser.add_argument_group("可视化")
        g_viz.add_argument("--top-heatmap-genes", type=int, default=30, help="热图展示基因数（默认 30）")

        return parser

    sys.exit(run_cli(main, parser=build_parser(), skill="bulk-rnaseq"))
