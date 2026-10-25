#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""多肽序列聚类去冗余入口。

优先调用外部 MMseqs2（``mmseqs easy-cluster``），若不可用则回退到自实现的
序列 identity 聚类：

1. **k-mer 签名**：每条序列提取 k-mer（默认 k=3）集合
2. **相似度矩阵**：用 Jaccard 相似度（|A∩B| / |A∪B|）作为距离度量
3. **层次聚类**：自实现的 average-linkage（UPGMA）凝聚聚类
4. **代表序列选择**：簇内选择与簇心（medoid）identity 最高的序列，
   即"最居中"的序列，保证代表序列对簇内所有成员都具代表性

用法示例::

    # 用 MMseqs2 聚类（若已安装）
    python cluster_peptides.py --input peptides.fasta --output-dir clusters/

    # 强制用内置算法（不依赖 MMseqs2）
    python cluster_peptides.py --input peptides.fasta --method internal \\
        --min-identity 0.3 --output clusters/representatives.fasta

    # 输出代表序列 CSV
    python cluster_peptides.py --input peptides.fasta --method internal \\
        --output clusters/representatives.csv
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "packages"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from biopipe_core import (  # noqa: E402
    SkillError,
    cli_parser,
    get_config,
    get_logger,
    read_fasta,
    resolve_path,
    run_cli,
    write_csv,
    write_fasta,
)

LOGGER = get_logger("biopipe.peptide-design")


# ---------------------------------------------------------------- 序列工具


def kmer_set(sequence: str, k: int = 3) -> frozenset[str]:
    """提取序列的 k-mer 集合（去重）。

    序列短于 k 时返回整条序列作为唯一 k-mer。
    """
    if not sequence:
        return frozenset()
    if len(sequence) < k:
        return frozenset({sequence})
    return frozenset(sequence[i : i + k] for i in range(len(sequence) - k + 1))


def jaccard(a: frozenset[str], b: frozenset[str]) -> float:
    """Jaccard 相似度。"""
    if not a and not b:
        return 1.0
    if not a or not b:
        return 0.0
    inter = len(a & b)
    union = len(a | b)
    return inter / union if union else 0.0


def pairwise_identity(seq_a: str, seq_b: str) -> float:
    """两条序列的全局 identity（无 gap 的字符匹配比例）。

    使用最长公共子序列（Needleman-Wunsch 的简化版：只算匹配不计 gap 惩罚）。
    """
    if not seq_a or not seq_b:
        return 0.0
    if seq_a == seq_b:
        return 1.0
    # LCS 动态规划（滚动数组，O(n*m) 时间、O(min(n,m)) 空间）
    if len(seq_a) < len(seq_b):
        seq_a, seq_b = seq_b, seq_a
    prev = [0] * (len(seq_b) + 1)
    for i in range(1, len(seq_a) + 1):
        cur = [0] * (len(seq_b) + 1)
        char_a = seq_a[i - 1]
        for j in range(1, len(seq_b) + 1):
            if char_a == seq_b[j - 1]:
                cur[j] = prev[j - 1] + 1
            else:
                cur[j] = max(prev[j], cur[j - 1])
        prev = cur
    lcs_len = prev[len(seq_b)]
    return lcs_len / max(len(seq_a), len(seq_b))


# ---------------------------------------------------------------- 聚类


def average_linkage(
    distances: list[list[float]],
    threshold: float,
) -> list[list[int]]:
    """自实现的 average-linkage 凝聚层次聚类（UPGMA）。

    距离矩阵用 1 - similarity 表示；簇间距离取成员对距离的平均。
    合并时使用 Lance-Williams 更新公式::

        d(k, i∪j) = (n_i·d(k,i) + n_j·d(k,j)) / (n_i + n_j)

    因此每次合并只需 O(n) 更新而不是 O(n²) 全量重算。

    当下一次合并的簇间距离超过 ``threshold`` 时停止。

    Returns:
        簇列表，每簇是成员在原始输入中的下标。
    """
    n = len(distances)
    if n == 0:
        return []
    if n == 1:
        return [[0]]

    # 固定槽位：active[i] 为该槽位的成员列表，None 表示已被合并
    active: list[list[int] | None] = [[i] for i in range(n)]
    sizes = [1] * n
    # d[i][j] 为两个槽位之间的距离（i != j 时有效）
    d = [list(row) for row in distances]
    active_count = n

    while active_count > 1:
        best_a = -1
        best_b = -1
        best_dist = None
        for i in range(n):
            if active[i] is None:
                continue
            for j in range(i + 1, n):
                if active[j] is None:
                    continue
                value = d[i][j]
                if best_dist is None or value < best_dist:
                    best_dist = value
                    best_a, best_b = i, j
        if best_a < 0:
            break
        if best_dist is not None and best_dist > threshold:
            break

        # 把 b 合并进 a
        merged_size = sizes[best_a] + sizes[best_b]
        for k in range(n):
            if k == best_a or k == best_b or active[k] is None:
                continue
            new_value = (
                sizes[best_a] * d[k][best_a] + sizes[best_b] * d[k][best_b]
            ) / merged_size
            d[best_a][k] = new_value
            d[k][best_a] = new_value

        active[best_a] = active[best_a] + active[best_b]  # type: ignore[operator]
        active[best_b] = None
        sizes[best_a] = merged_size
        active_count -= 1

    return [members for members in active if members]


def cluster_by_identity(
    sequences: list[str],
    min_identity: float = 0.3,
    kmer_k: int = 3,
    similarity_weight: float = 0.5,
) -> tuple[list[list[int]], list[list[float]]]:
    """基于序列 identity 的聚类（MMseqs2 不可用时的回退实现）。

    步骤：
    1. k-mer Jaccard 相似度作为快速预筛
    2. 对通过预筛的候选对计算精确 identity
    3. 综合相似度 = similarity_weight×identity + (1-w)×Jaccard
    4. average-linkage 凝聚聚类，阈值 = 1 - min_identity

    Returns:
        ``(簇列表, 综合相似度矩阵)``
    """
    n = len(sequences)
    if n == 0:
        return [], []

    signatures = [kmer_set(seq, kmer_k) for seq in sequences]

    # 预筛：Jaccard 相似度下界
    jaccard_lower = max(0.0, (min_identity - similarity_weight) / max(1e-9, 1.0 - similarity_weight))

    combined = [[0.0] * n for _ in range(n)]
    for i in range(n):
        combined[i][i] = 1.0
    for i in range(n):
        for j in range(i + 1, n):
            j_sim = jaccard(signatures[i], signatures[j])
            if j_sim < jaccard_lower:
                continue
            identity = pairwise_identity(sequences[i], sequences[j])
            score = similarity_weight * identity + (1.0 - similarity_weight) * j_sim
            combined[i][j] = score
            combined[j][i] = score

    distances = [[1.0 - combined[i][j] for j in range(n)] for i in range(n)]
    threshold = 1.0 - min_identity
    clusters = average_linkage(distances, threshold)
    return clusters, combined


def select_representatives(
    sequences: list[str],
    clusters: list[list[int]],
    combined: list[list[float]],
    method: str = "medoid",
) -> dict[int, int]:
    """从每个簇中选出代表序列。

    - ``medoid``：与簇内所有成员相似度之和最大的序列（最居中）
    - ``longest``：簇内最长的序列
    - ``first``：簇内第一条序列（保持输入顺序，MMseqs2 风格）

    Returns:
        ``簇下标 -> 代表序列下标`` 映射。
    """
    representatives: dict[int, int] = {}
    for cluster_idx, members in enumerate(clusters):
        if not members:
            continue
        if method == "longest":
            best = max(members, key=lambda i: len(sequences[i]))
        elif method == "first":
            best = min(members)
        else:
            best = max(
                members,
                key=lambda i: sum(combined[i][j] for j in members),
            )
        representatives[cluster_idx] = best
    return representatives


# ---------------------------------------------------------------- MMseqs2


def find_mmseqs() -> str | None:
    """查找可用的 MMseqs2 可执行文件。"""
    for name in ("mmseqs", "mmseqs2", "mmseqs-linux-avx2"):
        path = shutil.which(name)
        if path:
            return path
    return None


def run_mmseqs(
    fasta_path: Path,
    output_dir: Path,
    min_identity: float,
    coverage: float,
    sensitivity: float,
    threads: int,
) -> tuple[dict[str, list[str]], list[list[str]]]:
    """调用 MMseqs2 easy-cluster 完成聚类。

    Returns:
        ``(簇名 -> 成员序列名列表, 代表序列名列表)``
    """
    binary = find_mmseqs()
    if binary is None:
        raise SkillError("未找到 MMseqs2 可执行文件")

    tmp_dir = Path(tempfile.mkdtemp(prefix="biopipe_mmseqs_"))
    try:
        cluster_prefix = tmp_dir / "clu"
        cmd = [
            binary, "easy-cluster",
            str(fasta_path),
            str(cluster_prefix),
            str(tmp_dir),
            "--min-seq-id", str(min_identity),
            "-c", str(coverage),
            "-s", str(sensitivity),
            "--threads", str(threads),
        ]
        LOGGER.info("运行 MMseqs2: %s", " ".join(cmd))
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        if result.returncode != 0:
            raise SkillError(
                f"MMseqs2 执行失败 (exit {result.returncode}): {result.stderr[-800:]}"
            )

        tsv_path = Path(f"{cluster_prefix}_cluster.tsv")
        if not tsv_path.exists():
            raise SkillError(f"MMseqs2 未生成聚类文件: {tsv_path}")

        clusters: dict[str, list[str]] = {}
        for line in tsv_path.read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            rep, member = parts[0], parts[1]
            clusters.setdefault(rep, []).append(member)

        # 复制产物到输出目录
        output_dir.mkdir(parents=True, exist_ok=True)
        for source in tmp_dir.glob("*"):
            if source.is_file():
                shutil.copy2(source, output_dir / source.name)

        representatives = list(clusters.keys())
        return clusters, representatives
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


# ---------------------------------------------------------------- CLI


def build_parser() -> argparse.ArgumentParser:
    """构建命令行解析器。"""
    parser = cli_parser(
        prog="peptide-cluster",
        description="多肽序列聚类去冗余（MMseqs2 优先，内置算法回退）",
    )
    parser.add_argument("--input", required=True, help="输入 FASTA 文件")
    parser.add_argument("--output", help="输出代表序列 FASTA 的路径")
    parser.add_argument("--output-csv", help="输出代表序列明细 CSV 的路径")
    parser.add_argument("--output-dir", help="输出目录（存放中间产物与报告）")
    parser.add_argument(
        "--method",
        default="auto",
        choices=["auto", "mmseqs", "internal"],
        help="聚类方法：auto 优先用 MMseqs2，失败时回退到内置算法",
    )
    parser.add_argument("--min-identity", type=float, default=0.3,
                        help="簇内最小序列 identity")
    parser.add_argument("--coverage", type=float, default=0.8,
                        help="MMseqs2 覆盖度阈值")
    parser.add_argument("--sensitivity", type=float, default=7.5,
                        help="MMseqs2 灵敏度")
    parser.add_argument("--threads", type=int, default=4, help="MMseqs2 线程数")
    parser.add_argument("--kmer-k", type=int, default=3, help="内置算法的 k-mer 长度")
    parser.add_argument("--similarity-weight", type=float, default=0.5,
                        help="内置算法中精确 identity 的权重（其余给 k-mer Jaccard）")
    parser.add_argument(
        "--representative",
        default="medoid",
        choices=["medoid", "longest", "first"],
        help="代表序列选择策略",
    )
    parser.add_argument("--report", help="输出聚类统计报告的路径")
    return parser


def main(args: argparse.Namespace) -> int:
    """主流程。"""
    config = get_config(args.config)
    input_path = resolve_path(args.input)
    if not input_path.exists():
        raise SkillError(f"输入文件不存在: {input_path}")

    records = list(read_fasta(input_path))
    if not records:
        LOGGER.warning("输入 FASTA 没有序列")
        return 0
    LOGGER.info("载入 %d 条序列", len(records))

    ids = [h for h, _ in records]
    sequences = [s for _, s in records]

    output_dir = (
        resolve_path(args.output_dir) if args.output_dir else config.output_dir
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    clusters: list[list[str]] = []
    representatives: list[str] = []
    method_used = ""

    # ---- MMseqs2 路径 --------------------------------------------------
    if args.method in ("auto", "mmseqs"):
        binary = find_mmseqs()
        if binary:
            try:
                cluster_map, reps = run_mmseqs(
                    input_path,
                    output_dir,
                    args.min_identity,
                    args.coverage,
                    args.sensitivity,
                    args.threads,
                )
                clusters = [cluster_map[rep] for rep in reps]
                representatives = reps
                method_used = f"mmseqs ({binary})"
            except SkillError as exc:
                LOGGER.warning("MMseqs2 聚类失败，回退到内置算法: %s", exc)
                if args.method == "mmseqs":
                    raise
        elif args.method == "mmseqs":
            raise SkillError("未找到 MMseqs2 可执行文件（--method mmseqs）")
        else:
            LOGGER.info("未检测到 MMseqs2，使用内置聚类算法")

    # ---- 内置算法路径 --------------------------------------------------
    if not clusters:
        LOGGER.info("使用内置 identity 聚类算法（k=%d）", args.kmer_k)
        index_clusters, combined = cluster_by_identity(
            sequences,
            min_identity=args.min_identity,
            kmer_k=args.kmer_k,
            similarity_weight=args.similarity_weight,
        )
        reps_map = select_representatives(sequences, index_clusters, combined, args.representative)
        clusters = [[ids[i] for i in members] for members in index_clusters]
        representatives = [ids[reps_map[ci]] for ci in sorted(reps_map)]
        method_used = f"internal (k-mer k={args.kmer_k}, identity>={args.min_identity})"

    # ---- 汇总 -----------------------------------------------------------
    id_to_seq = dict(zip(ids, sequences))
    LOGGER.info("聚类完成: %d 条序列 -> %d 个簇（方法: %s）", len(ids), len(clusters), method_used)

    csv_rows = []
    for idx, members in enumerate(clusters):
        rep = representatives[idx] if idx < len(representatives) else members[0]
        for member in members:
            csv_rows.append(
                {
                    "cluster_id": f"c{idx:05d}",
                    "representative": rep,
                    "is_representative": "Yes" if member == rep else "No",
                    "member_id": member,
                    "length": len(id_to_seq.get(member, "")),
                    "sequence": id_to_seq.get(member, ""),
                    "cluster_size": len(members),
                    "method": method_used,
                }
            )

    # ---- 输出 -----------------------------------------------------------
    if args.output:
        fasta_out = resolve_path(args.output)
        write_fasta([(rep, id_to_seq.get(rep, "")) for rep in representatives], fasta_out)
        LOGGER.info("已写出代表序列 FASTA: %s (%d 条)", fasta_out, len(representatives))

    csv_out = (
        resolve_path(args.output_csv)
        if args.output_csv
        else output_dir / "cluster_members.csv"
    )
    if csv_rows:
        write_csv(csv_rows, csv_out, fieldnames=list(csv_rows[0].keys()))
        LOGGER.info("已写出聚类明细 CSV: %s", csv_out)

    # ---- 报告 -----------------------------------------------------------
    report_path = (
        resolve_path(args.report) if args.report else output_dir / "cluster_report.txt"
    )
    lines = [
        "=" * 68,
        "多肽聚类去冗余报告",
        "=" * 68,
        f"聚类方法: {method_used}",
        f"输入序列数: {len(ids)}",
        f"簇数量: {len(clusters)}",
        f"唯一序列数: {len(set(sequences))}",
        f"冗余率: {(1 - len(clusters) / len(ids)) * 100:.1f}%" if ids else "冗余率: N/A",
        "",
    ]
    sizes = sorted((len(c) for c in clusters), reverse=True)
    if sizes:
        lines.append(f"最大簇大小: {sizes[0]}")
        lines.append(f"中位簇大小: {sizes[len(sizes) // 2]}")
        lines.append(f"单序列簇数: {sum(1 for s in sizes if s == 1)}")
        lines.append("")
        lines.append("--- Top 10 簇 ---")
        for idx in range(min(10, len(clusters))):
            lines.append(f"  c{idx:05d}: {len(clusters[idx])} 条 | {representatives[idx]}")
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    LOGGER.info("已写出聚类报告: %s", report_path)

    return 0


if __name__ == "__main__":
    sys.exit(
        run_cli(
            main,
            parser=build_parser(),
            skill="peptide-design",
        )
    )
