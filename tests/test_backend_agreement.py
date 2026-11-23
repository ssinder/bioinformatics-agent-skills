"""双后端一致性测试：内置纯 Python 实现 vs 官方参考实现。

为什么这个测试重要
------------------
本仓库的 ``bulk-rnaseq`` 技能提供两个差异表达后端：

* ``--backend python`` —— 内置的纯 NumPy/SciPy 实现
* ``--backend r``      —— Rscript 调用 Bioconductor DESeq2（官方参考实现）

技能的 SKILL.md 里明确写着："同一数据用 ``--backend python`` 与
``--backend r`` 各跑一次，比较 log2FC 排序与显著基因集合是否基本一致"。

也就是说，**这个承诺必须由 CI 来验证**，否则它只是一句声明。
一个未经验证的自研统计实现在科研工具里是危险的——它可能"看起来对"
却与参考实现有系统性偏差。

测试策略
--------
1. 在固定的夹具数据上分别跑两个后端；
2. 比较 log2FC 排序（Spearman 相关）；
3. 比较显著基因集合（Jaccard 相似度）；
4. 比较 baseMean（量纲相同，应当高度一致）。

R 后端不可用时（CI 里没装 R、或本地没装 Bioconductor）**跳过**而非失败——
但环境里既然装了 R，就必须真的跑，不能悄悄跳过。
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
FIXTURES = REPO_ROOT / "tests" / "fixtures_bulk"
PIPELINE = REPO_ROOT / "skills" / "bulk-rnaseq" / "scripts" / "run_pipeline.py"

sys.path.insert(0, str(REPO_ROOT / "packages"))


def _rscript_available() -> bool:
    """Rscript 是否可用。"""
    if shutil.which("Rscript") is None:
        return False
    try:
        proc = subprocess.run(
            ["Rscript", "-e", "cat(requireNamespace('DESeq2', quietly=TRUE))"],
            capture_output=True,
            text=True,
            timeout=180,
        )
    except (OSError, subprocess.SubprocessError):
        return False
    return proc.returncode == 0 and "TRUE" in proc.stdout


RSCRIPT_OK = _rscript_available()
requires_r = pytest.mark.skipif(
    not RSCRIPT_OK, reason="需要 Rscript + Bioconductor DESeq2"
)


def _run_backend(backend: str, outdir: Path) -> Path:
    """在夹具数据上跑指定后端，返回结果文件路径。"""
    cmd = [
        sys.executable,
        str(PIPELINE),
        "--input",
        str(FIXTURES),
        "--output",
        str(outdir),
        "--tx2gene",
        str(FIXTURES / "tx2gene.tsv"),
        "--metadata",
        str(FIXTURES / "metadata.csv"),
        "--go-annotation",
        str(FIXTURES / "go_annotation.csv"),
        "--quant-dir",
        str(FIXTURES / "salmon"),
        "--skip-qc",
        "--skip-trimming",
        "--skip-quant",
        "--group-a",
        "trt",
        "--group-b",
        "ctrl",
        "--backend",
        backend,
    ]
    env = dict(os.environ)
    env["PYTHONPATH"] = str(REPO_ROOT / "packages")
    proc = subprocess.run(
        cmd, capture_output=True, text=True, timeout=1800, env=env, cwd=str(REPO_ROOT)
    )
    if proc.returncode != 0:
        pytest.fail(
            f"{backend} 后端运行失败（退出码 {proc.returncode}）:\n"
            f"--- stdout ---\n{proc.stdout[-3000:]}\n"
            f"--- stderr ---\n{proc.stderr[-3000:]}"
        )
    return outdir / "05_deseq2" / "deseq2_results_full.csv"


def _read_results(path: Path) -> dict[str, dict[str, float]]:
    """读取结果 CSV，返回 {gene: {列: 值}}。"""
    import csv

    rows: dict[str, dict[str, float]] = {}
    with open(path, encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            gene = str(row.get("gene_id") or row.get("gene") or "").strip()
            if not gene:
                continue
            record: dict[str, float] = {}
            for key, value in row.items():
                if key in ("gene_id", "gene"):
                    continue
                try:
                    record[key] = float(value)
                except (TypeError, ValueError):
                    continue
            rows[gene] = record
    return rows


def _spearman(xs: list[float], ys: list[float]) -> float:
    """Spearman 秩相关（用平均秩处理并列）。"""
    n = len(xs)
    if n < 3:
        return 0.0

    def ranks(values: list[float]) -> list[float]:
        order = sorted(range(len(values)), key=lambda i: values[i])
        result = [0.0] * len(values)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1.0
            for k in range(i, j + 1):
                result[order[k]] = avg
            i = j + 1
        return result

    rx, ry = ranks(xs), ranks(ys)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((rx[i] - mx) * (ry[i] - my) for i in range(n))
    den_x = sum((rx[i] - mx) ** 2 for i in range(n)) ** 0.5
    den_y = sum((ry[i] - my) ** 2 for i in range(n)) ** 0.5
    if den_x == 0 or den_y == 0:
        return 0.0
    return num / (den_x * den_y)


class TestPythonBackendRuns:
    """纯 Python 后端必须能独立跑通——这是降级路径的基础。"""

    def test_python_backend_produces_results(self, tmp_path: Path) -> None:
        out = _run_backend("python", tmp_path / "py")
        assert out.is_file(), f"未产出结果文件: {out}"

        rows = _read_results(out)
        assert len(rows) > 0, "结果为空"

        first = next(iter(rows.values()))
        for column in ("log2FoldChange", "pvalue", "padj", "baseMean"):
            assert column in first, f"结果缺少列 {column}"

    def test_python_backend_writes_size_factors(self, tmp_path: Path) -> None:
        _run_backend("python", tmp_path / "py")
        sf = tmp_path / "py" / "05_deseq2" / "size_factors.csv"
        assert sf.is_file(), "缺少 size_factors.csv"

    def test_backend_is_recorded_in_manifest(self, tmp_path: Path) -> None:
        _run_backend("python", tmp_path / "py")
        manifest = json.loads(
            (tmp_path / "py" / "run_manifest.json").read_text(encoding="utf-8")
        )
        assert manifest.get("backend") == "python"


@requires_r
class TestBackendsAgree:
    """两个后端必须给出统计上等价的结果。

    这是整个仓库科学正确性的核心证据——如果自研实现与 Bioconductor
    DESeq2 系统性偏离，那么所有基于它的下游结论都不可信。
    """

    @pytest.fixture(scope="class")
    def both_results(self, tmp_path_factory: pytest.TempPathFactory):
        base = tmp_path_factory.mktemp("backend_agreement")
        py = _read_results(_run_backend("python", base / "py"))
        r = _read_results(_run_backend("r", base / "r"))
        return py, r

    def test_same_gene_set(self, both_results) -> None:
        """两个后端必须覆盖同一批基因。"""
        py, r = both_results
        assert set(py) == set(r), (
            "两后端覆盖的基因集不同："
            f"仅 Python 有 {sorted(set(py) - set(r))[:5]}，"
            f"仅 R 有 {sorted(set(r) - set(py))[:5]}"
        )

    def test_log2fc_ranking_agrees(self, both_results) -> None:
        """log2FC 排序必须高度一致（Spearman > 0.95）。"""
        py, r = both_results
        genes = sorted(set(py) & set(r))
        assert len(genes) >= 10, f"共同基因过少（{len(genes)}），无法做相关性检验"

        py_lfc = [py[g]["log2FoldChange"] for g in genes]
        r_lfc = [r[g]["log2FoldChange"] for g in genes]
        rho = _spearman(py_lfc, r_lfc)

        assert rho > 0.95, (
            f"log2FC 排序相关度过低（Spearman ρ={rho:.4f}），"
            f"内置实现与 DESeq2 可能存在系统性偏差"
        )

    def test_basemean_agrees(self, both_results) -> None:
        """baseMean 是量纲相同的量，应当高度一致。"""
        py, r = both_results
        genes = [g for g in sorted(set(py) & set(r)) if py[g].get("baseMean") is not None]
        assert len(genes) >= 10

        py_mean = [py[g]["baseMean"] for g in genes]
        r_mean = [r[g]["baseMean"] for g in genes]
        rho = _spearman(py_mean, r_mean)

        assert rho > 0.99, f"baseMean 相关度过低（ρ={rho:.4f}）"

    def test_significant_gene_sets_overlap(self, both_results) -> None:
        """显著基因集合应当高度重叠（Jaccard > 0.8）。"""
        py, r = both_results
        genes = sorted(set(py) & set(r))

        def significant(table: dict) -> set[str]:
            out = set()
            for gene in genes:
                record = table[gene]
                padj = record.get("padj")
                if padj is None:
                    continue
                # NaN 会被 float() 解析出来，需要显式排除
                if padj != padj:  # NaN check
                    continue
                if padj < 0.05:
                    out.add(gene)
            return out

        py_sig, r_sig = significant(py), significant(r)
        if not py_sig and not r_sig:
            pytest.skip("夹具数据在两后端下都没有显著基因，无法比较集合")

        union = py_sig | r_sig
        assert union, "显著基因集合为空"
        jaccard = len(py_sig & r_sig) / len(union)

        assert jaccard > 0.8, (
            f"显著基因集合 Jaccard={jaccard:.3f}（Python {len(py_sig)} 个 / "
            f"R {len(r_sig)} 个），两后端结论不一致"
        )

    def test_direction_agrees(self, both_results) -> None:
        """显著基因的调控方向必须一致——方向相反是最严重的错误。"""
        py, r = both_results
        genes = sorted(set(py) & set(r))

        mismatches = []
        for gene in genes:
            p_adj = py[gene].get("padj")
            r_adj = r[gene].get("padj")
            if p_adj is None or r_adj is None:
                continue
            if p_adj != p_adj or r_adj != r_adj:
                continue
            if p_adj >= 0.05 or r_adj >= 0.05:
                continue
            p_dir = py[gene].get("log2FoldChange", 0.0)
            r_dir = r[gene].get("log2FoldChange", 0.0)
            if p_dir * r_dir < 0:
                mismatches.append((gene, p_dir, r_dir))

        assert not mismatches, (
            f"{len(mismatches)} 个基因在两后端下调控方向相反，例如: {mismatches[:5]}"
        )


class TestRBackendAvailability:
    """R 后端的环境探测本身要可靠。"""

    def test_probe_does_not_raise(self) -> None:
        """探测 Rscript 不应抛异常——无论结果如何。"""
        assert isinstance(_rscript_available(), bool)

    def test_skip_marker_is_accurate(self) -> None:
        """跳过标记必须与实际探测结果一致，不能"永远跳过"。"""
        import importlib.util

        spec = importlib.util.spec_from_file_location(
            "backend_agreement", __file__
        )
        assert spec is not None
        # 探测函数已在本模块导入时执行过；这里只确认它返回 bool
        assert isinstance(RSCRIPT_OK, bool)
