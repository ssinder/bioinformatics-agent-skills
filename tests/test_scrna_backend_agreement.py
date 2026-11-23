"""单细胞后端一致性测试：内置纯 NumPy 实现 vs scanpy。

与 ``test_backend_agreement.py`` 同理——技能在 SKILL.md 里承诺
"同一数据用 ``--backend scanpy`` 与 ``--backend numpy`` 各跑一次，
聚类数与簇大小分布应大体一致（±1 个簇以内）"，那这个承诺必须由 CI 验证。

夹具里带了一份 ``truth_clusters.csv``（600 个细胞的人工真值聚类），
所以除了两后端互比，还能额外验证两个后端各自对真值的还原度——
这比"两个实现碰巧一致"更强：两个实现可能一起错，但对上真值就说明没错。
"""

from __future__ import annotations

import csv
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
FIXTURES = REPO_ROOT / "tests" / "fixtures_scrna"
PIPELINE = REPO_ROOT / "skills" / "single-cell-rnaseq" / "scripts" / "run_pipeline.py"

sys.path.insert(0, str(REPO_ROOT / "packages"))


def _scanpy_available() -> bool:
    try:
        import anndata  # noqa: F401
        import scanpy  # noqa: F401
    except ImportError:
        return False
    return True


SCANPY_OK = _scanpy_available()
requires_scanpy = pytest.mark.skipif(
    not SCANPY_OK, reason="需要 scanpy + anndata"
)


def _run_backend(backend: str, outdir: Path) -> dict:
    """在夹具数据上跑指定后端，返回 run_manifest 字典。"""
    cmd = [
        sys.executable,
        str(PIPELINE),
        "--input",
        str(FIXTURES / "test_matrix.h5"),
        "--output",
        str(outdir),
        "--marker-file",
        str(FIXTURES / "markers.csv"),
        "--backend",
        backend,
        "--seed",
        "0",
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
    return json.loads((outdir / "run_manifest.json").read_text(encoding="utf-8"))


def _read_clusters(outdir: Path) -> dict[str, int]:
    """读取聚类结果，返回 {cell_id: cluster}。

    产出布局（见 run_pipeline.py 的 DIR_* 常量）：

    - ``01_qc/obs_qc.csv`` —— 每个细胞的 obs（含 cluster / leiden 列）
    - ``03_dimensionality/embeddings.csv`` —— 嵌入坐标
    - ``05_annotated/cluster_annotations.csv`` —— 簇的注释汇总

    细胞级标签在 ``01_qc/obs_qc.csv`` 里——那是流程里第一个落盘的
    obs 快照，包含后续步骤写入的 cluster / leiden / cell_type 列。
    """
    candidates = [
        outdir / "01_qc" / "obs_qc.csv",
        outdir / "04_clustered" / "clusters.csv",
        outdir / "04_clustered" / "obs.csv",
    ]
    for path in candidates:
        if not path.is_file():
            continue
        with open(path, encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            if not reader.fieldnames:
                continue
            id_col = next(
                (c for c in reader.fieldnames if c in ("cell_id", "cell", "barcode", "index")),
                None,
            )
            cl_col = next(
                (c for c in reader.fieldnames if c in ("cluster", "leiden", "louvain", "labels")),
                None,
            )
            if id_col and cl_col:
                clusters: dict[str, int] = {}
                for row in reader:
                    raw_id = row.get(id_col)
                    raw_cl = row.get(cl_col)
                    if raw_id is None or raw_cl in (None, ""):
                        continue
                    try:
                        clusters[str(raw_id)] = int(float(raw_cl))
                    except ValueError:
                        continue
                if clusters:
                    return clusters

    # 兜底：在输出目录里搜任何含 cluster 列的 CSV
    for path in sorted(outdir.rglob("*.csv")):
        try:
            with open(path, encoding="utf-8-sig", newline="") as handle:
                reader = csv.DictReader(handle)
                if not reader.fieldnames:
                    continue
                id_col = next(
                    (c for c in reader.fieldnames if c in ("cell_id", "cell", "barcode")),
                    None,
                )
                cl_col = next(
                    (c for c in reader.fieldnames if c in ("cluster", "leiden", "louvain")),
                    None,
                )
                if not (id_col and cl_col):
                    continue
                clusters = {}
                for row in reader:
                    if row.get(id_col) and row.get(cl_col) not in (None, ""):
                        try:
                            clusters[str(row[id_col])] = int(float(row[cl_col]))
                        except ValueError:
                            continue
                if clusters:
                    return clusters
        except OSError:
            continue

    pytest.fail(f"在 {outdir} 下找不到聚类结果 CSV")
    return {}


def _read_truth() -> dict[str, int]:
    truth: dict[str, int] = {}
    with open(FIXTURES / "truth_clusters.csv", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            truth[row["cell_id"]] = int(row["true_cluster"])
    return truth


def _adjusted_rand_index(a: dict[str, int], b: dict[str, int]) -> float:
    """调整兰德指数（ARI），衡量两种聚类划分的一致性。"""
    common = sorted(set(a) & set(b))
    if not common:
        return 0.0


    contingency: dict[tuple[int, int], int] = {}
    for cell in common:
        key = (a[cell], b[cell])
        contingency[key] = contingency.get(key, 0) + 1

    def comb2(x: int) -> int:
        return x * (x - 1) // 2

    sum_ij = sum(comb2(v) for v in contingency.values())
    sum_i = sum(
        comb2(sum(v for (i, _), v in contingency.items() if i == ai))
        for ai in {k[0] for k in contingency}
    )
    sum_j = sum(
        comb2(sum(v for (_, j), v in contingency.items() if j == bj))
        for bj in {k[1] for k in contingency}
    )
    total = comb2(len(common))

    expected = sum_i * sum_j / total if total else 0.0
    maximum = (sum_i + sum_j) / 2.0
    if maximum == expected:
        return 1.0
    return (sum_ij - expected) / (maximum - expected)


class TestNumpyBackend:
    """内置 NumPy 后端必须能独立跑通。"""

    def test_numpy_backend_runs(self, tmp_path: Path) -> None:
        manifest = _run_backend("numpy", tmp_path / "np")
        assert manifest.get("backend") in ("numpy", "numpy 降级后端", "numpy_backend")
        assert manifest.get("n_cells", 0) > 0

    def test_numpy_backend_writes_markers(self, tmp_path: Path) -> None:
        _run_backend("numpy", tmp_path / "np")
        markers = tmp_path / "np" / "05_annotated" / "marker_genes.csv"
        assert markers.is_file(), "缺少 marker_genes.csv"


@requires_scanpy
class TestScanpyAgreement:
    """scanpy 与内置实现必须给出一致的聚类结构。"""

    @pytest.fixture(scope="class")
    def both(self, tmp_path_factory: pytest.TempPathFactory):
        base = tmp_path_factory.mktemp("scrna_agreement")
        _run_backend("numpy", base / "np")
        _run_backend("scanpy", base / "sc")
        return (
            _read_clusters(base / "np"),
            _read_clusters(base / "sc"),
            base,
        )

    def test_cluster_counts_are_close(self, both) -> None:
        """簇数差异不超过 1（SKILL.md 的承诺）。"""
        np_clusters, sc_clusters, _ = both
        n_np = len(set(np_clusters.values()))
        n_sc = len(set(sc_clusters.values()))
        assert abs(n_np - n_sc) <= 1, (
            f"簇数不一致：内置实现 {n_np} 个，scanpy {n_sc} 个"
        )

    def test_clustering_structure_agrees(self, both) -> None:
        """两个后端的聚类划分应当高度一致（ARI > 0.8）。"""
        np_clusters, sc_clusters, _ = both
        ari = _adjusted_rand_index(np_clusters, sc_clusters)
        assert ari > 0.8, f"两后端聚类结构一致度过低（ARI={ari:.3f}）"

    def test_both_backends_recover_truth(self, both) -> None:
        """两个后端都应对得上夹具里的真值聚类（ARI > 0.6）。"""
        np_clusters, sc_clusters, _ = both
        truth = _read_truth()

        for label, clusters in (("numpy", np_clusters), ("scanpy", sc_clusters)):
            ari = _adjusted_rand_index(truth, clusters)
            assert ari > 0.6, f"{label} 后端与真值聚类一致度过低（ARI={ari:.3f}）"
