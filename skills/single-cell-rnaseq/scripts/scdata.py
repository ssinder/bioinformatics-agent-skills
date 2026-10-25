#!/usr/bin/env python3
"""轻量 AnnData 容器与 10x / 矩阵文件读取。

``scanpy`` / ``anndata`` 未必在所有环境可用，而单细胞流程的第一步
（把原始数据读进来）必须有稳定的行为。因此本模块实现一个最小化的
AnnData 兼容子集：

- ``X``：表达矩阵（``scipy.sparse.csr_matrix`` 或稠密 ``ndarray``）
- ``obs``：细胞级元数据（``pandas.DataFrame``）
- ``var``：基因级元数据（``pandas.DataFrame``）
- ``obsm`` / ``varm`` / ``uns``：降维坐标、载荷与额外信息

支持的数据格式：

1. **10x Genomics HDF5**（``*.h5``）——  ``matrix/barcodes``、``matrix/data``、
   ``matrix/indices``、``matrix/indptr``、``matrix/shape``、``matrix/features/name``
2. **10x MTX 目录**（``matrix.mtx`` + ``barcodes.tsv`` + ``features.tsv``）
3. **CSV / TSV 稠密矩阵**（首行是基因名，首列是细胞名，或反过来）
4. **``.h5ad``** —— 需要 ``anndata``，缺失时给出明确报错

Usage:
    from scdata import read_matrix, SimpleAnnData
    adata = read_matrix(Path("data/filtered_feature_bc_matrix.h5"))
    print(adata.n_obs, adata.n_vars)
"""

from __future__ import annotations

import gzip
import io
import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import scipy.sparse as sp

__all__ = [
    "SimpleAnnData",
    "read_matrix",
    "read_10x_h5",
    "read_10x_mtx",
    "read_dense_table",
    "to_anndata",
    "from_anndata",
]


class SimpleAnnData:
    """最小化的 AnnData 兼容容器。

    Attributes:
        X: 表达矩阵（sparse CSR 或稠密 ndarray）。
        obs: 细胞级元数据 DataFrame。
        var: 基因级元数据 DataFrame。
        obsm: 细胞级嵌入字典（如 ``X_pca`` / ``X_umap``）。
        varm: 基因级载荷字典。
        uns: 其他信息字典。
        raw: 可选的原始计数副本。
    """

    def __init__(
        self,
        X: Any,
        obs: pd.DataFrame,
        var: pd.DataFrame,
        obsm: dict[str, Any] | None = None,
        varm: dict[str, Any] | None = None,
        uns: dict[str, Any] | None = None,
    ) -> None:
        self.X = X
        self.obs = obs
        self.var = var
        self.obsm = obsm or {}
        self.varm = varm or {}
        self.uns = uns or {}
        self.raw: SimpleAnnData | None = None

    # ------------------------------------------------------------------ #
    @property
    def n_obs(self) -> int:
        """细胞数。"""
        return int(self.obs.shape[0])

    @property
    def n_vars(self) -> int:
        """基因数。"""
        return int(self.var.shape[0])

    @property
    def shape(self) -> tuple[int, int]:
        """``(n_obs, n_vars)``。"""
        return (self.n_obs, self.n_vars)

    @property
    def obs_names(self) -> pd.Index:
        """细胞名。"""
        return self.obs.index

    @property
    def var_names(self) -> pd.Index:
        """基因名。"""
        return self.var.index

    # ------------------------------------------------------------------ #
    def copy(self) -> SimpleAnnData:
        """深拷贝（矩阵复制，元数据复制）。"""
        X = self.X.copy() if hasattr(self.X, "copy") else np.array(self.X, copy=True)
        return SimpleAnnData(
            X,
            self.obs.copy(),
            self.var.copy(),
            {k: (v.copy() if hasattr(v, "copy") else v) for k, v in self.obsm.items()},
            {k: (v.copy() if hasattr(v, "copy") else v) for k, v in self.varm.items()},
            dict(self.uns),
        )

    def subset_cells(self, mask: np.ndarray) -> SimpleAnnData:
        """按布尔掩码筛选细胞，返回新对象。"""
        mask = np.asarray(mask, dtype=bool)
        X = self.X[mask]
        return SimpleAnnData(
            X,
            self.obs.loc[mask].copy(),
            self.var.copy(),
            {k: (np.asarray(v)[mask] if len(v) == self.n_obs else v) for k, v in self.obsm.items()},
            {k: v for k, v in self.varm.items()},
            dict(self.uns),
        )

    def subset_genes(self, mask: np.ndarray) -> SimpleAnnData:
        """按布尔掩码筛选基因，返回新对象。"""
        mask = np.asarray(mask, dtype=bool)
        X = self.X[:, mask]
        return SimpleAnnData(
            X,
            self.obs.copy(),
            self.var.loc[mask].copy(),
            {k: v for k, v in self.obsm.items()},
            {k: (np.asarray(v)[mask] if len(v) == self.n_vars else v) for k, v in self.varm.items()},
            dict(self.uns),
        )

    def toarray(self) -> np.ndarray:
        """返回稠密表达矩阵。"""
        if sp.issparse(self.X):
            return self.X.toarray()
        return np.asarray(self.X)

    def var_names_make_unique(self) -> None:
        """把重复的基因名改写为唯一名（原名 + ``-1``、``-2`` …）。"""
        self.var.index = _make_unique(self.var.index)

    def write_h5ad(self, path: str | Path) -> Path:
        """写出为 h5ad（需要 anndata）。

        Args:
            path: 输出路径。

        Returns:
            输出文件路径。
        """
        try:
            import anndata as ad
        except ImportError as exc:
            raise RuntimeError(
                "写出 .h5ad 需要 anndata，请先安装: pip install anndata"
            ) from exc
        adata = ad.AnnData(
            X=self.X,
            obs=self.obs.copy(),
            var=self.var.copy(),
            obsm={k: np.asarray(v) for k, v in self.obsm.items()},
            uns={k: _sanitize_uns(v) for k, v in self.uns.items()},
        )
        out = Path(path)
        out.parent.mkdir(parents=True, exist_ok=True)
        adata.write_h5ad(out)
        return out


def _make_unique(index: pd.Index) -> pd.Index:
    """把重复标签改写为唯一标签。"""
    seen: dict[str, int] = {}
    new: list[str] = []
    for label in index:
        key = str(label)
        if key in seen:
            seen[key] += 1
            new.append(f"{key}-{seen[key]}")
        else:
            seen[key] = 0
            new.append(key)
    return pd.Index(new, name=index.name)


def _sanitize_uns(obj: Any) -> Any:
    """把 uns 中的值转成 anndata 可序列化类型。"""
    if isinstance(obj, dict):
        return {str(k): _sanitize_uns(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_sanitize_uns(v) for v in obj]
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj
    if isinstance(obj, (str, int, float, bool)) or obj is None:
        return obj
    return str(obj)


# --------------------------------------------------------------------------- #
# 读取
# --------------------------------------------------------------------------- #
def _decode(values: Any) -> list[str]:
    """把 h5py 中的字节串数组解码为 str 列表。"""
    out: list[str] = []
    for v in values:
        if isinstance(v, bytes):
            out.append(v.decode("utf-8"))
        else:
            out.append(str(v))
    return out


def read_10x_h5(path: str | Path) -> SimpleAnnData:
    """读取 10x Genomics HDF5 矩阵。

    Args:
        path: ``.h5`` 文件路径。

    Returns:
        :class:`SimpleAnnData` 实例。
    """
    import h5py

    path = Path(path)
    with h5py.File(path, "r") as f:
        if "matrix" not in f:
            raise ValueError(f"{path} 不是 10x HDF5 格式（缺少 /matrix 组）")
        matrix = f["matrix"]
        barcodes = _decode(matrix["barcodes"][:])
        data = matrix["data"][:]
        indices = matrix["indices"][:]
        indptr = matrix["indptr"][:]
        shape = matrix["shape"][:]

        if "features" in matrix and "name" in matrix["features"]:
            gene_names = _decode(matrix["features"]["name"][:])
        else:
            gene_names = [f"gene_{i}" for i in range(int(shape[0]))]

    # 10x 存的是 (n_genes, n_cells)
    csr = sp.csr_matrix((data, indices, indptr), shape=(int(shape[0]), int(shape[1])))
    csr = csr.T.tocsr()  # 转为 (n_cells, n_genes)

    obs = pd.DataFrame(index=pd.Index(barcodes, name="cell_id"))
    var = pd.DataFrame(index=pd.Index(gene_names, name="gene_id"))
    return SimpleAnnData(csr, obs, var)


def _read_lines(path: Path) -> list[str]:
    """读取文本文件的全部行（自动处理 gzip 与 BOM）。"""
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            return [ln.rstrip("\n") for ln in handle]
    return path.read_text(encoding="utf-8-sig").splitlines()


def read_10x_mtx(directory: str | Path) -> SimpleAnnData:
    """读取 10x MTX 目录。

    需要包含 ``matrix.mtx(.gz)``、``barcodes.tsv(.gz)`` 与
    ``features.tsv(.gz)``（或 ``genes.tsv``）三个文件。

    Args:
        directory: 目录路径。

    Returns:
        :class:`SimpleAnnData` 实例。
    """
    directory = Path(directory)
    if not directory.is_dir():
        raise ValueError(f"{directory} 不是目录")

    def _find(*names: str) -> Path | None:
        for name in names:
            for candidate in (name, name + ".gz"):
                p = directory / candidate
                if p.exists():
                    return p
        return None

    mtx_file = _find("matrix.mtx")
    barcode_file = _find("barcodes.tsv")
    feature_file = _find("features.tsv", "genes.tsv")
    if mtx_file is None or barcode_file is None or feature_file is None:
        raise ValueError(
            f"{directory} 缺少 10x MTX 必需文件（matrix.mtx / barcodes.tsv / features.tsv）"
        )

    barcodes = [ln.split("\t")[0] for ln in _read_lines(barcode_file) if ln.strip()]
    feature_lines = [ln for ln in _read_lines(feature_file) if ln.strip()]
    gene_names = [ln.split("\t")[0] for ln in feature_lines]

    with open(mtx_file, encoding="utf-8") as handle:
        header = handle.readline()
        if not header.startswith("%"):
            raise ValueError(f"{mtx_file} 不是合法的 MTX 文件")
        n_rows, n_cols, nnz = (int(x) for x in handle.readline().split())
        rows: list[int] = []
        cols: list[int] = []
        vals: list[float] = []
        for line in handle:
            if not line.strip():
                continue
            parts = line.split()
            rows.append(int(parts[0]) - 1)
            cols.append(int(parts[1]) - 1)
            vals.append(float(parts[2]))

    # MTX 行=基因，列=细胞
    mtx = sp.csr_matrix((vals, (rows, cols)), shape=(n_rows, n_cols)).T.tocsr()

    obs = pd.DataFrame(index=pd.Index(barcodes, name="cell_id"))
    var = pd.DataFrame(index=pd.Index(gene_names, name="gene_id"))
    return SimpleAnnData(mtx, obs, var)


def read_dense_table(path: str | Path, sep: str | None = None) -> SimpleAnnData:
    """读取稠密的 CSV/TSV 表达矩阵。

    自动判断方向：若行数明显多于列数（基因数 > 细胞数），则视为"行=基因"，
    转置为"行=细胞"。

    Args:
        path: 表格文件路径。
        sep: 分隔符，None 时自动嗅探（制表符或逗号）。

    Returns:
        :class:`SimpleAnnData` 实例。
    """
    path = Path(path)
    text = path.read_text(encoding="utf-8-sig", errors="replace")
    first_line = text.splitlines()[0] if text.strip() else ""
    if sep is None:
        sep = "\t" if first_line.count("\t") >= first_line.count(",") else ","

    df = pd.read_csv(path, sep=sep, index_col=0)
    df.index = df.index.astype(str)
    df.columns = df.columns.astype(str)

    # 去掉可能存在的非数值列（如基因名重复列）
    numeric = df.select_dtypes(include=[np.number])
    if numeric.shape[1] == 0:
        raise ValueError(f"{path} 没有数值列，无法作为表达矩阵读取")

    if numeric.shape[0] > numeric.shape[1]:
        # 行=基因，列=细胞 -> 转置
        numeric = numeric.T

    obs = pd.DataFrame(index=pd.Index(numeric.index.astype(str), name="cell_id"))
    var = pd.DataFrame(index=pd.Index(numeric.columns.astype(str), name="gene_id"))
    return SimpleAnnData(numeric.to_numpy(dtype=float), obs, var)


def read_h5ad(path: str | Path) -> SimpleAnnData:
    """读取 ``.h5ad`` 文件（需要 anndata）。"""
    try:
        import anndata as ad
    except ImportError as exc:
        raise RuntimeError(
            f"读取 {path} 需要 anndata，请先安装: pip install anndata；"
            "或改用 10x HDF5 / MTX / CSV 格式"
        ) from exc
    adata = ad.read_h5ad(path)
    obs = pd.DataFrame(adata.obs.copy())
    var = pd.DataFrame(adata.var.copy())
    obsm = {k: np.asarray(v) for k, v in adata.obsm.items()}
    uns = {k: _sanitize_uns(v) for k, v in adata.uns.items()}
    return SimpleAnnData(adata.X, obs, var, obsm=obsm, uns=uns)


def to_anndata(adata: SimpleAnnData) -> Any:
    """把 :class:`SimpleAnnData` 转换为官方的 ``anndata.AnnData``。

    scanpy 后端需要真正的 ``AnnData`` 对象；本函数负责在两种容器之间
    无损搬运 ``X`` / ``obs`` / ``var`` / ``obsm`` / ``uns``。

    Args:
        adata: 轻量容器。

    Returns:
        ``anndata.AnnData`` 实例。

    Raises:
        RuntimeError: 缺少 ``anndata`` 时抛出。
    """
    try:
        import anndata as ad
    except ImportError as exc:
        raise RuntimeError(
            "scanpy 后端需要 anndata，请先安装: pip install scanpy anndata"
        ) from exc

    out = ad.AnnData(
        X=adata.X,
        obs=adata.obs.copy(),
        var=adata.var.copy(),
        obsm={k: np.asarray(v) for k, v in adata.obsm.items() if k != "distances"},
    )
    for key, value in adata.uns.items():
        out.uns[key] = value
    if adata.raw is not None:
        out.raw = adata.raw.X
    return out


def from_anndata(adata: Any) -> SimpleAnnData:
    """把官方的 ``anndata.AnnData`` 转换回 :class:`SimpleAnnData`。

    Args:
        adata: ``anndata.AnnData`` 实例。

    Returns:
        :class:`SimpleAnnData` 实例。
    """
    obs = pd.DataFrame(adata.obs.copy())
    var = pd.DataFrame(adata.var.copy())
    obsm: dict[str, Any] = {}
    for key, value in adata.obsm.items():
        arr = np.asarray(value)
        # distances 是稀疏邻接矩阵，保留为对象引用即可
        obsm[key] = arr
    uns = {k: _sanitize_uns(v) for k, v in adata.uns.items()}
    out = SimpleAnnData(adata.X, obs, var, obsm=obsm, uns=uns)
    if adata.raw is not None:
        out.raw = SimpleAnnData(
            adata.raw.X,
            pd.DataFrame(index=adata.obs_names.copy()),
            pd.DataFrame(index=adata.var_names.copy()),
        )
    return out


def read_matrix(path: str | Path) -> SimpleAnnData:
    """按扩展名/内容自动选择读取方式。

    Args:
        path: 输入路径（文件或 10x MTX 目录）。

    Returns:
        :class:`SimpleAnnData` 实例。
    """
    path = Path(path)
    if path.is_dir():
        return read_10x_mtx(path)

    suffix = path.suffix.lower()
    if suffix == ".h5":
        return read_10x_h5(path)
    if suffix == ".h5ad":
        return read_h5ad(path)
    if suffix in (".csv", ".tsv", ".txt", ".csv.gz", ".tsv.gz"):
        return read_dense_table(path)
    # 尝试按内容嗅探
    try:
        return read_10x_h5(path)
    except Exception:  # noqa: BLE001
        return read_dense_table(path)


# --------------------------------------------------------------------------- #
# 工具
# --------------------------------------------------------------------------- #
def gene_is_mitochondrial(var_names: pd.Index, prefix: str = "MT-") -> np.ndarray:
    """识别线粒体基因。

    同时支持人类（``MT-``）与小鼠（``mt-``）命名，符号匹配不区分大小写。

    Args:
        var_names: 基因名索引。
        prefix: 线粒体基因前缀。

    Returns:
        布尔数组。
    """
    upper = var_names.astype(str).str.upper()
    return upper.str.startswith(prefix.upper()) | upper.str.startswith("MT_")


def gene_is_hb(var_names: pd.Index) -> np.ndarray:
    """识别血红蛋白基因（``HB*``，常见污染标志）。

    Args:
        var_names: 基因名索引。

    Returns:
        布尔数组。
    """
    upper = var_names.astype(str).str.upper()
    return upper.str.match(r"^HB[ABDEGQZ][0-9A-Z\-]*$")


def summarize_sparsity(X: Any) -> dict[str, float]:
    """统计表达矩阵的稀疏度信息。

    Args:
        X: 表达矩阵。

    Returns:
        含 ``shape`` / ``nnz`` / ``density`` / ``sparsity`` 的字典。
    """
    if sp.issparse(X):
        shape = X.shape
        nnz = int(X.nnz)
    else:
        arr = np.asarray(X)
        shape = arr.shape
        nnz = int(np.count_nonzero(arr))
    total = max(shape[0] * shape[1], 1)
    return {
        "shape": shape,
        "nnz": nnz,
        "density": nnz / total,
        "sparsity": 1.0 - nnz / total,
    }


def strip_version(gene_id: str) -> str:
    """去掉 Ensembl 风格的版本后缀（``ENSG000001.2`` -> ``ENSG000001``）。"""
    return re.sub(r"\.\d+$", "", str(gene_id))
