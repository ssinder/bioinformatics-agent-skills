"""分子动力学轨迹分析入口。

用法::

    python run_analysis.py --input traj.dcd --topology top.pdb --output results/

支持的分析项：RMSD（逐原子 / 逐残基）、RMSF、回旋半径、氢键、PCA、
二面角分布、Ramachandran 图、轨迹聚类、链间相互作用、SASA。

若 MDAnalysis 可用，优先使用它读取轨迹；否则使用内置的纯 numpy DCD / PDB
轨迹解析器。核心算法（Kabsch 叠合、SVD PCA、层次聚类、几何氢键判据）
为内置实现，不依赖 MDAnalysis。
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "packages"))

from biopipe_core.optional_deps import optional_import  # noqa: E402

np = optional_import('numpy')
np_HINT = '数值计算基础库'


def _require_np():
    """确保 numpy 可用，否则抛出带安装指引的异常。"""
    if np is None:
        from biopipe_core.optional_deps import MissingDependencyError

        raise MissingDependencyError('numpy', hint=np_HINT)
    return np



from biopipe_core import cli_parser, get_config, get_logger, write_csv  # noqa: E402
from biopipe_core.cli import SkillError, run_cli  # noqa: E402
from biopipe_core.optional_deps import is_available  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))

import md_analysis as md  # noqa: E402
import trajectory_io as tio  # noqa: E402

SKILL_NAME = "molecular-dynamics"

# 各元素的质量（用于质量加权 Rg），单位 amu
ATOMIC_MASS = {
    "H": 1.008, "D": 2.014, "C": 12.011, "N": 14.007, "O": 15.999,
    "S": 32.06, "P": 30.974, "F": 18.998, "CL": 35.45, "BR": 79.904,
    "I": 126.904, "SE": 78.971, "FE": 55.845, "ZN": 65.38, "MG": 24.305,
    "CA": 40.078, "NA": 22.990, "K": 39.098,
}


def _build_parser() -> argparse.ArgumentParser:
    """构造命令行解析器。"""
    parser = cli_parser(
        prog=SKILL_NAME,
        description="分子动力学轨迹分析：RMSD / RMSF / Rg / 氢键 / PCA / 聚类 / Ramachandran / SASA",
        skill=SKILL_NAME,
    )
    parser.add_argument("--input", required=True, metavar="PATH", help="轨迹文件路径（DCD / PDB / XYZ）")
    parser.add_argument("--topology", metavar="PATH", help="拓扑文件路径（PDB），提供原子/残基信息")
    parser.add_argument("--output", metavar="DIR", help="输出目录（默认取配置中的 results/）")
    parser.add_argument("--start", type=int, default=0, help="起始帧（默认 0）")
    parser.add_argument("--stop", type=int, default=None, help="结束帧（不含）")
    parser.add_argument("--stride", type=int, default=1, help="帧步长（默认 1）")
    parser.add_argument(
        "--selection",
        metavar="SEL",
        help="分析选择串，如 'protein' / 'name CA' / 'resname LIG'（默认 protein）",
    )
    parser.add_argument(
        "--analysis",
        default="all",
        help="要执行的分析项，逗号分隔：rmsd,rmsf,rg,hbond,pca,dihedral,ramachandran,cluster,interchain,sasa（默认 all）",
    )
    parser.add_argument("--fit-reference", type=int, default=0, help="RMSD 叠合的参考帧（默认 0）")
    parser.add_argument("--hbond-angle", type=float, default=120.0, help="氢键最小 D-H···A 角度（度，默认 120）")
    parser.add_argument("--hbond-h-dist", type=float, default=2.5, help="H···A 最大距离（Å，默认 2.5）")
    parser.add_argument("--hbond-da-dist", type=float, default=3.5, help="D···A 最大距离（Å，默认 3.5）")
    parser.add_argument("--n-clusters", type=int, default=None, help="聚类簇数（默认按阈值自动）")
    parser.add_argument(
        "--cluster-threshold",
        type=float,
        default=None,
        help="聚类 RMSD 阈值（Å），给定则优先于 --n-clusters",
    )
    parser.add_argument("--contact-dist", type=float, default=4.5, help="链间接触距离阈值（Å，默认 4.5）")
    parser.add_argument("--sasa-points", type=int, default=100, help="SASA 每原子球面采样点数（默认 100）")
    parser.add_argument("--max-cluster-frames", type=int, default=200, help="参与聚类的最大帧数（默认 200）")
    parser.add_argument(
        "--format", choices=("json", "csv", "both"), default="both", help="输出格式（默认 both）"
    )
    return parser


# ---------------------------------------------------------------------------
# 选择串解析（轻量子集）
# ---------------------------------------------------------------------------


def _parse_selection(selection: str | None, topology: tio.SimpleTopology) -> np.ndarray:
    """解析选择串，返回原子索引数组。

    支持：

    * ``protein`` / ``all``  -- 全部原子
    * ``name CA`` / ``name N`` -- 按原子名
    * ``resname LIG``          -- 按残基名
    * ``resid 1-20``          -- 按残基号区间
    * ``chain A``             -- 按链

    Args:
        selection: 选择串。
        topology: 拓扑对象。

    Returns:
        原子索引数组（升序）。空选择时返回全部原子。
    """
    if not selection:
        selection = "protein"
    sel = selection.strip()
    if sel.lower() in ("all", "protein", "*"):
        return np.arange(topology.n_atoms, dtype=int)

    tokens = sel.split()
    if not tokens:
        return np.arange(topology.n_atoms, dtype=int)

    keyword = tokens[0].lower()
    value = tokens[1] if len(tokens) > 1 else ""
    idx: list[int] = []
    for atom in topology.atoms:
        keep = False
        if keyword == "name":
            keep = atom.name.strip().upper() == value.upper()
        elif keyword == "resname":
            keep = atom.resname.strip().upper() == value.upper()
        elif keyword == "chain":
            keep = atom.chain.strip() == value
        elif keyword == "resid":
            keep = _match_resid(value, atom.resid)
        elif keyword == "element":
            keep = str(atom.element).strip().upper() == value.upper()
        elif keyword == "protein":
            keep = atom.resname.strip().upper() not in ("HOH", "WAT", "SOL", "TIP3", "NA")
        idx.append(atom.index) if keep else None
    if not idx:
        raise SkillError(f"选择串 '{selection}' 未匹配到任何原子")
    return np.asarray(idx, dtype=int)


def _match_resid(spec: str, resid: int) -> bool:
    """匹配残基号，支持 ``5``、``1-20``、``1,5,9``。"""
    for part in spec.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            lo, hi = part.split("-", 1)
            try:
                if int(lo) <= resid <= int(hi):
                    return True
            except ValueError:
                continue
        else:
            try:
                if int(part) == resid:
                    return True
            except ValueError:
                continue
    return False


# ---------------------------------------------------------------------------
# 分析主流程
# ---------------------------------------------------------------------------


def _analyze(args: argparse.Namespace) -> dict:
    """执行轨迹分析，返回结果字典。"""
    logger = get_logger(f"biopipe.{SKILL_NAME}")
    traj = tio.load_trajectory(
        args.input, start=args.start, stop=args.stop, stride=args.stride
    )
    logger.info(
        "已加载轨迹 %s：%d 帧 / %d 原子 / %d 残基（读取器 %s）",
        traj.source, traj.n_frames, traj.n_atoms, traj.topology.n_residues, traj.reader,
    )

    topology = traj.topology
    if args.topology:
        topo_path = Path(args.topology)
        if not topo_path.exists():
            raise SkillError(f"拓扑文件不存在: {topo_path}")
        try:
            topo_traj = tio.read_pdb_trajectory(topo_path)
            if topo_traj.n_atoms != traj.n_atoms:
                raise SkillError(
                    f"拓扑原子数 ({topo_traj.n_atoms}) 与轨迹原子数 ({traj.n_atoms}) 不一致"
                )
            topology = topo_traj.topology
            logger.info("已从拓扑文件读取原子信息：%s", topo_path)
        except SkillError:
            raise
        except Exception as exc:  # noqa: BLE001
            raise SkillError(f"解析拓扑文件失败: {exc}") from exc

    elements = tio.atom_elements(topology)
    chain_ids = tio.atom_chain_ids(topology)
    selection = _parse_selection(args.selection, topology)
    logger.info("分析选择包含 %d 个原子", selection.size)

    coords = traj.coords[:, selection, :]
    weights = np.array(
        [ATOMIC_MASS.get(e, 12.011) for e in [elements[i] for i in selection]],
        dtype=np.float64,
    )

    requested = {s.strip().lower() for s in args.analysis.split(",") if s.strip()}
    if "all" in requested:
        requested = {
            "rmsd", "rmsf", "rg", "hbond", "pca",
            "dihedral", "ramachandran", "cluster", "interchain", "sasa",
        }

    results: dict[str, Any] = {
        "input": str(args.input),
        "topology": str(args.topology) if args.topology else None,
        "n_frames": traj.n_frames,
        "n_atoms_total": traj.n_atoms,
        "n_atoms_selected": int(selection.size),
        "reader": traj.reader,
        "selection": args.selection or "protein",
        "analyses": sorted(requested),
    }
    frame_rows: list[dict] = []

    # --- RMSD ------------------------------------------------------------
    bb = tio.extract_backbone(topology)
    if "rmsd" in requested and traj.n_frames > 1:
        ref_idx = min(max(args.fit_reference, 0), traj.n_frames - 1)
        ref = coords[ref_idx]
        per_frame = np.zeros(traj.n_frames)
        aligned_all = np.empty_like(coords)
        for f in range(traj.n_frames):
            rot, trans, _ = md.kabsch_align(coords[f], ref, weights)
            aligned_all[f] = coords[f] @ rot.T + trans
            per_frame[f] = md.rmsd(coords[f], ref, align=True, weights=weights)
        results["rmsd"] = {
            "reference_frame": ref_idx,
            "per_frame": per_frame.tolist(),
            "mean": float(per_frame.mean()),
            "std": float(per_frame.std()),
            "min": float(per_frame.min()),
            "max": float(per_frame.max()),
            "final": float(per_frame[-1]),
        }
        # 逐残基 RMSD：先整体叠合，再统计每个 Cα 相对参考帧的偏差
        res_rmsd: list[dict] = []
        for k, res_idx in enumerate(bb["residue_index"]):
            ca_global = int(bb["ca"][k])
            ca_local = _to_local(ca_global, selection)
            if ca_local < 0:
                continue
            dev = np.linalg.norm(aligned_all[:, ca_local, :] - ref[ca_local], axis=1)
            res_rmsd.append(
                {
                    "residue_index": int(res_idx),
                    "resid": int(topology.residues[res_idx].resid),
                    "resname": topology.residues[res_idx].resname,
                    "chain": topology.residues[res_idx].chain,
                    "ca_index": ca_global,
                    "rmsd_mean": round(float(dev.mean()), 4),
                    "rmsd_max": round(float(dev.max()), 4),
                }
            )
        results["rmsd"]["per_residue_backbone"] = res_rmsd
        logger.info("RMSD: 均值 %.3f Å, 末帧 %.3f Å", per_frame.mean(), per_frame[-1])

    # --- RMSF ------------------------------------------------------------
    if "rmsf" in requested and traj.n_frames > 1:
        values = md.rmsf(coords, align_reference=min(max(args.fit_reference, 0), traj.n_frames - 1))
        results["rmsf"] = {
            "per_atom": values.tolist(),
            "mean": float(values.mean()),
            "max": float(values.max()),
            "max_atom": int(int(np.argmax(values))),
        }
        logger.info("RMSF: 均值 %.3f Å, 最大 %.3f Å", values.mean(), values.max())

    # --- 回旋半径 --------------------------------------------------------
    if "rg" in requested:
        rg = md.radius_of_gyration(coords, weights)
        results["radius_of_gyration"] = {
            "per_frame": rg.tolist(),
            "mean": float(rg.mean()),
            "std": float(rg.std()),
            "min": float(rg.min()),
            "max": float(rg.max()),
        }
        logger.info("Rg: 均值 %.3f Å", rg.mean())

    # --- 氢键 ------------------------------------------------------------
    if "hbond" in requested:
        # 氢键判据在选择子集的局部索引上运行，因此残基映射也要投影到局部索引
        local_res_map: dict[int, tuple[str, int]] = {}
        for res in topology.residues:
            for ai in res.atoms:
                loc = _to_local(int(ai), selection)
                if loc >= 0:
                    local_res_map[loc] = (res.chain, int(res.resid))
        md.register_residue_map(local_res_map)
        sub_elements = [elements[i] for i in selection]
        hb = md.hbond_analysis(
            coords,
            sub_elements,
            h_dist_max=args.hbond_h_dist,
            ha_dist_max=args.hbond_da_dist,
            angle_min=args.hbond_angle,
        )
        results["hbond"] = {
            "per_frame_count": hb["per_frame_count"].tolist(),
            "mean_count": hb["mean_count"],
            "criteria": hb["criteria"],
            "n_candidate_pairs": hb["n_candidate_pairs"],
            "occupancy": hb["occupancy"][:100],
        }
        logger.info("氢键: 平均 %.1f 个/帧", hb["mean_count"])

    # --- PCA -------------------------------------------------------------
    if "pca" in requested and traj.n_frames > 2:
        aligned = np.empty_like(coords)
        ref = coords[min(max(args.fit_reference, 0), traj.n_frames - 1)]
        for f in range(traj.n_frames):
            rot, trans, _ = md.kabsch_align(coords[f], ref, weights)
            aligned[f] = coords[f] @ rot.T + trans
        pca = md.pca_analysis(aligned, n_components=min(3, traj.n_frames - 1))
        results["pca"] = {
            "eigenvalues": pca["eigenvalues"].tolist(),
            "explained_variance_ratio": pca["explained_variance_ratio"].tolist(),
            "cumulative_variance": pca["cumulative_variance"].tolist(),
            "projection": pca["projection"].tolist(),
        }
        logger.info(
            "PCA: PC1/PC2 解释方差 %s",
            np.round(pca["explained_variance_ratio"][:2], 4).tolist(),
        )

    # --- 二面角 / Ramachandran -------------------------------------------
    if "dihedral" in requested or "ramachandran" in requested:
        phi_arr, psi_arr = _collect_phi_psi(coords, bb, selection)
        if "dihedral" in requested:
            results["dihedral"] = {
                "description": "主链 phi = C(i-1)-N(i)-CA(i)-C(i)，psi = N(i)-CA(i)-C(i)-N(i+1)",
                "phi_per_frame": np.nan_to_num(phi_arr, nan=0.0).tolist(),
                "psi_per_frame": np.nan_to_num(psi_arr, nan=0.0).tolist(),
                "phi_circular_mean": round(float(md._circular_mean(phi_arr)), 3),
                "psi_circular_mean": round(float(md._circular_mean(psi_arr)), 3),
            }
        if "ramachandran" in requested and phi_arr.size >= 2:
            rama = md.ramachandran(phi_arr, psi_arr)
            results["ramachandran"] = {
                "n_valid": rama["n_valid"],
                "fractions": rama["fractions"],
                "counts": rama["counts"],
                "outliers": rama["outliers"],
                "histogram": rama["histogram"],
            }
            logger.info(
                "Ramachandran: alpha %.1f%% / beta %.1f%% / coil %.1f%%",
                rama["fractions"]["alpha"] * 100,
                rama["fractions"]["beta"] * 100,
                rama["fractions"]["coil"] * 100,
            )

    # --- 聚类 ------------------------------------------------------------
    if "cluster" in requested and traj.n_frames > 2:
        matrix, frame_idx = md.rmsd_matrix(coords, max_frames=args.max_cluster_frames)
        cluster = md.hierarchical_cluster(
            matrix, n_clusters=args.n_clusters, threshold=args.cluster_threshold
        )
        results["clustering"] = {
            "n_clusters": cluster["n_clusters"],
            "frame_indices": frame_idx,
            "labels": cluster["labels"],
            "clusters": cluster["clusters"],
            "rmsd_matrix": matrix.tolist(),
        }
        logger.info("轨迹聚类: %d 个簇", cluster["n_clusters"])

    # --- 链间相互作用 ----------------------------------------------------
    if "interchain" in requested and len(set(chain_ids)) > 1:
        inter = md.interchain_contacts(
            coords, [chain_ids[i] for i in selection], contact_dist=args.contact_dist
        )
        results["interchain"] = inter
        logger.info("链间接触: %d 对链存在稳定接触", len(inter["pairs"]))

    # --- SASA ------------------------------------------------------------
    if "sasa" in requested:
        sasa = md.sasa_estimate(
            coords, [elements[i] for i in selection], n_points=args.sasa_points
        )
        results["sasa"] = {
            "per_frame_total": sasa["per_frame_total"].tolist(),
            "per_frame_polar": sasa["per_frame_polar"].tolist(),
            "per_frame_nonpolar": sasa["per_frame_nonpolar"].tolist(),
            "mean_total": sasa["mean_total"],
            "std_total": sasa["std_total"],
            "n_points": sasa["n_points"],
            "probe_radius": sasa["probe_radius"],
        }
        logger.info("SASA: 平均 %.1f Å²", sasa["mean_total"])

    # --- 逐帧汇总表 ------------------------------------------------------
    for f in range(traj.n_frames):
        row = {"frame": f, "time_ps": float(traj.times[f])}
        if "rmsd" in results:
            row["rmsd"] = round(results["rmsd"]["per_frame"][f], 4)
        if "rg" in results:
            row["rg"] = round(results["radius_of_gyration"]["per_frame"][f], 4)
        if "hbond" in results:
            row["n_hbond"] = int(results["hbond"]["per_frame_count"][f])
        if "sasa" in results:
            row["sasa"] = round(results["sasa"]["per_frame_total"][f], 2)
        frame_rows.append(row)

    results["frames"] = frame_rows
    return results


def _to_local(global_index: int, selection: np.ndarray) -> int:
    """全局原子索引 -> 选择子集内的局部索引（找不到返回 -1）。"""
    pos = np.searchsorted(selection, global_index)
    if pos < selection.size and selection[pos] == global_index:
        return int(pos)
    return -1


def _local_dihedral(coords: np.ndarray, a: int, b: int, c: int, d: int) -> np.ndarray:
    """逐帧计算局部索引的四原子二面角（索引非法时返回 NaN）。"""
    out = np.full(coords.shape[0], np.nan)
    if min(a, b, c, d) < 0:
        return out
    for f in range(coords.shape[0]):
        p = coords[f]
        out[f] = md._dihedral(p[a], p[b], p[c], p[d])
    return out


def _collect_phi_psi(
    coords: np.ndarray, bb: dict, selection: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """收集所有残基的 (phi, psi) 时间序列，展平为样本集合。"""
    n_res = len(bb["residue_index"])
    phis: list[np.ndarray] = []
    psis: list[np.ndarray] = []
    for k in range(n_res):
        n_i, ca_i, c_i = int(bb["n"][k]), int(bb["ca"][k]), int(bb["c"][k])
        if n_i < 0 or ca_i < 0 or c_i < 0:
            continue
        if k > 0:
            c_prev = int(bb["c"][k - 1])
            if c_prev >= 0:
                phi = _local_dihedral(
                    coords,
                    _to_local(c_prev, selection),
                    _to_local(n_i, selection),
                    _to_local(ca_i, selection),
                    _to_local(c_i, selection),
                )
                phis.append(phi)
        if k + 1 < n_res:
            n_next = int(bb["n"][k + 1])
            if n_next >= 0:
                psi = _local_dihedral(
                    coords,
                    _to_local(n_i, selection),
                    _to_local(ca_i, selection),
                    _to_local(c_i, selection),
                    _to_local(n_next, selection),
                )
                psis.append(psi)
    phi_arr = np.concatenate(phis) if phis else np.array([])
    psi_arr = np.concatenate(psis) if psis else np.array([])
    size = min(phi_arr.size, psi_arr.size)
    return phi_arr[:size], psi_arr[:size]


def _write_outputs(results: dict, output_dir: Path, fmt: str) -> list[str]:
    """写出 JSON / CSV 结果。"""
    logger = get_logger(f"biopipe.{SKILL_NAME}")
    written: list[str] = []
    if fmt in ("json", "both"):
        path = output_dir / "md_analysis.json"
        path.write_text(json.dumps(results, indent=2, ensure_ascii=False), encoding="utf-8")
        written.append(str(path))
        logger.info("已写出 %s", path)
    if fmt in ("csv", "both") and results.get("frames"):
        path = output_dir / "md_analysis_frames.csv"
        write_csv(results["frames"], path)
        written.append(str(path))
        logger.info("已写出 %s", path)
    return written


def _resolve_output_dir(output: str | None) -> Path:
    """解析输出目录。"""
    config = get_config()
    if output:
        path = Path(output).expanduser()
        if not path.is_absolute():
            path = Path.cwd() / path
    else:
        path = config.output_dir
    path.mkdir(parents=True, exist_ok=True)
    return path


def main(args: argparse.Namespace) -> int:
    """技能主入口。"""
    # 依赖守卫：numpy 缺失时给出明确安装指引，
    # 而不是让 --help 之类的元信息查询也一并崩溃。
    if not is_available("numpy"):
        raise SkillError(
            "缺少依赖 numpy，轨迹分析无法运行。"
            "\n请先安装: pip install numpy"
        )
    logger = get_logger(f"biopipe.{SKILL_NAME}")
    output_dir = _resolve_output_dir(args.output)
    results = _analyze(args)
    _write_outputs(results, output_dir, args.format)

    logger.info("完成：%d 帧 / %d 个选中原子", results["n_frames"], results["n_atoms_selected"])
    if "rmsd" in results:
        logger.info("  RMSD 均值 %.3f Å", results["rmsd"]["mean"])
    if "rmsf" in results:
        logger.info("  RMSF 均值 %.3f Å", results["rmsf"]["mean"])
    if "radius_of_gyration" in results:
        logger.info("  Rg 均值 %.3f Å", results["radius_of_gyration"]["mean"])
    if "hbond" in results:
        logger.info("  氢键平均 %.1f 个/帧", results["hbond"]["mean_count"])
    if "pca" in results:
        logger.info(
            "  PC1 解释方差 %.1f%%",
            results["pca"]["explained_variance_ratio"][0] * 100,
        )
    if "clustering" in results:
        logger.info("  聚类簇数 %d", results["clustering"]["n_clusters"])
    if "sasa" in results:
        logger.info("  SASA 平均 %.1f Å²", results["sasa"]["mean_total"])
    return 0


if __name__ == "__main__":
    sys.exit(run_cli(main, parser=_build_parser(), skill=SKILL_NAME, version="1.0.0"))
