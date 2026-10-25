#!/usr/bin/env python3
"""蛋白质结构分析入口：PDB 解析 → 几何计算 → 二级结构 → 元数据。

功能：

1. **纯 Python PDB 解析**（不依赖 BioPython）：解析 ATOM / HETATM 记录，
   建立 链-残基-原子 层级，统计元素组成；
2. **几何计算**：回旋半径 Rg（质量加权）、陀螺张量主轴分解、
   形状描述符（各向异性 / 棒状系数）、瞬时偶极矩、端到端距离、包围盒；
3. **二级结构分析**：优先调用外部 DSSP（mkdssp），不可用时回退到内置的
   简化 DSSP 实现（Kabsch-Sander 氢键能量 + φ/ψ 区间判定）；
4. **UniProt 元数据**：可选项，从 ``rest.uniprot.org`` 拉取 accession 对应的
   蛋白名称、物种、长度与功能注释（仅在显式开启时联网）。

用法示例::

    # 解析单个结构并输出 JSON 报告
    python analyze_structure.py --input 1crn.pdb --output results/

    # 只做二级结构，强制使用内置实现
    python analyze_structure.py --input 1crn.pdb --output results/ --dssp-method builtin

    # 批量处理目录并附带 UniProt 注释
    python analyze_structure.py --input pdb_dir/ --output results/ --query-uniprot
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "packages"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from typing import Optional

from biopipe_core import SkillError, cli_parser, get_config, get_logger, write_csv  # noqa: E402
from biopipe_core.cli import run_cli  # noqa: E402
from geometry import (  # noqa: E402
    bounding_box,
    end_to_end_distance,
    instantaneous_dipole,
    radius_of_gyration,
    shape_descriptors,
)
from pdb_parser import parse_pdb_file  # noqa: E402
from secondary_structure import (  # noqa: E402
    assign_secondary_structure,
    find_dssp_executable,
    run_external_dssp,
)

LOGGER = get_logger("biopipe.protein-structure")

SKILL_NAME = "protein-structure"

__all__ = ["analyze_structure", "main"]


def _build_parser() -> argparse.ArgumentParser:
    """构造带业务参数的命令行解析器。"""
    parser = cli_parser(
        prog=SKILL_NAME,
        description="蛋白质结构分析：PDB 解析、几何计算、二级结构判定与 UniProt 注释",
        skill=SKILL_NAME,
    )
    parser.add_argument("--input", required=True, metavar="PATH", help="PDB 文件或包含 PDB 的目录")
    parser.add_argument("--output", metavar="DIR", help="输出目录（默认取配置中的 results/）")
    parser.add_argument(
        "--dssp-method",
        choices=("auto", "external", "builtin"),
        default="auto",
        help="二级结构计算方式：auto 优先外部 DSSP，external 强制外部，builtin 用内置实现",
    )
    parser.add_argument("--dssp-executable", metavar="PATH", help="DSSP 可执行文件路径")
    parser.add_argument(
        "--hbond-energy",
        type=float,
        default=-0.5,
        help="内置 DSSP 判定氢键的最大能量（kcal/mol，默认 -0.5）",
    )
    parser.add_argument(
        "--model",
        type=int,
        default=1,
        help="使用第几个模型（1-based，默认 1）",
    )
    parser.add_argument("--query-uniprot", action="store_true", help="联网查询 UniProt 元数据")
    parser.add_argument(
        "--format",
        choices=("json", "csv", "both"),
        default="both",
        help="输出格式（默认 both）",
    )
    parser.add_argument("--max-files", type=int, default=0, help="最多处理的文件数（0 表示不限制）")
    parser.add_argument("--quiet-summary", action="store_true", help="不输出逐条结构的详细日志")
    return parser


def _resolve_output_dir(output: Optional[str]) -> Path:
    """解析输出目录（命令行 > 配置 > 当前目录）。"""
    config = get_config()
    if output:
        path = Path(output).expanduser()
        if not path.is_absolute():
            path = Path.cwd() / path
    else:
        path = config.output_dir
    path.mkdir(parents=True, exist_ok=True)
    return path


def _collect_pdb_files(input_path: Path, max_files: int) -> list[Path]:
    """收集待处理的 PDB 文件列表。"""
    if input_path.is_dir():
        files = sorted(
            [p for p in input_path.rglob("*") if p.suffix.lower() in (".pdb", ".ent")]
        )
    elif input_path.is_file():
        files = [input_path]
    else:
        raise SkillError(f"输入路径不存在: {input_path}")
    if max_files and max_files > 0:
        files = files[:max_files]
    return files


def _second_structure_for(
    pdb_path: Path,
    structure,
    method: str,
    executable: Optional[str],
    min_hbond_energy: float,
) -> Optional[dict]:
    """按指定方式计算二级结构。"""
    if method in ("auto", "external"):
        exe = executable or find_dssp_executable()
        if exe:
            result = run_external_dssp(pdb_path, executable=exe)
            if result:
                LOGGER.info("%s: 使用外部 DSSP (%s)", pdb_path.name, exe)
                return result
            if method == "external":
                LOGGER.warning("%s: 外部 DSSP 执行失败", pdb_path.name)
                return None
        elif method == "external":
            LOGGER.warning("%s: 未找到 DSSP 可执行文件", pdb_path.name)
            return None
        LOGGER.info("%s: 未找到外部 DSSP，回退到内置实现", pdb_path.name)

    result = assign_secondary_structure(
        structure, min_hbond_energy=min_hbond_energy, use_hbond_pattern=True
    )
    LOGGER.info(
        "%s: 内置 DSSP 完成（螺旋 %.1f%%, 折叠 %.1f%%, 卷曲 %.1f%%）",
        pdb_path.name,
        result["statistics"]["helix_pct"],
        result["statistics"]["sheet_pct"],
        result["statistics"]["coil_pct"],
    )
    return result


def analyze_structure(
    pdb_path: Path,
    dssp_method: str = "auto",
    dssp_executable: Optional[str] = None,
    min_hbond_energy: float = -0.5,
    model_index: int = 1,
    query_uniprot: bool = False,
) -> dict:
    """分析单个 PDB 结构，返回结构化的结果字典。

    Args:
        pdb_path: PDB 文件路径。
        dssp_method: ``auto`` / ``external`` / ``builtin``。
        dssp_executable: 外部 DSSP 可执行文件路径。
        min_hbond_energy: 内置 DSSP 的氢键能量阈值。
        model_index: 使用的模型序号（1-based）。
        query_uniprot: 是否联网查询 UniProt 元数据。

    Returns:
        含 ``structure``（解析概要）、``geometry``、``secondary_structure``、
        ``statistics`` 的字典。
    """
    structure = parse_pdb_file(pdb_path)
    if not structure.models or structure.n_atoms == 0:
        raise SkillError(f"{pdb_path.name}: 未解析到任何原子记录")

    if model_index > 1:
        if model_index > structure.n_models:
            LOGGER.warning(
                "%s: 模型数只有 %d，回退到模型 1", pdb_path.name, structure.n_models
            )
        else:
            structure.models = [structure.models[model_index - 1]]

    atoms = structure.atoms()
    residues = structure.residues()

    geometry = {
        "radius_of_gyration": radius_of_gyration(atoms),
        "shape": shape_descriptors(atoms),
        "dipole": instantaneous_dipole(atoms),
        "end_to_end_distance": end_to_end_distance(atoms),
        "bounding_box": bounding_box(atoms),
    }

    ss_result = _second_structure_for(
        pdb_path, structure, dssp_method, dssp_executable, min_hbond_energy
    )

    chains = []
    for chain in structure.model.chains:
        chains.append(
            {
                "chain_id": chain.chain_id,
                "n_residues": chain.n_residues,
                "n_atoms": chain.n_atoms,
                "sequence": chain.sequence(),
                "ca_count": sum(1 for r in chain if r.get("CA") is not None),
            }
        )

    result = {
        "structure_id": structure.structure_id,
        "source_file": str(pdb_path),
        "analyzed_at": datetime.now().isoformat(timespec="seconds"),
        "structure": structure.summary(),
        "chains": chains,
        "geometry": geometry,
        "secondary_structure": (ss_result or {}).get("secondary_structure", []),
        "ss_statistics": (ss_result or {}).get("statistics", {}),
        "dssp_method": (ss_result or {}).get("method", "none"),
        "residue_types": _residue_type_counts(residues),
    }

    if query_uniprot:
        result["uniprot"] = _query_uniprot(pdb_path.stem)
    else:
        result["uniprot"] = None

    return result


def _residue_type_counts(residues) -> dict:
    """统计残基类型分布。"""
    counts: dict[str, int] = {}
    for residue in residues:
        counts[residue.res_name] = counts.get(residue.res_name, 0) + 1
    return dict(sorted(counts.items(), key=lambda kv: (-kv[1], kv[0])))


def _query_uniprot(accession: str) -> Optional[dict]:
    """按文件名（通常就是 UniProt accession）查询元数据。"""
    try:
        from uniprot_client import fetch_uniprot_entries
    except ImportError as exc:  # pragma: no cover
        LOGGER.warning("无法导入 uniprot_client: %s", exc)
        return None
    rows = fetch_uniprot_entries([accession])
    return rows[0] if rows else None


def _write_outputs(results: list[dict], output_dir: Path, fmt: str) -> list[str]:
    """把结果写成 JSON / CSV 文件，返回文件路径列表。"""
    written: list[str] = []

    if fmt in ("json", "both"):
        json_path = output_dir / "structure_analysis.json"
        payload = {
            "generated_at": datetime.now().isoformat(timespec="seconds"),
            "n_structures": len(results),
            "structures": results,
        }
        json_path.write_text(
            json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8"
        )
        written.append(str(json_path))
        LOGGER.info("已写出 %s", json_path)

    if fmt in ("csv", "both"):
        rows = []
        for item in results:
            struct = item.get("structure", {})
            geom = item.get("geometry", {})
            rg = geom.get("radius_of_gyration", {})
            shape = geom.get("shape", {})
            stats = item.get("ss_statistics", {})
            rows.append(
                {
                    "structure_id": item.get("structure_id", ""),
                    "source_file": item.get("source_file", ""),
                    "n_models": struct.get("n_models", 0),
                    "n_chains": struct.get("n_chains", 0),
                    "n_residues": struct.get("n_residues", 0),
                    "n_atoms": struct.get("n_atoms", 0),
                    "seq_length": struct.get("seq_length", 0),
                    "rg_angstrom": round(rg.get("rg", 0.0), 3),
                    "rg_ca": round(rg["rg_ca"], 3) if rg.get("rg_ca") else "",
                    "shape_class": shape.get("shape_class", ""),
                    "asphericity": round(shape.get("asphericity", 0.0), 4),
                    "dipole_debye": round(geom.get("dipole", {}).get("magnitude", 0.0), 3),
                    "helix_pct": stats.get("helix_pct", 0.0),
                    "sheet_pct": stats.get("sheet_pct", 0.0),
                    "coil_pct": stats.get("coil_pct", 0.0),
                    "dssp_method": item.get("dssp_method", ""),
                    "protein_name": (item.get("uniprot") or {}).get("protein_name", ""),
                    "organism": (item.get("uniprot") or {}).get("organism_name", ""),
                }
            )
        if rows:
            csv_path = output_dir / "structure_analysis.csv"
            write_csv(rows, csv_path)
            written.append(str(csv_path))
            LOGGER.info("已写出 %s", csv_path)

    return written


def main(args: argparse.Namespace) -> int:
    """入口主函数，返回退出码。"""
    input_path = Path(args.input).expanduser()
    if not input_path.is_absolute():
        input_path = Path.cwd() / input_path
    output_dir = _resolve_output_dir(args.output)

    files = _collect_pdb_files(input_path, args.max_files)
    if not files:
        raise SkillError(f"在 {input_path} 下没有找到 PDB 文件")
    LOGGER.info("待处理 %d 个结构文件", len(files))

    results: list[dict] = []
    failures: list[tuple[str, str]] = []
    for pdb_file in files:
        try:
            result = analyze_structure(
                pdb_file,
                dssp_method=args.dssp_method,
                dssp_executable=args.dssp_executable,
                min_hbond_energy=args.hbond_energy,
                model_index=args.model,
                query_uniprot=args.query_uniprot,
            )
        except Exception as exc:  # noqa: BLE001
            LOGGER.warning("分析 %s 失败: %s", pdb_file.name, exc)
            failures.append((pdb_file.name, str(exc)))
            continue
        results.append(result)
        if not args.quiet_summary:
            stats = result.get("ss_statistics", {})
            rg = result["geometry"]["radius_of_gyration"]
            LOGGER.info(
                "%s: %d 链 / %d 残基 / %d 原子, Rg=%.2f Å, 螺旋 %.1f%%",
                result["structure_id"],
                result["structure"]["n_chains"],
                result["structure"]["n_residues"],
                result["structure"]["n_atoms"],
                rg.get("rg", 0.0),
                stats.get("helix_pct", 0.0),
            )

    if not results:
        raise SkillError("没有成功分析任何结构")

    written = _write_outputs(results, output_dir, args.format)

    LOGGER.info("完成：成功 %d / 失败 %d", len(results), len(failures))
    for name, reason in failures:
        LOGGER.warning("失败 %s: %s", name, reason)
    for path in written:
        LOGGER.info("输出: %s", path)
    return 0


if __name__ == "__main__":
    sys.exit(run_cli(main, parser=_build_parser(), skill=SKILL_NAME, version="1.0.0"))
