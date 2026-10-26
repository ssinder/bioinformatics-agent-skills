#!/usr/bin/env python3
"""从 AlphaFold 预测结果中选出最佳模型并转换为 PDB。

对应老仓库 ``dock_lightdock/scripts/extract_best_alphafold.py``。

老脚本的问题
------------
1. **硬编码个人路径**。下载目录与输入目录都写死在模块常量里。
2. **只认一种命名**。只匹配 ``fold_{design}_model_*.cif`` 与
   ``fold_{design}_summary_confidences_*.json``，且用字符串 ``replace``
   从模型文件名反推置信度文件名——一旦命名规则不同就彻底失效。
   实际上 AlphaFold Server / ColabFold / AlphaFold DB 有好几套命名约定。
3. **只取第一个 CIF**。``sorted(glob(...))[0]`` 在多个模型存在时取到的是
   **字典序最小**的那个，不是置信度最高的那个。
4. **转换失败只返回 0 原子**。调用方据此判定"失败"，但没有区分
   "文件本身是空的"与"解析器出错了"。
5. **没有汇总表**。批量处理多个设计时，只能看 stdout 才知道谁最好。

本脚本的做法
------------
- **不依赖具体命名**：在目录里找出所有结构文件（``.cif`` / ``.pdb``）与
  所有置信度文件（``*summary_confidences*.json``），按**排名编号**配对。
- 按 **ranking_score** 排序选最优，并给出完整的模型排名表。
- CIF → PDB 转换优先用 Biopython；没装 Biopython 时退回到**内置的
  mmCIF 原子记录解析器**（只提取 ``ATOM``/``HETATM`` 行，转成 PDB 格式），
  这样这个技能在零第三方依赖时也能用。

ranking_score 的含义（AlphaFold 官方）
--------------------------------------
::

    ranking_score = 0.8 * ipTM + 0.2 * pTM + 0.5 * fraction_disordered - 100 * has_clash

它把整体结构置信度（pTM、ipTM）与无序比例、原子冲突惩罚揉成一个标量，
**越高越好**。官方明确说明：这个分数**只应该用来给结构排序**，不能当
"精度"指标解读。

Examples:
    # 从 AlphaFold Server 下载目录里选最佳模型
    python select_best_model.py --input-dir folds_2026_07_03/ --output-dir best/

    # 只出排名表，不做转换
    python select_best_model.py --input-dir folds/ --rank-only

    # 强制按某个键排序
    python select_best_model.py --input-dir folds/ --rank-by ptm
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path
from typing import Optional

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "packages"))

from biopipe_core import ensure_dir, get_logger  # noqa: E402

LOGGER = get_logger("biopipe.alphafold-db")

#: 置信度 JSON 里可能出现的排序键，按优先级排列。
RANK_KEYS = ("ranking_score", "ptm", "iptm", "mean_plddt")

#: 结构文件的扩展名。
STRUCTURE_SUFFIXES = (".cif", ".pdb")

#: 置信度 JSON 的文件名模式。
CONFIDENCE_PATTERN = re.compile(r"summary_confidences[_-]?(\d+)", re.IGNORECASE)

#: 结构文件的文件名模式（用来提取排名编号）。
MODEL_PATTERN = re.compile(r"(?:model|cif)[_-]?(\d+)", re.IGNORECASE)


def build_parser() -> argparse.ArgumentParser:
    """构造命令行参数解析器。"""
    parser = argparse.ArgumentParser(
        prog="select_best_model",
        description="从 AlphaFold 预测结果中选出最佳模型并转换为 PDB",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "示例:\n"
            "  python select_best_model.py --input-dir folds/ --output-dir best/\n"
            "  python select_best_model.py --input-dir folds/ --rank-only\n"
            "  python select_best_model.py --input-dir folds/ --rank-by ptm\n"
            "\n"
            "说明:\n"
            "  · 不依赖具体命名，目录里的结构文件与置信度 JSON 会被自动配对。\n"
            "  · ranking_score 越高越好（AlphaFold 官方排序指标）。\n"
        ),
    )
    parser.add_argument(
        "--input-dir",
        required=True,
        help="AlphaFold 预测结果目录（含结构文件与 summary_confidences_*.json）",
    )
    parser.add_argument(
        "--output-dir",
        default="",
        help="最佳模型的输出目录（缺省则放在输入目录旁边）",
    )
    parser.add_argument(
        "--rank-only",
        action="store_true",
        help="只输出排名表，不做任何文件转换",
    )
    parser.add_argument(
        "--rank-by",
        default="ranking_score",
        choices=sorted(RANK_KEYS),
        help="按哪个指标排序（默认: ranking_score）",
    )
    parser.add_argument(
        "--keep-all",
        action="store_true",
        help="转换全部模型（不只最优的那个）",
    )
    parser.add_argument(
        "--no-convert",
        action="store_true",
        help="跳过 CIF → PDB 转换（只复制原始文件）",
    )
    return parser


def find_confidence_files(directory: Path) -> dict[int, Path]:
    """找出目录里的置信度 JSON，按排名编号索引。

    Args:
        directory: 预测结果目录。

    Returns:
        排名编号到文件路径的映射。
    """
    result: dict[int, Path] = {}
    for path in sorted(directory.glob("*.json")):
        if "summary_confidences" not in path.name.lower():
            continue
        match = CONFIDENCE_PATTERN.search(path.name)
        if match:
            result[int(match.group(1))] = path
    return result


def find_structure_files(directory: Path) -> dict[int, Path]:
    """找出目录里的结构文件，按排名编号索引。

    只认顶层目录——AlphaFold 的下载包里结构文件与置信度文件是平级的。

    文件名里**可能没有编号**（不同来源的命名约定差别很大）。这种情况
    下先按发现顺序占位，等 :func:`collect_models` 拿到置信度文件的编号
    后再对齐。

    Args:
        directory: 预测结果目录。

    Returns:
        排名编号到文件路径的映射。没有编号的文件用负数占位。
    """
    result: dict[int, Path] = {}
    anonymous: list[Path] = []
    for path in sorted(directory.iterdir()):
        if not path.is_file():
            continue
        suffix = path.suffix.lower()
        if suffix not in STRUCTURE_SUFFIXES:
            continue
        match = MODEL_PATTERN.search(path.stem)
        if match:
            result[int(match.group(1))] = path
        else:
            anonymous.append(path)
    # 没有编号的文件按发现顺序补在后面（负数占位）
    for index, path in enumerate(anonymous):
        result.setdefault(-1 - index, path)
    return result


def read_confidence(path: Path) -> dict:
    """读取置信度 JSON。

    Args:
        path: JSON 路径。

    Returns:
        解析后的字典；文件损坏时返回空字典。
    """
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return {}
    return data if isinstance(data, dict) else {}


def rank_value(confidence: dict, key: str) -> Optional[float]:
    """从置信度字典里取排序指标的数值。

    Args:
        confidence: 置信度字典。
        key: 指标名。

    Returns:
        数值；缺失或无法解析时返回 ``None``。
    """
    value = confidence.get(key)
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def collect_models(directory: Path, rank_by: str) -> list[dict]:
    """收集目录里所有模型及其指标。

    配对策略：

    1. 结构文件名里能提取到编号时，直接按编号配对置信度文件。
    2. 提取不到编号时（命名不规则），按**发现顺序**与置信度文件的编号
       顺序对齐——因为同一个预测包里的文件本来就是按排名顺序产出的。

    Args:
        directory: 预测结果目录。
        rank_by: 排序键。

    Returns:
        每项含 ``rank``、``structure``、``confidence_file``、``confidence``。
        已按 ``rank_by`` 降序排列。
    """
    structures = find_structure_files(directory)
    confidences = find_confidence_files(directory)

    # 编号明确的先占位
    paired: dict[int, tuple[Path, Optional[Path]]] = {}
    anonymous: list[Path] = []
    for rank, structure in structures.items():
        if rank < 0:
            anonymous.append(structure)
        else:
            paired[rank] = (structure, confidences.get(rank))

    # 编号不明确的按顺序与剩余置信度文件对齐
    remaining = sorted(
        (rank for rank in confidences if rank not in paired), reverse=True
    )
    for index, structure in enumerate(anonymous):
        confidence_rank = remaining[index] if index < len(remaining) else None
        key = confidence_rank if confidence_rank is not None else -1 - index
        paired[key] = (structure, confidences.get(confidence_rank) if confidence_rank is not None else None)

    models: list[dict] = []
    for rank, (structure, confidence_path) in paired.items():
        confidence = read_confidence(confidence_path) if confidence_path else {}
        models.append(
            {
                "rank": rank,
                "structure": structure,
                "confidence_file": confidence_path,
                "confidence": confidence,
            }
        )

    # 排序：有该指标的按数值降序在前，没有的排最后
    def sort_key(model: dict) -> tuple[int, float]:
        value = rank_value(model["confidence"], rank_by)
        if value is None:
            return (1, 0.0)
        return (0, -value)

    models.sort(key=sort_key)
    return models


def cif_to_pdb_biopython(cif_path: Path, pdb_path: Path) -> int:
    """用 Biopython 把 mmCIF 转成 PDB。

    Args:
        cif_path: 源 mmCIF 文件。
        pdb_path: 目标 PDB 文件。

    Returns:
        写出的原子数；失败返回 0。
    """
    try:
        from Bio.PDB import PDBIO, MMCIFParser  # noqa: PLC0415
    except ImportError:
        return 0

    try:
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("structure", str(cif_path))
        writer = PDBIO()
        writer.set_structure(structure)
        pdb_path.parent.mkdir(parents=True, exist_ok=True)
        writer.save(str(pdb_path))
        return sum(1 for _ in structure.get_atoms())
    except Exception:  # noqa: BLE001 - 解析器的异常类型不统一
        return 0


_CIF_ATOM_RE = re.compile(
    r"^_atom_site\.(?P<group>ATOM|HETATM)\s+"
    r"\d+\s+"
    r"(?P<atom>\S+)\s+"
    r"(?P<alt>\S+)\s+"
    r"(?P<resname>\S+)\s+"
    r"(?P<chain>\S+)\s+"
    r"(?P<resseq>\S+)\s+"
    r"(?P<icode>\S+)\s+"
    r"(?P<x>-?\d+\.?\d*)\s+"
    r"(?P<y>-?\d+\.?\d*)\s+"
    r"(?P<z>-?\d+\.?\d*)\s+"
    r"(?P<occ>-?\d+\.?\d*)\s+"
    r"(?P<bfac>-?\d+\.?\d*)",
    re.IGNORECASE,
)


def _format_pdb_atom(
    group: str,
    serial: int,
    atom: str,
    resname: str,
    chain: str,
    resseq: str,
    x: float,
    y: float,
    z: float,
    occupancy: float,
    bfactor: float,
) -> str:
    """按 PDB 固定列宽拼一行 ATOM/HETATM 记录。

    PDB 格式对列位置极其敏感。以真实文件的一行为准::

        ATOM      7  CB  CYS E   1      12.737  23.108  12.717
        ^^^^^^^^  ^^ ^^^^^^ ^^^^^ ^^ ^^^^    ^^^^^^^^ ^^^^^^^^ ^^^^^^^^
        1-6      11 13-16 18-20 21 22    31-38    39-46    47-54

    要点：

    - 序号占 6-10，索引 11 是**分隔空格**，索引 12 也是空格；
    - **原子名占索引 13-16**（左对齐，如 ``  CB``），索引 17 是空格；
    - **残基名占索引 18-20**（左对齐，如 ``CYS``），索引 21 是空格；
    - **链 ID 在索引 22**；
    - 残基号占 23-26，索引 27-30 是空格，坐标从索引 31 开始。

    注意原子名与残基名都是**左对齐**（不是右对齐）。

    手写很容易错位，所以严格按格式来。
    """
    chain_id = (chain or "A")[:1]
    # 原子名占索引 13-16（4 列，左对齐）
    atom_name = atom.strip()[:4].ljust(4)
    # 残基名占索引 18-20（3 列，左对齐）
    res_name = resname[:3].ljust(3)
    try:
        res_seq = f"{int(resseq):>4d}"
    except ValueError:
        res_seq = str(resseq)[:4].rjust(4)

    return (
        f"{group:<6}{serial:>5}  {atom_name}{res_name} {chain_id}"
        f"{res_seq}    "
        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}"
        f"{occupancy:>6.2f}{bfactor:>6.2f}          "
    )


def cif_to_pdb_fallback(cif_path: Path, pdb_path: Path) -> int:
    """不用 Biopython，把 mmCIF 的原子记录转成 PDB。

    只提取 ``_atom_site`` 的原子行——没有连接信息、没有 TER 记录，
    但对接、PyMOL 显示、几何统计这些场景够用。

    Args:
        cif_path: 源 mmCIF 文件。
        pdb_path: 目标 PDB 文件。

    Returns:
        写出的原子数；失败返回 0。
    """
    try:
        lines = cif_path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError:
        return 0

    # 找到 _atom_site 段落：它以 loop_ 开头，列名以 _atom_site. 标注
    in_loop = False
    columns: list[str] = []
    records: list[tuple[str, ...]] = []

    for line in lines:
        stripped = line.strip()
        if stripped == "loop_":
            in_loop = True
            columns = []
            continue
        if in_loop and stripped.startswith("_atom_site."):
            columns.append(stripped[len("_atom_site.") :])
            continue
        if in_loop and stripped.startswith("#"):
            break
        if in_loop and columns and stripped and not stripped.startswith("_"):
            records.append(tuple(stripped.split()))
        elif in_loop and not stripped:
            continue

    if not columns or not records:
        return 0

    index = {name: position for position, name in enumerate(columns)}

    def get(record: tuple[str, ...], name: str, default: str = "") -> str:
        position = index.get(name)
        if position is None or position >= len(record):
            return default
        return record[position]

    pdb_path.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    with pdb_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("REMARK   1 CONVERTED FROM mmCIF BY biopipe select_best_model\n")
        for serial, record in enumerate(records, 1):
            group = get(record, "group_PDB", "ATOM").upper()
            if group not in ("ATOM", "HETATM"):
                continue
            atom = get(record, "label_atom_id") or get(record, "auth_atom_id")
            if atom in (".", "?", ""):
                continue
            resname = get(record, "label_comp_id") or get(record, "auth_comp_id")
            chain = get(record, "auth_asym_id") or get(record, "label_asym_id")
            resseq = get(record, "auth_seq_id") or get(record, "label_seq_id")
            try:
                x = float(get(record, "Cartn_x", "0"))
                y = float(get(record, "Cartn_y", "0"))
                z = float(get(record, "Cartn_z", "0"))
                occ = float(get(record, "occupancy", "1.0"))
                bfac = float(get(record, "B_iso_or_equiv", "0"))
            except ValueError:
                continue
            handle.write(
                _format_pdb_atom(
                    group, serial, atom, resname, chain, resseq, x, y, z, occ, bfac
                )
                + "\n"
            )
            written += 1
        handle.write("END\n")
    return written


def convert_to_pdb(structure_path: Path, pdb_path: Path) -> tuple[bool, str]:
    """把结构文件转成 PDB（必要时先做 CIF → PDB）。

    Args:
        structure_path: 源结构文件。
        pdb_path: 目标 PDB 文件。

    Returns:
        ``(是否成功, 说明)``。
    """
    if structure_path.suffix.lower() == ".pdb":
        import shutil

        try:
            pdb_path.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(structure_path, pdb_path)
        except OSError as exc:
            return False, str(exc)
        return True, "已是 PDB，直接复制"

    atoms = cif_to_pdb_biopython(structure_path, pdb_path)
    if atoms:
        return True, f"Biopython 转换完成，{atoms} 个原子"

    LOGGER.debug("Biopython 不可用或转换失败，改用内置解析器")
    atoms = cif_to_pdb_fallback(structure_path, pdb_path)
    if atoms:
        return True, f"内置解析器转换完成，{atoms} 个原子"

    return False, "CIF 转换失败（Biopython 与内置解析器都没能解析出原子）"


def write_report(models: list[dict], rank_by: str, output_dir: Path) -> tuple[Path, Path]:
    """写出模型排名表。

    Args:
        models: 模型列表。
        rank_by: 排序键。
        output_dir: 输出目录。

    Returns:
        ``(csv 路径, json 路径)``。
    """
    fieldnames = [
        "position",
        "rank",
        "structure_file",
        "confidence_file",
        "ranking_score",
        "ptm",
        "iptm",
        "mean_plddt",
        "fraction_disordered",
        "has_clash",
    ]

    rows: list[dict] = []
    for position, model in enumerate(models, 1):
        confidence = model["confidence"]
        rows.append(
            {
                "position": position,
                "rank": model["rank"],
                "structure_file": model["structure"].name,
                "confidence_file": (
                    model["confidence_file"].name if model["confidence_file"] else ""
                ),
                "ranking_score": _fmt(rank_value(confidence, "ranking_score")),
                "ptm": _fmt(rank_value(confidence, "ptm")),
                "iptm": _fmt(rank_value(confidence, "iptm")),
                "mean_plddt": _fmt(rank_value(confidence, "mean_plddt")),
                "fraction_disordered": _fmt(
                    rank_value(confidence, "fraction_disordered")
                ),
                "has_clash": confidence.get("has_clash", ""),
            }
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = output_dir / "model_ranking.csv"
    with csv_path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    json_path = output_dir / "model_ranking.json"
    json_path.write_text(
        json.dumps(
            {
                "rank_by": rank_by,
                "models": [
                    {
                        "position": row["position"],
                        "rank": row["rank"],
                        "structure_file": row["structure_file"],
                        "confidence_file": row["confidence_file"],
                        "ranking_score": row["ranking_score"],
                        "ptm": row["ptm"],
                        "iptm": row["iptm"],
                    }
                    for row in rows
                ],
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )
    return csv_path, json_path


def _fmt(value: Optional[float]) -> str:
    """把浮点数格式化成便于阅读的字符串。"""
    if value is None:
        return ""
    return f"{value:.6f}"


def main(argv: Optional[list[str]] = None) -> int:
    """命令行主入口。"""
    args = build_parser().parse_args(argv)

    input_dir = Path(args.input_dir)
    if not input_dir.is_dir():
        LOGGER.error("输入目录不存在: %s", input_dir)
        return 2

    models = collect_models(input_dir, args.rank_by)
    if not models:
        LOGGER.error("在 %s 里没有找到任何结构文件", input_dir)
        return 1

    output_dir = (
        ensure_dir(args.output_dir)
        if args.output_dir
        else input_dir.parent / f"{input_dir.name}_selected"
    )

    csv_path, json_path = write_report(models, args.rank_by, output_dir)
    LOGGER.info("共找到 %d 个模型，排名表: %s", len(models), csv_path)

    best = models[0]
    best_score = rank_value(best["confidence"], args.rank_by)
    LOGGER.info(
        "最佳模型: %s（%s=%s）",
        best["structure"].name,
        args.rank_by,
        _fmt(best_score) if best_score is not None else "N/A",
    )

    if args.rank_only:
        return 0

    targets = models if args.keep_all else models[:1]
    succeeded = 0
    for model in targets:
        structure = model["structure"]
        destination = output_dir / f"{structure.stem}.pdb"
        if args.no_convert and structure.suffix.lower() != ".pdb":
            destination = output_dir / structure.name
            import shutil

            shutil.copy2(structure, destination)
            LOGGER.info("复制: %s", destination.name)
            succeeded += 1
            continue
        ok, detail = convert_to_pdb(structure, destination)
        if ok:
            LOGGER.info("已转换: %s —— %s", destination.name, detail)
            succeeded += 1
        else:
            LOGGER.error("转换失败: %s —— %s", structure.name, detail)

    LOGGER.info("-" * 56)
    LOGGER.info("完成：%d/%d 个模型处理成功，输出目录 %s", succeeded, len(targets), output_dir)
    return 0 if succeeded else 1


if __name__ == "__main__":
    sys.exit(main())
