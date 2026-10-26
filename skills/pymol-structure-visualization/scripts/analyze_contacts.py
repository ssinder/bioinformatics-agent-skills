#!/usr/bin/env python3
"""链间残基接触分析 —— 主入口。

从 PDB 文件算出链间残基接触，按化学类型分类，输出 CSV / JSON /
Mermaid 网络图 / PyMOL 可视化脚本。

**零依赖**：核心计算只用 Python 标准库，不需要 PyMOL、不需要 numpy、
不需要 PLIP。PyMOL 只在需要渲染脚本时才作为可选依赖出现。

用法::

    # 分析两条链之间的接触
    python analyze_contacts.py --input complex.pdb --chains A B --output out/

    # 所有链两两分析 + 生成 PyMOL 脚本
    python analyze_contacts.py --input complex.pdb --output out/ --make-pymol

    # 只要接触类型统计
    python analyze_contacts.py --input complex.pdb --chains A B --json out/contacts.json

从老仓库 ``pymol_show/`` 移植时修掉的问题
------------------------------------------
1. **必须启动 PyMOL GUI** —— 老脚本全部以 ``pymol.finish_launching(['pymol'])``
   开头，且在最后 ``cmd.zoom()``。这意味着：无法 headless 运行、
   无 PyMOL 的机器上完全跑不了、CI 里无法测试。本技能把几何计算
   与 PyMOL 渲染彻底分离 —— 分析部分纯 Python，PyMOL 脚本是**生成物**
   而不是运行依赖。
2. **硬编码个人路径** —— 老脚本写死 ``C:\\Users\\DELL\\Desktop\\1.pdb``。
   现在全部走命令行参数。
3. **分类表互相矛盾** —— ``donors``/``acceptors`` 在不同脚本里不一样，
   ASP/GLU 在一处被当供体。已统一到 :mod:`contact_chemistry`。
4. **文件名带空格** —— 老脚本名如 ``interaction count.py``、``hbond count.py``，
   在 shell 里必须加引号，且 PyMOL 的 ``run`` 命令对带空格路径的处理
   也容易出问题。现在文件名全部是合法的 Python 模块名。
5. **PLIP 结果汇总靠固定路径** —— 老 ``use plip find interact region.py``
   读 ``C:\\Users\\DELL\\Desktop\\output\\*\\interaction_residues.txt``。
   现在支持目录扫描，且 PLIP 报告解析由既有的 ``protein-interactions``
   技能负责（避免重复实现）。
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "packages"))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "protein-structure" / "scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from typing import Optional  # noqa: E402

from biopipe_core import cli_parser, get_logger, run_cli  # noqa: E402
from biopipe_core.cli import SkillError  # noqa: E402
from contact_chemistry import (  # noqa: E402
    INTERACTION_TYPES,
    LIGAND_CONTACT,
    classify_pair_multi,
    describe_rules,
)
from contact_geometry import (  # noqa: E402
    DEFAULT_MAX_DISTANCE,
    AtomRef,
    ContactDetector,
    ResiduePair,
)

LOGGER = get_logger("biopipe.pymol-structure-visualization")

SKILL_NAME = "pymol-structure-visualization"


# ---------------------------------------------------------------- PDB 读取


def _load_atoms(pdb_path: Path) -> tuple[list[AtomRef], dict]:
    """读取 PDB 文件，返回原子列表与结构概要。

    优先用 ``protein-structure`` 技能自带的 PDB 解析器（零依赖）；
    不可用时退回本模块的精简解析。
    """
    try:
        from pdb_parser import parse_pdb_file

        structure = parse_pdb_file(str(pdb_path))
        atoms: list[AtomRef] = []
        for atom in structure.atoms():
            atoms.append(AtomRef(
                serial=atom.serial,
                name=atom.name,
                res_name=atom.res_name,
                chain_id=atom.chain_id,
                res_seq=atom.res_seq,
                icode=atom.icode,
                x=atom.coord[0],
                y=atom.coord[1],
                z=atom.coord[2],
                element=atom.element,
                hetero=atom.hetatm,
            ))
        summary = structure.summary()
        summary["source"] = "protein-structure/pdb_parser"
        return atoms, summary

    except ImportError:
        LOGGER.debug("protein-structure 的 pdb_parser 不可用，使用内置精简解析")
        return _load_atoms_fallback(pdb_path)


def _load_atoms_fallback(pdb_path: Path) -> tuple[list[AtomRef], dict]:
    """内置的精简 PDB 解析器（只取坐标，不做完整校验）。"""
    atoms: list[AtomRef] = []
    chains: set[str] = set()
    residues: set[tuple[str, int, str]] = set()
    n_models = 0

    with pdb_path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            record = line[:6]
            if record == "MODEL ":
                n_models += 1
                continue
            if record == "ENDMDL":
                # 只用第一个 MODEL
                if n_models >= 1:
                    break
                continue
            if record not in ("ATOM  ", "HETATM"):
                continue

            try:
                serial = int(line[6:11].strip() or 0)
            except ValueError:
                serial = 0
            name = line[12:16].strip()
            res_name = line[17:20].strip()
            chain_id = line[21:22].strip() or "_"
            try:
                res_seq = int(line[22:26].strip() or 0)
            except ValueError:
                res_seq = 0
            icode = line[26:27].strip()
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            element = line[76:78].strip()

            chains.add(chain_id)
            residues.add((chain_id, res_seq, icode))
            atoms.append(AtomRef(
                serial=serial,
                name=name,
                res_name=res_name,
                chain_id=chain_id,
                res_seq=res_seq,
                icode=icode,
                x=x, y=y, z=z,
                element=element,
                hetero=(record == "HETATM"),
            ))

    summary = {
        "structure_id": pdb_path.stem,
        "n_models": max(n_models, 1),
        "n_chains": len(chains),
        "n_residues": len(residues),
        "n_atoms": len(atoms),
        "chain_ids": sorted(chains),
        "source": "builtin-fallback",
    }
    return atoms, summary


# ---------------------------------------------------------------- 分析


def _pair_record(
    pair: ResiduePair,
    include_histidine: bool,
    include_ambiguous_donors: bool,
) -> dict:
    """把一个残基对转成记录字典（可能含多种相互作用类型）。"""
    types = classify_pair_multi(
        pair.res_name_a, pair.res_name_b, pair.min_distance,
        include_histidine=include_histidine,
        include_ambiguous_donors=include_ambiguous_donors,
    )
    return {
        "chain_a": pair.chain_a,
        "res_name_a": pair.res_name_a,
        "res_seq_a": pair.res_seq_a,
        "label_a": pair.label_a,
        "chain_b": pair.chain_b,
        "res_name_b": pair.res_name_b,
        "res_seq_b": pair.res_seq_b,
        "label_b": pair.label_b,
        "min_distance": round(pair.min_distance, 3),
        "n_atom_pairs": pair.n_atom_pairs,
        "interaction_types": types,
        "interaction_type": types[0] if types else "",
    }


def _res_sort(label: str) -> tuple:
    """残基标签排序键（链 + 序号）。"""
    chain = label.split(":")[0]
    rest = label.split(":")[1] if ":" in label else label
    digits = "".join(ch for ch in rest if ch.isdigit())
    return (chain, int(digits or 0))


#: 配体接触类型的显示名与配色（单独处理，因为它不属于残基级化学类型）。
LIGAND_STYLE = {
    "name": LIGAND_CONTACT,
    "label": "配体接触",
    "color": "#7F7F7F",
    "pymol_color": "grey",
}


def _count_by_type(records: list[dict]) -> dict[str, int]:
    """按相互作用类型统计接触数量。

    一个残基对可能有多个类型（例如盐桥 + 氢键），分别计入各自类型。
    配体接触（``Ligand_contact``）单独占一个计数位。
    """
    by_type: dict[str, int] = {rule.name: 0 for rule in INTERACTION_TYPES}
    by_type[LIGAND_CONTACT] = 0
    for record in records:
        for type_name in record["interaction_types"]:
            by_type[type_name] = by_type.get(type_name, 0) + 1
    return by_type


def _type_label(type_name: str) -> str:
    """类型的中文显示名。"""
    for rule in INTERACTION_TYPES:
        if rule.name == type_name:
            return rule.label
    if type_name == LIGAND_CONTACT:
        return LIGAND_STYLE["label"]
    return type_name


def _analyze(
    atoms: list[AtomRef],
    chain_pairs: list[tuple[str, str]],
    max_distance: float,
    include_hetero: bool,
    ignore_hydrogen: bool,
    include_histidine: bool,
    include_ambiguous_donors: bool,
) -> dict:
    """执行接触分析，返回完整结果字典。"""
    detector = ContactDetector(
        max_distance=max_distance,
        include_hetero=include_hetero,
        ignore_hydrogen=ignore_hydrogen,
    )
    detector.add_atoms(atoms)

    all_records: list[dict] = []
    per_pair: dict[str, int] = {}

    for chain_a, chain_b in chain_pairs:
        pairs = detector.detect(chain_a, chain_b)
        per_pair[f"{chain_a}-{chain_b}"] = len(pairs)
        for pair in pairs:
            record = _pair_record(pair, include_histidine, include_ambiguous_donors)
            record["chain_pair"] = f"{chain_a}-{chain_b}"
            all_records.append(record)

    # 按类型统计（一个残基对可能有多个类型，分别计入）
    by_type = _count_by_type(all_records)

    # 热点残基（参与相互作用最多）
    hotspot_a: dict[str, int] = {}
    hotspot_b: dict[str, int] = {}
    for record in all_records:
        if not record["interaction_types"]:
            continue
        hotspot_a[record["label_a"]] = hotspot_a.get(record["label_a"], 0) + 1
        hotspot_b[record["label_b"]] = hotspot_b.get(record["label_b"], 0) + 1

    return {
        "contacts": all_records,
        "summary": {
            "n_contacts": len(all_records),
            "n_typed": sum(1 for r in all_records if r["interaction_types"]),
            "n_untyped": sum(1 for r in all_records if not r["interaction_types"]),
            "by_type": by_type,
            "per_chain_pair": per_pair,
            "hotspot_a": dict(sorted(hotspot_a.items(), key=lambda kv: (-kv[1], _res_sort(kv[0])))),
            "hotspot_b": dict(sorted(hotspot_b.items(), key=lambda kv: (-kv[1], _res_sort(kv[0])))),
        },
        "rules": describe_rules(),
    }


# ---------------------------------------------------------------- 输出


def _write_csv(records: list[dict], path: Path) -> Path:
    """写出接触明细 CSV。"""
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "Chain_Pair", "ChainA", "ResnA", "ResiA", "ChainB", "ResnB", "ResiB",
        "Min_Distance", "N_Atom_Pairs", "Interaction_Type", "All_Types",
    ]
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(columns)
        for record in records:
            writer.writerow([
                record.get("chain_pair", ""),
                record["chain_a"], record["res_name_a"], record["res_seq_a"],
                record["chain_b"], record["res_name_b"], record["res_seq_b"],
                f"{record['min_distance']:.3f}",
                record["n_atom_pairs"],
                record["interaction_type"],
                "; ".join(record["interaction_types"]),
            ])
    return path


def _write_summary_csv(result: dict, path: Path) -> Path:
    """写出按类型的统计 CSV。"""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Interaction_Type", "Label_CN", "Count"])
        for name, count in result["summary"]["by_type"].items():
            writer.writerow([name, _type_label(name), count])
    return path


#: Mermaid 节点配色（按主要类型）。
MERMAID_COLORS = {
    "Hydrogen_bond": "#4C78A8",
    "Salt_bridge": "#F58518",
    "Hydrophobic": "#54A24B",
    "Pi_Pi": "#B279A2",
    "Cation_Pi": "#E45756",
    LIGAND_CONTACT: "#7F7F7F",
}

#: PyMOL 配色（按主要类型）。
PYMOL_COLORS = {
    "Hydrogen_bond": "marine",
    "Salt_bridge": "orange",
    "Hydrophobic": "forest",
    "Pi_Pi": "purple",
    "Cation_Pi": "red",
    LIGAND_CONTACT: "grey",
}


def _write_mermaid(result: dict, path: Path, max_nodes: int = 30) -> Path:
    """生成 Mermaid 网络图。"""
    path.parent.mkdir(parents=True, exist_ok=True)

    colors = MERMAID_COLORS

    contacts = [r for r in result["contacts"] if r["interaction_types"]]
    # 优先展示多类型与近距离的接触
    contacts.sort(key=lambda r: (-len(r["interaction_types"]), r["min_distance"]))
    contacts = contacts[:max_nodes]

    lines = ["graph LR"]
    lines.append("    %% 链间残基接触网络（由 pymol-structure-visualization 生成）")

    for record in contacts:
        label_a = record["label_a"].replace(":", "_")
        label_b = record["label_b"].replace(":", "_")
        primary = record["interaction_types"][0]
        color = colors.get(primary, "#999999")
        lines.append(
            f'    {label_a}["{record["label_a"]}"]'
            f' --> |"{record["min_distance"]:.1f}Å {primary}"|'
            f' {label_b}["{record["label_b"]}"]'
        )
        lines.append(f"    style {label_a} fill:{color},stroke:#333,color:#fff")
        lines.append(f"    style {label_b} fill:{color},stroke:#333,color:#fff")

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _write_pymol(result: dict, path: Path) -> Path:
    """生成 PyMOL 交互式可视化脚本。

    这是**产物**而不是依赖 —— 用户在有 PyMOL 的机器上执行它即可。
    """
    path.parent.mkdir(parents=True, exist_ok=True)

    contacts = [r for r in result["contacts"] if r["interaction_types"]]
    by_type: dict[str, list[dict]] = {}
    for record in contacts:
        by_type.setdefault(record["interaction_types"][0], []).append(record)

    colors = PYMOL_COLORS

    lines = [
        "# 由 pymol-structure-visualization 生成 —— PyMOL 交互式可视化脚本",
        "#",
        "# 运行方式：",
        "#   pymol -r this_script.pml",
        "#   或在已打开的 PyMOL 里：run this_script.pml",
        "",
        "from pymol import cmd",
        "",
        "# 载入结构（按需改成你的路径）",
        "# cmd.load('complex.pdb', 'complex')",
        "",
        "# 选中并标注接触残基",
        "cmd.set('dash_gap', 0.3)",
        "cmd.set('dash_radius', 0.15)",
        "cmd.set('label_size', 14)",
        "cmd.set('label_color', black)",
        "cmd.set('label_outline_color', white)",
        "",
    ]

    for type_name, records in by_type.items():
        color = colors.get(type_name, "grey")
        labels = [r["label_a"] for r in records] + [r["label_b"] for r in records]
        unique = sorted(set(labels))
        selection = " or ".join(unique)
        lines.append(f"# --- {type_name}（{len(records)} 对）---")
        lines.append(f"cmd.select('{type_name}_res', '({selection})')")
        lines.append(f"cmd.color('{color}', '{type_name}_res')")
        lines.append(f"cmd.show('sticks', '{type_name}_res')")
        lines.append(f"cmd.label('{type_name}_res', '\"%s%s\"')")
        lines.append("")

    lines.extend([
        "# 显示所有接触残基",
        "cmd.show('sticks', 'all and not polymer.hydrophobic')",
        "cmd.zoom('all')",
        "",
        "# 导出图片（需要时取消注释）",
        "# cmd.png('contacts.png', width=1600, height=1200, dpi=300, ray=1)",
    ])

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _write_pymol_contacts(result: dict, path: Path) -> Path:
    """生成基于 ``distance`` 的接触绘制脚本（复刻老仓库的 dash 风格）。"""
    path.parent.mkdir(parents=True, exist_ok=True)

    lines = [
        "# 链间接触绘制脚本（pymol-structure-visualization）",
        "#",
        "# 用 distance 对象画出原子间的虚线，并按类型着色。",
        "",
        "from pymol import cmd",
        "",
        "cmd.set('dash_gap', 0.25)",
        "cmd.set('dash_radius', 0.12)",
        "",
    ]

    colors = PYMOL_COLORS

    for index, record in enumerate(result["contacts"]):
        if not record["interaction_types"]:
            continue
        primary = record["interaction_types"][0]
        color = colors.get(primary, "grey")
        name = f"contact_{index}"
        sel_a = f"chain {record['chain_a']} and resi {record['res_seq_a']}"
        sel_b = f"chain {record['chain_b']} and resi {record['res_seq_b']}"
        lines.append(f"cmd.distance('{name}', '{sel_a}', '{sel_b}', 6.0, mode=0)")
        lines.append(f"cmd.color('{color}', '{name}')")
        lines.append(f"cmd.set('dash_color', '{color}', '{name}')")

    lines.extend([
        "",
        "cmd.show('sticks', 'all and sidechain and not hydro')",
        "cmd.zoom('all')",
    ])

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


# ---------------------------------------------------------------- 主入口


def _build_parser() -> argparse.ArgumentParser:
    """构造命令行解析器。"""
    parser = cli_parser(
        prog=SKILL_NAME,
        description="从 PDB 文件分析链间残基接触，输出 CSV/JSON/Mermaid/PyMOL 脚本",
        skill=SKILL_NAME,
    )
    parser.add_argument("--input", required=True, metavar="PDB", help="PDB 结构文件")
    parser.add_argument("--output", metavar="DIR", help="输出目录（默认写到输入文件同级）")
    parser.add_argument("--chains", nargs="+", metavar="ID", help="要分析的链对，如 A B（默认全部链两两）")
    parser.add_argument("--max-distance", type=float, default=DEFAULT_MAX_DISTANCE,
                        help="最大接触距离（Å，默认 6.0）")
    parser.add_argument("--exclude-hetero", action="store_true",
                        help="忽略 HETATM（水、离子、配体）")
    parser.add_argument("--include-hydrogen", action="store_true",
                        help="把氢原子也算进去（默认忽略）")
    parser.add_argument("--include-histidine", action="store_true",
                        help="把 HIS 当作正电残基（影响盐桥与阳离子-π 计数）")
    parser.add_argument("--include-ambiguous-donors", action="store_true",
                        help="把 ASP/GLU 也当作弱氢键供体")
    parser.add_argument("--min-distance", type=float, default=0.0,
                        help="只报告距离不小于此值的接触（Å）")
    parser.add_argument("--make-mermaid", action="store_true", help="生成 Mermaid 网络图")
    parser.add_argument("--make-pymol", action="store_true", help="生成 PyMOL 可视化脚本")
    parser.add_argument("--make-pymol-contacts", action="store_true",
                        help="生成基于 distance 的接触绘制脚本")
    parser.add_argument("--json", metavar="PATH", help="把完整结果写成 JSON")
    return parser


def main(args: argparse.Namespace) -> int:
    """主入口。"""
    pdb_path = Path(args.input)
    if not pdb_path.exists():
        raise SkillError(f"输入文件不存在: {pdb_path}")

    output_dir = Path(args.output) if args.output else pdb_path.parent / f"{pdb_path.stem}_contacts"
    output_dir.mkdir(parents=True, exist_ok=True)

    LOGGER.info("读取结构: %s", pdb_path)
    atoms, summary = _load_atoms(pdb_path)
    LOGGER.info(
        "结构概要: %d 链 / %d 残基 / %d 原子（解析器: %s）",
        summary.get("n_chains", 0), summary.get("n_residues", 0),
        summary.get("n_atoms", 0), summary.get("source", "?"),
    )

    available_chains = sorted({atom.chain_id for atom in atoms})
    if not available_chains:
        raise SkillError(f"文件里没有找到任何原子: {pdb_path}")

    # 决定要分析的链对
    if args.chains:
        if len(args.chains) != 2:
            raise SkillError("--chains 需要恰好两个链标识（如 A B）")
        chain_a, chain_b = args.chains
        missing = [c for c in (chain_a, chain_b) if c not in available_chains]
        if missing:
            raise SkillError(f"结构里没有链 {missing}（可用: {available_chains}）")
        chain_pairs = [(chain_a, chain_b)]
    else:
        chain_pairs = [
            (available_chains[i], available_chains[j])
            for i in range(len(available_chains))
            for j in range(i + 1, len(available_chains))
        ]

    LOGGER.info("分析链对: %s", ", ".join(f"{a}-{b}" for a, b in chain_pairs))

    result = _analyze(
        atoms=atoms,
        chain_pairs=chain_pairs,
        max_distance=args.max_distance,
        include_hetero=not args.exclude_hetero,
        ignore_hydrogen=not args.include_hydrogen,
        include_histidine=args.include_histidine,
        include_ambiguous_donors=args.include_ambiguous_donors,
    )

    # 距离下限过滤
    if args.min_distance > 0:
        result["contacts"] = [
            r for r in result["contacts"] if r["min_distance"] >= args.min_distance
        ]
        result["summary"]["by_type"] = _count_by_type(result["contacts"])
        result["summary"]["n_contacts"] = len(result["contacts"])
        result["summary"]["n_typed"] = sum(1 for r in result["contacts"] if r["interaction_types"])
        result["summary"]["n_untyped"] = sum(1 for r in result["contacts"] if not r["interaction_types"])

    summary_out = result["summary"]
    LOGGER.info(
        "接触对 %d 个（已分类 %d，未分类 %d）",
        summary_out["n_contacts"], summary_out["n_typed"], summary_out["n_untyped"],
    )
    for type_name, count in summary_out["by_type"].items():
        LOGGER.info("  %-16s %d", type_name, count)

    # 写输出
    written: list[str] = []

    if result["contacts"]:
        csv_path = _write_csv(result["contacts"], output_dir / "contacts.csv")
        written.append(str(csv_path))
        LOGGER.info("接触明细 → %s", csv_path)

    summary_csv = _write_summary_csv(result, output_dir / "contact_summary.csv")
    written.append(str(summary_csv))

    if args.make_mermaid:
        path = _write_mermaid(result, output_dir / "contact_network.mmd")
        written.append(str(path))
        LOGGER.info("Mermaid 网络图 → %s", path)

    if args.make_pymol:
        path = _write_pymol(result, output_dir / "visualize_contacts.pml")
        written.append(str(path))
        LOGGER.info("PyMOL 标注脚本 → %s", path)

    if args.make_pymol_contacts:
        path = _write_pymol_contacts(result, output_dir / "draw_contacts.pml")
        written.append(str(path))
        LOGGER.info("PyMOL 接触绘制脚本 → %s", path)

    json_path = Path(args.json) if args.json else output_dir / "contacts.json"
    json_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    written.append(str(json_path))
    LOGGER.info("完整结果 → %s", json_path)

    LOGGER.info("共写出 %d 个文件", len(written))
    return 0


if __name__ == "__main__":
    sys.exit(run_cli(main, parser=_build_parser(), skill=SKILL_NAME, version="1.0.0"))
