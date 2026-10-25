#!/usr/bin/env python3
"""相互作用网络图生成（Mermaid）与 PyMOL 标注脚本生成。

- :func:`build_mermaid` —— 生成 Mermaid flowchart 描述的残基相互作用网络；
- :func:`build_pymol_script` —— 生成可直接在 PyMOL 中 ``run`` 的标注脚本，
  用棒状表示相互作用残基、用虚线表示相互作用连线并加标签；
- :func:`residue_statistics` —— 残基级统计（每个残基参与多少种相互作用）。

Usage:
    from interaction_network import build_mermaid, build_pymol_script
    md = build_mermaid(interactions)
    pml = build_pymol_script(interactions, structure="model.pdb")
"""

from __future__ import annotations

import re
from collections.abc import Iterable, Sequence
from typing import Optional

__all__ = [
    "build_mermaid",
    "build_pymol_script",
    "residue_statistics",
    "INTERACTION_COLORS",
    "INTERACTION_STYLES",
]

# Mermaid 边线颜色
INTERACTION_COLORS = {
    "salt_bridges": "#E74C3C",
    "hydrophobic_interactions": "#F39C12",
    "hydrogen_bonds": "#3498DB",
    "pi_cation_interactions": "#9B59B6",
    "pi_stacks": "#E91E63",
    "water_bridges": "#00BCD4",
    "halogen_bonds": "#FFC107",
    "metal_complexes": "#27AE60",
}

# PyMOL 连线样式
INTERACTION_STYLES = {
    "hydrogen_bonds": ("blue", 2.5, 0.7),
    "salt_bridges": ("yellow", 2.5, 0.8),
    "hydrophobic_interactions": ("grey", 2.0, 0.6),
    "pi_stacks": ("magenta", 2.5, 0.7),
    "pi_cation_interactions": ("violet", 2.5, 0.7),
    "water_bridges": ("cyan", 2.0, 0.5),
    "halogen_bonds": ("orange", 2.0, 0.6),
    "metal_complexes": ("green", 2.5, 0.8),
}

_SAFE = re.compile(r"[^0-9A-Za-z_]")


def _residue_key(record_part: dict) -> str:
    """把残基信息拼成 ``CHAIN:RESNAME+NUM`` 形式的稳定键。

    配体残基（无链标识）统一用 ``LIG`` 作为链标识，便于在网络图中分组。
    """
    chain = record_part.get("reschain") or record_part.get("chain") or ""
    resname = record_part.get("restype") or record_part.get("resname") or "UNK"
    resnr = record_part.get("resnr") or "?"
    if not chain:
        chain = "LIG"
    return f"{chain}:{resname}+{resnr}"


def _display(key: str) -> str:
    """把残基键转成人类可读标签（``TYR102``）。"""
    body = key.split(":", 1)[-1]
    if "+" in body:
        resname, _, resnr = body.partition("+")
        return f"{resname}{resnr}"
    return body


def _node_id(prefix: str, key: str) -> str:
    """生成 Mermaid / PyMOL 中合法的节点标识符。"""
    return f"{prefix}_{_SAFE.sub('_', key)}"


def _iter_residue_pairs(records: Iterable[dict]):
    """遍历记录，产出 (受体键, 供体键, 类型)。"""
    for record in records:
        receptor = _residue_key(record.get("receptor") or {})
        ligand = _residue_key(record.get("ligand") or {})
        if not receptor or not ligand:
            continue
        yield receptor, ligand, record.get("type", "")


def build_mermaid(records: Sequence[dict], title: str = "Residue Interaction Network") -> str:
    """生成 Mermaid 网络图。

    节点按链分组放在两个 ``subgraph`` 中（受体链 / 配体链），
    边按相互作用类型着色并带文字标签。

    Args:
        records: PLIP 解析出的相互作用记录列表。
        title: 图标题（写进图前的一行注释）。

    Returns:
        Mermaid markdown 文本（含 ```` ```mermaid ```` 围栏）。
    """
    receptor_nodes: set[str] = set()
    ligand_nodes: set[str] = set()
    edges: list[tuple[str, str, str]] = []

    for receptor, ligand, itype in _iter_residue_pairs(records):
        receptor_nodes.add(receptor)
        ligand_nodes.add(ligand)
        edges.append((receptor, ligand, itype))

    def _sort_key(key: str):
        digits = "".join(ch for ch in key.split("+")[-1] if ch.isdigit())
        return (int(digits) if digits else 0, key)

    lines: list[str] = []
    lines.append("```mermaid")
    lines.append("flowchart LR")
    lines.append(f"    %% {title}")
    lines.append("    classDef default fill:#fff,stroke:#000,stroke-width:1px,color:#000;")
    lines.append("    classDef receptor fill:#e8f0fe,stroke:#4285f4,stroke-width:1px,color:#000;")
    lines.append("    classDef ligand fill:#fef7e0,stroke:#f9ab00,stroke-width:1px,color:#000;")
    lines.append("")

    lines.append('    subgraph RECEPTOR["Receptor Residues"]')
    for node in sorted(receptor_nodes, key=_sort_key):
        lines.append(f'        {_node_id("R", node)}["{_display(node)}"]:::receptor')
    lines.append("    end")
    lines.append("")

    lines.append('    subgraph LIGAND["Ligand Residues"]')
    for node in sorted(ligand_nodes, key=_sort_key):
        lines.append(f'        {_node_id("L", node)}["{_display(node)}"]:::ligand')
    lines.append("    end")
    lines.append("")

    # 先写所有边，再统一用 linkStyle <index> 上色（Mermaid 的 linkStyle 按边序号索引）
    for receptor, ligand, itype in edges:
        label = itype.replace("_", " ").title()
        lines.append(
            f'    {_node_id("R", receptor)} -.->|{label}| {_node_id("L", ligand)}'
        )
    lines.append("")

    for index, (_receptor, _ligand, itype) in enumerate(edges):
        color = INTERACTION_COLORS.get(itype, "#95a5a6")
        lines.append(f"    linkStyle {index} stroke:{color},stroke-width:1.5px;")

    lines.append("```")
    return "\n".join(lines) + "\n"


def build_pymol_script(
    records: Sequence[dict],
    structure: str = "model.pdb",
    output_image: str = "interaction_network.png",
) -> str:
    """生成 PyMOL 标注脚本。

    脚本会：

    1. 载入结构并以 cartoon 显示；
    2. 把参与相互作用的残基显示为棒状并按类型着色；
    3. 在残基之间画虚线表示相互作用；
    4. 渲染输出图片。

    Args:
        records: 相互作用记录。
        structure: 结构文件路径（会被 ``load`` 进 PyMOL）。
        output_image: 输出图片路径。

    Returns:
        PyMOL 脚本文本（可 ``run`` 执行）。
    """
    lines: list[str] = []
    lines.append("# PyMOL 相互作用标注脚本（由 protein-interactions 技能生成）")
    lines.append(f"load {structure}, model")
    lines.append("hide everything, model")
    lines.append("show cartoon, model")
    lines.append("set cartoon_transparency, 0.75")
    lines.append("set label_font_id, 7")
    lines.append("set dash_gap, 0.35")
    lines.append("set dash_radius, 0.09")
    lines.append("set ray_opaque_background, 1")
    lines.append("bg_color white")
    lines.append("")

    # 收集需要标注的残基
    residue_selections: list[tuple[str, str, str]] = []
    for receptor, ligand, itype in _iter_residue_pairs(records):
        for key in (receptor, ligand):
            residue_selections.append((key, itype, "receptor" if key == receptor else "ligand"))

    seen: set[str] = set()
    unique: list[tuple[str, str]] = []
    for key, itype, _role in residue_selections:
        if key in seen:
            continue
        seen.add(key)
        unique.append((key, itype))

    lines.append("# --- 相互作用残基（棒状表示）---")
    for key, itype in unique:
        chain = key.split(":", 1)[0]
        body = key.split(":", 1)[-1]
        resname, _, resnr = body.partition("+")
        color = INTERACTION_STYLES.get(itype, ("grey", 2.0, 0.6))[0]
        sel = _node_id("res", key)
        lines.append(
            f"select {sel}, chain {chain} and resi {resnr} and resn {resname}"
        )
        lines.append(f"show sticks, {sel}")
        lines.append(f"color {color}, {sel}")
        lines.append(f'label {sel}, "{resname}{resnr}"')
    lines.append("")

    lines.append("# --- 相互作用连线 ---")
    edge_index = 0
    for receptor, ligand, itype in _iter_residue_pairs(records):
        edge_index += 1
        color = INTERACTION_STYLES.get(itype, ("grey", 2.0, 0.6))[0]
        rchain = receptor.split(":", 1)[0]
        rbody = receptor.split(":", 1)[-1]
        rname, _, rnum = rbody.partition("+")
        lchain = ligand.split(":", 1)[0]
        lbody = ligand.split(":", 1)[-1]
        lname, _, lnum = lbody.partition("+")
        lines.append(
            f"distance edge{edge_index}, "
            f"(chain {rchain} and resi {rnum} and resn {rname} and elem C), "
            f"(chain {lchain} and resi {lnum} and resn {lname} and elem C), "
            f"mode=2"
        )
        lines.append(f"set dash_color, {color}, edge{edge_index}")
        lines.append(f"hide labels, edge{edge_index}")
    lines.append("")

    lines.append("# --- 渲染 ---")
    lines.append(f"png {output_image}, width=2400, height=1800, dpi=300, ray=1")
    lines.append("print('相互作用网络图已保存: %s' % 'residual_interaction_network.png')")
    return "\n".join(lines) + "\n"


def residue_statistics(records: Sequence[dict]) -> dict:
    """统计残基级的相互作用参与情况。

    Returns:
        ``{"receptor": [...], "ligand": [...], "interaction_types": {...}}``，
        其中 receptor/ligand 是按参与次数降序排列的
        ``[{"key", "label", "count", "types"}]`` 列表。
    """
    receptor_counts: dict[str, dict] = {}
    ligand_counts: dict[str, dict] = {}
    type_counts: dict[str, int] = {}

    for record in records:
        itype = record.get("type", "unknown")
        type_counts[itype] = type_counts.get(itype, 0) + 1
        receptor = _residue_key(record.get("receptor") or {})
        ligand = _residue_key(record.get("ligand") or {})
        for key, bucket in ((receptor, receptor_counts), (ligand, ligand_counts)):
            if not key:
                continue
            entry = bucket.setdefault(key, {"count": 0, "types": set()})
            entry["count"] += 1
            entry["types"].add(itype)

    def _materialise(bucket: dict) -> list[dict]:
        items = [
            {
                "key": key,
                "label": _display(key),
                "count": entry["count"],
                "types": sorted(entry["types"]),
            }
            for key, entry in bucket.items()
        ]
        items.sort(key=lambda item: (-item["count"], item["label"]))
        return items

    return {
        "receptor": _materialise(receptor_counts),
        "ligand": _materialise(ligand_counts),
        "interaction_types": type_counts,
    }


def _main(argv: Optional[list[str]] = None) -> int:
    """命令行入口：解析 PLIP 报告并输出网络图。"""
    import argparse
    from pathlib import Path

    parser = argparse.ArgumentParser(
        prog="interaction_network", description="相互作用网络图生成（命令行自检工具）"
    )
    parser.add_argument("input", help="PLIP report.xml 路径")
    parser.add_argument("--mermaid", metavar="PATH", help="写出 Mermaid 图")
    parser.add_argument("--pymol", metavar="PATH", help="写出 PyMOL 脚本")
    args = parser.parse_args(argv)

    from plip_parser import parse_plip_xml

    report = parse_plip_xml(args.input)
    records = report["interactions"]

    stats = residue_statistics(records)
    logger.info("总相互作用: %d", report["summary"]["total"])
    for key, value in stats["interaction_types"].items():
        logger.info("  %-26s %d", key, value)

    if args.mermaid:
        Path(args.mermaid).write_text(build_mermaid(records), encoding="utf-8")
        logger.info("已写出 %s", args.mermaid)
    if args.pymol:
        Path(args.pymol).write_text(
            build_pymol_script(records, structure=report["source"]), encoding="utf-8"
        )
        logger.info("已写出 %s", args.pymol)
    return 0


if __name__ == "__main__":
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "packages"))
    from biopipe_core import get_logger

    logger = get_logger("biopipe.protein-interactions")
    sys.exit(_main())
