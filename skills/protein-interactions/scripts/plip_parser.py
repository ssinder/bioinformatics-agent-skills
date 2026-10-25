#!/usr/bin/env python3
"""PLIP 报告（XML）解析器 —— 纯标准库实现。

PLIP（Protein-Ligand Interaction Profiler）在 `--xml` 模式下会输出一个
描述全部相互作用的 XML 报告。本模块用 ``xml.etree.ElementTree`` 解析该报告，
把嵌套的 site / pocket / ligand / residue 结构抽取为扁平的相互作用记录。

PLIP XML 的典型结构::

    <report>
      <addhydrogen>
        <add hydrogen="1" atom="15" type="polar" charge_group="0"/>
      </addhydrogen>
      <interactions>
        <hydrogen_bonds>
          <hydrogen_bond id="1" onf="A:ARG102:NE" onp="A:TYR93:OH"
                         distance="2.9" distance_h_a="1.9" angle="163.5"
                         donor_type="Protein" acceptor_type="Protein" d_h_a="H"/>
          ...
        </hydrogen_bonds>
        <salt_bridges>
          <salt_bridge id="1" onf="A:ARG102:NE" onp="A:ASP23:OD1" distance="2.8" .../>
        </salt_bridges>
        <hydrophobic_interactions>
          <hydrophobic_interaction id="1" onf="A:LEU5" onp="LIG:1" distance="3.9" .../>
        </hydrophobic_interactions>
        <pi_stacks>
          <pi_stack id="1" onf="A:TYR93" onp="LIG:1" distance="3.5" type="parallel" .../>
        </pi_stacks>
        <pi_cation_interactions>...</pi_cation_interactions>
        <water_bridges>...</water_bridges>
        <halogen_bonds>...</halogen_bonds>
        <metal_complexes>...</metal_complexes>
      </interactions>
    </report>

PLIP 也有 **旧版"表格"式 XML**（`05_extract_residues.py` 处理的那种），
其残基信息放在子元素的 ``<resnr>`` / ``<restype>`` / ``<reschain>`` 节点里。
本模块**同时支持两种格式**：

- ``new``：新版属性式（``onf`` / ``onp`` 编码 "chain:resname+resnr:atom"）
- ``legacy``：旧版子元素式（``resnr`` / ``restype`` / ``reschain`` +
  ``resnr_lig`` / ``restype_lig`` / ``reschain_lig``）

Usage:
    from plip_parser import parse_plip_xml
    report = parse_plip_xml("report.xml")
    print(report["summary"])
"""

from __future__ import annotations

import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any, Optional, Union

__all__ = [
    "PLIP_INTERACTION_TYPES",
    "INTERACTION_LABELS",
    "parse_plip_xml",
    "parse_residue_label",
    "summarize_interactions",
]

# PLIP 报告中的相互作用类型（按 XML 标签名）
PLIP_INTERACTION_TYPES = (
    "salt_bridges",
    "hydrophobic_interactions",
    "hydrogen_bonds",
    "pi_cation_interactions",
    "pi_stacks",
    "water_bridges",
    "halogen_bonds",
    "metal_complexes",
)

# 标签 -> 人类可读名称
INTERACTION_LABELS = {
    "salt_bridges": "Salt Bridge",
    "hydrophobic_interactions": "Hydrophobic Interaction",
    "hydrogen_bonds": "Hydrogen Bond",
    "pi_cation_interactions": "Pi-Cation Interaction",
    "pi_stacks": "Pi-Stack Interaction",
    "water_bridges": "Water Bridge",
    "halogen_bonds": "Halogen Bond",
    "metal_complexes": "Metal Complex",
}

# 新版 XML 中每个类型下承载单条相互作用的子标签名
_ITEM_TAGS = {
    "salt_bridges": "salt_bridge",
    "hydrophobic_interactions": "hydrophobic_interaction",
    "hydrogen_bonds": "hydrogen_bond",
    "pi_cation_interactions": "pi_cation_interaction",
    "pi_stacks": "pi_stack",
    "water_bridges": "water_bridge",
    "halogen_bonds": "halogen_bond",
    "metal_complexes": "metal_complex",
}


def parse_residue_label(label: Optional[str]) -> dict:
    """解析 PLIP 的 ``onf`` / ``onp`` 残基标签。

    PLIP 使用两种标签格式：

    - **蛋白残基** ``chain:RESNAME+resnr:atom``，例如 ``A:ARG+102:NE``；
    - **配体 / 杂原子残基** ``RESNAME:resnr:atom``，例如 ``LIG:1:CL1``
      （配体没有链标识与加号分隔的残基名，残基号直接跟在冒号后）。

    Args:
        label: PLIP 原始标签字符串。

    Returns:
        ``{"chain", "resname", "resnr", "atom", "raw", "is_ligand"}``；
        无法解析时相应字段为 None。
    """
    result: dict[str, Any] = {
        "chain": None,
        "resname": None,
        "resnr": None,
        "atom": None,
        "is_ligand": False,
        "raw": label or "",
    }
    if not label:
        return result

    text = label.strip()
    # 去掉可能存在的距离后缀（如 "A:ARG+102:NE (3.2)"）
    if " (" in text:
        text = text.split(" (", 1)[0].strip()

    parts = [p for p in text.split(":")]
    # 去掉尾部空段（"A:ARG+102:" 形式）
    while parts and not parts[-1].strip():
        parts.pop()

    if not parts:
        return result

    head = parts[0].strip()
    middle = parts[1].strip() if len(parts) >= 2 else ""
    atom = parts[2].strip() if len(parts) >= 3 else None

    if "+" in middle:
        # 蛋白格式：A:ARG+102:NE
        result["chain"] = head or None
        resname, _, resnr = middle.partition("+")
        result["resname"] = resname.strip() or None
        result["resnr"] = resnr.strip() or None
    elif "+" in head:
        # 形如 A+ARG:102 的变体
        chain, _, resname = head.partition("+")
        result["chain"] = chain.strip() or None
        result["resname"] = resname.strip() or None
        result["resnr"] = middle or None
    else:
        # 配体格式：LIG:1[:CL1]
        result["is_ligand"] = True
        result["resname"] = head or None
        result["resnr"] = middle or None

    if atom:
        result["atom"] = atom

    return result


def _legacy_item_to_record(tag: str, element: ET.Element) -> dict:
    """把旧版 XML 的子元素式相互作用转成统一记录。"""
    def text_of(name: str) -> str:
        node = element.find(name)
        return (node.text or "").strip() if node is not None and node.text else ""

    receptor = {
        "resnr": text_of("resnr"),
        "restype": text_of("restype"),
        "reschain": text_of("reschain"),
        "atom": text_of("atomid") or text_of("atom"),
    }
    ligand = {
        "resnr": text_of("resnr_lig"),
        "restype": text_of("restype_lig"),
        "reschain": text_of("reschain_lig"),
        "atom": text_of("atomid_lig") or text_of("atom_lig"),
    }
    return {
        "type": tag,
        "type_label": INTERACTION_LABELS.get(tag, tag),
        "receptor": receptor,
        "ligand": ligand,
        "distance": text_of("distance"),
        "angle": text_of("angle"),
        "extra": {},
    }


def _new_item_to_record(tag: str, element: ET.Element) -> dict:
    """把新版 XML 的属性式相互作用转成统一记录。"""
    attrib = element.attrib
    receptor = parse_residue_label(attrib.get("onf"))
    ligand = parse_residue_label(attrib.get("onp"))

    extra = {
        key: value
        for key, value in attrib.items()
        if key
        not in (
            "id",
            "onf",
            "onp",
            "distance",
            "distance_h_a",
            "distance_h_b",
            "angle",
            "angle_h_a",
            "angle_h_b",
            "donor_type",
            "acceptor_type",
            "d_h_a",
            "protdist",
        )
    }
    return {
        "type": tag,
        "type_label": INTERACTION_LABELS.get(tag, tag),
        "receptor": receptor,
        "ligand": ligand,
        "distance": attrib.get("distance", ""),
        "angle": attrib.get("angle", ""),
        "distance_h_a": attrib.get("distance_h_a", ""),
        "donor_type": attrib.get("donor_type", ""),
        "acceptor_type": attrib.get("acceptor_type", ""),
        "extra": extra,
    }


def parse_plip_xml(path: Union[str, Path]) -> dict:
    """解析 PLIP XML 报告。

    Args:
        path: PLIP 生成的 ``report.xml`` 路径。

    Returns:
        包含以下键的字典：

        - ``interactions``：扁平的相互作用记录列表
        - ``by_type``：按类型分组的记录字典
        - ``summary``：各类型计数与总数
        - ``source``：源文件路径
        - ``format``：``new`` 或 ``legacy``

    Raises:
        FileNotFoundError: 报告文件不存在。
        ValueError: XML 无法解析。
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"PLIP 报告不存在: {p}")

    try:
        tree = ET.parse(p)
    except ET.ParseError as exc:
        raise ValueError(f"PLIP XML 解析失败: {p}: {exc}") from exc
    root = tree.getroot()

    records: list[dict] = []

    # 优先定位 <interactions> 容器；新版报告都有该节点
    container = root.find("interactions")
    if container is None:
        container = root

    fmt = "new"
    for tag in PLIP_INTERACTION_TYPES:
        group = container.find(tag)
        if group is None:
            continue
        item_tag = _ITEM_TAGS.get(tag)
        items = group.findall(item_tag) if item_tag else []
        if items:
            for element in items:
                records.append(_new_item_to_record(tag, element))
        else:
            # 旧版：直接子节点没有统一标签名，按属性判定
            for element in group:
                if not isinstance(element.tag, str):
                    continue
                record = _legacy_item_to_record(tag, element)
                if any(record["receptor"].values()) or any(record["ligand"].values()):
                    records.append(record)
                    fmt = "legacy"

    by_type: dict[str, list[dict]] = {tag: [] for tag in PLIP_INTERACTION_TYPES}
    for record in records:
        by_type.setdefault(record["type"], []).append(record)

    summary = summarize_interactions(records)
    return {
        "interactions": records,
        "by_type": by_type,
        "summary": summary,
        "source": str(p),
        "format": fmt,
    }


def summarize_interactions(records: list[dict]) -> dict:
    """统计相互作用数量。

    Returns:
        含各类型计数与 ``total`` 的字典。
    """
    summary = {tag: 0 for tag in PLIP_INTERACTION_TYPES}
    for record in records:
        summary[record["type"]] = summary.get(record["type"], 0) + 1
    summary["total"] = len(records)
    return summary


def _main(argv: Optional[list[str]] = None) -> int:
    """命令行入口：解析 PLIP 报告并打印摘要。"""
    import argparse
    import json

    parser = argparse.ArgumentParser(
        prog="plip_parser", description="PLIP XML 报告解析器（命令行自检工具）"
    )
    parser.add_argument("input", help="PLIP report.xml 路径")
    parser.add_argument("--json", action="store_true", help="以 JSON 输出完整记录")
    args = parser.parse_args(argv)

    report = parse_plip_xml(args.input)
    logger.info("格式: %s，总相互作用数: %d", report["format"], report["summary"]["total"])
    for tag, count in report["summary"].items():
        if tag == "total":
            continue
        logger.info("  %-26s %d", tag, count)
    if args.json:
        logger.info("%s", json.dumps(report["interactions"], ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "packages"))
    from biopipe_core import get_logger

    logger = get_logger("biopipe.protein-interactions")
    sys.exit(_main())
