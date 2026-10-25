#!/usr/bin/env python3
"""蛋白质相互作用分析入口。

流程：

1. **可选** —— 从 AlphaFold DB 的 zip 包中按置信度选出最佳模型
   （:mod:`alphafold_models`，支持 pTM / ipTM / 平均 pLDDT 排序）；
2. **可选** —— 调用 PLIP 计算蛋白-配体相互作用（需要安装 ``plip``）；
3. **解析** —— 用纯标准库 ``xml.etree`` 解析 PLIP 的 ``report.xml``，
   抽取 site / pocket / ligand / residue 级的相互作用记录
   （:mod:`plip_parser`）；
4. **统计** —— 按类型、按残基汇总相互作用分布；
5. **可视化** —— 生成 Mermaid 残基相互作用网络图与 PyMOL 标注脚本
   （:mod:`interaction_network`）。

用法示例::

    # 解析已有的 PLIP 报告
    python analyze_interactions.py --input report.xml --output results/

    # 从 AlphaFold zip 选最佳模型并运行 PLIP
    python analyze_interactions.py --input AF-P69905-F1.zip --output results/ --run-plip

    # 只做统计与出图
    python analyze_interactions.py --input report.xml --output results/ \
        --format both --make-mermaid --make-pymol
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "packages"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from alphafold_models import select_best_model, summarize_confidence  # noqa: E402
from biopipe_core import SkillError, cli_parser, get_logger, write_csv  # noqa: E402
from biopipe_core.cli import run_cli  # noqa: E402
from interaction_network import (  # noqa: E402
    build_mermaid,
    build_pymol_script,
    residue_statistics,
)
from plip_parser import INTERACTION_LABELS, parse_plip_xml  # noqa: E402

LOGGER = get_logger("biopipe.protein-interactions")
SKILL_NAME = "protein-interactions"


def _build_parser() -> argparse.ArgumentParser:
    """构造命令行解析器。"""
    parser = cli_parser(
        prog=SKILL_NAME,
        description="蛋白质相互作用分析：AlphaFold 模型选择、PLIP 报告解析与网络图生成",
        skill=SKILL_NAME,
    )
    parser.add_argument(
        "--input",
        required=True,
        metavar="PATH",
        help="PLIP report.xml、AlphaFold zip 包，或包含 report.xml 的目录",
    )
    parser.add_argument("--output", metavar="DIR", help="输出目录（默认取配置中的 results/）")
    parser.add_argument(
        "--run-plip",
        action="store_true",
        help="对输入的结构文件运行 PLIP（需要已安装 plip 可执行程序）",
    )
    parser.add_argument("--plip-structure", metavar="PATH", help="要运行 PLIP 的结构文件（PDB）")
    parser.add_argument("--plip-executable", metavar="PATH", help="PLIP 可执行文件路径")
    parser.add_argument("--plip-timeout", type=int, default=600, help="PLIP 超时秒数（默认 600）")
    parser.add_argument(
        "--make-mermaid", action="store_true", help="生成 Mermaid 相互作用网络图"
    )
    parser.add_argument("--make-pymol", action="store_true", help="生成 PyMOL 标注脚本")
    parser.add_argument(
        "--format",
        choices=("json", "csv", "both"),
        default="both",
        help="输出格式（默认 both）",
    )
    parser.add_argument("--top-residues", type=int, default=10, help="展示的热门残基数量")
    return parser


def _resolve_output_dir(output: str | None) -> Path:
    """解析输出目录。"""
    from biopipe_core import get_config

    config = get_config()
    if output:
        path = Path(output).expanduser()
        if not path.is_absolute():
            path = Path.cwd() / path
    else:
        path = config.output_dir
    path.mkdir(parents=True, exist_ok=True)
    return path


def _find_plip_executable(explicit: str | None = None) -> str | None:
    """查找 PLIP 可执行文件。"""
    if explicit:
        return explicit if Path(explicit).exists() else None
    return shutil.which("plip")


def _run_plip(structure: Path, output_dir: Path, timeout: int = 600) -> Path:
    """调用 PLIP 分析结构，返回 report.xml 路径。

    Raises:
        SkillError: PLIP 不可用或执行失败。
    """
    executable = _find_plip_executable()
    if not executable:
        raise SkillError(
            "未找到 plip 可执行文件。请先 `pip install plip` 或 `conda install -c conda-forge plip`，"
            "或直接提供已生成的 report.xml。"
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        executable,
        "-f",
        str(structure),
        "--inter",
        "A",
        "-v",
        "-t",
        "-x",
        "-y",
        "-o",
        str(output_dir),
    ]
    LOGGER.info("运行 PLIP: %s", " ".join(cmd))
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired as exc:
        raise SkillError(f"PLIP 执行超时（{timeout}s）") from exc
    except OSError as exc:
        raise SkillError(f"PLIP 执行失败: {exc}") from exc

    if proc.returncode != 0:
        raise SkillError(f"PLIP 返回非零退出码 {proc.returncode}: {proc.stderr[:500]}")

    report = output_dir / "report.xml"
    if not report.exists():
        candidates = sorted(output_dir.rglob("report.xml"))
        if not candidates:
            raise SkillError("PLIP 未生成 report.xml")
        report = candidates[0]
    return report


def _collect_reports(input_path: Path) -> list[Path]:
    """收集输入中的 PLIP 报告文件。"""
    if input_path.is_dir():
        reports = sorted(input_path.rglob("report.xml"))
        if not reports:
            raise SkillError(f"在 {input_path} 下没有找到 report.xml")
        return reports
    if input_path.is_file():
        if input_path.suffix.lower() == ".xml":
            return [input_path]
        raise SkillError(f"输入不是 XML 文件: {input_path}")
    raise SkillError(f"输入路径不存在: {input_path}")


def _records_to_rows(records: list[dict]) -> list[dict]:
    """把相互作用记录转成扁平表格行。"""
    rows = []
    for record in records:
        receptor = record.get("receptor") or {}
        ligand = record.get("ligand") or {}
        rows.append(
            {
                "interaction_type": record.get("type", ""),
                "interaction_label": record.get("type_label", ""),
                "receptor_chain": receptor.get("reschain") or receptor.get("chain") or "",
                "receptor_resname": receptor.get("restype") or receptor.get("resname") or "",
                "receptor_resnr": receptor.get("resnr") or "",
                "receptor_atom": receptor.get("atom") or "",
                "ligand_chain": ligand.get("reschain") or ligand.get("chain") or "",
                "ligand_resname": ligand.get("restype") or ligand.get("resname") or "",
                "ligand_resnr": ligand.get("resnr") or "",
                "ligand_atom": ligand.get("atom") or "",
                "distance": record.get("distance", ""),
                "angle": record.get("angle", ""),
            }
        )
    return rows


def _build_summary(report: dict, records: list[dict]) -> dict:
    """生成汇总统计。"""
    stats = residue_statistics(records)
    type_counter = Counter(record.get("type", "") for record in records)
    return {
        "total": len(records),
        "by_type": {label: type_counter.get(tag, 0) for tag, label in INTERACTION_LABELS.items()},
        "by_type_raw": dict(type_counter),
        "top_receptor_residues": stats["receptor"][:10],
        "top_ligand_residues": stats["ligand"][:10],
        "xml_format": report.get("format"),
        "source": report.get("source"),
    }


def analyze_report(
    xml_path: Path,
    output_dir: Path,
    make_mermaid: bool = False,
    make_pymol: bool = False,
    structure: str | None = None,
) -> dict:
    """解析单个 PLIP 报告并写出结果。

    Args:
        xml_path: report.xml 路径。
        output_dir: 输出目录。
        make_mermaid: 是否生成 Mermaid 图。
        make_pymol: 是否生成 PyMOL 脚本。
        structure: PyMOL 脚本中要载入的结构路径。

    Returns:
        含 ``summary`` 与 ``interactions`` 的结果字典。
    """
    report = parse_plip_xml(xml_path)
    records = report["interactions"]
    summary = _build_summary(report, records)

    stem = xml_path.parent.name or xml_path.stem
    written: list[str] = []

    if make_mermaid:
        md_path = output_dir / f"{stem}_interaction_network.md"
        md_path.write_text(build_mermaid(records, title=stem), encoding="utf-8")
        written.append(str(md_path))
        LOGGER.info("已写出 Mermaid 网络图 %s", md_path)

    if make_pymol:
        pml_path = output_dir / f"{stem}_annotation.pml"
        pml_path.write_text(
            build_pymol_script(records, structure=structure or str(xml_path)),
            encoding="utf-8",
        )
        written.append(str(pml_path))
        LOGGER.info("已写出 PyMOL 脚本 %s", pml_path)

    return {
        "structure_name": stem,
        "summary": summary,
        "interactions": records,
        "files": written,
    }


def _write_outputs(results: list[dict], output_dir: Path, fmt: str) -> list[str]:
    """写出汇总 JSON 与残基级 CSV。"""
    written: list[str] = []

    if fmt in ("json", "both"):
        json_path = output_dir / "interaction_analysis.json"
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
        all_rows: list[dict] = []
        for item in results:
            for row in _records_to_rows(item["interactions"]):
                row = {"structure": item["structure_name"], **row}
                all_rows.append(row)
        if all_rows:
            csv_path = output_dir / "interactions_residue_level.csv"
            write_csv(all_rows, csv_path)
            written.append(str(csv_path))
            LOGGER.info("已写出 %s", csv_path)

    return written


def main(args: argparse.Namespace) -> int:
    """入口主函数。"""
    input_path = Path(args.input).expanduser()
    if not input_path.is_absolute():
        input_path = Path.cwd() / input_path
    output_dir = _resolve_output_dir(args.output)

    results: list[dict] = []

    # --- 情形 1：输入是 AlphaFold zip，需要先选模型再跑 PLIP ---
    if input_path.is_file() and input_path.suffix.lower() == ".zip":
        LOGGER.info("检测到 AlphaFold zip 包，按置信度选择最佳模型")
        best = select_best_model(input_path)
        if best is None:
            raise SkillError("未在 zip 中找到可用的模型文件")
        LOGGER.info(
            "最佳模型: %s (model_%s, 分数=%.3f)",
            best["model_member"],
            best["model_index"],
            best["score"],
        )
        conf_summary = summarize_confidence(best["confidence"])

        structure_path: Path | None = None
        if args.run_plip:
            if not args.plip_structure:
                raise SkillError("使用 zip 输入时必须通过 --plip-structure 指定结构文件")
            structure_path = Path(args.plip_structure)
            plip_dir = output_dir / best.get("folder", "model")
            report_path = _run_plip(structure_path, plip_dir, timeout=args.plip_timeout)
            item = analyze_report(
                report_path,
                output_dir,
                make_mermaid=args.make_mermaid,
                make_pymol=args.make_pymol,
                structure=str(structure_path),
            )
            item["alphafold"] = {
                "model_member": best["model_member"],
                "model_index": best["model_index"],
                "score": best["score"],
                "confidence": conf_summary,
            }
            results.append(item)
        else:
            # 未跑 PLIP：只输出模型选择信息
            LOGGER.warning("未指定 --run-plip，仅输出最佳模型信息")
            results.append(
                {
                    "structure_name": best.get("folder", "model"),
                    "summary": {"total": 0, "by_type": {}, "by_type_raw": {}},
                    "interactions": [],
                    "files": [],
                    "alphafold": {
                        "model_member": best["model_member"],
                        "model_index": best["model_index"],
                        "score": best["score"],
                        "confidence": conf_summary,
                    },
                }
            )

    # --- 情形 2：输入是 XML 报告或包含报告的目录 ---
    else:
        if args.run_plip:
            if not args.plip_structure:
                raise SkillError("--run-plip 需要同时提供 --plip-structure")
            structure_path = Path(args.plip_structure)
            report_path = _run_plip(structure_path, output_dir, timeout=args.plip_timeout)
            reports = [report_path]
        else:
            reports = _collect_reports(input_path)

        for report_path in reports:
            item = analyze_report(
                report_path,
                output_dir,
                make_mermaid=args.make_mermaid,
                make_pymol=args.make_pymol,
                structure=str(args.plip_structure) if args.plip_structure else None,
            )
            results.append(item)

    if not results:
        raise SkillError("没有可分析的数据")

    written = _write_outputs(results, output_dir, args.format)

    total = sum(item["summary"]["total"] for item in results)
    LOGGER.info("完成：%d 个报告，共 %d 条相互作用", len(results), total)
    for item in results:
        by_type = item["summary"]["by_type"]
        active = {k: v for k, v in by_type.items() if v}
        LOGGER.info(
            "  %s: %s",
            item["structure_name"],
            ", ".join(f"{k}={v}" for k, v in active.items()) or "无相互作用",
        )
    for path in written:
        LOGGER.info("输出: %s", path)
    return 0


if __name__ == "__main__":
    sys.exit(run_cli(main, parser=_build_parser(), skill=SKILL_NAME, version="1.0.0"))
