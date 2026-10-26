"""生成 AlphaFold 最佳模型选择的测试夹具。

夹具模拟 AlphaFold Server / ColabFold 的真实输出布局：

- 多个 ``*_model_<rank>.cif`` 结构文件
- 多个 ``*_summary_confidences_<rank>.json`` 置信度文件
- 排名**故意不与文件编号一致**——rank 2 的文件置信度最高，
  这样"按字典序取第一个"（老脚本的做法）会挑错。

另外还生成一份**命名不规则**的目录（文件名里没有 model/cif 编号），
用来验证新脚本不依赖具体命名约定。

运行：``python tests/make_alphafold_select_fixtures.py``
"""

from __future__ import annotations

import json
from pathlib import Path

FIXTURES = Path(__file__).resolve().parent / "fixtures_alphafold_select"

#: 各排名对应的置信度指标。
#:
#: 关键设计：``rank`` 编号与 ``ranking_score`` **不成正比**——
#: 编号 1 的文件置信度反而最低。这正是"按文件名排序取第一个"会挑错的场景。
CONFIDENCES = {
    0: {"ptm": 0.72, "iptm": 0.61, "ranking_score": 0.551, "mean_plddt": 82.1,
        "fraction_disordered": 0.11, "has_clash": 0.0},
    1: {"ptm": 0.55, "iptm": 0.48, "ranking_score": 0.472, "mean_plddt": 74.3,
        "fraction_disordered": 0.19, "has_clash": 0.0},
    2: {"ptm": 0.93, "iptm": 0.88, "ranking_score": 0.884, "mean_plddt": 91.7,
        "fraction_disordered": 0.03, "has_clash": 0.0},
    3: {"ptm": 0.88, "iptm": 0.81, "ranking_score": 0.816, "mean_plddt": 89.2,
        "fraction_disordered": 0.05, "has_clash": 0.0},
}

#: 一个残基的 mmCIF 原子记录（3 个原子）。
ATOM_SITE_LOOP = """data_test
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1 N N 1 MET A 1 1 ? 11.104 12.220 10.383 1.00 88.12 ? 1 MET A N 1
ATOM 2 C CA 2 MET A 1 1 ? 12.560 12.660 10.900 1.00 88.12 ? 1 MET A CA 1
ATOM 3 C C 3 MET A 1 1 ? 13.010 13.400 9.910 1.00 88.12 ? 1 MET A C 1
#
"""


def build_directory(directory: Path, prefix: str, ranks: list[int]) -> None:
    """生成一个标准的预测结果目录。

    Args:
        directory: 目标目录。
        prefix: 文件名前缀（模拟 ``fold_<name>`` 之类）。
        ranks: 要生成的排名编号列表。
    """
    directory.mkdir(parents=True, exist_ok=True)
    for rank in ranks:
        (directory / f"{prefix}_model_{rank}.cif").write_text(
            ATOM_SITE_LOOP, encoding="utf-8"
        )
        confidence = dict(CONFIDENCES[rank])
        confidence["fraction_disordered"] = confidence.get("fraction_disordered", 0.0)
        (directory / f"{prefix}_summary_confidences_{rank}.json").write_text(
            json.dumps(confidence, indent=2), encoding="utf-8"
        )


def build_irregular_directory(directory: Path) -> None:
    """生成一个**命名不规则**的目录。

    文件名里既没有 ``model_<n>`` 也没有 ``cif_<n>``，只有置信度文件里
    带编号。用来验证新脚本不靠命名约定也能把结构与置信度配对。

    同时这里的置信度**故意缺** ``mean_plddt``，用来验证"缺失指标返回
    ``None`` 而不是 0"。
    """
    directory.mkdir(parents=True, exist_ok=True)
    for rank in (0, 2):
        (directory / f"structure_{rank}.pdb").write_text(
            "ATOM      1  N   MET A   1      11.104  12.220  10.383  1.00 88.12           N\n"
            "ATOM      2  CA  MET A   1      12.560  12.660  10.900  1.00 88.12           C\n"
            "END\n",
            encoding="utf-8",
        )
        confidence = {
            key: value
            for key, value in CONFIDENCES[rank].items()
            if key != "mean_plddt"
        }
        (directory / f"summary_confidences_{rank}.json").write_text(
            json.dumps(confidence, indent=2), encoding="utf-8"
        )


def build_empty_confidence(directory: Path) -> None:
    """生成一个置信度 JSON 损坏的目录。

    用来验证解析失败时不会让整个脚本崩溃。
    """
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "fold_model_0.cif").write_text(ATOM_SITE_LOOP, encoding="utf-8")
    (directory / "fold_summary_confidences_0.json").write_text(
        "{ broken json", encoding="utf-8"
    )


def main() -> None:
    """生成全部夹具。"""
    FIXTURES.mkdir(parents=True, exist_ok=True)

    # 标准布局：排名编号与置信度不成正比
    build_directory(FIXTURES / "standard", "fold", [0, 1, 2, 3])

    # 单模型目录
    build_directory(FIXTURES / "single", "fold", [0])

    # 命名不规则
    build_irregular_directory(FIXTURES / "irregular")

    # 置信度 JSON 损坏
    build_empty_confidence(FIXTURES / "broken")

    # 空目录
    (FIXTURES / "empty").mkdir(parents=True, exist_ok=True)

    print(f"夹具已生成到: {FIXTURES}")


if __name__ == "__main__":
    main()
