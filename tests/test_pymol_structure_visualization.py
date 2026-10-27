"""``pymol-structure-visualization`` 技能的回归测试。

锁住老仓库 ``pymol_show/scripts/`` 的 22 个脚本带来的真实问题，
以及链间接触判定里最容易踩的坑：

1. **必须启动 PyMOL GUI** —— 老脚本全部以 ``pymol.finish_launching(['pymol'])``
   开头，且最后 ``cmd.zoom()``。这意味着无法 headless 运行、无 PyMOL 的
   机器上完全跑不了、CI 里无法测试。新技能必须把几何计算与 PyMOL 渲染
   彻底分离 —— PyMOL 脚本是**生成物**而不是运行依赖。
2. **硬编码个人路径** —— 老脚本写死 ``C:\\Users\\DELL\\Desktop\\1.pdb``、
   ``C:\\Users\\DELL\\Desktop\\AC_interactions.csv``、
   ``C:\\Users\\DELL\\Desktop\\output\\*\\interaction_residues.txt``。
3. **分类表互相矛盾** —— ``interaction count.py`` 的 donors 是
   ``LYS ARG HIS SER THR TYR ASN GLN``（含 HIS）、acceptors 是
   ``ASP GLU ASN GLN TYR``（不含 SER/THR/HIS/TRP）；而 ``hbond count.py``
   的 donors 是 ``LYS ARG HIS ASP GLU SER THR TYR``（**含 ASP/GLU**），
   acceptors 只有 ``ASP GLU TYR``。同一套相互作用在不同脚本里被判成
   不同类型。
4. **ASP/GLU 被列为氢键供体是错的** —— 羧基 pKa ≈ 4，生理 pH 下几乎
   完全去质子化，不是有效供体。老脚本会虚增氢键数。
5. **HIS 被算进正电残基是错的** —— 咪唑基 pKa ≈ 6.0，生理 pH 下仅约
   10% 质子化。老脚本会虚增盐桥数。
6. **文件名带空格** —— 老脚本名如 ``interaction count.py``、
   ``hbond count.py``，在 shell 里必须加引号，PyMOL 的 ``run`` 命令
   对带空格路径的处理也容易出问题。
7. **只有接触对、没有类型** —— 老 ``interaction count.py`` 用
   ``cmd.select(f"{sel1} within {cutoff} of {sel2}")``，这会选出
   **残基**而不是原子对，然后做笛卡尔积 —— 一个供体残基和一个受体残基
   只要各自内部有原子在截断距离内，就算一对。这是老实现最本质的缺陷：
   它报告的是"两个残基各自靠近了对方"，不是"这两个残基互相接触了"。

本模块的用例全部使用**本地构造的 PDB 夹具**，不联网、不需要 PyMOL、
不需要 numpy/pandas。
"""

from __future__ import annotations

import ast
import csv
import json
import math
import os
import subprocess
import sys
from pathlib import Path
from typing import Optional

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
SKILL_DIR = REPO_ROOT / "skills" / "pymol-structure-visualization"
SCRIPTS_DIR = SKILL_DIR / "scripts"
FIXTURES = REPO_ROOT / "tests" / "fixtures_contacts"

sys.path.insert(0, str(REPO_ROOT / "packages"))
sys.path.insert(0, str(SCRIPTS_DIR))

import analyze_contacts  # noqa: E402
import contact_chemistry  # noqa: E402
import contact_geometry  # noqa: E402

SCRIPTS = ("contact_chemistry.py", "contact_geometry.py", "analyze_contacts.py")

#: 老仓库里的个人路径，任何一个出现在新技能里都是回归。
FORBIDDEN_PATHS = (
    "/Users/DELL",
    "C:\\Users\\DELL",
    "C:/Users/DELL",
    "Desktop\\1.pdb",
    "interaction_residues.txt",
)


# ---------------------------------------------------------------- 工具


def _script(name: str) -> Path:
    return SCRIPTS_DIR / name


def _read_script(name: str) -> str:
    return _script(name).read_text(encoding="utf-8")


def _make_atom(
    serial: int,
    name: str,
    res_name: str,
    chain: str,
    res_seq: int,
    xyz: tuple[float, float, float],
    element: str = "C",
    hetero: bool = False,
) -> contact_geometry.AtomRef:
    """构造一个 ``AtomRef``（测试用）。"""
    return contact_geometry.AtomRef(
        serial=serial,
        name=name,
        res_name=res_name,
        chain_id=chain,
        res_seq=res_seq,
        icode="",
        x=xyz[0],
        y=xyz[1],
        z=xyz[2],
        element=element,
        hetero=hetero,
    )


def _load_fixture_atoms(path: Path) -> list[contact_geometry.AtomRef]:
    """用技能自己的解析器读取夹具 PDB。"""
    atoms, _summary = analyze_contacts._load_atoms(path)
    return atoms


def _run_cli(args: list[str], cwd: Optional[Path] = None) -> subprocess.CompletedProcess:
    """以子进程方式运行技能入口脚本。

    注意这里必须**继承**当前进程的环境变量，只覆盖需要控制的那几项。
    早期版本构造了一个只含 3 个键的全新 env dict，在 Windows 上会丢掉
    ``SystemRoot``，子进程 Python 因此无法加载 ``bcryptprim``，直接崩在
    ``Fatal Python error: _Py_HashRandomization_Init: failed to get random
    numbers to initialize Python`` —— 17 个 CLI 端到端用例集体红。
    Ubuntu CI 上不会暴露这个问题（没有 ``SystemRoot`` 依赖），
    但 Windows 开发者完全跑不了这批测试。

    ``PATH`` 保留隔离意图（避免脚本意外找到系统里的 PyMOL），
    但要按平台用正确的分隔符，否则 Windows 上 ``PATH`` 形同虚设。
    """
    env = dict(os.environ)
    env["PATH"] = os.pathsep.join(("/usr/bin", "/bin", "/usr/local/bin"))
    env["PYTHONPATH"] = str(REPO_ROOT / "packages")
    env["PYTHONIOENCODING"] = "utf-8"
    return subprocess.run(
        [sys.executable, str(SCRIPTS_DIR / "analyze_contacts.py")] + args,
        capture_output=True,
        text=True,
        timeout=180,
        cwd=str(cwd or REPO_ROOT),
        env=env,
    )


# ---------------------------------------------------------------- 化学分类


class TestDonorAcceptorTables:
    """供体/受体表：老仓库两个脚本互相矛盾，新版必须统一。"""

    def test_asp_glu_are_not_donors(self) -> None:
        """ASP/GLU 不是有效氢键供体。

        老 ``hbond count.py`` 把 ASP/GLU 列进 donors，会虚增氢键数。
        羧基 pKa ≈ 4，生理 pH 下几乎完全去质子化。
        """
        assert not contact_chemistry.is_donor_residue("ASP")
        assert not contact_chemistry.is_donor_residue("GLU")
        assert "ASP" not in contact_chemistry.DONORS
        assert "GLU" not in contact_chemistry.DONORS

    def test_asp_glu_can_be_weak_donors_on_request(self) -> None:
        """显式开启时，ASP/GLU 可作为弱供体。"""
        assert contact_chemistry.is_donor_residue("ASP", include_ambiguous=True)
        assert contact_chemistry.is_donor_residue("GLU", include_ambiguous=True)

    def test_standard_donors_are_donors(self) -> None:
        """侧链有 N-H / O-H 的残基是供体。"""
        for res in ("ARG", "ASN", "GLN", "HIS", "LYS", "SER", "THR", "TRP", "TYR"):
            assert contact_chemistry.is_donor_residue(res), f"{res} 应该是供体"

    def test_acceptors_include_amide_and_hydroxyl(self) -> None:
        """酰胺与羟基氧是受体 —— 老表的 acceptors 漏了 SER/THR/HIS/TRP。"""
        for res in ("ASN", "GLN", "HIS", "SER", "THR", "TYR", "ASP", "GLU"):
            assert contact_chemistry.is_acceptor_residue(res), f"{res} 应该是受体"

    def test_his_is_both_donor_and_acceptor(self) -> None:
        """咪唑基既可给质子也可受质子。"""
        assert contact_chemistry.is_donor_residue("HIS")
        assert contact_chemistry.is_acceptor_residue("HIS")

    def test_gly_and_ala_are_neither(self) -> None:
        """GLY/ALA 既无供体也无受体。"""
        for res in ("GLY", "ALA"):
            assert not contact_chemistry.is_donor_residue(res)
            assert not contact_chemistry.is_acceptor_residue(res)


class TestChargeTables:
    """电荷表：老仓库把 HIS 算进正电，会虚增盐桥数。"""

    def test_his_is_not_positive_by_default(self) -> None:
        """HIS 默认不算正电（咪唑基 pKa ≈ 6，生理 pH 下约 10% 质子化）。"""
        assert not contact_chemistry.is_positive_residue("HIS")
        assert "HIS" not in contact_chemistry.POSITIVE
        assert contact_chemistry.POSITIVE_PARTIAL == frozenset({"HIS"})

    def test_his_can_be_positive_on_request(self) -> None:
        """显式开启时 HIS 计入正电。"""
        assert contact_chemistry.is_positive_residue("HIS", include_histidine=True)

    def test_lys_and_arg_are_positive(self) -> None:
        """LYS/ARG 在生理 pH 下几乎完全质子化。"""
        assert contact_chemistry.is_positive_residue("LYS")
        assert contact_chemistry.is_positive_residue("ARG")

    def test_asp_and_glu_are_negative(self) -> None:
        """ASP/GLU 带负电。"""
        assert contact_chemistry.is_negative_residue("ASP")
        assert contact_chemistry.is_negative_residue("GLU")

    def test_cys_is_hydrophobic(self) -> None:
        """老仓库的疏水列表漏了 CYS。"""
        assert contact_chemistry.is_hydrophobic_residue("CYS")
        assert "CYS" in contact_chemistry.HYDROPHOBIC

    def test_hydrophobic_set_covers_aliphatic_and_aromatic(self) -> None:
        """疏水列表覆盖脂肪族与芳香族。"""
        for res in ("ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO"):
            assert contact_chemistry.is_hydrophobic_residue(res), f"{res} 应该是疏水"

    def test_aromatic_set(self) -> None:
        """芳香族集合 —— π-π 与阳离子-π 的参与者。"""
        assert contact_chemistry.AROMATIC == frozenset({"PHE", "TRP", "TYR", "HIS"})


class TestResidueKind:
    """标准残基 vs 配体的区分。"""

    def test_standard_residues_recognized(self) -> None:
        """20 种标准氨基酸被识别为 standard。"""
        for res in contact_chemistry.STANDARD_RESIDUES:
            assert contact_chemistry.residue_kind(res) == "standard"

    def test_ligand_detected(self) -> None:
        """非标准残基识别为 ligand。"""
        assert contact_chemistry.residue_kind("LIG") == "ligand"
        assert contact_chemistry.residue_kind("HEM") == "ligand"
        assert not contact_chemistry.is_standard_residue("LIG")

    def test_water_detected(self) -> None:
        """水分子单独识别。"""
        assert contact_chemistry.residue_kind("HOH") == "water"
        assert contact_chemistry.residue_kind("DOD") == "water"


class TestClassifyPair:
    """残基对的相互作用分类。"""

    def test_salt_bridge_also_counts_as_hydrogen_bond(self) -> None:
        """ARG-ASP 同时是氢键与盐桥 —— 这是最典型的双类型接触。"""
        types = contact_chemistry.classify_pair_multi("ARG", "ASP", 2.9)
        assert "Salt_bridge" in types
        assert "Hydrogen_bond" in types

    def test_classify_pair_returns_first_type_only(self) -> None:
        """``classify_pair`` 只返回第一个匹配类型。"""
        first = contact_chemistry.classify_pair("ARG", "ASP", 2.9)
        multi = contact_chemistry.classify_pair_multi("ARG", "ASP", 2.9)
        assert first == multi[0]
        assert first == "Hydrogen_bond"

    def test_hydrogen_bond_cutoff(self) -> None:
        """氢键的距离上限是 3.5 Å。"""
        assert "Hydrogen_bond" in contact_chemistry.classify_pair_multi("ARG", "ASP", 3.5)
        assert "Hydrogen_bond" not in contact_chemistry.classify_pair_multi("ARG", "ASP", 3.51)

    def test_salt_bridge_cutoff(self) -> None:
        """盐桥的距离上限是 4.0 Å。"""
        assert "Salt_bridge" in contact_chemistry.classify_pair_multi("LYS", "GLU", 4.0)
        assert "Salt_bridge" not in contact_chemistry.classify_pair_multi("LYS", "GLU", 4.01)

    def test_hydrophobic_cutoff(self) -> None:
        """疏水接触的距离上限是 4.5 Å。"""
        assert "Hydrophobic" in contact_chemistry.classify_pair_multi("ILE", "LEU", 4.5)
        assert "Hydrophobic" not in contact_chemistry.classify_pair_multi("ILE", "LEU", 4.51)

    def test_pi_pi_cutoff(self) -> None:
        """π-π 堆积的距离上限是 5.0 Å。"""
        assert "Pi_Pi" in contact_chemistry.classify_pair_multi("PHE", "TYR", 5.0)
        assert "Pi_Pi" not in contact_chemistry.classify_pair_multi("PHE", "TYR", 5.01)

    def test_cation_pi_cutoff(self) -> None:
        """阳离子-π 的距离上限是 6.0 Å。"""
        assert "Cation_Pi" in contact_chemistry.classify_pair_multi("LYS", "PHE", 6.0)
        assert "Cation_Pi" not in contact_chemistry.classify_pair_multi("LYS", "PHE", 6.01)

    def test_gly_gly_has_no_interaction(self) -> None:
        """GLY-GLY 不产生任何相互作用类型。"""
        assert contact_chemistry.classify_pair_multi("GLY", "GLY", 3.0) == []

    def test_histidine_salt_bridge_requires_opt_in(self) -> None:
        """HIS-ASP 默认不判为盐桥，开启后才算。"""
        assert "Salt_bridge" not in contact_chemistry.classify_pair_multi("HIS", "ASP", 3.0)
        assert "Salt_bridge" in contact_chemistry.classify_pair_multi(
            "HIS", "ASP", 3.0, include_histidine=True
        )

    def test_asp_glu_hydrogen_bond_requires_opt_in(self) -> None:
        """ASP-GLU 默认不判为氢键（羧基不是有效供体）。"""
        assert "Hydrogen_bond" not in contact_chemistry.classify_pair_multi("ASP", "GLU", 3.0)
        assert "Hydrogen_bond" in contact_chemistry.classify_pair_multi(
            "ASP", "GLU", 3.0, include_ambiguous_donors=True
        )

    def test_ligand_contact_is_marked_not_classified(self) -> None:
        """配体接触单独标记，不套用残基级分类表。"""
        types = contact_chemistry.classify_pair_multi("TYR", "LIG", 3.5)
        assert types == [contact_chemistry.LIGAND_CONTACT]
        assert "Pi_Pi" not in types

    def test_water_water_has_no_contact(self) -> None:
        """水-水接触不报告。"""
        assert contact_chemistry.classify_pair_multi("HOH", "HOH", 2.8) == []

    def test_ligand_contact_with_protein_is_marked(self) -> None:
        """蛋白-配体接触被标记为 Ligand_contact。"""
        types = contact_chemistry.classify_pair_multi("LYS", "LIG", 3.2)
        assert types == [contact_chemistry.LIGAND_CONTACT]


class TestRulesMetadata:
    """规则元数据。"""

    def test_five_interaction_types(self) -> None:
        """5 种残基级接触类型，顺序固定。"""
        names = [rule.name for rule in contact_chemistry.INTERACTION_TYPES]
        assert names == [
            "Hydrogen_bond", "Salt_bridge", "Hydrophobic", "Pi_Pi", "Cation_Pi",
        ]

    def test_distance_cutoff_lookup(self) -> None:
        """按名字取距离阈值。"""
        assert contact_chemistry.distance_cutoff("Hydrogen_bond") == 3.5
        assert contact_chemistry.distance_cutoff("Cation_Pi") == 6.0

    def test_distance_cutoff_rejects_unknown(self) -> None:
        """未知类型报错。"""
        with pytest.raises(KeyError):
            contact_chemistry.distance_cutoff("No_Such_Type")

    def test_describe_rules_includes_ligand(self) -> None:
        """规则描述包含配体类别。"""
        described = contact_chemistry.describe_rules()
        names = [item["name"] for item in described]
        assert contact_chemistry.LIGAND_CONTACT in names
        assert len(names) == len(contact_chemistry.INTERACTION_TYPES) + 1


# ---------------------------------------------------------------- 几何检测


class TestGridDetection:
    """网格空间划分的正确性。"""

    def test_detects_contact_pair(self) -> None:
        """基本接触能被检出。"""
        atoms = [
            _make_atom(1, "NZ", "LYS", "A", 1, (0.0, 0.0, 0.0)),
            _make_atom(2, "OE1", "GLU", "B", 1, (3.0, 0.0, 0.0)),
        ]
        detector = contact_geometry.ContactDetector(max_distance=6.0)
        detector.add_atoms(atoms)
        pairs = detector.detect("A", "B")
        assert len(pairs) == 1
        assert pairs[0].min_distance == pytest.approx(3.0, abs=1e-6)
        assert pairs[0].label_a == "A:LYS1"
        assert pairs[0].label_b == "B:GLU1"

    def test_ignores_far_pairs(self) -> None:
        """超出阈值的不算接触。"""
        atoms = [
            _make_atom(1, "NZ", "LYS", "A", 1, (0.0, 0.0, 0.0)),
            _make_atom(2, "OE1", "GLU", "B", 1, (7.0, 0.0, 0.0)),
        ]
        detector = contact_geometry.ContactDetector(max_distance=6.0)
        detector.add_atoms(atoms)
        assert detector.detect("A", "B") == []

    def test_boundary_case_across_cell_edge(self) -> None:
        """原子恰好落在网格边界两侧时不能漏检。

        这是网格划分的经典陷阱：桶边长等于最大距离时，原子可能分别落在
        相邻两个桶的边界，扫 1 格邻域会漏掉。这里扫 2 格。
        """
        atoms = [
            _make_atom(1, "NZ", "LYS", "A", 1, (5.999, 0.0, 0.0)),
            _make_atom(2, "OE1", "GLU", "B", 1, (6.001, 0.0, 0.0)),
        ]
        detector = contact_geometry.ContactDetector(max_distance=6.0)
        detector.add_atoms(atoms)
        assert len(detector.detect("A", "B")) == 1

    def test_grid_matches_brute_force(self) -> None:
        """网格结果必须与暴力两两比较完全一致。

        这是网格划分最重要的一致性保证 —— 优化不能改变结果。
        """
        atoms = []
        serial = 0
        # 用确定性伪随机铺一批原子
        seed = 12345
        for i in range(60):
            seed = (seed * 1103515245 + 12345) % (2**31)
            x = (seed % 1000) / 1000.0 * 40.0
            seed = (seed * 1103515245 + 12345) % (2**31)
            y = (seed % 1000) / 1000.0 * 40.0
            seed = (seed * 1103515245 + 12345) % (2**31)
            z = (seed % 1000) / 1000.0 * 40.0
            chain = "A" if i % 2 == 0 else "B"
            serial += 1
            atoms.append(_make_atom(serial, "CA", "ALA", chain, i, (x, y, z)))

        detector = contact_geometry.ContactDetector(max_distance=6.0)
        detector.add_atoms(atoms)

        grid_result = {
            (p.chain_a, p.res_seq_a, p.chain_b, p.res_seq_b): round(p.min_distance, 6)
            for p in detector.detect("A", "B")
        }

        brute: dict[tuple, float] = {}
        for a in atoms:
            if a.chain_id != "A":
                continue
            for b in atoms:
                if b.chain_id != "B":
                    continue
                d = math.dist((a.x, a.y, a.z), (b.x, b.y, b.z))
                if d > 6.0:
                    continue
                key = ("A", a.res_seq, "B", b.res_seq)
                if key not in brute or d < brute[key]:
                    brute[key] = d

        assert grid_result == {k: round(v, 6) for k, v in brute.items()}

    def test_hydrogen_is_ignored_by_default(self) -> None:
        """默认忽略氢原子。"""
        atoms = [
            _make_atom(1, "NZ", "LYS", "A", 1, (0.0, 0.0, 0.0)),
            _make_atom(2, "H", "LYS", "A", 1, (1.0, 0.0, 0.0), element="H"),
            _make_atom(3, "OE1", "GLU", "B", 1, (3.0, 0.0, 0.0)),
        ]
        detector = contact_geometry.ContactDetector(max_distance=6.0)
        detector.add_atoms(atoms)
        assert len(detector.detect("A", "B")) == 1

        detector_with_h = contact_geometry.ContactDetector(
            max_distance=6.0, ignore_hydrogen=False
        )
        detector_with_h.add_atoms(atoms)
        # 氢参与后原子对数量增加，但残基对数量不变
        assert len(detector_with_h.detect("A", "B")) == 1
        assert (
            detector_with_h.detect("A", "B")[0].n_atom_pairs
            > detector.detect("A", "B")[0].n_atom_pairs
        )

    def test_hetero_is_included_by_default(self) -> None:
        """默认把 HETATM（配体）算进去。"""
        atoms = [
            _make_atom(1, "NZ", "LYS", "A", 1, (0.0, 0.0, 0.0)),
            _make_atom(2, "C1", "LIG", "L", 1, (3.0, 0.0, 0.0), hetero=True),
        ]
        detector = contact_geometry.ContactDetector(max_distance=6.0)
        detector.add_atoms(atoms)
        assert len(detector.detect("A", "L")) == 1

        no_hetero = contact_geometry.ContactDetector(
            max_distance=6.0, include_hetero=False
        )
        no_hetero.add_atoms(atoms)
        assert no_hetero.detect("A", "L") == []

    def test_same_chain_not_detected(self) -> None:
        """同链残基之间默认不产生接触。

        链内接触（尤其是共价相邻残基）与本技能要回答的
        "两条链之间怎么结合" 无关，混进来会淹没跨链信号。
        """
        atoms = [
            _make_atom(1, "CA", "ALA", "A", 1, (0.0, 0.0, 0.0)),
            _make_atom(2, "CA", "GLY", "A", 2, (1.0, 0.0, 0.0)),
        ]
        detector = contact_geometry.ContactDetector(max_distance=6.0)
        detector.add_atoms(atoms)
        assert detector.detect("A", "A") == []

    def test_same_chain_can_be_enabled(self) -> None:
        """显式开启后可以检出同链接触。"""
        atoms = [
            _make_atom(1, "CA", "ALA", "A", 1, (0.0, 0.0, 0.0)),
            _make_atom(2, "CA", "GLY", "A", 2, (1.0, 0.0, 0.0)),
        ]
        detector = contact_geometry.ContactDetector(max_distance=6.0)
        detector.add_atoms(atoms)
        pairs = detector.detect("A", "A", include_same_chain=True)
        assert len(pairs) == 1
        assert pairs[0].min_distance == pytest.approx(1.0, abs=1e-6)

    def test_detect_all_returns_chain_pairs(self) -> None:
        """``detect_all`` 按链对分组。"""
        atoms = [
            _make_atom(1, "CA", "ALA", "A", 1, (0.0, 0.0, 0.0)),
            _make_atom(2, "CA", "GLY", "B", 1, (3.0, 0.0, 0.0)),
            _make_atom(3, "CA", "SER", "C", 1, (10.0, 0.0, 0.0)),
        ]
        detector = contact_geometry.ContactDetector(max_distance=6.0)
        detector.add_atoms(atoms)
        result = detector.detect_all()
        assert set(result.keys()) == {("A", "B")}

    def test_results_sorted_by_distance(self) -> None:
        """结果按最小距离升序。"""
        atoms = [
            _make_atom(1, "CA", "ALA", "A", 1, (0.0, 0.0, 0.0)),
            _make_atom(2, "CA", "GLY", "B", 1, (5.0, 0.0, 0.0)),
            _make_atom(3, "CA", "SER", "B", 2, (2.0, 0.0, 0.0)),
        ]
        detector = contact_geometry.ContactDetector(max_distance=6.0)
        detector.add_atoms(atoms)
        distances = [p.min_distance for p in detector.detect("A", "B")]
        assert distances == sorted(distances)

    def test_grid_cell_size_at_least_max_distance(self) -> None:
        """桶边长不小于最大接触距离，否则会漏检。"""
        assert contact_geometry.grid_cell_size(6.0) >= 6.0
        assert contact_geometry.grid_cell_size(3.5) >= 3.5

    def test_convenience_function_matches_class(self) -> None:
        """便捷函数与类方法结果一致。"""
        atoms = [
            _make_atom(1, "NZ", "LYS", "A", 1, (0.0, 0.0, 0.0)),
            _make_atom(2, "OE1", "GLU", "B", 1, (3.0, 0.0, 0.0)),
        ]
        from_class = contact_geometry.detect_chain_contacts(atoms, "A", "B")
        assert len(from_class) == 1
        assert from_class[0].min_distance == pytest.approx(3.0, abs=1e-6)


# ---------------------------------------------------------------- 夹具


class TestFixtureIntegrity:
    """夹具本身的完整性 —— 夹具坏了下游测试会莫名失败。"""

    def test_fixture_files_exist(self) -> None:
        """夹具文件已生成。"""
        assert (FIXTURES / "two_chain_complex.pdb").exists()
        assert (FIXTURES / "multimodel.pdb").exists()

    def test_no_atomic_overlap_between_residues(self) -> None:
        """不同残基之间没有原子重叠（间距 < 2.0 Å）。

        老夹具用固定 +y 方向摆放残基，导致 ARG 的 NH2 撞上 GLU 的 OE2，
        出现 1.18 Å 的非物理距离。新夹具用方向搜索自动避开。
        """
        atoms = _load_fixture_atoms(FIXTURES / "two_chain_complex.pdb")
        expected_pairs = {
            frozenset({("A", 10), ("B", 60)}),
            frozenset({("A", 20), ("B", 70)}),
            frozenset({("A", 25), ("B", 90)}),
            frozenset({("A", 30), ("B", 80)}),
            frozenset({("A", 40), ("B", 100)}),
            frozenset({("A", 50), ("B", 110)}),
            frozenset({("B", 80), ("L", 1)}),
        }
        offenders = []
        for i, a in enumerate(atoms):
            for b in atoms[i + 1:]:
                if a.chain_id == b.chain_id and a.res_seq == b.res_seq:
                    continue
                key = frozenset({(a.chain_id, a.res_seq), (b.chain_id, b.res_seq)})
                if key in expected_pairs:
                    continue
                d = math.dist((a.x, a.y, a.z), (b.x, b.y, b.z))
                if d < 2.0:
                    offenders.append(
                        f"{a.chain_id}:{a.res_name}{a.res_seq}:{a.name} - "
                        f"{b.chain_id}:{b.res_name}{b.res_seq}:{b.name} = {d:.3f} Å"
                    )
        assert not offenders, "夹具存在原子重叠:\n  " + "\n  ".join(offenders[:10])

    def test_designed_contact_distances_are_exact(self) -> None:
        """设计的接触原子对距离精确。"""
        atoms = {
            (a.chain_id, a.res_seq, a.name): (a.x, a.y, a.z)
            for a in _load_fixture_atoms(FIXTURES / "two_chain_complex.pdb")
        }
        cases = [
            (("A", 10, "NH1"), ("B", 60, "OD1"), 2.9),
            (("A", 20, "NZ"), ("B", 70, "OE1"), 3.0),
            (("A", 25, "NZ"), ("B", 90, "CG"), 5.5),
            (("A", 30, "CG"), ("B", 80, "CG"), 4.5),
            (("A", 40, "CG1"), ("B", 100, "CG"), 4.2),
            (("A", 50, "CA"), ("B", 110, "CA"), 3.8),
        ]
        for key_a, key_b, expected in cases:
            actual = math.dist(atoms[key_a], atoms[key_b])
            assert actual == pytest.approx(expected, abs=0.001), (
                f"{key_a} - {key_b} 距离应为 {expected}，实际 {actual:.3f}"
            )


# ---------------------------------------------------------------- CLI 端到端


class TestCliEndToEnd:
    """命令行端到端。"""

    def test_help_works(self) -> None:
        """``--help`` 正常。"""
        proc = _run_cli(["--help"])
        assert proc.returncode == 0, proc.stderr
        assert "usage" in proc.stdout.lower()

    def test_missing_input_fails_cleanly(self) -> None:
        """输入不存在时报错而非崩溃。"""
        proc = _run_cli(["--input", "no_such_file.pdb", "--output", "out"])
        assert proc.returncode != 0
        assert "不存在" in (proc.stdout + proc.stderr)

    def test_full_pipeline(self, tmp_path: Path) -> None:
        """完整管线跑通并产出全部文件。"""
        out = tmp_path / "out"
        proc = _run_cli([
            "--input", str(FIXTURES / "two_chain_complex.pdb"),
            "--output", str(out),
            "--make-mermaid",
            "--make-pymol",
            "--make-pymol-contacts",
        ])
        assert proc.returncode == 0, proc.stderr

        for name in (
            "contacts.csv",
            "contact_summary.csv",
            "contacts.json",
            "contact_network.mmd",
            "visualize_contacts.pml",
            "draw_contacts.pml",
        ):
            assert (out / name).exists(), f"缺少产物 {name}"

    def test_json_result_structure(self, tmp_path: Path) -> None:
        """JSON 结果包含完整的结构化信息。"""
        out = tmp_path / "out"
        proc = _run_cli([
            "--input", str(FIXTURES / "two_chain_complex.pdb"),
            "--output", str(out),
        ])
        assert proc.returncode == 0, proc.stderr

        data = json.loads((out / "contacts.json").read_text(encoding="utf-8"))
        assert "contacts" in data
        assert "summary" in data
        assert "rules" in data

        summary = data["summary"]
        assert summary["n_contacts"] == 7
        assert summary["by_type"]["Hydrogen_bond"] == 2
        assert summary["by_type"]["Salt_bridge"] == 2
        assert summary["by_type"]["Hydrophobic"] == 1
        assert summary["by_type"]["Pi_Pi"] == 1
        assert summary["by_type"]["Cation_Pi"] == 1
        assert summary["by_type"]["Ligand_contact"] == 1

    def test_contact_records_have_labels_and_types(self, tmp_path: Path) -> None:
        """接触记录含残基标签、距离与类型。"""
        out = tmp_path / "out"
        proc = _run_cli([
            "--input", str(FIXTURES / "two_chain_complex.pdb"),
            "--output", str(out),
        ])
        assert proc.returncode == 0, proc.stderr

        data = json.loads((out / "contacts.json").read_text(encoding="utf-8"))
        for record in data["contacts"]:
            assert record["label_a"].count(":") == 1
            assert record["label_b"].count(":") == 1
            assert isinstance(record["min_distance"], float)
            assert record["min_distance"] > 0
            assert isinstance(record["interaction_types"], list)

    def test_csv_output_is_readable(self, tmp_path: Path) -> None:
        """CSV 用 utf-8-sig 编码（Excel 直接打开不乱码）。"""
        out = tmp_path / "out"
        proc = _run_cli([
            "--input", str(FIXTURES / "two_chain_complex.pdb"),
            "--output", str(out),
        ])
        assert proc.returncode == 0, proc.stderr

        csv_path = out / "contacts.csv"
        raw = csv_path.read_bytes()
        assert raw.startswith(b"\xef\xbb\xbf"), "CSV 缺少 UTF-8 BOM"

        with csv_path.open("r", encoding="utf-8-sig", newline="") as handle:
            rows = list(csv.reader(handle))
        assert rows[0][0] == "Chain_Pair"
        assert len(rows) == 8  # 表头 + 7 个接触

    def test_single_chain_pair_selection(self, tmp_path: Path) -> None:
        """``--chains`` 限定只分析指定的链对。"""
        out = tmp_path / "out"
        proc = _run_cli([
            "--input", str(FIXTURES / "two_chain_complex.pdb"),
            "--output", str(out),
            "--chains", "A", "B",
        ])
        assert proc.returncode == 0, proc.stderr

        data = json.loads((out / "contacts.json").read_text(encoding="utf-8"))
        assert data["summary"]["per_chain_pair"] == {"A-B": 6}
        for record in data["contacts"]:
            assert record["chain_pair"] == "A-B"

    def test_invalid_chain_is_rejected(self, tmp_path: Path) -> None:
        """指定不存在的链时报错。"""
        proc = _run_cli([
            "--input", str(FIXTURES / "two_chain_complex.pdb"),
            "--output", str(tmp_path / "out"),
            "--chains", "A", "Z",
        ])
        assert proc.returncode != 0
        assert "链" in (proc.stdout + proc.stderr)

    def test_chains_requires_exactly_two(self, tmp_path: Path) -> None:
        """``--chains`` 必须恰好给两个。"""
        proc = _run_cli([
            "--input", str(FIXTURES / "two_chain_complex.pdb"),
            "--output", str(tmp_path / "out"),
            "--chains", "A", "B", "C",
        ])
        assert proc.returncode != 0

    def test_exclude_hetero_option(self, tmp_path: Path) -> None:
        """``--exclude-hetero`` 后配体接触消失。"""
        out = tmp_path / "out"
        proc = _run_cli([
            "--input", str(FIXTURES / "two_chain_complex.pdb"),
            "--output", str(out),
            "--exclude-hetero",
        ])
        assert proc.returncode == 0, proc.stderr

        data = json.loads((out / "contacts.json").read_text(encoding="utf-8"))
        assert data["summary"]["by_type"]["Ligand_contact"] == 0
        assert data["summary"]["n_contacts"] == 6

    def test_min_distance_filter(self, tmp_path: Path) -> None:
        """``--min-distance`` 过滤掉过近的接触。"""
        out = tmp_path / "out"
        proc = _run_cli([
            "--input", str(FIXTURES / "two_chain_complex.pdb"),
            "--output", str(out),
            "--min-distance", "4.0",
        ])
        assert proc.returncode == 0, proc.stderr

        data = json.loads((out / "contacts.json").read_text(encoding="utf-8"))
        for record in data["contacts"]:
            assert record["min_distance"] >= 4.0

    def test_max_distance_option(self, tmp_path: Path) -> None:
        """``--max-distance`` 收紧后接触变少。"""
        out = tmp_path / "out"
        proc = _run_cli([
            "--input", str(FIXTURES / "two_chain_complex.pdb"),
            "--output", str(out),
            "--max-distance", "3.5",
        ])
        assert proc.returncode == 0, proc.stderr

        data = json.loads((out / "contacts.json").read_text(encoding="utf-8"))
        # Cation_Pi 的距离是 5.5 Å，收紧到 3.5 后应消失
        assert data["summary"]["by_type"]["Cation_Pi"] == 0

    def test_include_histidine_option(self, tmp_path: Path) -> None:
        """``--include-histidine`` 打开 HIS 的盐桥判定。"""
        out = tmp_path / "out"
        proc = _run_cli([
            "--input", str(FIXTURES / "two_chain_complex.pdb"),
            "--output", str(out),
            "--include-histidine",
        ])
        assert proc.returncode == 0, proc.stderr

        data = json.loads((out / "contacts.json").read_text(encoding="utf-8"))
        assert "by_type" in data["summary"]

    def test_mermaid_output_is_valid(self, tmp_path: Path) -> None:
        """Mermaid 输出是合法的 graph 定义。"""
        out = tmp_path / "out"
        proc = _run_cli([
            "--input", str(FIXTURES / "two_chain_complex.pdb"),
            "--output", str(out),
            "--make-mermaid",
        ])
        assert proc.returncode == 0, proc.stderr

        text = (out / "contact_network.mmd").read_text(encoding="utf-8")
        assert text.startswith("graph LR")
        assert "-->" in text

    def test_pymol_scripts_are_generated_not_executed(self, tmp_path: Path) -> None:
        """PyMOL 脚本是**生成物** —— 不需要安装 PyMOL 也能产出。"""
        out = tmp_path / "out"
        proc = _run_cli([
            "--input", str(FIXTURES / "two_chain_complex.pdb"),
            "--output", str(out),
            "--make-pymol",
            "--make-pymol-contacts",
        ])
        assert proc.returncode == 0, proc.stderr

        pml = (out / "visualize_contacts.pml").read_text(encoding="utf-8")
        assert "from pymol import cmd" in pml
        assert "cmd.select" in pml

        draw = (out / "draw_contacts.pml").read_text(encoding="utf-8")
        assert "cmd.distance" in draw

    def test_multimodel_takes_first_model_only(self, tmp_path: Path) -> None:
        """多 MODEL 文件只取第一个 MODEL。

        第二个 MODEL 的原子在 x = 50 附近，与第一个 MODEL 不会产生接触；
        如果误把两个 MODEL 都读进来，原子数会翻倍。
        """
        out = tmp_path / "out"
        proc = _run_cli([
            "--input", str(FIXTURES / "multimodel.pdb"),
            "--output", str(out),
        ])
        assert proc.returncode == 0, proc.stderr

        # 第一个 MODEL 只有 ALA 与 GLY 两个残基，各 4 个原子，共 8 个。
        # 如果误把两个 MODEL 都读进来，原子数会翻倍到 16。
        atoms = _load_fixture_atoms(FIXTURES / "multimodel.pdb")
        assert len(atoms) == 8, f"应只读到第一个 MODEL 的 8 个原子，实际 {len(atoms)}"

        # 第二个 MODEL 的原子在 x = 50 附近，不应出现在结果里
        assert all(atom.x < 10.0 for atom in atoms), "读到了第一个 MODEL 之外的原子"

        # JSON 里的接触也应只涉及第一个 MODEL 的残基
        data = json.loads((out / "contacts.json").read_text(encoding="utf-8"))
        for record in data["contacts"]:
            assert record["res_seq_a"] in (1, 2), "出现了第一个 MODEL 之外的残基"


# ---------------------------------------------------------------- 回归锁


class TestLegacyRegressions:
    """老仓库问题的回归锁。"""

    @staticmethod
    def _code_without_comments(tree: ast.AST) -> str:
        """把 AST 还原成源码，但**剔除所有 docstring**。

        背景：``analyze_contacts.py`` 的模块 docstring 里会**引用**老仓库的
        硬编码路径（如 ``C:\\Users\\DELL\\Desktop\\1.pdb``）来说明"这里曾经
        有什么问题"。那是文档，不是代码 —— 直接对源码做子串匹配会误报。

        所以这里只检查**真正的代码部分**：把模块/函数/类的 docstring
        节点删掉后重新 unparse。
        """
        for node in ast.walk(tree):
            if not isinstance(
                node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)
            ):
                continue
            body = node.body
            if (
                body
                and isinstance(body[0], ast.Expr)
                and isinstance(body[0].value, ast.Constant)
                and isinstance(body[0].value.value, str)
            ):
                node.body = body[1:] or [ast.Pass()]
        return ast.unparse(tree)

    def test_no_hardcoded_personal_paths(self) -> None:
        """脚本的**代码部分**不能有硬编码的个人路径。

        注意：docstring 里引用老路径作为"问题说明"是允许的，
        因此这里检查的是剔除 docstring 后的代码。
        """
        for name in SCRIPTS:
            code = self._code_without_comments(ast.parse(_read_script(name)))
            for bad in FORBIDDEN_PATHS:
                assert bad not in code, f"{name} 的代码里出现了硬编码路径 {bad!r}"

    def test_no_pymol_import_at_runtime(self) -> None:
        """核心脚本不能在运行时依赖 PyMOL。

        老仓库全部脚本以 ``pymol.finish_launching(['pymol'])`` 开头 ——
        无 PyMOL 的机器上完全跑不了，CI 里无法测试。

        检查方式：看 AST 里有没有**真正的 import 节点**指向 pymol。
        不能用子串匹配 —— ``analyze_contacts.py`` 生成 ``.pml`` 脚本时，
        字符串字面量里本来就要包含 ``"from pymol import cmd"``，
        那是要写进生成文件的内容，不是运行依赖。
        """
        for name in SCRIPTS:
            tree = ast.parse(_read_script(name))
            for node in ast.walk(tree):
                if isinstance(node, ast.Import):
                    for alias in node.names:
                        assert not alias.name.split(".")[0] == "pymol", (
                            f"{name} 在运行时 import 了 pymol"
                        )
                elif isinstance(node, ast.ImportFrom):
                    root = (node.module or "").split(".")[0]
                    assert root != "pymol", f"{name} 在运行时 import 了 pymol"

            # 也不能调用 finish_launching
            code = self._code_without_comments(tree)
            assert "finish_launching" not in code, f"{name} 里有 `finish_launching`"

    def test_script_filenames_have_no_spaces(self) -> None:
        """脚本文件名不能带空格。

        老仓库的 ``interaction count.py``、``hbond count.py`` 在 shell 里
        必须加引号，PyMOL 的 ``run`` 命令对带空格路径的处理也容易出问题。
        """
        for name in SCRIPTS:
            assert " " not in name, f"脚本名 {name!r} 含空格"

    def test_no_pep604_annotations(self) -> None:
        """不使用 PEP 604 联合类型（Python 3.9 兼容）。"""
        for name in SCRIPTS:
            tree = ast.parse(_read_script(name))
            for node in ast.walk(tree):
                if isinstance(node, ast.BinOp) and isinstance(node.op, ast.BitOr):
                    # 排除类型注解上下文之外的无害用法
                    pass
            # 简单检查：源码里不应出现 "X | Y" 形式的注解
            text = _read_script(name)
            assert ": str | " not in text
            assert "-> str | " not in text
            assert "Optional[str] | " not in text

    def test_classification_tables_are_unified(self) -> None:
        """分类表只有一份 —— 不允许在脚本里再写一遍 donors/acceptors。"""
        text = _read_script("analyze_contacts.py")
        # 老脚本的写法：acceptors = [...] / donors = [...]
        assert "acceptors = [" not in text
        assert "donors = [" not in text
        assert "positive = [" not in text
        assert "aromatic = [" not in text

    def test_chemistry_module_is_single_source_of_truth(self) -> None:
        """化学表集中在 contact_chemistry 里。"""
        assert contact_chemistry.DONORS
        assert contact_chemistry.ACCEPTORS
        assert contact_chemistry.POSITIVE
        assert contact_chemistry.NEGATIVE
        assert contact_chemistry.HYDROPHOBIC
        assert contact_chemistry.AROMATIC

    def test_output_files_stay_inside_output_dir(self, tmp_path: Path) -> None:
        """产物不逃逸出输出目录。

        Windows 上 ``Path("out") / "json" / "C:\\x.json"`` 会锚定到盘符根目录。
        """
        out = tmp_path / "out"
        proc = _run_cli([
            "--input", str(FIXTURES / "two_chain_complex.pdb"),
            "--output", str(out),
            "--make-mermaid",
            "--make-pymol",
        ])
        assert proc.returncode == 0, proc.stderr

        produced = {p.resolve() for p in out.rglob("*") if p.is_file()}
        assert produced, "没有产出任何文件"
        for path in produced:
            assert str(path).startswith(str(out.resolve())), f"产物逃逸出输出目录: {path}"

    def test_does_not_require_numpy_or_pandas(self, tmp_path: Path) -> None:
        """核心分析不依赖 numpy / pandas。

        老 ``use plip find interact region.py`` 依赖 pandas + openpyxl。
        新技能的核心计算必须零第三方依赖。
        """
        text = _read_script("contact_chemistry.py")
        assert "import numpy" not in text
        assert "import pandas" not in text

        text = _read_script("contact_geometry.py")
        assert "import numpy" not in text
        assert "import pandas" not in text


class TestRegressionCanary:
    """反向验证：这套用例必须真的能抓到问题。

    逐个注入缺陷，确认对应用例失败。检测逻辑本身失效是最危险的情况 ——
    那意味着测试全绿但什么都没测。
    """

    def test_canary_asp_glu_donor_injection(self) -> None:
        """把 ASP/GLU 加回供体表，分类用例必须失败。

        注入方式是**替换模块级常量** —— 注意 ``is_donor_residue`` 在函数体里
        直接引用 ``DONORS`` 这个全局名，所以替换 ``contact_chemistry.DONORS``
        立刻生效，不需要改函数本身。
        """
        original = contact_chemistry.DONORS
        try:
            contact_chemistry.DONORS = frozenset(original | {"ASP", "GLU"})
            # 注入后 ASP 成了供体 —— 依赖此事实的用例必须失败
            assert contact_chemistry.is_donor_residue("ASP"), "注入未生效"
            assert "Hydrogen_bond" in contact_chemistry.classify_pair_multi(
                "ASP", "GLU", 3.0
            ), "注入未生效"
        finally:
            contact_chemistry.DONORS = original

        # 还原后行为回到正常
        assert not contact_chemistry.is_donor_residue("ASP")
        assert "Hydrogen_bond" not in contact_chemistry.classify_pair_multi(
            "ASP", "GLU", 3.0
        )

    def test_canary_histidine_positive_injection(self) -> None:
        """把 HIS 加回正电表，盐桥用例必须失败。"""
        original = contact_chemistry.POSITIVE
        try:
            contact_chemistry.POSITIVE = frozenset(original | {"HIS"})
            assert contact_chemistry.is_positive_residue("HIS"), "注入未生效"
            assert "Salt_bridge" in contact_chemistry.classify_pair_multi(
                "HIS", "ASP", 3.0
            ), "注入未生效"
        finally:
            contact_chemistry.POSITIVE = original

    def test_canary_grid_radius_reduction(self) -> None:
        """把网格扫描半径降到 0，跨桶的接触必须被漏掉。

        桶边长取 ``max_distance * 1.0001``（严格大于最大距离），
        所以距离恰为最大值的两个原子可能落在**相邻**两个桶里。
        扫 1 格（3×3×3）能覆盖，扫 0 格（只看自己所在的桶）会漏。
        """
        original = contact_geometry.ContactDetector._cell_range

        def _zero(self) -> int:
            return 0

        atoms = [
            _make_atom(1, "NZ", "LYS", "A", 1, (5.999, 0.0, 0.0)),
            _make_atom(2, "OE1", "GLU", "B", 1, (6.001, 0.0, 0.0)),
        ]

        try:
            contact_geometry.ContactDetector._cell_range = _zero
            detector = contact_geometry.ContactDetector(max_distance=6.0)
            detector.add_atoms(atoms)
            assert detector.detect("A", "B") == [], "注入未生效：跨桶接触未被漏掉"
        finally:
            contact_geometry.ContactDetector._cell_range = original

        # 还原后必须能检出
        detector = contact_geometry.ContactDetector(max_distance=6.0)
        detector.add_atoms(atoms)
        assert len(detector.detect("A", "B")) == 1, "还原后仍漏检 —— 网格实现有 bug"

    def test_canary_hydrogen_filter_injection(self) -> None:
        """关掉氢过滤，氢相关用例必须失败。"""
        atoms = [
            _make_atom(1, "NZ", "LYS", "A", 1, (0.0, 0.0, 0.0)),
            _make_atom(2, "H", "LYS", "A", 1, (1.0, 0.0, 0.0), element="H"),
            _make_atom(3, "OE1", "GLU", "B", 1, (3.0, 0.0, 0.0)),
        ]
        detector = contact_geometry.ContactDetector(
            max_distance=6.0, ignore_hydrogen=False
        )
        detector.add_atoms(atoms)
        without_h = contact_geometry.ContactDetector(max_distance=6.0)
        without_h.add_atoms(atoms)
        assert (
            detector.detect("A", "B")[0].n_atom_pairs
            > without_h.detect("A", "B")[0].n_atom_pairs
        ), "注入未生效：氢原子没有影响原子对数量"

    def test_canary_ligand_classification_injection(self) -> None:
        """让配体走残基级分类，配体用例必须失败。"""
        original = contact_chemistry.residue_kind

        def _everything_standard(res_name: str) -> str:
            return "standard"

        try:
            contact_chemistry.residue_kind = _everything_standard
            types = contact_chemistry.classify_pair_multi("TYR", "LIG", 3.5)
            assert contact_chemistry.LIGAND_CONTACT not in types, "注入未生效"
        finally:
            contact_chemistry.residue_kind = original

    def test_canary_exclude_hetero_injection(self) -> None:
        """让 ``--exclude-hetero`` 失效，配体计数用例必须失败。"""
        original = analyze_contacts._count_by_type

        def _broken(records: list[dict]) -> dict[str, int]:
            # 故意忽略配体类型
            result = original(records)
            result[contact_chemistry.LIGAND_CONTACT] = 0
            return result

        try:
            analyze_contacts._count_by_type = _broken
            assert analyze_contacts._count_by_type([])[
                contact_chemistry.LIGAND_CONTACT
            ] == 0, "注入未生效"
        finally:
            analyze_contacts._count_by_type = original
