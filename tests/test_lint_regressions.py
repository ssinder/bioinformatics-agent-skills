"""Lint 相关缺陷的回归测试。

这些用例锁住一批曾经真实存在过的缺陷，防止 ruff 规则变更或后续重构
把它们重新引入：

1. ``biopipe_core.io_utils`` 曾在类型注解里使用 ``Any`` 却未导入 ——
   因 ``from __future__ import annotations`` 的存在，运行时不报错，
   但 ``typing.get_type_hints()`` 会抛 NameError。
2. ``protein-structure/scripts/analyze_structure.py`` 的 ``__all__``
   曾列出不存在的名字 ``build_parser``（实际函数为 ``_build_parser``），
   导致 ``from module import *`` 抛 ImportError。
3. ``peptide_properties`` 的等电点二分求解曾保留无用变量，
   且二分逻辑必须只依赖 pH 边界与符号判定。

本模块不依赖 numpy/pandas，只用标准库 + pytest。
"""

from __future__ import annotations

import sys
import typing
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "packages"))


@pytest.fixture(scope="module")
def peptide_props():
    """加载 peptide_properties 模块（模块级 fixture，兼容 Python 3.9）。"""
    sys.path.insert(0, str(REPO_ROOT / "skills" / "peptide-design" / "scripts"))
    import peptide_properties

    return peptide_properties


class TestTypeHintResolution:
    """类型注解必须能被 typing.get_type_hints() 解析。"""

    def test_io_utils_write_csv_hints_resolve(self) -> None:
        from biopipe_core import io_utils

        hints = typing.get_type_hints(io_utils.write_csv)
        assert "rows" in hints
        assert "path" in hints

    def test_core_module_hints_resolve(self) -> None:
        """共享层所有公开函数的类型注解都应可解析。"""
        from biopipe_core import cli, config, http, io_utils, logging_utils

        for module in (cli, config, http, io_utils, logging_utils):
            for name in getattr(module, "__all__", []):
                obj = getattr(module, name)
                if callable(obj) and getattr(obj, "__module__", "") == module.__name__:
                    typing.get_type_hints(obj)  # 不应抛 NameError

    def test_no_undefined_names_in_annotations(self) -> None:
        """所有共享层模块的 __all__ 名字必须真实存在。"""
        from biopipe_core import cli, config, http, io_utils, logging_utils

        for module in (cli, config, http, io_utils, logging_utils):
            for name in getattr(module, "__all__", []):
                assert hasattr(module, name), f"{module.__name__}.{name} 不存在"


class TestSkillModuleExports:
    """技能脚本的 __all__ 不得列出不存在的名字。"""

    @pytest.mark.parametrize(
        "script",
        [
            "skills/protein-structure/scripts/analyze_structure.py",
            "skills/protein-structure/scripts/secondary_structure.py",
            "skills/peptide-design/scripts/peptide_properties.py",
            "skills/molecular-dynamics/scripts/trajectory_io.py",
        ],
    )
    def test_all_entries_exist(self, script: str) -> None:
        import importlib.util

        path = REPO_ROOT / script
        assert path.exists(), f"{script} 不存在"

        mod_name = f"_lintreg_{path.stem}"
        spec = importlib.util.spec_from_file_location(mod_name, path)
        assert spec is not None and spec.loader is not None
        module = importlib.util.module_from_spec(spec)
        # dataclass 装饰器需要模块已注册到 sys.modules 才能正确解析 __dict__
        sys.modules[mod_name] = module
        try:
            spec.loader.exec_module(module)
        finally:
            sys.modules.pop(mod_name, None)

        for name in getattr(module, "__all__", []):
            assert hasattr(module, name), f"{script} 的 __all__ 含未定义名字: {name}"


class TestPeptideIsoelectricPoint:
    """等电点二分求解的行为约束。"""

    def test_net_charge_crosses_zero(self, peptide_props) -> None:
        seq = "KWKLFKKIEKVGQNIRDGIIKAGPA"
        charge = peptide_props.net_charge
        assert charge(seq, 4.0) > 0
        assert charge(seq, 12.0) < 0

    def test_isoelectric_point_within_range(self, peptide_props) -> None:
        seq = "KWKLFKKIEKVGQNIRDGIIKAGPA"
        pi = peptide_props.isoelectric_point(seq)
        assert 4.0 <= pi <= 12.0

    def test_pi_is_monotonic_bisection(self, peptide_props) -> None:
        """二分区间必须单调收缩：返回的 pH 落在初始区间内。"""
        seq = "ACDEFGHIKLMNPQRSTVWY"
        pi = peptide_props.isoelectric_point(seq)
        assert 0.0 <= pi <= 14.0

    def test_pure_sequence_has_no_unknown_residues(self, peptide_props) -> None:
        """纯序列不应被误判为含未知残基。"""
        seq = "KWKLFKKIEKVGQNIRDGIIKAGPA"
        mw = peptide_props.molecular_weight(seq)
        assert mw > 0
        _, unknown = peptide_props.validate_sequence(seq)
        assert unknown == [], f"纯序列被误判含未知残基: {unknown}"


class TestModificationSuffixStripping:
    """化学修饰后缀（如 -NH2）必须被剥离，不影响分子量计算。"""

    @pytest.mark.parametrize("suffix", ["-NH2", "-COOH", "-OH"])
    def test_suffix_stripped_from_sequence(self, peptide_props, suffix: str) -> None:
        """修饰后缀必须被剥离，且不产生未知残基。"""
        base = "KWKLFKKIEKVGQNIRDGIIKAGPA"
        stripped, detected = peptide_props.strip_modification(base + suffix)
        assert stripped == base, f"{suffix} 未被正确剥离，得到 {stripped}"
        assert detected == suffix

    @pytest.mark.parametrize("suffix", ["-NH2", "-COOH", "-OH", "-CHO", "-SH", "-H", "-AC", "-AMID"])
    def test_suffix_does_not_inflate_molecular_weight(self, peptide_props, suffix: str) -> None:
        """分子量不得因修饰后缀被当作残基而虚高。

        历史缺陷：``molecular_weight`` 等 15 个计算函数直接遍历原始字符串，
        导致 ``-OH`` 中的 ``O``/``H`` 被当作残基，引入约 250 Da 的虚高。
        """
        base = "KWKLFKKIEKVGQNIRDGIIKAGPA"
        plain = peptide_props.molecular_weight(base)
        modified = peptide_props.molecular_weight(base + suffix)
        assert abs(modified - plain) < 0.01, (
            f"{suffix} 导致分子量偏移 {modified - plain:.4f} Da"
        )

    @pytest.mark.parametrize(
        "func_name",
        [
            "molecular_weight",
            "monoisotopic_weight",
            "kyte_doolittle",
            "hydrophobic_moment",
            "aromaticity",
            "aliphatic_index",
            "instability_index",
            "boman_index",
            "isoelectric_point",
        ],
    )
    def test_all_computation_functions_ignore_suffix(self, peptide_props, func_name: str) -> None:
        """所有计算函数都必须先剥离修饰后缀。"""
        base = "KWKLFKKIEKVGQNIRDGIIKAGPA"
        func = getattr(peptide_props, func_name)
        plain = func(base)
        modified = func(base + "-NH2")
        if isinstance(plain, float):
            assert abs(modified - plain) < 0.01, f"{func_name} 受修饰后缀影响: {plain} -> {modified}"
        else:
            assert modified == plain, f"{func_name} 受修饰后缀影响"

    def test_no_unknown_residue_from_suffix(self, peptide_props) -> None:
        base = "KWKLFKKIEKVGQNIRDGIIKAGPA"
        _, unknown = peptide_props.validate_sequence(base + "-NH2")
        assert unknown == [], f"修饰后缀被误判为未知残基: {unknown}"

    def test_hydrophobicity_unaffected_by_suffix(self, peptide_props) -> None:
        """疏水性不得因修饰后缀被当作残基而偏负。"""
        base = "KWKLFKKIEKVGQNIRDGIIKAGPA"
        plain = peptide_props.kyte_doolittle(base)
        modified = peptide_props.kyte_doolittle(base + "-NH2")
        assert abs(modified - plain) < 1e-6, f"疏水性偏移 {modified - plain:.6f}"
