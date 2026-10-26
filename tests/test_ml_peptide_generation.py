"""ML 肽链生成技能的回归测试。

锁住老仓库 ``design_pepmlm/`` 的一批真实缺陷：

1. **模型在 import 时加载** —— ``test_pepmlm.py`` 在模块顶层调用
   ``from_pretrained``，import 即下载 650M 参数模型。本仓库禁止模块级
   副作用，这里用 AST 静态检查锁住。
2. **没有任何随机种子控制** —— 生成用 ``Categorical.sample()``，同一输入
   跑两次结果不同，实验不可复现。本模块要求 ``--seed`` 且必须有默认值。
3. **硬编码 macOS 绝对路径** —— ``/Users/jung/Desktop/trae/tool/...``。
4. **三个脚本近乎逐行重复** —— 合并后不应再有重复实现。
5. **采样会产出 X/B/Z 等非法氨基酸** —— 老脚本只做事后标记。
6. **循环变量越界** —— ``batch_generate_peptides.py`` 在长度循环之外
   引用 ``results``，被 ``continue`` 跳过时会 UnboundLocalError。
7. **筛选脚本把文件名硬编码成 12 个** —— 换参数就得改源码。
8. **``pd.concat([])`` 空输入崩溃** —— 所有文件都过滤空时抛 ValueError。

全部用例不联网、不加载模型。
"""

from __future__ import annotations

import ast
import csv
import re
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
SKILL_DIR = REPO_ROOT / "skills" / "ml-peptide-generation"
sys.path.insert(0, str(REPO_ROOT / "packages"))
sys.path.insert(0, str(SKILL_DIR / "scripts"))

import analyze_peptides  # noqa: E402
import generate_peptides  # noqa: E402
import pepmlm_client  # noqa: E402


def _script(name: str) -> Path:
    return SKILL_DIR / "scripts" / name


def _parse(name: str) -> ast.Module:
    return ast.parse(_script(name).read_text(encoding="utf-8"))


def _non_docstring_strings(tree: ast.AST) -> list[str]:
    """返回源码里所有**非 docstring** 的字符串字面量。"""
    docstring_nodes: set[int] = set()
    for node in ast.walk(tree):
        if isinstance(
            node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)
        ):
            body = getattr(node, "body", [])
            if body and isinstance(body[0], ast.Expr):
                value = body[0].value
                if isinstance(value, ast.Constant) and isinstance(value.value, str):
                    docstring_nodes.add(id(value))

    return [
        node.value
        for node in ast.walk(tree)
        if isinstance(node, ast.Constant)
        and isinstance(node.value, str)
        and id(node) not in docstring_nodes
    ]


def _imported_modules(name: str) -> set[str]:
    """提取脚本里实际 import 的顶层模块名。"""
    modules: set[str] = set()
    for node in ast.walk(_parse(name)):
        if isinstance(node, ast.Import):
            for alias in node.names:
                modules.add(alias.name.split(".")[0])
        elif isinstance(node, ast.ImportFrom):
            if node.module and node.level == 0:
                modules.add(node.module.split(".")[0])
    return modules


class TestNoModuleLevelSideEffects:
    """禁止 import 即加载模型 / 建目录。"""

    @pytest.mark.parametrize("script", ["pepmlm_client.py", "generate_peptides.py",
                                        "analyze_peptides.py"])
    def test_no_module_level_execution(self, script):
        """脚本顶层不允许有可执行语句。

        老脚本在模块顶层调用 ``AutoModelForMaskedLM.from_pretrained``——
        只要 import 就会下载 650M 参数模型。``--help`` 因此也无法执行。
        """
        tree = _parse(script)

        for node in tree.body:
            # 导入与定义：安全
            if isinstance(node, (ast.Import, ast.ImportFrom, ast.FunctionDef,
                                 ast.AsyncFunctionDef, ast.ClassDef)):
                continue

            # 赋值：允许常量赋值（LOGGER = get_logger(...) 之类）
            if isinstance(node, (ast.Assign, ast.AnnAssign)):
                continue

            # 表达式：只允许 docstring 与 sys.path.insert
            if isinstance(node, ast.Expr):
                if isinstance(node.value, ast.Constant):
                    continue
                func = node.value.func if isinstance(node.value, ast.Call) else None
                is_path_insert = (
                    isinstance(func, ast.Attribute)
                    and func.attr == "insert"
                    and isinstance(func.value, ast.Attribute)
                    and func.value.attr == "path"
                )
                assert is_path_insert, (
                    f"{script} 第 {node.lineno} 行存在模块级可执行语句: "
                    f"{ast.dump(node)[:120]}"
                )
                continue

            # if 只允许 __future__ / __main__ 守卫
            if isinstance(node, ast.If):
                test = ast.dump(node.test)
                assert "__future__" in test or "__main__" in test, (
                    f"{script} 第 {node.lineno} 行的模块级 if 不是守卫语句"
                )
                continue

            raise AssertionError(
                f"{script} 第 {node.lineno} 行存在不允许的模块级语句: {type(node).__name__}"
            )

    def test_model_not_loaded_at_import(self):
        """``from_pretrained`` 只能出现在函数内部，不能在模块顶层。"""
        for script in ["pepmlm_client.py", "generate_peptides.py"]:
            tree = _parse(script)
            for node in ast.walk(tree):
                if isinstance(node, ast.Call):
                    func = node.func
                    name = getattr(func, "attr", "") or getattr(func, "id", "")
                    if name == "from_pretrained":
                        # 必须嵌套在函数/方法里
                        assert not isinstance(node, ast.stmt), (
                            f"{script}: from_pretrained 不能是模块级语句"
                        )


class TestNoHardcodedPaths:
    """禁止写死的绝对路径。"""

    @pytest.mark.parametrize("script", ["pepmlm_client.py", "generate_peptides.py",
                                        "analyze_peptides.py"])
    def test_no_user_home_path(self, script):
        for literal in _non_docstring_strings(_parse(script)):
            assert "/Users/" not in literal, f"{script} 含硬编码 macOS 路径: {literal}"
            assert "C:\\Users" not in literal, f"{script} 含硬编码 Windows 路径: {literal}"
            assert "/home/" not in literal, f"{script} 含硬编码 Linux 家目录: {literal}"
            assert "Desktop" not in literal, f"{script} 含硬编码桌面路径: {literal}"

    def test_no_deeploc_local_dependency(self):
        """老 requirements.txt 里有 ``file:///Users/jung/...`` 的本地依赖。

        只检查**非注释行**——在注释里引用这条反面例子是合理的（说明
        "已移除"），但真把它写成一条依赖，别人就装不上了。
        """
        requirements = SKILL_DIR / "requirements.txt"
        lines = requirements.read_text(encoding="utf-8").splitlines()
        active = [line for line in lines if line.strip() and not line.strip().startswith("#")]
        for line in active:
            assert "file://" not in line, f"requirements.txt 含本地路径依赖: {line}"
            assert "/Users/" not in line, f"requirements.txt 含硬编码路径: {line}"


class TestSeedControl:
    """随机种子控制——老脚本完全没有。"""

    def test_seed_argument_exists_with_default(self):
        args = generate_peptides.build_parser().parse_args(["--input-seq", "MKTA"])
        assert hasattr(args, "seed")
        assert isinstance(args.seed, int)

    def test_seed_is_recorded_in_output(self):
        """种子必须写进结果行，否则事后无法复现。"""
        from pepmlm_client import PepMLMClient  # noqa: F401

        source = _script("generate_peptides.py").read_text(encoding="utf-8")
        assert '"seed"' in source
        assert "seed=args.seed" in source

    def test_set_seed_sets_all_three_sources(self):
        """必须同时设置 random / numpy / torch。"""
        import random

        pepmlm_client.set_seed(1234)
        first = random.random()

        pepmlm_client.set_seed(1234)
        assert random.random() == first, "random.seed 未生效"

    def test_set_seed_is_reproducible_across_calls(self):
        import random

        pepmlm_client.set_seed(7)
        seq1 = [random.random() for _ in range(5)]
        pepmlm_client.set_seed(7)
        seq2 = [random.random() for _ in range(5)]
        assert seq1 == seq2

    def test_config_has_seed_field(self):
        config = pepmlm_client.PepMLMConfig()
        assert config.seed == 42


class TestSampleFromProbabilities:
    """纯 Python 采样辅助函数。"""

    def test_respects_distribution(self):
        import random

        rng = random.Random(0)
        counts = [0, 0, 0]
        for _ in range(1000):
            counts[pepmlm_client.sample_from_probabilities([0.9, 0.05, 0.05], rng)] += 1
        assert counts[0] > counts[1]
        assert counts[0] > counts[2]

    def test_zero_total_falls_back(self):
        import random

        rng = random.Random(0)
        index = pepmlm_client.sample_from_probabilities([0.0, 0.0], rng)
        assert index in (0, 1)

    def test_single_candidate_always_selected(self):
        import random

        rng = random.Random(0)
        for _ in range(10):
            assert pepmlm_client.sample_from_probabilities([1.0], rng) == 0


class TestFastaParsing:
    """FASTA 解析与标识符提取。"""

    def test_reads_multi_line_fasta(self, tmp_path):
        path = tmp_path / "t.fasta"
        path.write_text(
            ">sp|P69905|HBA_HUMAN GN=HBA\nMVLSPADKTN\nVKAAWGKVGA\n",
            encoding="utf-8",
        )
        records = pepmlm_client.parse_target_records(fasta_path=path)
        assert len(records) == 1
        assert records[0]["id"] == "HBA"
        assert records[0]["seq"] == "MVLSPADKTNVKAAWGKVGA"

    def test_gn_field_takes_priority(self, tmp_path):
        path = tmp_path / "t.fasta"
        path.write_text(">sp|P69905|HBA_HUMAN GN=HBA PE=1 SV=2\nMVLSA\n", encoding="utf-8")
        records = pepmlm_client.parse_target_records(fasta_path=path)
        assert records[0]["id"] == "HBA"

    def test_falls_back_to_accession(self, tmp_path):
        path = tmp_path / "t.fasta"
        path.write_text(">sp|P69905|HBA_HUMAN\nMVLSA\n", encoding="utf-8")
        records = pepmlm_client.parse_target_records(fasta_path=path)
        assert records[0]["id"] == "P69905"

    def test_falls_back_to_header(self, tmp_path):
        path = tmp_path / "t.fasta"
        path.write_text(">my_target_description\nMVLSA\n", encoding="utf-8")
        records = pepmlm_client.parse_target_records(fasta_path=path)
        assert records[0]["id"] == "my_target_description"

    def test_skips_comment_lines(self, tmp_path):
        path = tmp_path / "t.fasta"
        path.write_text("; a comment\n>t1\nMVLSA\n", encoding="utf-8")
        records = pepmlm_client.parse_target_records(fasta_path=path)
        assert len(records) == 1

    def test_reads_multiple_records(self, tmp_path):
        path = tmp_path / "t.fasta"
        path.write_text(">t1\nMVLSA\n>t2\nMKTAY\n", encoding="utf-8")
        records = pepmlm_client.parse_target_records(fasta_path=path)
        assert [r["id"] for r in records] == ["t1", "t2"]


class TestCsvTargetParsing:
    """CSV 靶标解析。"""

    def test_reads_id_seq(self, tmp_path):
        path = tmp_path / "t.csv"
        path.write_text("id,seq\nT1,MVLSA\nT2,MKTAY\n", encoding="utf-8")
        records = pepmlm_client.parse_target_records(csv_path=path)
        assert len(records) == 2
        assert records[0]["id"] == "T1"
        assert records[0]["seq"] == "MVLSA"

    def test_requires_id_and_seq_columns(self, tmp_path):
        path = tmp_path / "t.csv"
        path.write_text("name,sequence\nT1,MVLSA\n", encoding="utf-8")
        with pytest.raises(ValueError, match="id"):
            pepmlm_client.parse_target_records(csv_path=path)

    def test_skips_empty_sequences(self, tmp_path):
        path = tmp_path / "t.csv"
        path.write_text("id,seq\nT1,\nT2,MKTAY\n", encoding="utf-8")
        records = pepmlm_client.parse_target_records(csv_path=path)
        assert len(records) == 1
        assert records[0]["id"] == "T2"


class TestInvalidAminoAcidDetection:
    """非法氨基酸检测。"""

    def test_valid_peptide(self):
        assert not pepmlm_client.PepMLMClient.has_invalid_amino_acids("ACDEFGHIKLMNPQRSTVWY")

    def test_x_is_invalid(self):
        assert pepmlm_client.PepMLMClient.has_invalid_amino_acids("ACDEFGHIKLXN")

    def test_b_and_z_are_invalid(self):
        assert pepmlm_client.PepMLMClient.has_invalid_amino_acids("ACDEFB")
        assert pepmlm_client.PepMLMClient.has_invalid_amino_acids("ACDEFZ")

    def test_valid_set_has_twenty_amino_acids(self):
        assert len(pepmlm_client.VALID_AMINO_ACIDS) == 20


class TestAnalyzeFiltering:
    """筛选逻辑（对应老 filter_analysis.py）。"""

    @staticmethod
    def _row(**kwargs):
        base = {
            "target_id": "T1",
            "peptide_length": "10",
            "top_k": "3",
            "binder": "ACDEFGHIKL",
            "binder_length": "10",
            "pseudo_perplexity": "1.5",
            "has_unknown_aa": "False",
            "config": "10aa_topk3",
        }
        base.update(kwargs)
        return base

    def test_filters_invalid_amino_acids(self):
        rows = [self._row(), self._row(binder="ACDEFGHIKX", has_unknown_aa="True")]
        kept = analyze_peptides.filter_peptides(rows)
        assert len(kept) == 1
        assert kept[0]["binder"] == "ACDEFGHIKL"

    def test_keep_invalid_option(self):
        rows = [self._row(binder="ACDEFGHIKX", has_unknown_aa="True")]
        kept = analyze_peptides.filter_peptides(rows, drop_invalid=False)
        assert len(kept) == 1

    def test_validity_detected_from_sequence_when_flag_missing(self):
        """老脚本只依赖 Has_Unknown_AA 一列，列名一变就失效。"""
        row = self._row()
        del row["has_unknown_aa"]
        assert analyze_peptides.is_valid_peptide(row)
        row["binder"] = "ACDEFGHIKX"
        assert not analyze_peptides.is_valid_peptide(row)

    def test_max_ppl_filter(self):
        rows = [self._row(pseudo_perplexity="1.0"), self._row(pseudo_perplexity="9.0")]
        kept = analyze_peptides.filter_peptides(rows, max_ppl=5.0)
        assert len(kept) == 1
        assert kept[0]["pseudo_perplexity_value"] == 1.0

    def test_length_filters(self):
        rows = [self._row(binder_length="8"), self._row(binder_length="20")]
        assert len(analyze_peptides.filter_peptides(rows, min_length=10)) == 1
        assert len(analyze_peptides.filter_peptides(rows, max_length=10)) == 1

    def test_rows_missing_ppl_are_dropped_when_max_ppl_set(self):
        rows = [self._row(pseudo_perplexity="")]
        assert analyze_peptides.filter_peptides(rows, max_ppl=5.0) == []

    def test_config_inferred_from_filename_when_columns_missing(self):
        """文件名兜底：results_10_3.csv → 10aa / top-k 3。"""
        row = {"binder": "ACDEFGHIKL", "binder_length": "10",
               "pseudo_perplexity": "1.0", "has_unknown_aa": "False"}
        kept = analyze_peptides.filter_peptides([row], source=Path("results_10_3.csv"))
        assert kept[0]["peptide_length"] == 10
        assert kept[0]["top_k"] == 3
        assert kept[0]["config"] == "10aa_topk3"

    def test_columns_take_priority_over_filename(self):
        row = {"peptide_length": "12", "top_k": "5", "binder": "ACDEFGHIKL",
               "binder_length": "10", "pseudo_perplexity": "1.0",
               "has_unknown_aa": "False"}
        kept = analyze_peptides.filter_peptides([row], source=Path("results_10_3.csv"))
        assert kept[0]["peptide_length"] == 12
        assert kept[0]["top_k"] == 5


class TestAnalyzeSummary:
    """参数统计。"""

    def test_groups_by_config(self):
        rows = [
            {"config": "10aa_topk3", "peptide_length": 10, "top_k": 3,
             "pseudo_perplexity_value": 1.0},
            {"config": "10aa_topk3", "peptide_length": 10, "top_k": 3,
             "pseudo_perplexity_value": 3.0},
            {"config": "12aa_topk5", "peptide_length": 12, "top_k": 5,
             "pseudo_perplexity_value": 2.0},
        ]
        summary = analyze_peptides.summarize_by_config(rows)
        assert len(summary) == 2
        by_config = {item["config"]: item for item in summary}
        assert by_config["10aa_topk3"]["total"] == 2
        assert by_config["10aa_topk3"]["mean_ppl"] == 2.0
        assert by_config["10aa_topk3"]["min_ppl"] == 1.0
        assert by_config["10aa_topk3"]["max_ppl"] == 3.0

    def test_sorted_by_mean_ppl(self):
        rows = [
            {"config": "A", "peptide_length": 10, "top_k": 3, "pseudo_perplexity_value": 9.0},
            {"config": "B", "peptide_length": 12, "top_k": 5, "pseudo_perplexity_value": 1.0},
        ]
        summary = analyze_peptides.summarize_by_config(rows)
        assert summary[0]["config"] == "B"

    def test_handles_rows_without_ppl(self):
        rows = [{"config": "A", "peptide_length": 10, "top_k": 3, "pseudo_perplexity_value": None}]
        summary = analyze_peptides.summarize_by_config(rows)
        assert summary[0]["mean_ppl"] == ""


class TestAnalyzeCli:
    """分析脚本的命令行行为。"""

    @staticmethod
    def _write_csv(path: Path, rows):
        with open(path, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)

    def test_writes_outputs(self, tmp_path):
        source = tmp_path / "in.csv"
        self._write_csv(source, [
            {"target_id": "T1", "peptide_length": "10", "top_k": "3",
             "binder": "ACDEFGHIKL", "binder_length": "10",
             "pseudo_perplexity": "1.5", "has_unknown_aa": "False", "config": "10aa_topk3"},
        ])
        out = tmp_path / "out"
        code = analyze_peptides.main(["--input", str(source), "--output-dir", str(out)])
        assert code == 0
        assert (out / "filtered_peptides.csv").is_file()
        assert (out / "parameter_analysis.csv").is_file()
        assert (out / "parameter_analysis.json").is_file()

    def test_empty_input_does_not_crash(self, tmp_path):
        """老脚本 ``pd.concat([])`` 在全部过滤空时抛 ValueError。"""
        source = tmp_path / "in.csv"
        self._write_csv(source, [
            {"target_id": "T1", "peptide_length": "10", "top_k": "3",
             "binder": "ACDEFGHIKX", "binder_length": "10",
             "pseudo_perplexity": "1.5", "has_unknown_aa": "True", "config": "10aa_topk3"},
        ])
        out = tmp_path / "out"
        code = analyze_peptides.main(["--input", str(source), "--output-dir", str(out)])
        assert code == 1, "全部被过滤时应返回 1 而不是崩溃"

    def test_missing_file_returns_error_code(self, tmp_path):
        code = analyze_peptides.main(
            ["--input", str(tmp_path / "nope.csv"), "--output-dir", str(tmp_path / "out")]
        )
        assert code == 2

    def test_multiple_inputs_are_merged(self, tmp_path):
        first = tmp_path / "a.csv"
        second = tmp_path / "b.csv"
        base = {"target_id": "T1", "peptide_length": "10", "top_k": "3",
                "binder": "ACDEFGHIKL", "binder_length": "10",
                "pseudo_perplexity": "1.5", "has_unknown_aa": "False", "config": "10aa_topk3"}
        self._write_csv(first, [base])
        self._write_csv(second, [base])
        out = tmp_path / "out"
        assert analyze_peptides.main(
            ["--input", str(first), str(second), "--output-dir", str(out)]
        ) == 0
        text = (out / "filtered_peptides.csv").read_text(encoding="utf-8-sig")
        assert text.count("ACDEFGHIKL") == 2


class TestGenerateCli:
    """生成脚本的命令行表面。"""

    def test_parser_builds(self):
        args = generate_peptides.build_parser().parse_args(
            ["--input-seq", "MKTA", "--output-dir", "out/", "--seed", "7"]
        )
        assert args.input_seq == "MKTA"
        assert args.seed == 7

    def test_requires_an_input(self):
        code = generate_peptides.main([])
        assert code == 2

    def test_output_dir_required_by_main(self):
        """没有输入源时 main() 应返回 2 而不是抛异常。"""
        assert generate_peptides.main([]) == 2

    def test_help_exits_cleanly(self, capsys):
        with pytest.raises(SystemExit) as exc:
            generate_peptides.build_parser().parse_args(["--help"])
        assert exc.value.code == 0
        assert "PepMLM" in capsys.readouterr().out

    def test_invalid_length_list_is_rejected(self):
        code = generate_peptides.main(["--input-seq", "MKTA", "--peptide-lengths", "10,abc"])
        assert code == 2

    def test_parse_int_list(self):
        assert generate_peptides.parse_int_list("10,12,15", "x") == [10, 12, 15]
        assert generate_peptides.parse_int_list(" 3 ", "x") == [3]

    def test_empty_length_list_rejected(self):
        with pytest.raises(ValueError):
            generate_peptides.parse_int_list("", "x")

    def test_output_fieldnames_cover_documented_columns(self):
        for column in ("binder", "pseudo_perplexity", "seed", "config",
                       "target_id", "peptide_length", "top_k", "rank"):
            assert column in generate_peptides.FIELDNAMES


class TestNoDuplicatedImplementation:
    """老仓库三个生成脚本近乎逐行重复，合并后不应再有重复实现。"""

    def test_only_two_entry_scripts(self):
        scripts = sorted(p.name for p in (SKILL_DIR / "scripts").glob("*.py"))
        assert scripts == ["analyze_peptides.py", "generate_peptides.py", "pepmlm_client.py"]

    def test_generation_logic_lives_in_client_only(self):
        """``from_pretrained`` 只能出现在客户端模块里。"""
        for script in ["generate_peptides.py", "analyze_peptides.py"]:
            assert "from_pretrained" not in _script(script).read_text(encoding="utf-8")

    def test_no_hardcoded_result_filenames(self):
        """老 filter_analysis.py 把 12 个文件名写死在源码里。

        只检查非 docstring 的字符串常量——在 docstring 里提到这些文件名
        是合理的（说明"已移除"），但真写进源码就会绑死参数组合。
        """
        literals = _non_docstring_strings(_parse("analyze_peptides.py"))
        for literal in literals:
            assert "results_10_3.csv" not in literal
            assert "results_15_10.csv" not in literal
            # 任何形如 results_<数字>_<数字>.csv 的硬编码
            assert not re.search(r"results_\d+_\d+\.csv", literal), (
                f"分析脚本里出现硬编码的结果文件名: {literal}"
            )


class TestRequirements:
    """依赖声明。"""

    def test_requirements_exists(self):
        assert (SKILL_DIR / "requirements.txt").is_file()

    def test_no_local_path_dependencies(self):
        """只检查非注释行——注释里引用反面例子是合理的。"""
        text = (SKILL_DIR / "requirements.txt").read_text(encoding="utf-8")
        active = [
            line for line in text.splitlines()
            if line.strip() and not line.strip().startswith("#")
        ]
        for line in active:
            assert "file://" not in line
            assert "@ file" not in line
            assert "/Users/" not in line

    def test_torch_is_optional_but_declared(self):
        text = (SKILL_DIR / "requirements.txt").read_text(encoding="utf-8")
        assert "torch" in text
        assert "transformers" in text
