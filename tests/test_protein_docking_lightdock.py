"""LightDock 对接技能的回归测试。

锁住老仓库 ``dock_lightdock/scripts/`` 的 5 个脚本带来的真实问题，
以及 LightDock 输出语义里最容易踩的坑：

1. **把 LigID 当成 glowworm 索引** —— 老正则的捕获组 1 注释写着
   "LigID"，却被当成 glowworm 编号使用。而单体对接时 LigID **恒为 0**。
   于是"最佳构象"永远是 glowworm #0，与真实得分最高的 glowworm 毫无关系。
   正确语义：**glowworm 编号 = 该行在 GSO 文件中的行号**（从 0 起）。
2. **GSO 文件按字典序取"最后一个"** —— LightDock 每 10 步存一次，
   文件名不补零，字典序下 ``gso_100.out`` 排在 ``gso_20.out`` 前面。
   老代码 ``sorted(glob(...))[-1]`` 会用上**未收敛的中间结果**。
3. **生成构象用原始 PDB** —— ``lightdock3_setup.py`` 会把两个结构平移到
   坐标原点并另存为 ``lightdock_*.pdb``，GSO 里记录的平移量正是相对这个
   原点。老代码传的是原始文件，复合物两个部分对不上。
4. **构象文件按字典序排序** —— ``lgd_generate_conformations.py`` 输出
   ``lightdock_<i>.pdb``，i 从 0 起不补零，字典序下 ``lightdock_10.pdb``
   排在 ``lightdock_2.pdb`` 前面。
5. **硬编码个人路径** —— ``/Users/jung/Desktop/trae/tool/lightdock/dockallss``
   与 ``/Users/jung/sinder/conda/anaconda3/envs/lightdock/bin`` 遍布全部脚本。
6. **三份逐字重复的脚本** —— ``dock_completess.py`` / ``dock_simple.py`` /
   ``dock_test_simple.py`` 各自复制了一份 run_docking / extract_best_score /
   generate_scores_csv。
7. **配体链 ID 原地改写 + 备份还原** —— 靠 ``finally`` 兜底，中途崩溃
   就留下被改坏的配体文件，后续所有对接都用错误的输入。
8. **每 500 个作业才写一次 CSV** —— 中途 Ctrl-C 丢掉前面全部结果。
9. **单作业失败中止整个批次** —— 一次失败毁掉全部。
10. **笛卡尔积没有预演** —— 直接开几千个作业。

全部用例不联网，也不真的跑 LightDock。
"""

from __future__ import annotations

import ast
import csv
import json
import shutil
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
SKILL_DIR = REPO_ROOT / "skills" / "protein-docking-lightdock"
FIXTURES = REPO_ROOT / "tests" / "fixtures_lightdock"
sys.path.insert(0, str(REPO_ROOT / "packages"))
sys.path.insert(0, str(SKILL_DIR / "scripts"))

import extract_scores  # noqa: E402
import generate_conformations  # noqa: E402
import lightdock_client  # noqa: E402
import run_docking  # noqa: E402

SCRIPTS = ("lightdock_client.py", "run_docking.py", "generate_conformations.py", "extract_scores.py")

#: 老仓库里的个人路径，任何一个出现在新技能里都是回归。
FORBIDDEN_PATHS = (
    "/Users/jung",
    "/Users/DELL",
    "C:\\Users\\DELL",
    "trae/tool/lightdock",
    "anaconda3/envs/lightdock",
)

#: 老仓库里被写死的作业目录名。
FORBIDDEN_LITERALS = ("dockallss", "completess", "ar781.pdb", "extracted_proteins")


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


def _all_strings() -> list[str]:
    """返回全部脚本的所有字符串字面量（含 docstring）。"""
    result: list[str] = []
    for name in SCRIPTS:
        result.extend(_non_docstring_strings(_parse(name)))
    return result


def _top_level_calls(name: str) -> list[str]:
    """返回**模块顶层**直接执行的调用名。

    只看 ``tree.body`` 的顶层语句与函数默认值——下钻到函数体会把
    正常业务逻辑里的 ``open()`` 误报成"模块级副作用"。
    """
    tree = _parse(name)
    calls: list[str] = []
    for node in tree.body:
        if isinstance(node, ast.ClassDef):
            # 类体里的方法体是运行时才执行的，不算模块级副作用；
            # 但方法**默认值**里的调用是在类定义时求值的，要看。
            for sub in node.body:
                if isinstance(sub, (ast.FunctionDef, ast.AsyncFunctionDef)):
                    for default in list(sub.args.defaults) + [
                        d for d in sub.args.kw_defaults if d is not None
                    ]:
                        for inner in ast.walk(default):
                            if isinstance(inner, ast.Call):
                                calls.append(getattr(inner.func, "id", "?"))
            continue
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            for default in list(node.args.defaults) + [
                d for d in node.args.kw_defaults if d is not None
            ]:
                for sub in ast.walk(default):
                    if isinstance(sub, ast.Call):
                        calls.append(getattr(sub.func, "id", "?"))
            continue
        if isinstance(node, ast.Expr) and isinstance(node.value, ast.Constant):
            continue  # docstring
        for sub in ast.walk(node):
            if isinstance(sub, ast.Call):
                func = sub.func
                calls.append(
                    getattr(func, "id", None) or getattr(func, "attr", None) or "?"
                )
    return calls


class TestNoHardcodedPersonalPaths:
    """老仓库把作者的个人路径写死在所有脚本里。"""

    @pytest.mark.parametrize("name", SCRIPTS)
    def test_no_personal_paths(self, name: str) -> None:
        content = _script(name).read_text(encoding="utf-8")
        for bad in FORBIDDEN_PATHS:
            assert bad not in content, f"{name} 残留硬编码路径: {bad}"

    @pytest.mark.parametrize("name", SCRIPTS)
    def test_no_old_project_literals(self, name: str) -> None:
        strings = _non_docstring_strings(_parse(name))
        for literal in FORBIDDEN_LITERALS:
            assert not any(literal in s for s in strings), (
                f"{name} 残留老仓库专用字面量: {literal}"
            )

    @pytest.mark.parametrize("name", SCRIPTS)
    def test_no_module_level_side_effects(self, name: str) -> None:
        """模块导入时不应有写文件/建目录等副作用。"""
        forbidden = {"open", "mkdir", "makedirs", "rmtree", "copytree", "system"}
        for call in _top_level_calls(name):
            assert call not in forbidden, f"{name} 顶层调用了 {call}()"


class TestGsoParsing:
    """GSO 文件的字段语义——最容易读错的地方。"""

    def test_parses_all_fields(self) -> None:
        line = (
            "(31.4171143, 1.8570079, -6.3956223, -0.1058407, -0.4849369, "
            "0.5997430, -0.6276482) 0 0 11.25395618 0 4.200 7.52800101"
        )
        entry = lightdock_client.parse_gso_line(line, glowworm=3)
        assert entry is not None
        assert entry.rec_id == 0
        assert entry.lig_id == 0
        assert entry.luciferin == pytest.approx(11.25395618)
        assert entry.neighbors == 0
        assert entry.vision_range == pytest.approx(4.200)
        assert entry.scoring == pytest.approx(7.52800101)
        assert entry.glowworm == 3

    def test_glowworm_id_comes_from_line_number_not_ligid(self) -> None:
        """核心回归：glowworm 编号必须来自行号，而不是 LigID 列。

        老脚本的正则捕获 LigID 当作 glowworm 索引。单体对接时 LigID
        恒为 0，所以老脚本挑出来的"最佳构象"永远是 0 号。
        """
        line = "(1.0, 2.0, 3.0, 0.5, 0.5, 0.5, 0.5) 0 0 12.5 3 4.2 -41.5"
        for line_index in range(10):
            entry = lightdock_client.parse_gso_line(line, glowworm=line_index)
            assert entry is not None
            # 无论 LigID 是什么，glowworm 编号必须等于传入的行号
            assert entry.glowworm == line_index
            # LigID 仍然是 0（这正是老脚本踩坑的原因）
            assert entry.lig_id == 0

    def test_skips_comment_lines(self) -> None:
        assert lightdock_client.parse_gso_line("#Coordinates RecID ...", 0) is None
        assert lightdock_client.parse_gso_line("", 0) is None
        assert lightdock_client.parse_gso_line("   \n", 0) is None

    def test_rejects_malformed_line(self) -> None:
        assert lightdock_client.parse_gso_line("total garbage", 0) is None

    def test_keeps_negative_scores(self) -> None:
        """DFIRE 打分为负，必须原样保留，不能被当成无效行丢弃。"""
        line = "(1.0, 2.0, 3.0, 0.5, 0.5, 0.5, 0.5) 0 0 9.5 2 4.2 -58.3"
        entry = lightdock_client.parse_gso_line(line, glowworm=0)
        assert entry is not None
        assert entry.scoring == pytest.approx(-58.3)

    def test_supports_anm_extended_coordinates(self) -> None:
        """开启 ANM 时坐标括号里会有更多数字，解析仍需正确。"""
        line = (
            "(1.0, 2.0, 3.0, 0.5, 0.5, 0.5, 0.5, 0.1, -0.2, 0.05) "
            "0 0 27.83 6 0.32 18.66"
        )
        entry = lightdock_client.parse_gso_line(line, glowworm=0)
        assert entry is not None
        assert entry.scoring == pytest.approx(18.66)
        assert entry.neighbors == 6

    def test_parse_gso_file(self) -> None:
        entries = lightdock_client.parse_gso_file(FIXTURES / "gso_sample.out")
        assert len(entries) == 5
        assert [e.glowworm for e in entries] == [0, 1, 2, 3, 4]
        assert entries[0].scoring == pytest.approx(7.52800101)

    def test_parse_gso_file_missing(self) -> None:
        with pytest.raises(lightdock_client.LightDockError):
            lightdock_client.parse_gso_file(FIXTURES / "does_not_exist.out")


class TestLatestGsoSelection:
    """GSO 文件必须按步数**数值**挑，而不是字典序。"""

    def test_picks_highest_numeric_step(self, tmp_path: Path) -> None:
        swarm = tmp_path / "swarm_0"
        swarm.mkdir()
        for step in (0, 10, 20, 100):
            (swarm / f"gso_{step}.out").write_text("x\n", encoding="utf-8")
        latest = lightdock_client.find_latest_gso(swarm)
        assert latest is not None
        assert latest.name == "gso_100.out"

    def test_dictionary_order_would_pick_wrong_file(self, tmp_path: Path) -> None:
        """夹具必须能区分"字典序最后"与"数值最大"。"""
        swarm = tmp_path / "swarm_0"
        swarm.mkdir()
        for step in (0, 10, 20, 100):
            (swarm / f"gso_{step}.out").write_text("x\n", encoding="utf-8")
        dictionary_last = sorted(swarm.glob("gso_*.out"))[-1]
        assert dictionary_last.name == "gso_20.out"  # 字典序的最后是 20
        latest = lightdock_client.find_latest_gso(swarm)
        assert latest is not None
        assert latest.name == "gso_100.out"  # 数值最大的是 100

    def test_ignores_non_gso_files(self, tmp_path: Path) -> None:
        swarm = tmp_path / "swarm_0"
        swarm.mkdir()
        (swarm / "gso_10.out").write_text("x\n", encoding="utf-8")
        (swarm / "gso_backup.out").write_text("x\n", encoding="utf-8")
        (swarm / "setup.json").write_text("{}", encoding="utf-8")
        latest = lightdock_client.find_latest_gso(swarm)
        assert latest is not None
        assert latest.name == "gso_10.out"

    def test_returns_none_without_gso(self, tmp_path: Path) -> None:
        assert lightdock_client.find_latest_gso(tmp_path) is None


class TestSwarmResults:
    """按 swarm 汇总与挑全局最优。"""

    @pytest.fixture()
    def job_dir(self) -> Path:
        return FIXTURES / "job_basic"

    def test_collects_all_swarms(self, job_dir: Path) -> None:
        results = lightdock_client.collect_swarm_results(job_dir)
        assert len(results) == 3
        assert [r.swarm for r in results] == [0, 1, 2]

    def test_picks_correct_best_glowworm(self, job_dir: Path) -> None:
        """swarm_0 的最优在 glowworm 7，不是 0 号。"""
        results = {r.swarm: r for r in lightdock_client.collect_swarm_results(job_dir)}
        best = results[0].best
        assert best is not None
        assert best.glowworm == 7
        assert best.scoring == pytest.approx(-58.3)

    def test_uses_highest_step_file(self, job_dir: Path) -> None:
        """swarm_0 同时有 20 步与 100 步，必须用 100 步那份。"""
        results = {r.swarm: r for r in lightdock_client.collect_swarm_results(job_dir)}
        assert results[0].gso_file.name == "gso_100.out"

    def test_global_best_is_lowest_score(self, job_dir: Path) -> None:
        """全局最优是所有 swarm 里分数最低（最负）的那个。"""
        best = lightdock_client.find_global_best(job_dir)
        assert best is not None
        swarm, entry = best
        assert swarm.swarm == 2
        assert entry.glowworm == 4
        assert entry.scoring == pytest.approx(-61.1)

    def test_lower_score_is_better(self, job_dir: Path) -> None:
        """确认方向：分数越负越好，不是越大越好。"""
        best = lightdock_client.find_global_best(job_dir)
        assert best is not None
        assert best[1].scoring < best[1].luciferin

    def test_empty_directory(self, tmp_path: Path) -> None:
        assert lightdock_client.collect_swarm_results(tmp_path) == []
        assert lightdock_client.find_global_best(tmp_path) is None


class TestConformationSorting:
    """构象文件必须按 glowworm 编号**数值**排序。"""

    def test_numeric_sort(self, tmp_path: Path) -> None:
        for index in (0, 1, 2, 10, 11):
            (tmp_path / f"lightdock_{index}.pdb").write_text("END\n", encoding="utf-8")
        files = lightdock_client.sort_conformation_files(tmp_path)
        names = [p.name for p in files]
        assert names == [
            "lightdock_0.pdb",
            "lightdock_1.pdb",
            "lightdock_2.pdb",
            "lightdock_10.pdb",
            "lightdock_11.pdb",
        ]

    def test_dictionary_order_would_be_wrong(self, tmp_path: Path) -> None:
        for index in (0, 1, 2, 10, 11):
            (tmp_path / f"lightdock_{index}.pdb").write_text("END\n", encoding="utf-8")
        dictionary = sorted(tmp_path.glob("lightdock_*.pdb"))
        assert dictionary[-1].name == "lightdock_2.pdb"  # 字典序的最后是 2
        numeric = lightdock_client.sort_conformation_files(tmp_path)
        assert numeric[-1].name == "lightdock_11.pdb"

    def test_moved_structure_path(self) -> None:
        moved = lightdock_client.moved_structure_path(Path("/tmp/job/receptor.pdb"))
        assert moved.name == "lightdock_receptor.pdb"
        assert moved.parent == Path("/tmp/job")


class TestExecutableLocation:
    """LightDock 的脚本不在 PATH 里，必须多找几处。"""

    def test_finds_in_sys_prefix(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        fake_prefix = tmp_path / "env"
        bin_dir = fake_prefix / "bin"
        bin_dir.mkdir(parents=True)
        script = bin_dir / lightdock_client.SETUP_SCRIPT
        script.write_text("#!/usr/bin/env python3\n", encoding="utf-8")

        monkeypatch.setattr(lightdock_client.sys, "prefix", str(fake_prefix))
        found = lightdock_client.locate_executable(lightdock_client.SETUP_SCRIPT)
        assert found == script

    def test_explicit_path_wins(self, tmp_path: Path) -> None:
        script = tmp_path / "custom_setup.py"
        script.write_text("#!/usr/bin/env python3\n", encoding="utf-8")
        found = lightdock_client.locate_executable(
            lightdock_client.SETUP_SCRIPT, str(script)
        )
        assert found == script

    def test_explicit_directory(self, tmp_path: Path) -> None:
        script = tmp_path / lightdock_client.RUN_SCRIPT
        script.write_text("#!/usr/bin/env python3\n", encoding="utf-8")
        found = lightdock_client.locate_executable(lightdock_client.RUN_SCRIPT, str(tmp_path))
        assert found == script

    def test_missing_returns_none(self, tmp_path: Path) -> None:
        found = lightdock_client.locate_executable(
            "definitely_not_a_real_script_xyz.py", str(tmp_path)
        )
        assert found is None

    def test_resolve_requires_setup_and_run(self) -> None:
        with pytest.raises(lightdock_client.LightDockError) as excinfo:
            lightdock_client.resolve_executables(
                setup="/nonexistent/setup.py", run="/nonexistent/run.py"
            )
        assert "lightdock3_setup.py" in str(excinfo.value)


class TestChainIdRewrite:
    """配体链 ID 必须改副本，不能动原文件。"""

    def test_creates_copy_and_leaves_original(self, tmp_path: Path) -> None:
        source = tmp_path / "lightdock_ligand.pdb"
        source.write_text(
            "ATOM      1  N   GLY B   1      50.000  50.000  50.000  1.00 20.00           N\n"
            "ATOM      2  CA  GLY B   1      51.000  50.000  50.000  1.00 20.00           C\n"
            "END\n",
            encoding="utf-8",
        )
        original = source.read_text(encoding="utf-8")

        copy = generate_conformations.rewrite_chain_id(source, "C")
        assert copy != source
        assert copy.is_file()
        # 原文件完全没动
        assert source.read_text(encoding="utf-8") == original
        # 副本的链 ID 改了
        assert "GLY C" in copy.read_text(encoding="utf-8")

    def test_preserves_non_atom_records(self, tmp_path: Path) -> None:
        source = tmp_path / "lig.pdb"
        source.write_text(
            "HEADER    PROTEIN\nATOM      1  N   GLY B   1      50.000  50.000  50.000  1.00 20.00           N\n"
            "END\n",
            encoding="utf-8",
        )
        copy = generate_conformations.rewrite_chain_id(source, "A")
        lines = copy.read_text(encoding="utf-8").splitlines()
        assert lines[0] == "HEADER    PROTEIN"
        assert lines[-1] == "END"
        assert "GLY A" in lines[1]


class TestExtractScores:
    """分数汇总。"""

    def test_global_summary(self, tmp_path: Path) -> None:
        job_root = tmp_path / "results"
        shutil.copytree(FIXTURES / "job_basic", job_root / "job_basic")
        output_dir = tmp_path / "scores"

        code = extract_scores.main(
            ["--jobs-root", str(job_root), "--output-dir", str(output_dir)]
        )
        assert code == 0

        csv_path = output_dir / "docking_scores.csv"
        assert csv_path.is_file()
        with csv_path.open(encoding="utf-8-sig") as handle:
            rows = list(csv.DictReader(handle))
        assert len(rows) == 1
        assert rows[0]["job_name"] == "job_basic"
        assert rows[0]["swarm"] == "2"
        assert rows[0]["glowworm"] == "4"
        assert float(rows[0]["score"]) == pytest.approx(-61.1)

    def test_per_swarm_summary(self, tmp_path: Path) -> None:
        job_root = tmp_path / "results"
        shutil.copytree(FIXTURES / "job_basic", job_root / "job_basic")
        output_dir = tmp_path / "scores"

        extract_scores.main(
            [
                "--jobs-root",
                str(job_root),
                "--output-dir",
                str(output_dir),
                "--per-swarm",
            ]
        )
        with (output_dir / "docking_scores.csv").open(encoding="utf-8-sig") as handle:
            rows = list(csv.DictReader(handle))
        assert len(rows) == 3
        # 按分数升序（最优在前）
        scores = [float(r["best_score"]) for r in rows]
        assert scores == sorted(scores)
        assert scores[0] == pytest.approx(-61.1)

    def test_json_written(self, tmp_path: Path) -> None:
        job_root = tmp_path / "results"
        shutil.copytree(FIXTURES / "job_basic", job_root / "job_basic")
        output_dir = tmp_path / "scores"
        extract_scores.main(["--jobs-root", str(job_root), "--output-dir", str(output_dir)])
        data = json.loads((output_dir / "docking_scores.json").read_text(encoding="utf-8"))
        assert len(data) == 1
        assert data[0]["glowworm"] == 4

    def test_missing_args_returns_2(self) -> None:
        assert extract_scores.main([]) == 2


class TestRunDockingCli:
    """批量对接 CLI 的参数与前置检查。"""

    def test_requires_receptor(self, tmp_path: Path) -> None:
        code = run_docking.main(
            ["--ligands", str(FIXTURES / "job_basic" / "ligand.pdb"), "--output-dir", str(tmp_path)]
        )
        assert code == 2

    def test_missing_file_returns_2(self, tmp_path: Path) -> None:
        code = run_docking.main(
            [
                "--receptor",
                str(tmp_path / "nope.pdb"),
                "--ligands",
                str(tmp_path / "nope2.pdb"),
                "--output-dir",
                str(tmp_path / "out"),
            ]
        )
        assert code == 2

    def test_missing_executables_returns_3(self, tmp_path: Path) -> None:
        code = run_docking.main(
            [
                "--receptor",
                str(FIXTURES / "job_basic" / "receptor.pdb"),
                "--ligands",
                str(FIXTURES / "job_basic" / "ligand.pdb"),
                "--output-dir",
                str(tmp_path / "out"),
                "--setup-exe",
                str(tmp_path / "no_setup.py"),
                "--run-exe",
                str(tmp_path / "no_run.py"),
            ]
        )
        assert code == 3

    def test_expand_ligands_from_directory(self, tmp_path: Path) -> None:
        (tmp_path / "a.pdb").write_text("END\n", encoding="utf-8")
        (tmp_path / "b.pdb").write_text("END\n", encoding="utf-8")
        (tmp_path / "notes.txt").write_text("x\n", encoding="utf-8")
        ligands = run_docking.expand_ligands([str(tmp_path)])
        assert [p.name for p in ligands] == ["a.pdb", "b.pdb"]

    def test_expand_ligands_keeps_explicit_paths(self, tmp_path: Path) -> None:
        path = tmp_path / "x.pdb"
        path.write_text("END\n", encoding="utf-8")
        assert run_docking.expand_ligands([str(path)]) == [path]


class TestSkillDocumentation:
    """SKILL.md 与参考文档必须把三个坑讲清楚。"""

    @pytest.fixture(scope="class")
    @classmethod
    def skill_md(cls) -> str:
        return (SKILL_DIR / "SKILL.md").read_text(encoding="utf-8")

    def test_has_workflow_section(self, skill_md: str) -> None:
        assert "## Workflow" in skill_md

    def test_documents_glowworm_line_number(self, skill_md: str) -> None:
        assert "line number" in skill_md.lower()

    def test_documents_lower_score_is_better(self, skill_md: str) -> None:
        assert "lower" in skill_md.lower() and "better" in skill_md.lower()

    def test_methodology_explains_gso_fields(self) -> None:
        text = (SKILL_DIR / "references" / "methodology.md").read_text(encoding="utf-8")
        assert "LigID" in text
        assert "行号" in text
        assert "gso_100.out" in text

    def test_methodology_explains_coordinate_frame(self) -> None:
        text = (SKILL_DIR / "references" / "methodology.md").read_text(encoding="utf-8")
        assert "move_to_origin" in text
        assert "lightdock_receptor.pdb" in text

    def test_troubleshooting_covers_common_failures(self) -> None:
        text = (SKILL_DIR / "references" / "troubleshooting.md").read_text(encoding="utf-8")
        assert "Windows" in text
        assert "already exists" in text
