"""生成 LightDock 测试夹具。

GSO 文件内容取自 LightDock 官方教程（lightdock.org/tutorials/2UUY）的
真实输出格式：

    (坐标) RecID LigID Luciferin Neighbors VisionRange Scoring

这里生成的夹具刻意包含以下**会让老脚本出错**的情形：

- ``gso_0.out`` / ``gso_10.out`` / ``gso_20.out`` / ``gso_100.out``
  并存——字典序下 ``gso_100.out`` 排在 ``gso_20.out`` 前面，
  老脚本的 ``sorted(glob(...))[-1]`` 会取到 100 步之前的文件。
- 不同 swarm 的最优分数不同，且最优 glowworm **不是 0 号**——
  老脚本把 LigID（恒为 0）当成 glowworm 编号，永远挑中 0 号。
- 最优 glowworm 编号为 7，用于验证"编号 7 的构象文件"能被正确选中。

运行：``python tests/make_lightdock_fixtures.py``
"""

from __future__ import annotations

from pathlib import Path

FIXTURES = Path(__file__).resolve().parent / "fixtures_lightdock"

#: 真实的 GSO 行模板：(坐标) RecID LigID Luciferin Neighbors VisionRange Scoring
GSO_TEMPLATE = (
    "({tx:.7f}, {ty:.7f}, {tz:.7f}, {qw:.7f}, {qx:.7f}, {qy:.7f}, {qz:.7f}) "
    "0 0 {luciferin:.8f} {neighbors} {vision:.3f} {scoring:.8f}"
)

#: 每个 swarm 的最优分数（越负越好）与对应 glowworm 编号。
#:
#: swarm_0 最优在 glowworm 7（-58.3），swarm_1 最优在 glowworm 2（-52.7），
#: swarm_2 最优在 glowworm 4（-61.1）。全局最优是 swarm_2 的 glowworm 4。
SWARM_BEST = {
    0: (7, -58.30000000),
    1: (2, -52.70000000),
    2: (4, -61.10000000),
}


def _swarm_lines(swarm: int, count: int, step: int) -> list[str]:
    """生成一个 swarm 在某一步的 GSO 行。

    打分随行号整体下降（模拟收敛过程），但最优行固定在指定编号，
    以便测试能断言"挑中的 glowworm 编号"是否正确。
    """
    best_index, best_score = SWARM_BEST[swarm]
    lines = ["#Coordinates RecID LigID Luciferin Neighbor's number Vision Range Scoring"]
    for index in range(count):
        if index == best_index:
            scoring = best_score
        else:
            # 其余行分数都比最优行差（数值更大）
            scoring = best_score + 3.0 + index * 0.37
        lines.append(
            GSO_TEMPLATE.format(
                tx=10.0 + index * 1.7 + step * 0.11,
                ty=-4.0 + index * 0.9 - step * 0.07,
                tz=-12.0 + index * 2.3 + step * 0.13,
                qw=0.7071 - index * 0.011,
                qx=0.0 + index * 0.003,
                qy=0.7071 - index * 0.009,
                qz=0.0 + index * 0.005,
                luciferin=abs(scoring) * 0.19,
                neighbors=index % 13,
                vision=4.2 + index * 0.31,
                scoring=scoring,
            )
        )
    return lines


def build_job(root: Path, job_name: str, swarms: int = 3, glowworms: int = 10) -> Path:
    """生成一个完整的模拟作业目录。

    结构::

        {job_name}/
            receptor.pdb
            ligand.pdb
            lightdock_receptor.pdb
            lightdock_ligand.pdb
            setup.json
            swarm_0/gso_0.out  gso_10.out  gso_20.out  gso_100.out
            swarm_1/gso_0.out  gso_20.out
            swarm_2/gso_0.out  gso_100.out

    其中 ``swarm_0`` 同时有 20 步与 100 步的结果，用于验证"按步数数值
    取最新"而不是"按字典序取最后"。

    Args:
        root: 父目录。
        job_name: 作业目录名。
        swarms: swarm 数量。
        glowworms: 每个 swarm 的 glowworm 数量。

    Returns:
        生成的作业目录路径。
    """
    job_dir = root / job_name
    job_dir.mkdir(parents=True, exist_ok=True)

    (job_dir / "receptor.pdb").write_text(
        "HEADER    PROTEIN                           01-JAN-00   0DOCK\n"
        "ATOM      1  N   MET A   1      10.000  20.000  30.000  1.00 20.00           N\n"
        "ATOM      2  CA  MET A   1      11.000  20.000  30.000  1.00 20.00           C\n"
        "END\n",
        encoding="utf-8",
    )
    (job_dir / "ligand.pdb").write_text(
        "HEADER    PROTEIN                           01-JAN-00   0DOCK\n"
        "ATOM      1  N   GLY B   1      50.000  50.000  50.000  1.00 20.00           N\n"
        "ATOM      2  CA  GLY B   1      51.000  50.000  50.000  1.00 20.00           C\n"
        "END\n",
        encoding="utf-8",
    )
    # setup 阶段平移后的结构
    (job_dir / "lightdock_receptor.pdb").write_text(
        "HEADER    PROTEIN                           01-JAN-00   0DOCK\n"
        "ATOM      1  N   MET A   1       0.000   0.000   0.000  1.00 20.00           N\n"
        "ATOM      2  CA  MET A   1       1.000   0.000   0.000  1.00 20.00           C\n"
        "END\n",
        encoding="utf-8",
    )
    (job_dir / "lightdock_ligand.pdb").write_text(
        "HEADER    PROTEIN                           01-JAN-00   0DOCK\n"
        "ATOM      1  N   GLY B   1       0.000   0.000   0.000  1.00 20.00           N\n"
        "ATOM      2  CA  GLY B   1       1.000   0.000   0.000  1.00 20.00           C\n"
        "END\n",
        encoding="utf-8",
    )
    (job_dir / "setup.json").write_text(
        '{\n'
        '    "swarms": 3,\n'
        '    "glowworms": 10,\n'
        '    "receptor_pdb": "receptor.pdb",\n'
        '    "ligand_pdb": "ligand.pdb",\n'
        '    "use_anm": false,\n'
        '    "scoring_function": "dfire"\n'
        '}\n',
        encoding="utf-8",
    )

    # 每个 swarm 存在的步数文件——刻意让字典序与数值序不一致
    steps_per_swarm = {
        0: (0, 10, 20, 100),
        1: (0, 20),
        2: (0, 100),
    }
    for swarm in range(swarms):
        swarm_dir = job_dir / f"swarm_{swarm}"
        swarm_dir.mkdir(parents=True, exist_ok=True)
        for step in steps_per_swarm.get(swarm, (0,)):
            lines = _swarm_lines(swarm, glowworms, step)
            (swarm_dir / f"gso_{step}.out").write_text(
                "\n".join(lines) + "\n", encoding="utf-8"
            )

    return job_dir


def main() -> None:
    """生成全部夹具。"""
    FIXTURES.mkdir(parents=True, exist_ok=True)

    # 单个作业（用于主测试）
    build_job(FIXTURES, "job_basic")

    # 多个作业（用于批量测试）
    for name in ("job_a", "job_b"):
        build_job(FIXTURES, name)

    print(f"夹具已生成到: {FIXTURES}")


if __name__ == "__main__":
    main()
