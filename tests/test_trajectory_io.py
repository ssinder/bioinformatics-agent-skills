"""分子动力学技能的轨迹解析回归测试。

重点覆盖 DCD 解析器——它是纯 Python 重写的高风险组件，
必须与 MDAnalysis 的权威实现逐字节一致。
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "skills" / "molecular-dynamics" / "scripts"))

mda = pytest.importorskip("MDAnalysis", reason="需要 MDAnalysis 作为参考实现")

from MDAnalysis.coordinates.DCD import DCDWriter  # noqa: E402


def _write_reference_dcd(path: Path, n_atoms: int, n_frames: int) -> np.ndarray:
    """用 MDAnalysis 写一个标准 DCD，并返回真值坐标。"""
    universe = mda.Universe.empty(
        n_atoms,
        n_residues=n_atoms,
        atom_resindex=np.arange(n_atoms),
        residue_segindex=np.zeros(n_atoms, dtype=int),
        trajectory=True,
    )
    truth = []
    with DCDWriter(str(path), n_atoms) as writer:
        for frame in range(n_frames):
            coords = np.arange(n_atoms * 3, dtype=float).reshape(n_atoms, 3) + frame * 10.0
            truth.append(coords)
            universe.atoms.positions = coords
            writer.write(universe)
    return np.array(truth)


class TestDcdParser:
    """DCD 解析器与 MDAnalysis 的一致性测试。"""

    @pytest.mark.parametrize("n_atoms", [7, 12, 30, 61])
    def test_matches_mdanalysis(self, tmp_path, n_atoms):
        """解析结果必须与 MDAnalysis 完全一致。"""
        from trajectory_io import read_dcd

        path = tmp_path / "ref.dcd"
        truth = _write_reference_dcd(path, n_atoms, 4)

        coords, info = read_dcd(path)

        assert coords.shape == truth.shape
        assert np.allclose(coords, truth, atol=1e-3), "坐标与 MDAnalysis 不一致"

    def test_header_fields(self, tmp_path):
        """头信息字段应正确解析。"""
        from trajectory_io import dcd_header_info

        path = tmp_path / "ref.dcd"
        _write_reference_dcd(path, 12, 5)

        info = dcd_header_info(path)

        assert info["n_frames"] == 5
        assert info["n_atoms"] == 12
        assert info["nsavc"] == 1

    def test_rejects_non_dcd(self, tmp_path):
        """非 DCD 文件应抛出清晰错误而非静默失败。"""
        from trajectory_io import read_dcd

        bad = tmp_path / "bad.dcd"
        bad.write_bytes(b"NOTADCD" * 20)

        with pytest.raises(ValueError):
            read_dcd(bad)

    def test_read_frames_limit(self, tmp_path):
        """read_frames 参数应限制读取帧数。"""
        from trajectory_io import read_dcd

        path = tmp_path / "ref.dcd"
        truth = _write_reference_dcd(path, 12, 6)

        coords, _ = read_dcd(path, read_frames=3)

        assert coords.shape[0] == 3
        assert np.allclose(coords, truth[:3], atol=1e-3)


class TestFallbackPaths:
    """纯 Python 回退路径（无 MDAnalysis 时）应可用。"""

    def test_pdb_trajectory_reader(self, tmp_path):
        """多模型 PDB 轨迹应能被解析。"""
        from trajectory_io import read_pdb_trajectory

        pdb = tmp_path / "multi.pdb"
        pdb.write_text(
            "MODEL        1\n"
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n"
            "ATOM      2  CA  ALA A   2       1.000   0.000   0.000  1.00  0.00           C\n"
            "ENDMDL\n"
            "MODEL        2\n"
            "ATOM      1  CA  ALA A   1       0.000   0.000   1.000  1.00  0.00           C\n"
            "ATOM      2  CA  ALA A   2       1.000   0.000   1.000  1.00  0.00           C\n"
            "ENDMDL\n"
            "END\n",
            encoding="utf-8",
        )

        traj = read_pdb_trajectory(pdb)

        assert traj.coords.shape == (2, 2, 3)
        assert np.allclose(traj.coords[1, 0], [0.0, 0.0, 1.0])

    def test_xyz_trajectory_reader(self, tmp_path):
        """XYZ 轨迹应能被解析。"""
        from trajectory_io import read_xyz_trajectory

        xyz = tmp_path / "traj.xyz"
        xyz.write_text(
            "2\nframe 1\n"
            "C 0.0 0.0 0.0\n"
            "C 1.0 0.0 0.0\n"
            "2\nframe 2\n"
            "C 0.0 0.0 1.0\n"
            "C 1.0 0.0 1.0\n",
            encoding="utf-8",
        )

        traj = read_xyz_trajectory(xyz)

        assert traj.coords.shape == (2, 2, 3)
        assert np.allclose(traj.coords[1, 0], [0.0, 0.0, 1.0])
