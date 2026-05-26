from pathlib import Path

from vhm.io_utils import read_xyz


def test_read_xyz_allows_blank_comment_line(tmp_path: Path):
    xyz_file = tmp_path / "blank_comment.xyz"
    xyz_file.write_text(
        "2\n"
        "\n"
        "H 0.0 0.0 0.0\n"
        "H 0.0 0.0 1.0\n",
        encoding="utf-8",
    )

    symbols, coords = read_xyz(xyz_file, "angstrom", 1.0)

    assert symbols == ["H", "H"]
    assert coords.tolist() == [[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]]


def test_read_xyz_ignores_blank_lines_between_atoms(tmp_path: Path):
    xyz_file = tmp_path / "blank_atom_lines.xyz"
    xyz_file.write_text(
        "2\n"
        "comment\n"
        "H 0.0 0.0 0.0\n"
        "\n"
        "H 0.0 0.0 1.0\n",
        encoding="utf-8",
    )

    symbols, coords = read_xyz(xyz_file, "angstrom", 1.0)

    assert symbols == ["H", "H"]
    assert coords.tolist() == [[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]]
