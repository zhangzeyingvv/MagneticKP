from __future__ import annotations

import pytest

from magnetickp import VaspEigenvalError, read_vasp_eigenval


EIGENVAL_TEXT = """    1    1    1    2
 header
 header
 header
 header
      4      2      3

  0.0  0.0  0.0  0.5
    1  -2.0  -1.5  1.0  1.0
    2  -1.0  -0.5  1.0  1.0
    3   1.0   1.5  0.0  0.0

  0.5  0.0  0.0  0.5
    1  -1.8  -1.3  1.0  1.0
    2  -0.8  -0.3  1.0  1.0
    3   1.2   1.7  0.0  0.0
"""


def eigenval_file(tmp_path):
    path = tmp_path / "EIGENVAL"
    path.write_text(EIGENVAL_TEXT, encoding="utf-8")
    return path


def test_read_vasp_eigenval_selects_one_based_bands_and_shifts_fermi_level(tmp_path):
    records = read_vasp_eigenval(
        eigenval_file(tmp_path),
        efermi=0.2,
        spin=1,
        start_band=2,
        end_band=3,
    )

    assert records == [
        ((0.0, 0.0, 0.0), [-1.2, 0.8]),
        ((0.5, 0.0, 0.0), [-1.0, 1.0]),
    ]


def test_read_vasp_eigenval_selects_the_second_spin_channel(tmp_path):
    records = read_vasp_eigenval(
        eigenval_file(tmp_path), spin=2, start_band=1, end_band=1
    )

    assert records[0] == ((0.0, 0.0, 0.0), [-1.5])
    assert records[1] == ((0.5, 0.0, 0.0), [-1.3])


def test_read_vasp_eigenval_rejects_invalid_selections(tmp_path):
    path = eigenval_file(tmp_path)

    with pytest.raises(VaspEigenvalError, match="positive one-based"):
        read_vasp_eigenval(path, spin=0)
    with pytest.raises(VaspEigenvalError, match="exceeds"):
        read_vasp_eigenval(path, start_band=3, end_band=4)
    with pytest.raises(VaspEigenvalError, match="cannot read"):
        read_vasp_eigenval(tmp_path / "missing-EIGENVAL")
