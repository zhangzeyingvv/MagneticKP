from __future__ import annotations

import pytest
from sympy import Matrix, symbols

from magnetickp import KpFittingError, fit_kp_to_vasp, fittingKP
from magnetickp.transforms import kx


np = pytest.importorskip("numpy")
pytest.importorskip("scipy")


def write_eigenval(path):
    path.write_text(
        """    1    1    1    1
 header
 header
 header
 header
      2      3      1

  0.0  0.0  0.0  0.333333
    1   1.0  1.0

  1.0  0.0  0.0  0.333333
    1   3.0  1.0

  2.0  0.0  0.0  0.333333
    1   5.0  0.0
""",
        encoding="utf-8",
    )
    return path


def test_fittingkp_fits_a_hamiltonian_from_one_file_and_kpoint_range(tmp_path):
    offset, slope = symbols("offset slope", real=True)
    hamiltonian = Matrix([[offset + slope * kx]])

    result = fittingKP(
        hamiltonian,
        write_eigenval(tmp_path / "EIGENVAL"),
        range(1, 4),
    )

    assert result["converged"] is True
    assert result["kpoint_indices"] == (1, 2, 3)
    assert result["band_range"] == (1, 1)
    assert result["fitted_parameters"][offset] == pytest.approx(1.0)
    assert result["fitted_parameters"][slope] == pytest.approx(2.0)
    assert result["fitted_least_squares"] == pytest.approx(0.0, abs=1e-20)


def test_fit_alias_and_named_initial_parameters(tmp_path):
    offset, slope = symbols("offset slope", real=True)
    result = fit_kp_to_vasp(
        Matrix([[offset + slope * kx]]),
        write_eigenval(tmp_path / "EIGENVAL"),
        [1, 2, 3],
        {"offset": 0, "slope": 0},
    )

    assert tuple(map(str, result["parameter_order"])) == ("offset", "slope")
    assert result["initial_least_squares"] == pytest.approx(35.0)


def test_fittingkp_shifts_eigenval_coordinates_to_the_expansion_center(tmp_path):
    offset, slope = symbols("offset slope", real=True)
    result = fittingKP(
        Matrix([[offset + slope * kx]]),
        write_eigenval(tmp_path / "EIGENVAL"),
        [1, 2, 3],
        kpoint_center=(1, 0, 0),
    )

    assert result["kpoint_center"] == (1.0, 0.0, 0.0)
    assert result["fitted_parameters"][offset] == pytest.approx(3.0)
    assert result["fitted_parameters"][slope] == pytest.approx(2.0)


def test_fittingkp_rejects_invalid_kpoint_indices(tmp_path):
    parameter = symbols("p", real=True)
    path = write_eigenval(tmp_path / "EIGENVAL")

    with pytest.raises(KpFittingError, match="positive one-based"):
        fittingKP(Matrix([[parameter]]), path, [0])
    with pytest.raises(KpFittingError, match="exceeds"):
        fittingKP(Matrix([[parameter]]), path, [4])
