from __future__ import annotations

import pytest
from sympy import Matrix, symbols

from magnetickp import BandPlotError, bandplot
from magnetickp.transforms import kx, ky


np = pytest.importorskip("numpy")
matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")


def test_bandplot_is_one_call_for_a_compact_path():
    mass = symbols("m")
    hamiltonian = Matrix([[mass, kx], [kx, -mass]])
    path = [
        ((0, 0, 0), (1, 0, 0), ("G", "X")),
        ((1, 0, 0), (1, 1, 0), ("X", "M")),
    ]

    figure, axis = bandplot(
        hamiltonian,
        path,
        {mass: 2},
        subdivisions=4,
        coordinate_basis="cartesian",
        show=False,
    )

    band_lines = axis.lines[:4]
    assert len(band_lines) == 4
    assert np.allclose(band_lines[0].get_ydata()[0], -2)
    assert np.allclose(band_lines[1].get_ydata()[0], 2)
    assert [tick.get_text() for tick in axis.get_xticklabels()] == ["G", "X", "M"]
    assert tuple(axis.get_xticks()) == (0, 4, 8)
    assert tuple(axis.get_xlim()) == (0, 8)
    matplotlib.pyplot.close(figure)


def test_bandplot_accepts_generated_path_records_and_string_parameters():
    mass = symbols("m")
    hamiltonian = Matrix([[mass, kx - 1j * ky], [kx + 1j * ky, -mass]])
    path = [
        {
            "start": {"coordinates": (0, 0, 0), "label": "K"},
            "end": {"coordinates": (1, 0, 0), "label": "L1"},
        }
    ]

    figure, axis = bandplot(
        hamiltonian,
        path,
        {"m": 1},
        subdivisions=2,
        coordinate_basis="cartesian",
        show=False,
    )

    assert [tick.get_text() for tick in axis.get_xticklabels()] == ["K", "L1"]
    assert np.allclose(axis.lines[0].get_ydata(), [-1, -np.sqrt(1.25), -np.sqrt(2)])
    matplotlib.pyplot.close(figure)


def test_bandplot_rejects_missing_parameters_and_nonhermitian_input():
    mass = symbols("m")
    path = [((0, 0, 0), (1, 0, 0), ("G", "X"))]

    with pytest.raises(BandPlotError, match="parameter values are required"):
        bandplot(Matrix([[mass]]), path, show=False)
    with pytest.raises(BandPlotError, match="not Hermitian"):
        bandplot(
            Matrix([[0, 1], [0, 0]]),
            path,
            coordinate_basis="cartesian",
            show=False,
        )
