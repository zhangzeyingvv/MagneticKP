from __future__ import annotations

import pytest
from sympy import Matrix, Rational, eye

from magnetickp import (
    HighSymmetryPathError,
    format_band_path,
    high_symmetry_k_path,
)
from magnetickp.transforms import kx, ky, kz


def test_c2_path_uses_the_largest_line_stabilizer_and_input_metadata():
    inputdict = {
        "Unitary": {"C2z": (eye(1), (-kx, -ky, kz))},
        "KPointLabel": "K",
        "CoordinateSystem": "ReciprocalFractional",
    }

    data = high_symmetry_k_path(inputdict, return_data=True)

    assert data["directions"] == ((0, 0, 1),)
    assert data["little_cogroup_order"] == 2
    assert data["enumerated_subgroup_count"] == 2
    assert data["center_label"] == "K"
    assert data["coordinate_basis"] == "fractional"
    assert data["path"] == [
        {
            "start": {
                "coordinates": (0, 0, -Rational(1, 10)),
                "label": "L1-",
            },
            "end": {"coordinates": (0, 0, 0), "label": "K"},
        },
        {
            "start": {"coordinates": (0, 0, 0), "label": "K"},
            "end": {
                "coordinates": (0, 0, Rational(1, 10)),
                "label": "L1+",
            },
        },
    ]


def test_d4h_returns_three_representatives_and_five_fixed_lines():
    inputdict = {
        "Unitary": {
            "C4z": (eye(1), (-ky, kx, kz)),
            "mxz": (eye(1), (kx, -ky, kz)),
            "mxy": (eye(1), (kx, ky, -kz)),
        }
    }

    representatives = high_symmetry_k_path(inputdict, return_data=True)
    all_lines = high_symmetry_k_path(
        inputdict,
        equivalent_directions="all",
        return_data=True,
    )

    assert representatives["little_cogroup_order"] == 16
    assert len(representatives["directions"]) == 3
    assert [line["subgroup_order"] for line in representatives["line_data"]] == [
        8,
        4,
        4,
    ]
    assert len(all_lines["directions"]) == 5
    assert len(high_symmetry_k_path(inputdict, include_negative_directions=False)) == 3


def test_msgcorep_fractional_example_matches_mathematica_directions():
    pytest.importorskip("msgcorep")
    from magnetickp.catalog import build_catalog_core_input

    core_input, _metadata, _group_input = build_catalog_core_input(
        space_group=156,
        magnetic_group="156.50",
        k_name="K",
        representation_family="double",
        representations=[{"index": 2}, {"index": 3}],
        coordinate_basis="fractional",
    )

    data = high_symmetry_k_path(core_input, return_data=True)

    assert data["directions"] == ((0, 0, 1), (1, 1, 0))
    assert data["little_cogroup_order"] == 6
    assert data["center_label"] == "K"
    assert format_band_path(data["path"]) == (
        "L1-:0,0,-1/10 -> K:0,0,0; K:0,0,0 -> L1+:0,0,1/10; "
        "L2-:-1/10,-1/10,0 -> K:0,0,0; K:0,0,0 -> L2+:1/10,1/10,0"
    )


def test_invalid_path_options_fail_explicitly():
    with pytest.raises(HighSymmetryPathError):
        high_symmetry_k_path({})
    with pytest.raises(HighSymmetryPathError):
        high_symmetry_k_path(
            {"Unitary": {"identity": (Matrix([[1]]), eye(3))}},
            radius=0,
        )
