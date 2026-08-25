from __future__ import annotations

import pytest
import sympy
from sympy import Matrix

msgcorep = pytest.importorskip("msgcorep")

from magnetickp.catalog import (
    CatalogInputError,
    build_catalog_core_input,
    catalog_status,
    corepresentation_catalog,
    exact_to_sympy,
    k_point_catalog,
    list_bz_variants,
    list_magnetic_groups,
    list_space_groups,
)
from magnetickp.web.app import (
    CatalogSolveRequest,
    app,
    compute_catalog_request,
    run_catalog_isolated,
)


def _default_request(**changes):
    values = {
        "space_group": 92,
        "magnetic_group": "92.112",
        "bz_variant": "a",
        "k_name": "Γ",
        "coordinate_basis": "cartesian",
        "representation_family": "double",
        "representations": [{"index": 1, "multiplicity": 1}],
        "order": 1,
        "method": "iterative",
        "max_group_order": 512,
    }
    values.update(changes)
    return CatalogSolveRequest(**values)


def test_catalog_exposes_all_space_groups_and_msg_family_members():
    status = catalog_status()
    assert status["available"] is True
    assert status["version"] == msgcorep.__version__

    space_groups = list_space_groups()
    assert len(space_groups) == 230
    assert space_groups[91]["number"] == 92
    assert space_groups[91]["symbol"] == "P41212"
    assert space_groups[91]["symbol_tex"] == "P4_12_12"
    assert space_groups[91]["bc_symbol_tex"] == "P4_12_12"

    magnetic_groups = list_magnetic_groups(92)
    assert [row["identifier"] for row in magnetic_groups] == [
        f"92.{number}" for number in range(111, 119)
    ]
    assert {row["type"] for row in magnetic_groups} == {1, 2, 3, 4}
    assert magnetic_groups[1]["bns_symbol_tex"] == "P4_{1}2_{1}21'"
    assert magnetic_groups[5]["bns_symbol_tex"] == "P_{c}4_{1}2_{1}2"


def test_k_point_and_corepresentation_catalogs_keep_exact_data():
    assert list_bz_variants(92) == [
        {"value": "a", "label": "TetrPrim", "full_bz_type": "TetrPrim"}
    ]
    points = k_point_catalog(92, magnetic_group="92.112")["points"]
    gamma = next(point for point in points if point["name"] == "Γ")
    delta = next(point for point in points if point["name"] == "Δ")
    assert gamma["coordinate"] == ["0", "0", "0"]
    assert r"\begin{matrix}0 & 0 & 0\end{matrix}" in gamma["coordinate_latex"]
    assert gamma["coordinate_latex"].endswith(r"^{\mathsf T}")
    assert gamma["parameters"] == []
    assert delta["coordinate"] == ["0", "u", "0"]
    assert delta["parameters"] == ["u"]

    coreps = corepresentation_catalog("92.112", "Γ")
    assert coreps["group_order"] == 16
    assert [row["label"] for row in coreps["double"]] == ["Γ6", "Γ7"]
    assert [row["dimension"] for row in coreps["double"]] == [2, 2]
    assert len(coreps["generators"]) == 3
    assert [row["momentum_image"] for row in coreps["generators"]] == [
        ["-ky", "kx", "kz"],
        ["kx", "-ky", "-kz"],
        ["-kx", "-ky", "-kz"],
    ]
    assert [row["antiunitary"] for row in coreps["generators"]] == [
        False,
        False,
        True,
    ]
    gamma6 = coreps["double"][0]
    assert [row["generator_index"] for row in gamma6["generator_matrices"]] == [
        5,
        2,
        9,
    ]
    assert [row["antiunitary"] for row in gamma6["generator_matrices"]] == [
        False,
        False,
        True,
    ]
    assert gamma6["generator_matrices"][0]["matrix"] == [
        ["sqrt(2)/2", "sqrt(2)/2"],
        ["-sqrt(2)/2", "sqrt(2)/2"],
    ]
    assert gamma6["generator_matrices"][2]["matrix"] == [
        ["0", "1"],
        ["-1", "0"],
    ]
    assert r"\begin{matrix}" in gamma6["generator_matrices"][0]["matrix_latex"]


def test_catalog_builds_generator_input_and_direct_sums():
    core_input, metadata, group_input = build_catalog_core_input(
        space_group=92,
        magnetic_group="92.112",
        k_name="Γ",
        representation_family="double",
        representations=[{"index": 1}, {"index": 2, "multiplicity": 2}],
        coordinate_basis="cartesian",
        method="iterative",
    )
    assert group_input == "generators"
    assert metadata["dimension"] == 6
    assert metadata["constraint_count"] < metadata["group_order"]
    assert metadata["momentum_basis"] == "cartesian"
    assert core_input["KPointLabel"] == "Γ"
    assert core_input["CoordinateSystem"] == "Cartesian"
    assert core_input["Group"] == (92, 112)
    matrices = [
        item[0]
        for category_name in ("Unitary", "Antiunitary")
        for item in core_input.get(category_name, {}).values()
    ]
    assert matrices
    assert all(matrix.shape == (6, 6) for matrix in matrices)
    assert "Antiunitary" in core_input


def test_catalog_resolves_exact_line_parameters():
    _core_input, metadata, _group_input = build_catalog_core_input(
        space_group=4,
        magnetic_group="4.10",
        k_name="Λ",
        parameter_values={"u": "1/5"},
        representation_family="double",
        representations=[{"index": 1}],
        coordinate_basis="fractional",
        method="iterative",
    )
    assert metadata["coordinate"]["coordinate"] == ["0", "0", "1/5"]
    assert metadata["coordinate"]["parameter_values"] == {"u": "1/5"}


def test_cartesian_and_fractional_momentum_actions_follow_msgcorep():
    cartesian, cartesian_metadata, _ = build_catalog_core_input(
        space_group=156,
        magnetic_group="156.50",
        k_name="K",
        representation_family="single",
        representations=[{"index": 1}],
        coordinate_basis="cartesian",
        method="iterative",
    )
    fractional, fractional_metadata, _ = build_catalog_core_input(
        space_group=156,
        magnetic_group="156.50",
        k_name="K",
        representation_family="single",
        representations=[{"index": 1}],
        coordinate_basis="fractional",
        method="iterative",
    )
    assert cartesian_metadata["momentum_basis"] == "cartesian"
    assert fractional_metadata["momentum_basis"] == "fractional"
    for category in ("Unitary", "Antiunitary"):
        operations = cartesian.get(category, {})
        for label, (_representation, momentum) in operations.items():
            operation = next(
                row for row in cartesian_metadata["operations"] if row["label"] == label
            )
            name = operation["rotation"].removeprefix("bar")
            sign = -1 if category == "Antiunitary" else 1
            assert momentum == sign * Matrix(exact_to_sympy(msgcorep.RotMatCart[name]))
    for category in ("Unitary", "Antiunitary"):
        operations = fractional.get(category, {})
        for label, (_representation, momentum) in operations.items():
            operation = next(
                row for row in fractional_metadata["operations"] if row["label"] == label
            )
            name = operation["rotation"].removeprefix("bar")
            sign = -1 if category == "Antiunitary" else 1
            assert momentum == sign * Matrix(
                exact_to_sympy(msgcorep.getRotMatOfK("HexaPrim", name))
            )


def test_catalog_request_produces_the_expected_gamma_model():
    result = compute_catalog_request(_default_request())
    assert result["status"] == "ok"
    assert result["dimension"] == 2
    assert result["parameters"] == ["e_1", "t_1", "t_2"]
    assert [row["invariant_dimension"] for row in result["diagnostics"]] == [1, 2]
    assert result["selection"]["magnetic_group"] == "92.112"
    assert result["selection"]["constraint_count"] < result["selection"]["group_order"]


def test_reynolds_catalog_input_uses_all_elements_and_agrees_in_dimension():
    iterative = compute_catalog_request(_default_request(method="iterative"))
    reynolds = compute_catalog_request(_default_request(method="reynolds"))
    assert reynolds["selection"]["group_input"] == "elements"
    assert reynolds["selection"]["constraint_count"] == reynolds["selection"]["group_order"]
    assert reynolds["number_of_parameters"] == iterative["number_of_parameters"]
    assert [row["invariant_dimension"] for row in reynolds["diagnostics"]] == [1, 2]


def test_catalog_errors_survive_the_worker_boundary_with_their_type():
    request = _default_request(k_name="not-a-k-point")
    data = request.model_dump() if hasattr(request, "model_dump") else request.dict()
    with pytest.raises(CatalogInputError):
        run_catalog_isolated(data, timeout=20)


def test_catalog_http_workflow():
    testclient = pytest.importorskip("fastapi.testclient")
    with testclient.TestClient(app) as client:
        groups = client.get("/api/catalog/space-groups")
        assert groups.status_code == 200
        assert groups.json()["count"] == 230

        points = client.get(
            "/api/catalog/k-points/92",
            params={"magnetic_group": "92.112", "variant": "a"},
        )
        assert points.status_code == 200
        assert any(row["name"] == "Γ" for row in points.json()["points"])

        coreps = client.get(
            "/api/catalog/corepresentations",
            params={
                "magnetic_group": "92.112",
                "k_name": "Γ",
                "variant": "a",
                "coordinate_basis": "cartesian",
            },
        )
        assert coreps.status_code == 200
        assert len(coreps.json()["generators"]) == 3

        result = client.post(
            "/api/catalog/solve",
            json=(
                _default_request().model_dump()
                if hasattr(_default_request(), "model_dump")
                else _default_request().dict()
            ),
        )
        assert result.status_code == 200
        assert result.json()["parameters"] == ["e_1", "t_1", "t_2"]


def test_catalog_rejects_inconsistent_selections():
    with pytest.raises(CatalogInputError, match="does not belong"):
        build_catalog_core_input(
            space_group=91,
            magnetic_group="92.112",
            k_name="Γ",
            representation_family="double",
            representations=[{"index": 1}],
        )
    with pytest.raises(CatalogInputError, match="coordinate basis"):
        build_catalog_core_input(
            space_group=92,
            magnetic_group="92.112",
            k_name="Γ",
            representation_family="double",
            representations=[{"index": 1}],
            coordinate_basis="oblique",
        )
