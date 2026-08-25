from __future__ import annotations

from sympy import Matrix, Poly, expand, symbols

from magnetickp.solvers import (
    iterative_invariant_basis,
    reynolds_invariant_basis,
    stacked_invariant_basis,
)


kx, ky, kz = symbols("kx ky kz")
K_VARIABLES = (kx, ky, kz)


def same_column_space(left, right):
    assert left.rows == right.rows
    return (
        left.rank()
        == right.rank()
        == Matrix.hstack(left, right).rank()
    )


def hamiltonian_coordinate_rows(hamiltonians):
    parameters = sorted(
        set().union(*(ham.free_symbols for ham in hamiltonians))
        - set(K_VARIABLES),
        key=str,
    )
    monomials = set()
    coefficient_matrices = []
    for hamiltonian in hamiltonians:
        matrices = [
            hamiltonian.applyfunc(
                lambda entry, parameter=parameter: expand(
                    entry.coeff(parameter)
                )
            )
            for parameter in parameters
            if parameter in hamiltonian.free_symbols
        ]
        coefficient_matrices.append(matrices)
        for matrix in matrices:
            for entry in matrix:
                monomials.update(
                    Poly(entry, *K_VARIABLES).monoms()
                )
    monomials = sorted(monomials)

    rows_by_hamiltonian = []
    for matrices in coefficient_matrices:
        rows = []
        for matrix in matrices:
            rows.append(
                [
                    Poly(entry, *K_VARIABLES).coeff_monomial(
                        monomial
                    )
                    for monomial in monomials
                    for entry in matrix
                ]
            )
        width = (
            len(monomials) * hamiltonians[0].rows**2
        )
        rows_by_hamiltonian.append(
            Matrix(rows) if rows else Matrix.zeros(0, width)
        )
    return rows_by_hamiltonian


def same_hamiltonian_subspace(left, right):
    left_rows, right_rows = hamiltonian_coordinate_rows(
        [left, right]
    )
    return (
        left_rows.rank()
        == right_rows.rank()
        == Matrix.vstack(left_rows, right_rows).rank()
    )


def test_stacked_iterative_reynolds_factor_subspaces_agree():
    elements = [
        {
            "F": Matrix.eye(3),
            "J": Matrix([[1]]),
            "eta": 1,
        },
        {
            "F": Matrix.diag(1, -1, -1),
            "J": Matrix([[1]]),
            "eta": 1,
        },
    ]
    iterative = iterative_invariant_basis(elements)
    stacked = stacked_invariant_basis(elements)
    reynolds, diagnostics = reynolds_invariant_basis(
        elements,
        group_input="elements",
    )
    assert diagnostics["fallback"] is False
    assert same_column_space(iterative, stacked)
    assert same_column_space(iterative, reynolds)


def test_type3_antisymmetry_direct_product_dispatch():
    import magnetickp

    sigma_z = Matrix.diag(1, -1)
    input_data = {
        "AntisymmetryUnitaryTest": {
            "type3": (sigma_z, Matrix.eye(3))
        }
    }
    model = magnetickp.kpHam(input_data)

    iterative_ham, iterative_count = model.IterativelySimplify(
        0,
        method="iterative",
    )
    direct_ham, direct_count = model.DirectProductDecomposition(0)
    alias_ham, alias_count = model.IterativelySimplify(
        0,
        method="stacked",
    )
    reynolds_ham, reynolds_count, diagnostics = (
        model.IterativelySimplify(
            0,
            method="reynolds",
            group_input="generators",
            return_diagnostics=True,
        )
    )

    assert {
        iterative_count,
        direct_count,
        alias_count,
        reynolds_count,
    } == {2}
    assert diagnostics["fallback"] is False
    assert same_hamiltonian_subspace(iterative_ham, direct_ham)
    assert same_hamiltonian_subspace(iterative_ham, alias_ham)
    assert same_hamiltonian_subspace(iterative_ham, reynolds_ham)


def test_automatic_large_ambient_selects_iterative():
    import magnetickp

    cyclic_five = Matrix(
        5,
        5,
        lambda row, column: (
            1 if column == (row + 1) % 5 else 0
        ),
    )
    model = magnetickp.kpHam(
        {
            "Unitary": {
                "C5": (cyclic_five, Matrix.eye(3))
            }
        }
    )
    automatic_ham, automatic_count, automatic_diagnostics = (
        model.IterativelySimplify(
            2,
            method="automatic",
            return_diagnostics=True,
        )
    )
    iterative_ham, iterative_count, iterative_diagnostics = (
        model.IterativelySimplify(
            2,
            method="iterative",
            return_diagnostics=True,
        )
    )

    assert automatic_diagnostics["ambient_dimension"] > 128
    assert automatic_diagnostics["solver"] == "iterative"
    assert iterative_diagnostics["solver"] == "iterative"
    assert automatic_count == iterative_count
    assert same_hamiltonian_subspace(
        automatic_ham,
        iterative_ham,
    )
