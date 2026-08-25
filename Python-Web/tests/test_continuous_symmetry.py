from __future__ import annotations

import pytest
from sympy import I, Matrix, exp, sqrt, symbols


def coefficient_space(hamiltonian, momentum_symbols):
    parameters = sorted(
        hamiltonian.free_symbols - set(momentum_symbols),
        key=str,
    )
    rows = []
    for parameter in parameters:
        coefficient = hamiltonian.applyfunc(
            lambda entry: entry.coeff(parameter)
        )
        rows.append(list(coefficient))
    return Matrix(rows) if rows else Matrix.zeros(0, hamiltonian.rows**2)


def same_coefficient_space(left, right, momentum_symbols):
    left_space = coefficient_space(left, momentum_symbols)
    right_space = coefficient_space(right, momentum_symbols)
    return (
        left_space.rank()
        == right_space.rank()
        == Matrix.vstack(left_space, right_space).rank()
    )


@pytest.fixture
def continuous_data():
    from magnetickp.magkp import kx, ky, kz

    theta = symbols("theta")
    representation = Matrix.diag(
        exp(-I * theta / 2),
        exp(I * theta / 2),
    )
    variables = (kx, ky, kz)
    return theta, representation, variables


def test_continuous_only_reduces_operator_basis(continuous_data):
    import magnetickp

    theta, representation, variables = continuous_data
    model = magnetickp.kpHam(
        {
            "Unitary": {
                "CInfinity": (representation, variables, theta)
            }
        }
    )
    result = model.getkpHam(1, return_diagnostics=True)
    assert result["NumberOfParameters"] == 8
    assert [
        diagnostic["invariant_dimension"]
        for diagnostic in result["Diagnostics"]
    ] == [2, 6]
    assert all(
        diagnostic["continuous_symmetry"] == "CInfinity"
        and diagnostic["original_operator_dimension"] == 4
        and diagnostic["reduced_operator_dimension"] == 2
        and diagnostic["solver"] == "none"
        for diagnostic in result["Diagnostics"]
    )

    constant, count = model.IterativelySimplify(0)
    assert count == 2
    assert constant[0, 1] == constant[1, 0] == 0
    assert coefficient_space(constant, variables).rank() == 2


def test_low_order_coefficients_use_named_one_based_families(continuous_data):
    import magnetickp

    theta, representation, variables = continuous_data
    model = magnetickp.kpHam(
        {
            "Unitary": {
                "CInfinity": (representation, variables, theta)
            }
        }
    )
    hamiltonian = model.getkpHam(3)["ham"]
    names = {
        str(parameter)
        for parameter in hamiltonian.free_symbols - set(variables)
    }
    expected = {
        *(f"e_{index}" for index in range(1, 3)),
        *(f"t_{index}" for index in range(1, 7)),
        *(f"r_{index}" for index in range(1, 13)),
        *(f"s_{index}" for index in range(1, 21)),
    }
    assert names == expected

    fourth_order, count = model.IterativelySimplify(4)
    fourth_order_names = {
        str(parameter)
        for parameter in fourth_order.free_symbols - set(variables)
    }
    assert count == 30
    assert fourth_order_names == {f"C_{index}_4" for index in range(count)}


def test_continuous_equals_dense_irrational_cyclic_subgroup():
    """Cross-check CInfinity against one exact dense cyclic generator."""
    import magnetickp
    from magnetickp.magkp import kx, ky, kz

    theta = symbols("theta_dense")
    variables = (kx, ky, kz)
    representation = Matrix.diag(exp(-I * theta), 1, exp(I * theta))
    continuous_model = magnetickp.kpHam(
        {
            "Unitary": {
                "CInfinity": (representation, variables, theta)
            }
        }
    )

    # This is theta0 = acos(3/5).  Niven's theorem implies theta0/pi is
    # irrational, while exact rational-complex entries keep the test fast.
    dense_rotation = Matrix.diag((3 - 4 * I) / 5, 1, (3 + 4 * I) / 5)
    dense_model = magnetickp.kpHam(
        {
            "Unitary": {
                "DenseCyclicGenerator": (dense_rotation, variables)
            }
        }
    )

    # At the rational angle pi the outer spin-one weights are degenerate.
    rational_model = magnetickp.kpHam(
        {
            "Unitary": {
                "OrderTwoRotation": (Matrix.diag(-1, 1, -1), variables)
            }
        }
    )

    # Python IterativelySimplify returns one homogeneous order at a time.
    expected = ((3, 3, 5), (9, 9, 15), (18, 18, 30))
    for order, expected_counts in enumerate(expected):
        continuous_ham, continuous_count = (
            continuous_model.IterativelySimplify(order)
        )
        dense_ham, dense_count = dense_model.IterativelySimplify(order)
        _, rational_count = rational_model.IterativelySimplify(order)
        assert (continuous_count, dense_count, rational_count) == (
            expected_counts
        )
        assert same_coefficient_space(
            continuous_ham,
            dense_ham,
            variables,
        )


def test_continuous_then_discrete_normalizer(continuous_data):
    import magnetickp

    theta, representation, variables = continuous_data
    model = magnetickp.kpHam(
        {
            "Unitary": {
                "Cinfty": (representation, variables, theta),
                "SpinFlip": (Matrix([[0, 1], [1, 0]]), variables),
            }
        }
    )
    iterative_ham, iterative_count = model.IterativelySimplify(
        0,
        method="iterative",
    )
    reynolds_ham, reynolds_count, diagnostics = model.IterativelySimplify(
        0,
        method="reynolds",
        group_input="generators",
        return_diagnostics=True,
    )
    assert iterative_count == reynolds_count == 1
    assert diagnostics["fallback"] is False
    assert same_coefficient_space(iterative_ham, reynolds_ham, variables)
    parameter = next(iter(iterative_ham.free_symbols))
    assert iterative_ham.diff(parameter).equals(Matrix.eye(2))


def test_continuous_rejects_invalid_inputs(continuous_data):
    import magnetickp

    theta, representation, variables = continuous_data
    bad_normalizer = Matrix([[1, 1], [1, -1]]) / sqrt(2)
    with pytest.raises(ValueError, match="does not normalize"):
        magnetickp.kpHam(
            {
                "Unitary": {
                    "CInfinity": (representation, variables, theta),
                    "Bad": (bad_normalizer, variables),
                }
            }
        )

    with pytest.raises(ValueError, match="leave momentum unchanged"):
        magnetickp.kpHam(
            {
                "Unitary": {
                    "CInfinity": (
                        representation,
                        (-variables[0], variables[1], variables[2]),
                        theta,
                    )
                }
            }
        )

    extra = symbols("m")
    representation_with_extra = Matrix.diag(
        exp(-I * extra * theta),
        exp(I * extra * theta),
    )
    with pytest.raises(ValueError, match="only the declared parameter"):
        magnetickp.kpHam(
            {
                "Unitary": {
                    "CInfinity": (
                        representation_with_extra,
                        variables,
                        theta,
                    )
                }
            }
        )

    with pytest.raises(ValueError, match="only in the Unitary"):
        magnetickp.kpHam(
            {
                "Antiunitary": {
                    "CInfinity": (representation, variables, theta)
                }
            }
        )

    nonrepresentation = Matrix.diag(
        exp(-I * theta**2 / 2),
        exp(I * theta**2 / 2),
    )
    with pytest.raises(ValueError, match="theta-independent generator"):
        magnetickp.kpHam(
            {
                "Unitary": {
                    "CInfinity": (nonrepresentation, variables, theta)
                }
            }
        )
