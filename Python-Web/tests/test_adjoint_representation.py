from __future__ import annotations

import pytest
from sympy import I, Matrix, sqrt

from conftest import (
    basis_list,
    call_adjoint,
    call_operator_basis,
)


def legacy_adjoint(representation, basis, antiunitary):
    inverse = representation.inv()
    images = [
        representation
        * (matrix.conjugate() if antiunitary else matrix)
        * inverse
        for matrix in basis
    ]
    return Matrix(
        len(basis),
        len(basis),
        lambda j, i: (
            (basis[j].H * images[i]).trace()
            / (basis[j].H * basis[j]).trace()
        ).simplify(),
    )


def exact_representation(n):
    return Matrix(
        n,
        n,
        lambda i, j: (-1) ** i if j == n - 1 - i else 0,
    )


def test_adjoint_old_loop_parity_n1_n10(
    operator_basis_callable,
    adjoint_callable,
):
    for n in range(1, 11):
        data = call_operator_basis(operator_basis_callable, n)
        basis = basis_list(data)
        representation = exact_representation(n)
        for antiunitary in (False, True):
            computed = Matrix(
                call_adjoint(
                    adjoint_callable,
                    representation,
                    data,
                    antiunitary,
                )
            )
            expected = legacy_adjoint(
                representation, basis, antiunitary
            )
            assert computed.equals(expected)


def test_adjoint_reconstruction_and_composition(
    operator_basis_callable,
    adjoint_callable,
):
    data = call_operator_basis(operator_basis_callable, 2)
    basis = basis_list(data)
    d1 = Matrix([[1, 0], [0, I]])
    d2 = Matrix([[1, 1], [-1, 1]]) / sqrt(2)

    j1 = Matrix(call_adjoint(adjoint_callable, d1, data, False))
    j2 = Matrix(call_adjoint(adjoint_callable, d2, data, False))
    j12 = Matrix(
        call_adjoint(adjoint_callable, d1 * d2, data, False)
    )
    assert j12.equals(j1 * j2)

    ja1 = Matrix(call_adjoint(adjoint_callable, d1, data, True))
    ja2 = Matrix(call_adjoint(adjoint_callable, d2, data, True))
    expected_aa = Matrix(
        call_adjoint(
            adjoint_callable,
            d1 * d2.conjugate(),
            data,
            False,
        )
    )
    assert (ja1 * ja2).equals(expected_aa)

    for antiunitary in (False, True):
        j = Matrix(
            call_adjoint(
                adjoint_callable, d2, data, antiunitary
            )
        )
        for i, matrix in enumerate(basis):
            image = (
                d2
                * (matrix.conjugate() if antiunitary else matrix)
                * d2.inv()
            )
            reconstruction = sum(
                (j[q, i] * basis[q] for q in range(len(basis))),
                Matrix.zeros(2),
            )
            assert image.equals(reconstruction)


def test_adjoint_projective_sign_invariance(
    operator_basis_callable,
    adjoint_callable,
):
    data = call_operator_basis(operator_basis_callable, 2)
    d = Matrix([[1, 0], [0, I]])
    positive = Matrix(call_adjoint(adjoint_callable, d, data, False))
    negative = Matrix(call_adjoint(adjoint_callable, -d, data, False))
    assert positive.equals(negative)
