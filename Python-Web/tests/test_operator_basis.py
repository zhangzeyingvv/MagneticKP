from __future__ import annotations

from sympy import Matrix

from conftest import (
    basis_list,
    call_operator_basis,
    gram_matrix,
)


def orthogonal_hilbert_schmidt_gram(basis):
    return Matrix.diag(
        *((matrix.H * matrix).trace() for matrix in basis)
    )


def test_operator_basis_dimensions_1_through_10(
    operator_basis_callable,
):
    for n in range(1, 11):
        data = call_operator_basis(operator_basis_callable, n)
        basis = basis_list(data)
        assert len(basis) == n * n
        assert all(matrix.shape == (n, n) for matrix in basis)
        assert all(
            matrix.equals(matrix.H)
            for matrix in basis
        )
        flattened = Matrix.hstack(
            *(Matrix(matrix).reshape(n * n, 1) for matrix in basis)
        )
        assert flattened.rank() == n * n
        gram = gram_matrix(data)
        expected = orthogonal_hilbert_schmidt_gram(basis)
        if gram is not None:
            assert Matrix(gram).equals(expected)
        assert expected.rank() == n * n


def test_operator_basis_generic_dimensions(operator_basis_callable):
    counts = {}
    for n in (5, 7, 9, 10):
        data = call_operator_basis(operator_basis_callable, n)
        counts[n] = len(basis_list(data))
    assert counts == {5: 25, 7: 49, 9: 81, 10: 100}
