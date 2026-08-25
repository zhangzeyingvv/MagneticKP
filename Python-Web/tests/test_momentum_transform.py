from __future__ import annotations

import io
import math
import pathlib
import pickle
import subprocess

from sympy import Matrix, Poly, symbols

from conftest import (
    REPOSITORY_ROOT,
    call_momentum_transform,
)


kx, ky, kz = symbols("kx ky kz")


def monomial_basis(order: int, dimension: int):
    variables = (kx, ky, kz)[:dimension]
    exponents = []

    def visit(prefix, remaining, slots):
        if slots == 1:
            exponents.append(tuple(prefix + [remaining]))
            return
        for value in range(remaining + 1):
            visit(prefix + [value], remaining - value, slots - 1)

    visit([], order, dimension)
    exponents.sort(
        key=lambda exponent: tuple(reversed(exponent))
    )
    return [
        Poly.from_dict({exponent: 1}, variables).as_expr()
        for exponent in exponents
    ], variables


def load_legacy_fd_pickle():
    result = subprocess.run(
        [
            "git",
            "show",
            "a319488:MagneticKP-python/magnetickp/fd.pkl",
        ],
        cwd=REPOSITORY_ROOT,
        check=True,
        stdout=subprocess.PIPE,
    )
    return pickle.load(io.BytesIO(result.stdout))


def matrix_exact_equal(left, right):
    left = Matrix(left)
    right = Matrix(right)
    return left.shape == right.shape and left.equals(right)


def test_all_legacy_fd_entries(momentum_transform_callable):
    legacy = load_legacy_fd_pickle()
    assert len(legacy) > 0
    mismatches = []
    for key, expected in legacy.items():
        order = int(key[-1])
        flat = key[:-1]
        dimension = math.isqrt(len(flat))
        assert dimension * dimension == len(flat)
        spatial = Matrix(dimension, dimension, flat)
        basis, variables = monomial_basis(order, dimension)
        image = spatial * Matrix(variables)
        computed = call_momentum_transform(
            momentum_transform_callable,
            order,
            image,
            variables,
            basis,
        )
        if not matrix_exact_equal(computed, expected):
            mismatches.append(key)
    assert mismatches == []


def test_momentum_transform_definition_and_homomorphism(
    momentum_transform_callable,
):
    a = Matrix([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    b = Matrix([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
    for order in range(0, 9):
        basis, variables = monomial_basis(order, 3)

        def transform(matrix):
            return Matrix(
                call_momentum_transform(
                    momentum_transform_callable,
                    order,
                    matrix * Matrix(variables),
                    variables,
                    basis,
                )
            )

        f_a = transform(a)
        f_b = transform(b)
        f_ab = transform(a * b)
        assert f_ab.equals(f_a * f_b)

        transformed_basis = Matrix(
            [
                value.subs(
                    dict(zip(variables, a * Matrix(variables))),
                    simultaneous=True,
                )
                for value in basis
            ]
        )
        assert Matrix([basis]).equals(transformed_basis.T * f_a)
