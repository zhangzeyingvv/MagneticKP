"""Homogeneous-monomial bases and direct momentum transformation matrices."""

from __future__ import annotations

from itertools import permutations
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from sympy import Matrix, Poly, SparseMatrix, expand, prod, symbols
from sympy.matrices import MatrixBase


kx, ky, kz = symbols("kx ky kz")
K_SYMBOLS = (kx, ky, kz)
_MOMENTUM_CACHE: Dict[tuple, MatrixBase] = {}


def _integer_partitions(n: int, length: int, largest: int):
    if n == 0:
        yield [0] * length
    elif length == 0 or n < 0:
        return
    else:
        for value in range(largest, 0, -1):
            for tail in _integer_partitions(n - value, length - 1, value):
                yield [value] + tail


def _powers(expr, variables: Sequence) -> Tuple[int, ...]:
    powers = expr.as_powers_dict()
    return tuple(int(powers.get(variable, 0)) for variable in variables)


def homogeneous_monomials(order: int, dimension: int = 3) -> List:
    """Return the legacy MagneticKP monomial ordering."""

    if not isinstance(order, int) or order < 0 or dimension not in (1, 2, 3):
        raise ValueError("order must be non-negative and dimension must be 1, 2 or 3")
    variables = (
        (kz,)
        if dimension == 1
        else (kx, ky)
        if dimension == 2
        else (kx, ky, kz)
    )
    if dimension == 1:
        return [kz**order]
    monomials = set()
    for partition in _integer_partitions(order, dimension, order):
        for exponent in set(permutations(partition)):
            monomials.add(prod(v**p for v, p in zip(variables, exponent)))
    monomials = list(monomials)
    if dimension == 2:
        monomials.sort(key=lambda expr: (_powers(expr, variables)[1], _powers(expr, variables)[0]))
    else:
        monomials.sort(
            key=lambda expr: (
                _powers(expr, variables)[2],
                _powers(expr, variables)[1],
                _powers(expr, variables)[0],
            )
        )
    return monomials


def _is_zero(value) -> bool:
    return value == 0 or value.is_zero is True


def _linear_matrix(k_image, variables: Sequence) -> Tuple[Matrix, Matrix]:
    dimension = len(variables)
    if isinstance(k_image, MatrixBase) and k_image.shape == (dimension, dimension):
        matrix = Matrix(k_image)
        image = matrix * Matrix(variables)
        return matrix, image
    image = Matrix(list(k_image))
    if image.shape != (dimension, 1):
        raise ValueError("the momentum image has the wrong dimension")
    matrix = Matrix(
        dimension,
        dimension,
        lambda i, j: expand(image[i]).coeff(variables[j]),
    )
    residual = image - matrix * Matrix(variables)
    if any(not _is_zero(expand(value)) for value in residual):
        raise ValueError("the momentum image must be linear and homogeneous")
    return matrix, image


def _monomial_linear_data(forms: Matrix, variables: Sequence):
    matrix = Matrix(
        len(variables),
        len(variables),
        lambda i, j: expand(forms[i]).coeff(variables[j]),
    )
    positions = []
    scales = []
    for row in range(matrix.rows):
        nonzero = [column for column in range(matrix.cols) if not _is_zero(matrix[row, column])]
        if len(nonzero) != 1:
            return None
        positions.append(nonzero[0])
        scales.append(matrix[row, nonzero[0]])
    if len(set(positions)) != len(positions):
        return None
    return positions, scales


def momentum_transform_matrix(
    order: int,
    k_image,
    variables: Sequence = K_SYMBOLS,
    klist: Optional[Sequence] = None,
) -> MatrixBase:
    """Construct ``F = T(A^-1).T`` directly from homogeneous coefficients."""

    variables = tuple(variables)
    if len(variables) not in (1, 2, 3):
        raise ValueError("the spatial dimension must be 1, 2 or 3")
    monomials = list(
        homogeneous_monomials(order, len(variables)) if klist is None else klist
    )
    matrix, _ = _linear_matrix(k_image, variables)
    key = (
        order,
        tuple(monomials),
        tuple(matrix),
        variables,
    )
    cached = _MOMENTUM_CACHE.get(key)
    if cached is not None:
        return cached
    if matrix.det() == 0:
        raise ValueError("the momentum transformation must be invertible")
    if matrix == SparseMatrix.eye(len(variables)) or matrix == Matrix.eye(len(variables)):
        result = SparseMatrix.eye(len(monomials))
        _MOMENTUM_CACHE[key] = result
        return result

    inverse_k = matrix.inv() * Matrix(variables)
    exponents = [_powers(monomial, variables) for monomial in monomials]
    position = {exponent: index for index, exponent in enumerate(exponents)}
    monomial_data = _monomial_linear_data(inverse_k, variables)
    entries = {}
    if monomial_data is not None:
        permutation, scales = monomial_data
        for row, alpha in enumerate(exponents):
            beta = [0] * len(variables)
            coefficient = 1
            for index, power in enumerate(alpha):
                beta[permutation[index]] += power
                coefficient *= scales[index] ** power
            entries[(row, position[tuple(beta)])] = coefficient
    else:
        for row, alpha in enumerate(exponents):
            polynomial = Poly(
                expand(prod(inverse_k[index] ** alpha[index] for index in range(len(variables)))),
                *variables,
            )
            for beta, coefficient in polynomial.terms():
                entries[(row, position[tuple(beta)])] = coefficient
    result = SparseMatrix(len(monomials), len(monomials), entries).T
    _MOMENTUM_CACHE[key] = result
    return result
