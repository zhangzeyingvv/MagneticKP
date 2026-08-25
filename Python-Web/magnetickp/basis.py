"""Operator bases and adjoint-representation matrices for MagneticKP."""

from __future__ import annotations

from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from sympy import (
    I,
    Matrix,
    SparseMatrix,
    Symbol,
    conjugate,
    eye,
    kronecker_product,
    simplify,
    sqrt,
)
from sympy.matrices import MatrixBase


_BASIS_CACHE: Dict[int, dict] = {}
_ADJOINT_CACHE: Dict[Tuple[tuple, tuple], Matrix] = {}


def _pauli_basis() -> List[Matrix]:
    return [
        Matrix([[1, 0], [0, 1]]),
        Matrix([[0, 1], [1, 0]]),
        Matrix([[0, -I], [I, 0]]),
        Matrix([[1, 0], [0, -1]]),
    ]


def _gell_mann_basis() -> List[Matrix]:
    return [
        Matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
        Matrix([[0, -I, 0], [I, 0, 0], [0, 0, 0]]),
        Matrix([[0, 0, -I], [0, 0, 0], [I, 0, 0]]),
        Matrix([[0, 0, 0], [0, 0, -I], [0, I, 0]]),
        Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]]),
        Matrix([[1, 0, 0], [0, -1, 0], [0, 0, 0]]),
        Matrix([[0, 0, 1], [0, 0, 0], [1, 0, 0]]),
        Matrix([[0, 0, 0], [0, 0, 1], [0, 1, 0]]),
        Matrix([[1, 0, 0], [0, 1, 0], [0, 0, -2]]) / sqrt(3),
    ]


def _matrix_unit(n: int, row: int, column: int) -> SparseMatrix:
    return SparseMatrix(n, n, {(row, column): 1})


def operator_basis(n: int) -> List[MatrixBase]:
    """Return the legacy-compatible Hermitian basis for dimension ``n``.

    Dimensions 1, 2, 3, 4, 6 and 8 retain the original basis and ordering.
    Other positive dimensions use the root-free sparse basis
    ``E_aa``, ``E_ab + E_ba`` and ``-I E_ab + I E_ba``.
    """

    if not isinstance(n, int) or n <= 0:
        raise ValueError("the representation dimension must be a positive integer")
    pauli = _pauli_basis()
    gm = _gell_mann_basis()
    if n == 1:
        return [eye(1)]
    if n == 2:
        return pauli
    if n == 3:
        return gm
    if n == 4:
        return [
            kronecker_product(pauli[i], pauli[j])
            for i in range(4)
            for j in range(4)
        ]
    if n == 6:
        return [
            kronecker_product(pauli[i], gm[j])
            for i in range(4)
            for j in range(9)
        ]
    if n == 8:
        return [
            kronecker_product(pauli[i], pauli[j], pauli[k])
            for i in range(4)
            for j in range(4)
            for k in range(4)
        ]

    diagonal = [_matrix_unit(n, a, a) for a in range(n)]
    symmetric = [
        _matrix_unit(n, a, b) + _matrix_unit(n, b, a)
        for a in range(n - 1)
        for b in range(a + 1, n)
    ]
    antisymmetric = [
        -I * _matrix_unit(n, a, b) + I * _matrix_unit(n, b, a)
        for a in range(n - 1)
        for b in range(a + 1, n)
    ]
    result = diagonal + symmetric + antisymmetric
    if len(result) != n * n:
        raise RuntimeError("internal operator-basis count mismatch")
    return result


def _row_flatten(matrix: MatrixBase) -> List:
    return [matrix[i, j] for i in range(matrix.rows) for j in range(matrix.cols)]


def _basis_cache_key(basis: Sequence[MatrixBase]) -> tuple:
    return tuple(tuple(_row_flatten(matrix)) for matrix in basis)


def basis_data(n: int, basis: Optional[Sequence[MatrixBase]] = None) -> dict:
    """Build Gram, dual and complex-conjugation data for an operator basis."""

    default_basis = basis is None
    if default_basis and n in _BASIS_CACHE:
        return _BASIS_CACHE[n]
    matrices = list(operator_basis(n) if default_basis else basis)
    if not matrices:
        raise ValueError("the operator basis cannot be empty")
    if any(matrix.shape != (n, n) for matrix in matrices):
        raise ValueError("all operator-basis matrices must have shape (n, n)")

    rows = Matrix([_row_flatten(matrix) for matrix in matrices])
    columns = rows.T
    if default_basis:
        # Every built-in basis is Hilbert--Schmidt orthogonal.  Constructing
        # its diagonal Gram matrix avoids a symbolic rank/inverse of size n^2.
        norms = [
            simplify((matrix.H * matrix).trace())
            for matrix in matrices
        ]
        gram = SparseMatrix.diag(*norms)
        gram_inverse = SparseMatrix.diag(
            *[simplify(1 / norm) for norm in norms]
        )
        dual = gram_inverse * columns.conjugate().T
        signs = []
        for matrix in matrices:
            conjugate = matrix.conjugate()
            if conjugate == matrix:
                signs.append(1)
            elif conjugate == -matrix:
                signs.append(-1)
            else:
                signs = []
                break
        conjugation = (
            SparseMatrix.diag(*signs)
            if signs
            else simplify(dual * columns.conjugate())
        )
    else:
        gram = columns.conjugate().T * columns
        if gram.rank() != len(matrices):
            raise ValueError("the supplied operator matrices are linearly dependent")
        if gram.is_diagonal():
            gram_inverse = SparseMatrix.diag(
                *[simplify(1 / gram[i, i]) for i in range(gram.rows)]
            )
            dual = gram_inverse * columns.conjugate().T
        else:
            dual = gram.inv() * columns.conjugate().T
        conjugation = simplify(dual * columns.conjugate())
    data = {
        "dimension": n,
        "count": len(matrices),
        "basis": matrices,
        "rows": rows,
        "columns": columns,
        "gram": gram,
        "dual": dual,
        "conjugation": conjugation,
        "cache_key": (
            ("default", n)
            if default_basis
            else _basis_cache_key(matrices)
        ),
    }
    if default_basis:
        _BASIS_CACHE[n] = data
    return data


def _first_nonzero(entries: Iterable):
    for value in entries:
        if value == 0 or value.is_zero is True:
            continue
        return value
    return None


def _projective_representative(matrix: MatrixBase) -> Matrix:
    pivot = _first_nonzero(_row_flatten(matrix))
    if pivot is None:
        raise ValueError("the representation matrix must be invertible")
    return matrix.applyfunc(lambda value: simplify(value / pivot))


def adjoint_representation_matrix(
    representation: MatrixBase,
    data: dict,
    antiunitary: bool = False,
) -> Matrix:
    """Return the induced matrix on the operator basis.

    Columns follow ``D M_i D^-1 = sum_j J[j,i] M_j``.  The unitary part
    is cached projectively, so ``D`` and any exact non-zero scalar multiple
    reuse the same result.  Antiunitary action is obtained by right
    multiplication with the precomputed complex-conjugation matrix.
    """

    representative = _projective_representative(representation)
    key = (
        data["cache_key"],
        tuple(_row_flatten(representative)),
    )
    unitary = _ADJOINT_CACHE.get(key)
    if unitary is None:
        inverse = representation.inv()
        images = [
            representation * matrix * inverse
            for matrix in data["basis"]
        ]
        image_columns = Matrix([_row_flatten(matrix) for matrix in images]).T
        unitary = simplify(data["dual"] * image_columns)
        _ADJOINT_CACHE[key] = unitary
    if antiunitary:
        return simplify(unitary * data["conjugation"])
    return unitary


def _simplify_matrix(matrix: MatrixBase) -> Matrix:
    return Matrix(matrix).applyfunc(simplify)


def _zero_matrix(matrix: MatrixBase) -> bool:
    simplified = _simplify_matrix(matrix)
    return all(
        value == 0 or value.is_zero is True
        for value in simplified
    )


def continuous_generator(
    representation: MatrixBase,
    parameter: Symbol,
) -> dict:
    """Extract the Hermitian generator of an internal one-parameter group.

    The declared ``parameter`` is treated as real and must be the only free
    user symbol in the representation matrix.  A scalar phase at the identity
    is accepted; the returned working generator is made traceless because
    scalar generators act trivially on the operator commutant.
    """

    representation = Matrix(representation)
    if representation.rows != representation.cols:
        raise ValueError("CInfinity representation must be square")
    if not isinstance(parameter, Symbol):
        raise ValueError("the CInfinity third entry must be a symbol")
    extra_symbols = representation.free_symbols - {parameter}
    if extra_symbols:
        raise ValueError(
            "CInfinity representation may contain only the declared parameter"
        )

    n = representation.rows
    identity = eye(n)
    d0 = _simplify_matrix(representation.subs(parameter, 0))
    scalar = simplify(d0.trace() / n)
    if scalar == 0 or scalar.is_zero is True:
        raise ValueError("CInfinity D(0) must be invertible")
    if not _zero_matrix(d0 - scalar * identity):
        raise ValueError(
            "CInfinity D(0) must be a scalar multiple of the identity"
        )

    real_parameter = {conjugate(parameter): parameter}
    unitary_residual = (
        representation.H * representation - identity
    ).xreplace(real_parameter)
    if not _zero_matrix(unitary_residual):
        raise ValueError("CInfinity D(theta) must be unitary for real theta")

    derivative = representation.diff(parameter)
    raw_generator = _simplify_matrix(
        I * derivative.subs(parameter, 0) * d0.inv()
    )
    logarithmic_generator = _simplify_matrix(
        I * derivative * representation.inv()
    )
    if not _zero_matrix(logarithmic_generator - raw_generator):
        raise ValueError(
            "CInfinity D(theta) must have a theta-independent generator"
        )

    generator = _simplify_matrix(
        raw_generator - simplify(raw_generator.trace() / n) * identity
    )
    if parameter in generator.free_symbols or not _zero_matrix(
        generator - generator.H
    ):
        raise ValueError(
            "the CInfinity generator must be parameter-independent and Hermitian"
        )
    return {
        "raw_generator": raw_generator,
        "generator": generator,
    }


def continuous_invariant_basis_data(generator: MatrixBase, data: dict) -> dict:
    """Restrict an Hermitian operator basis to ``[generator, M] == 0``."""

    generator = Matrix(generator)
    n = data["dimension"]
    if generator.shape != (n, n):
        raise ValueError("CInfinity generator and operator basis dimensions differ")
    images = [
        _simplify_matrix(I * (generator * matrix - matrix * generator))
        for matrix in data["basis"]
    ]
    image_columns = Matrix(
        [_row_flatten(matrix) for matrix in images]
    ).T
    action = _simplify_matrix(data["dual"] * image_columns)
    kernel = action.nullspace()
    if not kernel:
        raise ValueError("the CInfinity commutant could not be constructed")
    reduced_basis = []
    for coordinates in kernel:
        matrix = Matrix.zeros(n)
        for coefficient, basis_matrix in zip(coordinates, data["basis"]):
            matrix += coefficient * basis_matrix
        reduced_basis.append(_simplify_matrix(matrix))
    reduced_data = basis_data(n, reduced_basis)
    return {
        "basis_data": reduced_data,
        "generator_action": action,
        "kernel": kernel,
    }


def continuous_normalizer(
    representation: MatrixBase,
    generator: MatrixBase,
    antiunitary: bool = False,
) -> bool:
    """Whether a discrete operation sends the generator to ``+Q`` or ``-Q``."""

    representation = Matrix(representation)
    generator = Matrix(generator)
    if representation.shape != generator.shape:
        return False
    try:
        inverse = representation.inv()
    except Exception:
        return False
    transformed = (
        -representation * generator.conjugate() * inverse
        if antiunitary
        else representation * generator * inverse
    )
    return _zero_matrix(transformed - generator) or _zero_matrix(
        transformed + generator
    )
