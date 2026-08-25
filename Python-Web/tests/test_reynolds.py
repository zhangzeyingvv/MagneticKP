from __future__ import annotations

from sympy import Matrix, Rational, sqrt
from sympy.matrices import MatrixBase

from conftest import call_reynolds


def image_basis(projector):
    return Matrix(projector).columnspace()


def stacked_invariant_basis(representations):
    dimension = representations[0].rows
    identity = Matrix.eye(dimension)
    stacked = Matrix.vstack(
        *(Matrix(rep) - identity for rep in representations)
    )
    return stacked.nullspace()


def same_column_space(left, right, dimension):
    left_matrix = (
        Matrix.hstack(*left) if left else Matrix.zeros(dimension, 0)
    )
    right_matrix = (
        Matrix.hstack(*right) if right else Matrix.zeros(dimension, 0)
    )
    return (
        left_matrix.rank()
        == right_matrix.rank()
        == Matrix.hstack(left_matrix, right_matrix).rank()
    )


def normalize_basis_result(value, dimension):
    if (
        isinstance(value, tuple)
        and len(value) == 2
        and isinstance(value[1], dict)
    ):
        value = value[0]
    if isinstance(value, MatrixBase):
        if value.rows == dimension:
            return value.columnspace()
        if value.cols == dimension:
            return value.T.columnspace()
    if isinstance(value, (list, tuple)):
        vectors = [Matrix(vector) for vector in value]
        return [
            vector if vector.rows == dimension else vector.T
            for vector in vectors
        ]
    raise TypeError(f"Unsupported Reynolds basis result: {type(value)!r}")


def test_reynolds_projector_properties(reynolds_callable):
    c2 = [Matrix.eye(2), Matrix.diag(1, -1)]
    c3_generator = Matrix([[0, -1], [1, -1]])
    c3 = [Matrix.eye(2), c3_generator, c3_generator**2]
    for representations in (c2, c3):
        projector = sum(representations, Matrix.zeros(2)) / len(
            representations
        )
        assert (projector * projector).equals(projector)
        assert all((rep * projector).equals(projector)
                   for rep in representations)
        assert projector.rank() == projector.trace()
        result = normalize_basis_result(
            call_reynolds(reynolds_callable, representations),
            2,
        )
        assert same_column_space(
            result,
            image_basis(projector),
            2,
        )
        assert same_column_space(
            result,
            stacked_invariant_basis(representations),
            2,
        )


def test_twisted_reynolds(reynolds_callable):
    representations = [Matrix.eye(2), Matrix.diag(1, -1)]
    characters = [1, -1]
    projector = sum(
        (
            (
                character.conjugate() * representation
                if hasattr(character, "conjugate")
                else character * representation
            )
            for character, representation in zip(
                characters, representations
            )
        ),
        Matrix.zeros(2),
    ) / 2
    result = normalize_basis_result(
        call_reynolds(
            reynolds_callable,
            representations,
            characters,
        ),
        2,
    )
    assert same_column_space(result, image_basis(projector), 2)
