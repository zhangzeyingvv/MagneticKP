"""Invariant-space solvers used by the Python MagneticKP implementation."""

from __future__ import annotations

from typing import Dict, List, Sequence, Tuple

from sympy import (
    Float,
    Matrix,
    SparseMatrix,
    eye,
    kronecker_product,
    simplify,
    srepr,
    sympify,
    zeros,
)
from sympy.matrices import MatrixBase


def _representation(factor: dict) -> MatrixBase:
    return factor["eta"] * kronecker_product(factor["F"], factor["J"])


def iterative_invariant_basis(factors: Sequence[dict]) -> Matrix:
    """Return invariant vectors as columns using generator constraints."""

    if not factors:
        return Matrix([])
    ambient = factors[0]["F"].rows * factors[0]["J"].rows
    candidate = eye(ambient)
    identity = eye(ambient)
    for factor in factors:
        kernel = ((_representation(factor) - identity) * candidate).nullspace()
        if not kernel:
            return zeros(ambient, 0)
        candidate = candidate * Matrix.hstack(*kernel)
    return candidate


def stacked_invariant_basis(factors: Sequence[dict]) -> Matrix:
    """Return the common invariant space from one stacked constraint matrix."""

    if not factors:
        return Matrix([])
    ambient = factors[0]["F"].rows * factors[0]["J"].rows
    identity = eye(ambient)
    constraints = Matrix.vstack(
        *(_representation(factor) - identity for factor in factors)
    )
    columns = constraints.nullspace()
    return Matrix.hstack(*columns) if columns else zeros(ambient, 0)


def _contains_inexact(matrix: MatrixBase) -> bool:
    return any(value.has(Float) for value in matrix)


def _factor_inexact(factor: dict) -> bool:
    return (
        _contains_inexact(factor["F"])
        or _contains_inexact(factor["J"])
        or sympify(factor["eta"]).has(Float)
    )


def _canonical_value(value):
    return simplify(value)


def _matrix_key(matrix: MatrixBase) -> tuple:
    return matrix.shape, tuple(srepr(_canonical_value(value)) for value in matrix)


def _factor_key(factor: dict) -> tuple:
    return (
        _matrix_key(factor["F"]),
        _matrix_key(factor["J"]),
        srepr(_canonical_value(factor["eta"])),
    )


def _identity_factor(factor: dict) -> dict:
    return {
        "F": SparseMatrix.eye(factor["F"].rows),
        "J": SparseMatrix.eye(factor["J"].rows),
        "eta": simplify(1),
    }


def _factor_product(left: dict, right: dict) -> dict:
    return {
        "F": (left["F"] * right["F"]).applyfunc(simplify),
        "J": (left["J"] * right["J"]).applyfunc(simplify),
        "eta": simplify(left["eta"] * right["eta"]),
    }


def _factor_closure(generators: Sequence[dict], max_group_order: int):
    if not generators or any(_factor_inexact(factor) for factor in generators):
        return None
    identity = _identity_factor(generators[0])
    queue = [identity]
    seen = {_factor_key(identity): identity}
    head = 0
    while head < len(queue):
        current = queue[head]
        head += 1
        for generator in generators:
            new = _factor_product(current, generator)
            key = _factor_key(new)
            if key in seen:
                continue
            if len(queue) >= max_group_order:
                return None
            seen[key] = new
            queue.append(new)
    return list(seen.values())


def _zero_matrix(matrix: MatrixBase) -> bool:
    return all(simplify(value) == 0 for value in matrix)


def _projector(elements: Sequence[dict]) -> Matrix:
    ambient = elements[0]["F"].rows * elements[0]["J"].rows
    total = zeros(ambient)
    for factor in elements:
        total += _representation(factor)
    return total.applyfunc(lambda value: simplify(value / len(elements)))


def reynolds_invariant_basis(
    factors: Sequence[dict],
    group_input: str = "generators",
    max_group_order: int = 512,
):
    """Return ``(column_basis, diagnostics)`` for a finite factor group.

    Generator input is closed in paired ``{F,J,eta}`` form, so the large
    Kronecker-product group is never stored.  Inexact or unclosed input falls
    back to the iterative solver.
    """

    diagnostics = {
        "solver": "reynolds",
        "fallback": False,
        "element_closure_completed": False,
        "group_order": None,
        "reason": None,
    }
    if not factors:
        return Matrix([]), diagnostics
    if group_input.lower() == "elements":
        elements = None if any(_factor_inexact(factor) for factor in factors) else list(factors)
    else:
        elements = _factor_closure(factors, max_group_order)
    if elements is None:
        diagnostics.update(
            fallback=True,
            reason="exact finite factor closure could not be established",
        )
        return iterative_invariant_basis(factors), diagnostics

    projector = _projector(elements)
    valid = _zero_matrix(projector * projector - projector)
    if valid:
        valid = all(
            _zero_matrix(_representation(factor) * projector - projector)
            for factor in factors
        )
    if not valid and group_input.lower() == "elements":
        completed_elements = _factor_closure(factors, max_group_order)
        if completed_elements is not None:
            completed_projector = _projector(completed_elements)
            completed_valid = _zero_matrix(
                completed_projector * completed_projector
                - completed_projector
            )
            if completed_valid:
                completed_valid = all(
                    _zero_matrix(
                        _representation(factor) * completed_projector
                        - completed_projector
                    )
                    for factor in factors
                )
            if completed_valid:
                elements = completed_elements
                projector = completed_projector
                valid = True
                diagnostics["element_closure_completed"] = True
    if not valid:
        diagnostics.update(
            fallback=True,
            group_order=len(elements),
            reason="Reynolds idempotence or generator absorption failed",
        )
        return iterative_invariant_basis(factors), diagnostics

    columns = projector.columnspace()
    basis = Matrix.hstack(*columns) if columns else zeros(projector.rows, 0)
    diagnostics["group_order"] = len(elements)
    return basis, diagnostics
