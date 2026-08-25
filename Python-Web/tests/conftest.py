from __future__ import annotations

import importlib
import inspect
import pathlib
import sys
from typing import Any, Callable, Iterable

import pytest


PYTHON_ROOT = pathlib.Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = PYTHON_ROOT.parent
if str(PYTHON_ROOT) not in sys.path:
    sys.path.insert(0, str(PYTHON_ROOT))


MODULE_NAMES = (
    "magnetickp",
    "magnetickp.magkp",
    "magnetickp.basis",
    "magnetickp.transforms",
    "magnetickp.solvers",
)


def loaded_modules() -> Iterable[Any]:
    for name in MODULE_NAMES:
        try:
            yield importlib.import_module(name)
        except ModuleNotFoundError:
            continue


def resolve_callable(*names: str) -> Callable[..., Any]:
    for module in loaded_modules():
        for name in names:
            value = getattr(module, name, None)
            if callable(value):
                return value
    pytest.fail(f"Required MagneticKP callable not found: {names}")


def call_variants(
    function: Callable[..., Any],
    variants: Iterable[tuple[tuple[Any, ...], dict[str, Any]]],
) -> Any:
    errors: list[str] = []
    for args, kwargs in variants:
        try:
            value = function(*args, **kwargs)
        except (TypeError, ValueError) as exc:
            errors.append(f"{args!r} {kwargs!r}: {exc}")
            continue
        return value
    pytest.fail(
        f"No supported signature for {function!r}. Attempts:\n"
        + "\n".join(errors)
    )


@pytest.fixture(scope="session")
def operator_basis_callable() -> Callable[..., Any]:
    return resolve_callable(
        "basis_data",
        "operator_basis_data",
        "OperatorBasisData",
        "_operator_basis_data",
    )


@pytest.fixture(scope="session")
def momentum_transform_callable() -> Callable[..., Any]:
    return resolve_callable(
        "momentum_transform_matrix",
        "MomentumTransformMatrix",
        "_momentum_transform_matrix",
    )


@pytest.fixture(scope="session")
def adjoint_callable() -> Callable[..., Any]:
    return resolve_callable(
        "adjoint_representation_matrix",
        "AdjointRepresentationMatrix",
        "_adjoint_representation_matrix",
    )


@pytest.fixture(scope="session")
def reynolds_callable() -> Callable[..., Any]:
    return resolve_callable(
        "reynolds_invariant_basis",
        "ReynoldsInvariantBasis",
        "_reynolds_invariant_basis",
    )


def basis_list(data: Any) -> list[Any]:
    if isinstance(data, dict):
        for key in ("Basis", "basis", "matrices"):
            if key in data:
                return list(data[key])
    if isinstance(data, (list, tuple)):
        return list(data)
    pytest.fail(f"Cannot extract basis list from {type(data)!r}")


def gram_matrix(data: Any) -> Any | None:
    if isinstance(data, dict):
        for key in ("Gram", "gram"):
            if key in data:
                return data[key]
    return None


def call_operator_basis(function: Callable[..., Any], n: int) -> Any:
    return call_variants(function, (((n,), {}),))


def call_momentum_transform(
    function: Callable[..., Any],
    order: int,
    image: Any,
    variables: Any,
    basis: Any,
) -> Any:
    return call_variants(
        function,
        (
            (
                (order, image),
                {"variables": variables, "klist": basis},
            ),
            ((order, image), {}),
            ((order, image, variables), {}),
            ((order, basis, image, variables), {}),
            ((basis, image, variables), {}),
        ),
    )


def call_adjoint(
    function: Callable[..., Any],
    representation: Any,
    basis_data: Any,
    antiunitary: bool,
) -> Any:
    return call_variants(
        function,
        (
            ((representation, basis_data, antiunitary), {}),
            (
                (representation, basis_data),
                {"antiunitary": antiunitary},
            ),
            (
                (representation, basis_data),
                {"anti_unitary": antiunitary},
            ),
            (
                (representation, basis_list(basis_data), antiunitary),
                {},
            ),
        ),
    )


def call_reynolds(
    function: Callable[..., Any],
    representations: list[Any],
    characters: list[Any] | None = None,
) -> Any:
    from sympy import Matrix

    characters = (
        [1] * len(representations)
        if characters is None
        else characters
    )
    factors = [
        {
            "F": Matrix(representation),
            "J": Matrix([[1]]),
            "eta": character,
        }
        for representation, character in zip(
            representations, characters
        )
    ]
    variants = [
        ((factors,), {"group_input": "elements"}),
        ((representations, characters), {}),
        ((representations,), {"characters": characters}),
        ((representations,), {}),
    ]
    return call_variants(function, variants)
