from __future__ import annotations

import importlib.util
import pathlib
import subprocess
import sys

from sympy import Matrix

from conftest import REPOSITORY_ROOT


def git_file(revision_path: str) -> bytes:
    return subprocess.run(
        ["git", "show", revision_path],
        cwd=REPOSITORY_ROOT,
        check=True,
        stdout=subprocess.PIPE,
    ).stdout


def load_legacy_package(tmp_path: pathlib.Path):
    package = tmp_path / "legacy_magnetickp"
    package.mkdir()
    for name in ("__init__.py", "magkp.py", "fd.pkl"):
        (package / name).write_bytes(
            git_file(
                f"a319488:MagneticKP-python/magnetickp/{name}"
            )
        )
    spec = importlib.util.spec_from_file_location(
        "legacy_magnetickp",
        package / "__init__.py",
        submodule_search_locations=[str(package)],
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules["legacy_magnetickp"] = module
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def coefficient_subspace(hamiltonian):
    k_symbols = {
        symbol
        for symbol in hamiltonian.free_symbols
        if str(symbol) in {"kx", "ky", "kz"}
    }
    parameters = sorted(
        hamiltonian.free_symbols - k_symbols,
        key=str,
    )
    if not parameters:
        return Matrix.zeros(0, hamiltonian.rows * hamiltonian.cols)
    rows = []
    for parameter in parameters:
        coefficient = hamiltonian.applyfunc(
            lambda entry: entry.coeff(parameter)
        )
        rows.append(Matrix(list(coefficient)).T)
    return Matrix.vstack(*rows)


def same_subspace(left, right):
    left = coefficient_subspace(left)
    right = coefficient_subspace(right)
    return (
        left.rank()
        == right.rank()
        == Matrix.vstack(left, right).rank()
    )


def test_python_kpham_old_parity(tmp_path):
    import magnetickp

    legacy = load_legacy_package(tmp_path)
    input_data = {
        "Unitary": {
            "C2z": (
                Matrix([[1, 0], [0, -1]]),
                Matrix([[-1, 0, 0], [0, -1, 0], [0, 0, 1]]),
            )
        },
        "Anitunitary": {
            "T": (
                Matrix.eye(2),
                -Matrix.eye(3),
            )
        },
    }
    old_result = legacy.kpHam(input_data).getkpHam(2)
    new_result = magnetickp.kpHam(input_data).getkpHam(2)
    assert old_result["NumberOfParameters"] == new_result[
        "NumberOfParameters"
    ]
    assert old_result["dim"] == new_result["dim"] == 2
    assert same_subspace(old_result["ham"], new_result["ham"])
