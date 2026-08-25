"""Numerical fitting of symbolic k.p Hamiltonians to VASP bands."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from numbers import Integral

from sympy import Matrix, lambdify

from .transforms import K_SYMBOLS
from .vasp import read_vasp_eigenval


class KpFittingError(ValueError):
    """Raised when a k.p-to-VASP fit cannot be constructed."""


def _fitting_dependencies():
    try:
        import numpy as np
        from scipy.optimize import least_squares
    except ModuleNotFoundError as error:
        if error.name not in {"numpy", "scipy", "scipy.optimize"}:
            raise
        raise ImportError(
            "fittingKP requires NumPy and SciPy; reinstall the MagneticKP "
            "wheel so its required dependencies are installed"
        ) from error
    return np, least_squares


def _initial_values(parameters, initial_parameters):
    if initial_parameters is None:
        return [0.0] * len(parameters)
    if isinstance(initial_parameters, Mapping):
        by_name = {str(symbol): symbol for symbol in parameters}
        values = {}
        unknown = []
        for key, value in initial_parameters.items():
            symbol = key if key in parameters else by_name.get(str(key))
            if symbol is None:
                unknown.append(str(key))
            else:
                values[symbol] = value
        if unknown:
            raise KpFittingError("unknown initial parameters: " + ", ".join(unknown))
        missing = [symbol for symbol in parameters if symbol not in values]
        if missing:
            raise KpFittingError(
                "missing initial values for: " + ", ".join(map(str, missing))
            )
        ordered = [values[symbol] for symbol in parameters]
    elif isinstance(initial_parameters, Sequence) and not isinstance(
        initial_parameters, (str, bytes)
    ):
        if len(initial_parameters) != len(parameters):
            raise KpFittingError(
                f"expected {len(parameters)} initial values in this order: "
                + ", ".join(map(str, parameters))
            )
        ordered = list(initial_parameters)
    else:
        raise KpFittingError(
            "initial_parameters must be a mapping or a sequence of values"
        )
    try:
        return [float(value) for value in ordered]
    except (TypeError, ValueError) as error:
        raise KpFittingError("initial parameter values must be real numbers") from error


def _one_based_indices(kpoint_indices, record_count):
    if isinstance(kpoint_indices, (str, bytes)):
        raise KpFittingError("kpoint_indices must be a non-empty sequence")
    try:
        indices = tuple(kpoint_indices)
    except TypeError as error:
        raise KpFittingError("kpoint_indices must be a non-empty sequence") from error
    if not indices:
        raise KpFittingError("kpoint_indices must be a non-empty sequence")
    normalized = []
    for index in indices:
        if not isinstance(index, Integral) or isinstance(index, bool) or index < 1:
            raise KpFittingError("kpoint_indices use positive one-based numbering")
        if index > record_count:
            raise KpFittingError(
                f"k-point index {index} exceeds the {record_count} EIGENVAL records"
            )
        normalized.append(int(index))
    return tuple(normalized)


def _center_coordinates(kpoint_center):
    if kpoint_center is None:
        return (0.0, 0.0, 0.0)
    if isinstance(kpoint_center, (str, bytes)):
        raise KpFittingError("kpoint_center must contain three coordinates")
    try:
        center = tuple(float(value) for value in kpoint_center)
    except (TypeError, ValueError) as error:
        raise KpFittingError("kpoint_center must contain three real coordinates") from error
    if len(center) != 3:
        raise KpFittingError("kpoint_center must contain three coordinates")
    return center


def fittingKP(
    hamiltonian,
    eigenval_path,
    kpoint_indices,
    initial_parameters=None,
    *,
    efermi: float = 0.0,
    spin: int = 1,
    start_band: int = 1,
    end_band=None,
    kpoint_center=None,
    method: str = "lm",
    ftol: float = 1e-8,
    xtol: float = 1e-8,
    gtol: float = 1e-8,
    max_nfev=None,
):
    """Fit free Hamiltonian parameters to selected VASP band records.

    ``kpoint_indices`` uses one-based EIGENVAL record numbers.  For a local
    Hamiltonian expanded at K, pass its absolute fractional coordinates as
    ``kpoint_center``; EIGENVAL coordinates are shifted by this center before
    substitution.  No ``2*pi`` factor is applied.  ``initial_parameters``
    defaults to zero for every free parameter.
    """

    np, least_squares = _fitting_dependencies()
    try:
        matrix = Matrix(hamiltonian)
    except (TypeError, ValueError) as error:
        raise KpFittingError("hamiltonian must be a square SymPy-compatible matrix") from error
    if matrix.rows != matrix.cols or matrix.rows == 0:
        raise KpFittingError("hamiltonian must be a non-empty square matrix")

    parameters = sorted(matrix.free_symbols - set(K_SYMBOLS), key=str)
    if not parameters:
        raise KpFittingError("hamiltonian has no free model parameters to fit")
    initial = np.asarray(
        _initial_values(parameters, initial_parameters),
        dtype=float,
    )

    if end_band is None:
        end_band = start_band + matrix.rows - 1
    records = read_vasp_eigenval(
        eigenval_path,
        efermi=efermi,
        spin=spin,
        start_band=start_band,
        end_band=end_band,
    )
    indices = _one_based_indices(kpoint_indices, len(records))
    selected = [records[index - 1] for index in indices]
    center = _center_coordinates(kpoint_center)
    if any(len(energies) != matrix.rows for _point, energies in selected):
        raise KpFittingError(
            "the number of selected VASP bands must equal the Hamiltonian dimension"
        )

    numeric_hamiltonian = lambdify(
        (*K_SYMBOLS, *parameters),
        matrix,
        modules="numpy",
    )
    target = np.asarray(
        [energy for _point, energies in selected for energy in energies],
        dtype=float,
    )
    if method == "lm" and target.size < len(parameters):
        raise KpFittingError(
            "method='lm' requires at least as many fitted energies as parameters; "
            "select more k points or use method='trf'"
        )

    def predictions(values):
        result = []
        for record_index, (point, _energies) in zip(indices, selected):
            local_point = tuple(
                coordinate - center_coordinate
                for coordinate, center_coordinate in zip(point, center)
            )
            try:
                numeric = np.asarray(
                    numeric_hamiltonian(*local_point, *values),
                    dtype=complex,
                )
            except (TypeError, ValueError) as error:
                raise KpFittingError(
                    f"Hamiltonian is not numeric at EIGENVAL record {record_index}"
                ) from error
            if numeric.shape != (matrix.rows, matrix.cols):
                raise KpFittingError("numeric Hamiltonian has an unexpected shape")
            if not np.all(np.isfinite(numeric)):
                raise KpFittingError(
                    f"Hamiltonian is not finite at EIGENVAL record {record_index}"
                )
            if not np.allclose(numeric, numeric.conjugate().T, rtol=1e-9, atol=1e-10):
                raise KpFittingError(
                    f"Hamiltonian is not Hermitian at EIGENVAL record {record_index}"
                )
            result.extend(np.linalg.eigvalsh(numeric))
        return np.asarray(result, dtype=float)

    initial_residual = predictions(initial) - target
    result = least_squares(
        lambda values: predictions(values) - target,
        initial,
        method=method,
        ftol=ftol,
        xtol=xtol,
        gtol=gtol,
        max_nfev=max_nfev,
    )
    fitted_residual = predictions(result.x) - target
    return {
        "fitted_parameters": {
            symbol: float(value) for symbol, value in zip(parameters, result.x)
        },
        "initial_parameters": {
            symbol: float(value) for symbol, value in zip(parameters, initial)
        },
        "parameter_order": tuple(parameters),
        "initial_least_squares": float(initial_residual @ initial_residual),
        "fitted_least_squares": float(fitted_residual @ fitted_residual),
        "converged": bool(result.success),
        "iterations": int(result.nfev),
        "message": str(result.message),
        "kpoint_indices": indices,
        "kpoint_center": center,
        "band_range": (start_band, end_band),
    }


fit_kp_to_vasp = fittingKP


__all__ = ["KpFittingError", "fit_kp_to_vasp", "fittingKP"]
