"""Static band plotting for symbolic MagneticKP Hamiltonians."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Optional

from sympy import Matrix, sympify
from sympy.matrices import MatrixBase

from .transforms import K_SYMBOLS


class BandPlotError(ValueError):
    """Raised when a Hamiltonian or momentum path cannot be plotted."""


def _plotting_dependencies():
    try:
        import numpy as np
        import matplotlib.pyplot as plt
    except ModuleNotFoundError as error:
        if error.name not in {"numpy", "matplotlib", "matplotlib.pyplot"}:
            raise
        raise ImportError(
            "bandplot requires NumPy and Matplotlib; reinstall the MagneticKP "
            "wheel so its required dependencies are installed"
        ) from error
    return np, plt


def _coordinates(value, *, location: str):
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise BandPlotError(f"{location} coordinates must be a three-item sequence")
    if len(value) != 3:
        raise BandPlotError(f"{location} coordinates must contain exactly three values")
    try:
        return tuple(float(sympify(component)) for component in value)
    except (TypeError, ValueError) as error:
        raise BandPlotError(f"{location} coordinates must be real numbers") from error


def _endpoint(value, *, location: str):
    if not isinstance(value, Mapping):
        raise BandPlotError(f"{location} must contain coordinates and a label")
    if "coordinates" not in value or "label" not in value:
        raise BandPlotError(f"{location} must contain coordinates and a label")
    return _coordinates(value["coordinates"], location=location), str(value["label"])


def _normalize_path(path):
    if isinstance(path, (str, bytes)) or not isinstance(path, Sequence) or not path:
        raise BandPlotError("path must be a non-empty sequence of path segments")

    normalized = []
    for index, segment in enumerate(path, start=1):
        location = f"path segment {index}"
        if isinstance(segment, Mapping):
            if "start" not in segment or "end" not in segment:
                raise BandPlotError(f"{location} must contain start and end points")
            start, start_label = _endpoint(
                segment["start"], location=f"{location} start"
            )
            end, end_label = _endpoint(segment["end"], location=f"{location} end")
        else:
            if (
                isinstance(segment, (str, bytes))
                or not isinstance(segment, Sequence)
                or len(segment) != 3
            ):
                raise BandPlotError(
                    f"{location} must be (start, end, (start_label, end_label))"
                )
            start = _coordinates(segment[0], location=f"{location} start")
            end = _coordinates(segment[1], location=f"{location} end")
            labels = segment[2]
            if (
                isinstance(labels, (str, bytes))
                or not isinstance(labels, Sequence)
                or len(labels) != 2
            ):
                raise BandPlotError(f"{location} labels must contain two values")
            start_label, end_label = map(str, labels)
        normalized.append((start, end, (start_label, end_label)))
    return normalized


def _parameter_rules(hamiltonian: MatrixBase, parameters):
    model_parameters = sorted(
        hamiltonian.free_symbols - set(K_SYMBOLS),
        key=str,
    )
    if not model_parameters:
        if parameters is None:
            return {}
        if isinstance(parameters, Mapping) and not parameters:
            return {}
        if (
            isinstance(parameters, Sequence)
            and not isinstance(parameters, (str, bytes))
            and len(parameters) == 0
        ):
            return {}
        raise BandPlotError("this Hamiltonian has no free model parameters")
    if parameters is None:
        names = ", ".join(map(str, model_parameters))
        raise BandPlotError(f"parameter values are required for: {names}")

    if isinstance(parameters, Mapping):
        by_name = {str(symbol): symbol for symbol in model_parameters}
        rules = {}
        unknown = []
        for key, value in parameters.items():
            symbol = key if key in model_parameters else by_name.get(str(key))
            if symbol is None:
                unknown.append(str(key))
            else:
                rules[symbol] = value
        if unknown:
            raise BandPlotError("unknown model parameters: " + ", ".join(unknown))
        missing = [symbol for symbol in model_parameters if symbol not in rules]
        if missing:
            raise BandPlotError(
                "missing parameter values for: " + ", ".join(map(str, missing))
            )
        return rules

    if isinstance(parameters, (str, bytes)) or not isinstance(parameters, Sequence):
        raise BandPlotError("parameters must be a mapping or a sequence of values")
    if len(parameters) != len(model_parameters):
        names = ", ".join(map(str, model_parameters))
        raise BandPlotError(
            f"expected {len(model_parameters)} parameter values in this order: {names}"
        )
    return dict(zip(model_parameters, parameters))


def _tick_labels(path):
    labels = [path[0][2][0]]
    for index, segment in enumerate(path):
        end_label = segment[2][1]
        if index + 1 < len(path):
            next_start_label = path[index + 1][2][0]
            if next_start_label != end_label:
                end_label = f"{end_label}|{next_start_label}"
        labels.append(end_label)
    return labels


def bandplot(
    hamiltonian,
    path,
    parameters=None,
    *,
    subdivisions: int = 30,
    coordinate_basis: str = "fractional",
    ylim=None,
    yticks=None,
    fontsize: float = 12,
    fontfamily: Optional[str] = None,
    figsize=None,
    ax=None,
    show: bool = True,
):
    """Compute and plot the bands of a symbolic Hamiltonian along ``path``.

    ``path`` accepts both ``high_symmetry_k_path`` output and compact segments
    of the form ``(start, end, (start_label, end_label))``.  Model parameter
    values may be a mapping or a sequence ordered by parameter name.  Path
    coordinates are reciprocal fractional coordinates by default, so they are
    multiplied by ``2*pi`` before substitution.

    Returns the Matplotlib ``(figure, axes)`` pair.
    """

    np, plt = _plotting_dependencies()
    try:
        matrix = Matrix(hamiltonian)
    except (TypeError, ValueError) as error:
        raise BandPlotError("hamiltonian must be a square SymPy-compatible matrix") from error
    if matrix.rows != matrix.cols or matrix.rows == 0:
        raise BandPlotError("hamiltonian must be a non-empty square matrix")
    if not isinstance(subdivisions, int) or isinstance(subdivisions, bool) or subdivisions < 1:
        raise BandPlotError("subdivisions must be a positive integer")
    if coordinate_basis not in {"fractional", "cartesian"}:
        raise BandPlotError("coordinate_basis must be 'fractional' or 'cartesian'")

    normalized_path = _normalize_path(path)
    rules = _parameter_rules(matrix, parameters)
    matrix = matrix.subs(rules)
    unresolved = sorted(matrix.free_symbols - set(K_SYMBOLS), key=str)
    if unresolved:
        raise BandPlotError(
            "parameter values must be numeric; unresolved symbols: "
            + ", ".join(map(str, unresolved))
        )

    momentum_scale = 2 * np.pi if coordinate_basis == "fractional" else 1.0
    band_segments = []
    for segment_index, (start, end, _labels) in enumerate(normalized_path, start=1):
        points = np.linspace(start, end, subdivisions + 1)
        energies = []
        for point_index, point in enumerate(points):
            numeric = matrix.subs(
                dict(zip(K_SYMBOLS, momentum_scale * np.asarray(point)))
            )
            try:
                array = np.asarray(numeric.evalf(), dtype=complex)
            except (TypeError, ValueError) as error:
                raise BandPlotError(
                    f"Hamiltonian is not numeric at path segment {segment_index}, "
                    f"sample {point_index}"
                ) from error
            if not np.all(np.isfinite(array)):
                raise BandPlotError(
                    f"Hamiltonian is not finite at path segment {segment_index}, "
                    f"sample {point_index}"
                )
            if not np.allclose(array, array.conjugate().T, rtol=1e-9, atol=1e-10):
                raise BandPlotError(
                    f"Hamiltonian is not Hermitian at path segment {segment_index}, "
                    f"sample {point_index}"
                )
            energies.append(np.linalg.eigvalsh(array))
        band_segments.append(np.asarray(energies))

    if ax is None:
        figure, axis = plt.subplots(figsize=figsize)
    else:
        axis = ax
        try:
            figure = axis.figure
        except AttributeError as error:
            raise BandPlotError("ax must be a Matplotlib Axes object") from error
        if figsize is not None:
            figure.set_size_inches(*figsize)

    for segment_index, segment in enumerate(band_segments):
        x_values = np.arange(subdivisions + 1) + segment_index * subdivisions
        for band in segment.T:
            axis.plot(x_values, band, color="black", linewidth=1.5)

    boundaries = np.arange(len(normalized_path) + 1) * subdivisions
    axis.set_xlim(boundaries[0], boundaries[-1])
    axis.margins(x=0)
    axis.set_xticks(boundaries, _tick_labels(normalized_path))
    for boundary in boundaries:
        axis.axvline(boundary, color="0.8", linewidth=0.8)
    axis.axhline(0, color="0.88", linewidth=0.8)
    axis.set_ylabel("Energy")
    if ylim is not None:
        if len(ylim) != 2:
            raise BandPlotError("ylim must contain (minimum, maximum)")
        axis.set_ylim(*ylim)
    if yticks is not None:
        if (
            isinstance(yticks, Sequence)
            and len(yticks) == 2
            and isinstance(yticks[1], Sequence)
            and not isinstance(yticks[1], (str, bytes))
        ):
            axis.set_yticks(yticks[0], yticks[1])
        else:
            axis.set_yticks(yticks)
    axis.tick_params(labelsize=fontsize)
    if fontfamily is not None:
        for label in axis.get_xticklabels() + axis.get_yticklabels():
            label.set_fontfamily(fontfamily)

    if show:
        plt.show()
    return figure, axis


__all__ = ["BandPlotError", "bandplot"]
