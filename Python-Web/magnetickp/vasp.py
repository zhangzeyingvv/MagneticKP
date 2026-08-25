"""Readers for VASP band-structure output."""

from __future__ import annotations

from numbers import Integral
from pathlib import Path


class VaspEigenvalError(ValueError):
    """Raised when a VASP EIGENVAL file or selection is invalid."""


def _positive_index(value, *, name: str) -> int:
    if not isinstance(value, Integral) or isinstance(value, bool) or value < 1:
        raise VaspEigenvalError(f"{name} must be a positive one-based integer")
    return int(value)


def read_vasp_eigenval(
    path,
    efermi: float = 0.0,
    spin: int = 1,
    start_band: int = 1,
    end_band=None,
):
    """Read selected bands from a VASP ``EIGENVAL`` file.

    ``spin``, ``start_band``, and ``end_band`` use the one-based numbering from
    VASP and the Mathematica ``vaspEig`` interface.  The returned value is a
    list of ``((kx, ky, kz), [energy, ...])`` records with ``efermi``
    subtracted from every selected energy.
    """

    spin = _positive_index(spin, name="spin")
    start_band = _positive_index(start_band, name="start_band")
    if end_band is not None:
        end_band = _positive_index(end_band, name="end_band")
        if end_band < start_band:
            raise VaspEigenvalError("end_band must not be smaller than start_band")
    try:
        fermi_level = float(efermi)
    except (TypeError, ValueError) as error:
        raise VaspEigenvalError("efermi must be a real number") from error

    source = Path(path)
    try:
        lines = source.read_text(encoding="utf-8").splitlines()
    except OSError as error:
        raise VaspEigenvalError(f"cannot read EIGENVAL file: {source}") from error
    if len(lines) < 6:
        raise VaspEigenvalError("EIGENVAL header must contain at least six lines")

    try:
        spin_channels = int(lines[0].split()[-1])
        _electron_count, kpoint_count, band_count = map(
            int, lines[5].split()[-3:]
        )
    except (IndexError, ValueError) as error:
        raise VaspEigenvalError("invalid EIGENVAL header") from error
    if spin_channels not in (1, 2):
        raise VaspEigenvalError("EIGENVAL header has an invalid spin count")
    if spin > spin_channels:
        raise VaspEigenvalError(
            f"spin={spin} is unavailable; this EIGENVAL has {spin_channels} spin channel(s)"
        )
    if kpoint_count < 1 or band_count < 1:
        raise VaspEigenvalError("EIGENVAL must contain at least one k point and band")
    if end_band is None:
        end_band = band_count
    if start_band > band_count or end_band > band_count:
        raise VaspEigenvalError(
            f"band selection {start_band}-{end_band} exceeds the {band_count} available bands"
        )

    cursor = 6
    records = []
    for kpoint_index in range(1, kpoint_count + 1):
        while cursor < len(lines) and not lines[cursor].strip():
            cursor += 1
        if cursor >= len(lines):
            raise VaspEigenvalError(
                f"missing coordinates for k point {kpoint_index}"
            )
        try:
            coordinates = tuple(
                float(value.replace("D", "E").replace("d", "e"))
                for value in lines[cursor].split()[:3]
            )
        except ValueError as error:
            raise VaspEigenvalError(
                f"invalid coordinates for k point {kpoint_index}"
            ) from error
        if len(coordinates) != 3:
            raise VaspEigenvalError(
                f"invalid coordinates for k point {kpoint_index}"
            )
        cursor += 1

        selected_energies = []
        for band_index in range(1, band_count + 1):
            if cursor >= len(lines):
                raise VaspEigenvalError(
                    f"missing band {band_index} at k point {kpoint_index}"
                )
            row = lines[cursor].split()
            cursor += 1
            try:
                file_band_index = int(row[0])
                energy = float(
                    row[spin].replace("D", "E").replace("d", "e")
                )
            except (IndexError, ValueError) as error:
                raise VaspEigenvalError(
                    f"invalid band {band_index} at k point {kpoint_index}"
                ) from error
            if file_band_index != band_index:
                raise VaspEigenvalError(
                    f"expected band {band_index} at k point {kpoint_index}, "
                    f"found band {file_band_index}"
                )
            if start_band <= band_index <= end_band:
                selected_energies.append(energy - fermi_level)
        records.append((coordinates, selected_energies))
    return records


__all__ = ["VaspEigenvalError", "read_vasp_eigenval"]
