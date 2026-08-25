"""Symmetry-derived local momentum paths for MagneticKP models."""

from __future__ import annotations

from collections.abc import Mapping, Sequence

from sympy import Float, Matrix, eye, expand, simplify, srepr, sstr, sympify, together
from sympy.matrices import MatrixBase

from .transforms import K_SYMBOLS


class HighSymmetryPathError(ValueError):
    """Raised when a finite three-dimensional momentum group cannot yield a path."""


def _canonical_value(value):
    return simplify(together(value))


def _canonical_matrix(value) -> Matrix:
    return Matrix(value).applyfunc(_canonical_value)


def _matrix_key(matrix: MatrixBase) -> tuple:
    canonical = _canonical_matrix(matrix)
    return canonical.shape, tuple(srepr(value) for value in canonical)


def _zero_value(value) -> bool:
    canonical = _canonical_value(value)
    return canonical == 0 or canonical.is_zero is True


def _zero_vector(vector) -> bool:
    return all(_zero_value(value) for value in vector)


def _zero_matrix(matrix: MatrixBase) -> bool:
    return _zero_vector(matrix)


def _element_key(element) -> tuple:
    spatial, parity = element
    return _matrix_key(spatial), int(parity) % 2


def _canonical_element(element):
    spatial, parity = element
    return _canonical_matrix(spatial), int(parity) % 2


def _element_product(left, right):
    return _canonical_element((left[0] * right[0], left[1] + right[1]))


def _element_action(element) -> Matrix:
    spatial, parity = element
    return spatial if parity == 0 else -spatial


def _operation_values(category):
    if category is None:
        return []
    if isinstance(category, Mapping):
        return list(category.values())
    if isinstance(category, Sequence) and not isinstance(category, (str, bytes)):
        return list(category)
    raise HighSymmetryPathError("symmetry categories must be mappings or sequences")


def _momentum_action_matrix(operation) -> Matrix:
    if not isinstance(operation, Sequence) or isinstance(operation, (str, bytes)):
        raise HighSymmetryPathError("each symmetry operation must contain (D, R)")
    if len(operation) < 2:
        raise HighSymmetryPathError("each symmetry operation must contain (D, R)")
    image = operation[1]
    if isinstance(image, MatrixBase) and image.shape == (3, 3):
        matrix = _canonical_matrix(image)
    else:
        try:
            image_vector = Matrix(list(image))
        except (TypeError, ValueError) as error:
            raise HighSymmetryPathError(
                "momentum images must be three-vectors or 3 x 3 matrices"
            ) from error
        if image_vector.shape != (3, 1):
            raise HighSymmetryPathError(
                "high-symmetry paths require three-dimensional momentum images"
            )
        matrix = Matrix(
            3,
            3,
            lambda row, column: expand(image_vector[row]).coeff(K_SYMBOLS[column]),
        )
        residual = image_vector - matrix * Matrix(K_SYMBOLS)
        if not _zero_vector(residual):
            raise HighSymmetryPathError(
                "momentum images must be linear and homogeneous"
            )
        matrix = _canonical_matrix(matrix)
    if _zero_value(matrix.det()):
        raise HighSymmetryPathError("momentum transformations must be invertible")
    if any(value.has(Float) for value in matrix):
        raise HighSymmetryPathError(
            "high-symmetry paths require exact momentum transformations"
        )
    return matrix


def _symmetry_generators(inputdict: Mapping):
    if not isinstance(inputdict, Mapping):
        raise HighSymmetryPathError("inputdict must be a mapping")
    unitary = _operation_values(inputdict.get("Unitary"))
    antiunitary = _operation_values(
        inputdict.get("Antiunitary", inputdict.get("Anitunitary"))
    )
    generators = [(_momentum_action_matrix(value), 0) for value in unitary]
    generators.extend(
        (-_momentum_action_matrix(value), 1) for value in antiunitary
    )
    canonical = []
    seen = set()
    for generator in generators:
        value = _canonical_element(generator)
        key = _element_key(value)
        if key not in seen:
            seen.add(key)
            canonical.append(value)
    if not canonical:
        raise HighSymmetryPathError(
            "at least one finite three-dimensional unitary or antiunitary operation is required"
        )
    return canonical


def _group_closure(generators, max_group_order: int):
    identity = (eye(3), 0)
    elements = [identity]
    index_by_key = {_element_key(identity): 0}
    head = 0
    while head < len(elements):
        current = elements[head]
        head += 1
        for generator in generators:
            product = _element_product(current, generator)
            key = _element_key(product)
            if key in index_by_key:
                continue
            if len(elements) >= max_group_order:
                raise HighSymmetryPathError(
                    f"little co-group exceeds max_group_order={max_group_order}"
                )
            index_by_key[key] = len(elements)
            elements.append(product)
    return elements, index_by_key


def _product_table(elements, index_by_key):
    table = []
    for left in elements:
        row = []
        for right in elements:
            key = _element_key(_element_product(left, right))
            try:
                row.append(index_by_key[key])
            except KeyError as error:
                raise HighSymmetryPathError(
                    "the completed momentum group is not closed"
                ) from error
        table.append(row)
    return table


def _generated_subgroup(table, generators):
    seen = {0}
    queue = [0]
    head = 0
    generators = tuple(generators)
    while head < len(queue):
        current = queue[head]
        head += 1
        products = [table[current][item] for item in generators]
        products.extend(table[item][current] for item in generators)
        for product in products:
            if product not in seen:
                seen.add(product)
                queue.append(product)
    return tuple(sorted(seen))


def _enumerate_subgroups(table, max_subgroups: int):
    subgroups = [(0,)]
    seen = {subgroups[0]}
    head = 0
    while head < len(subgroups):
        subgroup = subgroups[head]
        head += 1
        for element in range(len(table)):
            if element in subgroup:
                continue
            candidate = _generated_subgroup(table, subgroup + (element,))
            if candidate in seen:
                continue
            if len(subgroups) >= max_subgroups:
                raise HighSymmetryPathError(
                    f"subgroup enumeration exceeds max_subgroups={max_subgroups}"
                )
            seen.add(candidate)
            subgroups.append(candidate)
    return sorted(subgroups, key=lambda item: (len(item), item))


def _fixed_basis(elements, subgroup):
    identity = eye(3)
    constraints = Matrix.vstack(
        *(_element_action(elements[index]) - identity for index in subgroup)
    )
    return constraints.nullspace()


def _normalize_direction(direction):
    exact = [_canonical_value(value) for value in direction]
    magnitudes = [abs(complex(value.evalf())) for value in exact]
    maximum = max(magnitudes)
    if maximum == 0:
        raise HighSymmetryPathError("a fixed momentum direction vanished")
    pivot = magnitudes.index(maximum)
    return tuple(_canonical_value(value / exact[pivot]) for value in exact)


def _collinear(first, second) -> bool:
    return Matrix([list(first), list(second)]).rank() == 1


def _equivalent_line(first, second, elements) -> bool:
    first_vector = Matrix(first)
    return any(
        _collinear(tuple(_element_action(element) * first_vector), second)
        for element in elements
    )


def _distinct_matrices(matrices):
    result = []
    seen = set()
    for matrix in matrices:
        key = _matrix_key(matrix)
        if key not in seen:
            seen.add(key)
            result.append(matrix)
    return result


def _line_data(elements, subgroups):
    identity = eye(3)
    lines = []
    for subgroup in subgroups:
        basis = _fixed_basis(elements, subgroup)
        if len(basis) != 1:
            continue
        direction = _normalize_direction(basis[0])
        stabilizer = tuple(
            index
            for index, element in enumerate(elements)
            if _zero_vector(
                (_element_action(element) - identity) * Matrix(direction)
            )
        )
        actions = [_element_action(elements[index]) for index in stabilizer]
        if any(_collinear(line["direction"], direction) for line in lines):
            continue
        lines.append(
            {
                "direction": direction,
                "subgroup_order": len(stabilizer),
                "momentum_action_order": len(_distinct_matrices(actions)),
                "subgroup_indices": tuple(index + 1 for index in stabilizer),
                "subgroup_matrices": tuple(actions),
                "subgroup_elements": tuple(
                    {
                        "spatial_matrix": elements[index][0],
                        "antiunitary": bool(elements[index][1]),
                        "momentum_action_matrix": _element_action(elements[index]),
                    }
                    for index in stabilizer
                ),
            }
        )
    return sorted(
        lines,
        key=lambda line: (
            -line["subgroup_order"],
            tuple(-float(value.evalf()) for value in line["direction"]),
            tuple(sstr(value) for value in line["direction"]),
        ),
    )


def _validate_numeric_vector(values, name: str):
    if not isinstance(values, Sequence) or isinstance(values, (str, bytes)):
        raise HighSymmetryPathError(f"{name} must be a three-vector")
    if len(values) != 3:
        raise HighSymmetryPathError(f"{name} must be a three-vector")
    result = tuple(sympify(value) for value in values)
    if any(value.free_symbols or value.is_real is False for value in result):
        raise HighSymmetryPathError(f"{name} must contain real numeric values")
    return result


def high_symmetry_k_path(
    inputdict: Mapping,
    *,
    radius=sympify("1/10"),
    center=(0, 0, 0),
    center_label: str | None = None,
    equivalent_directions: str = "representatives",
    include_negative_directions: bool = True,
    return_data: bool = False,
    max_group_order: int = 512,
    max_subgroups: int = 4096,
    coordinate_basis: str | None = None,
):
    """Derive symmetry-inequivalent local momentum lines from ``inputdict``.

    The finite unitary and antiunitary momentum transformations are completed
    as a magnetic little co-group.  Each returned line is fixed by a maximal
    subgroup at a generic point.  Antiunitary parity is retained even when two
    group elements have the same action on momentum.
    """

    radius = sympify(radius)
    if radius.free_symbols or radius.is_real is False or radius.is_positive is not True:
        raise HighSymmetryPathError("radius must be a positive real number")
    center = _validate_numeric_vector(center, "center")
    if center_label is None:
        center_label = inputdict.get("KPointLabel", "k0")
    if not isinstance(center_label, str):
        raise HighSymmetryPathError(
            "center_label must be a string or inputdict must contain a string KPointLabel"
        )
    equivalent_directions = str(equivalent_directions).lower()
    if equivalent_directions not in {"representatives", "all"}:
        raise HighSymmetryPathError(
            "equivalent_directions must be 'representatives' or 'all'"
        )
    if not isinstance(include_negative_directions, bool):
        raise HighSymmetryPathError("include_negative_directions must be bool")
    if not isinstance(return_data, bool):
        raise HighSymmetryPathError("return_data must be bool")
    if not isinstance(max_group_order, int) or max_group_order <= 0:
        raise HighSymmetryPathError("max_group_order must be a positive integer")
    if not isinstance(max_subgroups, int) or max_subgroups <= 0:
        raise HighSymmetryPathError("max_subgroups must be a positive integer")

    if coordinate_basis is None:
        coordinate_system = inputdict.get("CoordinateSystem")
        coordinate_basis = {
            "Cartesian": "cartesian",
            "ReciprocalFractional": "fractional",
        }.get(coordinate_system)

    generators = _symmetry_generators(inputdict)
    elements, index_by_key = _group_closure(generators, max_group_order)
    table = _product_table(elements, index_by_key)
    subgroups = _enumerate_subgroups(table, max_subgroups)
    all_lines = _line_data(elements, subgroups)
    if not all_lines:
        raise HighSymmetryPathError(
            "no subgroup with a one-dimensional fixed momentum space was found"
        )
    if equivalent_directions == "all":
        selected_lines = all_lines
    else:
        selected_lines = []
        for candidate in all_lines:
            if any(
                _equivalent_line(
                    candidate["direction"], representative["direction"], elements
                )
                for representative in selected_lines
            ):
                continue
            selected_lines.append(candidate)

    path = []
    for index, line in enumerate(selected_lines, start=1):
        direction = line["direction"]
        negative = tuple(
            _canonical_value(value - radius * component)
            for value, component in zip(center, direction)
        )
        positive = tuple(
            _canonical_value(value + radius * component)
            for value, component in zip(center, direction)
        )
        if include_negative_directions:
            path.append(
                {
                    "start": {"coordinates": negative, "label": f"L{index}-"},
                    "end": {"coordinates": center, "label": center_label},
                }
            )
            path.append(
                {
                    "start": {"coordinates": center, "label": center_label},
                    "end": {"coordinates": positive, "label": f"L{index}+"},
                }
            )
        else:
            path.append(
                {
                    "start": {"coordinates": center, "label": center_label},
                    "end": {"coordinates": positive, "label": f"L{index}"},
                }
            )

    if not return_data:
        return path
    actions = [_element_action(element) for element in elements]
    return {
        "path": path,
        "directions": tuple(line["direction"] for line in selected_lines),
        "line_data": tuple(selected_lines),
        "little_cogroup_order": len(elements),
        "little_cogroup_action_order": len(_distinct_matrices(actions)),
        "little_cogroup_matrices": tuple(actions),
        "little_cogroup_elements": tuple(
            {
                "spatial_matrix": element[0],
                "antiunitary": bool(element[1]),
                "momentum_action_matrix": _element_action(element),
            }
            for element in elements
        ),
        "enumerated_subgroup_count": len(subgroups),
        "center": center,
        "center_label": center_label,
        "radius": radius,
        "coordinate_basis": coordinate_basis,
    }


def format_band_path(path) -> str:
    """Format a generated path for the Web band-path text field."""

    def point_text(point):
        coordinates = ",".join(sstr(value) for value in point["coordinates"])
        return f'{point["label"]}:{coordinates}'

    return "; ".join(
        f'{point_text(segment["start"])} -> {point_text(segment["end"])}'
        for segment in path
    )


__all__ = ["HighSymmetryPathError", "format_band_path", "high_symmetry_k_path"]
