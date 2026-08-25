"""Translate MSGCorep catalogs into native MagneticKP solver inputs.

This module is the only integration boundary between the two projects.  It
keeps MSGCorep's exact, database-backed crystallographic model out of the Web
layer and keeps MagneticKP's solver independent of MSGCorep's internal AST.
"""

from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
import importlib
import re
from typing import Mapping, Sequence

import sympy
from sympy import Matrix, diag
from sympy.matrices import MatrixBase


class CatalogDependencyError(RuntimeError):
    """Raised when a catalog workflow is requested without MSGCorep."""


class CatalogInputError(ValueError):
    """Raised when a crystallographic catalog selection is inconsistent."""


@lru_cache(maxsize=1)
def _msgcorep():
    try:
        return importlib.import_module("msgcorep")
    except ModuleNotFoundError as error:
        raise CatalogDependencyError(
            "MSGCorep is required for crystallographic catalog workflows; "
            "install MagneticKP with its web/catalog dependencies"
        ) from error


@lru_cache(maxsize=1)
def _repository():
    try:
        module = importlib.import_module("msgcorep.data")
    except ModuleNotFoundError as error:
        raise CatalogDependencyError(
            "MSGCorep is required for crystallographic catalog workflows"
        ) from error
    return module.get_repository()


def catalog_status():
    """Return availability and version without making package import global."""

    try:
        package = _msgcorep()
    except CatalogDependencyError as error:
        return {"available": False, "version": None, "detail": str(error)}
    return {
        "available": True,
        "version": getattr(package, "__version__", "unknown"),
        "detail": None,
    }


_CALLS = {
    "Abs": sympy.Abs,
    "ArcCos": sympy.acos,
    "ArcSin": sympy.asin,
    "ArcTan": sympy.atan,
    "Ceiling": sympy.ceiling,
    "Conjugate": sympy.conjugate,
    "Cos": sympy.cos,
    "Cot": sympy.cot,
    "Csc": sympy.csc,
    "Exp": sympy.exp,
    "Floor": sympy.floor,
    "Im": sympy.im,
    "Log": sympy.log,
    "Mod": sympy.Mod,
    "Re": sympy.re,
    "Sin": sympy.sin,
    "Sqrt": sympy.sqrt,
    "Tan": sympy.tan,
}


def exact_to_sympy(value, symbols: dict[str, sympy.Symbol] | None = None):
    """Convert MSGCorep's immutable exact AST into a SymPy value.

    Symbols other than ``E`` and ``Pi`` are real crystallographic parameters.
    The conversion is structural and never parses or evaluates source text.
    """

    symbol_table = {} if symbols is None else symbols
    if isinstance(value, sympy.Basic):
        return value
    if isinstance(value, bool):
        raise CatalogInputError("boolean values are not crystallographic scalars")
    if isinstance(value, int):
        return sympy.Integer(value)
    if isinstance(value, Fraction):
        return sympy.Rational(value.numerator, value.denominator)
    if isinstance(value, float):
        return sympy.Float(str(value))
    if isinstance(value, complex):
        return exact_to_sympy(value.real, symbol_table) + sympy.I * exact_to_sympy(
            value.imag, symbol_table
        )
    if isinstance(value, (tuple, list)):
        return tuple(exact_to_sympy(item, symbol_table) for item in value)

    head = getattr(value, "head", None)
    args = getattr(value, "args", None)
    if not isinstance(head, str) or not isinstance(args, tuple):
        raise CatalogInputError(
            f"unsupported MSGCorep exact value {type(value).__name__}"
        )
    if head == "Integer":
        return sympy.Integer(args[0])
    if head == "Rational":
        return sympy.Rational(args[0], args[1])
    if head == "RealDecimal":
        return sympy.Float(args[0])
    if head == "Symbol":
        name = str(args[0])
        if name == "E":
            return sympy.E
        if name == "Pi":
            return sympy.pi
        return symbol_table.setdefault(name, sympy.Symbol(name, real=True))
    if head == "Complex":
        return exact_to_sympy(args[0], symbol_table) + sympy.I * exact_to_sympy(
            args[1], symbol_table
        )
    if head == "Add":
        return sympy.Add(*(exact_to_sympy(item, symbol_table) for item in args))
    if head == "Mul":
        return sympy.Mul(*(exact_to_sympy(item, symbol_table) for item in args))
    if head == "Pow":
        return exact_to_sympy(args[0], symbol_table) ** exact_to_sympy(
            args[1], symbol_table
        )
    if head == "Call":
        name, arguments = args
        try:
            function = _CALLS[name]
        except KeyError as error:
            raise CatalogInputError(
                f"unsupported MSGCorep exact function {name!r}"
            ) from error
        return function(*(exact_to_sympy(item, symbol_table) for item in arguments))
    if head == "Matrix":
        return Matrix(
            [
                [exact_to_sympy(item, symbol_table) for item in row]
                for row in args[0]
            ]
        )
    if head == "List":
        return tuple(exact_to_sympy(item, symbol_table) for item in args[0])
    raise CatalogInputError(f"unsupported MSGCorep exact-expression head {head!r}")


def parse_magnetic_identifier(value) -> tuple[int, int]:
    """Normalize a BNS number written as ``92.112`` or ``(92, 112)``."""

    if isinstance(value, str):
        match = re.fullmatch(r"\s*(\d{1,3})\.(\d+)\s*", value)
        if not match:
            raise CatalogInputError(
                "magnetic group must be a BNS number such as 92.112"
            )
        identifier = int(match.group(1)), int(match.group(2))
    else:
        try:
            identifier = tuple(value)
        except TypeError as error:
            raise CatalogInputError(
                "magnetic group must be a BNS number such as 92.112"
            ) from error
        if len(identifier) != 2:
            raise CatalogInputError("magnetic group must contain two integers")
    if any(isinstance(item, bool) or not isinstance(item, int) for item in identifier):
        raise CatalogInputError("magnetic group must contain two integers")
    normalized = int(identifier[0]), int(identifier[1])
    _repository().get_magnetic_space_group(normalized)
    return normalized


def list_space_groups():
    """Return the compact 230-entry space-group catalog."""

    rows = _repository().get_space_group_catalog()
    return [
        {
            "number": row.number,
            "symbol": row.symbol_text,
            "symbol_tex": row.symbol_tex,
            "bc_symbol": row.symbol_text_bc,
            "bc_symbol_tex": row.symbol_tex_bc,
            "bravais": row.bravais_name,
            "schoenflies": f"{row.schoenflies_family}{row.schoenflies_index}",
            "orientation": row.orientation_label,
        }
        for row in rows
    ]


def list_magnetic_groups(space_group: int):
    """Return all BNS magnetic groups belonging to one SG family."""

    if isinstance(space_group, bool) or not isinstance(space_group, int):
        raise CatalogInputError("space group must be an integer")
    rows = _repository().get_magnetic_space_groups_for_family(space_group)
    return [
        {
            "identifier": f"{row.family_number}.{row.bns_number}",
            "family_number": row.family_number,
            "bns_number": row.bns_number,
            "type": row.msg_type,
            "bns_symbol": row.bns_symbol,
            "bns_symbol_tex": row.bns_symbol_tex,
            "bc_symbol": row.bc_symbol,
            "bc_symbol_tex": row.bc_symbol_tex,
            "og_number": row.og_number,
            "og_symbol": row.og_symbol,
            "og_symbol_tex": row.og_symbol_tex,
        }
        for row in rows
    ]


def _full_bz_type(space_group: int, variant: str = "a") -> str:
    package = _msgcorep()
    bravais = package.getSGLatt(space_group)
    if bravais in package.BCHighSymKpt:
        return bravais
    candidate = f"{bravais}({variant})"
    if candidate in package.BCHighSymKpt:
        return candidate
    variants = list_bz_variants(space_group)
    if not variants:
        raise CatalogInputError(f"no Brillouin-zone catalog for SG {space_group}")
    return str(variants[0]["full_bz_type"])


def list_bz_variants(space_group: int):
    package = _msgcorep()
    bravais = package.getSGLatt(space_group)
    matching = [
        name
        for name in package.BCHighSymKpt
        if name == bravais or name.startswith(f"{bravais}(")
    ]
    return [
        {
            "value": (
                full_name[len(bravais) + 1 : -1]
                if full_name != bravais
                else "a"
            ),
            "label": full_name,
            "full_bz_type": full_name,
        }
        for full_name in matching
    ]


def _coordinate_data(value):
    symbols: dict[str, sympy.Symbol] = {}
    converted = tuple(exact_to_sympy(item, symbols) for item in value)
    parameters = sorted(
        set().union(*(item.free_symbols for item in converted)), key=str
    )
    return {
        "coordinate": [sympy.sstr(item) for item in converted],
        "coordinate_latex": (
            sympy.latex(Matrix(1, len(converted), converted))
            + r"^{\mathsf T}"
        ),
        "parameters": [str(item) for item in parameters],
    }


def k_point_catalog(
    space_group: int,
    *,
    magnetic_group=None,
    variant: str = "a",
):
    """Return standard BZ points/lines with exact affine coordinates."""

    package = _msgcorep()
    full_name = _full_bz_type(space_group, variant)
    available = None
    if magnetic_group is not None:
        identifier = parse_magnetic_identifier(magnetic_group)
        if identifier[0] != space_group:
            raise CatalogInputError(
                "magnetic group does not belong to the selected space group"
            )
        available = set(_repository().get_mlg_corep_k_names(identifier))
    rows = []
    for name, path, coordinate, point_group, _operations in package.BCHighSymKpt[
        full_name
    ]:
        lookup_names = {str(name)}
        if str(name) == "F":
            lookup_names.update({"aF", "bF"})
        if available is not None and lookup_names.isdisjoint(available):
            continue
        data = _coordinate_data(coordinate)
        rows.append(
            {
                "name": str(name),
                "path": str(path) if path else "",
                "kind": "line" if data["parameters"] or path else "point",
                "point_group": str(point_group),
                **data,
            }
        )
    return {"bz_type": full_name, "points": rows}


def _fraction(value: str, name: str) -> Fraction:
    try:
        result = Fraction(str(value).strip())
    except (ValueError, ZeroDivisionError) as error:
        raise CatalogInputError(
            f"k-point parameter {name!r} must be an exact integer, decimal, or fraction"
        ) from error
    if result.denominator > 10**6 or abs(result.numerator) > 10**9:
        raise CatalogInputError(f"k-point parameter {name!r} is too large")
    return result


def _resolved_corep_info(
    identifier: tuple[int, int],
    k_name: str,
    variant: str,
    parameter_values: Mapping[str, str] | None,
):
    package = _msgcorep()
    try:
        symbolic = package.getMLGCorep(identifier, k_name, variant, format=False)
    except package.InputError as error:
        raise CatalogInputError(str(error)) from error
    coordinate_data = _coordinate_data(symbolic["k"][0])
    required = tuple(coordinate_data["parameters"])
    supplied = dict(parameter_values or {})
    if not required:
        if supplied:
            raise CatalogInputError("the selected k point has no free coordinate parameter")
        return symbolic, coordinate_data
    if not supplied:
        return symbolic, coordinate_data
    missing = [name for name in required if name not in supplied]
    extra = [name for name in supplied if name not in required]
    if missing or extra:
        details = []
        if missing:
            details.append("missing " + ", ".join(missing))
        if extra:
            details.append("unknown " + ", ".join(extra))
        raise CatalogInputError("invalid k-point parameters: " + "; ".join(details))
    substitutions = {
        sympy.Symbol(name, real=True): sympy.Rational(
            _fraction(supplied[name], name).numerator,
            _fraction(supplied[name], name).denominator,
        )
        for name in required
    }
    symbolic_coordinate = [
        exact_to_sympy(item) for item in symbolic["k"][0]
    ]
    coordinate = []
    for value in symbolic_coordinate:
        resolved = sympy.simplify(value.subs(substitutions))
        if not resolved.is_Rational:
            raise CatalogInputError(
                "resolved k-point coordinates must be exact rational numbers"
            )
        coordinate.append(Fraction(int(resolved.p), int(resolved.q)))
    try:
        numeric = package.getMLGCorep(
            identifier, tuple(coordinate), variant, format=False
        )
    except package.InputError as error:
        raise CatalogInputError(str(error)) from error
    return numeric, {
        "coordinate": [str(item) for item in coordinate],
        "coordinate_latex": sympy.latex(
            Matrix(
                1,
                len(coordinate),
                [
                    sympy.Rational(item.numerator, item.denominator)
                    for item in coordinate
                ],
            )
        ) + r"^{\mathsf T}",
        "parameters": list(required),
        "parameter_values": {name: str(_fraction(supplied[name], name)) for name in required},
    }


def _matrix(value, symbols=None) -> Matrix:
    if isinstance(value, MatrixBase):
        return Matrix(value)
    converted = exact_to_sympy(value, symbols)
    if isinstance(converted, MatrixBase):
        return Matrix(converted)
    if isinstance(converted, tuple) and converted and isinstance(converted[0], tuple):
        return Matrix(converted)
    return Matrix([[converted]])


def _representation_dimension(representation) -> int:
    return _matrix(representation[0]).rows


def corepresentation_catalog(
    magnetic_group,
    k_name: str,
    *,
    variant: str = "a",
    parameter_values: Mapping[str, str] | None = None,
    coordinate_basis: str = "cartesian",
):
    """Return selectable single- and double-valued MLG corepresentations."""

    identifier = parse_magnetic_identifier(magnetic_group)
    if coordinate_basis not in {"cartesian", "fractional"}:
        raise CatalogInputError(
            "coordinate basis must be cartesian or fractional"
        )
    info, coordinate = _resolved_corep_info(
        identifier, k_name, variant, parameter_values
    )

    def family(prefix: str):
        rows = []
        for index, (label, corep_type, representation) in enumerate(
            zip(
                info[prefix + "label"],
                info[prefix + "type"],
                info[prefix + "corep"],
            ),
            start=1,
        ):
            generator_matrices = []
            for operation_index in generator_indices:
                matrix = _matrix(representation[operation_index])
                generator_matrices.append(
                    {
                        "generator_index": operation_index + 1,
                        "antiunitary": bool(operations[operation_index][2]),
                        "matrix": [
                            [
                                sympy.sstr(matrix[row, column])
                                for column in range(matrix.cols)
                            ]
                            for row in range(matrix.rows)
                        ],
                        "matrix_latex": sympy.latex(matrix),
                    }
                )
            rows.append(
                {
                    "index": index,
                    "label": str(label),
                    "type": str(corep_type),
                    "dimension": _representation_dimension(representation),
                    "generator_matrices": generator_matrices,
                }
            )
        return rows

    operations = tuple(info["MLG"])
    generator_indices = _generator_indices(identifier, operations)
    package = _msgcorep()
    bravais = package.getSGLatt(identifier[0])
    momentum = Matrix(sympy.symbols("kx ky kz", real=True))
    generators = []
    for operation_index in generator_indices:
        operation = operations[operation_index]
        action = _momentum_action_matrix(
            package,
            bravais,
            operation,
            coordinate_basis,
        )
        image = action * momentum
        generators.append(
            {
                "index": operation_index + 1,
                "label": _operation_label(operation_index, operation),
                "rotation": str(operation[0]),
                "translation": [str(item) for item in operation[1]],
                "antiunitary": bool(operation[2]),
                "momentum_matrix": [
                    [sympy.sstr(action[row, column]) for column in range(3)]
                    for row in range(3)
                ],
                "momentum_image": [sympy.sstr(item) for item in image],
                "momentum_image_latex": (
                    sympy.latex(Matrix(1, 3, list(image)))
                    + r"^{\mathsf T}"
                ),
            }
        )

    return {
        "magnetic_group": f"{identifier[0]}.{identifier[1]}",
        "k_name": str(info["k"][1]),
        "coordinate": coordinate,
        "group_order": len(info["MLG"]),
        "coordinate_basis": coordinate_basis,
        "generators": generators,
        "single": family("s"),
        "double": family("d"),
    }


def _operation_key(operation) -> tuple[str, int]:
    return str(operation[0]), int(operation[2])


def _generator_indices(identifier: tuple[int, int], operations: Sequence) -> tuple[int, ...]:
    """Find a small generating subset of the magnetic little co-group."""

    if len(operations) <= 1:
        return (0,)
    package = _msgcorep()
    bravais = package.getSGLatt(identifier[0])
    multiply = package.MSGSeitzTimes(bravais)
    key_to_index = {_operation_key(operation): index for index, operation in enumerate(operations)}
    if len(key_to_index) != len(operations):
        return tuple(range(1, len(operations)))
    identity = _operation_key(operations[0])

    def product(left, right):
        value = multiply(
            operations[key_to_index[left]], operations[key_to_index[right]]
        )
        key = _operation_key(value)
        if key not in key_to_index:
            raise CatalogInputError("magnetic little-group catalog is not closed")
        return key

    def closure(generators):
        seen = {identity}
        queue = [identity]
        while queue:
            current = queue.pop()
            for generator in generators:
                value = product(current, generator)
                if value not in seen:
                    seen.add(value)
                    queue.append(value)
        return seen

    candidates = [key for key in key_to_index if key != identity]
    selected: list[tuple[str, int]] = []
    generated = {identity}
    while len(generated) < len(operations):
        best = None
        best_closure = generated
        for candidate in candidates:
            if candidate in generated:
                continue
            trial = closure(selected + [candidate])
            if len(trial) > len(best_closure):
                best = candidate
                best_closure = trial
        if best is None:
            return tuple(range(1, len(operations)))
        selected.append(best)
        generated = best_closure
    for candidate in tuple(reversed(selected)):
        trial_generators = [item for item in selected if item != candidate]
        if len(closure(trial_generators)) == len(operations):
            selected = trial_generators
    return tuple(key_to_index[key] for key in selected)


def _operation_label(index: int, operation) -> str:
    rotation, translation, antiunitary = operation
    translation_text = ",".join(str(item) for item in translation)
    suffix = "T" if antiunitary else ""
    return f"g{index + 1}:{rotation}|{translation_text}{suffix}"


def _momentum_action_matrix(package, bravais: str, operation, coordinate_basis: str):
    """Return the exact momentum action used by the MagneticKP constraint."""

    rotation_name = str(operation[0]).removeprefix("bar")
    if coordinate_basis == "cartesian":
        try:
            momentum_rotation = package.RotMatCart[rotation_name]
        except KeyError as error:
            raise CatalogInputError(
                f"Cartesian momentum rotation {rotation_name!r} is absent"
            ) from error
    else:
        try:
            momentum_rotation = package.getRotMatOfK(bravais, rotation_name)
        except (KeyError, ValueError) as error:
            raise CatalogInputError(
                f"fractional momentum rotation {rotation_name!r} is absent"
            ) from error
    rotation = _matrix(momentum_rotation)
    return -rotation if bool(operation[2]) else rotation


def build_catalog_core_input(
    *,
    space_group: int,
    magnetic_group,
    k_name: str,
    representation_family: str,
    representations: Sequence[Mapping[str, int]],
    variant: str = "a",
    parameter_values: Mapping[str, str] | None = None,
    coordinate_basis: str = "cartesian",
    method: str = "automatic",
    maximum_dimension: int = 12,
):
    """Build ``kpHam`` input and selection metadata from catalog identifiers."""

    identifier = parse_magnetic_identifier(magnetic_group)
    if identifier[0] != space_group:
        raise CatalogInputError(
            "magnetic group does not belong to the selected space group"
        )
    if representation_family not in {"single", "double"}:
        raise CatalogInputError("representation family must be single or double")
    if coordinate_basis not in {"cartesian", "fractional"}:
        raise CatalogInputError(
            "coordinate basis must be cartesian or fractional"
        )
    if not representations:
        raise CatalogInputError("select at least one corepresentation")
    info, coordinate = _resolved_corep_info(
        identifier, k_name, variant, parameter_values
    )
    prefix = "s" if representation_family == "single" else "d"
    available = tuple(info[prefix + "corep"])
    labels = tuple(str(item) for item in info[prefix + "label"])
    selected = []
    seen = set()
    for item in representations:
        try:
            index = int(item["index"])
            multiplicity = int(item.get("multiplicity", 1))
        except (KeyError, TypeError, ValueError) as error:
            raise CatalogInputError(
                "corepresentations require integer index and multiplicity"
            ) from error
        if index in seen:
            raise CatalogInputError(f"corepresentation index {index} is duplicated")
        if not 1 <= index <= len(available):
            raise CatalogInputError(
                f"corepresentation index {index} is out of range"
            )
        if not 1 <= multiplicity <= 8:
            raise CatalogInputError("corepresentation multiplicity must be from 1 to 8")
        seen.add(index)
        selected.append((index - 1, multiplicity))

    symbol_table: dict[str, sympy.Symbol] = {}
    dimensions = {
        index: _representation_dimension(available[index]) for index, _ in selected
    }
    total_dimension = sum(
        dimensions[index] * multiplicity for index, multiplicity in selected
    )
    if total_dimension > maximum_dimension:
        raise CatalogInputError(
            f"selected direct sum has dimension {total_dimension}; maximum is {maximum_dimension}"
        )
    operations = tuple(info["MLG"])
    if str(method).lower() == "reynolds":
        operation_indices = tuple(range(len(operations)))
        group_input = "elements"
    else:
        operation_indices = _generator_indices(identifier, operations)
        group_input = "generators"
    package = _msgcorep()
    bravais = package.getSGLatt(space_group)
    core_input: dict[str, dict[str, tuple[Matrix, Matrix]]] = {}
    used_operations = []
    for operation_index in operation_indices:
        operation = operations[operation_index]
        blocks = []
        for representation_index, multiplicity in selected:
            matrix = _matrix(
                available[representation_index][operation_index], symbol_table
            )
            blocks.extend(matrix for _ in range(multiplicity))
        representation_matrix = blocks[0] if len(blocks) == 1 else diag(*blocks)
        momentum_matrix = _momentum_action_matrix(
            package,
            bravais,
            operation,
            coordinate_basis,
        )
        antiunitary = bool(operation[2])
        category = "Antiunitary" if antiunitary else "Unitary"
        label = _operation_label(operation_index, operation)
        core_input.setdefault(category, {})[label] = (
            representation_matrix,
            momentum_matrix,
        )
        used_operations.append(
            {
                "index": operation_index + 1,
                "label": label,
                "rotation": str(operation[0]),
                "translation": [str(item) for item in operation[1]],
                "antiunitary": antiunitary,
            }
        )
    metadata = {
        "space_group": space_group,
        "magnetic_group": f"{identifier[0]}.{identifier[1]}",
        "k_name": str(info["k"][1]),
        "coordinate": coordinate,
        "representation_family": representation_family,
        "representations": [
            {
                "index": index + 1,
                "label": labels[index],
                "multiplicity": multiplicity,
                "dimension": dimensions[index],
            }
            for index, multiplicity in selected
        ],
        "dimension": total_dimension,
        "group_order": len(operations),
        "constraint_count": len(operation_indices),
        "group_input": group_input,
        "momentum_basis": coordinate_basis,
        "operations": used_operations,
    }
    core_input.update(
        {
            "KPoint": coordinate["coordinate"],
            "KPointLabel": str(info["k"][1]),
            "CoordinateSystem": (
                "Cartesian"
                if coordinate_basis == "cartesian"
                else "ReciprocalFractional"
            ),
            "Group": (identifier[0], identifier[1]),
        }
    )
    return core_input, metadata, group_input
