"""Python interface for constructing symmetry-constrained k.p Hamiltonians."""

from __future__ import annotations

from typing import Dict, Iterable, List, Mapping, MutableMapping, Optional, Sequence

from sympy import Matrix, Symbol, eye, symbols, zeros
from sympy.matrices import MatrixBase

from .basis import (
    adjoint_representation_matrix,
    basis_data,
    continuous_generator,
    continuous_invariant_basis_data,
    continuous_normalizer,
)
from .solvers import (
    iterative_invariant_basis,
    reynolds_invariant_basis,
    stacked_invariant_basis,
)
from .transforms import K_SYMBOLS, homogeneous_monomials, kx, ky, kz, momentum_transform_matrix


__author__ = "Zhang Zeying"


_COEFFICIENT_PREFIXES = ("e", "t", "r", "s")


def _coefficient_symbols(order: int, count: int):
    """Return stable, one-based coefficient names for low momentum orders."""

    if order < len(_COEFFICIENT_PREFIXES):
        prefix = _COEFFICIENT_PREFIXES[order]
        return [symbols(f"{prefix}_{line}") for line in range(1, count + 1)]
    return [symbols(f"C_{line}_{order}") for line in range(count)]


class kpHam(object):
    r"""Construct a k.p Hamiltonian from representation and momentum matrices.

    ``inputdict`` accepts the legacy key ``"Anitunitary"`` and the corrected
    spelling ``"Antiunitary"``.  Each operation maps to ``(D, R)`` where ``D``
    is the (co)representation matrix and ``R`` is either a momentum matrix or
    a list of linear momentum images.  A single internal continuous symmetry
    may be placed in ``"Unitary"`` as
    ``"CInfinity": (D_theta, (kx, ky, kz), theta)``.  Antisymmetry keys are
    accepted with the same spelling conventions as the Wolfram package.
    """

    _CATEGORY_NAMES = (
        ("Unitary",),
        ("Anitunitary", "Antiunitary"),
        ("AntisymmetryUnitaryTest",),
        ("AntisymmetryAnitunitaryTest", "AntisymmetryAntiunitaryTest"),
    )

    def __init__(self, inputdict: Mapping):
        self.input = dict(inputdict)
        self.debug = False
        continuous_record = self._continuous_record()
        raw_records = self._raw_records()
        if not raw_records and continuous_record is None:
            raise ValueError("at least one valid symmetry operation is required")
        first_record = continuous_record or raw_records[0]
        self.dim = first_record["D"].rows
        all_records = (
            ([continuous_record] if continuous_record is not None else [])
            + raw_records
        )
        if any(record["D"].shape != (self.dim, self.dim) for record in all_records):
            raise ValueError("all representation matrices must have the same dimension")
        self.spatial_dim = self._spatial_dimension(first_record["k_image"])
        if self.spatial_dim not in (1, 2, 3):
            raise ValueError("the momentum-space dimension must be 1, 2 or 3")
        if any(
            self._spatial_dimension(record["k_image"]) != self.spatial_dim
            for record in all_records
        ):
            raise ValueError("all momentum transformations must have the same dimension")

        full_basis_data = basis_data(self.dim)
        if continuous_record is not None:
            expected_k_image = Matrix(self._variables())
            actual_k_image = self._continuous_k_image_vector(
                continuous_record["k_image"]
            )
            if actual_k_image is None or not actual_k_image.equals(expected_k_image):
                raise ValueError(
                    "internal CInfinity symmetry must leave momentum unchanged"
                )
            generator_data = continuous_generator(
                continuous_record["D"],
                continuous_record["parameter"],
            )
            reduced = continuous_invariant_basis_data(
                generator_data["generator"],
                full_basis_data,
            )
            for record in raw_records:
                if not continuous_normalizer(
                    record["D"],
                    generator_data["generator"],
                    record["antiunitary"],
                ):
                    raise ValueError(
                        f"discrete operation {record['label']!r} does not "
                        "normalize the CInfinity generator"
                    )
            self._basis_data = reduced["basis_data"]
            self._continuous = {
                **continuous_record,
                **generator_data,
                **reduced,
                "original_operator_dimension": full_basis_data["count"],
                "reduced_operator_dimension": reduced["basis_data"]["count"],
            }
        else:
            self._basis_data = full_basis_data
            self._continuous = None
        self._basis = self._basis_data["basis"]
        self._records = []
        for record in raw_records:
            prepared = dict(record)
            prepared["J"] = adjoint_representation_matrix(
                record["D"],
                self._basis_data,
                record["antiunitary"],
            )
            self._records.append(prepared)

    @staticmethod
    def _spatial_dimension(k_image) -> int:
        if isinstance(k_image, MatrixBase):
            if k_image.cols == k_image.rows:
                return k_image.rows
            if k_image.cols == 1:
                return k_image.rows
        return len(k_image)

    @staticmethod
    def _continuous_label(label) -> bool:
        return str(label).lower() in {"cinfinity", "cinfty"}

    @staticmethod
    def _continuous_k_image_vector(k_image):
        if isinstance(k_image, MatrixBase):
            if k_image.cols == 1:
                return Matrix(k_image)
            if k_image.rows == 1:
                return Matrix(k_image).T
            return None
        try:
            return Matrix(list(k_image))
        except (TypeError, ValueError):
            return None

    def _continuous_record(self):
        entries = []
        for category, names in enumerate(self._CATEGORY_NAMES, start=1):
            for name in names:
                values = self.input.get(name, {})
                if not isinstance(values, Mapping):
                    continue
                for label, value in values.items():
                    if self._continuous_label(label):
                        entries.append((category, label, value))
        if not entries:
            return None
        if any(category != 1 for category, _, _ in entries):
            raise ValueError("CInfinity is allowed only in the Unitary category")
        if len(entries) != 1:
            raise ValueError("exactly one CInfinity entry is supported")
        _, label, value = entries[0]
        if not isinstance(value, (tuple, list)) or len(value) != 3:
            raise ValueError(
                "CInfinity entry must be (D_theta, (kx, ky, kz), theta)"
            )
        representation, k_image, parameter = value
        if not isinstance(parameter, Symbol):
            raise ValueError("the CInfinity third entry must be a symbol")
        return {
            "label": label,
            "D": Matrix(representation),
            "k_image": k_image,
            "parameter": parameter,
            "type": 1,
            "antiunitary": False,
            "eta": 1,
        }

    def _raw_records(self) -> List[dict]:
        records = []
        for category, names in enumerate(self._CATEGORY_NAMES, start=1):
            merged = {}
            for name in names:
                value = self.input.get(name, {})
                if isinstance(value, Mapping):
                    merged.update(value)
            antiunitary = category in (2, 4)
            eta = 1 if category in (1, 2) else -1
            for label, value in merged.items():
                if self._continuous_label(label):
                    continue
                if not isinstance(value, (tuple, list)) or len(value) < 2:
                    continue
                representation = Matrix(value[0])
                records.append(
                    {
                        "label": label,
                        "D": representation,
                        "k_image": value[1],
                        "type": category,
                        "antiunitary": antiunitary,
                        "eta": eta,
                    }
                )
        return records

    def _variables(self):
        if self.spatial_dim == 1:
            return (kz,)
        if self.spatial_dim == 2:
            return (kx, ky)
        return K_SYMBOLS

    def getklist(self, order):
        """Generate the homogeneous momentum basis in legacy ordering."""

        return homogeneous_monomials(order, self.spatial_dim)

    def ISAAlg(self, constraints):
        """Compatibility implementation of the iterative common-nullspace solver."""

        if not constraints:
            return Matrix([])
        ambient = constraints[0].cols
        candidate = eye(ambient)
        for constraint in constraints:
            kernel = (constraint * candidate).nullspace()
            if not kernel:
                return zeros(ambient, 0)
            candidate = candidate * Matrix.hstack(*kernel)
        return candidate

    def to_Ham(self, order, kbases, hbases, column_basis):
        """Convert an invariant column basis to a symbolic Hamiltonian."""

        ambient, count = column_basis.shape
        coefficients = _coefficient_symbols(order, count)
        product_basis = [
            monomial * operator
            for monomial in kbases
            for operator in hbases
        ]
        if len(product_basis) != ambient:
            raise RuntimeError("coefficient-space and product-basis dimensions disagree")
        hamiltonian = zeros(self.dim)
        for column in range(count):
            term = zeros(self.dim)
            for row, basis_matrix in enumerate(product_basis):
                term += column_basis[row, column] * basis_matrix
            hamiltonian += coefficients[column] * term
        return hamiltonian

    def _factors(self, order):
        kbases = self.getklist(order)
        variables = self._variables()
        factors = []
        for record in self._records:
            factors.append(
                {
                    "F": momentum_transform_matrix(
                        order,
                        record["k_image"],
                        variables=variables,
                        klist=kbases,
                    ),
                    "J": record["J"],
                    "eta": record["eta"],
                }
            )
        return kbases, factors

    def _resolved_group_input(self, group_input):
        if group_input is None or str(group_input).lower() == "automatic":
            group_input = self.input.get("GroupInput", "Generators")
        return "elements" if str(group_input).lower() == "elements" else "generators"

    @staticmethod
    def _normalise_method(method):
        text = str(method).lower()
        if text in ("reynolds",):
            return "reynolds"
        if text in ("automatic", "auto"):
            return "automatic"
        if text in (
            "directproduct",
            "direct_product",
            "directproductdecomposition",
            "stacked",
        ):
            return "direct_product"
        if text in (
            "iterative",
            "iterativelysimplify",
            "iterativesimplification",
            "iterativesimplification2016",
        ):
            return "iterative"
        raise ValueError(
            "method must be iterative, direct_product, reynolds, or automatic"
        )

    def IterativelySimplify(
        self,
        korder,
        method="iterative",
        group_input=None,
        max_group_order=512,
        return_diagnostics=False,
    ):
        """Construct the Hamiltonian at one homogeneous momentum order.

        The original ``IterativelySimplify(korder) -> (ham, count)`` call is
        unchanged.  The ``direct_product`` and ``reynolds`` methods select
        those backends explicitly.  ``automatic`` conservatively uses the
        generator-based iterative solver.
        """

        if not isinstance(korder, int) or korder < 0:
            raise ValueError("korder must be a non-negative integer")
        method = self._normalise_method(method)
        kbases, factors = self._factors(korder)
        ambient = len(kbases) * self._basis_data["count"]
        chosen = "iterative" if method == "automatic" else method
        if not factors:
            column_basis = eye(ambient)
            diagnostics = {
                "solver": "none",
                "fallback": False,
                "group_order": 0,
                "reason": None,
            }
        elif chosen == "reynolds":
            column_basis, diagnostics = reynolds_invariant_basis(
                factors,
                group_input=self._resolved_group_input(group_input),
                max_group_order=max_group_order,
            )
        elif chosen == "direct_product":
            column_basis = stacked_invariant_basis(factors)
            diagnostics = {
                "solver": "direct_product",
                "fallback": False,
                "group_order": None,
                "reason": None,
            }
        else:
            column_basis = iterative_invariant_basis(factors)
            diagnostics = {
                "solver": "iterative",
                "fallback": False,
                "group_order": None,
                "reason": None,
            }
        diagnostics.update(
            order=korder,
            ambient_dimension=ambient,
            invariant_dimension=column_basis.cols,
            continuous_symmetry=(
                "CInfinity" if self._continuous is not None else None
            ),
            original_operator_dimension=(
                self._continuous["original_operator_dimension"]
                if self._continuous is not None
                else self._basis_data["count"]
            ),
            reduced_operator_dimension=self._basis_data["count"],
        )
        hamiltonian = self.to_Ham(korder, kbases, self._basis, column_basis)
        if return_diagnostics:
            return hamiltonian, column_basis.cols, diagnostics
        return hamiltonian, column_basis.cols

    def DirectProductDecomposition(
        self,
        korder,
        return_diagnostics=False,
    ):
        """Construct one homogeneous order using stacked constraints."""

        return self.IterativelySimplify(
            korder,
            method="direct_product",
            return_diagnostics=return_diagnostics,
        )

    def getkpHam(
        self,
        order,
        method="iterative",
        group_input=None,
        max_group_order=512,
        return_diagnostics=False,
    ):
        """Construct the k.p model through the requested momentum orders.

        ``method="automatic"`` conservatively selects the iterative backend;
        request ``method="reynolds"`` explicitly to use group averaging.
        """

        if isinstance(order, int):
            if order < 0:
                raise ValueError("order must be non-negative")
            orders = list(range(order + 1))
        else:
            orders = list(order)
            if any(not isinstance(value, int) or value < 0 for value in orders):
                raise ValueError("all requested orders must be non-negative integers")

        hamiltonian = zeros(self.dim)
        number_of_parameters = 0
        all_diagnostics = []
        for current_order in orders:
            result = self.IterativelySimplify(
                current_order,
                method=method,
                group_input=group_input,
                max_group_order=max_group_order,
                return_diagnostics=return_diagnostics,
            )
            if return_diagnostics:
                term, count, diagnostics = result
                all_diagnostics.append(diagnostics)
            else:
                term, count = result
            hamiltonian += term
            number_of_parameters += count

        output = {
            "ham": hamiltonian,
            "NumberOfParameters": number_of_parameters,
            "korder": orders,
            "dim": self.dim,
        }
        if return_diagnostics:
            output["Diagnostics"] = all_diagnostics
        return output


if __name__ == "__main__":
    example = kpHam(
        {
            "Unitary": {
                "C4": (
                    Matrix([[1, 0], [0, -1]]),
                    Matrix([[0, -1, 0], [1, 0, 0], [0, 0, 1]]),
                )
            },
            "Antiunitary": {
                "IT": (
                    eye(2),
                    eye(3),
                )
            },
        }
    )
    print(example.getkpHam(2))
