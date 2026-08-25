"""Crystallographic catalogs backed by the bundled MSGCorep dependency.

The public MagneticKP solver does not require MSGCorep.  Catalog-driven
workflows import it lazily so that the ordinary matrix-input API remains a
small SymPy-only dependency.
"""

from .msgcorep_adapter import (
    CatalogDependencyError,
    CatalogInputError,
    build_catalog_core_input,
    catalog_status,
    corepresentation_catalog,
    exact_to_sympy,
    k_point_catalog,
    list_bz_variants,
    list_magnetic_groups,
    list_space_groups,
    parse_magnetic_identifier,
)

__all__ = [
    "CatalogDependencyError",
    "CatalogInputError",
    "build_catalog_core_input",
    "catalog_status",
    "corepresentation_catalog",
    "exact_to_sympy",
    "k_point_catalog",
    "list_bz_variants",
    "list_magnetic_groups",
    "list_space_groups",
    "parse_magnetic_identifier",
]
