# Python catalog / Mathematica parity

`PythonCatalogParity.wls` and
`Python-Web/tests/test_mathematica_catalog_parity.py` compare the new
MSGCorep-backed Python catalog workflow with the established Mathematica path

```wl
interfaceRep[msg, k, reps,
  "CartesianCoordinates" -> basis,
  "CalculateGenerators" -> True]
kpHam[order, input,
  "Method" -> "IterativeSimplification",
  "ReturnDiagnostics" -> True]
```

The comparison is exact.  It checks the matrix dimension, parameter count,
per-order invariant dimensions, momentum exponent set, and the full linear
space spanned by all Hamiltonian coefficient matrices.  Comparing coefficient
spaces instead of printed expressions allows Mathematica and Python to choose
different, but equivalent, bases for the free parameters.

## Covered cases

| Case | Dimension | Parameters | Per-order invariant dimensions | Exact subspace |
| --- | ---: | ---: | --- | --- |
| type-I MSG 13.65, B, double, Cartesian | 2 | 4 | 1, 3 | equal |
| type-I MSG 13.65, B, double, fractional | 2 | 4 | 1, 3 | equal |
| type-II MSG 92.112, Γ, double | 2 | 5 | 1, 2, 2 | equal |
| type-III MSG 92.113, Γ, double | 2 | 7 | 1, 2, 4 | equal |
| type-IV MSG 92.116, Γ, double | 2 | 5 | 1, 2, 2 | equal |
| type-IV MSG 4.10, Λ with `u=1/5`, fractional | 1 | 6 | 1, 1, 4 | equal |
| MSG 97.155, P, single direct sum 1+2 | 3 | 9 | 2, 7 | equal |
| MSG 92.112, Γ, double corep 1 with multiplicity 2 | 4 | 9 | 3, 6 | equal |
| MSG 156.50, K, single, hexagonal fractional coordinates | 1 | 3 | 1, 0, 2 | equal |
| MSG 211.57, Γ, single corep 3, order-3 C-4 Weyl model | 2 | 4 | 1, 0, 2, 1 | equal |

All ten cases also have identical momentum exponent sets and fully resolved
Hamiltonians.  The cases deliberately cover all four MSG types, single- and
double-valued corepresentations, direct sums, multiplicities, an exact
parameterized k line, both momentum-coordinate conventions, a non-orthogonal
hexagonal lattice, and orders one through three.

## Running the comparison

From the Python/Web module directory being tested:

```console
PYTHONPYCACHEPREFIX=/tmp/mkp-pyc pytest -q \
  tests/test_mathematica_catalog_parity.py
```

The test invokes a direct `WolframKernel -noinit -script` process for every
case.  This avoids the standalone `wolframscript` launcher/SharedMemory WSTP
failure mode seen when its output streams are captured by Python.
