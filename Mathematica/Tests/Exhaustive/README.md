# Exhaustive v1.04 comparison

These scripts compare the frozen Mathematica v1.04 implementation in
`../MagneticKP-Mathmatica` with the current package in independent Wolfram
kernels.

The first layer covers every magnetic space group at Gamma and every
single-valued and double-valued corepresentation stored by MSGCorep. Results
are resumable and are written outside the Git repository under
`../tmp/exhaustive_gamma_comparison/`.

From either Mathematica module directory, run a range of the 1651 BNS
magnetic groups:

```console
wolframscript -file Tests/Exhaustive/RunGammaBatch.wls legacy 1 100 1
wolframscript -file Tests/Exhaustive/RunGammaBatch.wls current 1 100 1
wolframscript -file Tests/Exhaustive/CompareGammaResults.wls 1
```

Arguments are the engine, inclusive magnetic-group index range, and maximum
momentum order. Add `force` as the last argument to overwrite completed group
files.

For every corepresentation the comparison records:

- representation dimension;
- total reported parameter count through the requested order;
- momentum exponent support;
- a canonical row basis of the Hamiltonian space;
- pure `kpHam` wall time.

The comparison succeeds only if both engines return resolved Hamiltonians
with the same representation dimension, momentum exponent support, and
Hamiltonian coefficient subspace. The coefficient subspace, rather than the
raw reported parameter count, is the mathematical equivalence criterion.

Reported parameter-count differences are retained as informational metadata
in `summary.wxf` under `ReportedParameterCountDifferences` and
`ReportedParameterCountDifferenceDetails`. They are printed separately and
do not by themselves fail the comparison. This distinction is necessary
because the frozen legacy implementation can over-report a parameter while
returning a rank-deficient parametrization. For example, MSG `{200,14}` at
Gamma, corepresentation index 12, reports two legacy parameters versus one
current parameter at order 1, although both canonical Hamiltonian spaces are
the same one-dimensional span of the identity matrix. The JSON summary keeps
the difference count but intentionally omits the detailed case list.
