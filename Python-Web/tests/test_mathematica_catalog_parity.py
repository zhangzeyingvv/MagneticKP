from __future__ import annotations

import json
import os
import pathlib
import shutil
import subprocess

import pytest
import sympy
from sympy.printing.mathematica import mathematica_code

pytest.importorskip("msgcorep")

from conftest import REPOSITORY_ROOT
from magnetickp import kpHam
from magnetickp.catalog import build_catalog_core_input
from magnetickp.web.parsing import momentum_symbols


STABLE_ROOT = (
    REPOSITORY_ROOT
    if (REPOSITORY_ROOT / "Mathematica").is_dir()
    else REPOSITORY_ROOT / "MagneticKP"
)
MATHEMATICA_PARITY_SCRIPT = (
    STABLE_ROOT
    / "Mathematica"
    / "Tests"
    / "Integration"
    / "PythonCatalogParity.wls"
)


CASES = (
    {
        "name": "type1_B_double_cartesian",
        "space_group": 13,
        "magnetic_group": "13.65",
        "python_k": "B",
        "mma_k": '"B"',
        "family": "double",
        "representations": ({"index": 1},),
        "mma_reps": (2,),
        "order": 1,
        "basis": "cartesian",
    },
    {
        "name": "type1_B_double_fractional",
        "space_group": 13,
        "magnetic_group": "13.65",
        "python_k": "B",
        "mma_k": '"B"',
        "family": "double",
        "representations": ({"index": 1},),
        "mma_reps": (2,),
        "order": 1,
        "basis": "fractional",
    },
    {
        "name": "type2_Gamma_double_cartesian",
        "space_group": 92,
        "magnetic_group": "92.112",
        "python_k": "Γ",
        "mma_k": "{0,0,0}",
        "family": "double",
        "representations": ({"index": 1},),
        "mma_reps": (6,),
        "order": 2,
        "basis": "cartesian",
    },
    {
        "name": "type3_Gamma_double_cartesian",
        "space_group": 92,
        "magnetic_group": "92.113",
        "python_k": "Γ",
        "mma_k": "{0,0,0}",
        "family": "double",
        "representations": ({"index": 1},),
        "mma_reps": (4,),
        "order": 2,
        "basis": "cartesian",
    },
    {
        "name": "type4_Gamma_double_cartesian",
        "space_group": 92,
        "magnetic_group": "92.116",
        "python_k": "Γ",
        "mma_k": "{0,0,0}",
        "family": "double",
        "representations": ({"index": 1},),
        "mma_reps": (6,),
        "order": 2,
        "basis": "cartesian",
    },
    {
        "name": "type4_exact_Lambda_fractional",
        "space_group": 4,
        "magnetic_group": "4.10",
        "python_k": "Λ",
        "mma_k": "{0,0,1/5}",
        "k_parameters": {"u": "1/5"},
        "family": "double",
        "representations": ({"index": 1},),
        "mma_reps": (3,),
        "order": 2,
        "basis": "fractional",
    },
    {
        "name": "single_corep_direct_sum",
        "space_group": 97,
        "magnetic_group": "97.155",
        "python_k": "P",
        "mma_k": '"P"',
        "family": "single",
        "representations": ({"index": 1}, {"index": 2}),
        "mma_reps": (1, 2),
        "order": 1,
        "basis": "cartesian",
    },
    {
        "name": "double_corep_multiplicity_two",
        "space_group": 92,
        "magnetic_group": "92.112",
        "python_k": "Γ",
        "mma_k": "{0,0,0}",
        "family": "double",
        "representations": ({"index": 1, "multiplicity": 2},),
        "mma_reps": (6, 6),
        "order": 1,
        "basis": "cartesian",
    },
    {
        "name": "hexagonal_K_fractional",
        "space_group": 156,
        "magnetic_group": "156.50",
        "python_k": "K",
        "mma_k": '"K"',
        "family": "single",
        "representations": ({"index": 1},),
        "mma_reps": (1,),
        "order": 2,
        "basis": "fractional",
    },
    {
        "name": "C4_Weyl_order_three",
        "space_group": 211,
        "magnetic_group": "211.57",
        "python_k": "Γ",
        "mma_k": "{0,0,0}",
        "family": "single",
        "representations": ({"index": 3},),
        "mma_reps": (3,),
        "order": 3,
        "basis": "cartesian",
    },
)


def _python_case_payload(spec):
    core_input, metadata, group_input = build_catalog_core_input(
        space_group=spec["space_group"],
        magnetic_group=spec["magnetic_group"],
        k_name=spec["python_k"],
        variant="a",
        parameter_values=spec.get("k_parameters"),
        coordinate_basis=spec["basis"],
        representation_family=spec["family"],
        representations=spec["representations"],
        method="iterative",
    )
    result = kpHam(core_input).getkpHam(
        spec["order"],
        method="iterative",
        group_input=group_input,
        return_diagnostics=True,
    )
    hamiltonian = result["ham"]
    parameters = sorted(
        hamiltonian.free_symbols - set(momentum_symbols(3)), key=str
    )
    renamed = hamiltonian.xreplace(
        {
            parameter: sympy.Symbol(f"pyP{index}", real=True)
            for index, parameter in enumerate(parameters, start=1)
        }
    )
    family, number = (int(value) for value in spec["magnetic_group"].split("."))
    return {
        "Name": spec["name"],
        "MSG": [family, number],
        "KExpression": spec["mma_k"],
        "MMAReps": list(spec["mma_reps"]),
        "Order": spec["order"],
        "CartesianCoordinates": spec["basis"] == "cartesian",
        "PythonHamiltonian": mathematica_code(renamed),
        "PythonParameterCount": result["NumberOfParameters"],
        "PythonInvariantDimensions": [
            row["invariant_dimension"] for row in result["Diagnostics"]
        ],
        "PythonConstraintCount": metadata["constraint_count"],
        "PythonMatrixDimension": metadata["dimension"],
    }


def _extract_report(output: str):
    start_marker = "MAGNETICKP_PYTHON_PARITY_BEGIN"
    end_marker = "MAGNETICKP_PYTHON_PARITY_END"
    start = output.rfind(start_marker)
    end = output.rfind(end_marker)
    if start < 0 or end < start:
        return None
    return json.loads(output[start + len(start_marker) : end].strip())


def _wolfram_kernel():
    configured = os.environ.get("WOLFRAM_KERNEL")
    candidates = (
        configured,
        "/Applications/Mathematica.app/Contents/MacOS/WolframKernel",
        shutil.which("WolframKernel"),
    )
    return next(
        (candidate for candidate in candidates if candidate and pathlib.Path(candidate).is_file()),
        None,
    )


@pytest.mark.parametrize("spec", CASES, ids=lambda case: case["name"])
def test_python_catalog_hamiltonian_spaces_match_mathematica_interface(
    tmp_path, spec
):
    kernel = _wolfram_kernel()
    if kernel is None:
        pytest.skip("a direct WolframKernel executable is unavailable")
    assert MATHEMATICA_PARITY_SCRIPT.is_file()
    payload_path = pathlib.Path(tmp_path) / "python_catalog_parity.json"
    payload_path.write_text(
        json.dumps(
            {"SchemaVersion": 1, "Cases": [_python_case_payload(spec)]},
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )
    environment = dict(os.environ)
    environment["MAGNETICKP_PARITY_PAYLOAD"] = str(payload_path)
    process = subprocess.run(
        [kernel, "-noinit", "-script", str(MATHEMATICA_PARITY_SCRIPT)],
        cwd=REPOSITORY_ROOT,
        env=environment,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=90,
        check=False,
    )
    report = _extract_report(process.stdout)
    assert report is not None, (
        f"Mathematica parity report was not produced (exit {process.returncode}).\n"
        f"stdout:\n{process.stdout}\nstderr:\n{process.stderr}"
    )
    failures = {name: row for name, row in report.items() if not row.get("Pass")}
    assert process.returncode == 0 and not failures, (
        f"Mathematica/Python catalog parity failures: {failures}\n"
        f"stderr:\n{process.stderr}"
    )
    assert set(report) == {spec["name"]}
