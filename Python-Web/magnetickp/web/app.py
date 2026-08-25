"""FastAPI application serving MagneticKP and its browser interface."""

from __future__ import annotations

import asyncio
import multiprocessing
import os
import queue
import re
from pathlib import Path
from typing import Dict, List, Literal, Optional

import sympy
from fastapi import FastAPI, HTTPException, Request
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field
from starlette.concurrency import run_in_threadpool
from sympy import E, I, Add, Mul, Pow, Symbol, cos, exp, latex, mathml, pi, sin, sstr, tan

from .. import __version__
from ..catalog import (
    CatalogDependencyError,
    CatalogInputError,
    build_catalog_core_input,
    catalog_status,
    corepresentation_catalog,
    k_point_catalog,
    list_bz_variants,
    list_magnetic_groups,
    list_space_groups,
)
from ..magkp import kpHam
from ..kpath import HighSymmetryPathError, format_band_path, high_symmetry_k_path
from .parsing import (
    ExpressionSyntaxError,
    momentum_symbols,
    parse_momentum_image,
    parse_representation,
)


STATIC_DIRECTORY = Path(__file__).with_name("static")
DEFAULT_TIMEOUT = float(os.environ.get("MAGNETICKP_WEB_TIMEOUT", "30"))
MAX_CONCURRENT_JOBS = max(
    1,
    int(os.environ.get("MAGNETICKP_WEB_CONCURRENCY", "2")),
)

_PARAMETER_ORDER = {"e": 0, "t": 1, "r": 2, "s": 3}


def parameter_sort_key(parameter):
    """Sort coefficient symbols by momentum order and numeric index."""

    text = str(parameter)
    low_order = re.fullmatch(r"([etrs])_(\d+)", text)
    if low_order:
        return (_PARAMETER_ORDER[low_order.group(1)], int(low_order.group(2)), text)
    legacy = re.fullmatch(r"C_(\d+)_(\d+)", text)
    if legacy:
        return (int(legacy.group(2)), int(legacy.group(1)), text)
    return (10**9, 10**9, text)


def normalize_base_path(value: str) -> str:
    base_path = value.strip()
    if not base_path or base_path == "/":
        return ""
    if not base_path.startswith("/"):
        raise RuntimeError("MAGNETICKP_BASE_PATH must start with '/'")
    if "?" in base_path or "#" in base_path or "//" in base_path:
        raise RuntimeError("MAGNETICKP_BASE_PATH is not a valid URL path")
    segments = base_path.strip("/").split("/")
    if any(segment in {"", ".", ".."} for segment in segments):
        raise RuntimeError("MAGNETICKP_BASE_PATH is not a valid URL path")
    return "/" + "/".join(segments)


BASE_PATH = normalize_base_path(os.environ.get("MAGNETICKP_BASE_PATH", ""))


class OperationRequest(BaseModel):
    kind: Literal[
        "unitary",
        "antiunitary",
        "antisymmetry_unitary",
        "antisymmetry_antiunitary",
        "continuous",
    ]
    label: str = Field(min_length=1, max_length=64)
    representation: List[List[str]] = Field(min_length=1, max_length=12)
    momentum_image: List[str] = Field(default_factory=list, max_length=3)
    parameter: Optional[Literal["theta"]] = None


class SolveRequest(BaseModel):
    order: int = Field(default=1, ge=0, le=6)
    spatial_dimension: Literal[1, 2, 3] = 3
    method: Literal[
        "automatic",
        "iterative",
        "direct_product",
        "reynolds",
    ] = "automatic"
    group_input: Literal["generators", "elements"] = "generators"
    max_group_order: int = Field(default=512, ge=1, le=4096)
    operations: List[OperationRequest] = Field(min_length=1, max_length=32)


class CatalogRepresentationRequest(BaseModel):
    index: int = Field(ge=1)
    multiplicity: int = Field(default=1, ge=1, le=8)


class CatalogSolveRequest(BaseModel):
    space_group: int = Field(ge=1, le=230)
    magnetic_group: str = Field(min_length=3, max_length=16)
    bz_variant: str = Field(default="a", min_length=1, max_length=8)
    k_name: str = Field(min_length=1, max_length=32)
    k_parameters: Dict[str, str] = Field(default_factory=dict)
    coordinate_basis: Literal["cartesian", "fractional"] = "cartesian"
    representation_family: Literal["single", "double"]
    representations: List[CatalogRepresentationRequest] = Field(
        min_length=1, max_length=12
    )
    order: int = Field(default=1, ge=0, le=6)
    method: Literal["automatic", "iterative", "reynolds"] = "automatic"
    max_group_order: int = Field(default=512, ge=1, le=4096)


_CATEGORY_NAMES = {
    "unitary": "Unitary",
    "antiunitary": "Antiunitary",
    "antisymmetry_unitary": "AntisymmetryUnitaryTest",
    "antisymmetry_antiunitary": "AntisymmetryAntiunitaryTest",
}


def _payload_dict(payload: SolveRequest):
    if hasattr(payload, "model_dump"):
        return payload.model_dump()
    return payload.dict()


def build_core_input(payload: SolveRequest):
    """Convert a validated Web request into the MagneticKP Python API."""

    result: Dict[str, Dict[str, tuple]] = {}
    representation_dimension = None
    continuous_count = 0
    variables = momentum_symbols(payload.spatial_dimension)

    for operation in payload.operations:
        if operation.kind == "continuous":
            continuous_count += 1
            theta = Symbol(operation.parameter or "theta", real=True)
            representation = parse_representation(
                operation.representation,
                continuous_parameter=theta,
            )
            category = "Unitary"
            label = "CInfinity"
            value = (representation, variables, theta)
        else:
            representation = parse_representation(operation.representation)
            category = _CATEGORY_NAMES[operation.kind]
            label = operation.label
            momentum_image = parse_momentum_image(
                operation.momentum_image,
                payload.spatial_dimension,
            )
            value = (representation, momentum_image)

        if representation_dimension is None:
            representation_dimension = representation.rows
        elif representation.rows != representation_dimension:
            raise ValueError(
                "all representation matrices must have the same dimension"
            )
        category_entries = result.setdefault(category, {})
        if label in category_entries:
            raise ValueError(f"duplicate operation label {label!r}")
        category_entries[label] = value

    if continuous_count > 1:
        raise ValueError("only one CInfinity operation is supported")
    return result


def _plain_value(value):
    if value is None or isinstance(value, (bool, int, float, str)):
        return value
    if isinstance(value, dict):
        return {str(key): _plain_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain_value(item) for item in value]
    if getattr(value, "is_Integer", False):
        return int(value)
    return sstr(value)


def _expression_tree(expression):
    """Return a JSON expression tree for safe in-browser numerical plots."""

    if expression == I:
        return {"type": "complex", "real": 0.0, "imaginary": 1.0}
    if expression == pi:
        return {"type": "number", "value": float(pi)}
    if expression == E:
        return {"type": "number", "value": float(E)}
    if isinstance(expression, Symbol):
        return {"type": "symbol", "name": str(expression)}
    if expression.is_Number:
        numeric = complex(expression.evalf())
        if abs(numeric.imag) < 1e-14:
            return {"type": "number", "value": numeric.real}
        return {
            "type": "complex",
            "real": numeric.real,
            "imaginary": numeric.imag,
        }
    if isinstance(expression, Add):
        return {
            "type": "add",
            "arguments": [_expression_tree(item) for item in expression.args],
        }
    if isinstance(expression, Mul):
        return {
            "type": "multiply",
            "arguments": [_expression_tree(item) for item in expression.args],
        }
    if isinstance(expression, Pow):
        return {
            "type": "power",
            "base": _expression_tree(expression.base),
            "exponent": _expression_tree(expression.exp),
        }
    function_names = {exp: "exp", sin: "sin", cos: "cos", tan: "tan"}
    if expression.func in function_names and len(expression.args) == 1:
        return {
            "type": function_names[expression.func],
            "argument": _expression_tree(expression.args[0]),
        }
    raise ValueError(
        f"Hamiltonian expression {sstr(expression)!r} cannot be plotted in the browser"
    )


def _high_symmetry_path_payload(core_input, spatial_dimension):
    """Return compact path data suitable for the browser result payload."""

    if spatial_dimension != 3:
        return {
            "available": False,
            "detail": "automatic high-symmetry paths require three momentum dimensions",
        }
    try:
        data = high_symmetry_k_path(core_input, return_data=True)
    except HighSymmetryPathError as error:
        return {"available": False, "detail": str(error)}
    return {
        "available": True,
        "path": _plain_value(data["path"]),
        "text": format_band_path(data["path"]),
        "directions": _plain_value(data["directions"]),
        "line_data": [
            {
                "direction": _plain_value(line["direction"]),
                "subgroup_order": line["subgroup_order"],
                "momentum_action_order": line["momentum_action_order"],
            }
            for line in data["line_data"]
        ],
        "little_cogroup_order": data["little_cogroup_order"],
        "little_cogroup_action_order": data["little_cogroup_action_order"],
        "enumerated_subgroup_count": data["enumerated_subgroup_count"],
        "center_label": data["center_label"],
        "radius": _plain_value(data["radius"]),
        "coordinate_basis": data["coordinate_basis"],
    }


def compute_core_input(
    core_input,
    *,
    order,
    method,
    group_input,
    max_group_order,
    spatial_dimension=3,
    selection=None,
):
    """Execute one already validated MagneticKP input dictionary."""

    model = kpHam(core_input)
    result = model.getkpHam(
        order,
        method=method,
        group_input=group_input,
        max_group_order=max_group_order,
        return_diagnostics=True,
    )
    hamiltonian = result["ham"]
    parameters = sorted(
        hamiltonian.free_symbols - set(momentum_symbols(spatial_dimension)),
        key=parameter_sort_key,
    )
    response = {
        "status": "ok",
        "dimension": result["dim"],
        "orders": result["korder"],
        "number_of_parameters": result["NumberOfParameters"],
        "parameters": [str(parameter) for parameter in parameters],
        "hamiltonian": {
            "text": sstr(hamiltonian),
            "pretty": sympy.pretty(hamiltonian, use_unicode=True),
            "latex": latex(hamiltonian),
            "mathml": mathml(hamiltonian, printer="presentation"),
            "entries": [
                [sstr(hamiltonian[row, column]) for column in range(hamiltonian.cols)]
                for row in range(hamiltonian.rows)
            ],
            "expression_tree": [
                [
                    _expression_tree(hamiltonian[row, column])
                    for column in range(hamiltonian.cols)
                ]
                for row in range(hamiltonian.rows)
            ],
        },
        "diagnostics": _plain_value(result["Diagnostics"]),
        "high_symmetry_path": _high_symmetry_path_payload(
            core_input,
            spatial_dimension,
        ),
    }
    if selection is not None:
        response["selection"] = _plain_value(selection)
    return response


def compute_request(payload: SolveRequest):
    """Execute a manual matrix request synchronously."""

    return compute_core_input(
        build_core_input(payload),
        order=payload.order,
        method=payload.method,
        group_input=payload.group_input,
        max_group_order=payload.max_group_order,
        spatial_dimension=payload.spatial_dimension,
    )


def compute_catalog_request(payload: CatalogSolveRequest):
    """Resolve one catalog selection and execute it synchronously."""

    representation_rows = [
        (
            item.model_dump()
            if hasattr(item, "model_dump")
            else item.dict()
        )
        for item in payload.representations
    ]
    core_input, metadata, group_input = build_catalog_core_input(
        space_group=payload.space_group,
        magnetic_group=payload.magnetic_group,
        k_name=payload.k_name,
        variant=payload.bz_variant,
        parameter_values=payload.k_parameters,
        coordinate_basis=payload.coordinate_basis,
        representation_family=payload.representation_family,
        representations=representation_rows,
        method=payload.method,
    )
    return compute_core_input(
        core_input,
        order=payload.order,
        method=payload.method,
        group_input=group_input,
        max_group_order=payload.max_group_order,
        spatial_dimension=3,
        selection=metadata,
    )


def _worker(payload_data, output_queue, request_kind="manual"):
    try:
        if request_kind == "catalog":
            payload = CatalogSolveRequest(**payload_data)
            result = compute_catalog_request(payload)
        else:
            payload = SolveRequest(**payload_data)
            result = compute_request(payload)
        output_queue.put(("ok", result))
    except Exception as error:  # pragma: no cover - exercised through parent
        output_queue.put(("error", type(error).__name__, str(error)))


def run_isolated(payload_data, timeout=DEFAULT_TIMEOUT, request_kind="manual"):
    """Run SymPy in a disposable process so timeouts actually stop work."""

    context = multiprocessing.get_context("spawn")
    output_queue = context.Queue(maxsize=1)
    process = context.Process(
        target=_worker, args=(payload_data, output_queue, request_kind)
    )
    process.start()
    process.join(timeout)
    if process.is_alive():
        process.terminate()
        process.join(2)
        if process.is_alive():
            process.kill()
            process.join()
        raise TimeoutError(f"calculation exceeded the {timeout:g} second limit")
    try:
        message = output_queue.get(timeout=1)
    except queue.Empty as error:
        raise RuntimeError(
            f"calculation process exited with code {process.exitcode}"
        ) from error
    finally:
        output_queue.close()
    if message[0] == "error":
        error_name, detail = message[1], message[2]
        error_types = {
            "CatalogDependencyError": CatalogDependencyError,
            "CatalogInputError": CatalogInputError,
            "ExpressionSyntaxError": ExpressionSyntaxError,
            "ValueError": ValueError,
        }
        raise error_types.get(error_name, RuntimeError)(detail)
    return message[1]


def run_catalog_isolated(payload_data, timeout=DEFAULT_TIMEOUT):
    return run_isolated(payload_data, timeout=timeout, request_kind="catalog")


app = FastAPI(
    title="MagneticKP",
    description="Construct symmetry-constrained k.p Hamiltonians",
    version=__version__,
    root_path=BASE_PATH,
    docs_url="/api/docs",
    redoc_url=None,
)
app.mount("/static", StaticFiles(directory=STATIC_DIRECTORY), name="static")
app.mount(
    "/vendor",
    StaticFiles(directory=STATIC_DIRECTORY / "vendor"),
    name="vendor",
)
app.mount(
    "/help",
    StaticFiles(directory=STATIC_DIRECTORY / "help", html=True),
    name="help",
)
_job_slots = asyncio.Semaphore(MAX_CONCURRENT_JOBS)


@app.middleware("http")
async def security_headers(request: Request, call_next):
    response = await call_next(request)
    response.headers["X-Content-Type-Options"] = "nosniff"
    response.headers["X-Frame-Options"] = "DENY"
    response.headers["Referrer-Policy"] = "no-referrer"
    response.headers["Permissions-Policy"] = (
        "camera=(), microphone=(), geolocation=()"
    )
    return response


@app.api_route("/", methods=["GET", "HEAD"], include_in_schema=False)
async def index():
    return FileResponse(
        STATIC_DIRECTORY / "index.html",
        headers={"Cache-Control": "no-store"},
    )


@app.api_route("/styles.css", methods=["GET", "HEAD"], include_in_schema=False)
async def root_stylesheet():
    """Serve the stylesheet for both direct and server-hosted index paths."""

    return FileResponse(STATIC_DIRECTORY / "styles.css", media_type="text/css")


@app.api_route("/app.js", methods=["GET", "HEAD"], include_in_schema=False)
async def root_script():
    """Serve the browser application beside the root index document."""

    return FileResponse(
        STATIC_DIRECTORY / "app.js",
        media_type="text/javascript",
    )


@app.get("/api/health")
async def health():
    catalog = catalog_status()
    return {
        "status": "ok",
        "magnetickp_version": __version__,
        "sympy_version": sympy.__version__,
        "timeout_seconds": DEFAULT_TIMEOUT,
        "msgcorep": catalog,
    }


@app.get("/api/catalog/space-groups")
async def catalog_space_groups():
    try:
        groups = await run_in_threadpool(list_space_groups)
        return {"count": len(groups), "groups": groups}
    except CatalogDependencyError as error:
        raise HTTPException(status_code=503, detail=str(error)) from error


@app.get("/api/catalog/magnetic-groups/{space_group}")
async def catalog_magnetic_groups(space_group: int):
    try:
        groups = await run_in_threadpool(list_magnetic_groups, space_group)
        return {"space_group": space_group, "count": len(groups), "groups": groups}
    except CatalogDependencyError as error:
        raise HTTPException(status_code=503, detail=str(error)) from error
    except (CatalogInputError, ValueError) as error:
        raise HTTPException(status_code=400, detail=str(error)) from error


@app.get("/api/catalog/bz-variants/{space_group}")
async def catalog_bz_variants(space_group: int):
    try:
        variants = await run_in_threadpool(list_bz_variants, space_group)
        return {"space_group": space_group, "variants": variants}
    except CatalogDependencyError as error:
        raise HTTPException(status_code=503, detail=str(error)) from error
    except (CatalogInputError, ValueError) as error:
        raise HTTPException(status_code=400, detail=str(error)) from error


@app.get("/api/catalog/k-points/{space_group}")
async def catalog_k_points(
    space_group: int,
    magnetic_group: Optional[str] = None,
    variant: str = "a",
):
    try:
        return await run_in_threadpool(
            k_point_catalog,
            space_group,
            magnetic_group=magnetic_group,
            variant=variant,
        )
    except CatalogDependencyError as error:
        raise HTTPException(status_code=503, detail=str(error)) from error
    except (CatalogInputError, ValueError) as error:
        raise HTTPException(status_code=400, detail=str(error)) from error


@app.get("/api/catalog/corepresentations")
async def catalog_corepresentations(
    magnetic_group: str,
    k_name: str,
    variant: str = "a",
    coordinate_basis: Literal["cartesian", "fractional"] = "cartesian",
):
    try:
        return await run_in_threadpool(
            corepresentation_catalog,
            magnetic_group,
            k_name,
            variant=variant,
            coordinate_basis=coordinate_basis,
        )
    except CatalogDependencyError as error:
        raise HTTPException(status_code=503, detail=str(error)) from error
    except (CatalogInputError, ValueError) as error:
        raise HTTPException(status_code=400, detail=str(error)) from error


@app.post("/api/solve")
async def solve(payload: SolveRequest):
    try:
        async with _job_slots:
            return await run_in_threadpool(
                run_isolated,
                _payload_dict(payload),
                DEFAULT_TIMEOUT,
            )
    except TimeoutError as error:
        raise HTTPException(status_code=408, detail=str(error)) from error
    except (ExpressionSyntaxError, ValueError) as error:
        raise HTTPException(status_code=400, detail=str(error)) from error
    except Exception as error:
        raise HTTPException(
            status_code=500,
            detail="the calculation process failed unexpectedly",
        ) from error


@app.post("/api/catalog/solve")
async def solve_catalog(payload: CatalogSolveRequest):
    try:
        async with _job_slots:
            return await run_in_threadpool(
                run_catalog_isolated,
                _payload_dict(payload),
                DEFAULT_TIMEOUT,
            )
    except TimeoutError as error:
        raise HTTPException(status_code=408, detail=str(error)) from error
    except CatalogDependencyError as error:
        raise HTTPException(status_code=503, detail=str(error)) from error
    except (CatalogInputError, ValueError) as error:
        raise HTTPException(status_code=400, detail=str(error)) from error
    except Exception as error:
        raise HTTPException(
            status_code=500,
            detail="the catalog calculation process failed unexpectedly",
        ) from error
