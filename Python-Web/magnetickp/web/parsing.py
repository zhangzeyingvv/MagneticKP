"""Strict expression parsing for the public MagneticKP Web API."""

from __future__ import annotations

import ast
from typing import Dict, Iterable, Mapping, Sequence

from sympy import E, I, Float, Integer, Matrix, Rational, cos, exp, pi, sin, sqrt, tan
from sympy.core.expr import Expr

from ..transforms import kx, ky, kz


MAX_EXPRESSION_LENGTH = 160
MAX_AST_NODES = 128

_FUNCTIONS = {
    "sqrt": sqrt,
    "exp": exp,
    "sin": sin,
    "cos": cos,
    "tan": tan,
}
_CONSTANTS = {"I": I, "pi": pi, "E": E}


class ExpressionSyntaxError(ValueError):
    """Raised when a Web expression is outside the supported grammar."""


class _ExpressionBuilder(ast.NodeVisitor):
    def __init__(self, symbols: Mapping[str, Expr]):
        self.symbols = dict(symbols)

    def generic_visit(self, node):
        raise ExpressionSyntaxError(
            f"unsupported expression element: {type(node).__name__}"
        )

    def visit_Expression(self, node: ast.Expression):
        return self.visit(node.body)

    def visit_Constant(self, node: ast.Constant):
        value = node.value
        if isinstance(value, bool):
            raise ExpressionSyntaxError("boolean values are not allowed")
        if isinstance(value, int):
            if abs(value) > 10**9:
                raise ExpressionSyntaxError("integer literal is too large")
            return Integer(value)
        if isinstance(value, float):
            return Float(str(value))
        raise ExpressionSyntaxError("only real numeric literals are allowed")

    def visit_Name(self, node: ast.Name):
        if node.id in self.symbols:
            return self.symbols[node.id]
        if node.id in _CONSTANTS:
            return _CONSTANTS[node.id]
        allowed = sorted(set(self.symbols) | set(_CONSTANTS) | set(_FUNCTIONS))
        raise ExpressionSyntaxError(
            f"name {node.id!r} is not allowed; use one of {', '.join(allowed)}"
        )

    def visit_UnaryOp(self, node: ast.UnaryOp):
        value = self.visit(node.operand)
        if isinstance(node.op, ast.UAdd):
            return value
        if isinstance(node.op, ast.USub):
            return -value
        raise ExpressionSyntaxError("only unary + and - are allowed")

    def visit_BinOp(self, node: ast.BinOp):
        left = self.visit(node.left)
        right = self.visit(node.right)
        if isinstance(node.op, ast.Add):
            return left + right
        if isinstance(node.op, ast.Sub):
            return left - right
        if isinstance(node.op, ast.Mult):
            return left * right
        if isinstance(node.op, ast.Div):
            if right == 0:
                raise ExpressionSyntaxError("division by zero")
            return left / right
        if isinstance(node.op, ast.Pow):
            if right.is_Integer and abs(int(right)) > 32:
                raise ExpressionSyntaxError("integer exponent is too large")
            return left**right
        raise ExpressionSyntaxError("only +, -, *, / and powers are allowed")

    def visit_Call(self, node: ast.Call):
        if (
            not isinstance(node.func, ast.Name)
            or node.func.id not in _FUNCTIONS
            or len(node.args) != 1
            or node.keywords
        ):
            raise ExpressionSyntaxError(
                "only one-argument sqrt, exp, sin, cos and tan calls are allowed"
            )
        return _FUNCTIONS[node.func.id](self.visit(node.args[0]))


def parse_expression(text: str, symbols: Mapping[str, Expr] | None = None):
    """Parse a small mathematical grammar without calling eval or sympify."""

    if not isinstance(text, str):
        raise ExpressionSyntaxError("matrix entries must be strings")
    source = text.strip()
    if source == "":
        return Integer(0)
    if len(source) > MAX_EXPRESSION_LENGTH:
        raise ExpressionSyntaxError("expression is too long")
    source = source.replace("^", "**")
    try:
        tree = ast.parse(source, mode="eval")
    except SyntaxError as error:
        raise ExpressionSyntaxError("invalid mathematical expression") from error
    if sum(1 for _ in ast.walk(tree)) > MAX_AST_NODES:
        raise ExpressionSyntaxError("expression is too complex")
    return _ExpressionBuilder(symbols or {}).visit(tree)


def parse_square_matrix(
    rows: Sequence[Sequence[str]],
    *,
    symbols: Mapping[str, Expr] | None = None,
    maximum_dimension: int = 12,
):
    """Parse and validate one non-empty square representation matrix."""

    if not rows or len(rows) > maximum_dimension:
        raise ValueError(
            f"representation dimension must be between 1 and {maximum_dimension}"
        )
    dimension = len(rows)
    if any(len(row) != dimension for row in rows):
        raise ValueError("representation matrix must be square")
    return Matrix(
        [
            [parse_expression(entry, symbols) for entry in row]
            for row in rows
        ]
    )


def momentum_symbols(spatial_dimension: int):
    if spatial_dimension == 1:
        return (kz,)
    if spatial_dimension == 2:
        return (kx, ky)
    if spatial_dimension == 3:
        return (kx, ky, kz)
    raise ValueError("spatial dimension must be 1, 2 or 3")


def parse_momentum_image(entries: Sequence[str], spatial_dimension: int):
    variables = momentum_symbols(spatial_dimension)
    if len(entries) != spatial_dimension:
        raise ValueError(
            "momentum image length must equal the spatial dimension"
        )
    symbols = {str(variable): variable for variable in variables}
    parsed = [parse_expression(entry, symbols) for entry in entries]
    allowed = set(variables)
    if any(not value.free_symbols.issubset(allowed) for value in parsed):
        raise ValueError("momentum images may contain only kx, ky and kz")
    return tuple(parsed)


def parse_representation(
    rows: Sequence[Sequence[str]],
    *,
    continuous_parameter=None,
):
    symbols: Dict[str, Expr] = {}
    if continuous_parameter is not None:
        symbols[str(continuous_parameter)] = continuous_parameter
    matrix = parse_square_matrix(rows, symbols=symbols)
    allowed = set(symbols.values())
    if any(not entry.free_symbols.issubset(allowed) for entry in matrix):
        if continuous_parameter is None:
            raise ValueError("discrete representation matrices must be parameter-free")
        raise ValueError(
            "CInfinity representation may contain only its declared parameter"
        )
    return matrix
