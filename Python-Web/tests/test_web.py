from __future__ import annotations

import pytest
from sympy import I, Matrix, Symbol, exp, eye, simplify, sqrt, zeros

pytest.importorskip("fastapi")

from magnetickp.web.app import (
    STATIC_DIRECTORY,
    SolveRequest,
    app,
    build_core_input,
    compute_request,
    normalize_base_path,
    run_isolated,
)
from magnetickp.web.parsing import ExpressionSyntaxError, parse_expression


def continuous_payload(order=0):
    return SolveRequest(
        order=order,
        spatial_dimension=3,
        method="automatic",
        operations=[
            {
                "kind": "continuous",
                "label": "CInfinity",
                "parameter": "theta",
                "representation": [
                    ["exp(-I*theta/2)", "0"],
                    ["0", "exp(I*theta/2)"],
                ],
                "momentum_image": ["kx", "ky", "kz"],
            }
        ],
    )


def test_direct_and_semidirect_continuous_examples_obey_their_group_relations():
    theta = Symbol("theta", real=True)
    direct_continuous = Matrix.diag(
        exp(I * theta / 2),
        exp(I * theta / 2),
        exp(I * theta / 2),
        exp(I * theta / 2),
    )
    semidirect_continuous = Matrix.diag(
        exp(I * theta / 2),
        exp(I * theta / 2),
        exp(I * theta / 2),
        exp(I * theta / 2),
    )
    independent_c2 = Matrix(
        [
            [0, 1, 0, 0],
            [1, 0, 0, 0],
            [0, 0, 0, 1],
            [0, 0, 1, 0],
        ]
    )
    antiunitary_c2 = eye(4)

    assert independent_c2**2 == eye(4)
    assert antiunitary_c2 * antiunitary_c2.conjugate() == eye(4)
    direct_commutator = independent_c2 * direct_continuous * independent_c2 - direct_continuous
    assert direct_commutator.applyfunc(simplify) == zeros(4)
    reversed_continuous = semidirect_continuous.subs(theta, -theta)
    flipped_continuous = antiunitary_c2 * semidirect_continuous.conjugate() * antiunitary_c2
    assert (flipped_continuous - reversed_continuous).applyfunc(simplify) == zeros(4)


def test_safe_expression_parser_accepts_only_the_web_grammar():
    theta = Symbol("theta", real=True)
    assert parse_expression("sqrt(2)/2 + I/3") == sqrt(2) / 2 + I / 3
    assert parse_expression("exp(-I*theta/2)", {"theta": theta}).has(theta)

    for malicious in (
        "__import__('os').system('id')",
        "(1).__class__",
        "open('/tmp/x')",
        "[x for x in (1, 2)]",
    ):
        with pytest.raises(ExpressionSyntaxError):
            parse_expression(malicious, {"theta": theta})


def test_web_payload_builds_the_continuous_core_input():
    core_input = build_core_input(continuous_payload())
    representation, momentum, theta = core_input["Unitary"]["CInfinity"]
    assert representation.shape == (2, 2)
    assert [str(value) for value in momentum] == ["kx", "ky", "kz"]
    assert str(theta) == "theta"


def test_web_computation_returns_json_math_and_plot_tree():
    result = compute_request(continuous_payload())
    assert result["status"] == "ok"
    assert result["dimension"] == 2
    assert result["number_of_parameters"] == 2
    assert result["parameters"] == ["e_1", "e_2"]
    assert result["orders"] == [0]
    assert len(result["hamiltonian"]["expression_tree"]) == 2
    assert "mtable" in result["hamiltonian"]["mathml"]
    assert result["diagnostics"][0]["continuous_symmetry"] == "CInfinity"
    assert result["high_symmetry_path"]["available"] is False


def test_disposable_worker_matches_direct_computation():
    payload = continuous_payload()
    payload_data = payload.model_dump() if hasattr(payload, "model_dump") else payload.dict()
    result = run_isolated(payload_data, timeout=20)
    assert result["number_of_parameters"] == 2
    assert result["hamiltonian"]["expression_tree"]


@pytest.mark.parametrize(
    "payload",
    [
        {
            "order": 1,
            "spatial_dimension": 3,
            "operations": [
                {
                    "kind": "unitary",
                    "label": "C4z",
                    "representation": [
                        ["sqrt(2)/2", "sqrt(2)/2"],
                        ["-sqrt(2)/2", "sqrt(2)/2"],
                    ],
                    "momentum_image": ["-ky", "kx", "kz"],
                },
                {
                    "kind": "unitary",
                    "label": "C2x",
                    "representation": [["0", "I"], ["I", "0"]],
                    "momentum_image": ["kx", "-ky", "-kz"],
                },
                {
                    "kind": "antiunitary",
                    "label": "T",
                    "representation": [["0", "-1"], ["1", "0"]],
                    "momentum_image": ["-kx", "-ky", "-kz"],
                },
            ],
        },
        {
            "order": 2,
            "spatial_dimension": 3,
            "operations": [
                {
                    "kind": "antiunitary",
                    "label": "PT",
                    "representation": [["1", "0"], ["0", "1"]],
                    "momentum_image": ["kx", "ky", "kz"],
                }
            ],
        },
        {
            "order": 1,
            "spatial_dimension": 3,
            "operations": [
                {
                    "kind": "continuous",
                    "label": "CInfinity",
                    "parameter": "theta",
                    "representation": [
                        ["exp(I*theta/2)", "0", "0", "0"],
                        ["0", "exp(I*theta/2)", "0", "0"],
                        ["0", "0", "exp(I*theta/2)", "0"],
                        ["0", "0", "0", "exp(I*theta/2)"],
                    ],
                    "momentum_image": ["kx", "ky", "kz"],
                },
                {
                    "kind": "unitary",
                    "label": "IndependentC2",
                    "representation": [
                        ["0", "1", "0", "0"],
                        ["1", "0", "0", "0"],
                        ["0", "0", "0", "1"],
                        ["0", "0", "1", "0"],
                    ],
                    "momentum_image": ["kx", "ky", "kz"],
                },
            ],
        },
        {
            "order": 1,
            "spatial_dimension": 3,
            "operations": [
                {
                    "kind": "continuous",
                    "label": "CInfinity",
                    "parameter": "theta",
                    "representation": [
                        ["exp(I*theta/2)", "0", "0", "0"],
                        ["0", "exp(I*theta/2)", "0", "0"],
                        ["0", "0", "exp(I*theta/2)", "0"],
                        ["0", "0", "0", "exp(I*theta/2)"],
                    ],
                    "momentum_image": ["kx", "ky", "kz"],
                },
                {
                    "kind": "antiunitary",
                    "label": "ThetaFlip",
                    "representation": [
                        ["1", "0", "0", "0"],
                        ["0", "1", "0", "0"],
                        ["0", "0", "1", "0"],
                        ["0", "0", "0", "1"],
                    ],
                    "momentum_image": ["kx", "ky", "kz"],
                },
            ],
        },
        {
            "order": 3,
            "spatial_dimension": 3,
            "operations": [
                {
                    "kind": "unitary",
                    "label": "C3m",
                    "representation": [
                        ["-1/2", "sqrt(3)/2"],
                        ["-sqrt(3)/2", "-1/2"],
                    ],
                    "momentum_image": ["ky", "kz", "kx"],
                },
                {
                    "kind": "unitary",
                    "label": "C2z",
                    "representation": [["1", "0"], ["0", "1"]],
                    "momentum_image": ["-kx", "-ky", "kz"],
                },
                {
                    "kind": "unitary",
                    "label": "C2x",
                    "representation": [["1", "0"], ["0", "1"]],
                    "momentum_image": ["kx", "-ky", "-kz"],
                },
                {
                    "kind": "unitary",
                    "label": "C2a",
                    "representation": [["1", "0"], ["0", "-1"]],
                    "momentum_image": ["ky", "kx", "-kz"],
                },
                {
                    "kind": "antiunitary",
                    "label": "T",
                    "representation": [["1", "0"], ["0", "1"]],
                    "momentum_image": ["-kx", "-ky", "-kz"],
                },
            ],
        },
        {
            "order": 1,
            "spatial_dimension": 3,
            "operations": [
                {
                    "kind": "antiunitary",
                    "label": "T",
                    "representation": [["0", "1"], ["-1", "0"]],
                    "momentum_image": ["-kx", "-ky", "-kz"],
                }
            ],
        },
    ],
)
def test_verified_browser_examples_compute_and_serialize(payload):
    result = compute_request(SolveRequest(**payload))
    assert result["status"] == "ok"
    assert result["number_of_parameters"] > 0
    assert result["hamiltonian"]["expression_tree"]


def test_topological_examples_have_the_expected_constraints():
    # Yu et al., Table S5 and Eq. (S75): type-II MSG 92 Gamma R6 C-1 WP.
    c1_result = compute_request(
        SolveRequest(
            order=1,
            spatial_dimension=3,
            operations=[
                {
                    "kind": "unitary",
                    "label": "C4z",
                    "representation": [
                        ["sqrt(2)/2", "sqrt(2)/2"],
                        ["-sqrt(2)/2", "sqrt(2)/2"],
                    ],
                    "momentum_image": ["-ky", "kx", "kz"],
                },
                {
                    "kind": "unitary",
                    "label": "C2x",
                    "representation": [["0", "I"], ["I", "0"]],
                    "momentum_image": ["kx", "-ky", "-kz"],
                },
                {
                    "kind": "antiunitary",
                    "label": "T",
                    "representation": [["0", "-1"], ["1", "0"]],
                    "momentum_image": ["-kx", "-ky", "-kz"],
                },
            ],
        )
    )
    assert c1_result["number_of_parameters"] == 3
    assert c1_result["parameters"] == ["e_1", "t_1", "t_2"]
    assert [item["invariant_dimension"] for item in c1_result["diagnostics"]] == [1, 2]
    c1_entries = " ".join(sum(c1_result["hamiltonian"]["entries"], []))
    assert "ky*t_1" in c1_entries
    assert "kx*t_1" in c1_entries
    assert "I*kz*t_2" in c1_entries

    # Spinless PT=K protects a real two-band Weyl nodal line with Berry phase pi.
    pt_result = compute_request(
        SolveRequest(
            order=2,
            spatial_dimension=3,
            operations=[
                {
                    "kind": "antiunitary",
                    "label": "PT",
                    "representation": [["1", "0"], ["0", "1"]],
                    "momentum_image": ["kx", "ky", "kz"],
                }
            ],
        )
    )
    assert pt_result["number_of_parameters"] == 30
    assert all(
        "I" not in entry
        for row in pt_result["hamiltonian"]["entries"]
        for entry in row
    )

    # Yu et al., Table II and Eq. (S5): SG 211 Gamma-point R3 C-4 WP.
    c4_result = compute_request(
        SolveRequest(
            order=3,
            spatial_dimension=3,
            operations=[
                {
                    "kind": "unitary",
                    "label": "C3m",
                    "representation": [
                        ["-1/2", "sqrt(3)/2"],
                        ["-sqrt(3)/2", "-1/2"],
                    ],
                    "momentum_image": ["ky", "kz", "kx"],
                },
                {
                    "kind": "unitary",
                    "label": "C2z",
                    "representation": [["1", "0"], ["0", "1"]],
                    "momentum_image": ["-kx", "-ky", "kz"],
                },
                {
                    "kind": "unitary",
                    "label": "C2x",
                    "representation": [["1", "0"], ["0", "1"]],
                    "momentum_image": ["kx", "-ky", "-kz"],
                },
                {
                    "kind": "unitary",
                    "label": "C2a",
                    "representation": [["1", "0"], ["0", "-1"]],
                    "momentum_image": ["ky", "kx", "-kz"],
                },
                {
                    "kind": "antiunitary",
                    "label": "T",
                    "representation": [["1", "0"], ["0", "1"]],
                    "momentum_image": ["-kx", "-ky", "-kz"],
                },
            ],
        )
    )
    assert c4_result["number_of_parameters"] == 4
    assert c4_result["parameters"] == ["e_1", "r_1", "r_2", "s_1"]
    assert [item["invariant_dimension"] for item in c4_result["diagnostics"]] == [1, 0, 2, 1]
    c4_entries = " ".join(sum(c4_result["hamiltonian"]["entries"], []))
    assert "sqrt(3)" in c4_entries
    assert "kx*ky*kz*s_1" in c4_entries
    assert "t_" not in c4_entries


def test_application_exposes_api_and_hilbertbalance_style_static_ui():
    paths = {route.path for route in app.routes}
    assert {
        "/",
        "/api/health",
        "/api/solve",
        "/api/catalog/space-groups",
        "/api/catalog/corepresentations",
        "/api/catalog/solve",
        "/styles.css",
        "/app.js",
        "/help",
        "/vendor",
        "/static",
    }.issubset(paths)
    page = (STATIC_DIRECTORY / "index.html").read_text(encoding="utf-8")
    script = (STATIC_DIRECTORY / "app.js").read_text(encoding="utf-8")
    style = (STATIC_DIRECTORY / "styles.css").read_text(encoding="utf-8")
    help_directory = STATIC_DIRECTORY / "help"
    help_page = (help_directory / "index.html").read_text(encoding="utf-8")
    help_script = (help_directory / "help.js").read_text(encoding="utf-8")
    help_style = (help_directory / "styles.css").read_text(encoding="utf-8")
    katex_directory = STATIC_DIRECTORY / "vendor" / "katex"
    katex_css = (katex_directory / "katex.min.css").read_text(encoding="utf-8")
    assert 'id="operationList"' in page
    assert 'id="bandPlot"' in page
    assert 'id="bandPathHint"' in page
    assert 'class="help-link" href="help/"' in page
    assert "api/docs" not in page
    assert "direct_product" not in page
    assert "methodDirect" not in script
    assert "Direct product" not in page
    assert 'vendor/katex/katex.min.css?v=0.18.4' in page
    assert 'vendor/katex/katex.min.js?v=0.18.4' in page
    assert page.index("katex.min.js") < page.index('src="app.js')
    assert 'href="styles.css?v=' in page
    assert 'src="app.js?v=' in page
    assert 'window.location.protocol === "file:"' in script
    assert "result.high_symmetry_path" in script
    assert "bandPathGeneratedHint" in script
    assert (katex_directory / "katex.min.js").stat().st_size > 100_000
    assert (katex_directory / "LICENSE").is_file()
    assert (katex_directory / "fonts" / "KaTeX_Main-Regular.woff2").is_file()
    assert "fonts/KaTeX_Main-Regular.woff2" in katex_css
    assert page.index('class="hero-copy"') < page.index('id="solverForm"')
    assert "theory-note" not in page
    assert "feature-strip" not in page
    assert 'inputMode === "catalog" ? "api/catalog/solve" : "api/solve"' in script
    assert 'id="catalogSpaceGroup"' in page
    assert 'id="catalogMagneticGroup"' in page
    assert 'id="catalogKPoint"' in page
    assert 'id="catalogRepresentationFamily"' in page
    assert 'id="catalogCoordinateBasis"' in page
    assert 'id="corepresentationList"' in page
    assert 'id="catalogGeneratorList"' in page
    assert 'id="catalogGeneratorCount"' in page
    assert 'id="catalogRepresentationMatrices"' in page
    assert 'id="convertCatalogToManualButton"' in page
    assert page.index("H(k)") < page.index('id="copyLatexButton"') < page.index('id="hamiltonianMath"')
    assert 'id="formulaViewButton"' in page
    assert 'id="texViewButton"' in page
    assert '<details class="latex-panel">' not in page
    assert 'value="cartesian"' in page
    assert 'value="fractional"' in page
    assert "MMA 默认" not in page
    assert "MMA default" not in script
    assert 'setInputMode("catalog");' in script
    assert "function renderCatalogGenerators" in script
    assert "generator.momentum_image_latex" in script
    assert "function renderRepresentationMatrices" in script
    assert "matrix.matrix_latex" in script
    assert "antiunitaryOperatorHint" not in script
    assert "\\widehat{g}" not in script
    assert "function convertCatalogToManual" in script
    assert "function blockDiagonalMatrix" in script
    assert "function hermitianEigenvalues" in script
    assert "function renderHamiltonian" in script
    assert "function setHamiltonianView" in script
    assert 'setHamiltonianView("formula")' in script
    assert 'setHamiltonianView("tex")' in script
    assert "function renderSubscriptedText" in script
    assert "function unicodeCrystallographicSymbol" in script
    assert "unicodeCrystallographicSymbol(group.symbol_tex, group.symbol)" in script
    assert "unicodeCrystallographicSymbol(group.bns_symbol_tex, group.bns_symbol)" in script
    assert "function renderInlineMath" in script
    assert "function renderSpatialDimensionTabs" in script
    assert "function syncSpatialDimensionTabs" in script
    assert "function displayLatex" in script
    assert '.replace(/\\bkx\\b/g, "k_{x}")' in script
    assert "renderInlineMath(symbol, displayLatex(name), name);" in script
    assert "window.katex.render" in script
    assert 'output: "htmlAndMathml"' in script
    assert 'hamiltonianMath.dataset.renderer = "mathml"' in script
    assert "const initial = 0;" in script
    assert 'addContinuous: "添加 C∞ 操作"' in script
    assert 'id="spatialDimensionTabs"' in page
    assert r'data-latex="k_x,\,k_y,\,k_z"' in page
    assert '<option value="3">3D</option>' in page
    assert 'hilbertDimension: "哈密顿量阶数"' in script
    assert "Hilbert 空间维数" not in page
    assert '<span aria-hidden="true">∞</span>' not in page
    assert 'id="viewToggle"' not in page
    assert 'id="hamiltonianPretty"' not in page
    assert "showText" not in script
    assert "showMath" not in script
    assert 'return "P:-1/2,0,0 -> O:0,0,0; O:0,0,0 -> Q:0,1/2,0";' in script
    assert "applyExample(examples[0], false);" in script
    assert "scroll-margin-top: 24px" in style
    assert "C-1 Weyl 点（SG 92 · Γ）" in script
    assert "C-4 Weyl 点（SG 211 · Γ）" in script
    assert "C 表示 Chern 电荷而非四重旋转" not in script
    assert "无自旋 PT=K" not in script
    assert "localized(example.description)" not in script
    assert 'title: { zh: "连续对称性直积"' in script
    assert 'meta: "C∞ × C₂ · 4 × 4 · ORDER 1"' in script
    assert 'title: { zh: "连续对称性半直积"' in script
    assert 'meta: "C∞ ⋊ C₂ · 4 × 4 · ORDER 1"' in script
    assert "--cyan: #24b8bd" in style
    assert "--blue: #5269e8" in style
    assert "font-size:1.2em" in style
    assert "font-size:17px" in style
    assert '<span class="lang-zh">帮助文档</span>' in help_page
    assert '<span class="lang-zh">打开 Web 版 MagneticKP</span>' in help_page
    assert '<span class="lang-en">Open MagneticKP Web</span>' in help_page
    assert "打开计算器" not in help_page
    assert "教程" in help_page
    assert "指南" in help_page
    assert "参考" in help_page
    assert "MAGNETICKP_WEB_HOST=127.0.0.1" in help_page
    assert "MAGNETICKP_WEB_PORT=18765" in help_page
    assert "dist/magnetickp-1.2.0-py3-none-any.whl[" not in help_page
    assert "先在项目目录创建虚拟环境" in help_page
    assert "python3 -m venv env" in help_page
    assert "source env/bin/activate" in help_page
    assert "python -m pip install --upgrade pip" in help_page
    assert "deactivate" in help_page
    assert "dist/msgcorep-0.1.0a0-py3-none-any.whl" in help_page
    assert "pip 会在这一次安装中自动找到它" in help_page
    assert "import msgcorep" in help_page
    assert "magnetickp.catalog" in help_page
    assert "build_catalog_core_input" in help_page
    assert "自动根据 K 点的小陪群计算高对称线" in help_page
    assert "high_symmetry_k_path(" in help_page
    assert "((0, 0, 1), (1, 1, 0))" in help_page
    assert "L2-:-1/10,-1/10,0" in help_page
    reference_ids = {
        "ref-kpham": "kpHam",
        "ref-getkpham": "GETkpHam",
        "ref-interface": "interfaceRep",
        "ref-high-symmetry": "highSymmetryKPath",
        "ref-band-manipulate": "bandManipulate",
        "ref-bandplot": "bandplot",
        "ref-vasp": "vaspEig",
        "ref-fitting": "fittingKP",
        "ref-generator": "getGenerator",
    }
    for section_id, title in reference_ids.items():
        assert f'id="{section_id}"' in help_page
        assert f'href="#{section_id}"' in help_page
        assert title in help_page
    assert 'id="getting-started"' in help_page
    assert 'id="guide-construction"' in help_page
    assert 'id="guide-bands"' in help_page
    assert 'id="guide-utilities"' in help_page
    assert "G:0,0,0 -> X:1/3,0,0" in help_page
    assert 'magnetic_group="13.65"' in help_page
    assert 'representation_family="double"' in help_page
    assert 'magnetic_group="191.234"' in help_page
    assert 'representation_family="single"' in help_page
    assert "from magnetickp import read_vasp_eigenval" in help_page
    assert "def read_vasp_eigenval" not in help_page
    assert "from magnetickp import fittingKP" in help_page
    assert 'fit = fittingKP(' in help_page
    assert 'ham = kpHam(core_input).getkpHam(2, group_input=group_input)["ham"]' in help_page
    assert 'vasp_file = "EIGENVAL"' in help_page
    assert "kpoint_center=(1/3, 1/3, 0)" in help_page
    assert "hamshift" not in help_page
    assert "def predictions" not in help_page
    assert "scipy.optimize import least_squares" not in help_page
    assert "3.6511886991089995" in help_page
    assert "0.001692172786575689" in help_page
    assert "from magnetickp import bandplot, kpHam" in help_page
    assert "bandplot(ham, path, rules, subdivisions=30)" in help_page
    assert "e_1, t_1, t_2, t_3 = sorted" in help_page
    assert "e_1: 0.2" in help_page
    assert "t_3: 0.8" in help_page
    assert "rules = dict(zip(params, (0.2, 0.4, 0.6, 0.8)))" not in help_page
    assert "plot_band_segments" not in help_page
    assert "def band_energies" not in help_page
    assert "倒易分数坐标并乘以 2π" in help_page
    assert "蓝点：VASP；灰色虚线：默认零初值；红线：拟合结果" in help_page
    assert "从 1 开始的 k 点编号" in help_page
    assert "Blue points: VASP; gray dashed: zero initial values; red lines: fitted result" in help_page
    assert 'src="assets/bandplot-examples.svg"' in help_page
    assert 'src="assets/fittingkp-comparison.svg"' in help_page
    assert "get_generators(range(4)" in help_page
    assert "bandManipulateWithVASP" not in help_page
    assert "ref-band-vasp" not in help_page
    assert "directSum" not in help_page
    assert "ref-direct-sum" not in help_page
    assert "IterativelySimplify" not in help_page
    assert "ref-iterative" not in help_page
    assert "COMPUTED EXAMPLE" not in help_page
    assert "已计算范例" not in help_page
    assert "正确性与结果解释" not in help_page
    assert "\\u0393" not in help_page
    assert 'if row["name"] == "Γ"' not in help_page
    assert 'href="../"' in help_page
    assert '../vendor/katex/katex.min.css?v=0.18.4' in help_page
    assert "from sympy import Matrix, latex" in help_page
    assert 'print(latex(ham))' in help_page
    assert "<code>inputdict</code>" in help_page
    assert '"CInfinity": (D_theta, (kz,), theta)' in help_page
    assert "e_1 + e_2 + kz*t_1 + kz*t_2" in help_page
    assert "例 2：两种反对称输入" not in help_page
    assert '"AntisymmetryAntiunitaryTest"' in help_page
    assert "kx*t_1 + ky*t_2 - I*kz*t_3" in help_page
    assert "return_diagnostics" not in help_page
    assert "navigator.clipboard" in help_script
    assert "function renderCodeComments()" in help_script
    assert 'line.indexOf(" || ")' in help_script
    assert "code.dataset.bilingualSource" in help_script
    assert "IntersectionObserver" in help_script
    assert 'id="languageToggle"' in help_page
    assert "Python + Web Documentation" in help_script
    assert 'localStorage.getItem("magnetickp-language")' in help_script
    assert 'languageToggle.addEventListener("click"' in help_script
    assert 'document.body.dataset.language = currentLanguage' in help_script
    assert ".language-toggle" in help_style
    assert 'body[data-language="en"] .lang-zh' in help_style
    assert "--doc-blue" in help_style


def test_source_distribution_manifest_includes_nested_web_assets():
    manifest = (STATIC_DIRECTORY.parents[2] / "MANIFEST.in").read_text(encoding="utf-8")
    assert "recursive-include magnetickp/web/static *" in manifest
    setup_source = (STATIC_DIRECTORY.parents[2] / "setup.py").read_text(encoding="utf-8")
    assert '"static/help/assets/*"' in setup_source
    assert '"web":' not in setup_source
    assert '"plot":' not in setup_source
    assert '"fit":' not in setup_source
    assert '"catalog":' not in setup_source
    for dependency in (
        "MSGCOREP_DEPENDENCY",
        "FASTAPI_DEPENDENCY",
        "NUMPY_DEPENDENCY",
        "MATPLOTLIB_DEPENDENCY",
        "SCIPY_DEPENDENCY",
        '"uvicorn[standard]>=0.23,<1"',
    ):
        assert dependency in setup_source


def test_root_relative_assets_work_when_index_is_server_hosted():
    testclient = pytest.importorskip("fastapi.testclient")
    with testclient.TestClient(app) as client:
        index = client.get("/")
        assert index.status_code == 200
        assert index.headers["cache-control"] == "no-store"
        stylesheet = client.get("/styles.css")
        script = client.get("/app.js")
        katex = client.get("/vendor/katex/katex.min.css")
        help_index = client.get("/help/")
        help_stylesheet = client.get("/help/styles.css")
        bandplot_figure = client.get("/help/assets/bandplot-examples.svg")
        fitting_figure = client.get("/help/assets/fittingkp-comparison.svg")
    assert stylesheet.status_code == 200
    assert stylesheet.headers["content-type"].startswith("text/css")
    assert script.status_code == 200
    assert "javascript" in script.headers["content-type"]
    assert katex.status_code == 200
    assert help_index.status_code == 200
    assert '<span class="lang-zh">帮助文档</span>' in help_index.text
    assert help_stylesheet.status_code == 200
    assert help_stylesheet.headers["content-type"].startswith("text/css")
    assert bandplot_figure.status_code == 200
    assert "svg" in bandplot_figure.headers["content-type"]
    assert fitting_figure.status_code == 200
    assert "svg" in fitting_figure.headers["content-type"]


def test_base_path_is_reverse_proxy_safe():
    assert normalize_base_path("") == ""
    assert normalize_base_path("/") == ""
    assert normalize_base_path("/apps/magnetickp/") == "/apps/magnetickp"
    for invalid in ("apps/magnetickp", "/apps//magnetickp", "/../x"):
        with pytest.raises(RuntimeError):
            normalize_base_path(invalid)
