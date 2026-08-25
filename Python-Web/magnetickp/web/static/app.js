const translations = {
  zh: {
    pageTitle: "MagneticKP · 对称性约束的 k·p 模型",
    metaDescription: "MagneticKP：从幺正、反幺正和连续自旋对称性构造 k·p 哈密顿量。",
    brandAria: "MagneticKP 首页",
    brandTagline: "SYMMETRY MODEL BUILDER",
    helpLink: "帮助",
    switchLanguage: "Switch to English",
    runtimeConnecting: "正在连接计算引擎",
    runtimeReady: "v{version} · SymPy {sympy}",
    runtimeFailed: "计算引擎未连接",
    eyebrow: "REPRESENTATION → HAMILTONIAN",
    heroLead: "输入对称性，",
    heroEmphasis: "构造完整的 k·p 模型。",
    heroDescription: "从表示矩阵和动量变换直接求公共不变子空间；支持幺正、反幺正、反对称以及内部连续 C∞ 自旋旋转。",
    defineStep: "01 · 定义",
    defineHeading: "设置模型空间",
    reset: "恢复示例",
    catalogMode: "空间群目录",
    manualMode: "手动输入矩阵",
    catalogIntroduction: "从 MSGCorep 选择空间群、磁群、k 点与小群共表示，自动构造表示矩阵和动量作用。",
    spaceGroupLabel: "空间群",
    magneticGroupLabel: "磁空间群",
    bzVariantLabel: "Brillouin 区",
    kPointLabel: "k 点 / k 线",
    representationFamilyLabel: "表示类型",
    singleValued: "单值表示",
    doubleValued: "双值表示（SOC）",
    coordinateBasisLabel: "动量坐标",
    cartesianBasis: "笛卡尔坐标",
    fractionalBasis: "分数坐标",
    coordinateLabel: "选中的 k 坐标",
    generatorsLabel: "生成元对 k 的作用",
    generatorsHint: "与求解器实际使用的最小小群生成元一致。",
    generatorCount: "{count} 个生成元",
    unitaryShort: "幺正",
    antiunitaryShort: "反幺正",
    corepresentationsLabel: "小群共表示",
    corepresentationsHint: "可多选并设置重数，程序自动构造直和表示。",
    multiplicityLabel: "重数",
    representationMatricesLabel: "选中共表示的生成元矩阵",
    representationMatricesHint: "展示各不可约块；重数表示该块在直和中重复的次数。",
    multiplicityShort: "重数 {multiplicity}",
    noSelectedCorepresentation: "勾选共表示后，这里会显示其生成元矩阵。",
    convertToManual: "转为手动输入矩阵",
    invalidCatalogMultiplicity: "共表示重数必须是 1 到 8 的整数。",
    catalogManualDimensionTooLarge: "直和后的矩阵阶数为 {dimension}，手动编辑器最多支持 {maximum} 阶。",
    catalogLoading: "正在读取 MSGCorep 目录…",
    catalogUnavailable: "MSGCorep 目录不可用，已切换到手动矩阵输入。",
    staticPreview: "当前是本地静态预览，样式和示例可正常显示；计算功能请通过 MagneticKP Web 服务器打开。",
    noCorepresentations: "这个 k 点没有可选的此类共表示。",
    selectCorepresentation: "请至少选择一个共表示。",
    catalogDimensionLabel: "H 阶数: {dimension}",
    hilbertDimension: "哈密顿量阶数",
    spatialDimension: "动量空间维数",
    orderLabel: "截止阶数",
    methodLabel: "求解方法",
    methodAutomatic: "Automatic · 推荐",
    methodIterative: "Iterative · 生成元迭代",
    methodReynolds: "Reynolds · 有限群平均",
    advancedTitle: "高级求解设置",
    advancedDescription: "Reynolds 群输入与闭包上限",
    groupInputLabel: "离散群输入",
    groupGenerators: "Generators · 生成元",
    groupElements: "Elements · 全部群元",
    maxGroupOrder: "最大群阶",
    operationsStep: "02 · 对称性",
    operationsHeading: "表示与动量变换",
    operationCount: "{count} 个操作",
    operationType: "操作类型",
    operationLabel: "标签",
    typeUnitary: "幺正",
    typeAntiunitary: "反幺正",
    typeAntisymmetryUnitary: "幺正反对称",
    typeAntisymmetryAntiunitary: "反幺正反对称",
    typeContinuous: "内部连续 C∞",
    representationMatrix: "表示矩阵 D",
    matrixDimension: "{dimension} × {dimension}",
    momentumImage: "动量映射",
    momentumFixed: "连续内部旋转固定 k",
    removeOperation: "删除这个操作",
    continuousParameter: "连续参数",
    addDiscrete: "添加离散操作",
    addContinuous: "添加 C∞ 操作",
    syntaxTitle: "输入语法",
    syntaxDescription: "允许整数、有理数、I、pi、sqrt、exp、sin、cos、tan；乘法请写 *。连续表示只允许参数 theta。",
    solve: "构造 Hamiltonian",
    solving: "正在约化并求公共不变量…",
    examplesStep: "示例库",
    examplesHeading: "从已验证模型开始",
    examplesDescription: "点击卡片会替换上面的输入；拓扑粒子示例按《Encyclopedia of emergent particles in three-dimensional crystals》的小群表示构造。",
    loadExample: "载入示例",
    errorHeading: "无法完成计算",
    connectionError: "无法连接计算服务。",
    resultsStep: "03 · 结果",
    resultsHeading: "对称性允许的 Hamiltonian",
    copyLatex: "复制 LaTeX",
    copied: "已复制",
    downloadJson: "下载 JSON",
    metricParameters: "独立参数",
    metricDimension: "矩阵维数",
    metricOrders: "包含阶数",
    metricTime: "计算时间",
    hamiltonianDescription: "累计到指定阶数的完整结果",
    formulaMode: "公式",
    parametersHeading: "自由参数",
    parametersDescription: "每个参数对应一个不变量方向",
    noParameters: "没有允许的独立参数",
    diagnosticsHeading: "逐阶诊断",
    diagnosticsDescription: "环境维数与不变子空间维数",
    diagnosticLine: "环境 {ambient} → 不变量 {invariant}",
    bandHeading: "交互能带",
    bandDescription: "浏览器本地计算；拖动参数不会重新调用服务器",
    browserComputed: "LOCAL · 即时",
    bandPathLabel: "高对称路径",
    bandPathHint: "按路径段输入：起点:坐标 -> 终点:坐标；支持 1/3 等精确分数。坐标乘以 2π 后代入哈密顿量。",
    bandPathGeneratedHint: "已根据完整磁小余群自动生成局域高对称路径；可以直接编辑，坐标乘以 2π 后代入哈密顿量。",
    bandPointsLabel: "每段采样点",
    resetBand: "恢复默认路径",
    bandParameterLimit: "仅显示前 {count} 个参数滑块，其余参数使用默认值。",
    bandPathError: "路径段格式应为“起点:坐标 -> 终点:坐标”，并且每个点的坐标维数必须正确。",
    footerNote: "Exact symbolic symmetry constraints · Powered by SymPy",
    invalidDimension: "Hilbert 空间维数必须是 1 到 12 的整数。",
    invalidOrder: "截止阶数必须是 0 到 6 的整数。",
    continuousExists: "目前只支持一个 CInfinity 连续操作。",
  },
  en: {
    pageTitle: "MagneticKP · Symmetry-constrained k·p models",
    metaDescription: "MagneticKP constructs k·p Hamiltonians from unitary, antiunitary and continuous spin symmetries.",
    brandAria: "MagneticKP home",
    brandTagline: "SYMMETRY MODEL BUILDER",
    helpLink: "Docs",
    switchLanguage: "切换至中文",
    runtimeConnecting: "Connecting to the calculation engine",
    runtimeReady: "v{version} · SymPy {sympy}",
    runtimeFailed: "Calculation engine unavailable",
    eyebrow: "REPRESENTATION → HAMILTONIAN",
    heroLead: "Enter symmetries.",
    heroEmphasis: "Build the complete k·p model.",
    heroDescription: "Compute the common invariant space directly from representation matrices and momentum maps, including unitary, antiunitary, antisymmetry and internal continuous C∞ spin rotations.",
    defineStep: "01 · DEFINE",
    defineHeading: "Set the model space",
    reset: "Restore example",
    catalogMode: "Space-group catalog",
    manualMode: "Enter matrices manually",
    catalogIntroduction: "Select a space group, magnetic group, k point and little-group corepresentation from MSGCorep; representation matrices and momentum actions are constructed automatically.",
    spaceGroupLabel: "Space group",
    magneticGroupLabel: "Magnetic space group",
    bzVariantLabel: "Brillouin zone",
    kPointLabel: "k point / line",
    representationFamilyLabel: "Representation family",
    singleValued: "Single-valued",
    doubleValued: "Double-valued (SOC)",
    coordinateBasisLabel: "Momentum coordinates",
    cartesianBasis: "Cartesian",
    fractionalBasis: "Fractional",
    coordinateLabel: "Selected k coordinate",
    generatorsLabel: "Generator actions on k",
    generatorsHint: "The same minimal little-group generators used by the solver.",
    generatorCount: "{count} generator(s)",
    unitaryShort: "unitary",
    antiunitaryShort: "antiunitary",
    corepresentationsLabel: "Little-group corepresentations",
    corepresentationsHint: "Select several entries and multiplicities to construct a direct sum.",
    multiplicityLabel: "Multiplicity",
    representationMatricesLabel: "Generator matrices of selected corepresentations",
    representationMatricesHint: "Each irreducible block is shown once; its multiplicity gives the number of copies in the direct sum.",
    multiplicityShort: "multiplicity {multiplicity}",
    noSelectedCorepresentation: "Select a corepresentation to display its generator matrices.",
    convertToManual: "Convert to manual matrices",
    invalidCatalogMultiplicity: "Corepresentation multiplicities must be integers from 1 to 8.",
    catalogManualDimensionTooLarge: "The direct sum has dimension {dimension}; the manual editor supports at most {maximum}.",
    catalogLoading: "Loading the MSGCorep catalog…",
    catalogUnavailable: "The MSGCorep catalog is unavailable; switched to manual matrix input.",
    staticPreview: "This is a local static preview: styling and examples work, but calculations require the MagneticKP Web server.",
    noCorepresentations: "No corepresentations of this family are available at this k point.",
    selectCorepresentation: "Select at least one corepresentation.",
    catalogDimensionLabel: "H dimension: {dimension}",
    hilbertDimension: "Hamiltonian dimension",
    spatialDimension: "Momentum dimension",
    orderLabel: "Maximum order",
    methodLabel: "Solver",
    methodAutomatic: "Automatic · Recommended",
    methodIterative: "Iterative · Generators",
    methodReynolds: "Reynolds · Finite-group average",
    advancedTitle: "Advanced solver settings",
    advancedDescription: "Reynolds group input and closure limit",
    groupInputLabel: "Discrete group input",
    groupGenerators: "Generators",
    groupElements: "All elements",
    maxGroupOrder: "Maximum group order",
    operationsStep: "02 · SYMMETRIES",
    operationsHeading: "Representations and momentum maps",
    operationCount: "{count} operation(s)",
    operationType: "Operation type",
    operationLabel: "Label",
    typeUnitary: "Unitary",
    typeAntiunitary: "Antiunitary",
    typeAntisymmetryUnitary: "Unitary antisymmetry",
    typeAntisymmetryAntiunitary: "Antiunitary antisymmetry",
    typeContinuous: "Internal continuous C∞",
    representationMatrix: "Representation matrix D",
    matrixDimension: "{dimension} × {dimension}",
    momentumImage: "Momentum image",
    momentumFixed: "internal rotation fixes k",
    removeOperation: "Remove this operation",
    continuousParameter: "Continuous parameter",
    addDiscrete: "Add discrete operation",
    addContinuous: "Add C∞ operation",
    syntaxTitle: "Expression syntax",
    syntaxDescription: "Use integers, rationals, I, pi, sqrt, exp, sin, cos and tan; write multiplication as *. A continuous representation may contain theta only.",
    solve: "Construct Hamiltonian",
    solving: "Reducing and solving the invariant space…",
    examplesStep: "EXAMPLE LIBRARY",
    examplesHeading: "Start from a verified model",
    examplesDescription: "The topological-particle examples follow the little-group representations in Encyclopedia of emergent particles in three-dimensional crystals.",
    loadExample: "Load example",
    errorHeading: "Unable to complete the calculation",
    connectionError: "Unable to connect to the calculation service.",
    resultsStep: "03 · RESULT",
    resultsHeading: "Symmetry-allowed Hamiltonian",
    copyLatex: "Copy LaTeX",
    copied: "Copied",
    downloadJson: "Download JSON",
    metricParameters: "Independent parameters",
    metricDimension: "Matrix dimension",
    metricOrders: "Included orders",
    metricTime: "Calculation time",
    hamiltonianDescription: "Complete result through the selected order",
    formulaMode: "Formula",
    parametersHeading: "Free parameters",
    parametersDescription: "One parameter for each invariant direction",
    noParameters: "No independent parameters are allowed",
    diagnosticsHeading: "Per-order diagnostics",
    diagnosticsDescription: "Ambient and invariant-space dimensions",
    diagnosticLine: "ambient {ambient} → invariant {invariant}",
    bandHeading: "Interactive bands",
    bandDescription: "Evaluated locally in the browser; sliders do not call the server",
    browserComputed: "LOCAL · LIVE",
    bandPathLabel: "High-symmetry path",
    bandPathHint: "Enter start:coordinates -> end:coordinates segments; exact fractions such as 1/3 are supported. Coordinates are multiplied by 2π before substitution.",
    bandPathGeneratedHint: "A local high-symmetry path was generated from the complete magnetic little co-group. It remains editable; coordinates are multiplied by 2π before substitution.",
    bandPointsLabel: "Samples per segment",
    resetBand: "Reset path",
    bandParameterLimit: "Only the first {count} parameter sliders are shown; the rest use defaults.",
    bandPathError: "Use segments of the form start:coordinates -> end:coordinates with the correct coordinate dimension.",
    footerNote: "Exact symbolic symmetry constraints · Powered by SymPy",
    invalidDimension: "The Hilbert-space dimension must be an integer from 1 to 12.",
    invalidOrder: "The maximum order must be an integer from 0 to 6.",
    continuousExists: "Only one CInfinity operation is currently supported.",
  },
};

const examples = [
  {
    id: "charge-1-weyl-sg92",
    icon: "C−1",
    title: { zh: "C-1 Weyl 点（SG 92 · Γ）", en: "C-1 Weyl point (SG 92 · Γ)" },
    meta: "|C| = 1 · MSG 92 · Γ · ORDER 1",
    dimension: 2,
    spatialDimension: 3,
    order: 1,
    method: "automatic",
    operations: [
      {
        kind: "unitary", label: "C4z",
        representation: [["sqrt(2)/2", "sqrt(2)/2"], ["-sqrt(2)/2", "sqrt(2)/2"]],
        momentum_image: ["-ky", "kx", "kz"],
      },
      {
        kind: "unitary", label: "C2x",
        representation: [["0", "I"], ["I", "0"]],
        momentum_image: ["kx", "-ky", "-kz"],
      },
      {
        kind: "antiunitary", label: "T",
        representation: [["0", "-1"], ["1", "0"]],
        momentum_image: ["-kx", "-ky", "-kz"],
      },
    ],
  },
  {
    id: "weyl-nodal-line-pt",
    icon: "WNL",
    title: { zh: "Weyl 节线（PT）", en: "Weyl nodal line (PT)" },
    meta: "WNL · BERRY PHASE π · ORDER 2",
    dimension: 2,
    spatialDimension: 3,
    order: 2,
    method: "automatic",
    operations: [{
      kind: "antiunitary", label: "PT",
      representation: [["1", "0"], ["0", "1"]],
      momentum_image: ["kx", "ky", "kz"],
    }],
  },
  {
    id: "charge-4-weyl-sg211",
    icon: "C−4",
    title: { zh: "C-4 Weyl 点（SG 211 · Γ）", en: "C-4 Weyl point (SG 211 · Γ)" },
    meta: "|C| = 4 · SG 211 · Γ · ORDER 3",
    dimension: 2,
    spatialDimension: 3,
    order: 3,
    method: "automatic",
    operations: [
      {
        kind: "unitary", label: "C3m",
        representation: [["-1/2", "sqrt(3)/2"], ["-sqrt(3)/2", "-1/2"]],
        momentum_image: ["ky", "kz", "kx"],
      },
      {
        kind: "unitary", label: "C2z",
        representation: [["1", "0"], ["0", "1"]],
        momentum_image: ["-kx", "-ky", "kz"],
      },
      {
        kind: "unitary", label: "C2x",
        representation: [["1", "0"], ["0", "1"]],
        momentum_image: ["kx", "-ky", "-kz"],
      },
      {
        kind: "unitary", label: "C2a",
        representation: [["1", "0"], ["0", "-1"]],
        momentum_image: ["ky", "kx", "-kz"],
      },
      {
        kind: "antiunitary", label: "T",
        representation: [["1", "0"], ["0", "1"]],
        momentum_image: ["-kx", "-ky", "-kz"],
      },
    ],
  },
  {
    id: "continuous-direct-product",
    icon: "×",
    title: { zh: "连续对称性直积", en: "Continuous-symmetry direct product" },
    meta: "C∞ × C₂ · 4 × 4 · ORDER 1",
    dimension: 4,
    spatialDimension: 3,
    order: 1,
    method: "automatic",
    operations: [
      {
        kind: "continuous", label: "CInfinity",
        representation: [
          ["exp(I*theta/2)", "0", "0", "0"],
          ["0", "exp(I*theta/2)", "0", "0"],
          ["0", "0", "exp(I*theta/2)", "0"],
          ["0", "0", "0", "exp(I*theta/2)"],
        ],
        momentum_image: ["kx", "ky", "kz"], parameter: "theta",
      },
      {
        kind: "unitary", label: "IndependentC2",
        representation: [
          ["0", "1", "0", "0"],
          ["1", "0", "0", "0"],
          ["0", "0", "0", "1"],
          ["0", "0", "1", "0"],
        ],
        momentum_image: ["kx", "ky", "kz"],
      },
    ],
  },
  {
    id: "continuous-semidirect-product",
    icon: "⋊",
    title: { zh: "连续对称性半直积", en: "Continuous-symmetry semidirect product" },
    meta: "C∞ ⋊ C₂ · 4 × 4 · ORDER 1",
    dimension: 4,
    spatialDimension: 3,
    order: 1,
    method: "automatic",
    operations: [
      {
        kind: "continuous", label: "CInfinity",
        representation: [
          ["exp(I*theta/2)", "0", "0", "0"],
          ["0", "exp(I*theta/2)", "0", "0"],
          ["0", "0", "exp(I*theta/2)", "0"],
          ["0", "0", "0", "exp(I*theta/2)"],
        ],
        momentum_image: ["kx", "ky", "kz"], parameter: "theta",
      },
      {
        kind: "antiunitary", label: "ThetaFlip",
        representation: [
          ["1", "0", "0", "0"],
          ["0", "1", "0", "0"],
          ["0", "0", "1", "0"],
          ["0", "0", "0", "1"],
        ],
        momentum_image: ["kx", "ky", "kz"],
      },
    ],
  },
  {
    id: "time-reversal",
    icon: "𝒯",
    title: { zh: "自旋时间反演", en: "Spinful time reversal" },
    meta: "T² = −1 · 2 × 2 · ORDER 1",
    dimension: 2,
    spatialDimension: 3,
    order: 1,
    method: "automatic",
    operations: [{
      kind: "antiunitary", label: "T",
      representation: [["0", "1"], ["-1", "0"]],
      momentum_image: ["-kx", "-ky", "-kz"],
    }],
  },
];

const metaDescription = document.getElementById("metaDescription");
const languageToggle = document.getElementById("languageToggle");
const runtimeBadge = document.getElementById("runtimeBadge");
const runtimeText = document.getElementById("runtimeText");
const solverForm = document.getElementById("solverForm");
const inputModeTabs = document.getElementById("inputModeTabs");
const catalogModeButton = document.getElementById("catalogModeButton");
const manualModeButton = document.getElementById("manualModeButton");
const catalogPanel = document.getElementById("catalogPanel");
const catalogSpaceGroup = document.getElementById("catalogSpaceGroup");
const catalogMagneticGroup = document.getElementById("catalogMagneticGroup");
const catalogBzVariant = document.getElementById("catalogBzVariant");
const catalogKPoint = document.getElementById("catalogKPoint");
const catalogRepresentationFamily = document.getElementById("catalogRepresentationFamily");
const catalogCoordinateBasis = document.getElementById("catalogCoordinateBasis");
const catalogCoordinateMath = document.getElementById("catalogCoordinateMath");
const catalogParameterFields = document.getElementById("catalogParameterFields");
const catalogGeneratorList = document.getElementById("catalogGeneratorList");
const catalogGeneratorCount = document.getElementById("catalogGeneratorCount");
const corepresentationList = document.getElementById("corepresentationList");
const catalogRepresentationMatrices = document.getElementById("catalogRepresentationMatrices");
const convertCatalogToManualButton = document.getElementById("convertCatalogToManualButton");
const catalogDimension = document.getElementById("catalogDimension");
const manualSymmetryPanel = document.getElementById("manualSymmetryPanel");
const hilbertDimension = document.getElementById("hilbertDimension");
const spatialDimension = document.getElementById("spatialDimension");
const spatialDimensionTabs = document.getElementById("spatialDimensionTabs");
const orderInput = document.getElementById("orderInput");
const methodInput = document.getElementById("methodInput");
const groupInput = document.getElementById("groupInput");
const maxGroupOrder = document.getElementById("maxGroupOrder");
const operationList = document.getElementById("operationList");
const operationCount = document.getElementById("operationCount");
const addDiscreteButton = document.getElementById("addDiscreteButton");
const addContinuousButton = document.getElementById("addContinuousButton");
const resetButton = document.getElementById("resetButton");
const solveButton = document.getElementById("solveButton");
const exampleGrid = document.getElementById("exampleGrid");
const errorPanel = document.getElementById("errorPanel");
const errorMessage = document.getElementById("errorMessage");
const resultsSection = document.getElementById("resultsSection");
const parameterCount = document.getElementById("parameterCount");
const resultDimension = document.getElementById("resultDimension");
const resultOrders = document.getElementById("resultOrders");
const elapsedTime = document.getElementById("elapsedTime");
const hamiltonianMath = document.getElementById("hamiltonianMath");
const formulaViewButton = document.getElementById("formulaViewButton");
const texViewButton = document.getElementById("texViewButton");
const parameterList = document.getElementById("parameterList");
const diagnosticList = document.getElementById("diagnosticList");
const latexOutput = document.getElementById("latexOutput");
const copyLatexButton = document.getElementById("copyLatexButton");
const downloadButton = document.getElementById("downloadButton");
const footerVersion = document.getElementById("footerVersion");
const bandPathInput = document.getElementById("bandPathInput");
const bandPathHint = document.getElementById("bandPathHint");
const bandPointsInput = document.getElementById("bandPointsInput");
const resetBandButton = document.getElementById("resetBandButton");
const parameterSliders = document.getElementById("parameterSliders");
const bandPlot = document.getElementById("bandPlot");
const bandError = document.getElementById("bandError");

let currentLanguage = "zh";
let runtimeState = { status: "connecting", version: "", sympy: "" };
let latestResult = null;
let operationSerial = 0;
let bandParameterValues = {};
let bandFrame = null;
let generatedBandPath = null;
let inputMode = "catalog";
let catalogLoaded = false;
let catalogAvailable = true;
let catalogPoints = [];
let catalogCorepresentations = null;
let catalogRequestSerial = 0;

try {
  currentLanguage = window.localStorage.getItem("magnetickp-language") === "en" ? "en" : "zh";
} catch (_) {}

function t(key, values = {}) {
  let message = translations[currentLanguage][key] || translations.zh[key] || key;
  Object.entries(values).forEach(([name, value]) => {
    message = message.replaceAll(`{${name}}`, value);
  });
  return message;
}

function localized(copy) { return copy[currentLanguage] || copy.zh; }

function setLanguage(language) {
  currentLanguage = language === "en" ? "en" : "zh";
  document.documentElement.lang = currentLanguage === "zh" ? "zh-CN" : "en";
  document.title = t("pageTitle");
  metaDescription.content = t("metaDescription");
  document.querySelectorAll("[data-i18n]").forEach((node) => { node.textContent = t(node.dataset.i18n); });
  document.querySelectorAll("[data-i18n-aria]").forEach((node) => { node.setAttribute("aria-label", t(node.dataset.i18nAria)); });
  languageToggle.textContent = currentLanguage === "zh" ? "EN" : "中文";
  languageToggle.setAttribute("aria-label", t("switchLanguage"));
  try { window.localStorage.setItem("magnetickp-language", currentLanguage); } catch (_) {}
  renderRuntime();
  refreshDynamicLabels();
  if (catalogCorepresentations) {
    renderCatalogGenerators();
    updateCatalogDimension();
  }
  renderExamples();
  setBusy(solveButton.disabled);
  if (latestResult) renderResult(latestResult.result, latestResult.elapsed, false);
}

function renderRuntime() {
  runtimeBadge.classList.toggle("ready", runtimeState.status === "ready");
  runtimeBadge.classList.toggle("failed", runtimeState.status === "failed");
  if (runtimeState.status === "ready") {
    runtimeText.textContent = t("runtimeReady", { version: runtimeState.version, sympy: runtimeState.sympy });
  } else if (runtimeState.status === "failed") {
    runtimeText.textContent = t("runtimeFailed");
  } else {
    runtimeText.textContent = t("runtimeConnecting");
  }
}

function setInputMode(mode) {
  inputMode = mode === "manual" ? "manual" : "catalog";
  const catalogActive = inputMode === "catalog";
  catalogPanel.hidden = !catalogActive;
  manualSymmetryPanel.hidden = catalogActive;
  document.querySelectorAll(".manual-only").forEach((element) => {
    element.hidden = catalogActive;
  });
  catalogModeButton.classList.toggle("active", catalogActive);
  manualModeButton.classList.toggle("active", !catalogActive);
  catalogModeButton.setAttribute("aria-selected", String(catalogActive));
  manualModeButton.setAttribute("aria-selected", String(!catalogActive));
  if (catalogActive) {
    spatialDimension.value = "3";
    syncSpatialDimensionTabs();
  }
  hideError();
  resultsSection.hidden = true;
}

async function apiJson(url, options) {
  const response = await fetch(url, options);
  let data = null;
  try { data = await response.json(); } catch (_) {}
  if (!response.ok) throw new Error(data?.detail || t("connectionError"));
  return data;
}

function replaceSelectOptions(select, rows, valueOf, labelOf, preferred) {
  select.replaceChildren();
  rows.forEach((row) => {
    const option = document.createElement("option");
    option.value = String(valueOf(row));
    option.textContent = labelOf(row);
    select.append(option);
  });
  if (preferred !== undefined && rows.some((row) => String(valueOf(row)) === String(preferred))) {
    select.value = String(preferred);
  }
}

function unicodeCrystallographicSymbol(tex, fallback = "") {
  const subscriptDigits = {
    0: "₀", 1: "₁", 2: "₂", 3: "₃", 4: "₄",
    5: "₅", 6: "₆", 7: "₇", 8: "₈", 9: "₉",
  };
  const source = String(tex || fallback);
  const subscript = (value) => /^[0-9]+$/.test(value)
    ? Array.from(value, (digit) => subscriptDigits[digit]).join("")
    : `₍${value}₎`;
  return source
    .replace(/\\bar\{([^{}]+)\}/g, (_, value) => Array.from(value, (character) => `${character}̅`).join(""))
    .replace(/_\{([^{}]+)\}|_(.)/g, (_, braced, single) => subscript(braced || single))
    .replaceAll("'", "′");
}

function selectedCatalogPoint() {
  return catalogPoints.find((point) => point.name === catalogKPoint.value) || null;
}

function renderCatalogCoordinate() {
  const point = selectedCatalogPoint();
  catalogParameterFields.replaceChildren();
  if (!point) {
    catalogCoordinateMath.textContent = "—";
    return;
  }
  renderInlineMath(
    catalogCoordinateMath,
    String.raw`\mathbf{k}=${point.coordinate_latex}`,
    `k = (${point.coordinate.join(", ")})`,
  );
  point.parameters.forEach((name) => {
    const field = document.createElement("label");
    field.className = "catalog-parameter-field";
    const symbol = document.createElement("span");
    renderInlineMath(symbol, name, name);
    const input = document.createElement("input");
    input.type = "text";
    input.spellcheck = false;
    input.value = "1/5";
    input.dataset.parameter = name;
    input.setAttribute("aria-label", `k parameter ${name}`);
    field.append(symbol, input);
    catalogParameterFields.append(field);
  });
}

function selectedCatalogRepresentations() {
  return Array.from(corepresentationList.querySelectorAll(".corepresentation-option"))
    .filter((row) => row.querySelector('input[type="checkbox"]').checked)
    .map((row) => ({
      index: Number(row.dataset.index),
      multiplicity: Number(row.querySelector('input[type="number"]').value),
    }));
}

function updateCatalogDimension() {
  let dimension = 0;
  corepresentationList.querySelectorAll(".corepresentation-option").forEach((row) => {
    if (row.querySelector('input[type="checkbox"]').checked) {
      dimension += Number(row.dataset.dimension) * Number(row.querySelector('input[type="number"]').value);
    }
  });
  catalogDimension.textContent = dimension ? t("catalogDimensionLabel", { dimension }) : "H: —";
  if (dimension) hilbertDimension.value = dimension;
  renderRepresentationMatrices();
}

function renderRepresentationMatrices() {
  catalogRepresentationMatrices.replaceChildren();
  const family = catalogRepresentationFamily.value;
  const rows = catalogCorepresentations?.[family] || [];
  const generators = catalogCorepresentations?.generators || [];
  const selected = selectedCatalogRepresentations();
  if (!selected.length) {
    const message = document.createElement("p");
    message.className = "catalog-empty";
    message.textContent = t("noSelectedCorepresentation");
    catalogRepresentationMatrices.append(message);
    return;
  }
  selected.forEach((selection) => {
    const item = rows.find((row) => Number(row.index) === selection.index);
    if (!item) return;
    const block = document.createElement("section");
    block.className = "catalog-representation-block";
    const heading = document.createElement("div");
    heading.className = "catalog-representation-block-heading";
    const identity = document.createElement("strong");
    renderSubscriptedText(identity, item.label);
    const detail = document.createElement("span");
    detail.textContent = `${item.dimension}D · type ${item.type}`;
    const multiplicity = document.createElement("b");
    multiplicity.textContent = t("multiplicityShort", { multiplicity: selection.multiplicity });
    heading.append(identity, detail, multiplicity);

    const matrixGrid = document.createElement("div");
    matrixGrid.className = "catalog-representation-matrix-grid";
    (item.generator_matrices || []).forEach((matrix) => {
      const generator = generators.find((row) => Number(row.index) === Number(matrix.generator_index));
      const card = document.createElement("article");
      card.className = "catalog-representation-matrix";
      const cardHeading = document.createElement("div");
      cardHeading.className = "catalog-representation-matrix-heading";
      const name = document.createElement("strong");
      name.textContent = `D(g${matrix.generator_index})${generator ? ` · ${generator.rotation}` : ""}`;
      const badge = document.createElement("span");
      badge.className = matrix.antiunitary ? "antiunitary" : "unitary";
      badge.textContent = t(matrix.antiunitary ? "antiunitaryShort" : "unitaryShort");
      cardHeading.append(name, badge);
      const math = document.createElement("div");
      math.className = "catalog-representation-matrix-math";
      renderInlineMath(
        math,
        String.raw`D(g_{${matrix.generator_index}})=${matrix.matrix_latex}`,
        `D(g${matrix.generator_index}) = ${JSON.stringify(matrix.matrix)}`,
      );
      card.append(cardHeading, math);
      matrixGrid.append(card);
    });
    block.append(heading, matrixGrid);
    catalogRepresentationMatrices.append(block);
  });
}

function blockDiagonalMatrix(blocks) {
  const dimension = blocks.reduce((total, block) => total + block.length, 0);
  const result = Array.from({ length: dimension }, () => Array(dimension).fill("0"));
  let offset = 0;
  blocks.forEach((block) => {
    block.forEach((row, rowIndex) => {
      row.forEach((value, columnIndex) => {
        result[offset + rowIndex][offset + columnIndex] = value;
      });
    });
    offset += block.length;
  });
  return result;
}

function convertCatalogToManual() {
  try {
    const family = catalogRepresentationFamily.value;
    const rows = catalogCorepresentations?.[family] || [];
    const selected = selectedCatalogRepresentations();
    if (!selected.length) throw new Error(t("selectCorepresentation"));

    const expanded = [];
    selected.forEach((selection) => {
      if (!Number.isInteger(selection.multiplicity) || selection.multiplicity < 1 || selection.multiplicity > 8) {
        throw new Error(t("invalidCatalogMultiplicity"));
      }
      const item = rows.find((row) => Number(row.index) === selection.index);
      if (!item) throw new Error(t("selectCorepresentation"));
      for (let copy = 0; copy < selection.multiplicity; copy += 1) expanded.push(item);
    });

    const dimension = expanded.reduce((total, item) => total + Number(item.dimension), 0);
    const maximum = Number(hilbertDimension.max) || 12;
    if (dimension > maximum) {
      throw new Error(t("catalogManualDimensionTooLarge", { dimension, maximum }));
    }
    const generators = catalogCorepresentations?.generators || [];
    const operations = generators.map((generator) => {
      const blocks = expanded.map((item) => {
        const entry = (item.generator_matrices || []).find(
          (matrix) => Number(matrix.generator_index) === Number(generator.index),
        );
        if (!entry) throw new Error(`Missing D(g${generator.index}) for ${item.label}.`);
        return entry.matrix;
      });
      return {
        kind: generator.antiunitary ? "antiunitary" : "unitary",
        label: generator.label || `g${generator.index}`,
        representation: blockDiagonalMatrix(blocks),
        momentum_image: generator.momentum_image,
      };
    });

    hilbertDimension.value = dimension;
    spatialDimension.value = "3";
    syncSpatialDimensionTabs();
    groupInput.value = "generators";
    operationList.replaceChildren();
    operations.forEach((operation) => createOperation(operation));
    setInputMode("manual");
    hideError();
    manualSymmetryPanel.scrollIntoView({ behavior: "smooth", block: "start" });
  } catch (error) {
    showError(error.message);
  }
}

function renderCatalogGenerators() {
  catalogGeneratorList.replaceChildren();
  const generators = catalogCorepresentations?.generators || [];
  catalogGeneratorCount.textContent = generators.length
    ? t("generatorCount", { count: generators.length })
    : "—";
  generators.forEach((generator) => {
    const card = document.createElement("article");
    card.className = "catalog-generator";
    const heading = document.createElement("div");
    heading.className = "catalog-generator-heading";
    const badge = document.createElement("span");
    badge.className = generator.antiunitary ? "antiunitary" : "unitary";
    badge.textContent = t(generator.antiunitary ? "antiunitaryShort" : "unitaryShort");
    const identity = document.createElement("strong");
    identity.textContent = `g${generator.index} · ${generator.rotation}`;
    heading.append(badge, identity);
    const translation = document.createElement("small");
    translation.textContent = `t = (${generator.translation.join(", ")})`;
    const action = document.createElement("div");
    action.className = "catalog-generator-action";
    renderInlineMath(
      action,
      String.raw`\mathbf{k}\mapsto ${generator.momentum_image_latex}`,
      `k -> (${generator.momentum_image.join(", ")})`,
    );
    card.append(heading, translation, action);
    catalogGeneratorList.append(card);
  });
}

function renderCorepresentations() {
  corepresentationList.replaceChildren();
  const family = catalogRepresentationFamily.value;
  const rows = catalogCorepresentations?.[family] || [];
  if (!rows.length) {
    const message = document.createElement("p");
    message.className = "catalog-empty";
    message.textContent = t("noCorepresentations");
    corepresentationList.append(message);
    updateCatalogDimension();
    return;
  }
  rows.forEach((item, position) => {
    const row = document.createElement("label");
    row.className = "corepresentation-option";
    row.dataset.index = item.index;
    row.dataset.dimension = item.dimension;
    const checkbox = document.createElement("input");
    checkbox.type = "checkbox";
    checkbox.checked = position === 0;
    const description = document.createElement("span");
    description.className = "corepresentation-description";
    const name = document.createElement("strong");
    renderSubscriptedText(name, item.label);
    const detail = document.createElement("small");
    detail.textContent = `${item.dimension}D · type ${item.type}`;
    description.append(name, detail);
    const multiplicity = document.createElement("input");
    multiplicity.type = "number";
    multiplicity.min = "1";
    multiplicity.max = "8";
    multiplicity.value = "1";
    multiplicity.setAttribute("aria-label", `${item.label} ${t("multiplicityLabel")}`);
    checkbox.addEventListener("change", updateCatalogDimension);
    multiplicity.addEventListener("input", updateCatalogDimension);
    row.append(checkbox, description, multiplicity);
    corepresentationList.append(row);
  });
  updateCatalogDimension();
}

async function loadCorepresentations() {
  const magneticGroup = catalogMagneticGroup.value;
  const kName = catalogKPoint.value;
  if (!magneticGroup || !kName) return;
  const serial = ++catalogRequestSerial;
  corepresentationList.textContent = t("catalogLoading");
  catalogGeneratorList.textContent = t("catalogLoading");
  catalogGeneratorCount.textContent = "—";
  const query = new URLSearchParams({
    magnetic_group: magneticGroup,
    k_name: kName,
    variant: catalogBzVariant.value || "a",
    coordinate_basis: catalogCoordinateBasis.value,
  });
  const data = await apiJson(`api/catalog/corepresentations?${query}`);
  if (serial !== catalogRequestSerial) return;
  catalogCorepresentations = data;
  renderCatalogGenerators();
  renderCorepresentations();
}

async function loadKPoints(preferred = "Γ") {
  const spaceGroup = catalogSpaceGroup.value;
  const magneticGroup = catalogMagneticGroup.value;
  if (!spaceGroup || !magneticGroup) return;
  const query = new URLSearchParams({
    magnetic_group: magneticGroup,
    variant: catalogBzVariant.value || "a",
  });
  const data = await apiJson(`api/catalog/k-points/${spaceGroup}?${query}`);
  catalogPoints = data.points;
  replaceSelectOptions(
    catalogKPoint,
    catalogPoints,
    (point) => point.name,
    (point) => `${point.name}  ·  (${point.coordinate.join(", ")})`,
    preferred,
  );
  renderCatalogCoordinate();
  await loadCorepresentations();
}

async function loadBzVariants(preferred = "a") {
  const data = await apiJson(`api/catalog/bz-variants/${catalogSpaceGroup.value}`);
  replaceSelectOptions(
    catalogBzVariant,
    data.variants,
    (variant) => variant.value,
    (variant) => variant.label,
    preferred,
  );
}

async function loadMagneticGroups(preferred = "92.112", preferredK = "Γ") {
  const spaceGroup = catalogSpaceGroup.value;
  const [groups] = await Promise.all([
    apiJson(`api/catalog/magnetic-groups/${spaceGroup}`),
    loadBzVariants("a"),
  ]);
  replaceSelectOptions(
    catalogMagneticGroup,
    groups.groups,
    (group) => group.identifier,
    (group) => `${group.identifier}  ·  type ${group.type}  ·  ${unicodeCrystallographicSymbol(group.bns_symbol_tex, group.bns_symbol)}`,
    preferred,
  );
  await loadKPoints(preferredK);
}

async function loadCatalog() {
  if (catalogLoaded || !catalogAvailable) return;
  catalogSpaceGroup.replaceChildren(new Option(t("catalogLoading"), ""));
  try {
    const data = await apiJson("api/catalog/space-groups");
    replaceSelectOptions(
      catalogSpaceGroup,
      data.groups,
      (group) => group.number,
      (group) => `${group.number}  ·  ${unicodeCrystallographicSymbol(group.symbol_tex, group.symbol)}`,
      92,
    );
    await loadMagneticGroups("92.112", "Γ");
    catalogLoaded = true;
  } catch (error) {
    catalogAvailable = false;
    catalogModeButton.disabled = true;
    setInputMode("manual");
    showError(error.message || t("catalogUnavailable"));
  }
}

function collectCatalogRequest() {
  const order = Number(orderInput.value);
  if (!Number.isInteger(order) || order < 0 || order > 6) throw new Error(t("invalidOrder"));
  const representations = selectedCatalogRepresentations();
  if (!representations.length) throw new Error(t("selectCorepresentation"));
  const kParameters = {};
  catalogParameterFields.querySelectorAll("input").forEach((input) => {
    kParameters[input.dataset.parameter] = input.value.trim();
  });
  return {
    space_group: Number(catalogSpaceGroup.value),
    magnetic_group: catalogMagneticGroup.value,
    bz_variant: catalogBzVariant.value || "a",
    k_name: catalogKPoint.value,
    k_parameters: kParameters,
    coordinate_basis: catalogCoordinateBasis.value,
    representation_family: catalogRepresentationFamily.value,
    representations,
    order,
    method: methodInput.value,
    max_group_order: Number(maxGroupOrder.value),
  };
}

function dimensionValue() {
  const value = Number(hilbertDimension.value);
  return Number.isInteger(value) && value >= 1 && value <= 12 ? value : null;
}

function spatialValue() { return Number(spatialDimension.value); }
function momentumNames() { return spatialValue() === 1 ? ["kz"] : spatialValue() === 2 ? ["kx", "ky"] : ["kx", "ky", "kz"]; }

function renderSubscriptedText(element, text) {
  const pattern = /\bk([xyz])\b|\b([A-Za-z]+)_([0-9]+(?:_[0-9]+)*)\b/g;
  const fragments = [];
  let cursor = 0;
  for (const match of text.matchAll(pattern)) {
    fragments.push(document.createTextNode(text.slice(cursor, match.index)));
    fragments.push(document.createTextNode(match[1] ? "k" : match[2]));
    const subscript = document.createElement("sub");
    subscript.textContent = match[1] || match[3].replaceAll("_", ",");
    fragments.push(subscript);
    cursor = match.index + match[0].length;
  }
  fragments.push(document.createTextNode(text.slice(cursor)));
  element.replaceChildren(...fragments);
}

function displayLatex(latex) {
  return latex
    .replace(/\bkx\b/g, "k_{x}")
    .replace(/\bky\b/g, "k_{y}")
    .replace(/\bkz\b/g, "k_{z}");
}

function renderInlineMath(element, latex, fallback) {
  if (window.katex?.render) {
    window.katex.render(latex, element, {
      displayMode: false,
      throwOnError: false,
      trust: false,
      output: "htmlAndMathml",
    });
  } else {
    renderSubscriptedText(element, fallback);
  }
}

function syncSpatialDimensionTabs() {
  spatialDimensionTabs.querySelectorAll(".dimension-tab").forEach((button) => {
    const active = button.dataset.dimension === spatialDimension.value;
    button.classList.toggle("active", active);
    button.setAttribute("aria-checked", String(active));
  });
}

function renderSpatialDimensionTabs() {
  spatialDimensionTabs.querySelectorAll(".dimension-tab-math").forEach((element) => {
    renderInlineMath(element, element.dataset.latex, element.dataset.fallback);
  });
  syncSpatialDimensionTabs();
}

function identityMatrix(dimension) {
  return Array.from({ length: dimension }, (_, row) =>
    Array.from({ length: dimension }, (_, column) => row === column ? "1" : "0"),
  );
}

function readMatrix(card) {
  const dimension = dimensionValue() || 1;
  const rows = Array.from({ length: dimension }, () => Array(dimension).fill("0"));
  card.querySelectorAll(".matrix-grid input").forEach((input) => {
    const row = Number(input.dataset.row);
    const column = Number(input.dataset.column);
    if (row < dimension && column < dimension) rows[row][column] = input.value.trim() || "0";
  });
  return rows;
}

function renderMatrix(card, values) {
  const dimension = dimensionValue() || 1;
  const matrix = values || identityMatrix(dimension);
  const grid = card.querySelector(".matrix-grid");
  grid.replaceChildren();
  grid.style.gridTemplateColumns = `repeat(${dimension}, minmax(58px, 1fr))`;
  for (let row = 0; row < dimension; row += 1) {
    for (let column = 0; column < dimension; column += 1) {
      const input = document.createElement("input");
      input.type = "text";
      input.spellcheck = false;
      input.dataset.row = row;
      input.dataset.column = column;
      input.value = matrix[row]?.[column] ?? (row === column ? "1" : "0");
      input.setAttribute("aria-label", `D ${row + 1}, ${column + 1}`);
      grid.append(input);
    }
  }
  card.querySelector(".matrix-size").textContent = t("matrixDimension", { dimension });
}

function readMomentum(card) {
  return Array.from(card.querySelectorAll(".momentum-grid input"), (input) => input.value.trim());
}

function renderMomentum(card, values) {
  const names = momentumNames();
  const continuous = card.dataset.kind === "continuous";
  const grid = card.querySelector(".momentum-grid");
  grid.replaceChildren();
  names.forEach((name, index) => {
    const row = document.createElement("div");
    row.className = "momentum-row";
    const symbol = document.createElement("span");
    symbol.className = "symbol-label momentum-symbol";
    renderInlineMath(symbol, displayLatex(name), name);
    const input = document.createElement("input");
    input.type = "text";
    input.spellcheck = false;
    input.value = continuous ? name : (values?.[index] || name);
    input.disabled = continuous;
    input.setAttribute("aria-label", `${name} image`);
    row.append(symbol, input);
    grid.append(row);
  });
  card.querySelector(".momentum-note").textContent = continuous ? t("momentumFixed") : "Rk";
}

function operationTypeOptions(select) {
  const options = [
    ["unitary", "typeUnitary"],
    ["antiunitary", "typeAntiunitary"],
    ["antisymmetry_unitary", "typeAntisymmetryUnitary"],
    ["antisymmetry_antiunitary", "typeAntisymmetryAntiunitary"],
    ["continuous", "typeContinuous"],
  ];
  const selected = select.value;
  select.replaceChildren();
  options.forEach(([value, key]) => {
    const option = document.createElement("option");
    option.value = value;
    option.textContent = t(key);
    select.append(option);
  });
  select.value = selected;
}

function createOperation(operation = {}) {
  operationSerial += 1;
  const card = document.createElement("article");
  const kind = operation.kind || "unitary";
  card.className = `operation-card${kind === "continuous" ? " continuous" : ""}`;
  card.dataset.kind = kind;
  card.dataset.operationId = operationSerial;

  const top = document.createElement("div");
  top.className = "operation-top";
  const typeField = document.createElement("div");
  const typeLabel = document.createElement("label");
  typeLabel.className = "type-caption";
  typeLabel.textContent = t("operationType");
  const typeSelect = document.createElement("select");
  typeSelect.className = "kind-input";
  typeSelect.value = kind;
  operationTypeOptions(typeSelect);
  typeSelect.value = kind;
  typeField.append(typeLabel, typeSelect);

  const labelField = document.createElement("div");
  labelField.className = "label-field";
  const labelCaption = document.createElement("label");
  labelCaption.className = "label-caption";
  labelCaption.textContent = t("operationLabel");
  const labelInput = document.createElement("input");
  labelInput.className = "label-input";
  labelInput.type = "text";
  labelInput.maxLength = 64;
  labelInput.value = kind === "continuous" ? "CInfinity" : (operation.label || `g${operationSerial}`);
  labelInput.disabled = kind === "continuous";
  labelField.append(labelCaption, labelInput);

  const remove = document.createElement("button");
  remove.className = "remove-operation";
  remove.type = "button";
  remove.textContent = "×";
  remove.setAttribute("aria-label", t("removeOperation"));
  top.append(typeField, labelField, remove);

  const body = document.createElement("div");
  body.className = "operation-body";
  const matrixArea = document.createElement("div");
  matrixArea.innerHTML = `<p class="input-heading"><span class="matrix-caption"></span><small class="matrix-size"></small></p><div class="matrix-scroll"><div class="matrix-grid"></div></div>`;
  const momentumArea = document.createElement("div");
  momentumArea.innerHTML = `<p class="input-heading"><span class="momentum-caption"></span><small class="momentum-note"></small></p><div class="momentum-grid"></div><div class="continuous-parameter" hidden><span></span><code>theta</code></div>`;
  body.append(matrixArea, momentumArea);
  card.append(top, body);
  operationList.append(card);

  renderMatrix(card, operation.representation);
  renderMomentum(card, operation.momentum_image);
  syncOperationKind(card, kind, false);

  remove.addEventListener("click", () => {
    card.remove();
    refreshOperationCount();
  });
  typeSelect.addEventListener("change", () => {
    if (typeSelect.value === "continuous" && document.querySelector('.operation-card[data-kind="continuous"]')) {
      typeSelect.value = card.dataset.kind;
      showError(t("continuousExists"));
      return;
    }
    syncOperationKind(card, typeSelect.value, true);
  });
  refreshOperationCount();
  return card;
}

function syncOperationKind(card, kind, resetValues) {
  const wasContinuous = card.dataset.kind === "continuous";
  const matrix = readMatrix(card);
  const momentum = readMomentum(card);
  card.dataset.kind = kind;
  card.classList.toggle("continuous", kind === "continuous");
  const label = card.querySelector(".label-input");
  label.disabled = kind === "continuous";
  if (kind === "continuous") label.value = "CInfinity";
  else if (wasContinuous && label.value === "CInfinity") label.value = `g${card.dataset.operationId}`;
  card.querySelector(".continuous-parameter").hidden = kind !== "continuous";
  renderMatrix(card, resetValues ? matrix : matrix);
  renderMomentum(card, kind === "continuous" ? null : momentum);
  refreshOperationCount();
}

function refreshOperationCount() {
  operationCount.textContent = t("operationCount", { count: operationList.children.length });
}

function refreshDynamicLabels() {
  document.querySelectorAll(".operation-card").forEach((card) => {
    const select = card.querySelector(".kind-input");
    operationTypeOptions(select);
    card.querySelector(".type-caption").textContent = t("operationType");
    card.querySelector(".label-caption").textContent = t("operationLabel");
    card.querySelector(".matrix-caption").textContent = t("representationMatrix");
    card.querySelector(".matrix-size").textContent = t("matrixDimension", { dimension: dimensionValue() || 1 });
    card.querySelector(".momentum-caption").textContent = t("momentumImage");
    card.querySelector(".continuous-parameter span").textContent = t("continuousParameter");
    card.querySelector(".remove-operation").setAttribute("aria-label", t("removeOperation"));
    card.querySelector(".momentum-note").textContent = card.dataset.kind === "continuous" ? t("momentumFixed") : "Rk";
  });
  refreshOperationCount();
}

function operationData(card) {
  const kind = card.dataset.kind;
  return {
    kind,
    label: kind === "continuous" ? "CInfinity" : card.querySelector(".label-input").value.trim(),
    representation: readMatrix(card),
    momentum_image: readMomentum(card),
    parameter: kind === "continuous" ? "theta" : null,
  };
}

function applyExample(example, scroll = true) {
  setInputMode("manual");
  hilbertDimension.value = example.dimension;
  spatialDimension.value = example.spatialDimension;
  syncSpatialDimensionTabs();
  orderInput.value = example.order;
  methodInput.value = example.method;
  groupInput.value = "generators";
  maxGroupOrder.value = 512;
  operationList.replaceChildren();
  example.operations.forEach((operation) => createOperation(operation));
  hideError();
  resultsSection.hidden = true;
  if (scroll) solverForm.scrollIntoView({ behavior: "smooth", block: "start" });
}

function renderExamples() {
  exampleGrid.replaceChildren();
  examples.forEach((example) => {
    const button = document.createElement("button");
    button.type = "button";
    button.className = "example-card";
    const icon = document.createElement("span");
    icon.className = "example-icon";
    icon.textContent = example.icon;
    const title = document.createElement("strong");
    title.textContent = localized(example.title);
    const meta = document.createElement("small");
    meta.textContent = `${example.meta} · ${t("loadExample")}`;
    button.append(icon, title, meta);
    button.addEventListener("click", () => applyExample(example));
    exampleGrid.append(button);
  });
}

function collectRequest() {
  const dimension = dimensionValue();
  if (!dimension) throw new Error(t("invalidDimension"));
  const order = Number(orderInput.value);
  if (!Number.isInteger(order) || order < 0 || order > 6) throw new Error(t("invalidOrder"));
  const operations = Array.from(operationList.children, operationData);
  if (!operations.length) throw new Error(currentLanguage === "zh" ? "至少需要一个对称操作。" : "At least one symmetry operation is required.");
  return {
    order,
    spatial_dimension: spatialValue(),
    method: methodInput.value,
    group_input: groupInput.value,
    max_group_order: Number(maxGroupOrder.value),
    operations,
  };
}

function setBusy(busy) {
  solveButton.disabled = busy;
  solveButton.firstElementChild.textContent = t(busy ? "solving" : "solve");
}

function showError(message) {
  errorMessage.textContent = message;
  errorPanel.hidden = false;
}
function hideError() { errorPanel.hidden = true; errorMessage.textContent = ""; }

function hamiltonianLatex(result) {
  return String.raw`H(\mathbf{k}) = ` + displayLatex(result.hamiltonian.latex);
}

function setHamiltonianView(view) {
  const texActive = view === "tex";
  hamiltonianMath.hidden = texActive;
  latexOutput.hidden = !texActive;
  formulaViewButton.classList.toggle("active", !texActive);
  texViewButton.classList.toggle("active", texActive);
  formulaViewButton.setAttribute("aria-selected", String(!texActive));
  texViewButton.setAttribute("aria-selected", String(texActive));
}

function renderHamiltonian(result) {
  const latex = hamiltonianLatex(result);
  if (window.katex?.render) {
    try {
      window.katex.render(latex, hamiltonianMath, {
        displayMode: true,
        throwOnError: false,
        trust: false,
        output: "htmlAndMathml",
      });
      hamiltonianMath.dataset.renderer = "katex";
      return;
    } catch (_) {
      // Fall through to the server-provided MathML when KaTeX cannot render.
    }
  }
  hamiltonianMath.innerHTML = result.hamiltonian.mathml;
  hamiltonianMath.dataset.renderer = "mathml";
}

function renderResult(result, elapsed, scroll = true) {
  parameterCount.textContent = result.number_of_parameters;
  resultDimension.textContent = `${result.dimension} × ${result.dimension}`;
  resultOrders.textContent = result.orders.join(", ");
  elapsedTime.textContent = `${elapsed.toFixed(0)} ms`;
  renderHamiltonian(result);
  latexOutput.textContent = hamiltonianLatex(result);
  setHamiltonianView("formula");
  parameterList.replaceChildren();
  if (result.parameters.length) {
    result.parameters.forEach((parameter) => {
      const chip = document.createElement("span");
      chip.className = "parameter-chip symbol-label";
      renderSubscriptedText(chip, parameter);
      parameterList.append(chip);
    });
  } else {
    const empty = document.createElement("p");
    empty.className = "empty-note";
    empty.textContent = t("noParameters");
    parameterList.append(empty);
  }
  diagnosticList.replaceChildren();
  result.diagnostics.forEach((diagnostic) => {
    const row = document.createElement("div");
    row.className = "diagnostic-row";
    const order = document.createElement("b");
    order.textContent = `k^${diagnostic.order}`;
    const detail = document.createElement("span");
    detail.textContent = t("diagnosticLine", {
      ambient: diagnostic.ambient_dimension,
      invariant: diagnostic.invariant_dimension,
    });
    const solver = document.createElement("strong");
    solver.textContent = diagnostic.solver;
    row.append(order, detail, solver);
    diagnosticList.append(row);
  });
  prepareBandControls(result);
  resultsSection.hidden = false;
  if (scroll) resultsSection.scrollIntoView({ behavior: "smooth", block: "start" });
}

function defaultBandPath() {
  if (spatialValue() === 1) return "−K:-1/2 -> Γ:0; Γ:0 -> K:1/2";
  if (spatialValue() === 2) return "Γ:-2/3,1/3 -> K:0,0; K:0,0 -> M:0,1/2";
  return "P:-1/2,0,0 -> O:0,0,0; O:0,0,0 -> Q:0,1/2,0";
}

function prepareBandControls(result) {
  bandParameterValues = {};
  parameterSliders.replaceChildren();
  result.parameters.forEach((name, index) => {
    const initial = 0;
    bandParameterValues[name] = initial;
    if (index >= 40) return;
    const wrapper = document.createElement("div");
    wrapper.className = "parameter-slider";
    const heading = document.createElement("div");
    heading.className = "parameter-slider-heading";
    const label = document.createElement("span");
    label.className = "symbol-label";
    renderSubscriptedText(label, name);
    const output = document.createElement("output");
    output.textContent = initial.toFixed(2);
    const input = document.createElement("input");
    input.type = "range";
    input.min = "-1";
    input.max = "1";
    input.step = "0.02";
    input.value = initial;
    input.setAttribute("aria-label", name);
    input.addEventListener("input", () => {
      bandParameterValues[name] = Number(input.value);
      output.textContent = Number(input.value).toFixed(2);
      scheduleBandPlot();
    });
    heading.append(label, output);
    wrapper.append(heading, input);
    parameterSliders.append(wrapper);
  });
  if (result.parameters.length > 40) {
    const note = document.createElement("p");
    note.className = "empty-note";
    note.textContent = t("bandParameterLimit", { count: 40 });
    parameterSliders.append(note);
  }
  const automaticPath = result.high_symmetry_path;
  generatedBandPath = automaticPath?.available ? automaticPath.text : null;
  bandPathInput.value = generatedBandPath || defaultBandPath();
  bandPathHint.dataset.i18n = generatedBandPath
    ? "bandPathGeneratedHint"
    : "bandPathHint";
  bandPathHint.textContent = t(bandPathHint.dataset.i18n);
  scheduleBandPlot();
}

function complex(real = 0, imaginary = 0) { return { real, imaginary }; }
function cadd(a, b) { return complex(a.real + b.real, a.imaginary + b.imaginary); }
function cmultiply(a, b) { return complex(a.real * b.real - a.imaginary * b.imaginary, a.real * b.imaginary + a.imaginary * b.real); }
function cdivide(a, b) {
  const denominator = b.real * b.real + b.imaginary * b.imaginary;
  return complex((a.real * b.real + a.imaginary * b.imaginary) / denominator, (a.imaginary * b.real - a.real * b.imaginary) / denominator);
}
function cexp(a) {
  const scale = Math.exp(a.real);
  return complex(scale * Math.cos(a.imaginary), scale * Math.sin(a.imaginary));
}
function clog(a) { return complex(Math.log(Math.hypot(a.real, a.imaginary)), Math.atan2(a.imaginary, a.real)); }
function cpow(a, b) {
  if (a.real === 0 && a.imaginary === 0 && b.imaginary === 0) {
    if (b.real > 0) return complex(0, 0);
    if (b.real === 0) return complex(1, 0);
  }
  return cexp(cmultiply(b, clog(a)));
}
function csin(a) { return complex(Math.sin(a.real) * Math.cosh(a.imaginary), Math.cos(a.real) * Math.sinh(a.imaginary)); }
function ccos(a) { return complex(Math.cos(a.real) * Math.cosh(a.imaginary), -Math.sin(a.real) * Math.sinh(a.imaginary)); }

function evaluateTree(node, environment) {
  if (node.type === "number") return complex(node.value, 0);
  if (node.type === "complex") return complex(node.real, node.imaginary);
  if (node.type === "symbol") {
    if (!(node.name in environment)) throw new Error(`Missing numerical value for ${node.name}`);
    return complex(environment[node.name], 0);
  }
  if (node.type === "add") return node.arguments.reduce((value, item) => cadd(value, evaluateTree(item, environment)), complex());
  if (node.type === "multiply") return node.arguments.reduce((value, item) => cmultiply(value, evaluateTree(item, environment)), complex(1, 0));
  if (node.type === "power") return cpow(evaluateTree(node.base, environment), evaluateTree(node.exponent, environment));
  const value = evaluateTree(node.argument, environment);
  if (node.type === "exp") return cexp(value);
  if (node.type === "sin") return csin(value);
  if (node.type === "cos") return ccos(value);
  if (node.type === "tan") return cdivide(csin(value), ccos(value));
  throw new Error(`Unsupported numerical node ${node.type}`);
}

function hermitianEigenvalues(matrix) {
  const n = matrix.length;
  const size = n * 2;
  const real = Array.from({ length: size }, () => Array(size).fill(0));
  for (let row = 0; row < n; row += 1) {
    for (let column = 0; column < n; column += 1) {
      const value = matrix[row][column];
      real[row][column] = value.real;
      real[row][column + n] = -value.imaginary;
      real[row + n][column] = value.imaginary;
      real[row + n][column + n] = value.real;
    }
  }
  for (let row = 0; row < size; row += 1) {
    for (let column = row + 1; column < size; column += 1) {
      const average = (real[row][column] + real[column][row]) / 2;
      real[row][column] = average;
      real[column][row] = average;
    }
  }
  const maximumIterations = Math.max(40, 60 * size * size);
  for (let iteration = 0; iteration < maximumIterations; iteration += 1) {
    let p = 0;
    let q = 1;
    let maximum = 0;
    for (let row = 0; row < size; row += 1) {
      for (let column = row + 1; column < size; column += 1) {
        const magnitude = Math.abs(real[row][column]);
        if (magnitude > maximum) { maximum = magnitude; p = row; q = column; }
      }
    }
    if (maximum < 1e-10) break;
    const app = real[p][p];
    const aqq = real[q][q];
    const apq = real[p][q];
    const tau = (aqq - app) / (2 * apq);
    const tangent = tau === 0 ? 1 : Math.sign(tau) / (Math.abs(tau) + Math.sqrt(1 + tau * tau));
    const cosine = 1 / Math.sqrt(1 + tangent * tangent);
    const sine = tangent * cosine;
    for (let index = 0; index < size; index += 1) {
      if (index === p || index === q) continue;
      const aip = real[index][p];
      const aiq = real[index][q];
      real[index][p] = cosine * aip - sine * aiq;
      real[p][index] = real[index][p];
      real[index][q] = sine * aip + cosine * aiq;
      real[q][index] = real[index][q];
    }
    real[p][p] = app - tangent * apq;
    real[q][q] = aqq + tangent * apq;
    real[p][q] = 0;
    real[q][p] = 0;
  }
  const doubled = real.map((row, index) => row[index]).sort((a, b) => a - b);
  return Array.from({ length: n }, (_, index) => (doubled[2 * index] + doubled[2 * index + 1]) / 2);
}

function parseFractionalCoordinate(source) {
  const normalized = source.trim().replaceAll("−", "-");
  const parts = normalized.split("/");
  const decimalPattern = /^[+-]?(?:\d+(?:\.\d*)?|\.\d+)$/;
  if (parts.length < 1 || parts.length > 2 || parts.some((part) => !decimalPattern.test(part.trim()))) return NaN;
  const numerator = Number(parts[0]);
  if (parts.length === 1) return numerator;
  const denominator = Number(parts[1]);
  return denominator === 0 ? NaN : numerator / denominator;
}

function parseBandPoint(source, dimension) {
  const separator = source.indexOf(":");
  if (separator < 1) throw new Error(t("bandPathError"));
  const label = source.slice(0, separator).trim();
  const coordinates = source.slice(separator + 1).split(",").map(parseFractionalCoordinate);
  if (!label || coordinates.length !== dimension || coordinates.some((value) => !Number.isFinite(value))) {
    throw new Error(t("bandPathError"));
  }
  return { label, coordinates };
}

function parseBandPathValue(source, dimension) {
  const parts = source.split(";").map((part) => part.trim()).filter(Boolean);
  if (!parts.length) throw new Error(t("bandPathError"));

  // Keep accepting the first web prototype's point-list syntax, while all
  // generated paths use the Mathematica segment-list semantics.
  if (!parts.some((part) => /->|→/.test(part))) {
    const points = parts.map((part) => parseBandPoint(part, dimension));
    if (points.length < 2) throw new Error(t("bandPathError"));
    return points.slice(0, -1).map((start, index) => ({ start, end: points[index + 1] }));
  }

  return parts.map((part) => {
    const endpoints = part.split(/\s*(?:->|→)\s*/);
    if (endpoints.length !== 2) throw new Error(t("bandPathError"));
    return {
      start: parseBandPoint(endpoints[0], dimension),
      end: parseBandPoint(endpoints[1], dimension),
    };
  });
}

function parseBandPath() {
  return parseBandPathValue(bandPathInput.value, spatialValue());
}

function sampledBandPath(path, samplesPerSegment) {
  const segments = path.map((segment, segmentIndex) => {
    const samples = [];
    for (let step = 0; step <= samplesPerSegment; step += 1) {
      const fraction = step / samplesPerSegment;
      samples.push({
        distance: segmentIndex + fraction,
        coordinates: segment.start.coordinates.map(
          (value, index) => value + fraction * (segment.end.coordinates[index] - value),
        ),
      });
    }
    return samples;
  });
  const ticks = [{ label: path[0].start.label, distance: 0 }];
  path.forEach((segment, index) => {
    const next = path[index + 1];
    const label = next && next.start.label !== segment.end.label
      ? `${segment.end.label} | ${next.start.label}`
      : segment.end.label;
    ticks.push({ label, distance: index + 1 });
  });
  return { segments, ticks, totalDistance: path.length };
}

function svgElement(name, attributes = {}) {
  const element = document.createElementNS("http://www.w3.org/2000/svg", name);
  Object.entries(attributes).forEach(([key, value]) => element.setAttribute(key, value));
  return element;
}

function scheduleBandPlot() {
  if (bandFrame !== null) cancelAnimationFrame(bandFrame);
  bandFrame = requestAnimationFrame(() => { bandFrame = null; renderBandPlot(); });
}

function renderBandPlot() {
  if (!latestResult) return;
  try {
    const path = parseBandPath();
    const samplesPerSegment = Math.max(8, Math.min(240, Number(bandPointsInput.value) || 100));
    const sampled = sampledBandPath(path, samplesPerSegment);
    const tree = latestResult.result.hamiltonian.expression_tree;
    const names = momentumNames();
    const bands = Array.from(
      { length: tree.length },
      () => Array.from({ length: sampled.segments.length }, () => []),
    );
    sampled.segments.forEach((segment, segmentIndex) => segment.forEach((sample) => {
      const environment = { ...bandParameterValues };
      names.forEach((name, index) => { environment[name] = 2 * Math.PI * sample.coordinates[index]; });
      const matrix = tree.map((row) => row.map((entry) => evaluateTree(entry, environment)));
      hermitianEigenvalues(matrix).forEach((energy, band) => bands[band][segmentIndex].push({ distance: sample.distance, energy }));
    }));
    const energies = bands.flat(2).map((point) => point.energy);
    if (!energies.length || energies.some((value) => !Number.isFinite(value))) throw new Error("Band evaluation produced a non-finite value.");
    let minimum = Math.min(...energies);
    let maximum = Math.max(...energies);
    if (Math.abs(maximum - minimum) < 1e-10) { minimum -= 1; maximum += 1; }
    const padding = 0.08 * (maximum - minimum);
    minimum -= padding;
    maximum += padding;
    const width = 760;
    const height = 430;
    const margin = { left: 58, right: 18, top: 20, bottom: 42 };
    const plotWidth = width - margin.left - margin.right;
    const plotHeight = height - margin.top - margin.bottom;
    const x = (distance) => margin.left + plotWidth * distance / sampled.totalDistance;
    const y = (energy) => margin.top + plotHeight * (maximum - energy) / (maximum - minimum);
    bandPlot.replaceChildren();
    for (let index = 0; index <= 4; index += 1) {
      const energy = minimum + (maximum - minimum) * index / 4;
      const line = svgElement("line", { x1: margin.left, x2: width - margin.right, y1: y(energy), y2: y(energy), class: Math.abs(energy) < (maximum - minimum) / 100 ? "band-zero-line" : "band-grid-line" });
      const label = svgElement("text", { x: margin.left - 8, y: y(energy) + 3, "text-anchor": "end", class: "band-energy-label" });
      label.textContent = energy.toFixed(2);
      bandPlot.append(line, label);
    }
    sampled.ticks.forEach((tick) => {
      const line = svgElement("line", { x1: x(tick.distance), x2: x(tick.distance), y1: margin.top, y2: height - margin.bottom, class: "band-axis" });
      const label = svgElement("text", { x: x(tick.distance), y: height - 17, "text-anchor": "middle", class: "band-label" });
      label.textContent = tick.label;
      bandPlot.append(line, label);
    });
    bands.forEach((bandSegments, bandIndex) => {
      bandSegments.forEach((band) => {
        const polyline = svgElement("polyline", {
          class: `band-line band-color-${bandIndex % 3}`,
          points: band.map((point) => `${x(point.distance).toFixed(2)},${y(point.energy).toFixed(2)}`).join(" "),
        });
        bandPlot.append(polyline);
      });
    });
    const energyLabel = svgElement("text", { x: 15, y: height / 2, transform: `rotate(-90 15 ${height / 2})`, "text-anchor": "middle", class: "band-label" });
    energyLabel.textContent = "E";
    bandPlot.append(energyLabel);
    bandError.hidden = true;
  } catch (error) {
    bandError.textContent = error.message || t("bandPathError");
    bandError.hidden = false;
  }
}

async function solveModel(event) {
  event.preventDefault();
  hideError();
  let payload;
  try {
    payload = inputMode === "catalog" ? collectCatalogRequest() : collectRequest();
  } catch (error) {
    showError(error.message);
    return;
  }
  setBusy(true);
  const started = performance.now();
  try {
    const endpoint = inputMode === "catalog" ? "api/catalog/solve" : "api/solve";
    const response = await fetch(endpoint, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(payload),
    });
    const data = await response.json();
    if (!response.ok) throw new Error(data.detail || t("connectionError"));
    const elapsed = performance.now() - started;
    latestResult = { result: data, elapsed };
    renderResult(data, elapsed);
  } catch (error) {
    showError(error.message || t("connectionError"));
  } finally {
    setBusy(false);
  }
}

async function checkHealth() {
  if (window.location.protocol === "file:") {
    runtimeState = { status: "failed", version: "", sympy: "" };
    catalogAvailable = false;
    catalogModeButton.disabled = true;
    setInputMode("manual");
    showError(t("staticPreview"));
    renderRuntime();
    return;
  }
  try {
    const response = await fetch("api/health");
    if (!response.ok) throw new Error("health");
    const data = await response.json();
    runtimeState = { status: "ready", version: data.magnetickp_version, sympy: data.sympy_version };
    footerVersion.textContent = data.magnetickp_version;
    catalogAvailable = Boolean(data.msgcorep?.available);
    if (catalogAvailable) {
      await loadCatalog();
    } else {
      catalogModeButton.disabled = true;
      setInputMode("manual");
      showError(data.msgcorep?.detail || t("catalogUnavailable"));
    }
  } catch (_) {
    runtimeState = { status: "failed", version: "", sympy: "" };
  }
  renderRuntime();
}

languageToggle.addEventListener("click", () => setLanguage(currentLanguage === "zh" ? "en" : "zh"));
inputModeTabs.addEventListener("click", (event) => {
  const button = event.target.closest(".input-mode-tab");
  if (!button || button.disabled) return;
  setInputMode(button.dataset.mode);
});
catalogSpaceGroup.addEventListener("change", async () => {
  try { await loadMagneticGroups("", "Γ"); } catch (error) { showError(error.message); }
});
catalogMagneticGroup.addEventListener("change", async () => {
  try { await loadKPoints("Γ"); } catch (error) { showError(error.message); }
});
catalogBzVariant.addEventListener("change", async () => {
  try { await loadKPoints("Γ"); } catch (error) { showError(error.message); }
});
catalogKPoint.addEventListener("change", async () => {
  renderCatalogCoordinate();
  try { await loadCorepresentations(); } catch (error) { showError(error.message); }
});
catalogRepresentationFamily.addEventListener("change", renderCorepresentations);
catalogCoordinateBasis.addEventListener("change", async () => {
  try { await loadCorepresentations(); } catch (error) { showError(error.message); }
});
convertCatalogToManualButton.addEventListener("click", convertCatalogToManual);
addDiscreteButton.addEventListener("click", () => createOperation({ kind: "unitary", label: `g${operationSerial + 1}`, representation: identityMatrix(dimensionValue() || 2), momentum_image: momentumNames() }));
addContinuousButton.addEventListener("click", () => {
  if (document.querySelector('.operation-card[data-kind="continuous"]')) { showError(t("continuousExists")); return; }
  createOperation({ kind: "continuous", label: "CInfinity", representation: identityMatrix(dimensionValue() || 2), momentum_image: momentumNames(), parameter: "theta" });
});
resetButton.addEventListener("click", async () => {
  if (inputMode === "manual") {
    applyExample(examples[0]);
    return;
  }
  catalogRepresentationFamily.value = "double";
  catalogCoordinateBasis.value = "cartesian";
  orderInput.value = 1;
  methodInput.value = "automatic";
  maxGroupOrder.value = 512;
  if (catalogLoaded) {
    catalogSpaceGroup.value = "92";
    try { await loadMagneticGroups("92.112", "Γ"); } catch (error) { showError(error.message); }
  }
});
hilbertDimension.addEventListener("change", () => {
  if (!dimensionValue()) { showError(t("invalidDimension")); return; }
  document.querySelectorAll(".operation-card").forEach((card) => renderMatrix(card, readMatrix(card)));
  hideError();
});
spatialDimension.addEventListener("change", () => {
  syncSpatialDimensionTabs();
  document.querySelectorAll(".operation-card").forEach((card) => renderMomentum(card, readMomentum(card)));
});
spatialDimensionTabs.addEventListener("click", (event) => {
  const button = event.target.closest(".dimension-tab");
  if (!button) return;
  spatialDimension.value = button.dataset.dimension;
  spatialDimension.dispatchEvent(new Event("change"));
});
solverForm.addEventListener("submit", solveModel);
formulaViewButton.addEventListener("click", () => setHamiltonianView("formula"));
texViewButton.addEventListener("click", () => setHamiltonianView("tex"));
copyLatexButton.addEventListener("click", async () => {
  if (!latestResult) return;
  await navigator.clipboard.writeText(hamiltonianLatex(latestResult.result));
  const original = t("copyLatex");
  copyLatexButton.textContent = t("copied");
  window.setTimeout(() => { copyLatexButton.textContent = original; }, 1200);
});
downloadButton.addEventListener("click", () => {
  if (!latestResult) return;
  const blob = new Blob([JSON.stringify(latestResult.result, null, 2)], { type: "application/json" });
  const url = URL.createObjectURL(blob);
  const link = document.createElement("a");
  link.href = url;
  link.download = "magnetickp-result.json";
  link.click();
  URL.revokeObjectURL(url);
});
bandPathInput.addEventListener("change", scheduleBandPlot);
bandPointsInput.addEventListener("change", scheduleBandPlot);
resetBandButton.addEventListener("click", () => { bandPathInput.value = generatedBandPath || defaultBandPath(); scheduleBandPlot(); });

renderSpatialDimensionTabs();
applyExample(examples[0], false);
setInputMode("catalog");
setLanguage(currentLanguage);
checkHealth();
