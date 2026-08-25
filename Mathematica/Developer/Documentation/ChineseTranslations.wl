(* Exact user-facing translations used by GenerateChineseDocumentationSources.wls.
   Code strings, symbol names, option values, contexts, and paclet URIs are
   intentionally left unchanged. *)

<|
  "History" -> "历史",
  "Categorization" -> "分类",
  "Keywords" -> "关键词",
  "New in: " -> "引入版本：",
  " | Modified in: " -> " | 修改版本：",
  " | Obsolete in: " -> " | 废弃版本：",
  "Examples Initialization" -> "范例初始化",
  "Basic Examples" -> "基本范例",
  "Related Functions" -> "相关函数",
  "Related Guides" -> "相关指南",
  "Related Tutorials" -> "相关教程",
  "Tutorials" -> "教程",
  "Getting Started with MagneticKP" -> "MagneticKP 入门",
  "MagneticKP constructs symmetry-allowed k.p Hamiltonians, connects to representation packages, plots model bands, and fits model parameters to VASP data." ->
    "MagneticKP 用于构造满足对称性要求的 k.p 哈密顿量，可连接表示数据库软件包、绘制模型能带，并根据 VASP 数据拟合模型参数。",
  "MagneticKP constructs symmetry-allowed k.p Hamiltonians" ->
    "MagneticKP 用于构造满足对称性要求的 k.p 哈密顿量",
  "Constructing k.p models" -> "构造 k.p 模型",
  "Bands and parameter fitting" -> "能带与参数拟合",
  "Utilities" -> "实用函数",
  "construct a symmetry-allowed k.p Hamiltonian." -> "构造满足对称性要求的 k.p 哈密顿量。",
  "prepare kpHam input from SpaceGroupIrep or MSGCorep." -> "根据 SpaceGroupIrep 或 MSGCorep 数据生成 kpHam 输入。",
  "call the default construction method directly." -> "直接调用默认构造方法。",
  "call the direct-product construction method directly." -> "直接调用直积构造方法。",
  "use the legacy homogeneous-order interface." -> "使用旧版齐次阶数接口。",
  "explore a model with interactive parameter controls." -> "通过交互式参数控件查看模型能带。",
  "plot a numerical model along a labeled momentum path." -> "沿带标签的动量路径绘制数值模型能带。",
  "import selected bands from a VASP EIGENVAL file." -> "从 VASP EIGENVAL 文件读取指定能带。",
  "compare model bands with imported VASP data." -> "比较模型能带和导入的 VASP 数据。",
  "fit model parameters to selected band data." -> "根据选定的能带数据拟合模型参数。",
  "find a generating subset of a finite group." -> "寻找有限群的一个生成子集。",
  "form a block-diagonal direct sum." -> "构造分块对角直和矩阵。",
  "represent the x momentum component." -> "表示动量的 x 分量。",
  "represent the y momentum component." -> "表示动量的 y 分量。",
  "represent the z momentum component." -> "表示动量的 z 分量。",

  "constructs the symmetry-allowed k.p Hamiltonian through korder and returns its data in an Association." ->
    "构造截至 korder 阶、满足对称性要求的 k.p 哈密顿量，并以 Association 返回结果。",
  "constructs only the listed homogeneous orders." -> "仅构造列表中指定的齐次阶数。",
  "input is an Association whose \"Unitary\" and optional \"Antiunitary\" entries associate symmetry labels with representation matrices and transformed momentum vectors." ->
    "input 是一个 Association，其中 \"Unitary\" 和可选的 \"Antiunitary\" 将对称操作标签与表示矩阵及变换后的动量向量关联起来。",
  "The result contains \"ham\", \"korder\", \"dim\", and \"NumberOfParameters\". Use result[\"ham\"] to obtain the Hamiltonian matrix." ->
    "返回结果包含 \"ham\"、\"korder\"、\"dim\" 和 \"NumberOfParameters\"。使用 result[\"ham\"] 取得哈密顿量矩阵。",
  "The historical input key \"Anitunitary\" remains accepted, but new input should use \"Antiunitary\"." ->
    "仍兼容历史拼写 \"Anitunitary\"，但新输入应使用 \"Antiunitary\"。",
  "selects \"IterativeSimplification\", \"DirectProductDecomposition\", \"Reynolds\", \"Automatic\", or the compatibility method \"IterativeSimplification2016\"." ->
    "选择 \"IterativeSimplification\"、\"DirectProductDecomposition\"、\"Reynolds\"、\"Automatic\" 或兼容方法 \"IterativeSimplification2016\"。",
  "specifies whether the supplied symmetry data contains generators or all elements; Automatic uses metadata in input when available." ->
    "指定输入的对称数据是生成元还是全部群元；Automatic 会优先读取 input 中的元数据。",
  "sets the largest finite group order accepted while completing generator input." ->
    "设置补全生成元输入时允许的最大有限群阶数。",
  "specifies whether the returned Association also includes per-order diagnostic information." ->
    "指定返回的 Association 是否同时包含逐阶诊断信息。",
  "Construct a first-order two-band model from explicitly supplied symmetry data:" ->
    "根据手动给出的对称数据构造一阶双带模型：",

  "converts SpaceGroupIrep or MSGCorep representation data into input accepted by kpHam." ->
    "将 SpaceGroupIrep 或 MSGCorep 的表示数据转换为 kpHam 可接受的输入。",
  "Use an integer group number with SpaceGroupIrep and a BNS pair {spaceGroup, serial} with MSGCorep." ->
    "SpaceGroupIrep 使用整数空间群编号；MSGCorep 使用 BNS 编号对 {spaceGroup, serial}。",
  "k may be a high-symmetry-point label or exact coordinates. rep may select one representation or a list whose direct sum is required." ->
    "k 可以是高对称点标签或精确坐标。rep 可以选择一个表示，也可以给出需要取直和的表示列表。",
  "For MSGCorep, each integer in rep selects the row with that Index in showMLGCorep; here {2} selects B_2." ->
    "对于 MSGCorep，rep 中的每个整数对应 showMLGCorep 表格第一列的 Index；此处 {2} 选择 B_2。",
  "Load SpaceGroupIrep or MSGCorep before calling interfaceRep." ->
    "调用 interfaceRep 前需要先载入 SpaceGroupIrep 或 MSGCorep。",
  "The returned Association records the selected k-point label and coordinate system so highSymmetryKPath can label a local path automatically." ->
    "返回的 Association 会记录所选 k 点标签和坐标系，使 highSymmetryKPath 可以自动标记局域路径。",
  "specifies whether symmetry operations are converted to Cartesian momentum coordinates." ->
    "指定是否把对称操作转换到笛卡尔动量坐标。",
  "specifies whether the returned input contains a generating subset instead of every group element." ->
    "指定返回的输入是生成子集还是全部群元。",
  "Obtain input for the second corepresentation at the B point of magnetic space group 13.65:" ->
    "取得磁空间群 13.65 的 B 点第二个共表示所对应的输入：",
  "List the corepresentations at the B point of magnetic space group 13.65:" ->
    "列出磁空间群 13.65 在 B 点的共表示：",
  "The second data row is labeled B_2, so select the second corepresentation with {2}. interfaceRep returns the complete input Association:" ->
    "表格中第二个数据行标记为 B_2，因此用 {2} 选择第二个共表示。interfaceRep 会直接返回完整的输入 Association：",
  "Pass the returned Association directly to kpHam to construct the first-order Hamiltonian:" ->
    "将返回的 Association 直接传给 kpHam，构造一阶哈密顿量：",

  "generate a local momentum path from little-co-group subgroups." ->
    "根据小余群的子群生成局域动量路径。",
  "returns a local labeled momentum path obtained from the complete little magnetic co-group contained in interfaceRep input." ->
    "根据 interfaceRep 输入中的完整磁小余群返回一条带标签的局域动量路径。",
  "input is normally the Association returned by interfaceRep. A manually prepared input with three-dimensional momentum transformations is also accepted." ->
    "input 通常使用 interfaceRep 返回的 Association；也可以使用手动准备的、包含三维动量变换的输入。",
  "The default result is an ordinary path list. highSymmetryKPath does not call bandplot; pass the returned path to bandplot or bandManipulate in a separate step." ->
    "默认结果是普通的 path 列表。highSymmetryKPath 不会调用 bandplot；请在下一步把返回的 path 单独传给 bandplot 或 bandManipulate。",
  "The path is centered at the local k.p expansion origin {0, 0, 0}. Each retained line has the complete largest magnetic subgroup that fixes a generic point on that line, and symmetry-equivalent lines are represented once by default." ->
    "路径以局域 k.p 展开原点 {0, 0, 0} 为中心。每条保留的直线都对应固定该线上一般点的完整最大磁子群；默认只保留每类对称等价直线的一个代表。",
  "Unitary and antiunitary little-co-group elements remain distinct even when they act on momentum in the same way. In returned diagnostics, SubgroupOrder is the full magnetic subgroup order, while MomentumActionOrder counts distinct momentum-action matrices." ->
    "即使幺正与反幺正小余群元素对动量的作用相同，程序也会把它们作为不同群元保留。在返回的诊断信息中，SubgroupOrder 是完整磁子群阶数，MomentumActionOrder 是不同动量作用矩阵的数量。",
  "highSymmetryKPath identifies symmetry lines; whether bands stay degenerate along a line depends on the selected representations and Hamiltonian." ->
    "highSymmetryKPath 负责识别对称线；能带是否沿某条线保持简并，取决于所选表示和哈密顿量。",
  "Use \"ReturnData\" -> True to inspect the directions, subgroup orders, little co-group order, and the generated path together." ->
    "使用 \"ReturnData\" -> True 可同时查看方向、子群阶数、小余群阶数和生成的 path。",
  "sets the extent of each path arm in the coordinate system used by input." ->
    "设置每条路径分支在 input 所用坐标系中的延伸范围。",
  "sets the local center of the returned path." ->
    "设置返回路径的局域中心。",
  "uses the k-point label recorded by interfaceRep, or sets an explicit label string." ->
    "使用 interfaceRep 记录的 k 点标签，或设置明确的标签字符串。",
  "returns one representative of each symmetry-equivalent line class; \"All\" keeps every distinct line." ->
    "每类对称等价直线返回一个代表；\"All\" 保留所有不同的直线。",
  "includes both arms of every line through the local center." ->
    "包含每条直线穿过局域中心的两个方向分支。",
  "returns only the path, or an Association containing the path and symmetry diagnostics." ->
    "仅返回 path，或返回同时包含 path 和对称性诊断信息的 Association。",
  "sets the largest little co-group order accepted while completing the supplied transformations." ->
    "设置补全所给变换时允许的最大小余群阶数。",
  "sets the largest number of subgroups enumerated while finding fixed lines." ->
    "设置寻找不变直线时允许枚举的最大子群数量。",
  "Create a separate local path at K from the MSGCorep representation input:" ->
    "根据 MSGCorep 表示输入，在 K 点单独生成局域 path：",
  "Request diagnostics and compare each retained direction with both the full magnetic subgroup order and the number of distinct momentum actions:" ->
    "请求诊断信息，并同时查看每个保留方向的完整磁子群阶数及不同动量作用的数量：",
  "Use the previously generated path as an ordinary argument in a separate bandplot call:" ->
    "把前面生成的 path 作为普通参数，单独调用 bandplot：",
  "Use the original MagneticKP K-point example: select K5 and K6 of magnetic space group 191.234 (P6/mmm1') and inspect the complete K-point magnetic symmetry and its retained lines:" ->
    "使用原版 MagneticKP 的 K 点示例：选择磁空间群 191.234（P6/mmm1'）的 K5 与 K6，并查看 K 点的完整磁对称性和保留的直线：",
  "At K, the unitary little co-group is D3h (order 12), and its antiunitary coset gives a complete magnetic little co-group of order 24. On the first automatic line, q parallel to (0,0,1), the unitary stabilizer is C3v (order 6) and the full magnetic line subgroup has order 12. Construct the first-order four-band Hamiltonian; K5 and K6 each remain twofold degenerate on this line, so the four bands appear as two coincident pairs:" ->
    "在 K 点，幺正小余群是 D3h（12 阶），加入反幺正陪集后得到 24 阶的完整磁小余群。自动找到的第一条直线满足 q 平行于 (0,0,1)，其幺正稳定子群为 C3v（6 阶），完整磁线子群为 12 阶。构造一阶四带哈密顿量后，K5 与 K6 沿这条线各自保持二重简并，因此四条能带显示为两组完全重合的曲线：",
  "Use the same Hamiltonian along a generic direction. Away from the protected line, the two degenerate pairs split into four branches:" ->
    "沿一般方向使用同一个哈密顿量。离开受保护的高对称线后，两组二重简并劈裂为四支：",
  "One complete example: for magnetic space group 191.234 (P6/mmm1') at K, use K5 and K6 to generate the local path and plot the protected q parallel to (0,0,1) line. The unitary K-point little co-group is D3h (order 12), the complete magnetic little co-group has order 24, and the full magnetic line subgroup has order 12. The four bands therefore appear as two coincident, twofold-degenerate curves:" ->
    "只保留一个完整示例：对磁空间群 191.234（P6/mmm1'）的 K 点，使用 K5 与 K6 生成局域路径，并绘制受保护的 q 平行于 (0,0,1) 直线。K 点幺正小余群为 D3h（12 阶），完整磁小余群为 24 阶，完整磁线子群为 12 阶。因此四条能带显示为两条完全重合的二重简并曲线：",
  "Generate the path for K5 and K6 at K of magnetic space group 191.234, construct the first-order Hamiltonian, and pass the returned path directly to bandplot:" ->
    "对磁空间群 191.234 的 K 点选择 K5 与 K6，生成路径、构造一阶哈密顿量，然后把返回的 path 直接传给 bandplot：",
  "Generate and display the path for K5 and K6 at K of magnetic space group 191.234:" ->
    "对磁空间群 191.234 的 K 点选择 K5 与 K6，生成并直接显示 path：",
  "Construct the first-order Hamiltonian and pass the displayed path directly to bandplot:" ->
    "构造一阶哈密顿量，并把上面显示的 path 直接传给 bandplot：",

  "constructs a symmetry-allowed k.p Hamiltonian with the default solver." ->
    "使用默认求解方法构造满足对称性要求的 k.p 哈密顿量。",
  "constructs a symmetry-allowed k.p Hamiltonian with the direct-product solver." ->
    "使用直积求解方法构造满足对称性要求的 k.p 哈密顿量。",
  "An integer korder includes every order from zero through korder; a list selects only the listed homogeneous orders." ->
    "整数 korder 包含从零阶到 korder 阶的所有项；列表只选择其中列出的齐次阶数。",
  "The return format is the same Association used by kpHam." ->
    "返回格式与 kpHam 使用的 Association 相同。",
  "Construct all terms through first order:" -> "构造从零阶到一阶的全部项：",

  "constructs only the homogeneous order-n part of a k.p Hamiltonian using the legacy calling convention." ->
    "按照旧版调用方式，仅构造 k.p 哈密顿量的 n 阶齐次部分。",
  "GETkpHam is retained for notebooks written for earlier MagneticKP releases." ->
    "GETkpHam 用于兼容早期 MagneticKP 版本编写的 Notebook。",
  "Use kpHam for new work, especially when multiple orders or diagnostic output are required." ->
    "新工作应使用 kpHam，特别是需要多个阶数或诊断输出时。",
  "Construct only the first-order contribution:" -> "仅构造一阶贡献：",

  "creates an interactive band plot with controls for the free parameters in ham." ->
    "创建交互式能带图，并为 ham 中的自由参数提供控件。",
  "Each path segment has the form {{k1, k2}, {label1, label2}} and n gives the sampling count per segment." ->
    "每段路径的格式为 {{k1, k2}, {label1, label2}}，n 是每段路径的采样数。",
  "The momentum coordinates in path are multiplied by 2 Pi before substitution into ham." ->
    "path 中的动量坐标在代入 ham 前会乘以 2 Pi。",
  "sets the displayed energy range." -> "设置显示的能量范围。",
  "sets a symmetric slider range, or may be given explicitly as {min, max}." ->
    "设置对称的滑块范围，也可以明确给出 {min, max}。",
  "Construct a two-band Hamiltonian with kpHam, then explore it with interactive parameter controls:" ->
    "先用 kpHam 构造双带哈密顿量，再通过交互式参数控件查看其能带：",

  "draws the numerical eigenvalues of ham along a labeled reciprocal-fractional path after applying parameter rules." ->
    "应用参数规则后，沿带标签的倒易分数坐标路径绘制 ham 的数值本征值。",
  "Each path segment has the form {{k1, k2}, {label1, label2}}. The two momentum points contain three reciprocal-fractional coordinates." ->
    "每段路径的格式为 {{k1, k2}, {label1, label2}}，两个动量点均包含三个倒易分数坐标。",
  "n is the number of subdivisions per segment, so each segment contains n + 1 sampled points. Momentum coordinates are multiplied by 2 Pi before substitution into ham." ->
    "n 是每段路径的细分数，因此每段包含 n + 1 个采样点。动量坐标乘以 2 Pi 后再代入 ham。",
  "rules must supply numerical values for every free parameter in ham." ->
    "rules 必须为 ham 中的每个自由参数提供数值。",
  "At a discontinuity between two consecutive segments, the boundary tick combines the two labels as end|start." ->
    "相邻两段路径不连续时，边界刻度会按 end|start 合并两个标签。",
  "The result is a framed Graphics object. Invalid paths, sampling counts, rules, tick specifications, font settings, or nonnumeric Hamiltonians return $Failed with a message." ->
    "结果为带边框的 Graphics 对象。路径、采样数、参数规则、刻度、字体设置无效或哈密顿量不能数值化时，会返回 $Failed 并给出消息。",
  "sets the vertical frame ticks to Automatic, None, numeric positions, or standard tick specifications." ->
    "将纵轴边框刻度设置为 Automatic、None、数值位置或标准刻度规格。",
  "sets the size of path labels and frame ticks." ->
    "设置路径标签和边框刻度的字号。",
  "sets the font family of path labels and frame ticks." ->
    "设置路径标签和边框刻度的字体。",
  "sets the final Graphics image size." ->
    "设置最终 Graphics 的图像尺寸。",
  "Construct a first-order two-band Hamiltonian with kpHam and assign a numerical value to every model parameter:" ->
    "使用 kpHam 构造一阶双带哈密顿量，并为每个模型参数指定数值：",
  "Define a two-segment path and draw the bands with explicit energy ticks and display settings:" ->
    "定义包含两段的路径，并使用明确的能量刻度和显示设置绘制能带：",
  "\"PlotRange\" - restrict the displayed energy interval:" ->
    "\"PlotRange\"——限制图中显示的能量区间：",
  "\"YTicks\" - place explicit labeled ticks on the energy axis:" ->
    "\"YTicks\"——在能量轴上设置明确的带标签刻度：",
  "FontSize - enlarge the path labels and frame ticks:" ->
    "FontSize——放大路径标签和边框刻度：",
  "FontFamily - use Helvetica for the path labels and frame ticks:" ->
    "FontFamily——路径标签和边框刻度使用 Helvetica 字体：",
  "ImageSize - set the final plot width:" ->
    "ImageSize——设置最终图形宽度：",

  "imports selected bands from a VASP EIGENVAL file and shifts their energies by efermi." ->
    "从 VASP EIGENVAL 文件导入指定能带，并用 efermi 平移能量。",
  "spin selects the energy column to import, and startband through endband select an inclusive band range." ->
    "spin 选择需要导入的能量列，startband 到 endband 指定包含端点的能带范围。",
  "The result is a list of {k, energies} pairs suitable for bandManipulateWithVASP and fittingKP." ->
    "结果为 {k, energies} 对的列表，可直接交给 bandManipulateWithVASP 和 fittingKP。",
  "Import bands 15 through 18 from an EIGENVAL file selected by the user:" ->
    "从用户选择的 EIGENVAL 文件导入第 15 到第 18 条能带：",
  "Import the first two bands from the EIGENVAL sample stored beside this help page:" ->
    "从此帮助页面随附的 EIGENVAL 样例中导入前两条能带：",

  "interactively compares the eigenvalues of ham with imported VASP eigenvalues." ->
    "交互比较 ham 的本征值与导入的 VASP 本征值。",
  "eigdata should be the result returned by vaspEig or data in the same {k, energies} form." ->
    "eigdata 应为 vaspEig 的返回结果，或采用相同 {k, energies} 格式的数据。",
  "The model bands are drawn in black and imported bands in red." ->
    "模型能带显示为黑色，导入的能带显示为红色。",
  "When kpHam is constructed in fractional coordinates about a nonzero expansion point, shift its momentum variables to the absolute fractional coordinates stored in EIGENVAL before comparing the bands." ->
    "如果 kpHam 使用非零展开点附近的分数坐标构造，比较能带前必须平移动量变量，使其对应 EIGENVAL 中保存的绝对分数坐标。",
  "Inspect the corepresentations at K used by the original MagneticKP fitting example:" ->
    "查看 MagneticKP 原始拟合示例在 K 点使用的共表示：",
  "Select corepresentations 2 and 3 in fractional coordinates and construct the second-order k.p Hamiltonian:" ->
    "在分数坐标下选择第 2、3 个共表示，并构造二阶 k.p 哈密顿量：",
  "Read bands 21 and 22 from the original VASP EIGENVAL file using the Fermi energy from the original example:" ->
    "使用原始示例给出的费米能，从 VASP 原始 EIGENVAL 文件读取第 21、22 条能带：",
  "Shift the local k.p expansion from K to the absolute fractional coordinates in EIGENVAL, then tune the model against the VASP bands:" ->
    "将 K 点附近的局域 k.p 展开平移到 EIGENVAL 的绝对分数坐标，再对照 VASP 能带调整模型参数：",

  "fits the free parameters of ham to selected imported eigenvalues." ->
    "根据选定的导入本征值拟合 ham 中的自由参数。",
  "fits the free parameters of ham to selected reference bands." ->
    "根据选定的参考能带拟合 ham 中的自由参数。",
  "krange is a list of k-point indices used in the fit, and initialRules supplies a starting value for every free parameter." ->
    "krange 是参与拟合的 k 点序号列表，initialRules 为每个自由参数提供初始值。",
  "The result contains \"FittedParams\", \"LSQForInitParams\", and \"LSQForFittedParams\"." ->
    "结果包含 \"FittedParams\"、\"LSQForInitParams\" 和 \"LSQForFittedParams\"。",
  "eigdata uses the {{kx, ky, kz}, {energy1, energy2, ...}} records returned by vaspEig. krange is a nonempty list of 1-based record indices." ->
    "eigdata 使用 vaspEig 返回的 {{kx, ky, kz}, {energy1, energy2, ...}} 记录格式。krange 是从 1 开始计数的非空记录序号列表。",
  "initialRules may be a rule list or Association and must give a real numeric starting value for every free parameter in ham." ->
    "initialRules 可以是规则列表或 Association，并且必须为 ham 中的每个自由参数给出实数初值。",
  "fittingKP substitutes the coordinates in eigdata directly into ham; it does not multiply them by 2 Pi. Shift a local k.p Hamiltonian to the absolute fractional coordinates stored in EIGENVAL before fitting." ->
    "fittingKP 会把 eigdata 中的坐标直接代入 ham，不会再乘以 2 Pi。拟合前，应将局域 k.p 哈密顿量平移到 EIGENVAL 中保存的绝对分数坐标。",
  "At each selected k point, sorted model eigenvalue j is compared with reference energy j. Inspect \"ComparisonPlot\" after fitting, especially near band crossings." ->
    "在每个选定的 k 点，排序后的第 j 个模型本征值会与第 j 个参考能量比较。拟合后应查看 \"ComparisonPlot\"，尤其需要留意能带交叉附近。",
  "Use result[\"FittedParams\"] as the parameter rules for bandplot. \"LSQForInitParams\" and \"LSQForFittedParams\" report the sum of squared band residuals before and after fitting." ->
    "把 result[\"FittedParams\"] 作为 bandplot 的参数规则。\"LSQForInitParams\" 和 \"LSQForFittedParams\" 分别给出拟合前后的能带残差平方和。",
  "The result also contains \"ComparisonPlot\", \"Converged\", \"Iterations\", the effective data selection, and the number of residuals used." ->
    "结果还包含 \"ComparisonPlot\"、\"Converged\"、\"Iterations\"、实际使用的数据范围和残差数量。",
  "sets the maximum number of fitting iterations." ->
    "设置拟合的最大迭代次数。",
  "sets the fitting time limit in seconds, or may be Infinity." ->
    "设置拟合的时间上限（秒），也可以设为 Infinity。",
  "uses every selected reference energy, or only energies in the inclusive interval {emin, emax}." ->
    "使用全部选定的参考能量，或只使用闭区间 {emin, emax} 内的能量。",
  "uses every index in krange, or selects a Radius or componentwise Range around a fractional-coordinate Center." ->
    "使用 krange 中的全部序号，或在分数坐标 Center 周围按 Radius 或逐分量 Range 选择 k 点。",
  "uses every reference band position, or selects 1-based positions with an integer, increasing list, or positive-step Span." ->
    "使用全部参考能带位置，或用整数、递增列表、正步长 Span 选择从 1 开始计数的位置。",
  "uses equal weights, or a nonnegative record-by-band matrix; a zero entry excludes that residual." ->
    "使用相同权重，或给出非负的“记录数×能带数”权重矩阵；零权重会排除相应残差。",
  "sets the relative numerical-difference step used during fitting." ->
    "设置拟合中数值差分使用的相对步长。",
  "sets the initial damping used by the fitter." ->
    "设置拟合器使用的初始阻尼。",
  "sets the stopping tolerance for the fit." ->
    "设置拟合的停止容差。",
  "Inspect the K-point corepresentations used by the original MagneticKP fitting example:" ->
    "查看 MagneticKP 原始拟合示例使用的 K 点共表示：",
  "Construct the second-order model in fractional coordinates, shift its local K-point momenta to the absolute coordinates in EIGENVAL, and prepare starting values for every generated parameter:" ->
    "在分数坐标下构造二阶模型，把 K 点附近的局域动量平移到 EIGENVAL 中的绝对坐标，并为所有生成参数准备初值：",
  "Read VASP bands 21 and 22 from the supplied original EIGENVAL file:" ->
    "从随附的原始 EIGENVAL 文件读取 VASP 第 21、22 条能带：",
  "Fit five points near K and inspect the fitted rules, error reduction, and convergence status:" ->
    "拟合 K 点附近的五个点，并查看拟合规则、误差下降和收敛状态：",
  "Compare the selected VASP points with the initial and fitted k.p bands:" ->
    "比较选定的 VASP 数据点、初始 k.p 能带和拟合后的 k.p 能带：",
  "Fit a model after importing eigdata with vaspEig:" ->
    "使用 vaspEig 导入 eigdata 后拟合模型：",
  "Fit a one-parameter two-band model to a small data set:" ->
    "将单参数双带模型拟合到一个小型数据集：",

  "returns a generating subset for a finite group described by its elements and multiplication function." ->
    "根据有限群的全部元素和乘法函数返回一个生成子集。",
  "elements should contain the complete finite group and identity should be its identity element." ->
    "elements 应包含完整的有限群，identity 应为其单位元。",
  "The returned subset generates the group, but it is not guaranteed to have the minimum possible size." ->
    "返回的子集能够生成该群，但不保证其大小为最小。",
  "Find a generator for integers modulo four under addition:" ->
    "求模四整数加法群的生成元：",

  "returns the block-diagonal direct sum of matrices a and b." ->
    "返回矩阵 a 和 b 的分块对角直和。",
  "directSum is useful when combining representation matrices for several irreducible representations." ->
    "组合多个不可约表示的表示矩阵时可以使用 directSum。",
  "Form the direct sum of a two-dimensional and a one-dimensional matrix:" ->
    "构造一个二维矩阵与一个一维矩阵的直和：",

  "is the x component of momentum used in MagneticKP Hamiltonians." ->
    "是 MagneticKP 哈密顿量中使用的动量 x 分量。",
  "is the y component of momentum used in MagneticKP Hamiltonians." ->
    "是 MagneticKP 哈密顿量中使用的动量 y 分量。",
  "is the z component of momentum used in MagneticKP Hamiltonians." ->
    "是 MagneticKP 哈密顿量中使用的动量 z 分量。",
  "MagneticKP uses kx, ky, and kz as symbolic momentum components in symmetry data and Hamiltonians." ->
    "MagneticKP 在对称数据和哈密顿量中使用 kx、ky 和 kz 作为符号动量分量。",
  "Form a rotationally symmetric quadratic momentum expression:" ->
    "构造旋转对称的二次动量表达式：",

  "This tutorial shows the shortest complete workflow: load MagneticKP, enter symmetry data, construct a Hamiltonian, and inspect its bands." ->
    "本教程给出最短的完整流程：载入 MagneticKP、输入对称数据、构造哈密顿量并查看能带。",
  "This tutorial shows the shortest complete workflow" ->
    "本教程给出最短的完整流程",
  "Load MagneticKP" -> "载入 MagneticKP",
  "Load the package after installing the MagneticKP paclet:" ->
    "安装 MagneticKP Paclet 后载入软件包：",
  "Enter symmetry data" -> "输入对称数据",
  "For each symmetry operation, supply its representation matrix and its action on momentum. Put unitary operations under \"Unitary\" and add \"Antiunitary\" only when needed:" ->
    "为每个对称操作给出表示矩阵及其对动量的作用。幺正操作放在 \"Unitary\" 中；只有需要时才添加 \"Antiunitary\"：",
  "Construct the Hamiltonian" -> "构造哈密顿量",
  "Use an integer order to include every order from zero through that value. The returned Association keeps the Hamiltonian and basic summary data together:" ->
    "使用整数阶数时会包含从零阶到该阶的全部项。返回的 Association 同时保存哈密顿量和基本摘要信息：",
  "Extract the matrix from the \"ham\" entry:" ->
    "从 \"ham\" 条目中取得矩阵：",
  "Inspect model bands" -> "查看模型能带",
  "Define a labeled momentum path and use bandManipulate while choosing parameter values:" ->
    "定义带标签的动量路径，并使用 bandManipulate 调整参数值：",
  "Use representation databases" -> "使用表示数据库",
  "When representation matrices are available from MSGCorep or SpaceGroupIrep, interfaceRep prepares the same input format automatically:" ->
    "如果表示矩阵可由 MSGCorep 或 SpaceGroupIrep 获得，interfaceRep 会自动生成相同格式的输入：",
  "Next steps" -> "下一步",
  "See " -> "参阅 ",
  " for all construction options, " -> " 了解全部构造选项；参阅 ",
  " for a fixed publication-style band plot, and " -> " 绘制固定参数的能带图；结合 ",
  " with " -> " 与 ",
  " for fitting to VASP data." -> " 拟合 VASP 数据。",

  "k.p Hamiltonian" -> "k.p 哈密顿量",
  "symmetry" -> "对称性",
  "effective model" -> "有效模型",
  "band theory" -> "能带理论",
  "representation" -> "表示",
  "magnetic space group" -> "磁空间群",
  "symmetry constraints" -> "对称约束",
  "direct product" -> "直积",
  "legacy interface" -> "旧版接口",
  "homogeneous order" -> "齐次阶数",
  "band plot" -> "能带图",
  "interactive" -> "交互式",
  "parameters" -> "参数",
  "momentum path" -> "动量路径",
  "little co-group" -> "小余群",
  "high-symmetry line" -> "高对称线",
  "subgroup" -> "子群",
  "dispersion" -> "色散",
  "band data" -> "能带数据",
  "import" -> "导入",
  "band comparison" -> "能带比较",
  "fitting" -> "拟合",
  "parameter fitting" -> "参数拟合",
  "least squares" -> "最小二乘",
  "finite group" -> "有限群",
  "generator" -> "生成元",
  "group elements" -> "群元",
  "direct sum" -> "直和",
  "block diagonal" -> "分块对角",
  "representation matrix" -> "表示矩阵",
  "momentum" -> "动量",
  "x component" -> "x 分量",
  "y component" -> "y 分量",
  "z component" -> "z 分量",
  "getting started" -> "入门",
  "symmetry input" -> "对称输入"
|>
