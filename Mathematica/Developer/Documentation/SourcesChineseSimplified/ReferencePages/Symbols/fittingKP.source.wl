Notebook[{Cell[CellGroupData[{Cell["fittingKP", "ObjectName"],
     Cell[TextData[{Cell["   ", "ModInfo"], Cell[BoxData["fittingKP[ham, eigdata, krange, initialRules]"],
         "InlineFormula"], "  ", "根据选定的参考能带拟合 ham 中的自由参数。"}], "Usage"],
     Cell[
      "eigdata 使用 vaspEig 返回的 {{kx, ky, kz}, {energy1, energy2, ...}} 记录格式。krange 是从 1 开始计数的非空记录序号列表。"\
, "Notes"],
     Cell[
      "initialRules 可以是规则列表或 Association，并且必须为 ham 中的每个自由参数给出实数初值。"\
, "Notes"],
     Cell[
      "fittingKP 会把 eigdata 中的坐标直接代入 ham，不会再乘以 2 Pi。拟合前，应将局域 k.p 哈密顿量平移到 EIGENVAL 中保存的绝对分数坐标。"\
, "Notes"],
     Cell[
      "在每个选定的 k 点，排序后的第 j 个模型本征值会与第 j 个参考能量比较。拟合后应查看 \"ComparisonPlot\"，尤其需要留意能带交叉附近。"\
, "Notes"], Cell["把 result[\"FittedParams\"] 作为 bandplot 的参数规则。\"LSQForInitParams\" 和 \"LSQForFittedParams\" \
分别给出拟合前后的能带残差平方和。", "Notes"],
     Cell[
      "结果还包含 \"ComparisonPlot\"、\"Converged\"、\"Iterations\"、实际使用的数据范围和残差数量。"\
, "Notes"], Cell[BoxData[GridBox[{{Cell["      ", "TableRowIcon"], "MaxIterations", "100",
          Cell["设置拟合的最大迭代次数。", "TableText"]}, {Cell["      ", "TableRowIcon"],
          "TimeConstraint", "60", Cell["设置拟合的时间上限（秒），也可以设为 Infinity。",
           "TableText"]}, {Cell["      ", "TableRowIcon"], "\"EnergyWindow\"", "All",
          Cell["使用全部选定的参考能量，或只使用闭区间 {emin, emax} 内的能量。", "TableText"]},
         {Cell["      ", "TableRowIcon"], "\"KPointNeighborhood\"", "All",
          Cell[
           "使用 krange 中的全部序号，或在分数坐标 Center 周围按 Radius 或逐分量 Range 选择 k 点。"\
, "TableText"]}, {Cell["      ", "TableRowIcon"], "\"BandSelection\"", "All",
          Cell[
           "使用全部参考能带位置，或用整数、递增列表、正步长 Span 选择从 1 开始计数的位置。"\
, "TableText"]}, {Cell["      ", "TableRowIcon"], "\"ResidualWeights\"", "Automatic",
          Cell[
           "使用相同权重，或给出非负的“记录数×能带数”权重矩阵；零权重会排除相应残差。"\
, "TableText"]}, {Cell["      ", "TableRowIcon"], "\"FiniteDifferenceStep\"", "10^-5",
          Cell["设置拟合中数值差分使用的相对步长。", "TableText"]}, {Cell["      ", "TableRowIcon"],
          "\"InitialDamping\"", "10^-3", Cell["设置拟合器使用的初始阻尼。", "TableText"]},
         {Cell["      ", "TableRowIcon"], "\"FitTolerance\"", "10^-8", Cell["设置拟合的停止容差。",
           "TableText"]}}, GridBoxAlignment -> {"Columns" -> {{Left}}}]], "3ColumnTableMod"]}, Open]],
  Cell[CellGroupData[{Cell["范例初始化", "ExamplesInitializationSection"],
     Cell[BoxData["CompoundExpression[Needs[\"MagneticKP`\"], Needs[\"MSGCorep`\"]]"], "ExampleInitialization"]},
    Open]], Cell[CellGroupData[{Cell["基本范例", "PrimaryExamplesSection"],
     Cell["查看 MagneticKP 原始拟合示例使用的 K 点共表示：", "ExampleText"],
     Cell[BoxData["showMLGCorep[{156, 50}, {1/3, 1/3, 0}]"], "Input"],
     Cell[BoxData[TagBox[GridBox[{{TagBox[GridBox[{{TemplateBox[{" ", "\" \"", "\"MSG \"", StyleBox[
                  TemplateBox[{"156", "\".\"", "50"}, "RowDefault"], Bold, RGBColor[0, 0, 1], StripOnInput -> False],
                 "\"(\"", StyleBox[TemplateBox[{"\"P\"", "\"3\"", "\"m\"", "\"1\"", TemplateBox[{"\"1\"", "\"′\""},
                     "Superscript", SyntaxForm -> SuperscriptBox]}, "RowDefault"], Bold, RGBColor[0, 0, 1],
                  StripOnInput -> False], "\"):  k-point name is \"", StyleBox["\"K\"", Bold, RGBColor[0, 0, 1],
                  StripOnInput -> False], "\", for \"", TemplateBox[{",", "\",\"", "\"HexaPrim\""},
                  "RowWithSeparators"]}, "RowWithSeparators"]}, {TemplateBox[{TemplateBox[{"\"▶  k=(\"",
                   TemplateBox[{",", "\",\"", FractionBox["1", "3"], FractionBox["1", "3"], "0"}, "RowWithSeparators"],
                   "\")\""}, "RowDefault"], "\" non-standard, \"", "\"\\!\\(\\*SubscriptBox[\\(k\\), \\(BC\\)]\\)\"",
                 "\"=(\"", TemplateBox[{",", "\",\"", RowBox[{"-", FractionBox["1", "3"]}], FractionBox["2", "3"],
                   "0"}, "RowWithSeparators"], "\")  [\"", SubscriptBox["\"σ\"", "\"v1\""], "\"]\"",
                 "\"\\!\\(\\*SubscriptBox[\\(k\\), \\(BC\\)]\\)\"", "\" ⇔ k\""}, "RowDefault"]}},
             GridBoxAlignment -> {"Columns" -> {{Left}}}, DefaultBaseStyle -> "Column", GridBoxItemSize ->
              {"Columns" -> {{Automatic}}, "Rows" -> {{Automatic}}}], "Column"]},
          {TagBox[GridBox[{{"\"Index\"", "", "", "1", "2", "3", "4", "5", "6"}, {"\"Element\"", "", "",
               TemplateBox[{"\"{\"", "\"E\"", "\"|\"", TemplateBox[{"0", "0", "0"}, "RowDefault"], "\"}\""},
                "RowDefault"], TemplateBox[{"\"{\"", TemplateBox[{"\"C\"", "\"3\"", "\"+\""}, "Subsuperscript",
                  SyntaxForm -> SubsuperscriptBox], "\"|\"", TemplateBox[{"0", "0", "0"}, "RowDefault"], "\"}\""},
                "RowDefault"], TemplateBox[{"\"{\"", TemplateBox[{"\"C\"", "\"3\"", "\"-\""}, "Subsuperscript",
                  SyntaxForm -> SubsuperscriptBox], "\"|\"", TemplateBox[{"0", "0", "0"}, "RowDefault"], "\"}\""},
                "RowDefault"], StyleBox[TemplateBox[{TemplateBox[{"\"{\"", SubscriptBox["\"σ\"", "\"v1\""], "\"|\"",
                    TemplateBox[{"0", "0", "0"}, "RowDefault"], "\"}\""}, "RowDefault"], "\"′\""}, "Superscript",
                 SyntaxForm -> SuperscriptBox], RGBColor[1, 0, 0], StripOnInput -> False], StyleBox[
                TemplateBox[{TemplateBox[{"\"{\"", SubscriptBox["\"σ\"", "\"v2\""], "\"|\"", TemplateBox[
                     {"0", "0", "0"}, "RowDefault"], "\"}\""}, "RowDefault"], "\"′\""}, "Superscript",
                 SyntaxForm -> SuperscriptBox], RGBColor[1, 0, 0], StripOnInput -> False], StyleBox[
                TemplateBox[{TemplateBox[{"\"{\"", SubscriptBox["\"σ\"", "\"v3\""], "\"|\"", TemplateBox[
                     {"0", "0", "0"}, "RowDefault"], "\"}\""}, "RowDefault"], "\"′\""}, "Superscript",
                 SyntaxForm -> SuperscriptBox], RGBColor[1, 0, 0], StripOnInput -> False]},
              {TagBox[GridBox[{{"\"Rotation\""}, {"\"matrix\""}}, GridBoxAlignment -> {"Columns" -> {{Left}}},
                 DefaultBaseStyle -> "Column", GridBoxItemSize -> {"Columns" -> {{All}}, "Rows" -> {{All}}}],
                "Column"], "", "", TagBox[RowBox[{"(", "⁠", GridBox[{{"1", "0", "0"}, {"0", "1", "0"},
                    {"0", "0", "1"}}, RowSpacings -> 1, ColumnSpacings -> 1, RowAlignments -> Baseline,
                   ColumnAlignments -> Center], "⁠", ")"}], Function[BoxForm`e$, MatrixForm[BoxForm`e$]]], TagBox[
                RowBox[{"(", "⁠", GridBox[{{"0", RowBox[{"-", "1"}], "0"}, {"1", RowBox[{"-", "1"}], "0"},
                    {"0", "0", "1"}}, RowSpacings -> 1, ColumnSpacings -> 1, RowAlignments -> Baseline,
                   ColumnAlignments -> Center], "⁠", ")"}], Function[BoxForm`e$, MatrixForm[BoxForm`e$]]], TagBox[
                RowBox[{"(", "⁠", GridBox[{{RowBox[{"-", "1"}], "1", "0"}, {RowBox[{"-", "1"}], "0", "0"},
                    {"0", "0", "1"}}, RowSpacings -> 1, ColumnSpacings -> 1, RowAlignments -> Baseline,
                   ColumnAlignments -> Center], "⁠", ")"}], Function[BoxForm`e$, MatrixForm[BoxForm`e$]]], TagBox[
                RowBox[{"(", "⁠", GridBox[{{RowBox[{"-", "1"}], "1", "0"}, {"0", "1", "0"}, {"0", "0", "1"}},
                   RowSpacings -> 1, ColumnSpacings -> 1, RowAlignments -> Baseline, ColumnAlignments -> Center],
                  "⁠", ")"}], Function[BoxForm`e$, MatrixForm[BoxForm`e$]]], TagBox[RowBox[{"(", "⁠",
                  GridBox[{{"1", "0", "0"}, {"1", RowBox[{"-", "1"}], "0"}, {"0", "0", "1"}}, RowSpacings -> 1,
                   ColumnSpacings -> 1, RowAlignments -> Baseline, ColumnAlignments -> Center], "⁠", ")"}],
                Function[BoxForm`e$, MatrixForm[BoxForm`e$]]], TagBox[RowBox[{"(", "⁠", GridBox[
                   {{"0", RowBox[{"-", "1"}], "0"}, {RowBox[{"-", "1"}], "0", "0"}, {"0", "0", "1"}}, RowSpacings -> 1,
                   ColumnSpacings -> 1, RowAlignments -> Baseline, ColumnAlignments -> Center], "⁠", ")"}],
                Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]}, {TagBox[GridBox[{{"\"Spin(↓↑)\""},
                  {"\"rotation\""}, {"\"matrix\""}}, GridBoxAlignment -> {"Columns" -> {{Left}}}, DefaultBaseStyle ->
                  "Column", GridBoxItemSize -> {"Columns" -> {{All}}, "Rows" -> {{All}}}], "Column"], "", "",
               TagBox[RowBox[{"(", "⁠", GridBox[{{"1", "0"}, {"0", "1"}}, RowSpacings -> 1, ColumnSpacings -> 1,
                   RowAlignments -> Baseline, ColumnAlignments -> Center], "⁠", ")"}], Function[BoxForm`e$,
                 MatrixForm[BoxForm`e$]]], TagBox[RowBox[{"(", "⁠", GridBox[{{RowBox[{FractionBox["1", "2"], "+",
                       FractionBox[RowBox[{"", " ", SqrtBox["3"]}], "2"]}], "0"}, {"0", RowBox[{FractionBox["1",
                        "2"], "-", FractionBox[RowBox[{"", " ", SqrtBox["3"]}], "2"]}]}}, RowSpacings -> 1,
                   ColumnSpacings -> 1, RowAlignments -> Baseline, ColumnAlignments -> Center], "⁠", ")"}],
                Function[BoxForm`e$, MatrixForm[BoxForm`e$]]], TagBox[RowBox[{"(", "⁠", GridBox[
                   {{RowBox[{FractionBox["1", "2"], "-", FractionBox[RowBox[{"", " ", SqrtBox["3"]}], "2"]}], "0"},
                    {"0", RowBox[{FractionBox["1", "2"], "+", FractionBox[RowBox[{"", " ", SqrtBox["3"]}], "2"]}]}},
                   RowSpacings -> 1, ColumnSpacings -> 1, RowAlignments -> Baseline, ColumnAlignments -> Center],
                  "⁠", ")"}], Function[BoxForm`e$, MatrixForm[BoxForm`e$]]], TagBox[RowBox[{"(", "⁠",
                  GridBox[{{"0", RowBox[{"-", "1"}]}, {"1", "0"}}, RowSpacings -> 1, ColumnSpacings -> 1,
                   RowAlignments -> Baseline, ColumnAlignments -> Center], "⁠", ")"}], Function[BoxForm`e$,
                 MatrixForm[BoxForm`e$]]], TagBox[RowBox[{"(", "⁠", GridBox[{{"0", RowBox[{RowBox[{"-", FractionBox[
                          "1", "2"]}], "+", FractionBox[RowBox[{"", " ", SqrtBox["3"]}], "2"]}]},
                    {RowBox[{FractionBox["1", "2"], "+", FractionBox[RowBox[{"", " ", SqrtBox["3"]}], "2"]}], "0"}},
                   RowSpacings -> 1, ColumnSpacings -> 1, RowAlignments -> Baseline, ColumnAlignments -> Center],
                  "⁠", ")"}], Function[BoxForm`e$, MatrixForm[BoxForm`e$]]], TagBox[RowBox[{"(", "⁠",
                  GridBox[{{"0", RowBox[{RowBox[{"-", FractionBox["1", "2"]}], "-", FractionBox[RowBox[{"", " ",
                          SqrtBox["3"]}], "2"]}]}, {RowBox[{FractionBox["1", "2"], "-", FractionBox[RowBox[
                         {"", " ", SqrtBox["3"]}], "2"]}], "0"}}, RowSpacings -> 1, ColumnSpacings -> 1,
                   RowAlignments -> Baseline, ColumnAlignments -> Center], "⁠", ")"}], Function[BoxForm`e$,
                 MatrixForm[BoxForm`e$]]]}, {"1", "\" \\!\\(\\*SubscriptBox[\\(K\\), \\(1\\)]\\) \"", "\"a\"", "1",
               "1", "1", "1", "1", "1"}, {"2", "\" \\!\\(\\*SubscriptBox[\\(K\\), \\(2\\)]\\) \"", "\"a\"", "1",
               SuperscriptBox["", RowBox[{"-", FractionBox[RowBox[{"2", " ", "", " ", "π"}], "3"]}]],
               SuperscriptBox["", FractionBox[RowBox[{"2", " ", "", " ", "π"}], "3"]], "1", SuperscriptBox[
                "", FractionBox[RowBox[{"2", " ", "", " ", "π"}], "3"]], SuperscriptBox["",
                RowBox[{"-", FractionBox[RowBox[{"2", " ", "", " ", "π"}], "3"]}]]},
              {"3", "\" \\!\\(\\*SubscriptBox[\\(K\\), \\(3\\)]\\) \"", "\"a\"", "1", SuperscriptBox["",
                FractionBox[RowBox[{"2", " ", "", " ", "π"}], "3"]], SuperscriptBox["",
                RowBox[{"-", FractionBox[RowBox[{"2", " ", "", " ", "π"}], "3"]}]], "1", SuperscriptBox["",
                RowBox[{"-", FractionBox[RowBox[{"2", " ", "", " ", "π"}], "3"]}]], SuperscriptBox["",
                FractionBox[RowBox[{"2", " ", "", " ", "π"}], "3"]]},
              {"4", "\" \\!\\(\\*SubscriptBox[\\(K\\), \\(4\\)]\\) \"", "\"a\"", "1", SuperscriptBox["",
                RowBox[{"-", FractionBox[RowBox[{"", " ", "π"}], "3"]}]], SuperscriptBox["",
                FractionBox[RowBox[{"", " ", "π"}], "3"]], "1", SuperscriptBox["", FractionBox[
                 RowBox[{"", " ", "π"}], "3"]], SuperscriptBox["", RowBox[{"-", FractionBox[
                   RowBox[{"", " ", "π"}], "3"]}]]}, {"5", "\" \\!\\(\\*SubscriptBox[\\(K\\), \\(5\\)]\\) \"",
               "\"a\"", "1", SuperscriptBox["", FractionBox[RowBox[{"", " ", "π"}], "3"]], SuperscriptBox["",
                RowBox[{"-", FractionBox[RowBox[{"", " ", "π"}], "3"]}]], "1", SuperscriptBox["",
                RowBox[{"-", FractionBox[RowBox[{"", " ", "π"}], "3"]}]], SuperscriptBox["",
                FractionBox[RowBox[{"", " ", "π"}], "3"]]},
              {"6", "\" \\!\\(\\*SubscriptBox[\\(K\\), \\(6\\)]\\) \"", "\"a\"", "1", RowBox[{"-", "1"}], RowBox[
                {"-", "1"}], "1", RowBox[{"-", "1"}], RowBox[{"-", "1"}]}}, GridBoxAlignment ->
              {"Columns" -> {{Center}}, "Rows" -> {{Center}}}, AutoDelete -> False, GridBoxBackground ->
              {"Columns" -> {{None}}, "Rows" -> {{None}}, "ItemsIndexed" -> {{{1, 4}, {1, -1}} -> RGBColor[1., 1.,
                   0.9], {{2, 2}, {1, -1}} -> RGBColor[1., 1., 0.9], {}, {{5, 7}, {1, -1}} -> RGBColor[0.95, 1., 0.95],
                 {}, {{8, -1}, {1, -1}} -> RGBColor[0.95, 0.95, 1.]}}, GridBoxDividers -> {"Columns" ->
                {{Directive[Thickness[Tiny], GrayLevel[0.8]]}}, "ColumnsIndexed" -> {1 -> Directive[GrayLevel[0],
                   Thickness[2]], 4 -> Directive[GrayLevel[0], Thickness[2]], -1 -> Directive[GrayLevel[0],
                   Thickness[2]], 2 -> Directive[Thickness[Tiny], GrayLevel[0.8]], 3 -> Directive[Thickness[Tiny],
                   GrayLevel[0.8]]}, "Rows" -> {{Directive[Thickness[Tiny], GrayLevel[0.8]]}}, "RowsIndexed" ->
                {1 -> Directive[GrayLevel[0], Thickness[2]], 5 -> Directive[GrayLevel[0], Thickness[2]],
                 8 -> Directive[GrayLevel[0], Thickness[2]], -1 -> Directive[GrayLevel[0], Thickness[2]],
                 2 -> Directive[Thickness[Tiny], GrayLevel[0.8]]}}, GridBoxFrame -> {"Columns" -> {{True}}, "Rows" ->
                {{True}}}, GridBoxItemSize -> {"Columns" -> {{All}}, "Rows" -> {{All}}}], "Grid"]}},
         GridBoxAlignment -> {"Columns" -> {{Left}}}, DefaultBaseStyle -> "Column",
         GridBoxItemSize -> {"Columns" -> {{Automatic}}, "Rows" -> {{Automatic}}}], "Column"]], "Output"],
     Cell[
      "在分数坐标下构造二阶模型，把 K 点附近的局域动量平移到 EIGENVAL 中的绝对坐标，并为所有生成参数准备初值："\
, "ExampleText"], Cell[BoxData["input = interfaceRep[{156, 50}, \"K\", {2, 3}, \"CartesianCoordinates\" -> False];\nham \
= kpHam[2, input][\"ham\"];\nhamshift = ham /. {kx -> kx - 1/3, ky -> ky - 1/3};\nparams = \
Complement[Variables[hamshift], {kx, ky, kz}];\ninitialRules = Thread[params -> {0., 0.19, 4.33, 4.33, 0., 0., 0., 0., \
0.}];\nMatrixForm[hamshift]"], "Input"],
     Cell[BoxData[TagBox[RowBox[{"(", "⁠", GridBox[{{RowBox[{SubscriptBox["e", "1"], "+", SubscriptBox["e", "2"],
               "+", RowBox[{RowBox[{"(", RowBox[{SuperscriptBox[RowBox[{"(", RowBox[{RowBox[{"-", FractionBox["1",
                             "3"]}], "+", "kx"}], ")"}], "2"], "+", RowBox[{RowBox[{"(", RowBox[{RowBox[{"-",
                             FractionBox["1", "3"]}], "+", "kx"}], ")"}], " ", RowBox[{"(", RowBox[
                          {RowBox[{"-", FractionBox["1", "3"]}], "+", "ky"}], ")"}]}], "+", SuperscriptBox[
                      RowBox[{"(", RowBox[{RowBox[{"-", FractionBox["1", "3"]}], "+", "ky"}], ")"}], "2"]}], ")"}],
                 " ", SubscriptBox["r", "1"]}], "+", RowBox[{RowBox[{"(", RowBox[{SuperscriptBox[RowBox[
                       {"(", RowBox[{RowBox[{"-", FractionBox["1", "3"]}], "+", "kx"}], ")"}], "2"], "+",
                     RowBox[{RowBox[{"(", RowBox[{RowBox[{"-", FractionBox["1", "3"]}], "+", "kx"}], ")"}], " ",
                       RowBox[{"(", RowBox[{RowBox[{"-", FractionBox["1", "3"]}], "+", "ky"}], ")"}]}], "+",
                     SuperscriptBox[RowBox[{"(", RowBox[{RowBox[{"-", FractionBox["1", "3"]}], "+", "ky"}], ")"}],
                      "2"]}], ")"}], " ", SubscriptBox["r", "3"]}], "+", RowBox[{SuperscriptBox["kz", "2"], " ",
                 SubscriptBox["r", "5"]}], "+", RowBox[{SuperscriptBox["kz", "2"], " ", SubscriptBox["r", "6"]}]}],
             RowBox[{RowBox[{RowBox[{"(", RowBox[{RowBox[{RowBox[{"-", FractionBox["1", "2"]}], " ", SuperscriptBox[
                        RowBox[{"(", RowBox[{RowBox[{"-", FractionBox["1", "3"]}], "+", "kx"}], ")"}], "2"]}], "-",
                     RowBox[{FractionBox["1", "2"], " ", "", " ", SqrtBox["3"], " ", SuperscriptBox[RowBox[
                         {"(", RowBox[{RowBox[{"-", FractionBox["1", "3"]}], "+", "kx"}], ")"}], "2"]}], "+",
                     RowBox[{RowBox[{"(", RowBox[{RowBox[{"-", FractionBox["1", "3"]}], "+", "kx"}], ")"}], " ",
                       RowBox[{"(", RowBox[{RowBox[{"-", FractionBox["1", "3"]}], "+", "ky"}], ")"}]}], "-",
                     RowBox[{"", " ", SqrtBox["3"], " ", RowBox[{"(", RowBox[{RowBox[{"-", FractionBox["1", "3"]}],
                           "+", "kx"}], ")"}], " ", RowBox[{"(", RowBox[{RowBox[{"-", FractionBox["1", "3"]}], "+",
                           "ky"}], ")"}]}], "+", SuperscriptBox[RowBox[{"(", RowBox[{RowBox[{"-", FractionBox["1",
                             "3"]}], "+", "ky"}], ")"}], "2"]}], ")"}], " ", SubscriptBox["r", "2"]}], "+", RowBox[
                {RowBox[{"(", RowBox[{RowBox[{FractionBox["1", "2"], " ", RowBox[{"(", RowBox[{RowBox[{"-", FractionBox[
                              "1", "3"]}], "+", "kx"}], ")"}], " ", "kz"}], "+", RowBox[{FractionBox["1", "2"], " ",
                       "", " ", SqrtBox["3"], " ", RowBox[{"(", RowBox[{RowBox[{"-", FractionBox["1", "3"]}], "+",
                           "kx"}], ")"}], " ", "kz"}], "+", RowBox[{RowBox[{"(", RowBox[{RowBox[{"-", FractionBox["1",
                              "3"]}], "+", "ky"}], ")"}], " ", "kz"}]}], ")"}], " ", SubscriptBox["r", "4"]}], "+",
               RowBox[{RowBox[{"(", RowBox[{RowBox[{RowBox[{"-", FractionBox["1", "2"]}], " ", "", " ",
                       RowBox[{"(", RowBox[{RowBox[{"-", FractionBox["1", "3"]}], "+", "kx"}], ")"}]}], "+",
                     RowBox[{FractionBox["1", "2"], " ", SqrtBox["3"], " ", RowBox[{"(", RowBox[{RowBox[{"-",
                             FractionBox["1", "3"]}], "+", "kx"}], ")"}]}], "-", RowBox[{"", " ", RowBox[
                        {"(", RowBox[{RowBox[{"-", FractionBox["1", "3"]}], "+", "ky"}], ")"}]}]}], ")"}], " ",
                 SubscriptBox["t", "1"]}]}]}, {RowBox[{RowBox[{RowBox[{"(", RowBox[{RowBox[{RowBox[{"-", FractionBox[
                          "1", "2"]}], " ", SuperscriptBox[RowBox[{"(", RowBox[{RowBox[{"-", FractionBox["1", "3"]}],
                            "+", "kx"}], ")"}], "2"]}], "+", RowBox[{FractionBox["1", "2"], " ", "", " ",
                       SqrtBox["3"], " ", SuperscriptBox[RowBox[{"(", RowBox[{RowBox[{"-", FractionBox["1", "3"]}],
                            "+", "kx"}], ")"}], "2"]}], "+", RowBox[{RowBox[{"(", RowBox[{RowBox[{"-", FractionBox["1",
                              "3"]}], "+", "kx"}], ")"}], " ", RowBox[{"(", RowBox[{RowBox[{"-", FractionBox["1",
                              "3"]}], "+", "ky"}], ")"}]}], "+", RowBox[{"", " ", SqrtBox["3"], " ",
                       RowBox[{"(", RowBox[{RowBox[{"-", FractionBox["1", "3"]}], "+", "kx"}], ")"}], " ",
                       RowBox[{"(", RowBox[{RowBox[{"-", FractionBox["1", "3"]}], "+", "ky"}], ")"}]}], "+",
                     SuperscriptBox[RowBox[{"(", RowBox[{RowBox[{"-", FractionBox["1", "3"]}], "+", "ky"}], ")"}],
                      "2"]}], ")"}], " ", SubscriptBox["r", "2"]}], "+", RowBox[
                {RowBox[{"(", RowBox[{RowBox[{FractionBox["1", "2"], " ", RowBox[{"(", RowBox[{RowBox[{"-", FractionBox[
                              "1", "3"]}], "+", "kx"}], ")"}], " ", "kz"}], "-", RowBox[{FractionBox["1", "2"], " ",
                       "", " ", SqrtBox["3"], " ", RowBox[{"(", RowBox[{RowBox[{"-", FractionBox["1", "3"]}], "+",
                           "kx"}], ")"}], " ", "kz"}], "+", RowBox[{RowBox[{"(", RowBox[{RowBox[{"-", FractionBox["1",
                              "3"]}], "+", "ky"}], ")"}], " ", "kz"}]}], ")"}], " ", SubscriptBox["r", "4"]}], "+",
               RowBox[{RowBox[{"(", RowBox[{RowBox[{FractionBox["1", "2"], " ", "", " ", RowBox[{"(",
                         RowBox[{RowBox[{"-", FractionBox["1", "3"]}], "+", "kx"}], ")"}]}], "+", RowBox[
                      {FractionBox["1", "2"], " ", SqrtBox["3"], " ", RowBox[{"(", RowBox[{RowBox[{"-", FractionBox[
                              "1", "3"]}], "+", "kx"}], ")"}]}], "+", RowBox[{"", " ", RowBox[{"(", RowBox[
                          {RowBox[{"-", FractionBox["1", "3"]}], "+", "ky"}], ")"}]}]}], ")"}], " ",
                 SubscriptBox["t", "1"]}]}], RowBox[{SubscriptBox["e", "1"], "-", SubscriptBox["e", "2"], "+", RowBox[
                {RowBox[{"(", RowBox[{SuperscriptBox[RowBox[{"(", RowBox[{RowBox[{"-", FractionBox["1", "3"]}], "+",
                          "kx"}], ")"}], "2"], "+", RowBox[{RowBox[{"(", RowBox[{RowBox[{"-", FractionBox["1", "3"]}],
                           "+", "kx"}], ")"}], " ", RowBox[{"(", RowBox[{RowBox[{"-", FractionBox["1", "3"]}], "+",
                           "ky"}], ")"}]}], "+", SuperscriptBox[RowBox[{"(", RowBox[{RowBox[{"-", FractionBox["1",
                             "3"]}], "+", "ky"}], ")"}], "2"]}], ")"}], " ", SubscriptBox["r", "1"]}], "+", RowBox[
                {RowBox[{"(", RowBox[{RowBox[{"-", SuperscriptBox[RowBox[{"(", RowBox[{RowBox[{"-", FractionBox["1",
                               "3"]}], "+", "kx"}], ")"}], "2"]}], "-", RowBox[{RowBox[{"(", RowBox[
                          {RowBox[{"-", FractionBox["1", "3"]}], "+", "kx"}], ")"}], " ", RowBox[{"(", RowBox[
                          {RowBox[{"-", FractionBox["1", "3"]}], "+", "ky"}], ")"}]}], "-", SuperscriptBox[
                      RowBox[{"(", RowBox[{RowBox[{"-", FractionBox["1", "3"]}], "+", "ky"}], ")"}], "2"]}], ")"}],
                 " ", SubscriptBox["r", "3"]}], "+", RowBox[{SuperscriptBox["kz", "2"], " ", SubscriptBox["r", "5"]}],
               "-", RowBox[{SuperscriptBox["kz", "2"], " ", SubscriptBox["r", "6"]}]}]}}, RowSpacings -> 1,
           ColumnSpacings -> 1, RowAlignments -> Baseline, ColumnAlignments -> Center], "⁠", ")"}],
        Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]], "Output"],
     Cell["从随附的原始 EIGENVAL 文件读取 VASP 第 21、22 条能带：", "ExampleText"],
     Cell[BoxData["file = FileNameJoin[{NotebookDirectory[], \"EIGENVAL\"}];\neigdata = vaspEig[file, -0.0072, 1, 21, \
22];\n{Length[eigdata], First[eigdata]}"], "Input"],
     Cell[BoxData[RowBox[{"{", RowBox[{"90", ",", RowBox[{"{", RowBox[{RowBox[{"{", RowBox[{"0.`", ",", "0.`", ",",
                   "0.`"}], "}"}], ",", RowBox[{"{", RowBox[{RowBox[{"-", "0.337773`"}], ",", "2.612089`"}], "}"}]}],
             "}"}]}], "}"}]], "Output"],
     Cell["拟合 K 点附近的五个点，并查看拟合规则、误差下降和收敛状态：", "ExampleText"],
     Cell[BoxData["fit = fittingKP[\n  hamshift, eigdata, Range[50, 70, 5], initialRules\n];\nKeyTake[\n  fit,\n  \
{\"FittedParams\", \"LSQForInitParams\", \"LSQForFittedParams\", \"Converged\", \"Iterations\"}\n]"], "Input"],
     Cell[BoxData[RowBox[{"", RowBox[{RowBox[{"\"FittedParams\"", "",
             RowBox[{"{", RowBox[{RowBox[{SubscriptBox["e", "1"], "", "0.002835954779085878`"}], ",",
                 RowBox[{SubscriptBox["e", "2"], "", "0.13242722920568564`"}], ",", RowBox[{SubscriptBox["r", "1"],
                   "", "2.373071416840521`"}], ",", RowBox[{SubscriptBox["r", "2"], "", "2.499386631580799`"}],
                 ",", RowBox[{SubscriptBox["r", "3"], "", RowBox[{"-", "0.03396194828064344`"}]}], ",",
                 RowBox[{SubscriptBox["r", "4"], "", "0.`"}], ",", RowBox[{SubscriptBox["r", "5"], "", "0.`"}],
                 ",", RowBox[{SubscriptBox["r", "6"], "", "0.`"}], ",", RowBox[{SubscriptBox["t", "1"], "",
                   RowBox[{"-", "5.574114096277786`"}]}]}], "}"}]}], ",", RowBox[{"\"LSQForInitParams\"", "",
             "1.8981578637374088`"}], ",", RowBox[{"\"LSQForFittedParams\"", "", "0.0016921734401136772`"}], ",",
           RowBox[{"\"Converged\"", "", "True"}], ",", RowBox[{"\"Iterations\"", "", "13"}]}], ""}]],
      "Output"], Cell["比较选定的 VASP 数据点、初始 k.p 能带和拟合后的 k.p 能带：",
      "ExampleText"], Cell[BoxData["fit[\"ComparisonPlot\"]"], "Input"],
     Cell[BoxData[TemplateBox[{GraphicsBox[{{{}, InterpretationBox[
             {TagBox[{{TagBox[{Directive[PointSize[0.012833333333333334], AbsoluteThickness[2], RGBColor[0, 0, 1],
                    AbsolutePointSize[5]], GeometricTransformationBox[InsetBox[BoxData[FormBox[StyleBox[GraphicsBox[
                         {EdgeForm[], DiskBox[{0, 0}]}, PlotRangePadding -> Scaled[0.15]], GraphicsBoxOptions ->
                         {DefaultBaseStyle -> Directive[PointSize[0.012833333333333334], AbsoluteThickness[2],
                            RGBColor[0, 0, 1], AbsolutePointSize[5]]}, StripOnInput -> False], TraditionalForm]],
                     {0., 0.}, Automatic, Offset[10.976]], {{{50., -0.506383}}, {{55., -0.296921}},
                     {{60., -0.12791799999999998}}, {{65., -0.484715}}, {{70., -1.0025769999999998}}}]},
                  Annotation[#1, "Charting`Private`Tag#1"] & ], TagBox[{Directive[PointSize[0.012833333333333334],
                    AbsoluteThickness[2], RGBColor[0, 0, 1], AbsolutePointSize[5]], GeometricTransformationBox[
                    InsetBox[BoxData[FormBox[StyleBox[GraphicsBox[{EdgeForm[], PolygonBox[{{-1, -1}, {1, -1}, {1, 1},
                            {-1, 1}}]}, PlotRangePadding -> Scaled[0.15]], GraphicsBoxOptions -> {DefaultBaseStyle ->
                           Directive[PointSize[0.012833333333333334], AbsoluteThickness[2], RGBColor[0, 0, 1],
                            AbsolutePointSize[5]]}, StripOnInput -> False], TraditionalForm]], {0., 0.}, Automatic,
                     Offset[10.5]], {{{50., 0.569574}}, {{55., 0.31285399999999997}}, {{60., 0.126005}},
                     {{65., 0.52262}}, {{70., 1.157078}}}]}, Annotation[#1, "Charting`Private`Tag#2"] & ]}}, {
                "WolframDynamicHighlight", Association["Label" -> {"XYLabel"}, "Ball" -> {"IndicatedBall"}]}],
              StyleBox[DynamicBox[(Charting`HighlightActionBox["DynamicHighlight", {}, #HighlightElements,
                   #LayoutOptions, #Meta, Charting`HighlightActionFunction["DynamicHighlight",
                    {{Annotation[{Directive[PointSize[0.012833333333333334], AbsoluteThickness[2], RGBColor[0, 0, 1],
                         AbsolutePointSize[5]], GeometricTransformation[Inset[Style[Graphics[{EdgeForm[],
                             Disk[{0, 0}]}, PlotRangePadding -> Scaled[0.15]], GraphicsBoxOptions ->
                            {DefaultBaseStyle -> Directive[PointSize[0.012833333333333334], AbsoluteThickness[2],
                               RGBColor[0, 0, 1], AbsolutePointSize[5]]}], {0., 0.}, Automatic, Offset[10.976]],
                         {{{50., -0.506383}}, {{55., -0.296921}}, {{60., -0.12791799999999998}}, {{65., -0.484715}},
                          {{70., -1.0025769999999998}}}]}, "Charting`Private`Tag#1"], Annotation[
                       {Directive[PointSize[0.012833333333333334], AbsoluteThickness[2], RGBColor[0, 0, 1],
                         AbsolutePointSize[5]], GeometricTransformation[Inset[Style[Graphics[{EdgeForm[], Polygon[
                              {{-1, -1}, {1, -1}, {1, 1}, {-1, 1}}]}, PlotRangePadding -> Scaled[0.15]],
                           GraphicsBoxOptions -> {DefaultBaseStyle -> Directive[PointSize[0.012833333333333334],
                               AbsoluteThickness[2], RGBColor[0, 0, 1], AbsolutePointSize[5]]}], {0., 0.}, Automatic,
                          Offset[10.5]], {{{50., 0.569574}}, {{55., 0.31285399999999997}}, {{60., 0.126005}},
                          {{65., 0.52262}}, {{70., 1.157078}}}]}, "Charting`Private`Tag#2"]}}, Association[
                     "HighlightElements" -> Association["Label" -> {"XYLabel"}, "Ball" -> {"IndicatedBall"}],
                     "LayoutOptions" -> Association["PanelPlotLayout" -> Association[], "PlotRange" ->
                        {{49.58333333333344, 70.}, {-1.0025769999999998, 1.157078}}, "Frame" -> {{False, False},
                         {False, False}}, "AxesOrigin" -> {49.58333333333344, 0}, "ImageSize" -> {360,
                         360/GoldenRatio}, "Axes" -> {True, True}, "LabelStyle" -> {}, "AspectRatio" -> GoldenRatio^
                         (-1), "DefaultStyle" -> {Directive[PointSize[0.012833333333333334], AbsoluteThickness[2],
                          RGBColor[0, 0, 1], AbsolutePointSize[5]], Directive[PointSize[0.012833333333333334],
                          AbsoluteThickness[2], RGBColor[0, 0, 1], AbsolutePointSize[5]]},
                       "HighlightLabelingFunctions" -> Association["CoordinatesToolOptions" -> Identity,
                         "ScalingFunctions" -> {{Identity, Identity}, {Identity, Identity}}], "Primitives" -> {},
                       "GCFlag" -> False], "Meta" -> Association["DefaultHighlight" -> {"Dynamic", None},
                       "Index" -> {}, "Function" -> ListPlot, "GroupHighlight" -> False]]]] & )[Association[
                  "HighlightElements" -> Association["Label" -> {"XYLabel"}, "Ball" -> {"IndicatedBall"}],
                  "LayoutOptions" -> Association["PanelPlotLayout" -> Association[], "PlotRange" ->
                     {{49.58333333333344, 70.}, {-1.0025769999999998, 1.157078}}, "Frame" -> {{False, False},
                      {False, False}}, "AxesOrigin" -> {49.58333333333344, 0}, "ImageSize" -> {360, 360/GoldenRatio},
                    "Axes" -> {True, True}, "LabelStyle" -> {}, "AspectRatio" -> GoldenRatio^(-1), "DefaultStyle" ->
                     {Directive[PointSize[0.012833333333333334], AbsoluteThickness[2], RGBColor[0, 0, 1],
                       AbsolutePointSize[5]], Directive[PointSize[0.012833333333333334], AbsoluteThickness[2],
                       RGBColor[0, 0, 1], AbsolutePointSize[5]]}, "HighlightLabelingFunctions" -> Association[
                      "CoordinatesToolOptions" -> Identity, "ScalingFunctions" -> {{Identity, Identity}, {Identity,
                         Identity}}], "Primitives" -> {}, "GCFlag" -> False], "Meta" -> Association[
                    "DefaultHighlight" -> {"Dynamic", None}, "Index" -> {}, "Function" -> ListPlot, "GroupHighlight" ->
                     False]]]], Selectable -> False]}, Annotation[{{Annotation[{Directive[PointSize[
                    0.012833333333333334], AbsoluteThickness[2], RGBColor[0, 0, 1], AbsolutePointSize[5]],
                  GeometricTransformation[Inset[Style[Graphics[{EdgeForm[], Disk[{0, 0}]}, PlotRangePadding ->
                       Scaled[0.15]], GraphicsBoxOptions -> {DefaultBaseStyle -> Directive[PointSize[
                          0.012833333333333334], AbsoluteThickness[2], RGBColor[0, 0, 1], AbsolutePointSize[5]]}],
                    {0., 0.}, Automatic, Offset[10.976]], {{{50., -0.506383}}, {{55., -0.296921}},
                    {{60., -0.12791799999999998}}, {{65., -0.484715}}, {{70., -1.0025769999999998}}}]},
                 "Charting`Private`Tag#1"], Annotation[{Directive[PointSize[0.012833333333333334], AbsoluteThickness[
                    2], RGBColor[0, 0, 1], AbsolutePointSize[5]], GeometricTransformation[Inset[
                    Style[Graphics[{EdgeForm[], Polygon[{{-1, -1}, {1, -1}, {1, 1}, {-1, 1}}]}, PlotRangePadding ->
                       Scaled[0.15]], GraphicsBoxOptions -> {DefaultBaseStyle -> Directive[PointSize[
                          0.012833333333333334], AbsoluteThickness[2], RGBColor[0, 0, 1], AbsolutePointSize[5]]}],
                    {0., 0.}, Automatic, Offset[10.5]], {{{50., 0.569574}}, {{55., 0.31285399999999997}},
                    {{60., 0.126005}}, {{65., 0.52262}}, {{70., 1.157078}}}]}, "Charting`Private`Tag#2"]}},
              Association["HighlightElements" -> Association["Label" -> {"XYLabel"}, "Ball" -> {"IndicatedBall"}],
               "LayoutOptions" -> Association["PanelPlotLayout" -> Association[], "PlotRange" ->
                  {{49.58333333333344, 70.}, {-1.0025769999999998, 1.157078}}, "Frame" -> {{False, False},
                   {False, False}}, "AxesOrigin" -> {49.58333333333344, 0}, "ImageSize" -> {360, 360/GoldenRatio},
                 "Axes" -> {True, True}, "LabelStyle" -> {}, "AspectRatio" -> GoldenRatio^(-1), "DefaultStyle" ->
                  {Directive[PointSize[0.012833333333333334], AbsoluteThickness[2], RGBColor[0, 0, 1],
                    AbsolutePointSize[5]], Directive[PointSize[0.012833333333333334], AbsoluteThickness[2],
                    RGBColor[0, 0, 1], AbsolutePointSize[5]]}, "HighlightLabelingFunctions" -> Association[
                   "CoordinatesToolOptions" -> Identity, "ScalingFunctions" -> {{Identity, Identity}, {Identity,
                      Identity}}], "Primitives" -> {}, "GCFlag" -> False], "Meta" -> Association["DefaultHighlight" ->
                  {"Dynamic", None}, "Index" -> {}, "Function" -> ListPlot, "GroupHighlight" -> False]],
              "DynamicHighlight"]], {{}, {}}}, {{}, InterpretationBox[{TagBox[{{{}, {}, TagBox[{Hue[0.67, 0.6, 0.6],
                   Directive[PointSize[NCache[7/360, 0.019444444444444445]], AbsoluteThickness[2], GrayLevel[0.5],
                    Dashing[{Small, Small}], Thickness[Tiny]], LineBox[{{50., -0.15187886871475645},
                     {55., -0.1795762256065017}, {60., -0.1899999999999857}, {65., -0.164514605267474},
                     {70., -0.09641134208667265}}]}, Annotation[#1, "Charting`Private`Tag#1"] & ],
                 TagBox[{Hue[0.9060679774997897, 0.6, 0.6], Directive[PointSize[NCache[7/360, 0.019444444444444445]],
                    AbsoluteThickness[2], GrayLevel[0.5], Dashing[{Small, Small}], Thickness[Tiny]],
                   LineBox[{{50., 0.23768941858395973}, {55., 0.2010288381887729}, {60., 0.1900000000000143},
                     {65., 0.21943340496308672}, {70., 0.3744372728215587}}]}, Annotation[#1,
                    "Charting`Private`Tag#2"] & ]}}, {"WolframDynamicHighlight", Association["Label" -> {"XYLabel"},
                 "Ball" -> {"InterpolatedBall"}]}], StyleBox[DynamicBox[(Charting`HighlightActionBox[
                   "DynamicHighlight", {}, #HighlightElements, #LayoutOptions, #Meta, Charting`HighlightActionFunction[
                    "DynamicHighlight", {{{}, {}, Annotation[{Hue[0.67, 0.6, 0.6], Directive[PointSize[7/360],
                         AbsoluteThickness[2], GrayLevel[0.5], Dashing[{Small, Small}], Thickness[Tiny]],
                        Line[{{50., -0.15187886871475645}, {55., -0.1795762256065017}, {60., -0.1899999999999857},
                          {65., -0.164514605267474}, {70., -0.09641134208667265}}]}, "Charting`Private`Tag#1"],
                      Annotation[{Hue[0.9060679774997897, 0.6, 0.6], Directive[PointSize[7/360], AbsoluteThickness[2],
                         GrayLevel[0.5], Dashing[{Small, Small}], Thickness[Tiny]], Line[{{50., 0.23768941858395973},
                          {55., 0.2010288381887729}, {60., 0.1900000000000143}, {65., 0.21943340496308672},
                          {70., 0.3744372728215587}}]}, "Charting`Private`Tag#2"]}}, Association["HighlightElements" ->
                      Association["Label" -> {"XYLabel"}, "Ball" -> {"InterpolatedBall"}], "LayoutOptions" ->
                      Association["PanelPlotLayout" -> Association[], "PlotRange" -> {{49.58333333333344, 70.},
                         {-0.1899999999999857, 0.3744372728215587}}, "Frame" -> {{False, False}, {False, False}},
                       "AxesOrigin" -> {49.58333333333344, 0}, "ImageSize" -> {360, 360/GoldenRatio}, "Axes" ->
                        {True, True}, "LabelStyle" -> {}, "AspectRatio" -> GoldenRatio^(-1), "DefaultStyle" ->
                        {Directive[PointSize[7/360], AbsoluteThickness[2], GrayLevel[0.5], Dashing[{Small, Small}],
                          Thickness[Tiny]], Directive[PointSize[7/360], AbsoluteThickness[2], GrayLevel[0.5],
                          Dashing[{Small, Small}], Thickness[Tiny]]}, "HighlightLabelingFunctions" -> Association[
                         "CoordinatesToolOptions" -> Identity, "ScalingFunctions" -> {{Identity, Identity},
                           {Identity, Identity}}], "Primitives" -> {}, "GCFlag" -> False], "Meta" -> Association[
                       "DefaultHighlight" -> {"Dynamic", None}, "Index" -> {}, "Function" -> ListLinePlot,
                       "GroupHighlight" -> False]]]] & )[Association["HighlightElements" -> Association[
                    "Label" -> {"XYLabel"}, "Ball" -> {"InterpolatedBall"}], "LayoutOptions" -> Association[
                    "PanelPlotLayout" -> Association[], "PlotRange" -> {{49.58333333333344, 70.}, {-0.1899999999999857,
                       0.3744372728215587}}, "Frame" -> {{False, False}, {False, False}}, "AxesOrigin" ->
                     {49.58333333333344, 0}, "ImageSize" -> {360, 360/GoldenRatio}, "Axes" -> {True, True},
                    "LabelStyle" -> {}, "AspectRatio" -> GoldenRatio^(-1), "DefaultStyle" -> {Directive[PointSize[
                        7/360], AbsoluteThickness[2], GrayLevel[0.5], Dashing[{Small, Small}], Thickness[Tiny]],
                      Directive[PointSize[7/360], AbsoluteThickness[2], GrayLevel[0.5], Dashing[{Small, Small}],
                       Thickness[Tiny]]}, "HighlightLabelingFunctions" -> Association["CoordinatesToolOptions" ->
                       Identity, "ScalingFunctions" -> {{Identity, Identity}, {Identity, Identity}}], "Primitives" ->
                     {}, "GCFlag" -> False], "Meta" -> Association["DefaultHighlight" -> {"Dynamic", None},
                    "Index" -> {}, "Function" -> ListLinePlot, "GroupHighlight" -> False]]]], Selectable -> False]},
             Annotation[{{{}, {}, Annotation[{Hue[0.67, 0.6, 0.6], Directive[PointSize[7/360], AbsoluteThickness[2],
                   GrayLevel[0.5], Dashing[{Small, Small}], Thickness[Tiny]], Line[{{50., -0.15187886871475645},
                    {55., -0.1795762256065017}, {60., -0.1899999999999857}, {65., -0.164514605267474},
                    {70., -0.09641134208667265}}]}, "Charting`Private`Tag#1"], Annotation[{Hue[0.9060679774997897, 0.6,
                   0.6], Directive[PointSize[7/360], AbsoluteThickness[2], GrayLevel[0.5], Dashing[{Small, Small}],
                   Thickness[Tiny]], Line[{{50., 0.23768941858395973}, {55., 0.2010288381887729},
                    {60., 0.1900000000000143}, {65., 0.21943340496308672}, {70., 0.3744372728215587}}]},
                 "Charting`Private`Tag#2"]}}, Association["HighlightElements" -> Association["Label" -> {"XYLabel"},
                 "Ball" -> {"InterpolatedBall"}], "LayoutOptions" -> Association["PanelPlotLayout" -> Association[],
                 "PlotRange" -> {{49.58333333333344, 70.}, {-0.1899999999999857, 0.3744372728215587}},
                 "Frame" -> {{False, False}, {False, False}}, "AxesOrigin" -> {49.58333333333344, 0},
                 "ImageSize" -> {360, 360/GoldenRatio}, "Axes" -> {True, True}, "LabelStyle" -> {}, "AspectRatio" ->
                  GoldenRatio^(-1), "DefaultStyle" -> {Directive[PointSize[7/360], AbsoluteThickness[2],
                    GrayLevel[0.5], Dashing[{Small, Small}], Thickness[Tiny]], Directive[PointSize[7/360],
                    AbsoluteThickness[2], GrayLevel[0.5], Dashing[{Small, Small}], Thickness[Tiny]]},
                 "HighlightLabelingFunctions" -> Association["CoordinatesToolOptions" -> Identity,
                   "ScalingFunctions" -> {{Identity, Identity}, {Identity, Identity}}], "Primitives" -> {},
                 "GCFlag" -> False], "Meta" -> Association["DefaultHighlight" -> {"Dynamic", None}, "Index" -> {},
                 "Function" -> ListLinePlot, "GroupHighlight" -> False]], "DynamicHighlight"]], {{}, {}}},
           {{}, InterpretationBox[{TagBox[{{{}, {}, TagBox[{Hue[0.67, 0.6, 0.6], Directive[PointSize[NCache[7/360,
                      0.019444444444444445]], AbsoluteThickness[2], RGBColor[1, 0, 0], Thickness[Large]],
                   LineBox[{{50., -0.5199579745948327}, {55., -0.2930904620581086}, {60., -0.1295912744269829},
                     {65., -0.46049023179894427}, {70., -1.0079395588481899}}]}, Annotation[#1,
                    "Charting`Private`Tag#1"] & ], TagBox[{Hue[0.9060679774997897, 0.6, 0.6], Directive[
                    PointSize[NCache[7/360, 0.019444444444444445]], AbsoluteThickness[2], RGBColor[1, 0, 0],
                    Thickness[Large]], LineBox[{{50., 0.5726586516259777}, {55., 0.31051954984619734},
                     {60., 0.1352631839851702}, {65., 0.4962605787254412}, {70., 1.165984537544339}}]},
                  Annotation[#1, "Charting`Private`Tag#2"] & ]}}, {"WolframDynamicHighlight", Association[
                 "Label" -> {"XYLabel"}, "Ball" -> {"InterpolatedBall"}]}], StyleBox[DynamicBox[
                (Charting`HighlightActionBox["DynamicHighlight", {}, #HighlightElements, #LayoutOptions, #Meta,
                   Charting`HighlightActionFunction["DynamicHighlight", {{{}, {}, Annotation[{Hue[0.67, 0.6, 0.6],
                        Directive[PointSize[7/360], AbsoluteThickness[2], RGBColor[1, 0, 0], Thickness[Large]],
                        Line[{{50., -0.5199579745948327}, {55., -0.2930904620581086}, {60., -0.1295912744269829},
                          {65., -0.46049023179894427}, {70., -1.0079395588481899}}]}, "Charting`Private`Tag#1"],
                      Annotation[{Hue[0.9060679774997897, 0.6, 0.6], Directive[PointSize[7/360], AbsoluteThickness[2],
                         RGBColor[1, 0, 0], Thickness[Large]], Line[{{50., 0.5726586516259777}, {55.,
                           0.31051954984619734}, {60., 0.1352631839851702}, {65., 0.4962605787254412}, {70.,
                           1.165984537544339}}]}, "Charting`Private`Tag#2"]}}, Association["HighlightElements" ->
                      Association["Label" -> {"XYLabel"}, "Ball" -> {"InterpolatedBall"}], "LayoutOptions" ->
                      Association["PanelPlotLayout" -> Association[], "PlotRange" -> {{49.58333333333344, 70.},
                         {-1.0079395588481899, 1.165984537544339}}, "Frame" -> {{False, False}, {False, False}},
                       "AxesOrigin" -> {49.58333333333344, 0}, "ImageSize" -> {360, 360/GoldenRatio}, "Axes" ->
                        {True, True}, "LabelStyle" -> {}, "AspectRatio" -> GoldenRatio^(-1), "DefaultStyle" ->
                        {Directive[PointSize[7/360], AbsoluteThickness[2], RGBColor[1, 0, 0], Thickness[Large]],
                         Directive[PointSize[7/360], AbsoluteThickness[2], RGBColor[1, 0, 0], Thickness[Large]]},
                       "HighlightLabelingFunctions" -> Association["CoordinatesToolOptions" -> Identity,
                         "ScalingFunctions" -> {{Identity, Identity}, {Identity, Identity}}], "Primitives" -> {},
                       "GCFlag" -> False], "Meta" -> Association["DefaultHighlight" -> {"Dynamic", None},
                       "Index" -> {}, "Function" -> ListLinePlot, "GroupHighlight" -> False]]]] & )[
                 Association["HighlightElements" -> Association["Label" -> {"XYLabel"}, "Ball" ->
                     {"InterpolatedBall"}], "LayoutOptions" -> Association["PanelPlotLayout" -> Association[],
                    "PlotRange" -> {{49.58333333333344, 70.}, {-1.0079395588481899, 1.165984537544339}},
                    "Frame" -> {{False, False}, {False, False}}, "AxesOrigin" -> {49.58333333333344, 0},
                    "ImageSize" -> {360, 360/GoldenRatio}, "Axes" -> {True, True}, "LabelStyle" -> {}, "AspectRatio" ->
                     GoldenRatio^(-1), "DefaultStyle" -> {Directive[PointSize[7/360], AbsoluteThickness[2],
                       RGBColor[1, 0, 0], Thickness[Large]], Directive[PointSize[7/360], AbsoluteThickness[2],
                       RGBColor[1, 0, 0], Thickness[Large]]}, "HighlightLabelingFunctions" -> Association[
                      "CoordinatesToolOptions" -> Identity, "ScalingFunctions" -> {{Identity, Identity}, {Identity,
                         Identity}}], "Primitives" -> {}, "GCFlag" -> False], "Meta" -> Association[
                    "DefaultHighlight" -> {"Dynamic", None}, "Index" -> {}, "Function" -> ListLinePlot,
                    "GroupHighlight" -> False]]]], Selectable -> False]}, Annotation[
              {{{}, {}, Annotation[{Hue[0.67, 0.6, 0.6], Directive[PointSize[7/360], AbsoluteThickness[2],
                   RGBColor[1, 0, 0], Thickness[Large]], Line[{{50., -0.5199579745948327}, {55., -0.2930904620581086},
                    {60., -0.1295912744269829}, {65., -0.46049023179894427}, {70., -1.0079395588481899}}]},
                 "Charting`Private`Tag#1"], Annotation[{Hue[0.9060679774997897, 0.6, 0.6], Directive[PointSize[7/360],
                   AbsoluteThickness[2], RGBColor[1, 0, 0], Thickness[Large]], Line[{{50., 0.5726586516259777},
                    {55., 0.31051954984619734}, {60., 0.1352631839851702}, {65., 0.4962605787254412},
                    {70., 1.165984537544339}}]}, "Charting`Private`Tag#2"]}}, Association["HighlightElements" ->
                Association["Label" -> {"XYLabel"}, "Ball" -> {"InterpolatedBall"}], "LayoutOptions" ->
                Association["PanelPlotLayout" -> Association[], "PlotRange" -> {{49.58333333333344, 70.},
                   {-1.0079395588481899, 1.165984537544339}}, "Frame" -> {{False, False}, {False, False}},
                 "AxesOrigin" -> {49.58333333333344, 0}, "ImageSize" -> {360, 360/GoldenRatio}, "Axes" -> {True, True},
                 "LabelStyle" -> {}, "AspectRatio" -> GoldenRatio^(-1), "DefaultStyle" -> {Directive[PointSize[7/360],
                    AbsoluteThickness[2], RGBColor[1, 0, 0], Thickness[Large]], Directive[PointSize[7/360],
                    AbsoluteThickness[2], RGBColor[1, 0, 0], Thickness[Large]]}, "HighlightLabelingFunctions" ->
                  Association["CoordinatesToolOptions" -> Identity, "ScalingFunctions" -> {{Identity, Identity},
                     {Identity, Identity}}], "Primitives" -> {}, "GCFlag" -> False], "Meta" -> Association[
                 "DefaultHighlight" -> {"Dynamic", None}, "Index" -> {}, "Function" -> ListLinePlot,
                 "GroupHighlight" -> False]], "DynamicHighlight"]], {{}, {}}}}, {Frame -> True,
           FrameLabel -> {FormBox["\"k-point index\"", TraditionalForm], FormBox["\"Energy\"", TraditionalForm]},
           PlotRange -> All, ImageSize -> Large, DisplayFunction -> Identity, DisplayFunction -> Identity,
           DisplayFunction -> Identity, AspectRatio -> NCache[GoldenRatio^(-1), 0.6180339887498948],
           Axes -> {True, True}, AxesLabel -> {None, None}, AxesOrigin -> {49.58333333333344, 0},
           DisplayFunction :> Identity, Frame -> {{False, False}, {False, False}},
           FrameLabel -> {{None, None}, {None, None}}, FrameTicks -> {{Automatic, Automatic}, {Automatic, Automatic}},
           GridLines -> {None, None}, GridLinesStyle -> Directive[GrayLevel[0.5, 0.4]],
           Method -> {"AxisPadding" -> Scaled[0.02], "DefaultBoundaryStyle" -> Automatic,
             "DefaultGraphicsInteraction" -> {"Version" -> 1.2, "TrackMousePosition" -> {True, False}, "Effects" ->
                {"Highlight" -> {"ratio" -> 2}, "HighlightPoint" -> {"ratio" -> 2}, "Droplines" ->
                  {"freeformCursorMode" -> True, "placement" -> {"x" -> "All", "y" -> "None"}}}},
             "DefaultMeshStyle" -> AbsolutePointSize[6], "DefaultPlotStyle" -> {Directive[RGBColor[0.368417, 0.506779,
                 0.709798], AbsoluteThickness[2]], Directive[RGBColor[0.880722, 0.611041, 0.142051],
                AbsoluteThickness[2]], Directive[RGBColor[0.560181, 0.691569, 0.194885], AbsoluteThickness[2]],
               Directive[RGBColor[0.922526, 0.385626, 0.209179], AbsoluteThickness[2]], Directive[RGBColor[0.528488,
                 0.470624, 0.701351], AbsoluteThickness[2]], Directive[RGBColor[0.772079, 0.431554, 0.102387],
                AbsoluteThickness[2]], Directive[RGBColor[0.363898, 0.618501, 0.782349], AbsoluteThickness[2]],
               Directive[RGBColor[1, 0.75, 0], AbsoluteThickness[2]], Directive[RGBColor[0.647624, 0.37816, 0.614037],
                AbsoluteThickness[2]], Directive[RGBColor[0.571589, 0.586483, 0.], AbsoluteThickness[2]], Directive[
                RGBColor[0.915, 0.3325, 0.2125], AbsoluteThickness[2]], Directive[RGBColor[0.40082222609352647,
                 0.5220066643438841, 0.85], AbsoluteThickness[2]], Directive[RGBColor[0.9728288904374106,
                 0.621644452187053, 0.07336199581899142], AbsoluteThickness[2]], Directive[RGBColor[0.736782672705901,
                 0.358, 0.5030266573755369], AbsoluteThickness[2]], Directive[RGBColor[0.28026441037696703, 0.715,
                 0.4292089322474965], AbsoluteThickness[2]]}, "DomainPadding" -> Scaled[0.02],
             "PointSizeFunction" -> "SmallPointSize", "RangePadding" -> Scaled[0.05], "OptimizePlotMarkers" -> True,
             "IncludeHighlighting" -> "CurrentPoint", "HighlightStyle" -> Automatic, "OptimizePlotMarkers" -> True,
             "CoordinatesToolOptions" -> {"DisplayFunction" -> ({Identity[#1[[1]]], Identity[#1[[2]]]} & ),
               "CopiedValueFunction" -> ({Identity[#1[[1]]], Identity[#1[[2]]]} & )}},
           PlotRange -> {{49.58333333333344, 70.}, {-1.0025769999999998, 1.157078}}, PlotRangeClipping -> True,
           PlotRangePadding -> {{Scaled[0.02], Scaled[0.02]}, {Scaled[0.05], Scaled[0.05]}},
           Ticks -> {Automatic, Automatic}}], FormBox[FormBox[TemplateBox[{"\"VASP\"", "\"Initial k.p\"",
             "\"Fitted k.p\""}, "LineLegend", DisplayFunction -> (FormBox[StyleBox[StyleBox[PaneBox[
                  TagBox[GridBox[{{TagBox[GridBox[{{GraphicsBox[{{Directive[EdgeForm[Directive[Opacity[0.3], GrayLevel[
                                 0]]], PointSize[0.5], AbsoluteThickness[1.6], RGBColor[0, 0, 1], Thickness[Large]],
                             {LineBox[{{0, 12.5}, {40, 12.5}}]}}, {Directive[EdgeForm[Directive[Opacity[0.3],
                                GrayLevel[0]]], PointSize[0.5], AbsoluteThickness[1.6], RGBColor[0, 0, 1], Thickness[
                               Large]], {}}}, AspectRatio -> Full, ImageSize -> {40, 12.5}, PlotRangePadding -> None,
                           ImagePadding -> Automatic, BaselinePosition -> Scaled[0.18000000000000002] -> Baseline],
                          #1}, {GraphicsBox[{{Directive[EdgeForm[Directive[Opacity[0.3], GrayLevel[0]]], PointSize[
                               0.5], AbsoluteThickness[1.6], GrayLevel[0.5], Dashing[{Small, Small}]],
                             {LineBox[{{0, 12.5}, {40, 12.5}}]}}, {Directive[EdgeForm[Directive[Opacity[0.3],
                                GrayLevel[0]]], PointSize[0.5], AbsoluteThickness[1.6], GrayLevel[0.5], Dashing[
                               {Small, Small}]], {}}}, AspectRatio -> Full, ImageSize -> {40, 12.5},
                           PlotRangePadding -> None, ImagePadding -> Automatic, BaselinePosition ->
                            Scaled[0.18000000000000002] -> Baseline], #2}, {GraphicsBox[{{Directive[EdgeForm[Directive[
                                Opacity[0.3], GrayLevel[0]]], PointSize[0.5], AbsoluteThickness[1.6], RGBColor[1, 0,
                               0], Thickness[Large]], {LineBox[{{0, 12.5}, {40, 12.5}}]}}, {Directive[EdgeForm[
                               Directive[Opacity[0.3], GrayLevel[0]]], PointSize[0.5], AbsoluteThickness[1.6],
                              RGBColor[1, 0, 0], Thickness[Large]], {}}}, AspectRatio -> Full, ImageSize -> {40, 12.5},
                           PlotRangePadding -> None, ImagePadding -> Automatic, BaselinePosition ->
                            Scaled[0.18000000000000002] -> Baseline], #3}}, GridBoxAlignment -> {"Columns" ->
                           {Center, Left}, "Rows" -> {{Baseline}}}, AutoDelete -> False, GridBoxDividers ->
                         {"Columns" -> {{False}}, "Rows" -> {{False}}}, GridBoxItemSize -> {"Columns" -> {{All}},
                          "Rows" -> {{All}}}, GridBoxSpacings -> {"Columns" -> {{0.5}}, "Rows" -> {{0.8}}}], "Grid"]}},
                    GridBoxAlignment -> {"Columns" -> {{Left}}, "Rows" -> {{Top}}}, AutoDelete -> False,
                    GridBoxItemSize -> {"Columns" -> {{Automatic}}, "Rows" -> {{Automatic}}}, GridBoxSpacings ->
                     {"Columns" -> {{1}}, "Rows" -> {{0}}}], "Grid"], Alignment -> Left, AppearanceElements -> None,
                  ImageMargins -> {{5, 5}, {5, 5}}, ImageSizeAction -> "ResizeToFit"], LineIndent -> 0,
                 StripOnInput -> False], {FontFamily -> "Arial"}, Background -> Automatic, StripOnInput -> False],
               TraditionalForm] & ), InterpretationFunction :> (RowBox[{"LineLegend", "[",
                RowBox[{RowBox[{"{", RowBox[{RowBox[{"Directive", "[", RowBox[{TemplateBox[Association["color" ->
                             RGBColor[0, 0, 1]], "RGBColorSwatchTemplate"], ",", RowBox[{"Thickness", "[", "Large",
                            "]"}]}], "]"}], ",", RowBox[{"Directive", "[", RowBox[{TemplateBox[Association["color" ->
                             GrayLevel[0.5]], "GrayLevelColorSwatchTemplate"], ",", RowBox[{"Dashing", "[",
                            RowBox[{"{", RowBox[{"Small", ",", "Small"}], "}"}], "]"}]}], "]"}], ",",
                      RowBox[{"Directive", "[", RowBox[{TemplateBox[Association["color" -> RGBColor[1, 0, 0]],
                           "RGBColorSwatchTemplate"], ",", RowBox[{"Thickness", "[", "Large", "]"}]}], "]"}]}], "}"}],
                  ",", RowBox[{"{", RowBox[{#1, ",", #2, ",", #3}], "}"}]}], "]"}] & ), Editable -> True],
           TraditionalForm], TraditionalForm]}, "Legended", DisplayFunction ->
         (GridBox[{{ItemBox[#2, DefaultBaseStyle -> "LabeledLabel"]},
            {TagBox[ItemBox[PaneBox[TagBox[#1, "SkipImageSizeLevel"], Alignment -> {Center, Baseline},
                BaselinePosition -> Baseline], DefaultBaseStyle -> "Labeled"], "SkipImageSizeLevel"]}},
           GridBoxAlignment -> {"Columns" -> {{Center}}, "Rows" -> {{Center}}}, AutoDelete -> False,
           GridBoxItemSize -> Automatic, BaselinePosition -> {2, 1}] & ), InterpretationFunction ->
         (RowBox[{"Legended", "[", RowBox[{#1, ",", RowBox[{"Placed", "[", RowBox[{#2, ",", "Above"}], "]"}]}],
            "]"}] & ), Editable -> True]], "Output"]}, Open]],
  Cell[CellGroupData[{Cell["相关函数", "SeeAlsoSection"],
     Cell[TextData[{ButtonBox["vaspEig", BaseStyle -> "Link", ButtonData -> "paclet:MagneticKP/ref/vaspEig"], " ◼ ",
        ButtonBox["bandManipulateWithVASP", BaseStyle -> "Link", ButtonData ->
          "paclet:MagneticKP/ref/bandManipulateWithVASP"], " ◼ ", ButtonBox["bandplot", BaseStyle -> "Link",
         ButtonData -> "paclet:MagneticKP/ref/bandplot"]}], "SeeAlso"]}, Open]],
  Cell[CellGroupData[{Cell["相关指南", "MoreAboutSection"],
     Cell[TextData[{ButtonBox["MagneticKP", BaseStyle -> "Link", ButtonData -> "paclet:MagneticKP/guide/MagneticKP"]}],
      "MoreAbout"]}, Open]], Cell[CellGroupData[{Cell["相关教程", "TutorialsSection"],
     Cell[TextData[{ButtonBox["MagneticKP 入门", BaseStyle -> "Link", ButtonData ->
          "paclet:MagneticKP/tutorial/GettingStarted"]}], "Tutorials"]}, Open]],
  Cell[CellGroupData[{Cell["历史", "HistorySection"],
     Cell[TextData[{"引入版本：", Cell["1.1.0", "HistoryData", CellTags -> "New"], " | 修改版本：",
        Cell["1.2.0", "HistoryData", CellTags -> "Modified"], " | 废弃版本：",
        Cell[" ", "HistoryData", CellTags -> "Obsolete"]}], "History"]}, Closed]],
  Cell[CellGroupData[{Cell["分类", "CategorizationSection"], Cell["Symbol", "Categorization",
      CellLabel -> "Entity Type"], Cell["MagneticKP", "Categorization", CellLabel -> "Paclet Name"],
     Cell["MagneticKP`", "Categorization", CellLabel -> "Context"], Cell["MagneticKP/ref/fittingKP", "Categorization",
      CellLabel -> "URI"]}, Closed]], Cell[CellGroupData[{Cell["关键词", "KeywordsSection"],
     Cell["参数拟合", "Keywords"], Cell["VASP", "Keywords"], Cell["最小二乘", "Keywords"],
     Cell["k.p 哈密顿量", "Keywords"]}, Closed]]}, TaggingRules -> <|"Paclet" -> "MagneticKP"|>,
 WindowTitle -> "fittingKP", StyleDefinitions -> FrontEnd`FileName[{"Wolfram"}, "FunctionPageStylesExt.nb",
   CharacterEncoding -> "UTF-8"]]