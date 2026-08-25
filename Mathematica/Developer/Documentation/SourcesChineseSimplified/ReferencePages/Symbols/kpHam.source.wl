Notebook[{Cell[CellGroupData[{Cell["kpHam", "ObjectName"],
     Cell[TextData[{Cell["   ", "ModInfo"], Cell[BoxData["kpHam[korder, input]"], "InlineFormula"], "  ",
        "构造截至 korder 阶、满足对称性要求的 k.p 哈密顿量，并以 Association 返回结果。"}],
      "Usage"], Cell[TextData[{Cell["   ", "ModInfo"], Cell[BoxData["kpHam[{n1, n2, ...}, input]"], "InlineFormula"],
        "  ", "仅构造列表中指定的齐次阶数。"}], "Usage"], Cell["input 是一个 Association，其中 \"Unitary\" \
和可选的 \"Antiunitary\" 将对称操作标签与表示矩阵及变换后的动量向量关联起来。", "Notes"],
     Cell[
      "返回结果包含 \"ham\"、\"korder\"、\"dim\" 和 \"NumberOfParameters\"。使用 result[\"ham\"] 取得哈密顿量矩阵。"\
, "Notes"], Cell["仍兼容历史拼写 \"Anitunitary\"，但新输入应使用 \"Antiunitary\"。", "Notes"],
     Cell[BoxData[GridBox[{{Cell["      ", "TableRowIcon"], "\"Method\"", "\"IterativeSimplification\"",
          Cell["选择 \"IterativeSimplification\"、\"DirectProductDecomposition\"、\"Reynolds\"、\"Automatic\" 或兼容方法 \
\"IterativeSimplification2016\"。", "TableText"]}, {Cell["      ", "TableRowIcon"], "\"GroupInput\"", "Automatic",
          Cell[
           "指定输入的对称数据是生成元还是全部群元；Automatic 会优先读取 input 中的元数据。"\
, "TableText"]}, {Cell["      ", "TableRowIcon"], "\"MaxGroupOrder\"", "512",
          Cell["设置补全生成元输入时允许的最大有限群阶数。", "TableText"]},
         {Cell["      ", "TableRowIcon"], "\"ReturnDiagnostics\"", "False",
          Cell["指定返回的 Association 是否同时包含逐阶诊断信息。", "TableText"]}},
        GridBoxAlignment -> {"Columns" -> {{Left}}}]], "3ColumnTableMod"]}, Open]],
  Cell[CellGroupData[{Cell["范例初始化", "ExamplesInitializationSection"],
     Cell[BoxData["Needs[\"MagneticKP`\"]"], "ExampleInitialization"]}, Open]],
  Cell[CellGroupData[{Cell["基本范例", "PrimaryExamplesSection"],
     Cell["根据手动给出的对称数据构造一阶双带模型：", "ExampleText"],
     Cell[BoxData["input = <|\n  \"Unitary\" -> <|\n    \"C2z\" -> {\n      {{0, 1}, {-1, 0}},\n      {-kx, -ky, kz}\n  \
  },\n    \"Inversion\" -> {\n      {{1, 0}, {0, -1}},\n      {-kx, -ky, -kz}\n    }\n  |>\n|>;\nresult = kpHam[1, \
input];\nMatrixForm[result[\"ham\"]]"], "Input"],
     Cell[BoxData[TagBox[RowBox[{"(", "⁠", GridBox[{{SubscriptBox["e", "1"],
             RowBox[{RowBox[{RowBox[{"-", "kx"}], " ", SubscriptBox["t", "1"]}], "-", RowBox[{"ky", " ",
                 SubscriptBox["t", "2"]}], "+", RowBox[{"", " ", "kz", " ", SubscriptBox["t", "3"]}]}]},
            {RowBox[{RowBox[{RowBox[{"-", "kx"}], " ", SubscriptBox["t", "1"]}], "-", RowBox[{"ky", " ",
                 SubscriptBox["t", "2"]}], "-", RowBox[{"", " ", "kz", " ", SubscriptBox["t", "3"]}]}],
             SubscriptBox["e", "1"]}}, RowSpacings -> 1, ColumnSpacings -> 1, RowAlignments -> Baseline,
           ColumnAlignments -> Center], "⁠", ")"}], Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]], "Output"]},
    Open]], Cell[CellGroupData[{Cell["相关函数", "SeeAlsoSection"],
     Cell[TextData[{ButtonBox["interfaceRep", BaseStyle -> "Link", ButtonData -> "paclet:MagneticKP/ref/interfaceRep"],
        " ◼ ", ButtonBox["IterativelySimplify", BaseStyle -> "Link", ButtonData ->
          "paclet:MagneticKP/ref/IterativelySimplify"], " ◼ ", ButtonBox["GETkpHam", BaseStyle -> "Link",
         ButtonData -> "paclet:MagneticKP/ref/GETkpHam"]}], "SeeAlso"]}, Open]],
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
     Cell["MagneticKP`", "Categorization", CellLabel -> "Context"], Cell["MagneticKP/ref/kpHam", "Categorization",
      CellLabel -> "URI"]}, Closed]], Cell[CellGroupData[{Cell["关键词", "KeywordsSection"],
     Cell["k.p 哈密顿量", "Keywords"], Cell["对称性", "Keywords"], Cell["有效模型", "Keywords"],
     Cell["能带理论", "Keywords"]}, Closed]]}, TaggingRules -> <|"Paclet" -> "MagneticKP"|>,
 WindowTitle -> "kpHam", StyleDefinitions -> FrontEnd`FileName[{"Wolfram"}, "FunctionPageStylesExt.nb",
   CharacterEncoding -> "UTF-8"]]