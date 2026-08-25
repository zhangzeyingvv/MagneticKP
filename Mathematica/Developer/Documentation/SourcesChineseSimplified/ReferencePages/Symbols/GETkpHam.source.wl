Notebook[{Cell[CellGroupData[{Cell["GETkpHam", "ObjectName"],
     Cell[TextData[{Cell["   ", "ModInfo"], Cell[BoxData["GETkpHam[n, input]"], "InlineFormula"], "  ",
        "按照旧版调用方式，仅构造 k.p 哈密顿量的 n 阶齐次部分。"}], "Usage"],
     Cell["GETkpHam 用于兼容早期 MagneticKP 版本编写的 Notebook。", "Notes"],
     Cell["新工作应使用 kpHam，特别是需要多个阶数或诊断输出时。", "Notes"]}, Open]],
  Cell[CellGroupData[{Cell["范例初始化", "ExamplesInitializationSection"],
     Cell[BoxData["Needs[\"MagneticKP`\"]"], "ExampleInitialization"]}, Open]],
  Cell[CellGroupData[{Cell["基本范例", "PrimaryExamplesSection"], Cell["仅构造一阶贡献：", "ExampleText"],
     Cell[BoxData["input = <|\n  \"Unitary\" -> <|\n    \"C2z\" -> {\n      {{0, 1}, {-1, 0}},\n      {-kx, -ky, kz}\n  \
  },\n    \"Inversion\" -> {\n      {{1, 0}, {0, -1}},\n      {-kx, -ky, -kz}\n    }\n  |>\n|>;\nGETkpHam[1, \
input][\"ham\"] // MatrixForm"], "Input"],
     Cell[BoxData[TagBox[RowBox[{"(", "⁠", GridBox[{{"0", RowBox[{RowBox[{RowBox[{"-", "kx"}], " ",
                 SubscriptBox["t", "1"]}], "-", RowBox[{"ky", " ", SubscriptBox["t", "2"]}], "+", RowBox[
                {"", " ", "kz", " ", SubscriptBox["t", "3"]}]}]},
            {RowBox[{RowBox[{RowBox[{"-", "kx"}], " ", SubscriptBox["t", "1"]}], "-", RowBox[{"ky", " ",
                 SubscriptBox["t", "2"]}], "-", RowBox[{"", " ", "kz", " ", SubscriptBox["t", "3"]}]}], "0"}},
           RowSpacings -> 1, ColumnSpacings -> 1, RowAlignments -> Baseline, ColumnAlignments -> Center], "⁠", ")"}],
        Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]], "Output"]}, Open]],
  Cell[CellGroupData[{Cell["相关函数", "SeeAlsoSection"],
     Cell[TextData[{ButtonBox["kpHam", BaseStyle -> "Link", ButtonData -> "paclet:MagneticKP/ref/kpHam"], " ◼ ",
        ButtonBox["IterativelySimplify", BaseStyle -> "Link", ButtonData ->
          "paclet:MagneticKP/ref/IterativelySimplify"]}], "SeeAlso"]}, Open]],
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
     Cell["MagneticKP`", "Categorization", CellLabel -> "Context"], Cell["MagneticKP/ref/GETkpHam", "Categorization",
      CellLabel -> "URI"]}, Closed]], Cell[CellGroupData[{Cell["关键词", "KeywordsSection"],
     Cell["旧版接口", "Keywords"], Cell["齐次阶数", "Keywords"], Cell["k.p 哈密顿量", "Keywords"]},
    Closed]]}, TaggingRules -> <|"Paclet" -> "MagneticKP"|>, WindowTitle -> "GETkpHam",
 StyleDefinitions -> FrontEnd`FileName[{"Wolfram"}, "FunctionPageStylesExt.nb", CharacterEncoding -> "UTF-8"]]