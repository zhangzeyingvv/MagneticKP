Notebook[{Cell["MagneticKP 入门", "Title"],
  Cell[
   "本教程给出最短的完整流程：载入 MagneticKP、输入对称数据、构造哈密顿量并查看能带。"\
, "Text"], Cell["载入 MagneticKP", "Section"], Cell["安装 MagneticKP Paclet 后载入软件包：", "Text"],
  Cell[BoxData["Needs[\"MagneticKP`\"]"], "Input"], Cell["输入对称数据", "Section"],
  Cell[
   "为每个对称操作给出表示矩阵及其对动量的作用。幺正操作放在 \"Unitary\" 中；只有需要时才添加 \"Antiunitary\"："\
, "Text"], Cell[BoxData["input = <|\n  \"Unitary\" -> <|\n    \"C2z\" -> {\n      {{0, 1}, {-1, 0}},\n      {-kx, -ky, \
kz}\n    },\n    \"Inversion\" -> {\n      {{1, 0}, {0, -1}},\n      {-kx, -ky, -kz}\n    }\n  |>\n|>;"], "Input"],
  Cell["构造哈密顿量", "Section"],
  Cell[
   "使用整数阶数时会包含从零阶到该阶的全部项。返回的 Association 同时保存哈密顿量和基本摘要信息："\
, "Text"], Cell[BoxData["result = kpHam[1, input];\nKeys[result]"], "Input"],
  Cell[BoxData[RowBox[{"{", RowBox[{"\"ham\"", ",", "\"korder\"", ",", "\"dim\"", ",", "\"NumberOfParameters\""}],
      "}"}]], "Output"], Cell["从 \"ham\" 条目中取得矩阵：", "Text"],
  Cell[BoxData["ham = result[\"ham\"];\nMatrixForm[ham]"], "Input"],
  Cell[BoxData[TagBox[RowBox[{"(", "⁠", GridBox[{{SubscriptBox["e", "1"],
          RowBox[{RowBox[{RowBox[{"-", "kx"}], " ", SubscriptBox["t", "1"]}], "-",
            RowBox[{"ky", " ", SubscriptBox["t", "2"]}], "+", RowBox[{"", " ", "kz", " ",
              SubscriptBox["t", "3"]}]}]}, {RowBox[{RowBox[{RowBox[{"-", "kx"}], " ", SubscriptBox["t", "1"]}], "-",
            RowBox[{"ky", " ", SubscriptBox["t", "2"]}], "-", RowBox[{"", " ", "kz", " ",
              SubscriptBox["t", "3"]}]}], SubscriptBox["e", "1"]}}, RowSpacings -> 1, ColumnSpacings -> 1,
        RowAlignments -> Baseline, ColumnAlignments -> Center], "⁠", ")"}],
     Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]], "Output"], Cell["查看模型能带", "Section"],
  Cell["定义带标签的动量路径，并使用 bandManipulate 调整参数值：", "Text"],
  Cell[BoxData["path = {\n  {{{0, 0, 0}, {1/3, 0, 0}}, {\"G\", \"X\"}},\n  {{{0, 0, 0}, {0, 1/3, 0}}, {\"G\", \
\"Y\"}},\n  {{{0, 0, 0}, {0, 0, 1/3}}, {\"G\", \"Z\"}}\n};\nbandManipulate[path, 40, ham, \"PlotRange\" -> {-2, 2}]"],
   "Input"], Cell[BoxData[DynamicBox[Uncompress["1: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"], SingleEvaluation -> True,
     TrackedSymbols :> {}]], "Output"], Cell["使用表示数据库", "Section"],
  Cell[
   "如果表示矩阵可由 MSGCorep 或 SpaceGroupIrep 获得，interfaceRep 会自动生成相同格式的输入："\
, "Text"],
  Cell[
   BoxData[
    "Needs[\"MSGCorep`\"];\ninput = interfaceRep[{13, 65}, \"B\", {2}];\n{input[\"lable\"], input[\"GroupInput\"]}"],
   "Input"], Cell[BoxData[RowBox[{"{", RowBox[{RowBox[{"{", "\"\\!\\(\\*SubscriptBox[\\(B\\), \\(2\\)]\\)\"", "}"}],
        ",", "\"Generators\""}], "}"}]], "Output"], Cell["下一步", "Section"],
  Cell[TextData[{"参阅 ", ButtonBox["kpHam", BaseStyle -> "Link", ButtonData -> "paclet:MagneticKP/ref/kpHam"],
     " 了解全部构造选项；参阅 ", ButtonBox["bandplot", BaseStyle -> "Link",
      ButtonData -> "paclet:MagneticKP/ref/bandplot"], " 绘制固定参数的能带图；结合 ",
     ButtonBox["vaspEig", BaseStyle -> "Link", ButtonData -> "paclet:MagneticKP/ref/vaspEig"], " 与 ",
     ButtonBox["fittingKP", BaseStyle -> "Link", ButtonData -> "paclet:MagneticKP/ref/fittingKP"],
     " 拟合 VASP 数据。"}], "Text"],
  Cell[CellGroupData[{Cell["历史", "HistorySection"],
     Cell[TextData[{"引入版本：", Cell["1.1.0", "HistoryData", CellTags -> "New"], " | 修改版本：",
        Cell["1.2.0", "HistoryData", CellTags -> "Modified"], " | 废弃版本：",
        Cell[" ", "HistoryData", CellTags -> "Obsolete"]}], "History"]}, Closed]],
  Cell[CellGroupData[{Cell["分类", "CategorizationSection"], Cell["Tech Note", "Categorization",
      CellLabel -> "Entity Type"], Cell["MagneticKP", "Categorization", CellLabel -> "Paclet Name"],
     Cell["MagneticKP`", "Categorization", CellLabel -> "Context"], Cell["MagneticKP/tutorial/GettingStarted",
      "Categorization", CellLabel -> "URI"]}, Closed]],
  Cell[CellGroupData[{Cell["关键词", "KeywordsSection"], Cell["入门", "Keywords"],
     Cell["k.p 哈密顿量", "Keywords"], Cell["对称输入", "Keywords"], Cell["能带图", "Keywords"]}, Closed]]},
 TaggingRules -> <|"Paclet" -> "MagneticKP"|>, WindowTitle -> "MagneticKP 入门",
 StyleDefinitions -> FrontEnd`FileName[{"Wolfram"}, "TechNotePageStylesExt.nb", CharacterEncoding -> "UTF-8"]]