Notebook[{Cell[CellGroupData[{Cell["bandManipulate", "ObjectName"],
     Cell[TextData[{Cell["   ", "ModInfo"], Cell[BoxData["bandManipulate[path, n, ham]"], "InlineFormula"], "  ",
        "创建交互式能带图，并为 ham 中的自由参数提供控件。"}], "Usage"],
     Cell["每段路径的格式为 {{k1, k2}, {label1, label2}}，n 是每段路径的采样数。", "Notes"],
     Cell["path 中的动量坐标在代入 ham 前会乘以 2 Pi。", "Notes"],
     Cell[BoxData[GridBox[{{Cell["      ", "TableRowIcon"], "\"PlotRange\"", "All",
          Cell["设置显示的能量范围。", "TableText"]}, {Cell["      ", "TableRowIcon"], "\"ManipulateRange\"",
          "1", Cell["设置对称的滑块范围，也可以明确给出 {min, max}。", "TableText"]}},
        GridBoxAlignment -> {"Columns" -> {{Left}}}]], "3ColumnTableMod"]}, Open]],
  Cell[CellGroupData[{Cell["范例初始化", "ExamplesInitializationSection"],
     Cell[BoxData["Needs[\"MagneticKP`\"]"], "ExampleInitialization"]}, Open]],
  Cell[CellGroupData[{Cell["基本范例", "PrimaryExamplesSection"],
     Cell["先用 kpHam 构造双带哈密顿量，再通过交互式参数控件查看其能带：", "ExampleText"],
     Cell[BoxData["input = <|\n  \"Unitary\" -> <|\n    \"C2z\" -> {\n      {{0, 1}, {-1, 0}},\n      {-kx, -ky, kz}\n  \
  },\n    \"Inversion\" -> {\n      {{1, 0}, {0, -1}},\n      {-kx, -ky, -kz}\n    }\n  |>\n|>;\nham = kpHam[1, \
input][\"ham\"];\npath = {{{{0, 0, 0}, {1/3, 0, 0}}, {\"G\", \"X\"}}};\nbandManipulate[path, 20, ham, \"PlotRange\" -> \
{-2, 2}]"], "Input"], Cell[BoxData[DynamicBox[Uncompress["1:eJztWE2IHEUU7p1JdOMGdzVG1EVhQx+ChxA35JSA7mY2ZIKbXWYmgdy2Zrp\
mt5iaqqaret1ZQU96ETx48aaeFDyo6EXiQXLyJxAQzEUUPCj+IGERRZSQtV737HZXT3XX7AaWgO5htqvre6++9+rVe/V6qslr7ZLjOOIe9dNAy7N8rb0Pxq\
Pqpy57FMMbGIgJ9VPpMdQlrXnuhfHMOMyAwLNEyFhTGSSxFPvVf9x6ynWJenAG52Th3HTB3Anj3P1gQc/HAsslscKfc92IWCMIsQ58OAVsct7poqAD9F237\
eSbY5JSfsCuWweBecxC8WAKgxjpIskDpVXT80AGs4oUhIwUmqOcDuvcp14dCZkkkmLviA4fT1vv45ZQ6+7XrNmnjcraaGRrdI5Tz7hzB2CdsClaAfGlADyO\
SJNN9ZeiD0pqKjjEmHo4w5kMOI14wXheGeyHFNyyQzrZYBmkI/eWzrSVTmkv6Zyw0infCZ0CAkn6mA2l5CyK0rk1nweygiSKwUByMSAsG4Uai+JAS6JyOCm\
ZSMmdSJUSqenhpcqJlNoKMZX2H16qhF1/qbFCRFVc4HKm72nLHsDEBc5wQZYj6yotRHxmQslVuiEtHT2RQoe+p7gYU2c6bxFGZEUtqxSDPWcRFQX5U3SIX0\
0k9IQLkdCvGFArgKd4TPdM8qgQWMR7dhAWkoh5KPDO8qAbaxvrO6qCKephrw7qLqGAoKZynimu9roQZQkeSjaU4iCh6gxSzYrCLi/IlbTUWFrKVOzyy1pB7\
copWeYqZahJprpjilFTIOaFXX58mZ0FTpnlXi/OUqN9Fy1S3s9UB6KgREz4fCuUYfsa4NcYARrnyDJmq4iGytcjW4rm4IUim8RzDfsUtfAMTR1fPZUlo7xU\
lqTTRRqKFCPSxdlhlLIFvOusFac6k2DmXaynV5z82mVdBmaqseR6cQLM+ODuttKk5w4tz93vbb+c54T9N0pg5sqgl7r0RkcBgoJYS8d0OUlvZwZdykGvG9G\
DMRpPtA+b4xZGgdP/G3zIAH/86PT4zetXnh4S+KUd+Okb718+Nvm1HfjxtRe+++X6N3bg5Vtvf1C7+r0deOO1V4Pjkz/YgadOvj55MPxpWI6/2oEnm+03Pz\
t60w7kGy990bi6YQe+90nv3TPsdzvwt5e7r0xP/mEHPllDdOrzP+1A8vjF+YnwLzvwnZ/PHduc+nvY6PnHDnRfPL7x7fO37MDmefera0dv24FvHX7kwys3N\
ofdGecZA1A/rJBHoFDXEFvG+tEkt6FXMRzwLZno60SUcGZVRe7kX6Tq6lpC2qSlrseciUxXXM4f7bYPtqncRS9rV7njftSucnc95f+doFYGjVF5rxov+HE4\
OoNuMAVxBbdRSGWB1PYBgbt2nbBlivvXV5DQO7PMCgJaFe3b3qX4ImwnF38WxHB9J4iS9WS1pHfNisC9/6Lqi9SpxV7S4YykFzM0ERkzHwUze6y1EnDGQ1E\
dZJBv70TEQKDV/PXNzYlZ3UOROo83siaV8k3Kb9r6XxXh9CTtccb40cjrPuWwvG7rNiZuo1VoLrAq80OZ6ucjU5/QFliqMokDP8AyeRcdnH8BTioCCQ=="]\
, SingleEvaluation -> True, TrackedSymbols :> {}]], "Output"]}, Open]],
  Cell[CellGroupData[{Cell["相关函数", "SeeAlsoSection"],
     Cell[TextData[{ButtonBox["bandplot", BaseStyle -> "Link", ButtonData -> "paclet:MagneticKP/ref/bandplot"],
        " ◼ ", ButtonBox["bandManipulateWithVASP", BaseStyle -> "Link",
         ButtonData -> "paclet:MagneticKP/ref/bandManipulateWithVASP"], " ◼ ",
        ButtonBox["kpHam", BaseStyle -> "Link", ButtonData -> "paclet:MagneticKP/ref/kpHam"]}], "SeeAlso"]}, Open]],
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
     Cell["MagneticKP`", "Categorization", CellLabel -> "Context"], Cell["MagneticKP/ref/bandManipulate",
      "Categorization", CellLabel -> "URI"]}, Closed]],
  Cell[CellGroupData[{Cell["关键词", "KeywordsSection"], Cell["能带图", "Keywords"],
     Cell["交互式", "Keywords"], Cell["参数", "Keywords"], Cell["动量路径", "Keywords"]}, Closed]]},
 TaggingRules -> <|"Paclet" -> "MagneticKP"|>, WindowTitle -> "bandManipulate",
 StyleDefinitions -> FrontEnd`FileName[{"Wolfram"}, "FunctionPageStylesExt.nb", CharacterEncoding -> "UTF-8"]]