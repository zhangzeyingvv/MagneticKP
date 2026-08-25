Notebook[{Cell[CellGroupData[{Cell["kpHam", "ObjectName"],
     Cell[TextData[{Cell["   ", "ModInfo"], Cell[BoxData["kpHam[korder, input]"], "InlineFormula"], "  ",
        "constructs the symmetry-allowed k.p Hamiltonian through korder and returns its data in an Association."}],
      "Usage"], Cell[TextData[{Cell["   ", "ModInfo"], Cell[BoxData["kpHam[{n1, n2, ...}, input]"], "InlineFormula"],
        "  ", "constructs only the listed homogeneous orders."}], "Usage"], Cell["input is an Association whose \
\"Unitary\" and optional \"Antiunitary\" entries associate symmetry labels with representation matrices and transformed \
momentum vectors.", "Notes"], Cell["The result contains \"ham\", \"korder\", \"dim\", and \"NumberOfParameters\". Use \
result[\"ham\"] to obtain the Hamiltonian matrix.", "Notes"],
     Cell["The historical input key \"Anitunitary\" remains accepted, but new input should use \"Antiunitary\".",
      "Notes"], Cell[BoxData[GridBox[{{Cell["      ", "TableRowIcon"], "\"Method\"", "\"IterativeSimplification\"",
          Cell["selects \"IterativeSimplification\", \"DirectProductDecomposition\", \"Reynolds\", \"Automatic\", or \
the compatibility method \"IterativeSimplification2016\".", "TableText"]}, {Cell["      ", "TableRowIcon"],
          "\"GroupInput\"", "Automatic", Cell["specifies whether the supplied symmetry data contains generators or all \
elements; Automatic uses metadata in input when available.", "TableText"]}, {Cell["      ", "TableRowIcon"],
          "\"MaxGroupOrder\"", "512",
          Cell["sets the largest finite group order accepted while completing generator input.", "TableText"]},
         {Cell["      ", "TableRowIcon"], "\"ReturnDiagnostics\"", "False",
          Cell["specifies whether the returned Association also includes per-order diagnostic information.",
           "TableText"]}}, GridBoxAlignment -> {"Columns" -> {{Left}}}]], "3ColumnTableMod"]}, Open]],
  Cell[CellGroupData[{Cell["Examples Initialization", "ExamplesInitializationSection"],
     Cell[BoxData["Needs[\"MagneticKP`\"]"], "ExampleInitialization"]}, Open]],
  Cell[CellGroupData[{Cell["Basic Examples", "PrimaryExamplesSection"],
     Cell["Construct a first-order two-band model from explicitly supplied symmetry data:", "ExampleText"],
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
    Open]], Cell[CellGroupData[{Cell["Related Functions", "SeeAlsoSection"],
     Cell[TextData[{ButtonBox["interfaceRep", BaseStyle -> "Link", ButtonData -> "paclet:MagneticKP/ref/interfaceRep"],
        " ◼ ", ButtonBox["IterativelySimplify", BaseStyle -> "Link", ButtonData ->
          "paclet:MagneticKP/ref/IterativelySimplify"], " ◼ ", ButtonBox["GETkpHam", BaseStyle -> "Link",
         ButtonData -> "paclet:MagneticKP/ref/GETkpHam"]}], "SeeAlso"]}, Open]],
  Cell[CellGroupData[{Cell["Related Guides", "MoreAboutSection"],
     Cell[TextData[{ButtonBox["MagneticKP", BaseStyle -> "Link", ButtonData -> "paclet:MagneticKP/guide/MagneticKP"]}],
      "MoreAbout"]}, Open]], Cell[CellGroupData[{Cell["Related Tutorials", "TutorialsSection"],
     Cell[TextData[{ButtonBox["Getting Started with MagneticKP", BaseStyle -> "Link",
         ButtonData -> "paclet:MagneticKP/tutorial/GettingStarted"]}], "Tutorials"]}, Open]],
  Cell[CellGroupData[{Cell["History", "HistorySection"],
     Cell[TextData[{"New in: ", Cell["1.1.0", "HistoryData", CellTags -> "New"], " | Modified in: ",
        Cell["1.2.0", "HistoryData", CellTags -> "Modified"], " | Obsolete in: ", Cell[" ", "HistoryData",
         CellTags -> "Obsolete"]}], "History"]}, Closed]],
  Cell[CellGroupData[{Cell["Categorization", "CategorizationSection"], Cell["Symbol", "Categorization",
      CellLabel -> "Entity Type"], Cell["MagneticKP", "Categorization", CellLabel -> "Paclet Name"],
     Cell["MagneticKP`", "Categorization", CellLabel -> "Context"], Cell["MagneticKP/ref/kpHam", "Categorization",
      CellLabel -> "URI"]}, Closed]], Cell[CellGroupData[{Cell["Keywords", "KeywordsSection"],
     Cell["k.p Hamiltonian", "Keywords"], Cell["symmetry", "Keywords"], Cell["effective model", "Keywords"],
     Cell["band theory", "Keywords"]}, Closed]]}, TaggingRules -> <|"Paclet" -> "MagneticKP"|>, WindowTitle -> "kpHam",
 StyleDefinitions -> FrontEnd`FileName[{"Wolfram"}, "FunctionPageStylesExt.nb", CharacterEncoding -> "UTF-8"]]