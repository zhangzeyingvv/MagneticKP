(* Run the frozen Mathematica implementation in an isolated Wolfram kernel. *)

$HistoryLength = 0;

testsRoot = ExpandFileName[
  ParentDirectory[DirectoryName[$InputFileName]]
];
mathematicaRoot = ParentDirectory[testsRoot];
containerRoot = ParentDirectory[mathematicaRoot];
repositoryRoot = If[
  DirectoryQ[FileNameJoin[{containerRoot, "Python-Web"}]],
  ParentDirectory[containerRoot],
  containerRoot
];
legacyRoot = FileNameJoin[{
  ParentDirectory[repositoryRoot],
  "MagneticKP-Mathmatica"
}];

If[! DirectoryQ[legacyRoot],
  Print["Legacy root not found: ", legacyRoot];
  Exit[2]
];

Needs["MSGCorep`"];
Needs["SpaceGroupIrep`"];

Scan[
  Get[FileNameJoin[{legacyRoot, #}]] &,
  {
    "Usage.wl",
    "Utilities.wl",
    "magkp.wl",
    "IO.m",
    "Fitting.wl"
  }
];

ClearAll[runCase];
runCase[spec_Association] := Module[
  {generatorInput, elementInput, order},
  order = spec["Order"];
  generatorInput = MagneticKP`interfaceRep[
    spec["MSG"],
    spec["K"],
    spec["Reps"],
    "CalculateGenerators" -> True
  ];
  elementInput = MagneticKP`interfaceRep[
    spec["MSG"],
    spec["K"],
    spec["Reps"],
    "CalculateGenerators" -> False
  ];
  <|
    "GeneratorInputSummary" -> <|
      "UnitaryCount" -> Length[generatorInput["Unitary"]],
      "AntiunitaryCount" ->
        Length[generatorInput["Anitunitary"]]
    |>,
    "ElementInputSummary" -> <|
      "UnitaryCount" -> Length[elementInput["Unitary"]],
      "AntiunitaryCount" ->
        Length[elementInput["Anitunitary"]]
    |>,
    "IterativeGenerators" ->
      MagneticKP`kpHam[
        order,
        generatorInput,
        "Method" -> "IterativeSimplification"
      ],
    "IterativeElements" ->
      MagneticKP`kpHam[
        order,
        elementInput,
        "Method" -> "IterativeSimplification"
      ]
  |>
];

caseSpecifications = <|
  "MSG13.65-B-2" -> <|
    "MSG" -> {13, 65}, "K" -> "B", "Reps" -> {2},
    "Order" -> 1
  |>,
  "MSG1.2-Gamma-1" -> <|
    "MSG" -> {1, 2}, "K" -> {0, 0, 0}, "Reps" -> {1},
    "Order" -> 2
  |>,
  "MSG1.2-Gamma-2" -> <|
    "MSG" -> {1, 2}, "K" -> {0, 0, 0}, "Reps" -> {2},
    "Order" -> 1
  |>,
  "MSG79.27-Gamma-1" -> <|
    "MSG" -> {79, 27}, "K" -> {0, 0, 0}, "Reps" -> {1},
    "Order" -> 1
  |>,
  "MSG97.155-P-1+2" -> <|
    "MSG" -> {97, 155}, "K" -> "P", "Reps" -> {1, 2},
    "Order" -> 1
  |>,
  "MSG80.32-Gamma-1" -> <|
    "MSG" -> {80, 32}, "K" -> {0, 0, 0}, "Reps" -> {1},
    "Order" -> 1
  |>,
  "MSG75.4-R-1+2" -> <|
    "MSG" -> {75, 4}, "K" -> "R", "Reps" -> {1, 2},
    "Order" -> 1
  |>,
  "MSG164.90-Gamma-5+6+7+8" -> <|
    "MSG" -> {164, 90}, "K" -> {0, 0, 0},
    "Reps" -> {5, 6, 7, 8}, "Order" -> 0
  |>
|>;

msgCases = Association @ KeyValueMap[
  #1 -> runCase[#2] &,
  caseSpecifications
];

dimensionReferences = Association @ Table[
  input = MagneticKP`interfaceRep[
    {1, 1},
    {0, 0, 0},
    ConstantArray[1, n],
    "CalculateGenerators" -> False
  ];
  ToString[n] -> MagneticKP`GETkpHam[0, input],
  {n, 1, 10}
];

type3Input = <|
  "Unitary" -> <|
    "E" -> {
      IdentityMatrix[2],
      {MagneticKP`kx, MagneticKP`ky, MagneticKP`kz}
    }
  |>,
  "AntisymmetryUnitaryTest" -> <|
    "A" -> {
      PauliMatrix[3],
      {MagneticKP`kx, MagneticKP`ky, MagneticKP`kz}
    }
  |>
|>;

payload = <|
  "LegacyRoot" -> legacyRoot,
  "WolframVersion" -> $Version,
  "Cases" -> msgCases,
  "DimensionReferencesGETkpHam" -> dimensionReferences,
  "DirectProductType3IterativeReference" ->
    MagneticKP`IterativelySimplify[0, type3Input]
|>;

Print["MAGNETICKP_ORACLE_BEGIN"];
Print[Compress[payload]];
Print["MAGNETICKP_ORACLE_END"];

Exit[0]
