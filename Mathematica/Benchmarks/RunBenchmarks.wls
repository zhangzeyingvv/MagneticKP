$HistoryLength = 0;

benchmarkRoot = ExpandFileName[DirectoryName[$InputFileName]];
mathematicaRoot = ParentDirectory[benchmarkRoot];
Get[FileNameJoin[{
  mathematicaRoot,
  "Tests",
  "Support",
  "TestHelpers.wl"
}]];
LoadMagneticKPForTests[];

ClearAll[MedianTimed];
SetAttributes[MedianTimed, HoldFirst];
MedianTimed[expr_, repeats_Integer: 5] := Module[
  {times, value},
  value = expr;
  times = Table[
    First @ AbsoluteTiming[value = expr],
    {repeats}
  ];
  <|"Seconds" -> Median[times], "Value" -> value|>
];

ClearAll[MedianTimedWithSetup];
SetAttributes[MedianTimedWithSetup, HoldAll];
MedianTimedWithSetup[
  setup_,
  expr_,
  repeats_Integer: 5
] := Module[{times, value},
  setup;
  value = expr;
  times = Table[
    setup;
    First @ AbsoluteTiming[value = expr],
    {repeats}
  ];
  <|"Seconds" -> Median[times], "Value" -> value|>
];

ClearAll[benchmarkF];
benchmarkF[] := Module[
  {a, orders, rows, warm},
  a = {{1, 1, 0}, {0, 1, 1}, {1, 0, 1}};
  orders = {2, 4, 8};
  rows = Table[
    old = MedianTimed[
      LegacySolveMomentumTransform[order, LinearKImage[a]]
    ];
    new = MedianTimedWithSetup[
      MagneticKP`Private`$MomentumTransformCache = <||>,
      CallMomentumTransform[order, LinearKImage[a]]
    ];
    warm = MedianTimed[
      CallMomentumTransform[order, LinearKImage[a]]
    ];
    <|
      "Component" -> "F",
      "Case" -> "dense-exact-order-" <> ToString[order],
      "Order" -> order,
      "Dimension" -> Length[KMonomialBasis[order, 3]],
      "LegacySeconds" -> old["Seconds"],
      "NewColdSeconds" -> new["Seconds"],
      "NewWarmSeconds" -> warm["Seconds"],
      "ColdSpeedup" -> old["Seconds"]/new["Seconds"],
      "WarmSpeedup" -> old["Seconds"]/warm["Seconds"],
      "EqualQ" -> MatrixExactEqualQ[
        old["Value"], new["Value"]
      ]
    |>,
    {order, orders}
  ];
  rows
];

ClearAll[LegacyAdjointMatrixWithSimplify];
LegacyAdjointMatrixWithSimplify[
  d_?MatrixQ,
  basis_List,
  antiunitary_: False
] := Module[{dinv = Inverse[d], action, inner},
  action[m_] := If[
    TrueQ[antiunitary],
    d . Conjugate[m] . dinv,
    d . m . dinv
  ];
  inner[a_, b_] := FullSimplify[
    Tr[ConjugateTranspose[a] . b]
  ];
  Table[
    inner[action[basis[[i]]], basis[[j]]] /
      inner[basis[[j]], basis[[j]]],
    {j, Length[basis]}, {i, Length[basis]}
  ]
];

ClearAll[benchmarkGg];
benchmarkGg[] := Module[{rows, dimensions, d, warm},
  dimensions = {2, 4, 8};
  rows = Table[
    data = CallOperatorBasisData[n];
    d = N[FourierMatrix[n]];
    old = MedianTimed[
      LegacyAdjointMatrixWithSimplify[
        d, data["Basis"], False
      ]
    ];
    new = MedianTimedWithSetup[
      MagneticKP`Private`$AdjointRepresentationCache = <||>,
      CallAdjointRepresentation[d, data, False]
    ];
    warm = MedianTimed[
      CallAdjointRepresentation[d, data, False]
    ];
    <|
      "Component" -> "Gg",
      "Case" -> "Fourier-unitary-n-" <> ToString[n],
      "StateDimension" -> n,
      "AdjointDimension" -> n^2,
      "LegacySeconds" -> old["Seconds"],
      "NewColdSeconds" -> new["Seconds"],
      "NewWarmSeconds" -> warm["Seconds"],
      "ColdSpeedup" -> old["Seconds"]/new["Seconds"],
      "WarmSpeedup" -> old["Seconds"]/warm["Seconds"],
      "EqualQ" -> MatrixNearQ[
        old["Value"], new["Value"], 10^-11
      ]
    |>,
    {n, dimensions}
  ];
  rows
];

c4z = {{0, -1, 0}, {1, 0, 0}, {0, 0, 1}};
c3Diagonal = {{0, 0, 1}, {1, 0, 0}, {0, 1, 0}};
inversion = -IdentityMatrix[3];
spatialGenerators = {c4z, c3Diagonal, inversion};
matrixGroup =
  MagneticKP`Private`generateMatrixGroup[spatialGenerators];

ClearAll[coefficientRepresentation];
coefficientRepresentation[g_, order_, copies_] := Module[
  {d, basisData, f, j},
  d = KroneckerProduct[g, IdentityMatrix[copies]];
  basisData = CallOperatorBasisData[Length[d]];
  f = CallMomentumTransform[order, LinearKImage[g]];
  j = CallAdjointRepresentation[d, basisData, False];
  SparseArray @ KroneckerProduct[f, j]
];

ClearAll[benchmarkSolverCase];
benchmarkSolverCase[copies_Integer, order_Integer] := Module[
  {
    allBuild, generatorBuild, allReps, generatorReps,
    iterative, reynolds, endToEndIterative, endToEndReynolds,
    ambient
  },
  allBuild = MedianTimed[
    coefficientRepresentation[#, order, copies] & /@ matrixGroup
  ];
  generatorBuild = MedianTimed[
    coefficientRepresentation[#, order, copies] & /@
      spatialGenerators
  ];
  allReps = allBuild["Value"];
  generatorReps = generatorBuild["Value"];
  ambient = Length[First[allReps]];
  iterative = MedianTimed[
    IterativeInvariantBasis[generatorReps]
  ];
  reynolds = MedianTimed[
    CallReynoldsInvariantBasis[allReps]
  ];
  endToEndIterative = MedianTimed[
    IterativeInvariantBasis[
      coefficientRepresentation[#, order, copies] & /@
        spatialGenerators
    ]
  ];
  endToEndReynolds = MedianTimed[
    CallReynoldsInvariantBasis[
      coefficientRepresentation[#, order, copies] & /@
        matrixGroup
    ]
  ];
  {
    <|
      "Component" -> "SolverOnly",
      "Case" -> "Oh-n-" <> ToString[3 copies] <>
        "-order-" <> ToString[order],
      "StateDimension" -> 3 copies,
      "AmbientDimension" -> ambient,
      "GroupOrder" -> Length[matrixGroup],
      "GeneratorCount" -> Length[spatialGenerators],
      "IterativeSeconds" -> iterative["Seconds"],
      "ReynoldsSeconds" -> reynolds["Seconds"],
      "ReynoldsVsIterativeSpeedup" ->
        iterative["Seconds"]/reynolds["Seconds"],
      "EqualQ" -> SubspaceEqualQ[
        iterative["Value"], reynolds["Value"]
      ]
    |>,
    <|
      "Component" -> "EndToEnd",
      "Case" -> "Oh-n-" <> ToString[3 copies] <>
        "-order-" <> ToString[order],
      "StateDimension" -> 3 copies,
      "AmbientDimension" -> ambient,
      "IterativeSeconds" -> endToEndIterative["Seconds"],
      "ReynoldsSeconds" -> endToEndReynolds["Seconds"],
      "ReynoldsVsIterativeSpeedup" ->
        endToEndIterative["Seconds"]/
          endToEndReynolds["Seconds"],
      "EqualQ" -> SubspaceEqualQ[
        endToEndIterative["Value"],
        endToEndReynolds["Value"]
      ]
    |>
  }
];

fResults = benchmarkF[];
gResults = benchmarkGg[];
solverResults = Flatten[
  benchmarkSolverCase[#[[1]], #[[2]]] & /@
    {{1, 2}, {2, 2}, {3, 2}},
  1
];

result = <|
  "WolframVersion" -> $Version,
  "Timestamp" -> DateString["ISODateTime"],
  "Warmup" -> 1,
  "Repeats" -> 5,
  "Results" -> Join[fResults, gResults, solverResults]
|>;

arguments = Rest[$ScriptCommandLine];
jsonPath = If[Length[arguments] >= 1,
  First[arguments],
  FileNameJoin[{
    $TemporaryDirectory,
    "magnetickp-benchmarks.json"
  }]
];
tsvPath = FileNameJoin[{
  DirectoryName[jsonPath],
  FileBaseName[jsonPath] <> ".tsv"
}];

Export[jsonPath, result, "RawJSON"];

allKeys = DeleteDuplicates @ Flatten[Keys /@ result["Results"]];
table = Prepend[
  Lookup[#, allKeys, ""] & /@ result["Results"],
  allKeys
];
Export[tsvPath, table, "TSV"];

Print[Dataset[result["Results"]]];
Print["JSON: ", jsonPath];
Print["TSV: ", tsvPath];
