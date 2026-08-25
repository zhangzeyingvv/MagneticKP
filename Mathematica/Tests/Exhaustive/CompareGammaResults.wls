$HistoryLength = 0;

scriptRoot = ExpandFileName[DirectoryName[$InputFileName]];
testsRoot = ParentDirectory[scriptRoot];
mathematicaRoot = ParentDirectory[testsRoot];
containerRoot = ParentDirectory[mathematicaRoot];
repositoryRoot = If[
  DirectoryQ[FileNameJoin[{containerRoot, "Python-Web"}]],
  ParentDirectory[containerRoot],
  containerRoot
];
arguments = Rest[$ScriptCommandLine];
order = If[arguments === {}, 1, ToExpression[First[arguments]]];

resultRoot = FileNameJoin[{
  ParentDirectory[repositoryRoot],
  "tmp",
  "exhaustive_gamma_comparison",
  "order-" <> ToString[order]
}];
legacyRoot = FileNameJoin[{resultRoot, "legacy"}];
currentRoot = FileNameJoin[{resultRoot, "current"}];

If[! DirectoryQ[legacyRoot] || ! DirectoryQ[currentRoot],
  Print["Missing result directories under ", resultRoot];
  Exit[2]
];

ClearAll[zeroScalarQ, subspaceEqualQ];
zeroScalarQ[value_] := TrueQ[PossibleZeroQ[value]] ||
  TrueQ[FullSimplify[value == 0]];

subspaceEqualQ[a_List, b_List] := Module[
  {rankA, rankB, joinedRank},
  rankA = If[a === {}, 0, MatrixRank[a]];
  rankB = If[b === {}, 0, MatrixRank[b]];
  joinedRank = If[
    a === {} && b === {},
    0,
    MatrixRank[Join[a, b]]
  ];
  TrueQ[rankA == rankB == joinedRank]
];

legacyFiles = FileNames["*.wxf", legacyRoot];
currentFiles = FileNames["*.wxf", currentRoot];
commonNames = Intersection[
  FileNameTake /@ legacyFiles,
  FileNameTake /@ currentFiles
];

comparisons = Reap[
  Do[
    legacy = Import[
      FileNameJoin[{legacyRoot, fileName}],
      "WXF"
    ];
    current = Import[
      FileNameJoin[{currentRoot, fileName}],
      "WXF"
    ];
    legacyCases = Association @ Map[
      ToString[#["RepIndex"], InputForm] -> # &,
      Lookup[legacy, "Cases", {}]
    ];
    currentCases = Association @ Map[
      ToString[#["RepIndex"], InputForm] -> # &,
      Lookup[current, "Cases", {}]
    ];
    repKeys = Union[Keys[legacyCases], Keys[currentCases]];
    Do[
      old = Lookup[legacyCases, repKey, Missing["Legacy"]];
      new = Lookup[currentCases, repKey, Missing["Current"]];
      If[MissingQ[old] || MissingQ[new],
        Sow @ <|
          "MSG" -> Lookup[legacy, "MSG", Lookup[current, "MSG", {}]],
          "RepIndex" -> ToExpression[repKey],
          "Pass" -> False,
          "Reason" -> "MissingCase"
        |>;
        Continue[]
      ];
      statusPass =
        Lookup[old, "Status", ""] === "OK" &&
        Lookup[new, "Status", ""] === "OK";
      dimensionPass = statusPass &&
        Lookup[old, "Dimension", Missing[]] ===
          Lookup[new, "Dimension", Missing[]];
      exponentPass = statusPass &&
        Lookup[old, "Exponents", Missing[]] ===
          Lookup[new, "Exponents", Missing[]];
      metadataPass = dimensionPass && exponentPass;
      parameterCountDifference = statusPass &&
        Lookup[old, "NumberOfParameters", Missing[]] =!=
          Lookup[new, "NumberOfParameters", Missing[]];
      subspacePass = metadataPass &&
        subspaceEqualQ[
          Lookup[old, "CanonicalRows", {}],
          Lookup[new, "CanonicalRows", {}]
        ];
      Sow @ <|
        "MSG" -> Lookup[legacy, "MSG", {}],
        "RepIndex" -> ToExpression[repKey],
        "Pass" -> TrueQ[subspacePass],
        "StatusPass" -> statusPass,
        "DimensionPass" -> dimensionPass,
        "ExponentPass" -> exponentPass,
        "MetadataPass" -> metadataPass,
        "SubspacePass" -> subspacePass,
        "ReportedParameterCountDifference" ->
          parameterCountDifference,
        "LegacyParameters" ->
          Lookup[old, "NumberOfParameters", Missing[]],
        "CurrentParameters" ->
          Lookup[new, "NumberOfParameters", Missing[]],
        "LegacySeconds" -> Lookup[old, "KPSeconds", 0.],
        "CurrentSeconds" -> Lookup[new, "KPSeconds", 0.]
      |>,
      {repKey, repKeys}
    ],
    {fileName, commonNames}
  ]
][[2]];
comparisons = If[comparisons === {}, {}, First[comparisons]];

passed = Count[Lookup[comparisons, "Pass", False], True];
failedCases = Select[comparisons, ! TrueQ[#["Pass"]] &];
parameterCountDifferences = Select[
  comparisons,
  TrueQ[#["ReportedParameterCountDifference"]] &
];
parameterCountDifferenceDetails = KeyTake[
  #,
  {
    "MSG",
    "RepIndex",
    "LegacyParameters",
    "CurrentParameters",
    "SubspacePass"
  }
] & /@ parameterCountDifferences;
legacySeconds = Total @ Lookup[comparisons, "LegacySeconds", 0.];
currentSeconds = Total @ Lookup[comparisons, "CurrentSeconds", 0.];

summary = <|
  "Order" -> order,
  "ComparedMSGs" -> Length[commonNames],
  "ComparedCoreps" -> Length[comparisons],
  "PassedCoreps" -> passed,
  "FailedCoreps" -> Length[failedCases],
  "ReportedParameterCountDifferences" ->
    Length[parameterCountDifferences],
  "ReportedParameterCountDifferenceDetails" ->
    parameterCountDifferenceDetails,
  "LegacyKPSeconds" -> legacySeconds,
  "CurrentKPSeconds" -> currentSeconds,
  "Speedup" -> If[currentSeconds > 0, legacySeconds/currentSeconds, Indeterminate],
  "Failures" -> failedCases
|>;

summaryFile = FileNameJoin[{
  resultRoot,
  "summary.wxf"
}];
Export[summaryFile, summary, "WXF"];
Export[
  FileNameJoin[{resultRoot, "summary.json"}],
  KeyDrop[
    summary,
    {
      "Failures",
      "ReportedParameterCountDifferenceDetails"
    }
  ],
  "RawJSON"
];

Print[
  Dataset @ KeyDrop[
    summary,
    {
      "Failures",
      "ReportedParameterCountDifferenceDetails"
    }
  ]
];
If[failedCases =!= {},
  Print["First failures:"];
  Print[Dataset @ Take[failedCases, UpTo[20]]]
];
If[parameterCountDifferenceDetails =!= {},
  Print[
    "Reported parameter-count differences ",
    "(informational; subspace comparison decides pass/fail):"
  ];
  Print[
    Dataset @ Take[
      parameterCountDifferenceDetails,
      UpTo[20]
    ]
  ]
];
Print["Summary: ", summaryFile];

Exit[If[failedCases === {}, 0, 1]]
