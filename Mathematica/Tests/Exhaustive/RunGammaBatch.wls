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

If[Length[arguments] < 4,
  Print[
    "Usage: wolframscript -file ",
    $InputFileName,
    " ",
    "<legacy|current> <start-msg-index> <end-msg-index> <order> [force]"
  ];
  Exit[2]
];

engine = arguments[[1]];
startIndex = ToExpression[arguments[[2]]];
endIndex = ToExpression[arguments[[3]]];
order = ToExpression[arguments[[4]]];
force = Length[arguments] >= 5 && arguments[[5]] === "force";

If[
  ! MemberQ[{"legacy", "current"}, engine] ||
    ! IntegerQ[startIndex] || ! IntegerQ[endIndex] ||
    ! IntegerQ[order] || order < 0,
  Print["Invalid arguments: ", InputForm[arguments]];
  Exit[2]
];

Needs["MSGCorep`"];
Needs["SpaceGroupIrep`"];

legacyRoot = FileNameJoin[{
  ParentDirectory[repositoryRoot],
  "MagneticKP-Mathmatica"
}];

If[engine === "legacy",
  Scan[
    Get[FileNameJoin[{legacyRoot, #}]] &,
    {
      "Usage.wl",
      "Utilities.wl",
      "magkp.wl",
      "IO.m",
      "Fitting.wl"
    }
  ],
  Get[FileNameJoin[{
    mathematicaRoot,
    "MagneticKP",
    "Kernel",
    "init.m"
  }]]
];

outputRoot = FileNameJoin[{
  ParentDirectory[repositoryRoot],
  "tmp",
  "exhaustive_gamma_comparison",
  "order-" <> ToString[order],
  engine
}];
If[! DirectoryQ[outputRoot],
  CreateDirectory[outputRoot, CreateIntermediateDirectories -> True]
];

allMSGs = Sort @ Select[
  Keys[MSGCorep`libMLGCorep],
  ListQ
];
startIndex = Max[1, startIndex];
endIndex = Min[Length[allMSGs], endIndex];

ClearAll[normalMatrix];
normalMatrix[value_] := If[MatrixQ[value], value, {{value}}];

ClearAll[gammaData];
gammaData[msg_List] := Module[{name},
  name = SelectFirst[
    Keys[MSGCorep`libMLGCorep[msg]],
    Function[kname,
      MemberQ[
        Keys[MSGCorep`libMLGCorep[msg, kname]],
        {0, 0, 0}
      ]
    ],
    Missing["NotFound"]
  ];
  If[MissingQ[name],
    $Failed,
    MSGCorep`lookupMLGCorep[msg, name, {0, 0, 0}]
  ]
];

ClearAll[cartesianRotations];
cartesianRotations[msg_List, mlg_List] := Module[
  {lattice, basic, toCartesian},
  lattice = SpaceGroupIrep`getSGLatt[First[msg]];
  basic = SpaceGroupIrep`BasicVectors[lattice];
  toCartesian[matrix_] := FullSimplify[
    Inverse[basic] . matrix . basic
  ];
  toCartesian[
    SpaceGroupIrep`getRotMatOfK[lattice, #[[1]]]
  ] & /@ mlg
];

ClearAll[buildInput];
buildInput[mlg_List, rotations_List, corep_List] := Module[
  {vars, unitaryIndices, antiunitaryIndices, makeEntry},
  vars = {
    MagneticKP`kx,
    MagneticKP`ky,
    MagneticKP`kz
  };
  unitaryIndices = Flatten @ Position[Last /@ mlg, 0];
  antiunitaryIndices = Flatten @ Position[Last /@ mlg, 1];
  makeEntry[index_Integer, antiunitary_] := {
    normalMatrix[corep[[index]]],
    Expand[
      rotations[[index]] .
        If[TrueQ[antiunitary], -vars, vars]
    ]
  };
  <|
    "Unitary" -> Association @ Table[
      ToString[index] -> makeEntry[index, False],
      {index, unitaryIndices}
    ],
    "Anitunitary" -> Association @ Table[
      ToString[index] -> makeEntry[index, True],
      {index, antiunitaryIndices}
    ],
    "GroupInput" -> "Elements"
  |>
];

ClearAll[hamiltonianParameters];
hamiltonianParameters[ham_?MatrixQ] := Module[
  {subscriptParameters, variables},
  subscriptParameters = DeleteDuplicates @ Cases[
    Normal[ham],
    _Subscript,
    Infinity
  ];
  If[subscriptParameters =!= {}, Return[subscriptParameters]];
  variables = Quiet @ Check[
    Variables[Flatten[Normal[ham]]],
    {}
  ];
  DeleteCases[
    variables,
    MagneticKP`kx | MagneticKP`ky | MagneticKP`kz
  ]
];

ClearAll[hamiltonianExponentSet];
hamiltonianExponentSet[ham_?MatrixQ] := Module[
  {vars, rules},
  vars = {
    MagneticKP`kx,
    MagneticKP`ky,
    MagneticKP`kz
  };
  rules = CoefficientRules[Expand[#], vars] & /@
    Flatten[Normal[ham]];
  Sort @ DeleteDuplicates @ Flatten[Keys /@ rules, 1]
];

ClearAll[hamiltonianCoordinateRows];
hamiltonianCoordinateRows[
  ham_?MatrixQ,
  exponents_List
] := Module[
  {vars, params, coefficientMatrix, ruleTables},
  vars = {
    MagneticKP`kx,
    MagneticKP`ky,
    MagneticKP`kz
  };
  params = hamiltonianParameters[ham];
  Table[
    coefficientMatrix = Map[
      Coefficient[#, parameter] &,
      Normal[ham],
      {2}
    ];
    ruleTables = Association @ CoefficientRules[
        Expand[#],
        vars
    ] & /@ Flatten[coefficientMatrix];
    Flatten @ Table[
      Lookup[#, Key[exponent], 0] & /@ ruleTables,
      {exponent, exponents}
    ],
    {parameter, params}
  ]
];

ClearAll[caseSignature];
caseSignature[
  msg_List,
  repIndex_Integer,
  input_Association
] := Module[
  {
    timing, result, resultValue, messages,
    ham, exponents, rows, canonical
  },
  timing = AbsoluteTiming[
    result = CheckAbort[
      Block[{$MessageList = {}},
        resultValue = Quiet @ TimeConstrained[
          MagneticKP`kpHam[
            order,
            input,
            "Method" -> "IterativeSimplification"
          ],
          180,
          $TimedOut
        ];
        messages = $MessageList;
        resultValue
      ],
      $Aborted
    ]
  ][[1]];
  If[! AssociationQ[result],
    Return @ <|
      "MSG" -> msg,
      "RepIndex" -> repIndex,
      "Status" -> ToString[result, InputForm],
      "Messages" -> ToString[messages, InputForm],
      "KPSeconds" -> timing
    |>
  ];
  ham = Lookup[result, "ham", $Failed];
  If[! MatrixQ[ham],
    Return @ <|
      "MSG" -> msg,
      "RepIndex" -> repIndex,
      "Status" -> "NonMatrixHamiltonian",
      "Messages" -> ToString[messages, InputForm],
      "KPSeconds" -> timing
    |>
  ];
  exponents = hamiltonianExponentSet[ham];
  rows = hamiltonianCoordinateRows[ham, exponents];
  canonical = If[rows === {}, {}, FullSimplify @ RowReduce[rows]];
  <|
    "MSG" -> msg,
    "RepIndex" -> repIndex,
    "Status" -> "OK",
    "Dimension" -> Lookup[result, "dim", Length[ham]],
    "NumberOfParameters" ->
      Lookup[result, "NumberOfParameters", Missing["NotFound"]],
    "Exponents" -> exponents,
    "CanonicalRows" -> canonical,
    "Messages" -> ToString[messages, InputForm],
    "KPSeconds" -> timing
  |>
];

Print[
  "Running ", engine, " engine for MSG indices ",
  startIndex, "..", endIndex, " at order ", order
];

Do[
  msg = allMSGs[[msgIndex]];
  outputFile = FileNameJoin[{
    outputRoot,
    StringRiffle[ToString /@ msg, "."] <> ".wxf"
  }];
  If[FileExistsQ[outputFile] && ! force,
    Print[
      "[", msgIndex, "/", Length[allMSGs], "] ",
      InputForm[msg], " skipped"
    ];
    Continue[]
  ];
  data = gammaData[msg];
  If[data === $Failed,
    payload = <|
      "Engine" -> engine,
      "Order" -> order,
      "MSG" -> msg,
      "Status" -> "MissingGamma",
      "Cases" -> {}
    |>;
    Export[outputFile, payload, "WXF"];
    Continue[]
  ];
  mlg = data["MLG"];
  rotations = cartesianRotations[msg, mlg];
  coreps = Join @@ data["corep"];
  groupTiming = AbsoluteTiming[
    signatures = Table[
      input = buildInput[mlg, rotations, coreps[[repIndex]]];
      caseSignature[msg, repIndex, input],
      {repIndex, Length[coreps]}
    ];
  ][[1]];
  payload = <|
    "Engine" -> engine,
    "Order" -> order,
    "MSGIndex" -> msgIndex,
    "MSG" -> msg,
    "CorepCount" -> Length[coreps],
    "GroupSeconds" -> groupTiming,
    "Cases" -> signatures
  |>;
  Export[outputFile, payload, "WXF"];
  Print[
    "[", msgIndex, "/", Length[allMSGs], "] ",
    InputForm[msg], ": ", Length[coreps],
    " coreps, ", NumberForm[groupTiming, {Infinity, 3}], " s"
  ],
  {msgIndex, startIndex, endIndex}
];

Exit[0]
