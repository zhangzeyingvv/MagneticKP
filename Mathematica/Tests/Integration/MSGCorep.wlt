If[! ValueQ[$MagneticKPTestHelpersLoaded],
  Get[FileNameJoin[{
    ParentDirectory[DirectoryName[$InputFileName]],
    "Support",
    "TestHelpers.wl"
  }]]
];
LoadMagneticKPForTests[];

msgCorepAvailable = LoadMSGCorepForTests[];
legacyOracle = GetLegacyOracleData[];

VerificationTest[
  msgCorepAvailable,
  True,
  TestID -> "MSGCorep-and-SpaceGroupIrep-available"
]

msgCaseSpecifications = <|
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

ClearAll[EvaluateMSGCase];
EvaluateMSGCase[name_String] := Module[
  {
    spec, old, generatorInput, elementInput,
    newGenerator, newElement,
    reynoldsGenerator, reynoldsElement,
    results, referenceHam, subspaceChecks, parameterCounts,
    resolvedChecks
  },
  spec = msgCaseSpecifications[name];
  old = legacyOracle["Cases", name];
  generatorInput = MagneticKP`interfaceRep[
    spec["MSG"], spec["K"], spec["Reps"],
    "CalculateGenerators" -> True
  ];
  elementInput = MagneticKP`interfaceRep[
    spec["MSG"], spec["K"], spec["Reps"],
    "CalculateGenerators" -> False
  ];
  newGenerator = RunIterativeKpHam[spec["Order"], generatorInput];
  newElement = RunIterativeKpHam[spec["Order"], elementInput];
  reynoldsGenerator = RunReynoldsKpHam[
    spec["Order"], generatorInput, "Generators"
  ];
  reynoldsElement = RunReynoldsKpHam[
    spec["Order"], elementInput, "Elements"
  ];
  results = {
    old["IterativeGenerators"],
    old["IterativeElements"],
    newGenerator,
    newElement,
    reynoldsGenerator,
    reynoldsElement
  };
  If[! And @@ (AssociationQ /@ results),
    Return[<|
      "Pass" -> False,
      "Reason" -> "NonAssociationResult",
      "AssociationQ" -> (AssociationQ /@ results)
    |>]
  ];
  referenceHam = results[[1, "ham"]];
  subspaceChecks =
    HamiltonianSubspaceEqualQ[referenceHam, #["ham"]] & /@
      Rest[results];
  parameterCounts = Lookup[results, "NumberOfParameters"];
  resolvedChecks = ResolvedHamiltonianQ[#["ham"]] & /@ results;
  <|
    "Pass" ->
      SameQ @@ parameterCounts &&
      And @@ subspaceChecks &&
      And @@ resolvedChecks,
    "ParameterCounts" -> parameterCounts,
    "SubspaceChecks" -> subspaceChecks,
    "ResolvedHamiltonians" -> resolvedChecks,
    "GeneratorCounts" -> {
      Length[generatorInput["Unitary"]],
      Length[generatorInput["Anitunitary"]]
    },
    "ElementCounts" -> {
      Length[elementInput["Unitary"]],
      Length[elementInput["Anitunitary"]]
    }
  |>
];

$MSGCorepCaseDiagnostics = Association @ KeyValueMap[
  #1 -> EvaluateMSGCase[#1] &,
  msgCaseSpecifications
];

VerificationTest[
  $MSGCorepCaseDiagnostics["MSG13.65-B-2", "Pass"],
  True,
  TestID -> "MSGCorep-MSG13.65-B-rep2-all-solvers-parity"
]

VerificationTest[
  $MSGCorepCaseDiagnostics["MSG1.2-Gamma-1", "Pass"],
  True,
  TestID -> "MSGCorep-MSG1.2-Gamma-rep1-all-solvers-parity"
]

VerificationTest[
  $MSGCorepCaseDiagnostics["MSG1.2-Gamma-2", "Pass"],
  True,
  TestID -> "MSGCorep-MSG1.2-Gamma-rep2-all-solvers-parity"
]

VerificationTest[
  $MSGCorepCaseDiagnostics["MSG79.27-Gamma-1", "Pass"],
  True,
  TestID -> "MSGCorep-MSG79.27-Gamma-rep1-all-solvers-parity"
]

VerificationTest[
  $MSGCorepCaseDiagnostics["MSG97.155-P-1+2", "Pass"],
  True,
  TestID -> "MSGCorep-MSG97.155-P-reps1+2-all-solvers-parity"
]

VerificationTest[
  $MSGCorepCaseDiagnostics["MSG80.32-Gamma-1", "Pass"],
  True,
  TestID -> "MSGCorep-MSG80.32-Gamma-rep1-all-solvers-parity"
]

VerificationTest[
  $MSGCorepCaseDiagnostics["MSG75.4-R-1+2", "Pass"],
  True,
  TestID -> "MSGCorep-MSG75.4-R-reps1+2-all-solvers-parity"
]

VerificationTest[
  $MSGCorepCaseDiagnostics[
    "MSG164.90-Gamma-5+6+7+8",
    "Pass"
  ],
  True,
  TestID -> "MSGCorep-MSG164.90-Gamma-reps5+6+7+8-all-solvers-parity"
]

centralExtensionElements = MagneticKP`interfaceRep[
  {164, 90},
  {0, 0, 0},
  {5, 6, 7, 8},
  "CalculateGenerators" -> False
];
centralExtensionReynolds = Quiet[
  MagneticKP`kpHam[
    0,
    centralExtensionElements,
    "Method" -> "Reynolds",
    "GroupInput" -> "Elements",
    "ReturnDiagnostics" -> True
  ]
];

VerificationTest[
  Module[{diagnostics},
    If[! AssociationQ[centralExtensionReynolds],
      Return[False]
    ];
    diagnostics = Last[
      centralExtensionReynolds["Diagnostics"]
    ];
    diagnostics["ReynoldsFallback"] === False &&
      diagnostics["ElementClosureCompleted"] === True &&
      diagnostics["GroupOrder"] == 48
  ],
  True,
  TestID -> "MSGCorep-MSG164.90-elements-complete-to-48-factor-group"
]

msg200GeneratorInput = MagneticKP`interfaceRep[
  {200, 14},
  {0, 0, 0},
  {12},
  "CalculateGenerators" -> True
];
msg200ElementInput = MagneticKP`interfaceRep[
  {200, 14},
  {0, 0, 0},
  {12},
  "CalculateGenerators" -> False
];
msg200GeneratorResult = MagneticKP`kpHam[
  1,
  msg200GeneratorInput,
  "Method" -> "IterativeSimplification",
  "ReturnDiagnostics" -> True
];
msg200ElementResult = MagneticKP`kpHam[
  1,
  msg200ElementInput,
  "Method" -> "IterativeSimplification",
  "ReturnDiagnostics" -> True
];
msg200Results = {
  msg200GeneratorResult,
  msg200ElementResult
};

VerificationTest[
  And @@ (AssociationQ /@ msg200Results) &&
    And @@ (
      Lookup[#["Diagnostics"], "InvariantDimension"] === {1, 0} &
        /@ msg200Results
    ),
  True,
  TestID -> "MSGCorep-MSG200.14-Gamma-rep12-order-dimensions"
]

VerificationTest[
  And @@ (
    #["NumberOfParameters"] == 1 &&
      Length[HamiltonianParameters[#["ham"]]] == 1 &
      /@ msg200Results
  ),
  True,
  TestID -> "MSGCorep-MSG200.14-reported-count-equals-actual-rank"
]

msg200ReferenceHamiltonian =
  msg200ReferenceCoefficient IdentityMatrix[2];

VerificationTest[
  And @@ (
    HamiltonianSubspaceEqualQ[
      #["ham"],
      msg200ReferenceHamiltonian
    ] & /@ msg200Results
  ),
  True,
  TestID -> "MSGCorep-MSG200.14-Hamiltonian-space-is-identity"
]

VerificationTest[
  HamiltonianSubspaceEqualQ[
    msg200GeneratorResult["ham"],
    msg200ElementResult["ham"]
  ],
  True,
  TestID -> "MSGCorep-MSG200.14-generators-elements-subspace-parity"
]

ClearAll[EvaluateMSGDimensionCase];
EvaluateMSGDimensionCase[n_Integer] := Module[
  {
    input, oldReference, iterative, reynolds,
    expected = n^2
  },
  input = MagneticKP`interfaceRep[
    {1, 1},
    {0, 0, 0},
    ConstantArray[1, n],
    "CalculateGenerators" -> False
  ];
  oldReference =
    legacyOracle["DimensionReferencesGETkpHam", ToString[n]];
  iterative = RunIterativeKpHam[0, input];
  reynolds = RunReynoldsKpHam[0, input, "Elements"];
  If[! And @@ (AssociationQ /@
      {oldReference, iterative, reynolds}),
    Return[<|
      "Pass" -> False,
      "Reason" -> "NonAssociationResult"
    |>]
  ];
  <|
    "Pass" ->
      oldReference["NumberOfParameters"] == expected &&
      iterative["NumberOfParameters"] == expected &&
      reynolds["NumberOfParameters"] == expected &&
      HamiltonianSubspaceEqualQ[
        oldReference["ham"], iterative["ham"]
      ] &&
      HamiltonianSubspaceEqualQ[
        oldReference["ham"], reynolds["ham"]
      ],
    "Counts" -> {
      oldReference["NumberOfParameters"],
      iterative["NumberOfParameters"],
      reynolds["NumberOfParameters"]
    }
  |>
];

$MSGCorepDimensionDiagnostics = Association @ Table[
  ToString[n] -> EvaluateMSGDimensionCase[n],
  {n, 1, 10}
];

VerificationTest[
  And @@ Lookup[
    Values[$MSGCorepDimensionDiagnostics],
    "Pass"
  ],
  True,
  TestID -> "MSGCorep-identity-direct-sum-dimensions-1-through-10"
]
