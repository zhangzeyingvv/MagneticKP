If[! ValueQ[$MagneticKPTestHelpersLoaded],
  Get[FileNameJoin[{
    ParentDirectory[DirectoryName[$InputFileName]],
    "Support",
    "TestHelpers.wl"
  }]]
];
LoadMagneticKPForTests[];

continuousTheta = Unique["theta"];
continuousVariables = {
  MagneticKP`kx,
  MagneticKP`ky,
  MagneticKP`kz
};
continuousRepresentation = DiagonalMatrix[{
  Exp[-I continuousTheta/2],
  Exp[I continuousTheta/2]
}];
continuousOnlyInput = <|
  "Unitary" -> <|
    "CInfinity" -> {
      continuousRepresentation,
      continuousVariables,
      continuousTheta
    }
  |>
|>;

continuousOnlyResult = MagneticKP`kpHam[
  1,
  continuousOnlyInput,
  "ReturnDiagnostics" -> True
];

VerificationTest[
  AssociationQ[continuousOnlyResult] &&
    continuousOnlyResult["NumberOfParameters"] == 8 &&
    Lookup[
      continuousOnlyResult["Diagnostics"],
      "InvariantDimension"
    ] === {2, 6},
  True,
  TestID -> "CInfinity-continuous-only-order-zero-and-one-count"
]

VerificationTest[
  Lookup[
    continuousOnlyResult["Diagnostics"],
    {
      "ContinuousSymmetry",
      "OriginalOperatorDimension",
      "ReducedOperatorDimension",
      "Solver"
    }
  ] === {
    {"CInfinity", 4, 2, "None"},
    {"CInfinity", 4, 2, "None"}
  },
  True,
  TestID -> "CInfinity-diagnostics-report-operator-reduction"
]

continuousOrderZero = MagneticKP`kpHam[0, continuousOnlyInput];
continuousReferenceHamiltonian =
  continuousIdentityCoefficient IdentityMatrix[2] +
    continuousZCoefficient PauliMatrix[3];

VerificationTest[
  AssociationQ[continuousOrderZero] &&
    continuousOrderZero["NumberOfParameters"] == 2 &&
    HamiltonianSubspaceEqualQ[
      continuousOrderZero["ham"],
      continuousReferenceHamiltonian
    ],
  True,
  TestID -> "CInfinity-spin-half-commutant-is-I-and-sigma-z"
]

continuousNamedResult = MagneticKP`kpHam[3, continuousOnlyInput];
continuousNamedParameters = Complement[
  Variables[continuousNamedResult["ham"]],
  continuousVariables
];
continuousNamedPairs = SortBy[
  {SymbolName[#[[1]]], #[[2]]} & /@ continuousNamedParameters,
  {First[#], Last[#]} &
];
continuousNamedFamilySymbols = DeleteDuplicates @ Cases[
  continuousNamedParameters,
  Subscript[symbol_Symbol, ___] :> HoldComplete[symbol],
  Infinity
];
continuousExpectedNamedPairs = SortBy[
  Join[
    Table[{"e", index}, {index, 2}],
    Table[{"t", index}, {index, 6}],
    Table[{"r", index}, {index, 12}],
    Table[{"s", index}, {index, 20}]
  ],
  {First[#], Last[#]} &
];

VerificationTest[
  continuousNamedResult["NumberOfParameters"] == 40 &&
    continuousNamedPairs === continuousExpectedNamedPairs &&
    Sort[continuousNamedFamilySymbols] === Sort[{
      HoldComplete[MagneticKP`e],
      HoldComplete[MagneticKP`t],
      HoldComplete[MagneticKP`r],
      HoldComplete[MagneticKP`s]
    }] &&
    StringFreeQ[
      ToString[continuousNamedResult["ham"], InputForm],
      "MagneticKP`Private`"
    ],
  True,
  TestID -> "low-orders-use-public-e-t-r-s-one-based-coefficient-families"
]

(*
  A single rotation by theta0 = ArcCos[3/5] generates a dense cyclic
  subgroup of CInfinity: theta0/Pi is irrational by Niven's theorem.
  Its exact rational-complex matrix avoids transcendental simplification.
  The spin-one representation also supplies a negative control at Pi,
  where the two outer weights become accidentally degenerate.
*)
continuousDenseTheta = Unique["denseTheta"];
continuousDenseRepresentation = DiagonalMatrix[{
  Exp[-I continuousDenseTheta],
  1,
  Exp[I continuousDenseTheta]
}];
continuousDenseInput = <|
  "Unitary" -> <|
    "CInfinity" -> {
      continuousDenseRepresentation,
      continuousVariables,
      continuousDenseTheta
    }
  |>
|>;
continuousDenseCyclicInput = <|
  "Unitary" -> <|
    "DenseCyclicGenerator" -> {
      DiagonalMatrix[{(3 - 4 I)/5, 1, (3 + 4 I)/5}],
      continuousVariables
    }
  |>
|>;
continuousRationalPiInput = <|
  "Unitary" -> <|
    "OrderTwoRotation" -> {
      DiagonalMatrix[{-1, 1, -1}],
      continuousVariables
    }
  |>
|>;
continuousDenseComparisons = Table[
  Module[{continuousResult, denseResult, rationalResult},
    continuousResult = MagneticKP`kpHam[
      order,
      continuousDenseInput,
      "Method" -> "IterativeSimplification"
    ];
    denseResult = MagneticKP`kpHam[
      order,
      continuousDenseCyclicInput,
      "Method" -> "IterativeSimplification"
    ];
    rationalResult = MagneticKP`kpHam[
      order,
      continuousRationalPiInput,
      "Method" -> "IterativeSimplification"
    ];
    {
      continuousResult["NumberOfParameters"],
      denseResult["NumberOfParameters"],
      HamiltonianSubspaceEqualQ[
        continuousResult["ham"],
        denseResult["ham"]
      ],
      rationalResult["NumberOfParameters"]
    }
  ],
  {order, 0, 2}
];

VerificationTest[
  continuousDenseComparisons,
  {
    {3, 3, True, 5},
    {12, 12, True, 20},
    {30, 30, True, 50}
  },
  TestID -> "CInfinity-equals-dense-irrational-cyclic-subgroup"
]

continuousFlipInput = <|
  "Unitary" -> <|
    "Cinfty" -> {
      continuousRepresentation,
      continuousVariables,
      continuousTheta
    },
    "SpinFlip" -> {
      PauliMatrix[1],
      continuousVariables
    }
  |>
|>;
continuousFlipIterative = MagneticKP`kpHam[
  0,
  continuousFlipInput,
  "Method" -> "IterativeSimplification",
  "ReturnDiagnostics" -> True
];
continuousFlipReynolds = Quiet @ MagneticKP`kpHam[
  0,
  continuousFlipInput,
  "Method" -> "Reynolds",
  "GroupInput" -> "Generators",
  "ReturnDiagnostics" -> True
];
continuousIdentityHamiltonian =
  continuousScalarCoefficient IdentityMatrix[2];

VerificationTest[
  And @@ (AssociationQ /@
      {continuousFlipIterative, continuousFlipReynolds}) &&
    continuousFlipIterative["NumberOfParameters"] == 1 &&
    continuousFlipReynolds["NumberOfParameters"] == 1 &&
    HamiltonianSubspaceEqualQ[
      continuousFlipIterative["ham"],
      continuousIdentityHamiltonian
    ] &&
    HamiltonianSubspaceEqualQ[
      continuousFlipIterative["ham"],
      continuousFlipReynolds["ham"]
    ],
  True,
  TestID -> "CInfinity-then-discrete-normalizer-iterative-Reynolds-parity"
]

continuousBadNormalizerInput = <|
  "Unitary" -> <|
    "CInfinity" -> {
      continuousRepresentation,
      continuousVariables,
      continuousTheta
    },
    "Bad" -> {
      {{1, 1}, {1, -1}}/Sqrt[2],
      continuousVariables
    }
  |>
|>;

VerificationTest[
  Quiet[MagneticKP`kpHam[0, continuousBadNormalizerInput]],
  $Failed,
  SameTest -> SameQ,
  TestID -> "CInfinity-rejects-discrete-nonnormalizer"
]

continuousBadMomentumInput = <|
  "Unitary" -> <|
    "CInfinity" -> {
      continuousRepresentation,
      {-MagneticKP`kx, MagneticKP`ky, MagneticKP`kz},
      continuousTheta
    }
  |>
|>;

VerificationTest[
  Quiet[MagneticKP`kpHam[0, continuousBadMomentumInput]],
  $Failed,
  SameTest -> SameQ,
  TestID -> "CInfinity-requires-trivial-momentum-action"
]

continuousExtraSymbol = Unique["m"];
continuousExtraSymbolInput = <|
  "Unitary" -> <|
    "CInfinity" -> {
      DiagonalMatrix[{
        Exp[-I continuousExtraSymbol continuousTheta],
        Exp[I continuousExtraSymbol continuousTheta]
      }],
      continuousVariables,
      continuousTheta
    }
  |>
|>;

VerificationTest[
  Quiet[MagneticKP`kpHam[0, continuousExtraSymbolInput]],
  $Failed,
  SameTest -> SameQ,
  TestID -> "CInfinity-representation-allows-only-declared-symbol"
]

continuousWrongCategoryInput = <|
  "Antiunitary" -> <|
    "CInfinity" -> {
      continuousRepresentation,
      continuousVariables,
      continuousTheta
    }
  |>
|>;

VerificationTest[
  Quiet[MagneticKP`kpHam[0, continuousWrongCategoryInput]],
  $Failed,
  SameTest -> SameQ,
  TestID -> "CInfinity-is-allowed-only-in-Unitary"
]

continuousNonRepresentationInput = <|
  "Unitary" -> <|
    "CInfinity" -> {
      DiagonalMatrix[{
        Exp[-I continuousTheta^2/2],
        Exp[I continuousTheta^2/2]
      }],
      continuousVariables,
      continuousTheta
    }
  |>
|>;

VerificationTest[
  Quiet[MagneticKP`kpHam[0, continuousNonRepresentationInput]],
  $Failed,
  SameTest -> SameQ,
  TestID -> "CInfinity-rejects-theta-dependent-logarithmic-generator"
]
