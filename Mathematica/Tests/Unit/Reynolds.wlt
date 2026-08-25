If[! ValueQ[$MagneticKPTestHelpersLoaded],
  Get[FileNameJoin[{
    ParentDirectory[DirectoryName[$InputFileName]],
    "Support",
    "TestHelpers.wl"
  }]]
];
LoadMagneticKPForTests[];

VerificationTest[
  PrivateHookAvailableQ["ReynoldsInvariantBasis"],
  True,
  TestID -> "Reynolds-private-hook-exists"
]

c2 = {
  IdentityMatrix[2],
  DiagonalMatrix[{1, -1}]
};
c3generator = {
  {0, -1},
  {1, -1}
};
c3 = NestList[c3generator . # &, IdentityMatrix[2], 2];
s3 = {
  IdentityMatrix[3],
  {{0, 1, 0}, {1, 0, 0}, {0, 0, 1}},
  {{0, 0, 1}, {0, 1, 0}, {1, 0, 0}},
  {{1, 0, 0}, {0, 0, 1}, {0, 1, 0}},
  {{0, 1, 0}, {0, 0, 1}, {1, 0, 0}},
  {{0, 0, 1}, {1, 0, 0}, {0, 1, 0}}
};

reynoldsCases = <|"C2" -> c2, "C3" -> c3, "S3" -> s3|>;

VerificationTest[
  And @@ KeyValueMap[
    Function[{name, reps},
      p = Total[reps]/Length[reps];
      basis = CallReynoldsInvariantBasis[reps];
      MatrixExactEqualQ[p . p, p] &&
        MatrixRank[p] == Tr[p] &&
        And @@ (MatrixExactEqualQ[# . p, p] & /@ reps) &&
        SubspaceEqualQ[basis, ImageBasis[p]] &&
        SubspaceEqualQ[basis, IterativeInvariantBasis[reps]] &&
        SubspaceEqualQ[basis, StackedInvariantBasis[reps]]
    ],
    reynoldsCases
  ],
  True,
  TestID -> "Reynolds-idempotence-absorption-rank-trace-subspace"
]

twistedCharacters = {1, -1};
twistedProjector =
  Total[MapThread[Conjugate[#1] #2 &, {twistedCharacters, c2}]]/2;
twistedBasis =
  CallReynoldsInvariantBasis[c2, twistedCharacters];

VerificationTest[
  MatrixExactEqualQ[
    twistedProjector . twistedProjector,
    twistedProjector
  ] &&
  And @@ MapThread[
    MatrixExactEqualQ[#1 . twistedProjector, #2 twistedProjector] &,
    {c2, twistedCharacters}
  ] &&
  SubspaceEqualQ[twistedBasis, ImageBasis[twistedProjector]],
  True,
  TestID -> "Reynolds-twisted-antisymmetry-character"
]

VerificationTest[
  Module[
    {
      r1 = {{1, 0}, {0, -1}},
      r2 = 1/2 {{-1, Sqrt[3]}, {Sqrt[3], 1}},
      naive
    },
    naive =
      ((IdentityMatrix[2] + r1)/2) .
      ((IdentityMatrix[2] + r2)/2);
    ! MatrixExactEqualQ[naive . naive, naive]
  ],
  True,
  TestID -> "Reynolds-noncommuting-generator-averages-not-a-projector"
]

automaticLargeInput = <|
  "Unitary" -> <|
    "C5" -> {
      RotateRight[IdentityMatrix[5]],
      {MagneticKP`kx, MagneticKP`ky, MagneticKP`kz}
    }
  |>
|>;
automaticLargeResult = MagneticKP`kpHam[
  {2},
  automaticLargeInput,
  "Method" -> Automatic,
  "ReturnDiagnostics" -> True
];
iterativeLargeResult = MagneticKP`kpHam[
  {2},
  automaticLargeInput,
  "Method" -> "IterativeSimplification",
  "ReturnDiagnostics" -> True
];

VerificationTest[
  Module[{automaticDiagnostics, iterativeDiagnostics},
    automaticDiagnostics = First[
      automaticLargeResult["Diagnostics"]
    ];
    iterativeDiagnostics = First[
      iterativeLargeResult["Diagnostics"]
    ];
    automaticDiagnostics["AmbientDimension"] > 128 &&
      automaticDiagnostics["Solver"] === "Iterative" &&
      iterativeDiagnostics["Solver"] === "Iterative" &&
      automaticLargeResult["NumberOfParameters"] ==
        iterativeLargeResult["NumberOfParameters"] &&
      HamiltonianSubspaceEqualQ[
        automaticLargeResult["ham"],
        iterativeLargeResult["ham"]
      ]
  ],
  True,
  TestID -> "Automatic-large-ambient-conservatively-selects-iterative"
]
