If[! ValueQ[$MagneticKPTestHelpersLoaded],
  Get[FileNameJoin[{
    ParentDirectory[DirectoryName[$InputFileName]],
    "Support",
    "TestHelpers.wl"
  }]]
];
LoadMagneticKPForTests[];

VerificationTest[
  PrivateHookAvailableQ["OperatorBasisData"],
  True,
  TestID -> "OperatorBasis-private-hook-exists"
]

basisScan = Association @ Table[
  data = CallOperatorBasisData[n];
  basis = If[AssociationQ[data], Lookup[data, "Basis", {}], {}];
  gram = If[AssociationQ[data], Lookup[data, "Gram", {}], {}];
  ToString[n] -> <|
    "AssociationQ" -> AssociationQ[data],
    "Count" -> Length[basis],
    "DimensionsQ" -> And @@ (Dimensions[#] === {n, n} & /@ basis),
    "HermitianQ" ->
      And @@ (TrueQ[# === ConjugateTranspose[#]] & /@ basis),
    "DistinctQ" -> DuplicateFreeQ[basis],
    "VectorRank" -> If[basis === {}, 0,
      MatrixRank[Flatten /@ basis]
    ],
    "GramRank" -> If[MatrixQ[gram], MatrixRank[gram], 0],
    "GramHermitianQ" -> If[MatrixQ[gram],
      TrueQ[gram === ConjugateTranspose[gram]],
      False
    ]
  |>,
  {n, 1, 10}
];

VerificationTest[
  And @@ KeyValueMap[
    Function[{nString, record},
      With[{n = ToExpression[nString]},
        record["AssociationQ"] &&
        record["Count"] == n^2 &&
        record["DimensionsQ"] &&
        record["HermitianQ"] &&
        record["DistinctQ"] &&
        record["VectorRank"] == n^2 &&
        record["GramRank"] == n^2 &&
        record["GramHermitianQ"]
      ]
    ],
    basisScan
  ],
  True,
  TestID -> "OperatorBasis-n-1-through-10-structural-invariants"
]

VerificationTest[
  Lookup[basisScan, #]["Count"] & /@
    {"5", "7", "9", "10"},
  {25, 49, 81, 100},
  TestID -> "OperatorBasis-generic-special-dimensions"
]

VerificationTest[
  And @@ Table[
    data = CallOperatorBasisData[n];
    basis = data["Basis"];
    gram = data["Gram"];
    norms = Tr[ConjugateTranspose[#] . #] & /@ basis;
    MatrixExactEqualQ[
      gram,
      DiagonalMatrix[norms]
    ],
    {n, 1, 10}
  ],
  True,
  TestID -> "OperatorBasis-reported-Gram-is-Hilbert-Schmidt-Gram"
]

VerificationTest[
  And @@ (
    Quiet @ Check[CallOperatorBasisData[#] === $Failed, True] & /@
      {0, -1}
  ),
  True,
  TestID -> "OperatorBasis-rejects-nonpositive-dimensions"
]
