If[! ValueQ[$MagneticKPTestHelpersLoaded],
  Get[FileNameJoin[{
    ParentDirectory[DirectoryName[$InputFileName]],
    "Support",
    "TestHelpers.wl"
  }]]
];
LoadMagneticKPForTests[];

VerificationTest[
  PrivateHookAvailableQ["AdjointRepresentationMatrix"],
  True,
  TestID -> "AdjointRepresentation-private-hook-exists"
]

exactRepresentation[n_] := Which[
  n == 1,
    {{1}},
  n >= 2,
    SparseArray[
      Join[
        Table[{i, n + 1 - i} -> (-1)^(i - 1), {i, n}]
      ],
      {n, n}
    ]
];

adjointOldMismatches = Reap[
  Do[
    data = CallOperatorBasisData[n];
    d = exactRepresentation[n];
    Do[
      newJ = CallAdjointRepresentation[d, data, anti];
      oldJ = LegacyAdjointMatrix[d, data["Basis"], anti];
      If[
        newJ === $Failed ||
          ! MatrixExactEqualQ[newJ, oldJ],
        Sow[{n, anti}]
      ],
      {anti, {False, True}}
    ],
    {n, 1, 10}
  ]
][[2]];
adjointOldMismatches =
  If[adjointOldMismatches === {}, {}, First[adjointOldMismatches]];

VerificationTest[
  adjointOldMismatches,
  {},
  TestID -> "AdjointRepresentation-old-double-loop-parity-n1-n10"
]

VerificationTest[
  And @@ Table[
    data = CallOperatorBasisData[n];
    d = exactRepresentation[n];
    j = CallAdjointRepresentation[d, data, False];
    imageColumns = Transpose[
      Flatten /@
        (d . # . Inverse[d] & /@ data["Basis"])
    ];
    MatrixExactEqualQ[
      imageColumns,
      data["Columns"] . j
    ],
    {n, 1, 10}
  ],
  True,
  TestID -> "AdjointRepresentation-unitary-column-reconstruction"
]

VerificationTest[
  And @@ Table[
    data = CallOperatorBasisData[n];
    d = exactRepresentation[n];
    j = CallAdjointRepresentation[d, data, True];
    imageColumns = Transpose[
      Flatten /@
        (d . Conjugate[#] . Inverse[d] & /@
          data["Basis"])
    ];
    MatrixExactEqualQ[
      imageColumns,
      data["Columns"] . j
    ],
    {n, 1, 10}
  ],
  True,
  TestID -> "AdjointRepresentation-antiunitary-column-reconstruction"
]

VerificationTest[
  And @@ Table[
    data = CallOperatorBasisData[n];
    identityJ = CallAdjointRepresentation[
      IdentityMatrix[n],
      data,
      False
    ];
    MatrixExactEqualQ[identityJ, IdentityMatrix[n^2]],
    {n, 1, 10}
  ],
  True,
  TestID -> "AdjointRepresentation-identity"
]

d1 = {{1, 0}, {0, I}};
d2 = 1/Sqrt[2] {{1, 1}, {-1, 1}};
data2 = CallOperatorBasisData[2];

VerificationTest[
  MatrixExactEqualQ[
    CallAdjointRepresentation[d1 . d2, data2, False],
    CallAdjointRepresentation[d1, data2, False] .
      CallAdjointRepresentation[d2, data2, False]
  ],
  True,
  TestID -> "AdjointRepresentation-unitary-homomorphism"
]

VerificationTest[
  MatrixExactEqualQ[
    CallAdjointRepresentation[d1, data2, True] .
      CallAdjointRepresentation[d2, data2, True],
    CallAdjointRepresentation[d1 . Conjugate[d2], data2, False]
  ],
  True,
  TestID -> "AdjointRepresentation-two-antiunitary-composition"
]

VerificationTest[
  MatrixExactEqualQ[
    CallAdjointRepresentation[-d1, data2, False],
    CallAdjointRepresentation[d1, data2, False]
  ],
  True,
  TestID -> "AdjointRepresentation-projective-sign-invariance"
]

VerificationTest[
  Module[{j = CallAdjointRepresentation[d2, data2, False],
    gram = data2["Gram"]},
    MatrixExactEqualQ[
      ConjugateTranspose[j] . gram . j,
      gram
    ]
  ],
  True,
  TestID -> "AdjointRepresentation-Gram-metric-invariance"
]

nonorthogonalBasis = {
  PauliMatrix[0] + PauliMatrix[3],
  PauliMatrix[1],
  PauliMatrix[2],
  PauliMatrix[3]
};
nonorthogonalData = CallOperatorBasisData[nonorthogonalBasis];

VerificationTest[
  MatrixExactEqualQ[
    CallAdjointRepresentation[d2, nonorthogonalData, False],
    ReferenceAdjointMatrix[d2, nonorthogonalBasis, False]
  ],
  True,
  TestID -> "AdjointRepresentation-full-Gram-nonorthogonal-basis"
]

singularData = {
  {PauliMatrix[0], PauliMatrix[0], PauliMatrix[1], PauliMatrix[2]}
};

VerificationTest[
  Module[{value = Quiet @ Check[
      CallOperatorBasisData[First[singularData]],
      $Failed
    ]},
    value === $Failed || FailureQ[value]
  ],
  True,
  TestID -> "AdjointRepresentation-rejects-singular-Gram"
]
