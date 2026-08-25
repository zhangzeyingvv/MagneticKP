If[! ValueQ[$MagneticKPTestHelpersLoaded],
  Get[FileNameJoin[{
    ParentDirectory[DirectoryName[$InputFileName]],
    "Support",
    "TestHelpers.wl"
  }]]
];
LoadMagneticKPForTests[];

VerificationTest[
  PrivateHookAvailableQ["MomentumTransformMatrix"],
  True,
  TestID -> "MomentumTransform-private-hook-exists"
]

legacyFDPath = FileNameJoin[{$MagneticKPLegacyRoot, "fd.mx"}];
legacyFD = If[FileExistsQ[legacyFDPath], Import[legacyFDPath], $Failed];

VerificationTest[
  AssociationQ[legacyFD] && Length[legacyFD] == 1368,
  True,
  TestID -> "MomentumTransform-legacy-fd-oracle-has-1368-entries"
]

fdMismatches = If[AssociationQ[legacyFD],
  Reap[
    KeyValueMap[
      Function[{key, expected},
        computed = CallMomentumTransform[key[[1]], key[[2]]];
        If[
          computed === $Failed ||
            ! MatrixExactEqualQ[computed, expected],
          Sow[key]
        ]
      ],
      legacyFD
    ]
  ][[2]],
  {{{"legacy-fd-unavailable"}}}
];
fdMismatches = If[fdMismatches === {}, {}, First[fdMismatches]];

VerificationTest[
  fdMismatches,
  {},
  TestID -> "MomentumTransform-all-1368-fd-entries-exact"
]

countCases = Flatten[Table[
  {
    dim,
    order,
    Length[KMonomialBasis[order, dim]],
    Binomial[order + dim - 1, dim - 1]
  },
  {dim, 1, 3}, {order, 0, 12}
], 1];

VerificationTest[
  And @@ (#[[3]] == #[[4]] & /@ countCases),
  True,
  TestID -> "MomentumTransform-monomial-count-d1-d3-order0-order12"
]

definitionMatrices = {
  {{1}},
  {{0, -1}, {1, 0}},
  {{0, -1, 0}, {1, 0, 0}, {0, 0, 1}},
  {{1, 1, 0}, {0, 1, 1}, {0, 0, 1}},
  {{0, 0, 1}, {1, 0, 0}, {0, 1, 0}}
};

VerificationTest[
  And @@ Flatten @ Table[
    dim = Length[a];
    vars = KVariables[dim];
    basis = KMonomialBasis[order, dim];
    image = LinearKImage[a];
    f = CallMomentumTransform[order, image];
    f =!= $Failed &&
      TrueQ[FullSimplify[
        basis == (basis /. Thread[vars -> image]) . f
      ]],
    {a, definitionMatrices}, {order, 0, 5}
  ],
  True,
  TestID -> "MomentumTransform-variable-substitution-definition"
]

a3 = {{0, -1, 0}, {1, 0, 0}, {0, 0, 1}};
b3 = {{0, 0, 1}, {1, 0, 0}, {0, 1, 0}};

VerificationTest[
  And @@ Table[
    MatrixExactEqualQ[
      CallMomentumTransform[order, LinearKImage[a3 . b3]],
      CallMomentumTransform[order, LinearKImage[a3]] .
        CallMomentumTransform[order, LinearKImage[b3]]
    ],
    {order, 0, 8}
  ],
  True,
  TestID -> "MomentumTransform-homomorphism"
]

VerificationTest[
  And @@ Table[
    f = CallMomentumTransform[
      order,
      KVariables[3]
    ];
    MatrixExactEqualQ[
      f,
      IdentityMatrix[Binomial[order + 2, 2]]
    ],
    {order, 0, 12}
  ],
  True,
  TestID -> "MomentumTransform-identity-orders-0-through-12"
]

VerificationTest[
  And @@ Table[
    f = CallMomentumTransform[order, LinearKImage[a3]];
    TrueQ[FullSimplify[
      Det[f] ==
        Det[a3]^(-Binomial[order + 2, 3])
    ]],
    {order, 0, 8}
  ],
  True,
  TestID -> "MomentumTransform-determinant-formula"
]

VerificationTest[
  ! FileExistsQ[
    FileNameJoin[{$MagneticKPPackageRoot, "fd.mx"}]
  ],
  True,
  TestID -> "MomentumTransform-production-fd-mx-removed"
]
