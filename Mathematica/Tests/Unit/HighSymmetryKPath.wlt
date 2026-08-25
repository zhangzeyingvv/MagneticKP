If[! ValueQ[$MagneticKPTestHelpersLoaded],
  Get[FileNameJoin[{
    ParentDirectory[DirectoryName[$InputFileName]],
    "Support",
    "TestHelpers.wl"
  }]]
];
LoadMagneticKPForTests[];

kPathC2Input = <|
  "Unitary" -> <|
    "C2z" -> {
      IdentityMatrix[1],
      {-MagneticKP`kx, -MagneticKP`ky, MagneticKP`kz}
    }
  |>,
  "KPointLabel" -> "K",
  "CoordinateSystem" -> "ReciprocalFractional"
|>;
kPathC2Data = MagneticKP`highSymmetryKPath[
  kPathC2Input,
  "ReturnData" -> True
];

VerificationTest[
  KeyTake[
    kPathC2Data,
    {
      "Directions",
      "LittleCoGroupOrder",
      "EnumeratedSubgroupCount",
      "CenterLabel",
      "CoordinateSystem"
    }
  ],
  <|
    "Directions" -> {{0, 0, 1}},
    "LittleCoGroupOrder" -> 2,
    "EnumeratedSubgroupCount" -> 2,
    "CenterLabel" -> "K",
    "CoordinateSystem" -> "ReciprocalFractional"
  |>,
  TestID -> "HighSymmetryKPath-C2-largest-line-stabilizer"
]

VerificationTest[
  {
    kPathC2Data["Path"],
    PrivateHook["validBandPathQ"][kPathC2Data["Path"]]
  },
  {
    {
      {{{0, 0, -1/10}, {0, 0, 0}}, {"L1-", "K"}},
      {{{0, 0, 0}, {0, 0, 1/10}}, {"K", "L1+"}}
    },
    True
  },
  TestID -> "HighSymmetryKPath-C2-ordinary-band-path"
]

VerificationTest[
  MagneticKP`highSymmetryKPath[kPathC2Input],
  kPathC2Data["Path"],
  TestID -> "HighSymmetryKPath-default-returns-path-only"
]

kPathD4hInput = <|
  "Unitary" -> <|
    "C4z" -> {
      IdentityMatrix[1],
      {-MagneticKP`ky, MagneticKP`kx, MagneticKP`kz}
    },
    "mxz" -> {
      IdentityMatrix[1],
      {MagneticKP`kx, -MagneticKP`ky, MagneticKP`kz}
    },
    "mxy" -> {
      IdentityMatrix[1],
      {MagneticKP`kx, MagneticKP`ky, -MagneticKP`kz}
    }
  |>,
  "KPointLabel" -> "K"
|>;
kPathD4hRepresentatives = MagneticKP`highSymmetryKPath[
  kPathD4hInput,
  "ReturnData" -> True
];
kPathD4hAll = MagneticKP`highSymmetryKPath[
  kPathD4hInput,
  "EquivalentDirections" -> "All",
  "ReturnData" -> True
];

VerificationTest[
  {
    kPathD4hRepresentatives["LittleCoGroupOrder"],
    Length[kPathD4hRepresentatives["Directions"]],
    Lookup[kPathD4hRepresentatives["LineData"], "SubgroupOrder"],
    Length[kPathD4hAll["Directions"]]
  },
  {16, 3, {8, 4, 4}, 5},
  TestID -> "HighSymmetryKPath-D4h-equivalent-line-representatives"
]

VerificationTest[
  Length @ MagneticKP`highSymmetryKPath[
    kPathD4hInput,
    "IncludeNegativeDirections" -> False
  ],
  3,
  TestID -> "HighSymmetryKPath-one-arm-per-direction"
]

VerificationTest[
  Quiet @ MagneticKP`highSymmetryKPath[
    kPathC2Input,
    "Radius" -> 0
  ],
  $Failed,
  TestID -> "HighSymmetryKPath-invalid-option-fails"
]

VerificationTest[
  Quiet @ MagneticKP`highSymmetryKPath[<||>],
  $Failed,
  TestID -> "HighSymmetryKPath-invalid-input-fails"
]
