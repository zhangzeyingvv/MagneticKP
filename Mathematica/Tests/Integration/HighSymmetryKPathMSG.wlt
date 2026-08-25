If[! ValueQ[$MagneticKPTestHelpersLoaded],
  Get[FileNameJoin[{
    ParentDirectory[DirectoryName[$InputFileName]],
    "Support",
    "TestHelpers.wl"
  }]]
];
LoadMagneticKPForTests[];
LoadMSGCorepForTests[];

kPathMSGGeneratorInput = MagneticKP`interfaceRep[
  {156, 50},
  "K",
  {2, 3},
  "CartesianCoordinates" -> False,
  "CalculateGenerators" -> True
];
kPathMSGElementInput = MagneticKP`interfaceRep[
  {156, 50},
  "K",
  {2, 3},
  "CartesianCoordinates" -> False,
  "CalculateGenerators" -> False
];
kPathMSGGeneratorData = MagneticKP`highSymmetryKPath[
  kPathMSGGeneratorInput,
  "ReturnData" -> True
];
kPathMSGElementData = MagneticKP`highSymmetryKPath[
  kPathMSGElementInput,
  "ReturnData" -> True
];

VerificationTest[
  {
    kPathMSGGeneratorData["Directions"],
    kPathMSGGeneratorData["LittleCoGroupOrder"],
    kPathMSGGeneratorData["CenterLabel"],
    kPathMSGGeneratorData["CoordinateSystem"],
    PrivateHook["validBandPathQ"][kPathMSGGeneratorData["Path"]]
  },
  {
    {{0, 0, 1}, {1, 1, 0}},
    6,
    "K",
    "ReciprocalFractional",
    True
  },
  TestID -> "HighSymmetryKPath-MSG156.50-K-lines"
]

VerificationTest[
  {
    kPathMSGGeneratorData["Path"],
    kPathMSGGeneratorData["Directions"],
    kPathMSGGeneratorData["LittleCoGroupOrder"]
  },
  {
    kPathMSGElementData["Path"],
    kPathMSGElementData["Directions"],
    kPathMSGElementData["LittleCoGroupOrder"]
  },
  TestID -> "HighSymmetryKPath-generators-elements-parity"
]

VerificationTest[
  MatrixQ @ MagneticKP`kpHam[1, kPathMSGGeneratorInput]["ham"],
  True,
  TestID -> "HighSymmetryKPath-interface-metadata-preserves-kpHam"
]

kPathDegenerateInput = MagneticKP`interfaceRep[
  {164, 90},
  {0, 0, 0},
  {5, 6},
  "CartesianCoordinates" -> False
];
kPathDegenerateData = MagneticKP`highSymmetryKPath[
  kPathDegenerateInput,
  "ReturnData" -> True
];
kPathDegenerateHamiltonian =
  MagneticKP`kpHam[1, kPathDegenerateInput]["ham"];
kPathDegenerateParameters = Complement[
  Variables[kPathDegenerateHamiltonian],
  {MagneticKP`kx, MagneticKP`ky, MagneticKP`kz}
];
kPathDegenerateRules = Thread[
  kPathDegenerateParameters ->
    Prime[Range[Length[kPathDegenerateParameters]] + 10]/100
];
kPathDiscriminant[direction_List] := Together @ Discriminant[
  CharacteristicPolynomial[
    kPathDegenerateHamiltonian /.
      Thread[
        {MagneticKP`kx, MagneticKP`ky, MagneticKP`kz} ->
          direction/10
      ] /.
      kPathDegenerateRules,
    Global`kPathEnergy
  ],
  Global`kPathEnergy
];

VerificationTest[
  {
    Dimensions[kPathDegenerateHamiltonian],
    kPathDegenerateData["Directions"],
    Lookup[kPathDegenerateData["LineData"], "SubgroupOrder"],
    Lookup[kPathDegenerateData["LineData"], "MomentumActionOrder"],
    PossibleZeroQ[kPathDiscriminant[#]] & /@
      kPathDegenerateData["Directions"],
    PossibleZeroQ[kPathDiscriminant[{1, 2, 3}]]
  },
  {
    {4, 4},
    {{0, 0, 1}, {1, 1, 0}},
    {12, 4},
    {6, 2},
    {True, True},
    False
  },
  TestID -> "HighSymmetryKPath-MSG164.90-line-degeneracy"
]

kPathOriginalExampleInput = MagneticKP`interfaceRep[
  {191, 234},
  "K",
  {5, 6},
  "CartesianCoordinates" -> False
];
kPathOriginalExampleData = MagneticKP`highSymmetryKPath[
  kPathOriginalExampleInput,
  "ReturnData" -> True
];
kPathOriginalExampleHamiltonian =
  MagneticKP`kpHam[1, kPathOriginalExampleInput]["ham"];
kPathOriginalExampleParameters = Complement[
  Variables[kPathOriginalExampleHamiltonian],
  {MagneticKP`kx, MagneticKP`ky, MagneticKP`kz}
];
kPathOriginalExampleRules = Thread[
  kPathOriginalExampleParameters ->
    Prime[Range[Length[kPathOriginalExampleParameters]] + 10]/100
];
kPathOriginalExampleDiscriminant[direction_List] := Together @ Discriminant[
  CharacteristicPolynomial[
    kPathOriginalExampleHamiltonian /.
      Thread[
        {MagneticKP`kx, MagneticKP`ky, MagneticKP`kz} ->
          direction/10
      ] /.
      kPathOriginalExampleRules,
    Global`kPathOriginalEnergy
  ],
  Global`kPathOriginalEnergy
];

VerificationTest[
  {
    Dimensions[kPathOriginalExampleHamiltonian],
    kPathOriginalExampleData["LittleCoGroupOrder"],
    kPathOriginalExampleData["LittleCoGroupActionOrder"],
    kPathOriginalExampleData["Directions"],
    Lookup[kPathOriginalExampleData["LineData"], "SubgroupOrder"],
    Lookup[kPathOriginalExampleData["LineData"], "MomentumActionOrder"],
    PossibleZeroQ[kPathOriginalExampleDiscriminant[#]] & /@
      kPathOriginalExampleData["Directions"],
    PossibleZeroQ[kPathOriginalExampleDiscriminant[{1, 2, 3}]]
  },
  {
    {4, 4},
    24,
    12,
    {{0, 0, 1}, {1, 1, 0}},
    {12, 8},
    {6, 4},
    {True, False},
    False
  },
  TestID -> "HighSymmetryKPath-original-MSG191.234-K-line-degeneracy"
]
