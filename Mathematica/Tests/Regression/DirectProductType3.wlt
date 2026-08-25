If[! ValueQ[$MagneticKPTestHelpersLoaded],
  Get[FileNameJoin[{
    ParentDirectory[DirectoryName[$InputFileName]],
    "Support",
    "TestHelpers.wl"
  }]]
];
LoadMagneticKPForTests[];

type3Input = <|
  "Unitary" -> <|
    "E" -> {
      IdentityMatrix[2],
      {MagneticKP`kx, MagneticKP`ky, MagneticKP`kz}
    }
  |>,
  "AntisymmetryUnitaryTest" -> <|
    "A" -> {
      PauliMatrix[3],
      {MagneticKP`kx, MagneticKP`ky, MagneticKP`kz}
    }
  |>
|>;

type3Iterative = MagneticKP`IterativelySimplify[0, type3Input];
type3Direct = MagneticKP`DirectProductDecomposition[0, type3Input];

VerificationTest[
  type3Iterative["NumberOfParameters"],
  2,
  TestID -> "DirectProduct-type3-independent-reference-parameter-count"
]

VerificationTest[
  type3Direct["NumberOfParameters"],
  2,
  TestID -> "DirectProduct-type3-does-not-reuse-stale-Gg"
]

VerificationTest[
  HamiltonianSubspaceEqualQ[
    type3Iterative["ham"],
    type3Direct["ham"]
  ],
  True,
  TestID -> "DirectProduct-type3-subspace-equals-iterative"
]
