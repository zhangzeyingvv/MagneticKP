If[! ValueQ[$MagneticKPTestHelpersLoaded],
  Get[FileNameJoin[{
    ParentDirectory[DirectoryName[$InputFileName]],
    "Support",
    "TestHelpers.wl"
  }]]
];
LoadMagneticKPForTests[];

legacyOracle = GetLegacyOracleData[];

VerificationTest[
  AssociationQ[legacyOracle] &&
    ! KeyExistsQ[legacyOracle, "Failure"],
  True,
  TestID -> "LegacyOracle-independent-kernel-completed"
]

VerificationTest[
  Lookup[legacyOracle, "LegacyRoot", Missing[]],
  $MagneticKPLegacyRoot,
  TestID -> "LegacyOracle-used-frozen-MagneticKP-Mathmatica-root"
]

VerificationTest[
  Sort @ Keys[legacyOracle["Cases"]],
  Sort @ {
    "MSG13.65-B-2",
    "MSG1.2-Gamma-1",
    "MSG1.2-Gamma-2",
    "MSG79.27-Gamma-1",
    "MSG97.155-P-1+2",
    "MSG80.32-Gamma-1",
    "MSG75.4-R-1+2",
    "MSG164.90-Gamma-5+6+7+8"
  },
  TestID -> "LegacyOracle-all-required-MSGCorep-cases-present"
]

VerificationTest[
  legacyOracle["Cases", "MSG13.65-B-2",
    "IterativeGenerators", "NumberOfParameters"],
  4,
  TestID -> "LegacyOracle-MSG13.65-documentation-baseline"
]

VerificationTest[
  legacyOracle["Cases", "MSG1.2-Gamma-1",
    "IterativeGenerators", "NumberOfParameters"],
  7,
  TestID -> "LegacyOracle-MSG1.2-Gamma-order2-baseline"
]

VerificationTest[
  legacyOracle["Cases", "MSG97.155-P-1+2",
    "IterativeGenerators", "NumberOfParameters"],
  9,
  TestID -> "LegacyOracle-MSG97.155-P-dimension3-baseline"
]

VerificationTest[
  Sort[
    ToExpression /@ Keys[
      legacyOracle["DimensionReferencesGETkpHam"]
    ]
  ],
  Range[10],
  TestID -> "LegacyOracle-GETkpHam-special-dimension-references"
]
