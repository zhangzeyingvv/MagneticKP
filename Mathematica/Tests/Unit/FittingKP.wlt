If[! ValueQ[$MagneticKPTestHelpersLoaded],
  Get[FileNameJoin[{
    ParentDirectory[DirectoryName[$InputFileName]],
    "Support",
    "TestHelpers.wl"
  }]]
];
LoadMagneticKPForTests[];

fittingKPHamiltonian = {
  {
    Global`fittingKPEnergy +
      Global`fittingKPSlope MagneticKP`kx,
    0
  },
  {0, 4 + Global`fittingKPEnergy/5}
};
fittingKPTrueRules = {
  Global`fittingKPEnergy -> 0.3,
  Global`fittingKPSlope -> 0.8
};
fittingKPKPoints = Table[{index/8, 0, 0}, {index, 0, 8}];
fittingKPReferenceData = Map[
  Function[k,
    {
      k,
      Sort@Eigenvalues@N[
        fittingKPHamiltonian /. fittingKPTrueRules /.
          Thread[
            {
              MagneticKP`kx,
              MagneticKP`ky,
              MagneticKP`kz
            } -> k
          ]
      ]
    }
  ],
  fittingKPKPoints
];
fittingKPResult = MagneticKP`fittingKP[
  MatrixForm[fittingKPHamiltonian],
  fittingKPReferenceData,
  Range[Length[fittingKPReferenceData]],
  <|
    Global`fittingKPEnergy -> 0.,
    Global`fittingKPSlope -> 0.
  |>
];
fittingKPFittedRules = Association[fittingKPResult["FittedParams"]];

VerificationTest[
  AssociationQ[fittingKPResult],
  True,
  TestID -> "FittingKP-valid-fit-returns-Association"
]

VerificationTest[
  {
    fittingKPResult["Objective"],
    fittingKPResult["Optimizer"],
    fittingKPResult["ParameterModel"],
    fittingKPResult["MomentumConvention"]
  },
  {
    "SumSquaredBandResiduals",
    "LevenbergMarquardt",
    "CachedAffineCoefficientMatrices",
    "DirectEigdataCoordinates"
  },
  TestID -> "FittingKP-reports-optimized-method"
]

VerificationTest[
  And[
    TrueQ[fittingKPResult["Converged"]],
    Abs[fittingKPFittedRules[Global`fittingKPEnergy] - 0.3] < 10^-5,
    Abs[fittingKPFittedRules[Global`fittingKPSlope] - 0.8] < 10^-5,
    fittingKPResult["LSQForFittedParams"] < 10^-10,
    fittingKPResult["LSQForFittedParams"] <
      fittingKPResult["LSQForInitParams"]
  ],
  True,
  TestID -> "FittingKP-recovers-direct-coordinate-model"
]

VerificationTest[
  And[
    Head[fittingKPResult["ComparisonPlot"]] === Legended,
    !FreeQ[fittingKPResult["ComparisonPlot"], "VASP"],
    !FreeQ[fittingKPResult["ComparisonPlot"], "Initial k.p"],
    !FreeQ[fittingKPResult["ComparisonPlot"], "Fitted k.p"]
  ],
  True,
  TestID -> "FittingKP-returns-comparison-plot"
]

fittingKPSelectedResult = MagneticKP`fittingKP[
  fittingKPHamiltonian,
  fittingKPReferenceData,
  Range[Length[fittingKPReferenceData]],
  fittingKPTrueRules,
  "BandSelection" -> 1,
  "EnergyWindow" -> {0.2, 1.2},
  "KPointNeighborhood" -> <|
    "Center" -> {1/2, 0, 0},
    "Radius" -> 1/2,
    "Periodic" -> False
  |>
];

VerificationTest[
  And[
    AssociationQ[fittingKPSelectedResult],
    TrueQ[fittingKPSelectedResult["SelectionRestricted"]],
    fittingKPSelectedResult["BandSelection"] === {1},
    fittingKPSelectedResult["ResidualCount"] === 9
  ],
  True,
  TestID -> "FittingKP-selection-options"
]

fittingKPWeightedResult = MagneticKP`fittingKP[
  fittingKPHamiltonian,
  fittingKPReferenceData,
  Range[Length[fittingKPReferenceData]],
  fittingKPTrueRules,
  "ResidualWeights" -> ConstantArray[
    {1., 0.5},
    Length[fittingKPReferenceData]
  ]
];

VerificationTest[
  {
    fittingKPWeightedResult["Objective"],
    fittingKPWeightedResult["Weighting"]
  },
  {"WeightedSumSquaredBandResiduals", "ExplicitMatrix"},
  TestID -> "FittingKP-residual-weights"
]

VerificationTest[
  Quiet@MagneticKP`fittingKP[
    fittingKPHamiltonian,
    fittingKPReferenceData,
    Range[Length[fittingKPReferenceData]],
    {Global`fittingKPEnergy -> 0.}
  ],
  $Failed,
  TestID -> "FittingKP-missing-initial-value-fails"
]

VerificationTest[
  Quiet@MagneticKP`fittingKP[
    fittingKPHamiltonian,
    fittingKPReferenceData,
    {0},
    fittingKPTrueRules
  ],
  $Failed,
  TestID -> "FittingKP-invalid-range-fails"
]

VerificationTest[
  Quiet@MagneticKP`fittingKP[
    {{Global`fittingKPEnergy^2}},
    {{{0, 0, 0}, {1.}}},
    {1},
    {Global`fittingKPEnergy -> 0.}
  ],
  $Failed,
  TestID -> "FittingKP-nonlinear-parameter-model-fails"
]

VerificationTest[
  Quiet@MagneticKP`fittingKP[
    fittingKPHamiltonian,
    fittingKPReferenceData,
    Range[Length[fittingKPReferenceData]],
    fittingKPTrueRules,
    MaxIterations -> 0
  ],
  $Failed,
  TestID -> "FittingKP-invalid-option-fails"
]
