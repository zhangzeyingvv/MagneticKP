If[! ValueQ[$MagneticKPTestHelpersLoaded],
  Get[FileNameJoin[{
    ParentDirectory[DirectoryName[$InputFileName]],
    "Support",
    "TestHelpers.wl"
  }]]
];
LoadMagneticKPForTests[];

bandPlotPath = {
  {{{0, 0, 0}, {1/2, 0, 0}}, {"G", "X"}},
  {{{1/2, 0, 0}, {1/2, 1/2, 0}}, {"X", "M"}}
};
bandPlotHamiltonian = {
  {Global`bandMass + MagneticKP`kx, 0},
  {0, -Global`bandMass - MagneticKP`kx}
};
bandPlotRules = {Global`bandMass -> 1/5};

VerificationTest[
  Head@MagneticKP`bandplot[
    bandPlotPath,
    4,
    bandPlotHamiltonian,
    bandPlotRules
  ],
  Graphics,
  TestID -> "BandPlot-valid-path-returns-Graphics"
]

VerificationTest[
  Head@MagneticKP`bandplot[
    bandPlotPath,
    4,
    bandPlotHamiltonian,
    bandPlotRules,
    "PlotRange" -> {-2, 2},
    "YTicks" -> {{-1, "low"}, {0, "zero"}, {1, "high"}},
    FontSize -> 16,
    FontFamily -> "Helvetica",
    ImageSize -> 320
  ],
  Graphics,
  TestID -> "BandPlot-display-options-return-Graphics"
]

VerificationTest[
  PrivateHook["bandPathBoundaryLabels"]@{
    {{{0, 0, 0}, {1/2, 0, 0}}, {"G", "X"}},
    {{{0, 1/2, 0}, {1/2, 1/2, 0}}, {"M", "R"}}
  },
  {"G", "X|M", "R"},
  TestID -> "BandPlot-discontinuous-path-labels"
]

bandPlotDataResult = PrivateHook["bandPlotData"][
  bandPlotPath,
  4,
  bandPlotHamiltonian,
  bandPlotRules
];

VerificationTest[
  {
    Length[bandPlotDataResult],
    {First[#][[1]], Last[#][[1]]} & /@ bandPlotDataResult
  },
  {
    4,
    {{0, 4}, {0, 4}, {4, 8}, {4, 8}}
  },
  TestID -> "BandPlot-segment-data-boundaries"
]

VerificationTest[
  Quiet@MagneticKP`bandplot[
    bandPlotPath,
    0,
    bandPlotHamiltonian,
    bandPlotRules
  ],
  $Failed,
  TestID -> "BandPlot-invalid-sample-count-fails"
]

VerificationTest[
  Quiet@MagneticKP`bandplot[
    {{{0, 0, 0}, {1/2, 0, 0}}},
    4,
    bandPlotHamiltonian,
    bandPlotRules
  ],
  $Failed,
  TestID -> "BandPlot-invalid-path-fails"
]

VerificationTest[
  Quiet@MagneticKP`bandplot[
    bandPlotPath,
    4,
    bandPlotHamiltonian,
    {}
  ],
  $Failed,
  TestID -> "BandPlot-unresolved-parameter-fails"
]
