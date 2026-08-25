$HistoryLength = 0;

testsRoot = ExpandFileName[DirectoryName[$InputFileName]];
Get[FileNameJoin[{testsRoot, "Support", "TestHelpers.wl"}]];

If[! LoadMagneticKPForTests[],
  Print["Unable to load MagneticKP from ", $MagneticKPPackageRoot];
  Exit[2]
];

scriptArguments = Rest[$ScriptCommandLine];
suiteFilter = If[
  scriptArguments === {},
  "",
  First[scriptArguments]
];
testFiles = Sort @ FileNames[
  "*.wlt",
  {
    FileNameJoin[{testsRoot, "Unit"}],
    FileNameJoin[{testsRoot, "Regression"}],
    FileNameJoin[{testsRoot, "Integration"}]
  }
];

If[StringQ[suiteFilter] && StringLength[suiteFilter] > 0,
  testFiles = Select[
    testFiles,
    StringContainsQ[#, suiteFilter, IgnoreCase -> True] &
  ]
];

If[testFiles === {},
  Print["No .wlt files selected."];
  Exit[2]
];

Print["Running MagneticKP test files:"];
Scan[Print["  ", #] &, testFiles];

report = TestReport[testFiles];
failed = report["TestsFailedCount"];
succeeded = report["TestsSucceededCount"];
expected = Total[
  StringCount[Import[#, "Text"], "VerificationTest["] & /@
    testFiles
];
Print[
  "MagneticKP tests: ",
  succeeded, " succeeded, ",
  failed, " failed; ",
  expected, " expected."
];
If[IntegerQ[failed] && failed > 0,
  Print["Failed test objects:"];
  Scan[Print, report["TestsFailed"]];
  If[ValueQ[$MSGCorepCaseDiagnostics],
    Print[
      "MSGCorep case diagnostics: ",
      InputForm[$MSGCorepCaseDiagnostics]
    ]
  ];
  If[ValueQ[$MSGCorepDimensionDiagnostics],
    Print[
      "MSGCorep dimension diagnostics: ",
      InputForm[$MSGCorepDimensionDiagnostics]
    ]
  ]
];
Exit[
  If[
    IntegerQ[failed] && IntegerQ[succeeded] &&
      failed == 0 && succeeded == expected,
    0,
    1
  ]
]
