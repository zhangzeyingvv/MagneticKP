(* Compare the Python MSGCorep catalog adapter with Mathematica interfaceRep. *)

$HistoryLength = 0;

testsRoot = ExpandFileName[
  ParentDirectory[DirectoryName[$InputFileName]]
];
Get[FileNameJoin[{testsRoot, "Support", "TestHelpers.wl"}]];

If[! LoadMagneticKPForTests[] || ! LoadMSGCorepForTests[],
  Print["MAGNETICKP_PYTHON_PARITY_ERROR_BEGIN"];
  Print["Unable to load MagneticKP, MSGCorep, or SpaceGroupIrep."];
  Print["MAGNETICKP_PYTHON_PARITY_ERROR_END"];
  Exit[2]
];

arguments = If[ListQ[$ScriptCommandLine] && Length[$ScriptCommandLine] > 0,
  Rest[$ScriptCommandLine],
  {}
];
payloadPath = Environment["MAGNETICKP_PARITY_PAYLOAD"];
If[! StringQ[payloadPath] || StringLength[payloadPath] == 0,
  payloadPath = If[Length[arguments] == 1, First[arguments], ""]
];
If[! FileExistsQ[payloadPath],
  Print["MAGNETICKP_PYTHON_PARITY_ERROR_BEGIN"];
  Print["Expected one existing JSON payload path."];
  Print["MAGNETICKP_PYTHON_PARITY_ERROR_END"];
  Exit[2]
];

payload = Import[payloadPath, "RawJSON"];
Print["PARITY_STAGE:payload-loaded"];

ClearAll[ParsePythonHamiltonian];
ParsePythonHamiltonian[source_String] := Quiet @ Check[
  ToExpression[source, InputForm] /. {
    Global`kx -> MagneticKP`kx,
    Global`ky -> MagneticKP`ky,
    Global`kz -> MagneticKP`kz
  },
  $Failed
];

ClearAll[EvaluateParityCase];
EvaluateParityCase[spec_Association] := Module[
  {
    input, mmaResult, mmaHamiltonian, pythonHamiltonian,
    mmaInvariantDimensions, mmaParameterCount,
    pythonInvariantDimensions, pythonParameterCount,
    subspaceEqual, exponentSetsEqual, dimensionsEqual,
    resolved, pass
  },
  Print["PARITY_STAGE:", spec["Name"], ":interfaceRep"];
  input = Quiet @ Check[
    MagneticKP`interfaceRep[
      spec["MSG"],
      ToExpression[spec["KExpression"], InputForm],
      spec["MMAReps"],
      "CartesianCoordinates" -> spec["CartesianCoordinates"],
      "CalculateGenerators" -> True
    ],
    $Failed
  ];
  If[! AssociationQ[input],
    Return[<|"Pass" -> False, "Reason" -> "interfaceRepFailed"|>]
  ];
  Print["PARITY_STAGE:", spec["Name"], ":kpHam"];
  mmaResult = Quiet @ Check[
    MagneticKP`kpHam[
      spec["Order"],
      input,
      "Method" -> "IterativeSimplification",
      "ReturnDiagnostics" -> True
    ],
    $Failed
  ];
  If[! AssociationQ[mmaResult],
    Return[<|"Pass" -> False, "Reason" -> "kpHamFailed"|>]
  ];
  Print["PARITY_STAGE:", spec["Name"], ":parse-python"];
  pythonHamiltonian = ParsePythonHamiltonian[
    spec["PythonHamiltonian"]
  ];
  If[! MatrixQ[pythonHamiltonian],
    Return[<|
      "Pass" -> False,
      "Reason" -> "PythonHamiltonianParseFailed"
    |>]
  ];
  Print["PARITY_STAGE:", spec["Name"], ":compare"];
  mmaHamiltonian = mmaResult["ham"];
  mmaInvariantDimensions = Lookup[
    mmaResult["Diagnostics"],
    "InvariantDimension",
    Missing["NotAvailable"]
  ];
  pythonInvariantDimensions = spec["PythonInvariantDimensions"];
  mmaParameterCount = mmaResult["NumberOfParameters"];
  pythonParameterCount = spec["PythonParameterCount"];
  dimensionsEqual = Dimensions[mmaHamiltonian] ===
    Dimensions[pythonHamiltonian];
  resolved = ResolvedHamiltonianQ[mmaHamiltonian] &&
    ResolvedHamiltonianQ[pythonHamiltonian];
  subspaceEqual = dimensionsEqual &&
    HamiltonianSubspaceEqualQ[mmaHamiltonian, pythonHamiltonian];
  exponentSetsEqual = HamiltonianExponentSet[mmaHamiltonian] ===
    HamiltonianExponentSet[pythonHamiltonian];
  pass = TrueQ[
    dimensionsEqual && resolved && subspaceEqual && exponentSetsEqual &&
    mmaParameterCount == pythonParameterCount &&
    mmaInvariantDimensions === pythonInvariantDimensions
  ];
  Print["PARITY_STAGE:", spec["Name"], ":done"];
  <|
    "Pass" -> pass,
    "MatrixDimension" -> First[Dimensions[mmaHamiltonian]],
    "MMAParameterCount" -> mmaParameterCount,
    "PythonParameterCount" -> pythonParameterCount,
    "MMAInvariantDimensions" -> mmaInvariantDimensions,
    "PythonInvariantDimensions" -> pythonInvariantDimensions,
    "SubspaceEqual" -> subspaceEqual,
    "ExponentSetsEqual" -> exponentSetsEqual,
    "ResolvedHamiltonians" -> resolved,
    "UnitaryGeneratorCount" -> Length[Lookup[input, "Unitary", <||>]],
    "AntiunitaryGeneratorCount" -> Length[
      Lookup[input, "Anitunitary", <||>]
    ]
  |>
];

report = Association @ Map[
  #["Name"] -> EvaluateParityCase[#] &,
  payload["Cases"]
];

Print["MAGNETICKP_PYTHON_PARITY_BEGIN"];
Print[ExportString[report, "RawJSON", "Compact" -> True]];
Print["MAGNETICKP_PYTHON_PARITY_END"];

Exit[If[And @@ Lookup[Values[report], "Pass", False], 0, 1]];
