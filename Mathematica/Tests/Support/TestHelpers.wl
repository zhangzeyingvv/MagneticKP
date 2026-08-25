(* Shared helpers for the MagneticKP regression tests. *)

If[! ValueQ[$MagneticKPTestHelpersLoaded],
  $MagneticKPTestHelpersLoaded = True;

  $MagneticKPTestsRoot = ExpandFileName[
    ParentDirectory[DirectoryName[$InputFileName]]
  ];
  $MagneticKPMathematicaRoot =
    ParentDirectory[$MagneticKPTestsRoot];
  $MagneticKPContainerRoot =
    ParentDirectory[$MagneticKPMathematicaRoot];
  $MagneticKPStableRoot = If[
    DirectoryQ[FileNameJoin[{$MagneticKPContainerRoot, "Python-Web"}]],
    $MagneticKPContainerRoot,
    FileNameJoin[{$MagneticKPContainerRoot, "MagneticKP"}]
  ];
  $MagneticKPRepositoryRoot = If[
    $MagneticKPStableRoot === $MagneticKPContainerRoot,
    ParentDirectory[$MagneticKPContainerRoot],
    $MagneticKPContainerRoot
  ];
  $MagneticKPPackageRoot =
    FileNameJoin[{$MagneticKPMathematicaRoot, "MagneticKP"}];
  $MagneticKPPythonRoot =
    FileNameJoin[{$MagneticKPStableRoot, "Python-Web"}];
  $MagneticKPLegacyRoot =
    FileNameJoin[{
      ParentDirectory[$MagneticKPRepositoryRoot],
      "MagneticKP-Mathmatica"
    }];

  ClearAll[LoadMagneticKPForTests];
  LoadMagneticKPForTests[] := Module[{},
    If[! MemberQ[$Path, $MagneticKPMathematicaRoot],
      PrependTo[$Path, $MagneticKPMathematicaRoot]
    ];
    If[! MemberQ[$Packages, "MagneticKP`"],
      Get[FileNameJoin[{$MagneticKPPackageRoot, "Kernel", "init.m"}]]
    ];
    MemberQ[$Packages, "MagneticKP`"]
  ];

  ClearAll[LoadMSGCorepForTests];
  LoadMSGCorepForTests[] := Module[{},
    Quiet[
      Needs["MSGCorep`"];
      Needs["SpaceGroupIrep`"]
    ];
    MemberQ[$Packages, "MSGCorep`"] &&
      MemberQ[$Packages, "SpaceGroupIrep`"]
  ];

  ClearAll[PrivateHook, PrivateHookAvailableQ];
  PrivateHook[name_String] :=
    Symbol["MagneticKP`Private`" <> name];
  PrivateHookAvailableQ[name_String] := Module[{fullName, values},
    fullName = "MagneticKP`Private`" <> name;
    If[! NameQ[fullName], Return[False]];
    values = ToExpression[
      fullName,
      InputForm,
      Function[symbol, DownValues[symbol], HoldAll]
    ];
    Length[values] > 0
  ];

  ClearAll[KVariables, KMonomialBasis];
  KVariables[1] := {MagneticKP`kz};
  KVariables[2] := {MagneticKP`kx, MagneticKP`ky};
  KVariables[3] := {
    MagneticKP`kx, MagneticKP`ky, MagneticKP`kz
  };

  KMonomialBasis[order_Integer?NonNegative, 1] :=
    MagneticKP`kz^#[[1]] & /@
      Sort @ Flatten[
        Permutations /@
          IntegerPartitions[order, {1}, Range[0, order]],
        1
      ];
  KMonomialBasis[order_Integer?NonNegative, 2] :=
    MagneticKP`ky^#[[1]] MagneticKP`kx^#[[2]] & /@
      Sort @ Flatten[
        Permutations /@
          IntegerPartitions[order, {2}, Range[0, order]],
        1
      ];
  KMonomialBasis[order_Integer?NonNegative, 3] :=
    MagneticKP`kz^#[[1]] MagneticKP`ky^#[[2]]
      MagneticKP`kx^#[[3]] & /@
      Sort @ Flatten[
        Permutations /@
          IntegerPartitions[order, {3}, Range[0, order]],
        1
      ];

  ClearAll[LinearKImage];
  LinearKImage[a_?MatrixQ] := a . KVariables[Length[a]];

  ClearAll[TryEvaluatedMatrixCalls];
  TryEvaluatedMatrixCalls[calls_List] := Module[{value},
    Do[
      value = Quiet @ Check[ReleaseHold[call], $Failed];
      If[value =!= $Failed && MatrixQ[value], Return[value]],
      {call, calls}
    ];
    $Failed
  ];

  ClearAll[CallOperatorBasisData];
  CallOperatorBasisData[argument_] := Module[
    {hook = PrivateHook["OperatorBasisData"], result},
    result = Quiet @ Check[hook[argument], $Failed];
    If[AssociationQ[result], result, $Failed]
  ];

  ClearAll[CallMomentumTransform];
  CallMomentumTransform[
    order_Integer?NonNegative,
    kImage_List
  ] := Module[
    {hook = PrivateHook["MomentumTransformMatrix"], dim, basis, vars},
    dim = Length[kImage];
    basis = KMonomialBasis[order, dim];
    vars = KVariables[dim];
    Quiet @ Check[
      hook[order, basis, kImage, vars],
      $Failed
    ]
  ];

  ClearAll[CallAdjointRepresentation];
  CallAdjointRepresentation[
    d_?MatrixQ,
    basisData_Association,
    antiunitary_: False
  ] := Module[
    {hook = PrivateHook["AdjointRepresentationMatrix"]},
    Quiet @ Check[
      hook[d, basisData, antiunitary],
      $Failed
    ]
  ];

  ClearAll[CallReynoldsInvariantBasis];
  CallReynoldsInvariantBasis[
    reps_List,
    characters_: Automatic
  ] := Module[
    {hook = PrivateHook["ReynoldsInvariantBasis"], value, chars, factors},
    chars = If[characters === Automatic,
      ConstantArray[1, Length[reps]],
      characters
    ];
    factors = MapThread[
      <|"F" -> #1, "J" -> {{1}}, "Eta" -> #2|> &,
      {reps, chars}
    ];
    value = Quiet @ Check[
      hook[
        factors,
        "GroupInput" -> "Elements",
        "Verify" -> True
      ],
      $Failed
    ];
    If[AssociationQ[value],
      value = Lookup[value, "Basis", $Failed]
    ];
    value
  ];

  ClearAll[LegacyAdjointMatrix];
  LegacyAdjointMatrix[
    d_?MatrixQ,
    basis_List,
    antiunitary_: False
  ] := Module[{dinv = Inverse[d], action, inner},
    action[m_] := If[
      TrueQ[antiunitary],
      d . Conjugate[m] . dinv,
      d . m . dinv
    ];
    inner[a_, b_] := Tr[ConjugateTranspose[a] . b];
    Table[
      inner[action[basis[[i]]], basis[[j]]] /
        inner[basis[[j]], basis[[j]]],
      {j, Length[basis]}, {i, Length[basis]}
    ]
  ];

  ClearAll[ReferenceAdjointMatrix];
  ReferenceAdjointMatrix[
    d_?MatrixQ,
    basis_List,
    antiunitary_: False
  ] := Module[
    {rows, columns, gram, dual, dinv, images, imageRows},
    rows = Flatten /@ basis;
    columns = Transpose[rows];
    gram = Conjugate[rows] . columns;
    If[MatrixRank[gram] < Length[basis], Return[$Failed]];
    dual = LinearSolve[gram, Conjugate[rows]];
    dinv = Inverse[d];
    images = If[
      TrueQ[antiunitary],
      d . Conjugate[#] . dinv & /@ basis,
      d . # . dinv & /@ basis
    ];
    imageRows = Flatten /@ images;
    FullSimplify[dual . Transpose[imageRows]]
  ];

  ClearAll[ReferenceBasisData];
  ReferenceBasisData[basis_List, basisID_: "TestBasis"] := Module[
    {rows, columns, gram, dual, conjugation},
    rows = Flatten /@ basis;
    columns = Transpose[rows];
    gram = Conjugate[rows] . columns;
    If[MatrixRank[gram] < Length[basis],
      Return[<|
        "Basis" -> basis,
        "Gram" -> gram,
        "Singular" -> True,
        "BasisID" -> basisID
      |>]
    ];
    dual = LinearSolve[gram, Conjugate[rows]];
    conjugation = FullSimplify[dual . Conjugate[columns]];
    <|
      "Basis" -> basis,
      "Rows" -> rows,
      "Columns" -> columns,
      "Gram" -> gram,
      "Dual" -> dual,
      "Conjugation" -> SparseArray[conjugation],
      "BasisID" -> basisID
    |>
  ];

  ClearAll[LegacySolveMomentumTransform];
  LegacySolveMomentumTransform[
    order_Integer?NonNegative,
    kImage_List
  ] := Module[
    {dim, vars, basis, transformed, tmp, unknown},
    dim = Length[kImage];
    vars = KVariables[dim];
    basis = KMonomialBasis[order, dim];
    transformed = basis /. Thread[vars -> kImage];
    unknown = Table[
      tmp[i, j],
      {i, Length[basis]}, {j, Length[basis]}
    ];
    unknown /. First @ SolveAlways[
      basis == transformed . unknown,
      vars
    ]
  ];

  ClearAll[MatrixExactEqualQ, MatrixNearQ];
  MatrixExactEqualQ[a_, b_] :=
    Dimensions[a] === Dimensions[b] &&
      TrueQ[FullSimplify[a == b]];
  MatrixNearQ[a_, b_, tolerance_: 10^-11] := Module[
    {scale},
    If[Dimensions[a] =!= Dimensions[b], Return[False]];
    scale = Max[1., Max[Abs @ N @ Flatten @ b]];
    Max[Abs @ N @ Flatten[a - b]] <= tolerance scale
  ];

  ClearAll[NonzeroRows, ImageBasis, SubspaceEqualQ];
  NonzeroRows[m_] := Select[
    Normal[m],
    Not @ TrueQ[PossibleZeroQ[Total[Abs[#]]]] &
  ];
  ImageBasis[projector_?MatrixQ] :=
    NonzeroRows @ RowReduce @ Transpose[projector];
  SubspaceEqualQ[a_List, b_List] := Module[
    {ra, rb, joined},
    ra = If[a === {}, 0, MatrixRank[a]];
    rb = If[b === {}, 0, MatrixRank[b]];
    joined = Join[a, b];
    ra == rb == If[joined === {}, 0, MatrixRank[joined]]
  ];

  ClearAll[IterativeInvariantBasis, StackedInvariantBasis];
  IterativeInvariantBasis[reps_List] := Module[
    {ambient, identity},
    If[reps === {}, Return[{}]];
    ambient = Length[First[reps]];
    identity = IdentityMatrix[ambient, SparseArray];
    Fold[
      Function[{basis, rep},
        If[basis === {},
          {},
          With[{kernel = NullSpace[
            (rep - identity) . Transpose[basis]
          ]},
            If[kernel === {}, {}, kernel . basis]
          ]
        ]
      ],
      IdentityMatrix[ambient],
      reps
    ]
  ];
  StackedInvariantBasis[reps_List] := Module[
    {ambient, identity},
    If[reps === {}, Return[{}]];
    ambient = Length[First[reps]];
    identity = IdentityMatrix[ambient, SparseArray];
    NullSpace[Join @@ ((# - identity) & /@ reps)]
  ];

  ClearAll[HamiltonianParameters];
  ClearAll[ResolvedHamiltonianQ];
  ResolvedHamiltonianQ[ham_?MatrixQ] := FreeQ[
    ham,
    _Return |
      _NullSpace |
      _MagneticKP`Private`hamiltonianFromRowBasis
  ];

  HamiltonianParameters[ham_?MatrixQ] := Module[{variables},
    variables = Quiet @ Check[Variables[Flatten[Normal[ham]]], {}];
    DeleteCases[
      variables,
      MagneticKP`kx | MagneticKP`ky | MagneticKP`kz
    ]
  ];

  ClearAll[HamiltonianExponentSet];
  HamiltonianExponentSet[ham_?MatrixQ] := Module[
    {vars, rules},
    vars = {
      MagneticKP`kx, MagneticKP`ky, MagneticKP`kz
    };
    rules = CoefficientRules[Expand[#], vars] & /@
      Flatten[Normal[ham]];
    Sort @ DeleteDuplicates @ Flatten[Keys /@ rules, 1]
  ];

  ClearAll[HamiltonianCoordinateRows];
  HamiltonianCoordinateRows[
    ham_?MatrixQ,
    exponents_: Automatic
  ] := Module[
    {vars, params, exps, coefficientMatrix, ruleTables},
    vars = {
      MagneticKP`kx, MagneticKP`ky, MagneticKP`kz
    };
    params = HamiltonianParameters[ham];
    exps = If[exponents === Automatic,
      HamiltonianExponentSet[ham],
      exponents
    ];
    Table[
      coefficientMatrix = Map[
        Coefficient[#, parameter] &,
        Normal[ham],
        {2}
      ];
      ruleTables = Association @ CoefficientRules[
          Expand[#],
          vars
        ] & /@ Flatten[coefficientMatrix];
      Flatten @ Table[
        Lookup[ruleTables, exponent, 0],
        {exponent, exps}
      ],
      {parameter, params}
    ]
  ];

  ClearAll[HamiltonianSubspaceEqualQ];
  HamiltonianSubspaceEqualQ[hamA_?MatrixQ, hamB_?MatrixQ] :=
    Module[{exponents, a, b},
      exponents = Sort @ DeleteDuplicates @ Join[
        HamiltonianExponentSet[hamA],
        HamiltonianExponentSet[hamB]
      ];
      a = HamiltonianCoordinateRows[hamA, exponents];
      b = HamiltonianCoordinateRows[hamB, exponents];
      SubspaceEqualQ[a, b]
    ];

  ClearAll[ExtractCompressedPayload, RunLegacyOracle];
  ExtractCompressedPayload[text_String] := Module[
    {startMarker, endMarker, start, end, payload},
    startMarker = "MAGNETICKP_ORACLE_BEGIN";
    endMarker = "MAGNETICKP_ORACLE_END";
    start = StringPosition[text, startMarker];
    end = StringPosition[text, endMarker];
    If[start === {} || end === {}, Return[$Failed]];
    payload = StringTrim @ StringTake[
      text,
      {
        start[[-1, 2]] + 1,
        First[end][[1]] - 1
      }
    ];
    Quiet @ Check[Uncompress[payload], $Failed]
  ];

  RunLegacyOracle[] := Module[{script, process},
    script = FileNameJoin[{
      $MagneticKPTestsRoot,
      "Support",
      "LegacyOracle.wls"
    }];
    process = RunProcess[{"wolframscript", "-file", script}];
    If[process["ExitCode"] =!= 0,
      Return[<|
        "Failure" -> "LegacyOracleProcess",
        "ExitCode" -> process["ExitCode"],
        "StandardError" -> process["StandardError"],
        "StandardOutput" -> process["StandardOutput"]
      |>]
    ];
    ExtractCompressedPayload[process["StandardOutput"]]
  ];

  ClearAll[GetLegacyOracleData];
  GetLegacyOracleData[] := If[
    ValueQ[$MagneticKPLegacyOracleData],
    $MagneticKPLegacyOracleData,
    $MagneticKPLegacyOracleData = RunLegacyOracle[]
  ];

  ClearAll[RunIterativeKpHam, RunReynoldsKpHam];
  RunIterativeKpHam[order_, input_] :=
    MagneticKP`kpHam[
      order,
      input,
      "Method" -> "IterativeSimplification"
    ];
  RunReynoldsKpHam[order_, input_, groupInput_: Automatic] :=
    Module[{value},
      value = Quiet @ CheckAbort[
        MagneticKP`kpHam[
          order,
          input,
          "Method" -> "Reynolds",
          "GroupInput" -> groupInput
        ],
        $Failed
      ];
      If[! AssociationQ[value],
        value = Quiet @ CheckAbort[
          MagneticKP`kpHam[
            order,
            input,
            "Method" -> "Reynolds"
          ],
          $Failed
        ]
      ];
      value
    ];
]
