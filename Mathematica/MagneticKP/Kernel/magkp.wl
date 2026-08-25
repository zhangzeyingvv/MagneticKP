(* ::Package:: *)

BeginPackage["MagneticKP`"]

Begin["`Private`"]

validKOrders[ikOrder_] := Which[
  IntegerQ[ikOrder] && ikOrder >= 0,
    Range[0, ikOrder],
  ListQ[ikOrder] && And @@ (IntegerQ[#] && # >= 0 & /@ ikOrder),
    ikOrder,
  True,
    $Failed
];

firstSymmetryValue[input_] := Module[
  {categories, values},
  categories = rawSymmetryCategories[input];
  values = Flatten[
    If[AssociationQ[#], Values[#], #] & /@ categories,
    1
  ];
  SelectFirst[
    values,
    ListQ[#] && Length[#] >= 2 && MatrixQ[#[[1]]] && ListQ[#[[2]]] &,
    Missing["NotFound"]
  ]
];

resolvedGroupInput[input_, requested_] := Module[
  {value},
  value = If[
    requested === Automatic,
    If[AssociationQ[input], Lookup[input, "GroupInput", "Generators"], "Generators"],
    requested
  ];
  If[value === "Elements", "Elements", "Generators"]
];

factorsForOrder[
  records_List,
  order_Integer?NonNegative,
  klist_List,
  vars_List
] := Map[
  Function[record,
    Module[{f},
      f = MomentumTransformMatrix[
        order,
        klist,
        record["KImage"],
        vars
      ];
      If[
        f === $Failed,
        Nothing,
        <|"F" -> f, "J" -> record["J"], "Eta" -> record["Eta"]|>
      ]
    ]
  ],
  records
];

coefficientSymbolForOrder[
  order_Integer?NonNegative,
  line_Integer?Positive
] := Switch[
  order,
  0, Subscript[MagneticKP`e, line],
  1, Subscript[MagneticKP`t, line],
  2, Subscript[MagneticKP`r, line],
  3, Subscript[MagneticKP`s, line],
  _, Subscript[System`C, order, line]
];

hamiltonianFromRowBasis[
  rowBasis_List,
  order_Integer?NonNegative,
  klist_List,
  operatorBasis_List,
  matrixDimension_Integer?Positive
] := Module[
  {productBasis},
  If[rowBasis === {}, Return[0 IdentityMatrix[matrixDimension]]];
  productBasis = Flatten[(# operatorBasis) & /@ klist, 1];
  Total @ Table[
    coefficientSymbolForOrder[order, line] Total @ MapThread[
      Times,
      {rowBasis[[line]], productBasis}
    ],
    {line, Length[rowBasis]}
  ]
];

constructKP[
  ikOrder_,
  input_,
  solver_,
  groupInput_,
  maxGroupOrder_,
  returnDiagnostics_
] := Module[
  {
    kOrder, firstValue, ndim, spatialDim, vars, fullBasisData,
    continuousData, continuousPresent, discreteInput, basisData, records,
    symham, numberOfParameters, diagnostics, order, klist, factors,
    ambient, chosenSolver, solverResult, rowBasis, result, resolvedInput
  },
  kOrder = validKOrders[ikOrder];
  If[kOrder === $Failed, Return[$Failed]];
  firstValue = firstSymmetryValue[input];
  If[MissingQ[firstValue], Return[$Failed]];

  ndim = Length[firstValue[[1]]];
  spatialDim = Length[firstValue[[2]]];
  vars = Switch[
    spatialDim,
    1, {kz},
    2, {kx, ky},
    3, {kx, ky, kz},
    _, Return[$Failed]
  ];

  fullBasisData = OperatorBasisData[ndim];
  If[fullBasisData === $Failed, Return[$Failed]];
  continuousData = PrepareContinuousSymmetry[
    input,
    fullBasisData,
    vars
  ];
  If[continuousData === $Failed, Return[$Failed]];
  continuousPresent = TrueQ[continuousData["Present"]];
  discreteInput = continuousData["DiscreteInput"];
  basisData = continuousData["BasisData"];
  records = PrepareSymmetryRecords[
    discreteInput,
    basisData,
    continuousPresent,
    If[
      continuousPresent,
      continuousData["Generator"],
      None
    ]
  ];
  If[records === $Failed, Return[$Failed]];
  If[records === {} && Not[continuousPresent], Return[$Failed]];
  If[
    Not[And @@ (Length[#["KImage"]] == spatialDim & /@ records)],
    Return[$Failed]
  ];

  resolvedInput = resolvedGroupInput[discreteInput, groupInput];
  symham = 0 IdentityMatrix[ndim];
  numberOfParameters = 0;
  diagnostics = {};

  Do[
    klist = HomogeneousMonomialBasis[order, spatialDim];
    factors = factorsForOrder[records, order, klist, vars];
    If[Length[factors] =!= Length[records], Return[$Failed]];
    ambient = Length[klist] basisData["Count"];
    chosenSolver = Which[
      solver === "Automatic",
        "Iterative",
      True,
        solver
    ];
    solverResult = Which[
      factors === {},
        <|
          "Basis" -> IdentityMatrix[ambient],
          "Success" -> True,
          "FallbackUsed" -> False,
          "GroupOrder" -> 0
        |>,
      chosenSolver === "Reynolds",
        ReynoldsInvariantBasis[
          factors,
          "GroupInput" -> resolvedInput,
          "MaxGroupOrder" -> maxGroupOrder,
          "Verify" -> True
        ],
      chosenSolver === "DirectProduct",
        <|
          "Basis" -> StackedInvariantBasis[factors],
          "Success" -> True,
          "FallbackUsed" -> False,
          "GroupOrder" -> Missing["NotUsed"]
        |>,
      True,
        <|
          "Basis" -> IterativeInvariantBasis[factors],
          "Success" -> True,
          "FallbackUsed" -> False,
          "GroupOrder" -> Missing["NotUsed"]
        |>
    ];
    rowBasis = Simplify[solverResult["Basis"]];
    symham += hamiltonianFromRowBasis[
      rowBasis,
      order,
      klist,
      basisData["Basis"],
      ndim
    ];
    numberOfParameters += Length[rowBasis];
    AppendTo[
      diagnostics,
      <|
        "Order" -> order,
        "Solver" -> If[factors === {}, "None", chosenSolver],
        "AmbientDimension" -> ambient,
        "InvariantDimension" -> Length[rowBasis],
        "ContinuousSymmetry" -> If[
          continuousPresent,
          "CInfinity",
          None
        ],
        "OriginalOperatorDimension" ->
          continuousData["OriginalOperatorDimension"],
        "ReducedOperatorDimension" ->
          continuousData["ReducedOperatorDimension"],
        "GroupInput" -> resolvedInput,
        "GroupOrder" -> Lookup[
          solverResult,
          "GroupOrder",
          Missing["NotAvailable"]
        ],
        "ReynoldsFallback" -> Lookup[
          solverResult,
          "FallbackUsed",
          False
        ],
        "ElementClosureCompleted" -> Lookup[
          solverResult,
          "ElementClosureCompleted",
          False
        ],
        "Reason" -> Lookup[solverResult, "Reason", None]
      |>
    ],
    {order, kOrder}
  ];

  result = <|
    "ham" -> symham,
    "korder" -> kOrder,
    "dim" -> ndim,
    "NumberOfParameters" -> numberOfParameters
  |>;
  If[
    TrueQ[returnDiagnostics],
    AssociateTo[result, "Diagnostics" -> diagnostics]
  ];
  result
];

Options[IterativelySimplify] = {
  "GroupInput" -> Automatic,
  "MaxGroupOrder" -> 512,
  "ReturnDiagnostics" -> False
};

IterativelySimplify[ikOrder_, input_, OptionsPattern[]] := constructKP[
  ikOrder,
  input,
  "Iterative",
  OptionValue["GroupInput"],
  OptionValue["MaxGroupOrder"],
  OptionValue["ReturnDiagnostics"]
];

Options[DirectProductDecomposition] = Options[IterativelySimplify];

DirectProductDecomposition[ikOrder_, input_, OptionsPattern[]] := constructKP[
  ikOrder,
  input,
  "DirectProduct",
  OptionValue["GroupInput"],
  OptionValue["MaxGroupOrder"],
  OptionValue["ReturnDiagnostics"]
];

(* The original 2016 entry point computes one homogeneous order, whereas
   IterativelySimplify[n,...] computes every order from zero through n. *)
GETkpHam[ikOrder_Integer?NonNegative, input_] := Module[
  {result},
  result = constructKP[
    {ikOrder},
    input,
    "Iterative",
    Automatic,
    512,
    False
  ];
  If[
    AssociationQ[result],
    AssociateTo[result, "korder" -> ikOrder]
  ];
  result
];

GETkpHam[___] := $Failed;

legacy2016InputQ[input_] := AssociationQ[input] &&
  Not[KeyExistsQ[input, "AntisymmetryUnitaryTest"]] &&
  Not[KeyExistsQ[input, "AntisymmetryAnitunitaryTest"]] &&
  Not[KeyExistsQ[input, "AntisymmetryAntiunitaryTest"]];

Options[kpHam] = {
  "Method" -> "IterativeSimplification",
  "GroupInput" -> Automatic,
  "MaxGroupOrder" -> 512,
  "ReturnDiagnostics" -> False
};

kpHam[ikOrder_, input_, OptionsPattern[]] := Module[
  {method, groupInput, maxGroupOrder, returnDiagnostics},
  method = OptionValue["Method"];
  groupInput = OptionValue["GroupInput"];
  maxGroupOrder = OptionValue["MaxGroupOrder"];
  returnDiagnostics = OptionValue["ReturnDiagnostics"];
  Which[
    method === "IterativeSimplification2016",
      If[
        IntegerQ[ikOrder] && ikOrder >= 0 && legacy2016InputQ[input],
        GETkpHam[ikOrder, input],
        $Failed
      ],
    method === "IterativeSimplification",
      constructKP[
        ikOrder, input, "Iterative",
        groupInput, maxGroupOrder, returnDiagnostics
      ],
    method === "DirectProductDecomposition",
      constructKP[
        ikOrder, input, "DirectProduct",
        groupInput, maxGroupOrder, returnDiagnostics
      ],
    method === "Reynolds",
      constructKP[
        ikOrder, input, "Reynolds",
        groupInput, maxGroupOrder, returnDiagnostics
      ],
    method === Automatic || method === "Automatic",
      constructKP[
        ikOrder, input, "Automatic",
        groupInput, maxGroupOrder, returnDiagnostics
      ],
    True,
      $Failed
  ]
];

End[]

EndPackage[]
