(* ::Package:: *)

BeginPackage["MagneticKP`"]

Begin["`Private`"]

ClearAll[
  OperatorBasisData,
  HomogeneousMonomialBasis,
  MomentumTransformMatrix,
  AdjointRepresentationMatrix,
  PrepareContinuousSymmetry,
  ContinuousNormalizerQ,
  PrepareSymmetryRecords,
  IterativeInvariantBasis,
  StackedInvariantBasis,
  ReynoldsInvariantBasis
]

$OperatorBasisCache = <||>;
$MomentumTransformCache = <||>;
$AdjointRepresentationCache = <||>;

CoreAlgorithms::basis = "The supplied operator matrices do not form a linearly independent basis.";
CoreAlgorithms::dimension = "Expected a positive matrix dimension, a non-negative polynomial order, and a spatial dimension from 1 to 3.";
CoreAlgorithms::kimage = "The momentum image `1` is not an invertible linear transformation of `2`.";
CoreAlgorithms::symmetry = "No valid symmetry record was found in the input.";
CoreAlgorithms::continuous = "Invalid internal CInfinity symmetry: `1`.";
CoreAlgorithms::normalizer = "The discrete representation matrix `1` does not normalize the CInfinity generator.";
ReynoldsInvariantBasis::fallback = "Reynolds construction fell back to the iterative generator solver: `1`.";

zeroScalarQ[x_] := TrueQ[PossibleZeroQ[x]] ||
  (InexactNumberQ[x] && TrueQ[Chop[x] == 0]);

zeroMatrixQ[m_] := And @@ (zeroScalarQ /@ Flatten[Normal[m]]);

nonzeroRows[m_] := Select[
  Normal[m],
  Not[And @@ (zeroScalarQ /@ #)] &
];

inexactMatrixQ[m_] := AnyTrue[Flatten[Normal[m]], InexactNumberQ];

canonicalScalar[x_] := Quiet @ Check[
  If[FreeQ[x, _Symbol] && Not[InexactNumberQ[x]], RootReduce[x], Simplify[x]],
  x
];

matrixCacheKey[m_] := ToString[
  {
    Dimensions[m],
    ArrayRules[SparseArray[m]] /. Rule[p_, value_] :> Rule[p, canonicalScalar[value]]
  },
  InputForm
];

projectiveMatrixRepresentative[d_] := Module[
  {flat, pivot},
  flat = Flatten[Normal[d]];
  pivot = SelectFirst[flat, Not[zeroScalarQ[#]] &, Missing["NotFound"]];
  If[MissingQ[pivot], Return[$Failed]];
  Simplify[d/pivot]
];

matrixUnit[n_, a_, b_] := SparseArray[{{a, b} -> 1}, {n, n}];

legacyGellMannBasis[] := {
  IdentityMatrix[3],
  {{0, -I, 0}, {I, 0, 0}, {0, 0, 0}},
  {{0, 0, -I}, {0, 0, 0}, {I, 0, 0}},
  {{0, 0, 0}, {0, 0, -I}, {0, I, 0}},
  {{0, 1, 0}, {1, 0, 0}, {0, 0, 0}},
  {{1, 0, 0}, {0, -1, 0}, {0, 0, 0}},
  {{0, 0, 1}, {0, 0, 0}, {1, 0, 0}},
  {{0, 0, 0}, {0, 0, 1}, {0, 1, 0}},
  1/Sqrt[3] {{1, 0, 0}, {0, 1, 0}, {0, 0, -2}}
};

defaultOperatorBasis[n_Integer?Positive] := Module[
  {pauli, gm, diagonal, symmetric, antisymmetric},
  pauli = Table[PauliMatrix[i], {i, 0, 3}];
  gm = legacyGellMannBasis[];
  Which[
    n == 1,
      {{{1}}},
    n == 2,
      pauli,
    n == 3,
      gm,
    n == 4,
      Flatten[
        Table[KroneckerProduct[pauli[[i]], pauli[[j]]], {i, 4}, {j, 4}],
        1
      ],
    n == 6,
      Flatten[
        Table[KroneckerProduct[pauli[[i]], gm[[j]]], {i, 4}, {j, 9}],
        1
      ],
    n == 8,
      Flatten[
        Table[
          KroneckerProduct[pauli[[i]], pauli[[j]], pauli[[k]]],
          {i, 4}, {j, 4}, {k, 4}
        ],
        2
      ],
    True,
      diagonal = Table[matrixUnit[n, a, a], {a, n}];
      symmetric = Flatten[
        Table[
          matrixUnit[n, a, b] + matrixUnit[n, b, a],
          {a, 1, n - 1}, {b, a + 1, n}
        ],
        1
      ];
      antisymmetric = Flatten[
        Table[
          -I matrixUnit[n, a, b] + I matrixUnit[n, b, a],
          {a, 1, n - 1}, {b, a + 1, n}
        ],
        1
      ];
      Join[diagonal, symmetric, antisymmetric]
  ]
];

buildBasisData[basis_List, kind_] := Module[
  {n, count, rows, columns, gram, dual, conjugation, rank},
  If[basis === {} || Not[And @@ (MatrixQ /@ basis)], Return[$Failed]];
  n = Length[First[basis]];
  If[Not[And @@ (Dimensions[#] === {n, n} & /@ basis)], Return[$Failed]];
  count = Length[basis];
  rows = Flatten /@ basis;
  columns = Transpose[rows];
  gram = Simplify[Conjugate[rows] . Transpose[rows]];
  rank = Quiet @ Check[MatrixRank[gram], 0];
  If[rank =!= count,
    Message[CoreAlgorithms::basis];
    Return[$Failed]
  ];
  dual = Quiet @ Check[
    LinearSolve[gram, Conjugate[rows]],
    $Failed
  ];
  If[dual === $Failed,
    Message[CoreAlgorithms::basis];
    Return[$Failed]
  ];
  conjugation = Simplify[dual . Conjugate[columns]];
  <|
    "Dimension" -> n,
    "Count" -> count,
    "Basis" -> basis,
    "Rows" -> rows,
    "Columns" -> columns,
    "Gram" -> gram,
    "Dual" -> dual,
    "Conjugation" -> SparseArray[conjugation],
    "Kind" -> kind,
    "CacheKey" -> ToString[{kind, n, rows}, InputForm]
  |>
];

OperatorBasisData[n_Integer?Positive] := Module[
  {key, cached, value},
  key = ToString[{"Default", n}, InputForm];
  cached = Lookup[$OperatorBasisCache, key, Missing["NotFound"]];
  If[Not[MissingQ[cached]], Return[cached]];
  value = buildBasisData[defaultOperatorBasis[n], "Default"];
  If[value =!= $Failed, AssociateTo[$OperatorBasisCache, key -> value]];
  value
];

OperatorBasisData[basis_List] := buildBasisData[basis, "Custom"];

OperatorBasisData[___] := (
  Message[CoreAlgorithms::dimension];
  $Failed
);

HomogeneousMonomialBasis[order_Integer?NonNegative, dim_Integer?Positive] /;
    1 <= dim <= 3 := Which[
  dim == 3,
    kz^#[[1]] ky^#[[2]] kx^#[[3]] & /@
      Sort @ Flatten[
        Permutations /@ IntegerPartitions[order, {3}, Range[0, order]],
        1
      ],
  dim == 2,
    ky^#[[1]] kx^#[[2]] & /@
      Sort @ Flatten[
        Permutations /@ IntegerPartitions[order, {2}, Range[0, order]],
        1
      ],
  dim == 1,
    {kz^order}
];

HomogeneousMonomialBasis[___] := (
  Message[CoreAlgorithms::dimension];
  $Failed
);

monomialExponents[klist_List, vars_List] :=
  (Exponent[#, #2] & @@@ Thread[{ConstantArray[#, Length[vars]], vars}]) & /@ klist;

(* The explicit form avoids relying on Association keys that are themselves lists. *)
exponentKey[exponent_List] := StringRiffle[ToString[#, InputForm] & /@ exponent, ","];

linearTransformationMatrix[kImage_List, vars_List] := Table[
  Coefficient[Expand[kImage[[i]]], vars[[j]]],
  {i, Length[vars]}, {j, Length[vars]}
];

linearImageQ[kImage_List, vars_List, a_] := Module[
  {residual},
  residual = Expand[kImage - a . vars];
  And @@ (zeroScalarQ /@ residual)
];

monomialLinearData[linearForms_List, vars_List] := Module[
  {a, positions, scales},
  a = Table[
    Coefficient[Expand[linearForms[[i]]], vars[[j]]],
    {i, Length[vars]}, {j, Length[vars]}
  ];
  positions = (Flatten @ Position[#, value_ /; Not[zeroScalarQ[value]]] &) /@ a;
  If[Not[And @@ (Length[#] == 1 & /@ positions)], Return[$Failed]];
  positions = First /@ positions;
  If[DuplicateFreeQ[positions] =!= True, Return[$Failed]];
  scales = MapThread[a[[#1, #2]] &, {Range[Length[vars]], positions}];
  <|"Permutation" -> positions, "Scales" -> scales|>
];

monomialTransformMatrix[klist_List, vars_List, inverseK_List, exponents_List] := Module[
  {data, position, entries, beta, coefficient, row, i},
  data = monomialLinearData[inverseK, vars];
  If[data === $Failed, Return[$Failed]];
  position = AssociationThread[exponentKey /@ exponents, Range[Length[exponents]]];
  entries = Reap[
    Do[
      beta = ConstantArray[0, Length[vars]];
      coefficient = 1;
      Do[
        beta[[data["Permutation"][[i]]]] += exponents[[row, i]];
        coefficient *= data["Scales"][[i]]^exponents[[row, i]],
        {i, Length[vars]}
      ];
      Sow[
        {row, Lookup[position, exponentKey[beta]]} -> coefficient
      ],
      {row, Length[klist]}
    ]
  ][[2, 1]];
  Transpose @ SparseArray[entries, {Length[klist], Length[klist]}]
];

MomentumTransformMatrix[
  order_Integer?NonNegative,
  kImage_List,
  vars_List
] := MomentumTransformMatrix[
  order,
  HomogeneousMonomialBasis[order, Length[vars]],
  kImage,
  vars
];

MomentumTransformMatrix[
  order_Integer?NonNegative,
  klist_List,
  kImage_List,
  vars_List
] := Module[
  {
    key, cached, a, inverseK, exponents, fast, n, position,
    polynomials, coefficientRules, entries, tInverse, row, rule, column
  },
  If[
    Not[1 <= Length[vars] <= 3] ||
    Length[kImage] =!= Length[vars] ||
    Length[klist] == 0,
    Message[CoreAlgorithms::dimension];
    Return[$Failed]
  ];
  key = ToString[{order, klist, kImage, vars}, InputForm];
  cached = Lookup[$MomentumTransformCache, key, Missing["NotFound"]];
  If[Not[MissingQ[cached]], Return[cached]];

  a = linearTransformationMatrix[kImage, vars];
  If[Not[linearImageQ[kImage, vars, a]],
    Message[CoreAlgorithms::kimage, kImage, vars];
    Return[$Failed]
  ];
  inverseK = Quiet @ Check[LinearSolve[a, vars], $Failed];
  If[inverseK === $Failed,
    Message[CoreAlgorithms::kimage, kImage, vars];
    Return[$Failed]
  ];

  n = Length[klist];
  If[zeroMatrixQ[a - IdentityMatrix[Length[vars]]],
    fast = SparseArray[IdentityMatrix[n]];
    AssociateTo[$MomentumTransformCache, key -> fast];
    Return[fast]
  ];

  exponents = (Exponent[#, #2] & @@@
      Thread[{ConstantArray[#, Length[vars]], vars}]) & /@ klist;

  fast = monomialTransformMatrix[klist, vars, inverseK, exponents];
  If[fast =!= $Failed,
    AssociateTo[$MomentumTransformCache, key -> fast];
    Return[fast]
  ];

  position = AssociationThread[exponentKey /@ exponents, Range[n]];
  polynomials = Expand[
    Times @@ MapThread[Power, {inverseK, #}]
  ] & /@ exponents;
  coefficientRules = CoefficientRules[#, vars] & /@ polynomials;
  entries = Reap[
    Do[
      Do[
        column = Lookup[position, exponentKey[First[rule]], Missing["NotFound"]];
        If[Not[MissingQ[column]], Sow[{row, column} -> Last[rule]]],
        {rule, coefficientRules[[row]]}
      ],
      {row, n}
    ]
  ][[2]];
  entries = If[entries === {}, {}, First[entries]];
  tInverse = SparseArray[entries, {n, n}];
  fast = Transpose[tInverse];
  AssociateTo[$MomentumTransformCache, key -> fast];
  fast
];

MomentumTransformMatrix[___] := (
  Message[CoreAlgorithms::dimension];
  $Failed
);

AdjointRepresentationMatrix[
  d_?MatrixQ,
  basisData_Association,
  antiunitary_: False
] := Module[
  {
    key, cached, dAdjoint, dInverse, unitaryQ, flatTransform,
    imageColumns, result
  },
  key = ToString[
    {
      basisData["CacheKey"],
      TrueQ[antiunitary],
      Dimensions[d],
      Normal[d]
    },
    InputForm
  ];
  cached = Lookup[$AdjointRepresentationCache, key, Missing["NotFound"]];
  If[Not[MissingQ[cached]], Return[cached]];

  dAdjoint = ConjugateTranspose[d];
  unitaryQ = zeroMatrixQ[
    Simplify[
      d . dAdjoint - IdentityMatrix[Length[d]]
    ]
  ];
  flatTransform = If[
    unitaryQ,
    KroneckerProduct[d, Conjugate[d]],
    dInverse = Quiet @ Check[Inverse[d], $Failed];
    If[dInverse === $Failed, Return[$Failed]];
    KroneckerProduct[d, Transpose[dInverse]]
  ];
  imageColumns = flatTransform . basisData["Columns"];
  result = basisData["Dual"] . imageColumns;
  result = If[
    TrueQ[antiunitary],
    result . basisData["Conjugation"],
    result
  ];
  result = Map[Simplify, result, {2}];
  result = If[inexactMatrixQ[result], Chop[result], result];
  AssociateTo[$AdjointRepresentationCache, key -> result];
  result
];

AdjointRepresentationMatrix[___] := $Failed;

associationOrEmpty[value_] := If[AssociationQ[value], value, <||>];

lookupSymmetryAssociation[input_Association, names_List] := Fold[
  Join,
  <||>,
  associationOrEmpty[Lookup[input, #, <||>]] & /@ names
];

rawSymmetryCategories[input_Association] := {
  lookupSymmetryAssociation[input, {"Unitary"}],
  lookupSymmetryAssociation[input, {"Anitunitary", "Antiunitary"}],
  lookupSymmetryAssociation[input, {"AntisymmetryUnitaryTest"}],
  lookupSymmetryAssociation[
    input,
    {"AntisymmetryAnitunitaryTest", "AntisymmetryAntiunitaryTest"}
  ]
};

rawSymmetryCategories[input_List] := PadRight[Take[input, UpTo[4]], 4, {}];

continuousKeyName[key_] := ToLowerCase @ Which[
  StringQ[key], key,
  Head[key] === Symbol, SymbolName[key],
  True, ToString[key, InputForm]
];

continuousLabelQ[key_] := MemberQ[
  {"cinfinity", "cinfty"},
  continuousKeyName[key]
];

continuousRules[association_Association] := Select[
  Normal[association],
  continuousLabelQ[First[#]] &
];

nonSystemSymbols[expression_] := DeleteDuplicates @ Cases[
  expression,
  symbol_Symbol /; Context[symbol] =!= "System`",
  Infinity,
  Heads -> True
];

continuousFailure[reason_String] := (
  Message[CoreAlgorithms::continuous, reason];
  $Failed
);

continuousGeneratorData[d_?MatrixQ, parameter_Symbol] := Module[
  {
    n, assumptions, d0, scalar, inverse0, derivative0,
    rawGenerator, generator, logarithmicDerivative
  },
  n = Length[d];
  assumptions = Element[parameter, Reals];
  d0 = FullSimplify[d /. parameter -> 0, assumptions];
  scalar = FullSimplify[Tr[d0]/n, assumptions];
  If[
    zeroScalarQ[scalar] ||
      Not[zeroMatrixQ[FullSimplify[
        d0 - scalar IdentityMatrix[n],
        assumptions
      ]]],
    Return @ continuousFailure[
      "D(0) must be a nonzero scalar multiple of the identity"
    ]
  ];
  If[
    Not[zeroMatrixQ[FullSimplify[
      ConjugateTranspose[d] . d - IdentityMatrix[n],
      assumptions
    ]]],
    Return @ continuousFailure[
      "D(theta) must be unitary for real theta"
    ]
  ];
  inverse0 = Quiet @ Check[Inverse[d0], $Failed];
  If[inverse0 === $Failed,
    Return @ continuousFailure["D(0) must be invertible"]
  ];
  derivative0 = FullSimplify[
    D[d, parameter] /. parameter -> 0,
    assumptions
  ];
  rawGenerator = FullSimplify[
    I derivative0 . inverse0,
    assumptions
  ];
  logarithmicDerivative = Quiet @ Check[
    FullSimplify[
      I D[d, parameter] . Inverse[d],
      assumptions
    ],
    $Failed
  ];
  If[
    logarithmicDerivative === $Failed ||
      Not[zeroMatrixQ[FullSimplify[
        logarithmicDerivative - rawGenerator,
        assumptions
      ]]],
    Return @ continuousFailure[
      "D(theta) must be a one-parameter representation with a theta-independent generator"
    ]
  ];
  generator = FullSimplify[
    rawGenerator - Tr[rawGenerator]/n IdentityMatrix[n],
    assumptions
  ];
  If[
    Not[FreeQ[generator, parameter]] ||
      Not[zeroMatrixQ[FullSimplify[
        generator - ConjugateTranspose[generator],
        assumptions
      ]]],
    Return @ continuousFailure[
      "the extracted traceless generator must be theta-independent and Hermitian"
    ]
  ];
  <|
    "RawGenerator" -> rawGenerator,
    "Generator" -> generator
  |>
];

continuousInvariantBasisData[
  generator_?MatrixQ,
  basisData_Association
] := Module[
  {
    n, identity, flatGenerator, action, kernel,
    reducedBasis, reducedData
  },
  n = basisData["Dimension"];
  identity = IdentityMatrix[n];
  flatGenerator = I (
    KroneckerProduct[generator, identity] -
      KroneckerProduct[identity, Transpose[generator]]
  );
  action = Map[
    Simplify,
    basisData["Dual"] . flatGenerator . basisData["Columns"],
    {2}
  ];
  kernel = Quiet @ Check[NullSpace[action], $Failed];
  If[kernel === $Failed || kernel === {},
    Return @ continuousFailure[
      "the commutant of the extracted generator could not be constructed"
    ]
  ];
  reducedBasis = Map[
    Function[coordinates,
      Map[
        Simplify,
        Total @ MapThread[
          Times,
          {coordinates, basisData["Basis"]}
        ],
        {2}
      ]
    ],
    kernel
  ];
  reducedData = buildBasisData[reducedBasis, "CInfinity"];
  If[reducedData === $Failed,
    Return @ continuousFailure[
      "the reduced Hermitian operator basis is singular"
    ]
  ];
  <|
    "BasisData" -> reducedData,
    "GeneratorAction" -> action,
    "Kernel" -> kernel
  |>
];

ContinuousNormalizerQ[
  d_?MatrixQ,
  generator_?MatrixQ,
  antiunitary_: False
] := Module[{inverse, transformed},
  If[Dimensions[d] =!= Dimensions[generator], Return[False]];
  inverse = Quiet @ Check[Inverse[d], $Failed];
  If[inverse === $Failed, Return[False]];
  transformed = If[
    TrueQ[antiunitary],
    -d . Conjugate[generator] . inverse,
    d . generator . inverse
  ];
  zeroMatrixQ[Simplify[transformed - generator]] ||
    zeroMatrixQ[Simplify[transformed + generator]]
];

ContinuousNormalizerQ[___] := False;

PrepareContinuousSymmetry[
  input_Association,
  basisData_Association,
  vars_List
] := Module[
  {
    categories, unitaryRules, forbiddenRules, rule, key, value,
    d, kImage, parameter, extraSymbols, generatorData, reduced,
    discreteInput, originalUnitary
  },
  categories = rawSymmetryCategories[input];
  unitaryRules = continuousRules[categories[[1]]];
  forbiddenRules = Flatten[
    continuousRules /@ Select[Rest[categories], AssociationQ]
  ];
  If[forbiddenRules =!= {},
    Return @ continuousFailure[
      "CInfinity is allowed only in the Unitary category"
    ]
  ];
  If[unitaryRules === {},
    Return @ <|
      "Present" -> False,
      "BasisData" -> basisData,
      "DiscreteInput" -> input,
      "OriginalOperatorDimension" -> basisData["Count"],
      "ReducedOperatorDimension" -> basisData["Count"]
    |>
  ];
  If[Length[unitaryRules] =!= 1,
    Return @ continuousFailure[
      "exactly one CInfinity entry is supported"
    ]
  ];
  rule = First[unitaryRules];
  key = First[rule];
  value = Last[rule];
  If[Not[ListQ[value]] || Length[value] =!= 3,
    Return @ continuousFailure[
      "the entry must have the form {D(theta), {kx,ky,kz}, theta}"
    ]
  ];
  {d, kImage, parameter} = value;
  If[
    Not[MatrixQ[d]] ||
      Dimensions[d] =!= {
        basisData["Dimension"],
        basisData["Dimension"]
      },
    Return @ continuousFailure[
      "D(theta) must have the same square dimension as the discrete representations"
    ]
  ];
  If[Head[parameter] =!= Symbol,
    Return @ continuousFailure[
      "the third entry must be the continuous-parameter symbol"
    ]
  ];
  If[
    Not[ListQ[kImage]] || Length[kImage] =!= Length[vars] ||
      Not[zeroMatrixQ[Simplify[kImage - vars]]],
    Return @ continuousFailure[
      "the internal CInfinity symmetry must leave momentum unchanged"
    ]
  ];
  extraSymbols = DeleteCases[
    nonSystemSymbols[d],
    symbol_ /; SameQ[symbol, parameter]
  ];
  If[extraSymbols =!= {},
    Return @ continuousFailure[
      "D(theta) may contain no user symbols other than the declared continuous parameter"
    ]
  ];
  generatorData = continuousGeneratorData[d, parameter];
  If[generatorData === $Failed, Return[$Failed]];
  reduced = continuousInvariantBasisData[
    generatorData["Generator"],
    basisData
  ];
  If[reduced === $Failed, Return[$Failed]];
  originalUnitary = Lookup[input, "Unitary", <||>];
  discreteInput = input;
  AssociateTo[
    discreteInput,
    "Unitary" -> If[
      AssociationQ[originalUnitary],
      KeyDrop[originalUnitary, {key}],
      <||>
    ]
  ];
  <|
    "Present" -> True,
    "Type" -> "CInfinity",
    "Key" -> key,
    "Parameter" -> parameter,
    "Representation" -> d,
    "RawGenerator" -> generatorData["RawGenerator"],
    "Generator" -> generatorData["Generator"],
    "GeneratorAction" -> reduced["GeneratorAction"],
    "Kernel" -> reduced["Kernel"],
    "BasisData" -> reduced["BasisData"],
    "DiscreteInput" -> discreteInput,
    "OriginalOperatorDimension" -> basisData["Count"],
    "ReducedOperatorDimension" -> reduced["BasisData", "Count"]
  |>
];

PrepareContinuousSymmetry[input_, basisData_Association, _List] := <|
  "Present" -> False,
  "BasisData" -> basisData,
  "DiscreteInput" -> input,
  "OriginalOperatorDimension" -> basisData["Count"],
  "ReducedOperatorDimension" -> basisData["Count"]
|>;

PrepareSymmetryRecords[
  input_,
  basisData_Association,
  allowEmpty_: False,
  continuousGenerator_: None
] := Module[
  {
    categories, records, values, type, anti, eta, value,
    d, kImage, j, failureTag
  },
  categories = rawSymmetryCategories[input];
  failureTag = Unique["PrepareSymmetryRecordsFailure"];
  records = Catch[
    Reap[
    Do[
      values = If[AssociationQ[categories[[type]]],
        Values[categories[[type]]],
        categories[[type]]
      ];
      anti = MemberQ[{2, 4}, type];
      eta = If[MemberQ[{1, 2}, type], 1, -1];
      Do[
        If[Not[ListQ[value]] || Length[value] < 2, Continue[]];
        d = value[[1]];
        kImage = value[[2]];
        If[
          continuousGenerator =!= None && MatrixQ[d] &&
            Not[ContinuousNormalizerQ[d, continuousGenerator, anti]],
          Message[CoreAlgorithms::normalizer, d];
          Throw[$Failed, failureTag]
        ];
        j = AdjointRepresentationMatrix[d, basisData, anti];
        If[j === $Failed, Continue[]];
        Sow @ <|
          "D" -> d,
          "KImage" -> kImage,
          "Type" -> type,
          "Antiunitary" -> anti,
          "Eta" -> eta,
          "J" -> j
        |>,
        {value, values}
      ],
      {type, 4}
    ]
    ],
    failureTag
  ];
  If[records === $Failed, Return[$Failed]];
  records = records[[2]];
  records = If[records === {}, {}, First[records]];
  If[records === {} && Not[TrueQ[allowEmpty]],
    Message[CoreAlgorithms::symmetry]
  ];
  records
];

factorRepresentation[factor_Association] :=
  factor["Eta"] KroneckerProduct[factor["F"], factor["J"]];

IterativeInvariantBasis[factors_List] := Module[
  {ambient, identity, generators},
  If[factors === {}, Return[{}]];
  ambient = Length[factors[[1, "F"]]] Length[factors[[1, "J"]]];
  identity = IdentityMatrix[ambient, SparseArray];
  generators = factorRepresentation /@ factors;
  Fold[
    Function[{basis, representation},
      If[
        basis === {},
        {},
        Module[{kernel},
          kernel = Quiet @ Check[
            NullSpace[(representation - identity) . Transpose[basis]],
            {}
          ];
          If[kernel === {}, {}, Simplify[kernel . basis]]
        ]
      ]
    ],
    IdentityMatrix[ambient],
    generators
  ]
];

StackedInvariantBasis[factors_List] := Module[
  {ambient, identity, constraints},
  If[factors === {}, Return[{}]];
  ambient = Length[factors[[1, "F"]]] Length[factors[[1, "J"]]];
  identity = IdentityMatrix[ambient, SparseArray];
  constraints = (factorRepresentation[#] - identity) & /@ factors;
  Quiet @ Check[NullSpace[Join @@ constraints], {}]
];

factorIdentity[factor_Association] := <|
  "F" -> IdentityMatrix[Length[factor["F"]], SparseArray],
  "J" -> IdentityMatrix[Length[factor["J"]], SparseArray],
  "Eta" -> 1
|>;

factorProduct[left_Association, right_Association] := <|
  "F" -> SparseArray[Simplify[left["F"] . right["F"]]],
  "J" -> SparseArray[Simplify[left["J"] . right["J"]]],
  "Eta" -> Simplify[left["Eta"] right["Eta"]]
|>;

factorKey[factor_Association] := StringRiffle[
  {
    matrixCacheKey[factor["F"]],
    matrixCacheKey[factor["J"]],
    ToString[canonicalScalar[factor["Eta"]], InputForm]
  },
  "::"
];

inexactFactorQ[factor_Association] :=
  inexactMatrixQ[factor["F"]] ||
  inexactMatrixQ[factor["J"]] ||
  InexactNumberQ[factor["Eta"]];

factorGroupClosure[generators_List, maxOrder_Integer?Positive] := Module[
  {identity, seen, queue, head, current, next, key, generator},
  If[generators === {}, Return[$Failed]];
  If[AnyTrue[generators, inexactFactorQ], Return[$Failed]];
  identity = factorIdentity[First[generators]];
  seen = <|factorKey[identity] -> identity|>;
  queue = {identity};
  head = 1;
  While[head <= Length[queue],
    current = queue[[head]];
    head++;
    Do[
      next = factorProduct[current, generator];
      key = factorKey[next];
      If[Not[KeyExistsQ[seen, key]],
        If[Length[queue] >= maxOrder, Return[$Failed]];
        AssociateTo[seen, key -> next];
        AppendTo[queue, next]
      ],
      {generator, generators}
    ]
  ];
  Values[seen]
];

reynoldsProjector[factors_List] := Module[
  {ambient, sum},
  ambient = Length[factors[[1, "F"]]] Length[factors[[1, "J"]]];
  sum = SparseArray[{}, {ambient, ambient}];
  Do[
    sum = sum + factor["Eta"] KroneckerProduct[factor["F"], factor["J"]],
    {factor, factors}
  ];
  Simplify[sum/Length[factors]]
];

verifyReynoldsProjector[projector_, generators_List] := Module[
  {idempotent, absorbed},
  idempotent = zeroMatrixQ[Simplify[projector . projector - projector]];
  absorbed = And @@ (
    zeroMatrixQ[Simplify[factorRepresentation[#] . projector - projector]] & /@
      generators
  );
  TrueQ[idempotent && absorbed]
];

Options[ReynoldsInvariantBasis] = {
  "GroupInput" -> "Generators",
  "MaxGroupOrder" -> 512,
  "Verify" -> True
};

ReynoldsInvariantBasis[factors_List, OptionsPattern[]] := Module[
  {
    groupInput, maxOrder, verify, elements, projector, basis,
    success, reason = None, fallback = False, completedElements,
    elementClosureCompleted = False
  },
  If[factors === {},
    Return @ <|
      "Basis" -> {},
      "Success" -> True,
      "FallbackUsed" -> False,
      "ElementClosureCompleted" -> False,
      "GroupOrder" -> 0,
      "Projector" -> {}
    |>
  ];
  groupInput = OptionValue["GroupInput"];
  maxOrder = OptionValue["MaxGroupOrder"];
  verify = TrueQ[OptionValue["Verify"]];
  elements = Which[
    groupInput === "Elements",
      If[AnyTrue[factors, inexactFactorQ], $Failed, factors],
    groupInput === "Generators",
      factorGroupClosure[factors, maxOrder],
    True,
      factorGroupClosure[factors, maxOrder]
  ];
  If[elements === $Failed,
    reason = "the exact finite factor closure could not be established within the group-order limit";
    fallback = True;
    Message[ReynoldsInvariantBasis::fallback, reason];
    Return @ <|
      "Basis" -> IterativeInvariantBasis[factors],
      "Success" -> False,
      "FallbackUsed" -> True,
      "ElementClosureCompleted" -> False,
      "Reason" -> reason,
      "GroupOrder" -> Missing["Unknown"],
      "Projector" -> Missing["NotConstructed"]
    |>
  ];
  projector = reynoldsProjector[elements];
  success = Not[verify] || verifyReynoldsProjector[projector, factors];
  If[Not[success] && groupInput === "Elements",
    (* Some projective or mixed single/double-valued inputs list one
       representative per spatial operation.  Their induced factors can
       close to a finite central extension.  Complete that small extension
       only when the nominal element average is not already a projector. *)
    completedElements = factorGroupClosure[factors, maxOrder];
    If[completedElements =!= $Failed,
      projector = reynoldsProjector[completedElements];
      success = Not[verify] ||
        verifyReynoldsProjector[projector, factors];
      If[success,
        elements = completedElements;
        elementClosureCompleted = True
      ]
    ]
  ];
  If[Not[success],
    reason = "the averaged operator failed idempotence or generator-absorption verification";
    fallback = True;
    Message[ReynoldsInvariantBasis::fallback, reason];
    basis = IterativeInvariantBasis[factors],
    basis = nonzeroRows[RowReduce[Transpose[projector]]]
  ];
  <|
    "Basis" -> basis,
    "Success" -> success,
    "FallbackUsed" -> fallback,
    "ElementClosureCompleted" -> elementClosureCompleted,
    "Reason" -> reason,
    "GroupOrder" -> Length[elements],
    "Projector" -> If[fallback, Missing["VerificationFailed"], projector]
  |>
];

End[]

EndPackage[]
