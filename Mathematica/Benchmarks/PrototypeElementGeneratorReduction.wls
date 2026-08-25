$HistoryLength = 0;

benchmarkRoot = ExpandFileName[DirectoryName[$InputFileName]];
mathematicaRoot = ParentDirectory[benchmarkRoot];
arguments = Rest[$ScriptCommandLine];
msg = If[
  Length[arguments] >= 2,
  ToExpression /@ Take[arguments, 2],
  {221, 93}
];
repIndex = If[
  Length[arguments] >= 3,
  ToExpression[arguments[[3]]],
  16
];
order = If[
  Length[arguments] >= 4,
  ToExpression[arguments[[4]]],
  1
];

Needs["MSGCorep`"];
Needs["SpaceGroupIrep`"];
Get[FileNameJoin[{
  mathematicaRoot,
  "MagneticKP",
  "Kernel",
  "init.m"
}]];

ClearAll[normalMatrix];
normalMatrix[value_] := If[MatrixQ[value], value, {{value}}];

ClearAll[gammaData];
gammaData[group_List] := Module[{name},
  name = SelectFirst[
    Keys[MSGCorep`libMLGCorep[group]],
    Function[kname,
      MemberQ[
        Keys[MSGCorep`libMLGCorep[group, kname]],
        {0, 0, 0}
      ]
    ],
    Missing["NotFound"]
  ];
  If[
    MissingQ[name],
    $Failed,
    MSGCorep`lookupMLGCorep[group, name, {0, 0, 0}]
  ]
];

ClearAll[cartesianRotations];
cartesianRotations[group_List, mlg_List] := Module[
  {lattice, basic, toCartesian},
  lattice = SpaceGroupIrep`getSGLatt[First[group]];
  basic = SpaceGroupIrep`BasicVectors[lattice];
  toCartesian[matrix_] := FullSimplify[
    Inverse[basic] . matrix . basic
  ];
  toCartesian[
    SpaceGroupIrep`getRotMatOfK[lattice, #[[1]]]
  ] & /@ mlg
];

ClearAll[buildInput];
buildInput[mlg_List, rotations_List, corep_List] := Module[
  {vars, unitaryIndices, antiunitaryIndices, makeEntry},
  vars = {
    MagneticKP`kx,
    MagneticKP`ky,
    MagneticKP`kz
  };
  unitaryIndices = Flatten @ Position[Last /@ mlg, 0];
  antiunitaryIndices = Flatten @ Position[Last /@ mlg, 1];
  makeEntry[index_Integer, antiunitary_] := {
    normalMatrix[corep[[index]]],
    Expand[
      rotations[[index]] .
        If[TrueQ[antiunitary], -vars, vars]
    ]
  };
  <|
    "Unitary" -> Association @ Table[
      ToString[index] -> makeEntry[index, False],
      {index, unitaryIndices}
    ],
    "Anitunitary" -> Association @ Table[
      ToString[index] -> makeEntry[index, True],
      {index, antiunitaryIndices}
    ],
    "GroupInput" -> "Elements"
  |>
];

ClearAll[greedyFactorGenerators];
greedyFactorGenerators[factors_List, maxOrder_Integer?Positive] := Module[
  {generators = {}, closure = {}, closureKeys = <||>, key},
  Do[
    key = MagneticKP`Private`factorKey[factor];
    If[Not[KeyExistsQ[closureKeys, key]],
      AppendTo[generators, factor];
      closure = MagneticKP`Private`factorGroupClosure[
        generators,
        maxOrder
      ];
      If[closure === $Failed, Return[factors]];
      closureKeys = AssociationThread[
        MagneticKP`Private`factorKey /@ closure,
        ConstantArray[True, Length[closure]]
      ]
    ],
    {factor, factors}
  ];
  generators
];

ClearAll[vectorizedUnitaryAdjoint];
vectorizedUnitaryAdjoint[
  d_?MatrixQ,
  basisData_Association,
  antiunitary_,
  canonicalizer_
] := Module[{imageColumns, unitaryResult, result},
  imageColumns = KroneckerProduct[d, Conjugate[d]] .
    basisData["Columns"];
  unitaryResult = basisData["Dual"] . imageColumns;
  result = If[
    TrueQ[antiunitary],
    unitaryResult . basisData["Conjugation"],
    unitaryResult
  ];
  Map[canonicalizer, result, {2}]
];

data = gammaData[msg];
If[data === $Failed,
  Print["No Gamma data for ", InputForm[msg]];
  Exit[2]
];
mlg = data["MLG"];
rotations = cartesianRotations[msg, mlg];
coreps = Join @@ data["corep"];
If[repIndex < 1 || repIndex > Length[coreps],
  Print["Rep index out of range: ", repIndex];
  Exit[2]
];
input = buildInput[mlg, rotations, coreps[[repIndex]]];
dimension = Length @ normalMatrix[coreps[[repIndex, 1]]];
basisTiming = AbsoluteTiming[
  basisData = MagneticKP`Private`OperatorBasisData[dimension];
][[1]];
coldDPairs = Join[
  {#[[1]], False} & /@ Values[input["Unitary"]],
  {#[[1]], True} & /@ Values[input["Anitunitary"]]
];
coldUnitaritySimplifyTiming = AbsoluteTiming[
  coldUnitarityChecks = MagneticKP`Private`zeroMatrixQ[
      Simplify[
        #[[1]] . ConjugateTranspose[#[[1]]] -
          IdentityMatrix[dimension]
      ]
    ] & /@ coldDPairs;
][[1]];
coldVectorizedSimplifyTiming = AbsoluteTiming[
  coldVectorizedSimplified = vectorizedUnitaryAdjoint[
      #[[1]],
      basisData,
      #[[2]],
      Simplify
    ] & /@ coldDPairs;
][[1]];
recordsTiming = AbsoluteTiming[
  records = MagneticKP`Private`PrepareSymmetryRecords[
    input,
    basisData
  ];
][[1]];
klist = MagneticKP`Private`HomogeneousMonomialBasis[order, 3];
factorTiming = AbsoluteTiming[
  factors = Map[
    Function[record,
      <|
        "F" -> MagneticKP`Private`MomentumTransformMatrix[
          order,
          klist,
          record["KImage"],
          {
            MagneticKP`kx,
            MagneticKP`ky,
            MagneticKP`kz
          }
        ],
        "J" -> record["J"],
        "Eta" -> record["Eta"]
      |>
    ],
    records
  ];
][[1]];
kpTiming = AbsoluteTiming[
  kpResult = MagneticKP`kpHam[
    order,
    input,
    "Method" -> "IterativeSimplification"
  ];
][[1]];
vectorizedRawTiming = AbsoluteTiming[
  vectorizedRaw = vectorizedUnitaryAdjoint[
      #["D"],
      basisData,
      #["Antiunitary"],
      Identity
    ] & /@ records;
][[1]];
vectorizedRootReduceTiming = AbsoluteTiming[
  vectorizedRootReduced = vectorizedUnitaryAdjoint[
      #["D"],
      basisData,
      #["Antiunitary"],
      RootReduce
    ] & /@ records;
][[1]];
vectorizedSimplifyTiming = AbsoluteTiming[
  vectorizedSimplified = vectorizedUnitaryAdjoint[
      #["D"],
      basisData,
      #["Antiunitary"],
      Simplify
    ] & /@ records;
][[1]];
unitarityCheckTiming = AbsoluteTiming[
  unitarityChecks = MagneticKP`Private`zeroMatrixQ[
      Simplify[
        #["D"] . ConjugateTranspose[#["D"]] -
          IdentityMatrix[dimension]
      ]
    ] & /@ records;
][[1]];
structuralKeyTiming = AbsoluteTiming[
  structuralKeys = ToString[
      {
        basisData["CacheKey"],
        Dimensions[#["D"]],
        Normal[#["D"]]
      },
      InputForm
    ] & /@ records;
][[1]];
hashKeyTiming = AbsoluteTiming[
  hashKeys = Hash[
      HoldComplete[
        basisData["CacheKey"],
        Normal[#["D"]]
      ],
      "SHA256"
    ] & /@ records;
][[1]];
vectorizedRawPass = And @@ MapThread[
  MagneticKP`Private`zeroMatrixQ[
    Simplify[#1 - #2]
  ] &,
  {vectorizedRaw, Lookup[records, "J"]}
];
vectorizedRootReducePass = And @@ MapThread[
  MagneticKP`Private`zeroMatrixQ[
    Simplify[#1 - #2]
  ] &,
  {vectorizedRootReduced, Lookup[records, "J"]}
];
vectorizedSimplifyPass = And @@ MapThread[
  MagneticKP`Private`zeroMatrixQ[
    Simplify[#1 - #2]
  ] &,
  {vectorizedSimplified, Lookup[records, "J"]}
];

fullTiming = AbsoluteTiming[
  fullBasis = MagneticKP`Private`IterativeInvariantBasis[factors];
][[1]];
fullSimplifyTiming = AbsoluteTiming[
  simplifiedFullBasis = Simplify[fullBasis];
][[1]];
fullHamiltonianTiming = AbsoluteTiming[
  fullHamiltonian = MagneticKP`Private`hamiltonianFromRowBasis[
    simplifiedFullBasis,
    order,
    klist,
    basisData["Basis"],
    dimension
  ];
][[1]];
reductionTiming = AbsoluteTiming[
  reducedFactors = greedyFactorGenerators[factors, 512];
][[1]];
reducedTiming = AbsoluteTiming[
  reducedBasis = MagneticKP`Private`IterativeInvariantBasis[
    reducedFactors
  ];
][[1]];

sameSpace = TrueQ[
  MatrixRank[fullBasis] ==
    MatrixRank[reducedBasis] ==
    MatrixRank[Join[fullBasis, reducedBasis]]
];

Print @ Dataset @ <|
  "MSG" -> msg,
  "RepIndex" -> repIndex,
  "Dimension" -> dimension,
  "Order" -> order,
  "FullKPSeconds" -> kpTiming,
  "FullKPParameters" ->
    Lookup[kpResult, "NumberOfParameters", Missing["NotFound"]],
  "ElementCount" -> Length[factors],
  "GeneratorCount" -> Length[reducedFactors],
  "PrepareBasisSeconds" -> basisTiming,
  "ColdUnitaritySimplifySeconds" ->
    coldUnitaritySimplifyTiming,
  "ColdAllUnitary" -> And @@ coldUnitarityChecks,
  "ColdVectorizedSimplifySeconds" ->
    coldVectorizedSimplifyTiming,
  "PrepareRecordsSeconds" -> recordsTiming,
  "VectorizedRawSeconds" -> vectorizedRawTiming,
  "VectorizedRootReduceSeconds" -> vectorizedRootReduceTiming,
  "VectorizedSimplifySeconds" -> vectorizedSimplifyTiming,
  "UnitarityCheckSeconds" -> unitarityCheckTiming,
  "AllUnitary" -> And @@ unitarityChecks,
  "StructuralKeySeconds" -> structuralKeyTiming,
  "HashKeySeconds" -> hashKeyTiming,
  "StructuralKeyCount" -> Length[DeleteDuplicates[structuralKeys]],
  "HashKeyCount" -> Length[DeleteDuplicates[hashKeys]],
  "VectorizedRawPass" -> vectorizedRawPass,
  "VectorizedRootReducePass" -> vectorizedRootReducePass,
  "VectorizedSimplifyPass" -> vectorizedSimplifyPass,
  "PrepareFactorsSeconds" -> factorTiming,
  "FullIterativeSeconds" -> fullTiming,
  "FullSimplifySeconds" -> fullSimplifyTiming,
  "FullHamiltonianSeconds" -> fullHamiltonianTiming,
  "GeneratorReductionSeconds" -> reductionTiming,
  "ReducedIterativeSeconds" -> reducedTiming,
  "CombinedReducedSeconds" -> reductionTiming + reducedTiming,
  "SolverSpeedup" -> fullTiming/reducedTiming,
  "CombinedSpeedup" ->
    fullTiming/(reductionTiming + reducedTiming),
  "InvariantDimension" -> Length[fullBasis],
  "SameSpace" -> sameSpace
|>;

Exit[If[sameSpace, 0, 1]]
