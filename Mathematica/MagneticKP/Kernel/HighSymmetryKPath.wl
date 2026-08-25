(* ::Package:: *)

BeginPackage["MagneticKP`"]

Begin["`Private`"]

ClearAll[
  kPathActionMatrix,
  kPathCanonicalMatrix,
  kPathCanonicalElement,
  kPathCollinearQ,
  kPathElementAction,
  kPathElementData,
  kPathElementEqualQ,
  kPathElementIndex,
  kPathElementMultiply,
  kPathEnumerateSubgroups,
  kPathEquivalentLineQ,
  kPathFixedBasis,
  kPathGeneratedSubgroup,
  kPathGroupClosure,
  kPathLineData,
  kPathNormalizeDirection,
  kPathProductTable,
  kPathSymmetryGenerators,
  kPathZeroMatrixQ,
  kPathZeroVectorQ
];

kPathCanonicalMatrix[matrix_?MatrixQ] := Map[
  RootReduce[Together[#]] &,
  matrix,
  {2}
];

kPathZeroVectorQ[vector_List] := And @@ (
  TrueQ[PossibleZeroQ[RootReduce[Together[#]]]] & /@ vector
);

kPathZeroMatrixQ[matrix_?MatrixQ] :=
  kPathZeroVectorQ[Flatten[matrix]];

kPathCanonicalElement[{matrix_?MatrixQ, parity_Integer}] := {
  kPathCanonicalMatrix[matrix],
  Mod[parity, 2]
};

kPathElementEqualQ[
    {firstMatrix_?MatrixQ, firstParity_Integer},
    {secondMatrix_?MatrixQ, secondParity_Integer}
  ] := firstParity === secondParity &&
    kPathZeroMatrixQ[firstMatrix - secondMatrix];

kPathElementIndex[elements_List, element_List] := SelectFirst[
  Range[Length[elements]],
  kPathElementEqualQ[elements[[#]], element] &,
  Missing["NotFound"]
];

kPathElementMultiply[first_List, second_List] := kPathCanonicalElement[{
  first[[1]] . second[[1]],
  Mod[first[[2]] + second[[2]], 2]
}];

kPathElementAction[{spatial_?MatrixQ, parity_Integer}] :=
  (-1)^parity spatial;

kPathElementData[element_List] := <|
  "SpatialMatrix" -> element[[1]],
  "Antiunitary" -> TrueQ[element[[2]] === 1],
  "MomentumActionMatrix" -> kPathElementAction[element]
|>;

kPathSymmetryGenerators[input_Association] := Module[
  {
    unitary,
    antiunitary,
    unitaryValues,
    antiunitaryValues,
    unitaryMatrices,
    antiunitaryMatrices
  },
  unitary = Lookup[input, "Unitary", <||>];
  antiunitary = Lookup[
    input,
    "Antiunitary",
    Lookup[input, "Anitunitary", <||>]
  ];
  unitaryValues = If[
    AssociationQ[unitary],
    Values[unitary],
    If[ListQ[unitary], unitary, {}]
  ];
  antiunitaryValues = If[
    AssociationQ[antiunitary],
    Values[antiunitary],
    If[ListQ[antiunitary], antiunitary, {}]
  ];
  unitaryMatrices = kPathActionMatrix /@ unitaryValues;
  antiunitaryMatrices = kPathActionMatrix /@ antiunitaryValues;
  If[
    unitaryValues === {} && antiunitaryValues === {} ||
      MemberQ[Join[unitaryMatrices, antiunitaryMatrices], $Failed],
    Return[$Failed]
  ];
  Join[
    ({#, 0} & /@ unitaryMatrices),
    ({-#, 1} & /@ antiunitaryMatrices)
  ]
];

kPathActionMatrix[value_] := Module[
  {image, variables = {kx, ky, kz}, matrix, reconstructed},
  If[!ListQ[value] || Length[value] < 2, Return[$Failed]];
  image = value[[2]];
  If[!ListQ[image] || Length[image] =!= 3, Return[$Failed]];
  matrix = Table[
    Coefficient[Expand[image[[row]]], variables[[column]]],
    {row, 3},
    {column, 3}
  ];
  reconstructed = matrix . variables;
  If[
    !kPathZeroVectorQ[Expand[image - reconstructed]],
    Return[$Failed]
  ];
  kPathCanonicalMatrix[matrix]
];

(* This is the same closure strategy used by MagneticTB's abstract-group
   layer: start at the identity, multiply the queue by every generator, and
   record each previously unseen element. *)
kPathGroupClosure[generators_List, maxOrder_Integer?Positive] := Module[
  {
    identity = {IdentityMatrix[3], 0},
    canonicalGenerators,
    elements,
    queue,
    head = 1,
    current,
    product
  },
  canonicalGenerators = DeleteDuplicates[
    kPathCanonicalElement /@ generators,
    kPathElementEqualQ
  ];
  elements = {identity};
  queue = {identity};
  While[head <= Length[queue],
    current = queue[[head++]];
    Do[
      product = kPathElementMultiply[current, generator];
      If[MissingQ[kPathElementIndex[elements, product]],
        If[Length[elements] >= maxOrder,
          Return[Failure[
            "GroupOrderLimit",
            <|"MaxGroupOrder" -> maxOrder|>
          ]]
        ];
        AppendTo[elements, product];
        AppendTo[queue, product]
      ],
      {generator, canonicalGenerators}
    ]
  ];
  elements
];

kPathProductTable[elements_List] := Module[
  {table, index},
  table = Table[
    index = kPathElementIndex[
      elements,
      kPathElementMultiply[elements[[left]], elements[[right]]]
    ];
    If[MissingQ[index], Return[$Failed, Module]];
    index,
    {left, Length[elements]},
    {right, Length[elements]}
  ];
  table
];

(* Ported in structure from MagneticTB`GroupAlgebra`GeneratedSubgroupIndices:
   close a generator-index set under left and right multiplication in the
   verified product table. *)
kPathGeneratedSubgroup[table_List, generators_List] := Module[
  {seen = <|1 -> True|>, queue = {1}, head = 1, current, products},
  While[head <= Length[queue],
    current = queue[[head++]];
    products = Join[
      table[[current, generators]],
      table[[generators, current]]
    ];
    Scan[
      Function[product,
        If[!KeyExistsQ[seen, product],
          AssociateTo[seen, product -> True];
          AppendTo[queue, product]
        ]
      ],
      products
    ]
  ];
  Sort[Keys[seen]]
];

kPathEnumerateSubgroups[table_List, maxSubgroups_Integer?Positive] := Module[
  {
    order = Length[table],
    subgroups = {{1}},
    seen = <|"1" -> True|>,
    head = 1,
    subgroup,
    candidate,
    key
  },
  While[head <= Length[subgroups],
    subgroup = subgroups[[head++]];
    Do[
      If[MemberQ[subgroup, element], Continue[]];
      candidate = kPathGeneratedSubgroup[
        table,
        Append[subgroup, element]
      ];
      key = StringRiffle[ToString /@ candidate, ","];
      If[!KeyExistsQ[seen, key],
        If[Length[subgroups] >= maxSubgroups,
          Return[Failure[
            "SubgroupLimit",
            <|"MaxSubgroups" -> maxSubgroups|>
          ]]
        ];
        AssociateTo[seen, key -> True];
        AppendTo[subgroups, candidate]
      ],
      {element, order}
    ]
  ];
  SortBy[subgroups, {Length[#], #} &]
];

kPathFixedBasis[elements_List, subgroup_List] := NullSpace[
  Join @@ (
    kPathElementAction[elements[[#]]] - IdentityMatrix[3] & /@ subgroup
  )
];

kPathNormalizeDirection[direction_List] := Module[
  {exact, magnitudes, pivot, normalized},
  exact = RootReduce[Together[#]] & /@ direction;
  magnitudes = Abs[N[exact]];
  If[Max[magnitudes] == 0, Return[$Failed]];
  pivot = First@FirstPosition[magnitudes, Max[magnitudes]];
  normalized = RootReduce[Together[#]] & /@ (exact/exact[[pivot]]);
  normalized
];

kPathCollinearQ[first_List, second_List] := TrueQ[
  MatrixRank[kPathCanonicalMatrix[{first, second}]] === 1
];

kPathEquivalentLineQ[
    first_List,
    second_List,
    elements_List
  ] := AnyTrue[
  elements,
  kPathCollinearQ[kPathElementAction[#] . first, second] &
];

kPathLineData[elements_List, subgroups_List] := Module[
  {
    identity = IdentityMatrix[3],
    lines = {},
    basis,
    direction,
    stabilizer,
    stabilizerActions,
    distinctActions
  },
  Do[
    basis = kPathFixedBasis[elements, subgroup];
    If[Length[basis] =!= 1, Continue[]];
    direction = kPathNormalizeDirection[First[basis]];
    If[direction === $Failed, Continue[]];
    stabilizer = Select[
      Range[Length[elements]],
      kPathZeroVectorQ[
        (kPathElementAction[elements[[#]]] - identity) . direction
      ] &
    ];
    stabilizerActions = kPathElementAction /@ elements[[stabilizer]];
    distinctActions = DeleteDuplicates[
      stabilizerActions,
      kPathZeroMatrixQ[#1 - #2] &
    ];
    (* Replacing the candidate subgroup by the complete stabilizer makes this
       the largest subgroup that fixes the corresponding generic line point. *)
    If[
      !AnyTrue[
        lines,
        Function[line,
          kPathCollinearQ[line["Direction"], direction]
        ]
      ],
      AppendTo[
        lines,
        <|
          "Direction" -> direction,
          "SubgroupOrder" -> Length[stabilizer],
          "MomentumActionOrder" -> Length[distinctActions],
          "SubgroupIndices" -> stabilizer,
          "SubgroupMatrices" -> stabilizerActions,
          "SubgroupElements" -> (kPathElementData /@ elements[[stabilizer]])
        |>
      ]
    ],
    {subgroup, subgroups}
  ];
  SortBy[
    lines,
    {
      -Lookup[#, "SubgroupOrder", 0],
      ToString[Lookup[#, "Direction", {}], InputForm]
    } &
  ]
];

highSymmetryKPath::input =
  "The input must be an interfaceRep-style Association containing three-dimensional finite momentum transformations.";
highSymmetryKPath::option =
  "Radius must be positive; Center must be a real numeric three-vector; CenterLabel must be Automatic or a string; EquivalentDirections must be \"Representatives\" or \"All\"; IncludeNegativeDirections and ReturnData must be True or False; MaxGroupOrder and MaxSubgroups must be positive integers.";
highSymmetryKPath::group =
  "The finite little co-group could not be completed within MaxGroupOrder -> `1`.";
highSymmetryKPath::subgroups =
  "Subgroup enumeration exceeded MaxSubgroups -> `1`.";
highSymmetryKPath::lines =
  "No subgroup with a one-dimensional fixed momentum space was found.";

Options[highSymmetryKPath] = {
  "Radius" -> 1/10,
  "Center" -> {0, 0, 0},
  "CenterLabel" -> Automatic,
  "EquivalentDirections" -> "Representatives",
  "IncludeNegativeDirections" -> True,
  "ReturnData" -> False,
  "MaxGroupOrder" -> 512,
  "MaxSubgroups" -> 4096
};

highSymmetryKPath[input_Association, OptionsPattern[]] := Module[
  {
    radius = OptionValue["Radius"],
    center = OptionValue["Center"],
    centerLabel = OptionValue["CenterLabel"],
    equivalentDirections = OptionValue["EquivalentDirections"],
    includeNegative = OptionValue["IncludeNegativeDirections"],
    returnData = OptionValue["ReturnData"],
    maxGroupOrder = OptionValue["MaxGroupOrder"],
    maxSubgroups = OptionValue["MaxSubgroups"],
    generators,
    elements,
    table,
    subgroups,
    allLines,
    selectedLines,
    directions,
    path,
    label
  },
  If[
    !NumericQ[radius] || !TrueQ[radius > 0] ||
      !ListQ[center] || Length[center] =!= 3 ||
      !VectorQ[
        center,
        NumericQ[#] && TrueQ[Chop[Im[N[#]]] === 0] &
      ] ||
      !(centerLabel === Automatic || StringQ[centerLabel]) ||
      !MemberQ[{"Representatives", "All"}, equivalentDirections] ||
      !MemberQ[{True, False}, includeNegative] ||
      !MemberQ[{True, False}, returnData] ||
      !IntegerQ[maxGroupOrder] || maxGroupOrder <= 0 ||
      !IntegerQ[maxSubgroups] || maxSubgroups <= 0,
    Message[highSymmetryKPath::option];
    Return[$Failed]
  ];
  generators = kPathSymmetryGenerators[input];
  If[generators === $Failed || generators === {},
    Message[highSymmetryKPath::input];
    Return[$Failed]
  ];
  elements = kPathGroupClosure[generators, maxGroupOrder];
  If[MatchQ[elements, Failure["GroupOrderLimit", _Association]],
    Message[highSymmetryKPath::group, maxGroupOrder];
    Return[$Failed]
  ];
  table = kPathProductTable[elements];
  If[table === $Failed,
    Message[highSymmetryKPath::input];
    Return[$Failed]
  ];
  subgroups = kPathEnumerateSubgroups[table, maxSubgroups];
  If[MatchQ[subgroups, Failure["SubgroupLimit", _Association]],
    Message[highSymmetryKPath::subgroups, maxSubgroups];
    Return[$Failed]
  ];
  allLines = kPathLineData[elements, subgroups];
  If[allLines === {},
    Message[highSymmetryKPath::lines];
    Return[$Failed]
  ];
  selectedLines = If[equivalentDirections === "All",
    allLines,
    Fold[
      Function[{representatives, candidate},
        If[
          AnyTrue[
            representatives,
            kPathEquivalentLineQ[
              candidate["Direction"],
              #["Direction"],
              elements
            ] &
          ],
          representatives,
          Append[representatives, candidate]
        ]
      ],
      {},
      allLines
    ]
  ];
  centerLabel = If[
    centerLabel === Automatic,
    Lookup[input, "KPointLabel", "k0"],
    centerLabel
  ];
  If[!StringQ[centerLabel],
    Message[highSymmetryKPath::option];
    Return[$Failed]
  ];
  directions = Lookup[selectedLines, "Direction"];
  path = Flatten[
    MapIndexed[
      Function[{direction, index},
        label = "L" <> ToString[First[index]];
        If[includeNegative,
          {
            {
              {center - radius direction, center},
              {label <> "-", centerLabel}
            },
            {
              {center, center + radius direction},
              {centerLabel, label <> "+"}
            }
          },
          {{
            {center, center + radius direction},
            {centerLabel, label}
          }}
        ]
      ],
      directions
    ],
    1
  ];
  If[returnData,
    <|
      "Path" -> path,
      "Directions" -> directions,
      "LineData" -> selectedLines,
      "LittleCoGroupOrder" -> Length[elements],
      "LittleCoGroupActionOrder" -> Length @ DeleteDuplicates[
        kPathElementAction /@ elements,
        kPathZeroMatrixQ[#1 - #2] &
      ],
      "LittleCoGroupMatrices" -> (kPathElementAction /@ elements),
      "LittleCoGroupElements" -> (kPathElementData /@ elements),
      "EnumeratedSubgroupCount" -> Length[subgroups],
      "Center" -> center,
      "CenterLabel" -> centerLabel,
      "Radius" -> radius,
      "CoordinateSystem" -> Lookup[
        input,
        "CoordinateSystem",
        Missing["NotRecorded"]
      ]
    |>,
    path
  ]
];

highSymmetryKPath[___] := (
  Message[highSymmetryKPath::input];
  $Failed
);

End[]

EndPackage[]
