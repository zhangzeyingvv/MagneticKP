(* ::Package:: *)

BeginPackage["MagneticKP`"]


Begin["`Private`"]


vaspEig[filename_,efermi_,spin_,startband_,endband_]:=Module[{nband,nk,fermi,eigfile,todata,band},eigfile=StringSplit[Import[filename],"\n"];
eigfile=StringReplace[#,{"E+":>"*^","E-":>"*^-"}]&[eigfile];
fermi=efermi;
{nk,nband}=ToExpression[StringSplit[eigfile[[6]]][[-2;;-1]]];
todata=ToExpression[StringSplit/@#[[2;;]]]&;
(*Print[nk,Partition[eigfile[[7;;-1]],nband+2][[1]]];*)
band=Table[{#[[1]][[1;;3]],#[[2;;]][[;;,spin+1]]-fermi}&@todata@Partition[eigfile[[7;;-1]],nband+2][[i]],{i,1,nk,1}];
{#[[1]],#[[2]][[startband;;endband]]}&/@band];



Options[bandManipulateWithVASP]={"PlotRange"->All,"ManipulateRange"->2Pi};
bandManipulateWithVASP[h_,eigdata_,OptionsPattern[]]:=Module[{h0,params,mparams,m,rule,klist,mpstr,debug=False,Range=OptionValue["PlotRange"],MRange=OptionValue["ManipulateRange"]},params=Complement[Variables[h],{kx,ky,kz}];
If[debug,Print[params]];
klist=eigdata[[;;,1]];
mparams=ToExpression[StringRiffle[#,"c"]&/@(params/.{Subscript->List})];
If[debug,Print[mparams],ToString[eigdata[[;;,2]]]];
Print["Number of params:",Length@params];
Print["params:",params];
(*Print["mparams:",mparams];*)rule=ToString[Thread[params->mparams],StandardForm](*~Join~{kx\[Rule]k[[1]],ky\[Rule]k[[2]],kz\[Rule]k[[3]]}*);
(*m=StringTake[ToString[{{#,0,StringTake[ToString[#],2;;]},-1,1}&/@mparams],{2;;-2}];*)If[debug,Print[rule]];
MRange=If[ListQ[MRange],MRange,{-Abs[MRange],Abs[MRange]}];
m=StringTake[ToString[{{#1,0,#2},MRange[[1]],MRange[[2]]}&@@@Transpose[{mparams,params}],InputForm],{2;;-2}];
If[debug,Print[mparams,m]];
mpstr="Manipulate[\[IndentingNewLine]Show[ListPlot[Transpose@Table[Eigenvalues[Evaluate["<>ToString[h,StandardForm]<>"/."<>rule<>"~Join~{kx\[Rule]k\[LeftDoubleBracket]1\[RightDoubleBracket],ky\[Rule]k\[LeftDoubleBracket]2\[RightDoubleBracket],kz\[Rule]k\[LeftDoubleBracket]3\[RightDoubleBracket]}]],{k,"<>ToString[klist]<>"}],
     PlotRange->"<>ToString[Range]<>",PlotStyle->Black],
ListPlot[Transpose@"<>ToString[eigdata[[;;,2]]]<>",PlotStyle->Red]],"<>m<>"
     ,Button[\"ExportData\",Print["<>rule<>"]]
     ]";
If[debug,Print[m]];
If[debug,Print[mpstr]];
ToExpression[mpstr]];


ClearAll[
  fittingKPBandIndices,
  fittingKPBandSeries,
  fittingKPComparisonPlot,
  fittingKPKNeighborhoodSpec,
  fittingKPKPointSelectedQ,
  fittingKPLevenbergMarquardt,
  fittingKPLinearHamiltonianModel,
  fittingKPNormalizeHamiltonian,
  fittingKPNumericEigenvalues,
  fittingKPParameterExpressions,
  fittingKPRealNumericQ,
  fittingKPRulesList,
  fittingKPValidEigenvalueDataQ,
  fittingKPValidHamiltonianQ
];

fittingKPNormalizeHamiltonian[hamiltonian_] :=
  hamiltonian /. HoldPattern[MatrixForm[matrix_]] :> matrix;

fittingKPValidHamiltonianQ[hamiltonian_] := Module[{dimensions},
  dimensions = Quiet@Check[Dimensions[hamiltonian], Missing[]];
  MatrixQ[hamiltonian] &&
    MatchQ[dimensions, {_Integer, _Integer}] &&
    dimensions[[1]] > 0 &&
    Equal @@ dimensions
];

fittingKPRealNumericQ[value_] :=
  NumericQ[value] && TrueQ[Chop[Im[N[value]]] === 0];

fittingKPValidEigenvalueDataQ[data_] :=
  ListQ[data] && data =!= {} &&
    And @@ Map[
      Function[record,
        MatchQ[record, {_List, _List}] &&
          Length[record[[1]]] === 3 &&
          VectorQ[record[[1]], fittingKPRealNumericQ] &&
          record[[2]] =!= {} &&
          VectorQ[record[[2]], fittingKPRealNumericQ]
      ],
      data
    ] &&
    SameQ @@ (Length /@ data[[All, 2]]);

fittingKPParameterExpressions[hamiltonian_] :=
  Complement[Variables[hamiltonian], {kx, ky, kz}];

fittingKPNumericEigenvalues[matrix_] :=
  Sort@Chop@Eigenvalues[N[matrix]];

fittingKPRulesList[rules_] := Which[
  AssociationQ[rules], Normal[rules],
  MatchQ[rules, {(_Rule | _RuleDelayed) ...}], rules,
  True, $Failed
];

fittingKPBandIndices[selection_, bandCount_Integer] := Module[
  {indices, arguments},
  indices = Which[
    selection === All,
      Range[bandCount],
    IntegerQ[selection],
      {selection},
    ListQ[selection],
      selection,
    Head[selection] === Span,
      arguments = List @@ selection;
      Which[
        MatchQ[arguments, {_Integer, _Integer}],
          Range @@ arguments,
        MatchQ[arguments, {_Integer, _Integer, _Integer}] &&
          Last[arguments] > 0,
          Range @@ arguments,
        True,
          $Failed
      ],
    True,
      $Failed
  ];
  If[
    indices === $Failed || indices === {} ||
      !VectorQ[indices, IntegerQ] ||
      !And @@ (Between[#, {1, bandCount}] & /@ indices) ||
      Length[DeleteDuplicates[indices]] =!= Length[indices] ||
      !And @@ Thread[Rest[indices] > Most[indices]],
    $Failed,
    indices
  ]
];

fittingKPKNeighborhoodSpec[All] := <|"Mode" -> "All"|>;
fittingKPKNeighborhoodSpec[input_] := Module[
  {allowedKeys, center, periodic, hasRadius, hasRange, radius, range},
  If[!AssociationQ[input], Return[$Failed]];
  allowedKeys = {"Center", "Radius", "Range", "Periodic"};
  If[Complement[Keys[input], allowedKeys] =!= {}, Return[$Failed]];
  If[!KeyExistsQ[input, "Center"], Return[$Failed]];
  center = Lookup[input, "Center"];
  periodic = Lookup[input, "Periodic", True];
  hasRadius = KeyExistsQ[input, "Radius"];
  hasRange = KeyExistsQ[input, "Range"];
  If[
    Length[center] =!= 3 ||
      !VectorQ[center, fittingKPRealNumericQ] ||
      !MemberQ[{True, False}, periodic] ||
      Boole[hasRadius] + Boole[hasRange] =!= 1,
    Return[$Failed]
  ];
  If[hasRadius,
    radius = Lookup[input, "Radius"];
    If[!fittingKPRealNumericQ[radius] || !TrueQ[radius >= 0],
      Return[$Failed]
    ];
    <|
      "Mode" -> "Radius",
      "Center" -> N[center],
      "Radius" -> N[radius],
      "Periodic" -> periodic
    |>,
    range = Lookup[input, "Range"];
    If[
      Length[range] =!= 3 ||
        !VectorQ[range, fittingKPRealNumericQ] ||
        !And @@ (TrueQ[# >= 0] & /@ range),
      Return[$Failed]
    ];
    <|
      "Mode" -> "Range",
      "Center" -> N[center],
      "Range" -> N[range],
      "Periodic" -> periodic
    |>
  ]
];

fittingKPKPointSelectedQ[_, <|"Mode" -> "All"|>] := True;
fittingKPKPointSelectedQ[point_List, specification_Association] := Module[
  {delta},
  delta = N[point - specification["Center"]];
  If[TrueQ[specification["Periodic"]], delta -= Round[delta]];
  Switch[
    specification["Mode"],
    "Radius", TrueQ[Norm[delta] <= specification["Radius"]],
    "Range", And @@ Thread[Abs[delta] <= specification["Range"]],
    _, False
  ]
];

fittingKPLinearHamiltonianModel[
    hamiltonian_List,
    params_List,
    momenta_List
  ] := Module[
  {
    coefficientRuleMatrix,
    coefficientAssociationMatrix,
    zeroPowers,
    constantTemplate,
    coefficientTemplates,
    constantMatrices,
    coefficientMatrices
  },
  coefficientRuleMatrix = Quiet@Check[
    Map[CoefficientRules[#, params] &, hamiltonian, {2}],
    $Failed
  ];
  If[
    coefficientRuleMatrix === $Failed ||
      !And @@ Flatten@Map[
        AllTrue[First /@ #, Total[#] <= 1 &] &,
        coefficientRuleMatrix,
        {2}
      ],
    Return[Failure["NonlinearParameterModel", <||>]]
  ];
  coefficientAssociationMatrix = Map[Association, coefficientRuleMatrix, {2}];
  zeroPowers = ConstantArray[0, Length[params]];
  constantTemplate = Map[
    Lookup[#, Key[zeroPowers], 0] &,
    coefficientAssociationMatrix,
    {2}
  ];
  coefficientTemplates = Table[
    Map[
      Lookup[#, Key[UnitVector[Length[params], index]], 0] &,
      coefficientAssociationMatrix,
      {2}
    ],
    {index, Length[params]}
  ];
  constantMatrices = N@Map[
    Function[k,
      constantTemplate /. Thread[{kx, ky, kz} -> k]
    ],
    momenta
  ];
  coefficientMatrices = N@Map[
    Function[template,
      Map[
        Function[k,
          template /. Thread[{kx, ky, kz} -> k]
        ],
        momenta
      ]
    ],
    coefficientTemplates
  ];
  If[
    !And @@ (MatrixQ[#, NumericQ] & /@ constantMatrices) ||
      !And @@ Flatten[
        Map[MatrixQ[#, NumericQ] &, coefficientMatrices, {2}]
      ],
    Return[Failure["NonNumericLinearModel", <||>]]
  ];
  <|
    "ConstantMatrices" -> constantMatrices,
    "CoefficientMatrices" -> coefficientMatrices
  |>
];

fittingKPLevenbergMarquardt[
    residualFunction_,
    initialValues_List,
    maxIterations_Integer,
    timeConstraint_,
    finiteDifferenceStep_?NumericQ,
    initialDamping_?NumericQ,
    tolerance_?NumericQ
  ] := Module[
  {
    values = N[initialValues],
    residual,
    loss,
    damping = N[initialDamping],
    startTime = AbsoluteTime[],
    timedOutQ,
    iteration,
    stepSizes,
    unitVectors,
    columns,
    columnFailure,
    jacobian,
    normalMatrix,
    gradient,
    dampingScale,
    delta,
    candidateValues,
    candidateResidual,
    candidateLoss,
    improvement,
    converged = False
  },
  timedOutQ[] := timeConstraint =!= Infinity &&
    AbsoluteTime[] - startTime >= timeConstraint;
  residual = Quiet@Check[N@residualFunction[values], $Failed];
  If[!VectorQ[residual, NumericQ],
    Return[Failure["InvalidResidual", <||>]]
  ];
  loss = residual.residual;
  unitVectors = IdentityMatrix[Length[values]];

  For[iteration = 1, iteration <= maxIterations, iteration++,
    If[timedOutQ[],
      Return[Failure[
        "TimeConstraint",
        <|
          "Values" -> values,
          "Loss" -> loss,
          "Iterations" -> iteration - 1
        |>
      ]]
    ];
    stepSizes = finiteDifferenceStep (Max[1., Abs[#]] & /@ values);
    columns = Table[
      If[timedOutQ[],
        Return[Failure[
          "TimeConstraint",
          <|
            "Values" -> values,
            "Loss" -> loss,
            "Iterations" -> iteration - 1
          |>
        ]],
        Module[{plusResidual, minusResidual},
          plusResidual = Quiet@Check[
            N@residualFunction[
              values + stepSizes[[index]] unitVectors[[index]]
            ],
            $Failed
          ];
          minusResidual = Quiet@Check[
            N@residualFunction[
              values - stepSizes[[index]] unitVectors[[index]]
            ],
            $Failed
          ];
          If[
            !VectorQ[plusResidual, NumericQ] ||
              !VectorQ[minusResidual, NumericQ] ||
              Length[plusResidual] =!= Length[residual] ||
              Length[minusResidual] =!= Length[residual],
            Return[Failure["InvalidResidual", <||>], Module]
          ];
          (plusResidual - minusResidual)/(2 stepSizes[[index]])
        ]
      ],
      {index, Length[values]}
    ];
    columnFailure = SelectFirst[columns, FailureQ, Missing["NoFailure"]];
    If[!MissingQ[columnFailure], Return[columnFailure]];
    jacobian = Transpose[columns];
    normalMatrix = Transpose[jacobian].jacobian;
    gradient = Transpose[jacobian].residual;
    If[Norm[gradient, Infinity] <= tolerance,
      converged = True;
      Break[]
    ];
    dampingScale = Max[1., #] & /@ Diagonal[normalMatrix];
    delta = Quiet@Check[
      LinearSolve[
        normalMatrix + damping DiagonalMatrix[dampingScale],
        -gradient
      ],
      $Failed
    ];
    If[!VectorQ[delta, NumericQ],
      damping *= 10.;
      Continue[]
    ];
    candidateValues = values + delta;
    candidateResidual = Quiet@Check[
      N@residualFunction[candidateValues],
      $Failed
    ];
    If[!VectorQ[candidateResidual, NumericQ],
      damping *= 10.;
      Continue[]
    ];
    candidateLoss = candidateResidual.candidateResidual;
    If[candidateLoss < loss,
      improvement = loss - candidateLoss;
      values = candidateValues;
      residual = candidateResidual;
      loss = candidateLoss;
      damping = Max[damping/3., 10.^-15];
      If[
        Norm[delta] <= tolerance (1. + Norm[values]) ||
          improvement <= tolerance Max[1., loss],
        converged = True;
        Break[]
      ],
      damping = Min[damping*10., 10.^15]
    ]
  ];

  <|
    "Values" -> values,
    "Loss" -> loss,
    "Iterations" -> Min[iteration, maxIterations],
    "Converged" -> converged,
    "FinalDamping" -> damping
  |>
];

fittingKPBandSeries[indices_List, bands_List, selectionMask_List] :=
  Select[
    MapThread[
      Pick[MapThread[List, {indices, #1}], #2, True] &,
      {Transpose[bands], Transpose[selectionMask]}
    ],
    # =!= {} &
  ];

fittingKPComparisonPlot[
    indices_List,
    referenceBands_List,
    initialBands_List,
    fittedBands_List,
    selectionMask_List
  ] := Legended[
  Show[
    ListPlot[
      fittingKPBandSeries[indices, referenceBands, selectionMask],
      PlotStyle -> Directive[Blue, AbsolutePointSize[5]],
      PlotMarkers -> Automatic
    ],
    ListLinePlot[
      fittingKPBandSeries[indices, initialBands, selectionMask],
      PlotStyle -> Directive[Gray, Dashed, Thin]
    ],
    ListLinePlot[
      fittingKPBandSeries[indices, fittedBands, selectionMask],
      PlotStyle -> Directive[Red, Thick]
    ],
    Frame -> True,
    FrameLabel -> {"k-point index", "Energy"},
    PlotRange -> All,
    ImageSize -> Large
  ],
  Placed[
    LineLegend[
      {
        Directive[Blue, Thick],
        Directive[Gray, Dashed],
        Directive[Red, Thick]
      },
      {"VASP", "Initial k.p", "Fitted k.p"}
    ],
    Above
  ]
];

fittingKP::matrix =
  "The Hamiltonian must be a nonempty square matrix after display wrappers are removed.";
fittingKP::data =
  "Reference data must contain numeric three-component k points and exactly `1` energies per selected point.";
fittingKP::range =
  "krange must be a nonempty list of valid 1-based reference-data indices.";
fittingKP::energywindow =
  "EnergyWindow must be All or an inclusive real numeric interval {emin, emax} with emin <= emax.";
fittingKP::kneighborhood =
  "KPointNeighborhood must be All or an Association with Center, exactly one of Radius or Range, and optional Periodic -> True or False.";
fittingKP::bands =
  "BandSelection must be All, a valid 1-based band index, a strictly increasing list of unique indices, or a positive-step Span.";
fittingKP::weights =
  "ResidualWeights must be Automatic or a nonnegative real numeric matrix with the same record-by-band shape as eigdata.";
fittingKP::selection =
  "The combined k-point, band, energy-window, and positive-weight selection contains no residuals.";
fittingKP::insufficient =
  "The selection contains only `1` residuals for `2` fitting parameters.";
fittingKP::init =
  "initialRules must assign a real numeric initial value to every fitting parameter: `1`.";
fittingKP::nonlinear =
  "The Hamiltonian must be affine-linear in every fitting parameter.";
fittingKP::parameter =
  "The Hamiltonian remains symbolic after parameter and momentum processing.";
fittingKP::fit =
  "The Levenberg-Marquardt fit failed for the supplied Hamiltonian and initial values.";
fittingKP::option =
  "MaxIterations must be a positive integer; TimeConstraint, FiniteDifferenceStep, InitialDamping, and FitTolerance must be positive (TimeConstraint may also be Infinity).";
fittingKP::timeout =
  "The fit exceeded the requested TimeConstraint of `1` seconds.";

Options[fittingKP] = {
  MaxIterations -> 100,
  TimeConstraint -> 60,
  "EnergyWindow" -> All,
  "KPointNeighborhood" -> All,
  "BandSelection" -> All,
  "ResidualWeights" -> Automatic,
  "FiniteDifferenceStep" -> 10^-5,
  "InitialDamping" -> 10^-3,
  "FitTolerance" -> 10^-8
};

fittingKP[h_, eigdata_, krange_, initparms_, OptionsPattern[]] := Module[
  {
    hamiltonian,
    dimensions,
    bandCount,
    params,
    initialRules,
    initialKeys,
    initialValues,
    energyWindowOption,
    energyWindow,
    kNeighborhoodOption,
    kNeighborhood,
    bandSelectionOption,
    bandIndices,
    weightOption,
    allWeights,
    candidateKPointCount,
    selectedKIndices,
    selectedData,
    targetBands,
    selectedWeights,
    bandMask,
    selectionMask,
    usedKPointMask,
    selectionMaskFlat,
    residualWeights,
    residualCount,
    selectionRestrictedQ,
    selectedBandIndicesByKPoint,
    selectionSummary,
    objectiveName,
    momenta,
    linearModel,
    constantMatrices,
    coefficientMatrices,
    numericHamiltonians,
    residual,
    initialResidual,
    initialLoss,
    initialBands,
    maxIterations,
    timeConstraint,
    finiteDifferenceStep,
    initialDamping,
    fitTolerance,
    fitResult,
    fittedRules,
    fittedBands,
    comparisonPlot
  },
  hamiltonian = fittingKPNormalizeHamiltonian[h];
  If[!fittingKPValidHamiltonianQ[hamiltonian],
    Message[fittingKP::matrix];
    Return[$Failed]
  ];
  dimensions = Dimensions[hamiltonian];
  bandCount = dimensions[[1]];
  If[
    !fittingKPValidEigenvalueDataQ[eigdata] ||
      !And @@ (Length[#] === bandCount & /@ eigdata[[All, 2]]),
    Message[fittingKP::data, bandCount];
    Return[$Failed]
  ];
  If[
    !ListQ[krange] || krange === {} ||
      !VectorQ[krange, IntegerQ] ||
      !And @@ (Between[#, {1, Length[eigdata]}] & /@ krange),
    Message[fittingKP::range];
    Return[$Failed]
  ];

  params = fittingKPParameterExpressions[hamiltonian];
  initialRules = fittingKPRulesList[initparms];
  If[initialRules === $Failed,
    Message[fittingKP::init, params];
    Return[$Failed]
  ];
  initialKeys = First /@ initialRules;
  If[
    Length[initialKeys] =!= Length[params] ||
      Length[DeleteDuplicates[initialKeys]] =!= Length[initialKeys] ||
      !And @@ (MemberQ[initialKeys, #] & /@ params),
    Message[fittingKP::init, params];
    Return[$Failed]
  ];
  initialValues = params /. initialRules;
  If[!VectorQ[initialValues, fittingKPRealNumericQ],
    Message[fittingKP::init, params];
    Return[$Failed]
  ];

  maxIterations = OptionValue[MaxIterations];
  timeConstraint = OptionValue[TimeConstraint];
  finiteDifferenceStep = OptionValue["FiniteDifferenceStep"];
  initialDamping = OptionValue["InitialDamping"];
  fitTolerance = OptionValue["FitTolerance"];
  If[
    !IntegerQ[maxIterations] || maxIterations <= 0 ||
      !(timeConstraint === Infinity ||
        (NumericQ[timeConstraint] && TrueQ[timeConstraint > 0])) ||
      !NumericQ[finiteDifferenceStep] ||
      !TrueQ[finiteDifferenceStep > 0] ||
      !NumericQ[initialDamping] || !TrueQ[initialDamping > 0] ||
      !NumericQ[fitTolerance] || !TrueQ[fitTolerance > 0],
    Message[fittingKP::option];
    Return[$Failed]
  ];

  energyWindowOption = OptionValue["EnergyWindow"];
  energyWindow = Which[
    energyWindowOption === All,
      All,
    ListQ[energyWindowOption] && Length[energyWindowOption] === 2 &&
      VectorQ[energyWindowOption, fittingKPRealNumericQ] &&
      TrueQ[energyWindowOption[[1]] <= energyWindowOption[[2]]],
      N[energyWindowOption],
    True,
      $Failed
  ];
  If[energyWindow === $Failed,
    Message[fittingKP::energywindow];
    Return[$Failed]
  ];

  kNeighborhoodOption = OptionValue["KPointNeighborhood"];
  kNeighborhood = fittingKPKNeighborhoodSpec[kNeighborhoodOption];
  If[kNeighborhood === $Failed,
    Message[fittingKP::kneighborhood];
    Return[$Failed]
  ];

  bandSelectionOption = OptionValue["BandSelection"];
  bandIndices = fittingKPBandIndices[bandSelectionOption, bandCount];
  If[bandIndices === $Failed,
    Message[fittingKP::bands];
    Return[$Failed]
  ];

  weightOption = OptionValue["ResidualWeights"];
  allWeights = Which[
    weightOption === Automatic,
      ConstantArray[1., {Length[eigdata], bandCount}],
    Dimensions[weightOption] === {Length[eigdata], bandCount} &&
      MatrixQ[
        weightOption,
        fittingKPRealNumericQ[#] && TrueQ[# >= 0] &
      ],
      N[weightOption],
    True,
      $Failed
  ];
  If[allWeights === $Failed,
    Message[fittingKP::weights];
    Return[$Failed]
  ];

  selectedKIndices = Select[
    krange,
    fittingKPKPointSelectedQ[eigdata[[#, 1]], kNeighborhood] &
  ];
  If[selectedKIndices === {},
    Message[fittingKP::selection];
    Return[$Failed]
  ];
  candidateKPointCount = Length[selectedKIndices];
  selectedData = eigdata[[selectedKIndices]];
  targetBands = N[selectedData[[All, 2]]];
  selectedWeights = allWeights[[selectedKIndices]];
  bandMask = MemberQ[bandIndices, #] & /@ Range[bandCount];
  selectionMask = MapThread[
    Function[{energies, weights},
      MapThread[
        Function[{bandAllowed, energy, weight},
          TrueQ[
            bandAllowed && weight > 0 &&
              (energyWindow === All || Between[energy, energyWindow])
          ]
        ],
        {bandMask, energies, weights}
      ]
    ],
    {targetBands, selectedWeights}
  ];
  usedKPointMask = Or @@ # & /@ selectionMask;
  If[!MemberQ[usedKPointMask, True],
    Message[fittingKP::selection];
    Return[$Failed]
  ];
  selectedKIndices = Pick[selectedKIndices, usedKPointMask, True];
  selectedData = Pick[selectedData, usedKPointMask, True];
  targetBands = Pick[targetBands, usedKPointMask, True];
  selectedWeights = Pick[selectedWeights, usedKPointMask, True];
  selectionMask = Pick[selectionMask, usedKPointMask, True];
  selectionMaskFlat = Flatten[selectionMask];
  residualWeights = Pick[Flatten[selectedWeights], selectionMaskFlat, True];
  residualCount = Length[residualWeights];
  If[residualCount === 0,
    Message[fittingKP::selection];
    Return[$Failed]
  ];
  selectionRestrictedQ = !(
    energyWindowOption === All &&
      kNeighborhoodOption === All &&
      bandSelectionOption === All &&
      weightOption === Automatic
  );
  If[selectionRestrictedQ && residualCount < Length[params],
    Message[fittingKP::insufficient, residualCount, Length[params]];
    Return[$Failed]
  ];
  selectedBandIndicesByKPoint =
    Pick[Range[bandCount], #, True] & /@ selectionMask;
  objectiveName = If[
    weightOption === Automatic,
    "SumSquaredBandResiduals",
    "WeightedSumSquaredBandResiduals"
  ];
  selectionSummary = <|
    "CandidateKPointCount" -> candidateKPointCount,
    "UsedKPointCount" -> Length[selectedKIndices],
    "UniqueKPointCount" -> Length[DeleteDuplicates[selectedKIndices]],
    "UsedKPointIndices" -> selectedKIndices,
    "ResidualCount" -> residualCount,
    "SelectionRestricted" -> selectionRestrictedQ,
    "SelectedBandIndicesByKPoint" -> selectedBandIndicesByKPoint,
    "BandSelection" -> bandIndices,
    "EnergyWindow" -> energyWindow,
    "KPointNeighborhood" -> kNeighborhood,
    "Weighting" -> If[weightOption === Automatic, "Uniform", "ExplicitMatrix"],
    "BandPairing" -> "SortedModelEigenvalueToReferencePosition",
    "MomentumConvention" -> "DirectEigdataCoordinates"
  |>;

  (* fittingKP preserves the original MagneticKP convention and substitutes
     eigdata coordinates directly.  Unlike fittingTB, no factor of 2 Pi is
     inserted here. *)
  momenta = N[selectedData[[All, 1]]];

  If[params === {},
    constantMatrices = N@Map[
      Function[k,
        hamiltonian /. Thread[{kx, ky, kz} -> k]
      ],
      momenta
    ];
    If[!And @@ (MatrixQ[#, NumericQ] & /@ constantMatrices),
      Message[fittingKP::parameter];
      Return[$Failed]
    ];
    initialBands = fittingKPNumericEigenvalues /@ constantMatrices;
    initialResidual = Sqrt[residualWeights] Pick[
      Flatten[initialBands - targetBands],
      selectionMaskFlat,
      True
    ];
    initialLoss = initialResidual.initialResidual;
    comparisonPlot = fittingKPComparisonPlot[
      selectedKIndices,
      targetBands,
      initialBands,
      initialBands,
      selectionMask
    ];
    Return[Join[
      <|
        "FittedParams" -> {},
        "Objective" -> objectiveName,
        "Optimizer" -> "None",
        "ParameterModel" -> "Constant",
        "Converged" -> True,
        "Iterations" -> 0,
        "FinalDamping" -> 0.,
        "LSQForInitParams" -> initialLoss,
        "LSQForFittedParams" -> initialLoss,
        "ComparisonPlot" -> comparisonPlot
      |>,
      selectionSummary
    ]]
  ];

  linearModel = fittingKPLinearHamiltonianModel[
    hamiltonian,
    params,
    momenta
  ];
  If[MatchQ[linearModel, Failure["NonlinearParameterModel", _Association]],
    Message[fittingKP::nonlinear];
    Return[$Failed]
  ];
  If[!AssociationQ[linearModel],
    Message[fittingKP::parameter];
    Return[$Failed]
  ];
  constantMatrices = linearModel["ConstantMatrices"];
  coefficientMatrices = linearModel["CoefficientMatrices"];
  numericHamiltonians[values_List] /; VectorQ[values, NumericQ] :=
    constantMatrices + Total[
      MapThread[Times, {N[values], coefficientMatrices}]
    ];
  residual[values_List] /; VectorQ[values, NumericQ] := Module[
    {modelBands},
    modelBands = fittingKPNumericEigenvalues /@ numericHamiltonians[values];
    Sqrt[residualWeights] Pick[
      Flatten[modelBands - targetBands],
      selectionMaskFlat,
      True
    ]
  ];

  initialBands = fittingKPNumericEigenvalues /@
    numericHamiltonians[N[initialValues]];
  initialResidual = residual[N[initialValues]];
  initialLoss = initialResidual.initialResidual;
  fitResult = fittingKPLevenbergMarquardt[
    Function[values, residual[values]],
    N[initialValues],
    maxIterations,
    timeConstraint,
    finiteDifferenceStep,
    initialDamping,
    fitTolerance
  ];
  If[MatchQ[fitResult, Failure["TimeConstraint", _Association]],
    Message[fittingKP::timeout, timeConstraint];
    Return[$Failed]
  ];
  If[!AssociationQ[fitResult],
    Message[fittingKP::fit];
    Return[$Failed]
  ];
  fittedRules = Thread[params -> fitResult["Values"]];
  fittedBands = fittingKPNumericEigenvalues /@
    numericHamiltonians[fitResult["Values"]];
  comparisonPlot = fittingKPComparisonPlot[
    selectedKIndices,
    targetBands,
    initialBands,
    fittedBands,
    selectionMask
  ];

  Join[
    <|
      "FittedParams" -> fittedRules,
      "Objective" -> objectiveName,
      "Optimizer" -> "LevenbergMarquardt",
      "ParameterModel" -> "CachedAffineCoefficientMatrices",
      "Converged" -> fitResult["Converged"],
      "Iterations" -> fitResult["Iterations"],
      "FinalDamping" -> fitResult["FinalDamping"],
      "LSQForInitParams" -> initialLoss,
      "LSQForFittedParams" -> fitResult["Loss"],
      "ComparisonPlot" -> comparisonPlot
    |>,
    selectionSummary
  ]
]


End[]
EndPackage[]

