(* ::Package:: *)

BeginPackage["MagneticKP`"]


Begin["`Private`"]


texOutput[mat_]:=Module[{ptex},
ptex=ToString[TeXForm@mat];
StringReplace[ToString[ptex],{
RegularExpression["\\\\text\\{(\\D)(\\d)\\}"]->"$1"<>"_"<>"{$2}",
RegularExpression["\\\\text\\{(k)(.)\\}"]->"$1"<>"_"<>"$2",
RegularExpression["\\\\text\\{(\\D)(\\d+)(\\D)(\\d+)\\}"]->"$1"<>"_"<>"{$4}"<>"^"<>"{$2}"
}]];


(*GenerateGroup is form SpaceGroupIrep*)
(*Seteq[s1_,s2_]:=SubsetQ[s1,s2]&&SubsetQ[s2,s1];*)
GenerateGroup[gens_,identityElement_,multiply_]:=
 Module[{i,j,ng,MAXORDER=200,orders,subs,mlist,g1,g2,g3},
   ng=Length[gens];   orders=subs=Range[ng]*0;
   For[i=1,i<=ng,i++,
     mlist=FoldList[multiply,Table[gens[[i]],MAXORDER]];
     orders[[i]]=FirstPosition[mlist,identityElement][[1]];
     subs[[i]]=mlist[[;;orders[[i]]]];
   ];
   g1=Union@@subs;
   g2=Union@@Table[multiply[g1[[i]],g1[[j]]],{i,Length[g1]},{j,Length[g1]}];
   g3=Complement[g2,g1];  g1=g2;
   While[g3!={},
     g2=Union@@Table[multiply[g3[[i]],g1[[j]]],{i,Length[g3]},{j,Length[g1]}];
     g3=Complement[g2,g1];  g1=g2;
   ];
   g2
 ];
 
getGenerator[groupele_,identityElement_,multiply_]:=
Module[
{try,group,
tmptry,tmpgroup,
seteq,
groupeleDisorder,
gereratorlist
},
(*Greedy Alg, May not give the minimal Generator set*)
seteq[s1_,s2_]:=SubsetQ[s1,s2]&&SubsetQ[s2,s1];
If[groupele=={identityElement},Return[groupele]];
gereratorlist=Table[
try={};
(*groupeleDisorder=RandomSample[groupele];*)
groupeleDisorder=groupele;
Do[
If[ele==identityElement,Continue[]];
group=GenerateGroup[try,identityElement,multiply];
tmptry=Append[try,ele];
(*Print[group];*)
tmpgroup=GenerateGroup[tmptry,identityElement,multiply];
If[Not@seteq[group,tmpgroup],try=tmptry];
If[seteq[groupele,group],Break[]];
(*If[seteq[tmpgroup,groupele],Print[tmptry,tmpgroup];Break[],
try=tmptry
]*)
,{ele,groupeleDisorder}];
try,1];
First[SortBy[gereratorlist,Length[#]&]]
];
transformationInput[input_,U_]:=Association[{"Unitary"->MapAt[FullSimplify[U . # . Inverse[U]]&,input["Unitary"],{;;,1}],
"Anitunitary"->MapAt[FullSimplify[U . # . Conjugate@Inverse[U]]&,input["Anitunitary"],{;;,1}]
}];
ComplexRationalMatQ=And@@Flatten[Map[(Element[Re[#],Rationals]&&Element[Im[#],Rationals]&),#,{-1}]]&;
directSum =(*FullSimplify@*)ArrayFlatten[{{#1, 0}, {0, #2}}] &;
generateMatrixGroup[gen_]:=GenerateGroup[gen,IdentityMatrix@Length[First[gen]],Dot];


Options[bandManipulate]={"PlotRange"->All,"ManipulateRange"->1};
bandManipulate[pathstr_,npoint_,h_,OptionsPattern[]]:=Module[
{h0,params,mparams,m,rule,mpstr,klist,debug=False,Range=OptionValue["PlotRange"],
MRange=OptionValue["ManipulateRange"]
},
params=Complement[Variables[h],{kx,ky,kz}];
If[debug,Print[params]];

mparams=ToExpression[StringRiffle[#,"c"]&/@(params/.{Subscript->List})];
If[debug,Print[mparams]];
klist=2.Pi Flatten[Subdivide[#[[1]],#[[2]],npoint]&/@(Transpose[pathstr][[1]]),1];
Print["Number of params:",Length@params];
Print["params:",params];
(*  Print["mparams:",mparams];*)
rule=ToString[Thread[params->mparams],StandardForm](*~Join~{kx\[Rule]k[[1]],ky\[Rule]k[[2]],kz\[Rule]k[[3]]}*);
(*  m=StringTake[ToString[{{#,0,StringTake[ToString[#],2;;]},-1,1}&/@mparams],{2;;-2}];*)
If[debug,Print[rule]];
MRange=If[ListQ[MRange],MRange,{-Abs[MRange],Abs[MRange]}];
m=StringTake[ToString[{{#1,0,#2},MRange[[1]],MRange[[2]]}&@@@Transpose[{mparams,params}],InputForm],{2;;-2}];
If[debug,Print[mparams,m]];
mpstr="Manipulate[\[IndentingNewLine]ListPlot[Transpose@Table[Eigenvalues[Evaluate["<>ToString[h,StandardForm]<>"/."<>rule<>"~Join~{kx\[Rule]k\[LeftDoubleBracket]1\[RightDoubleBracket],ky\[Rule]k\[LeftDoubleBracket]2\[RightDoubleBracket],kz\[Rule]k\[LeftDoubleBracket]3\[RightDoubleBracket]}]],{k,"<>ToString[klist]<>"}],
PlotRange->"<>ToString[Range]<>",PlotStyle->Black],"<>m<>"
,Button[\"ExportData\",Print["<>rule<>"]]
]";
If[debug,Print[m]];
If[debug,Print[mpstr]];
ToExpression[mpstr]
];


ClearAll[
  validBandPathQ,
  bandPathBoundaryLabels,
  styledBandYTicks,
  bandEnergySegments,
  bandPlotData
];

bandplot::path =
  "The band path must contain one or more {{{k1},{k2}},{label1,label2}} segments in three-dimensional reciprocal fractional coordinates.";
bandplot::npoint =
  "The number of subdivisions per path segment must be a positive integer; received `1`.";
bandplot::ham = "The Hamiltonian must be a nonempty square matrix.";
bandplot::rules = "Parameter values must be supplied as a list of replacement rules.";
bandplot::numeric =
  "The Hamiltonian did not become a numeric matrix with real eigenvalues after applying the parameter and momentum rules.";
bandplot::yticks =
  "\"YTicks\" must be Automatic, None, a list of numeric positions, or a standard list of tick specifications.";
bandplot::font =
  "FontFamily must be a string and FontSize must be a positive number.";

validBandPathQ[path_] := MatchQ[
  path,
  {
    Repeated[
      {
        {
          {_?NumericQ, _?NumericQ, _?NumericQ},
          {_?NumericQ, _?NumericQ, _?NumericQ}
        },
        {_String, _String}
      }
    ]
  }
];

bandPathBoundaryLabels[path_List] := Module[{labels, segmentCount},
  segmentCount = Length[path];
  labels = ConstantArray["", segmentCount + 1];
  labels[[1]] = path[[1, 2, 1]];
  Do[
    labels[[index + 1]] = If[
      index < segmentCount &&
        path[[index, 2, 2]] =!= path[[index + 1, 2, 1]],
      path[[index, 2, 2]] <> "|" <> path[[index + 1, 2, 1]],
      path[[index, 2, 2]]
    ],
    {index, segmentCount}
  ];
  StringReplace[#, "\\Gamma" -> "\[CapitalGamma]"] & /@ labels
];

styledBandYTicks[Automatic, _, _] := Automatic;
styledBandYTicks[None, _, _] := None;
styledBandYTicks[ticks_List, fontFamily_, fontSize_] := Which[
  VectorQ[ticks, NumericQ],
    {
      #,
      Style[#, Black, FontFamily -> fontFamily, FontSize -> fontSize]
    } & /@ ticks,
  And @@ Map[
    ListQ[#] && Length[#] >= 2 && NumericQ[First[#]] &,
    ticks
  ],
    ReplacePart[
      #,
      2 -> Style[
        #[[2]],
        Black,
        FontFamily -> fontFamily,
        FontSize -> fontSize
      ]
    ] & /@ ticks,
  True,
    $Failed
];
styledBandYTicks[_, _, _] := $Failed;

bandEnergySegments[
    path_List,
    npoint_Integer?Positive,
    ham_,
    rules_List
  ] := Module[
  {segments, numericHamiltonian, eigenvalues},
  segments = Subdivide[#[[1]], #[[2]], npoint] & /@ path[[All, 1]];
  Catch[
    Map[
      Function[kpoints,
        Transpose@Map[
          Function[kpoint,
            numericHamiltonian = Quiet@Check[
              N[ham /. rules /. {
                kx -> 2 Pi kpoint[[1]],
                ky -> 2 Pi kpoint[[2]],
                kz -> 2 Pi kpoint[[3]]
              }],
              $Failed
            ];
            If[
              numericHamiltonian === $Failed ||
                !MatrixQ[numericHamiltonian, NumericQ] ||
                Length[numericHamiltonian] == 0 ||
                Length[numericHamiltonian] =!= Length[First[numericHamiltonian]],
              Throw[$Failed, "BandPlotNumeric"]
            ];
            eigenvalues = Quiet@Check[
              Chop@Eigenvalues[numericHamiltonian],
              $Failed
            ];
            If[
              eigenvalues === $Failed ||
                !VectorQ[eigenvalues, NumericQ] ||
                !And @@ (
                  TrueQ[Chop[Im[N[#]]] == 0] & /@ eigenvalues
                ),
              Throw[$Failed, "BandPlotNumeric"]
            ];
            Sort[Re@N[eigenvalues]]
          ],
          kpoints
        ]
      ],
      segments
    ],
    "BandPlotNumeric"
  ]
];

bandPlotData[
    path_List,
    npoint_Integer?Positive,
    ham_,
    rules_List
  ] := Module[{segmentEnergies},
  segmentEnergies = bandEnergySegments[path, npoint, ham, rules];
  If[segmentEnergies === $Failed, Return[$Failed]];
  Chop@Flatten[
    MapIndexed[
      Function[{segmentBands, segmentIndex},
        Map[
          MapIndexed[
            {
              npoint (First[segmentIndex] - 1) + First[#2] - 1,
              #1
            } &,
            #
          ] &,
          segmentBands
        ]
      ],
      segmentEnergies
    ],
    1
  ]
];

Options[bandplot] = {
  "PlotRange" -> All,
  "YTicks" -> Automatic,
  FontSize -> 24,
  FontFamily -> "Times",
  ImageSize -> Automatic
};

bandplot[
    path_,
    npoint_,
    ham_,
    rules_,
    OptionsPattern[]
  ] := Module[
  {
    labels, fontFamily, fontSize, xticks, yticks, data,
    segmentCount
  },
  If[!IntegerQ[npoint] || npoint < 1,
    Message[bandplot::npoint, npoint];
    Return[$Failed]
  ];
  If[!MatrixQ[ham] || Length[ham] == 0 || Length[ham] =!= Length[First[ham]],
    Message[bandplot::ham];
    Return[$Failed]
  ];
  If[!ListQ[rules],
    Message[bandplot::rules];
    Return[$Failed]
  ];
  If[!validBandPathQ[path],
    Message[bandplot::path];
    Return[$Failed]
  ];

  fontFamily = OptionValue[FontFamily];
  fontSize = OptionValue[FontSize];
  If[!StringQ[fontFamily] || !NumericQ[fontSize] || fontSize <= 0,
    Message[bandplot::font];
    Return[$Failed]
  ];
  yticks = styledBandYTicks[
    OptionValue["YTicks"],
    fontFamily,
    fontSize
  ];
  If[yticks === $Failed,
    Message[bandplot::yticks];
    Return[$Failed]
  ];

  data = bandPlotData[path, npoint, ham, rules];
  If[data === $Failed,
    Message[bandplot::numeric];
    Return[$Failed]
  ];

  segmentCount = Length[path];
  labels = bandPathBoundaryLabels[path];
  xticks = MapThread[
    {
      #1,
      Style[
        #2,
        Black,
        FontFamily -> fontFamily,
        FontSize -> fontSize
      ]
    } &,
    {npoint Range[0, segmentCount], labels}
  ];

  ListLinePlot[
    data,
    PlotRange -> {
      {0, npoint segmentCount},
      OptionValue["PlotRange"]
    },
    PlotStyle -> Black,
    GridLines -> {npoint Range[0, segmentCount], {0}},
    Frame -> True,
    FrameStyle -> Black,
    FrameTicks -> {{yticks, None}, {xticks, None}},
    FrameTicksStyle -> Directive[
      Black,
      FontFamily -> fontFamily,
      FontSize -> fontSize
    ],
    BaseStyle -> Directive[
      FontFamily -> fontFamily,
      FontSize -> fontSize
    ],
    GridLinesStyle -> Directive[Black],
    ImageSize -> OptionValue[ImageSize]
  ]
];


End[]
EndPackage[]

