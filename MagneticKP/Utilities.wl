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


Options[bandManipulate]={"PlotRange"->All};
bandManipulate[pathstr_,npoint_,h_,OptionsPattern[]]:=Module[
{h0,params,mparams,m,rule,mpstr,klist,debug=False,Range=OptionValue["PlotRange"]},
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
m=StringTake[ToString[{{#1,0,#2},-1,1}&@@@Transpose[{mparams,params}],InputForm],{2;;-2}];
If[debug,Print[mparams,m]];
mpstr="Manipulate[\[IndentingNewLine]ListPlot[Transpose@Table[Eigenvalues[Evaluate["<>ToString[h,StandardForm]<>"/."<>rule<>"~Join~{kx\[Rule]k\[LeftDoubleBracket]1\[RightDoubleBracket],ky\[Rule]k\[LeftDoubleBracket]2\[RightDoubleBracket],kz\[Rule]k\[LeftDoubleBracket]3\[RightDoubleBracket]}]],{k,"<>ToString[klist]<>"}],
PlotRange->"<>ToString[Range]<>",PlotStyle->Black],"<>m<>"
,Button[\"ExportData\",Print["<>rule<>"]]
]";
If[debug,Print[m]];
If[debug,Print[mpstr]];
ToExpression[mpstr]
];


End[]
EndPackage[]

