(* ::Package:: *)

(* ::Subsubsection:: *)
(*Start*)


BeginPackage["MagneticKP`"]


Begin["`Private`"]


(* ::Subsection:: *)
(*Interface with SpaceGroupIrep and MSGCorep*)


(* ::Subsubsection:: *)
(*Common*)


Options[interfacekp] = {"generator" -> None,"CalculateGenerators"->True};
interfacekp[kinfo_, OptionsPattern[]]:=Module[
{generator = OptionValue["generator"],
CalculateGenerators= OptionValue["CalculateGenerators"],
grepmats, gmlg, rotnum,rotK,
yuinput,
mlg,rep,dim,
Unitary,Generator,
Anitunitary,
showSeitz=SpaceGroupIrep`showSeitz,
showMSGSeitz=MSGCorep`showMSGSeitz
},
(*Print[kinfo];*)
mlg=kinfo["Gkin"];
rep=kinfo["rep"];
If[Not@MatrixQ[First[rep]],rep={{#}}&/@rep];
If[ListQ[generator], grepmats = rep[[generator]]; gmlg = mlg[[generator]],
grepmats = rep; gmlg = mlg];
rotnum = Length[gmlg];
(*Print[grepmats];*)
yuinput={{},{}};
Do[
Which[
 Last[gmlg[[iop]]] == 0,
rotK = #[[2]] . {kx, ky, kz} &@(gmlg[[iop]]);
(*Print[gmlg[[iop]],rotK];*)
AppendTo[yuinput[[1]],gmlg[[iop]]->{grepmats[[iop]],rotK(*,directSum[#2,1]&@@gmlg[[iop]]*)}]
,Last[gmlg[[iop]]] ==1,
rotK = #[[2]] . {-kx,- ky,- kz} &@(gmlg[[iop]]);
(*Print[rotK];*)
AppendTo[yuinput[[2]],gmlg[[iop]]->{grepmats[[iop]],rotK(*,directSum[#2,-1]&@@gmlg[[iop]]*)}]
]
,{iop, rotnum}];

yuinput=Association@Thread[{"Unitary","Anitunitary"}->Association/@yuinput];

If[Not@ListQ[generator]&&CalculateGenerators,
Unitary=yuinput["Unitary"];
(*Print[Unitary];*)
Generator=getGenerator[Keys[Unitary][[;;,2]],IdentityMatrix[3],Dot];
Unitary=Association@Select[Normal@Unitary,MemberQ[Generator,Keys[#][[2]]]&];
yuinput["Unitary"]=Unitary;
(*Print[Length[yuinput["Anitunitary"]]\[Equal]0];*)
If[Not@(Length[yuinput["Anitunitary"]]==0),
(*Print@First@Keys[yuinput["Anitunitary"]];*)
yuinput["Anitunitary"]=Association[Normal[yuinput["Anitunitary"]][[{1}]]]]
];
yuinput["Unitary"]=KeyMap[showSeitz[{#[[1]],#[[3]]}]&,yuinput["Unitary"]];
yuinput["Anitunitary"]=KeyMap[showMSGSeitz[{#[[1]],#[[3]],#[[4]]}]&,yuinput["Anitunitary"]];
yuinput["lable"]=kinfo["lable"];
Return[yuinput]
];

Options[TakeSymmetryInfo]={"CartesianCoordinates"->True};
TakeSymmetryInfo[msgno_,k_,number_,OptionsPattern[]]:=Module[
  {CartesianCoordinates=OptionValue["CartesianCoordinates"]},
  Which[
  IntegerQ[msgno],TakeSymmetryInfoSG[msgno,k,number,"CartesianCoordinates"->CartesianCoordinates],
  ListQ[msgno],TakeSymmetryInfoMSG[msgno,k,number,"CartesianCoordinates"->CartesianCoordinates],
  True,Print["Error! Please Check your input."]
  ]

];
Options[interfaceRep] = {"CartesianCoordinates" -> True,"CalculateGenerators"->True};
interfaceRep[msgno_,k_,rep_,OptionsPattern[]]:=Module[{
CartesianCoordinates=OptionValue["CartesianCoordinates"],
CalculateGenerators= OptionValue["CalculateGenerators"],
repinfo,
input,
Ratk
},

Which[StringQ[k],Ratk=k,
ListQ[k],Ratk=Rationalize[k];(*Print[Ratk];*)If[(Or@@(InexactNumberQ/@Ratk)),Print@"Plz do not input inexact number.";Abort[]]];
Which[
NumberQ[msgno],
If[MemberQ[$Packages,"SpaceGroupIrep`"],
  repinfo=TakeSymmetryInfo[msgno,Ratk,rep,"CartesianCoordinates"->CartesianCoordinates];
  
  input=interfacekp[repinfo,"CalculateGenerators"->CalculateGenerators];
  ,Print["Please Import SpaceGroupIrep package"];Abort[]],
ListQ[msgno],
If[MemberQ[$Packages,"MSGCorep`"],
  (*Print[k];*)
  repinfo=TakeSymmetryInfo[msgno,Ratk,rep,"CartesianCoordinates"->CartesianCoordinates];
  (*Print[repinfo];*)
  input=interfacekp[repinfo,"CalculateGenerators"->CalculateGenerators];
  ,Print["Please Import MSGCorep package"];Abort[]]
];
Return[input]
]


(* ::Subsubsection:: *)
(*For SpaceGroupIrep*)


Options[TakeSymmetryInfoSG] = {"CartesianCoordinates" -> True};
TakeSymmetryInfoSG[msgno_, k_, number_, OptionsPattern[]] := Module[
  {Tab, kBZs, Gkin, n, ns, nd, rep, reality, lable,
   CartesianCoordinates = OptionValue["CartesianCoordinates"],
   toCartesianCoordinates,
   selectIr,
   replist, lablelist, realitylist,
   getSGLatt=SpaceGroupIrep`getSGLatt,
   getLGIrepTab=SpaceGroupIrep`getLGIrepTab,
   getRotMatOfK=SpaceGroupIrep`getRotMatOfK,
   BasicVectors=SpaceGroupIrep`BasicVectors
   
   },
  (*If[Not@MemberQ[$Packages,"SpaceGroupIrep`"],Print[
  "Error, This Function need to load SpaceGroupIrep package"];
  Abort[]];*)
  toCartesianCoordinates = 
   FullSimplify[
     Inverse[BasicVectors[getSGLatt[msgno]]] . # . 
      BasicVectors[getSGLatt[msgno]]] &;
  Tab = First@getLGIrepTab[msgno, k];
 (* Print[Tab];*)
  (*Print[SpaceGroupIrep`getLGIrepTab[msgno, k]];*)
  Gkin = Insert[#, 
      If[CartesianCoordinates, 
       toCartesianCoordinates[getRotMatOfK[getSGLatt[msgno], #[[1]]]],
       getRotMatOfK[getSGLatt[msgno], #[[1]]]
       ], 2] & /@ (Append[#, 0] & /@ Tab["Gkin"]);
  kBZs = Tab["kBZs"];
  ns = Length[Tab["sreality"]];
  (*Print[Tab["sreality"]];*)
  nd = Length[Tab["dreality"]];
  selectIr = Which[
     0 < # <= ns, n = #; rep = Tab["sirep"][[n]]; 
     reality = {"Single", Tab["sreality"][[n]]}; 
     lable = Tab["slabel"][[n]],
     ns < # <= nd + ns, n = # - ns; rep = Tab["direp"][[n]]; 
     reality = {"Double", Tab["dreality"][[n]]}; 
     lable = Tab["dlabel"][[n]],
     True, Print["Please input correct irrep number."]; Abort[]] &;
  replist = {}; lablelist = {}; realitylist = {};
  If[ListQ@number,
   Do[
    selectIr[nr];
    replist = AppendTo[replist, rep];
    lablelist = AppendTo[lablelist, lable];
    realitylist = AppendTo[realitylist, reality],
    {nr, number}];
   rep = Fold[directSum, #] & /@ (Transpose[replist]);
   reality = realitylist;
   lable = lablelist;
   ,
   selectIr[number]
   
   ];
  (*Print[MatrixForm/@(Fold[directSum,#]&/@(Transpose[replist]))];*)
  
  Association@
   Thread[{"kBZs", "Gkin", "rep", "reality", "lable"} -> {kBZs, Gkin, 
      rep, reality, lable}]
  ]


(* ::Subsubsection:: *)
(*For MSGCorep*)


Options[TakeSymmetryInfoMSG]={"CartesianCoordinates"->True};
TakeSymmetryInfoMSG[msgno_,k_,number_,OptionsPattern[]]:=Module[
  {Tab,kBZs,Gkin,n,ns,nd,rep,reality,lable,
  CartesianCoordinates=OptionValue["CartesianCoordinates"],
  toCartesianCoordinates,
  directSum=(*FullSimplify@*)ArrayFlatten[{{#1,0},{0,#2}}]&,
  selectIr,replist,lablelist,realitylist,
   getSGLatt=SpaceGroupIrep`getSGLatt,
   getMLGCorep=MSGCorep`getMLGCorep,
   getRotMatOfK=SpaceGroupIrep`getRotMatOfK,
   BasicVectors=SpaceGroupIrep`BasicVectors},
  (*If[Not@MemberQ[$Packages,"SpaceGroupIrep`"],Print["Error, This Function need to load SpaceGroupIrep package"];
Abort[]];*)
toCartesianCoordinates=FullSimplify[Inverse[BasicVectors[getSGLatt[First@msgno]]] . # . BasicVectors[getSGLatt[First@msgno]]]&;
(*Print[k];*)
Tab=getMLGCorep[msgno,k];

Gkin=Insert[#,If[CartesianCoordinates,toCartesianCoordinates[getRotMatOfK[getSGLatt[First@msgno],#[[1]]]],getRotMatOfK[getSGLatt[First@msgno],#[[1]]]],2]&/@(Tab["MLG"]);
kBZs=Tab["k"];
ns=Length[Tab["stype"]];
(*Print[Tab["sreality"]];*)
nd=Length[Tab["dtype"]];
selectIr=Which[0<#<=ns,n=#;rep=Tab["scorep"][[n]];
reality={"Single",Tab["stype"][[n]]};
lable=Tab["slabel"][[n]],ns<#<=nd+ns,n=#-ns;rep=Tab["dcorep"][[n]];
reality={"Double",Tab["dtype"][[n]]};
lable=Tab["dlabel"][[n]],True,Print["Please input correct irrep number."];Abort[]]&;
replist={};lablelist={};realitylist={};
If[ListQ@number,Do[selectIr[nr];
replist=AppendTo[replist,rep];
lablelist=AppendTo[lablelist,lable];
realitylist=AppendTo[realitylist,reality],{nr,number}];
rep=Fold[directSum,#]&/@(Transpose[replist]);
reality=realitylist;
lable=lablelist;,selectIr[number]];
(*Print[MatrixForm/@(Fold[directSum,#]&/@(Transpose[replist]))];*)
Association@Thread[{"kBZs","Gkin","rep","reality","lable"}->{kBZs,Gkin,rep,reality,lable}]]


(* ::Subsection:: *)
(*Interface with MagneticTB *)


(*The current version of interface only work for \[CapitalGamma] point*)
interfaceTB[k_, symmcompile_]:=Module[{input,rotK},
If[k!={0,0,0},Print["Not implemented yet",Abort[]]];
input={{},{}};
Do[
Which[
Last[ns[[2]]] == "F",
rotK = #[[4]] . {kx, ky, kz} &@(ns);
AppendTo[input[[1]],{ns[[3]],rotK}]
,Last[ns[[2]]]=="T",
rotK = #[[4]] . {-kx,- ky,- kz} &@(ns);
AppendTo[input[[2]],{ns[[3]],rotK}]
]
,{ns,symmcompile}];
input
]


(* ::Subsubsection:: *)
(*End*)


End[]
EndPackage[]



