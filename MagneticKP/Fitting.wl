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


fittingKP[h_,eigdata_,krange_,initparms_]:=Module[
{params,mparams,klist,hp,modelhams,nk,objparams,res,lsq},
params=Complement[Variables[h],{kx,ky,kz}];
klist=eigdata[[;;,1]];
nk=Length[krange];
mparams=ToExpression[StringRiffle[#,"c"]&/@(params/.{Subscript->List})];
hp=h/.Thread[params->mparams];
modelhams=Table[hp/.{kx->k[[1]],ky->k[[2]],kz->k[[3]]},{k,klist[[krange]]}];
res=FindMinimum[Total@Table[EuclideanDistance[Sort@Eigenvalues[modelhams[[i]]],eigdata[[;;,2]][[krange[[i]]]]],{i,nk}],Transpose[{mparams,params/.initparms}]];
objparams=res[[2]]/.Thread[mparams->params];

lsq={Total@Table[EuclideanDistance[Sort@Eigenvalues[modelhams[[i]]],eigdata[[;;,2]][[krange[[i]]]]],{i,nk}]/.(initparms/.Thread[params->mparams]),res[[1]]};
Association["FittedParams"->objparams,"LSQForInitParams"->lsq[[1]],"LSQForFittedParams"->lsq[[2]]]
]


End[]
EndPackage[]

