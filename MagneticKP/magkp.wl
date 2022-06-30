(* ::Package:: *)

(* ::Subsection:: *)
(*Begin*)


BeginPackage["MagneticKP`"]


Begin["`Private`"]
installdir=DirectoryName[$InputFileName];
(*installdir="C:\\Users\\zhang\\AppData\\Roaming\\Mathematica\\Applications\\MagneticKP\\";*)
fd=Import[installdir<>"fd.mx"];
If[Not@AssociationQ[fd],fd=Association[]];


(* ::Subsection:: *)
(*Iteratively Simplify Method*)


(* ::Text:: *)
(*\:7528\:7ebf\:6027\:4ee3\:6570\:8bed\:8a00\:91cd\:5199\:ff0c\:7b97\:6cd5\:548c\:539f\:59cb\:7248\:672c\:4e00\:81f4\:ff0c\:7531\:4e8e\:63d0\:524d\:8ba1\:7b97\:4e86F(k)\:6548\:7387\:7565\:6709\:63d0\:9ad8 by Zeying Zhang*)


IterativelySimplify[ikOrder_,input_] := Module[
{
 rotK,kOrder,
fkg,F,tmp,ndim,X,basedim,Gg,S,V,FrobeniusInnerProduct,Gf,Gfc,nullspace,toHam,
symham,symHgk,gHMinusHg,klist,
check,bases,GM,
MRmsymc,
matop,res,np,end,null,
allsymm,
dim
},
GM={{{1,0,0},{0,1,0},{0,0,1}},{{0,-I,0},{I,0,0},{0,0,0}},{{0,0,-I},{0,0,0},{I,0,0}},{{0,0,0},{0,0,-I},{0,I,0}},{{0,1,0},{1,0,0},{0,0,0}},{{1,0,0},{0,-1,0},{0,0,0}},{{0,0,1},{0,0,0},{1,0,0}},{{0,0,0},{0,0,1},{0,1,0}},1/Sqrt[3] {{1,0,0},{0,1,0},{0,0,-2}}};

If[AssociationQ[input],MRmsymc=Values/@(If[MissingQ[#1],Association[],#1]&/@{input["Unitary"],input["Anitunitary"],
input["AntisymmetryUnitaryTest"],input["AntisymmetryAnitunitaryTest"]}),MRmsymc=input];
dim=Length[(MRmsymc/.{{}->Nothing})[[1,-1,-1]]];
If[dim!=3,Print[ToString[dim]<>" Dimensional kp"]];
ndim = Length[(MRmsymc/.{{}->Nothing})[[1,1,1]]];
(*Print[ndim];*)
allsymm={};
Do[
Do[AppendTo[allsymm,Append[i,s]]
,{i,MRmsymc[[s]]}]
,{s,Length[MRmsymc]}];
allsymm=SortBy[allsymm,Not@ComplexRationalMatQ[#[[1]]]&];

(*Print[allsymm];*)
bases[n_]:=Block[{},
Which[n==1,{{{1}}},
n==2,Table[PauliMatrix[i],{i,0,3}],
n==3,GM,
n==4,Flatten[Table[KroneckerProduct[PauliMatrix[i],PauliMatrix[j]],{i,0,3},{j,0,3}],1],
n==6,Flatten[Table[KroneckerProduct[PauliMatrix[i],gm],{i,0,3},{gm,GM}],1],
n==8,Flatten[Table[KroneckerProduct[PauliMatrix[i],PauliMatrix[j],PauliMatrix[k]],{i,0,3},{j,0,3},{k,0,3}],2],
True,{IdentityMatrix[n]}~Join~Normal/@Flatten[Table[SparseArray[{{j,k}->1,{k,j}->1},{n,n}],{k,2,n},{j,1,k-1}],1]~Join~Flatten[Table[SparseArray[{{j,k}->-I,{k,j}->+I},{n,n}],{k,2,n},{j,1,k-1}],1]~Join~Table[Sqrt[2/l/(l+1)] SparseArray[Table[{j,j}->1,{j,1,l}]~Join~{{l+1,l+1}->-l},{n,n}],{l,1,n-1}],
True,{IdentityMatrix[n]}~Join~Table[ SparseArray[Table[{j,j}->1,{j,1,l}]~Join~{{l+1,l+1}->-l},{n,n}],{l,1,n-1}]~Join~Flatten[Table[SparseArray[{{j,k}->1,{k,j}->1},{n,n}],{k,2,n},{j,1,k-1}],1]~Join~Flatten[Table[SparseArray[{{j,k}->-I,{k,j}->+I},{n,n}],{k,2,n},{j,1,k-1}],1];
]

];
X=bases[ndim];
FrobeniusInnerProduct=FullSimplify@Tr[ConjugateTranspose[#1] . #2]&;
toHam[nullspace_,order_]:=
Module[{c,her},
her=X;
(*her[[1]]=0IdentityMatrix[ndim];*)
Total@Table[Subscript[C,order,line]=1;
 Subscript[C,order,line]=.;
 Subscript[C,order,line] Total@(Flatten[her # &/@ klist,1]nullspace[[line]]),{line,Length[nullspace]}]];
Gf[g_,x_]:=g . x . Inverse[g];
Gfc[g_,x_]:=g . Conjugate[x] . Inverse[g];
If[IntegerQ[ikOrder],kOrder=Range[0,ikOrder],kOrder=ikOrder];
symham=0 IdentityMatrix[ndim];
np=0;

Do[
Which[
dim==3,klist=kz^#[[1]]ky^#[[2]]kx^#[[3]]&/@Sort@Flatten[Permutations/@IntegerPartitions[Order,{3},Range[0,Order]],1],
dim==2,klist=ky^#[[1]]kx^#[[2]]&/@Sort@Flatten[Permutations/@IntegerPartitions[Order,{2},Range[0,Order]],1],
dim==1,klist=kz^#[[1]]&/@Sort@Flatten[Permutations/@IntegerPartitions[Order,{1},Range[0,Order]],1],
True,Print["Dimension of kp should be 1,2 or 3"];Abort[]
];
basedim=Length[klist];
(*Print[klist];*)
V=IdentityMatrix[ndim^2 basedim];
Do[
F=fd[{Order,iop[[2]]}];
(*Print[F,MissingQ@fd[{Order,iop[[2]]}]];*)
If[MissingQ[F],
Which[
dim==3,rotK = Thread[{kx, ky, kz} -> iop[[2]]],
dim==2,rotK = Thread[{kx, ky} -> iop[[2]]],
dim==1,rotK = Thread[{kz} -> iop[[2]]]
]

;
(*Print[rotK];*)
fkg=klist/.rotK;
F=Table[tmp[i,j],{i,basedim},{j,basedim}]/.First@SolveAlways[klist ==fkg . Table[tmp[i,j],{i,basedim},{j,basedim}],{kx,ky,kz}];
];
Which[
iop[[3]]==1, Gg=Table[FrobeniusInnerProduct[Gf[iop[[1]],i],j]/FrobeniusInnerProduct[j,j],{j,X},{i,X}];S=(IdentityMatrix[ndim^2 basedim]-KroneckerProduct[F,Gg]);,
iop[[3]]==2, Gg=Table[FrobeniusInnerProduct[Gfc[iop[[1]],i],j]/FrobeniusInnerProduct[j,j],{j,X},{i,X}];S=(IdentityMatrix[ndim^2 basedim]-KroneckerProduct[F,Gg]);,
iop[[3]]==3, Gg=Table[FrobeniusInnerProduct[Gf[iop[[1]],i],j]/FrobeniusInnerProduct[j,j],{j,X},{i,X}];S=(IdentityMatrix[ndim^2 basedim]+KroneckerProduct[F,Gg]);,
iop[[3]]==4, Gg=Table[FrobeniusInnerProduct[Gfc[iop[[1]],i],j]/FrobeniusInnerProduct[j,j],{j,X},{i,X}];S=(IdentityMatrix[ndim^2 basedim]+KroneckerProduct[F,Gg]);,
True, Print["Error!"];Abort[]
];
(*Print[KroneckerProduct[F,Gg]];*)
null=Simplify@NullSpace[S . Transpose@V];
(*Print[Dimensions[null]];*)
(*Print[{S,V,S . Transpose@V,null}];*)
If[null=={},V={};Goto[end]];
V=(*FullSimplify@*)null . V;
(*Print["V:",V]*)
,{iop, allsymm}];
V=FullSimplify[V];
Label[end];
symham+=toHam[V,Order];
np+=Length[V];
,{Order,kOrder}];
(*Print[MatrixForm@symham];*)
Return[Association[{"ham"->symham,"korder"->kOrder,"dim"-> ndim,"NumberOfParameters"->np(*,"NullSpace"\[Rule]V,"kbase"->klist,"base"->X*)}]]
];


(* ::Text:: *)
(*\:539f\:59cb\:7248\:672c by Zhi-Ming Yu*)


GETkpHam[ikOrder_,input_]:=(*GETkpHam[ikOrder,input]=*)
Module[{ik,id,msym,idm,imatr,iTsym,k,g,kg,
imd,tdc,ia,tcoe,
ca,ma,
tkxy,mdr,mdi,Ham,HamConj,dc(*,gx,gy,gz*)
,coea,isov,coed,
coeH={},
coeH2,
isov3,
MRmsymc
},
MRmsymc=Values/@(If[MissingQ[#1],Association[],#1]&/@{input["Unitary"],input["Anitunitary"]});
(*MRmsymc=Values/@{#Unitary,#Anitunitary}&@input;*)
(*Print[MRmsymc];*)
(*Clear[tkxy,mdr,mdi,Ham,HamConj,dc,qx,qy,qz,gx,gy,gz];*)
ia=Length[MRmsymc];  imatr=Table[Length[MRmsymc[[i]]],{i,1,ia}];iTsym=imatr[[2]]; (* MRmsym=DeleteCases[MRmsymc,{}];*)

{id,ik}={Length[MRmsymc[[1,1,1]]],ikOrder};
idm=Length[MRmsymc[[1,1,2]]];  
k=Take[{kx,ky,kz},idm];  imd=If[ik==0,1,Length[Expand[(1+Total[k])^ik]]];g=Take[{kx,ky,kz},idm];

(*Hamiltonian--------------------------------------------------------*)
tkxy[q_]:=Block[{bta,btb,btc},
If[ik==1,
MonomialList[Expand[(1+Total[q])^ik],"NegativeDegreeReverseLexicographic"],
bta=MonomialList[Expand[(1+Total[q])^ik],"NegativeDegreeReverseLexicographic"];
btb=Table[If[i==1,bta[[i]],bta[[i,1]]],{i,1,Length[bta]}];
btc=Table[If[MemberQ[q,btb[[i]]],1,btb[[i]]],{i,1,Length[bta]}];
bta/btc]
];      
(*Print------------------------------------------------------------*)
(* Print["k=", tkxy[k] ]; 
Print["Number of  MR=",imatr ]; *)
(*Print------------------------------------------------------------*)
tdc=Array[dc,{(Sum[i,{i,id}]+Sum[i,{i,id-1}])*imd}];

mdr=
Block[{ta,im},ta=Table[IdentityMatrix[id],{j,1,imd}];
im=0;
Do[Do[im=im+1;ta[[i,ii,ij]]=dc[im],{ii,1,id},{ij,ii,id}],{i,1,imd}];
Do[Do[ta[[j,ij,ii]]=ta[[j,ii,ij]],{ii,1,id},{ij,ii,id}],{j,1,imd}];
ta];

mdi=
Block[{ta,im},ta=Table[IdentityMatrix[id],{j,1,imd}];
im=Sum[i,{i,id}]*imd;
Do[Do[If[ii==ij,ta[[i,ii,ij]]=0,im=im+1;ta[[i,ii,ij]]=I*dc[im]],{ii,1,id},{ij,ii,id}],{i,1,imd}];

Do[Do[ta[[j,ij,ii]]=- ta[[j,ii,ij]],{ii,1,id},{ij,ii,id}],{j,1,imd}];
ta];

Ham[Q_]:=Block[{ta,tqx,tqy,tqz,tq},tq=Take[{tqx,tqy,tqz},Length[Q]];ta=tkxy[tq];  
{tqx,tqy,tqz}=PadRight[Q,3];ta . (mdr+mdi) ];
HamConj[Q_]:=Block[{ta,tqx,tqy,tqz,tq},tq=Take[{tqx,tqy,tqz},Length[Q]];ta=tkxy[tq];  
{tqx,tqy,tqz}=PadRight[Q,3];ta . (mdr-mdi) ];
(*Hamiltonian--------------------------------------------------------*)

(*symmetric constraints----------------------------------------------*)
Off[Solve::svars];  
msym=MRmsymc[[1]];
Do[
Clear[ma,coea,isov,coed];
coeH={};
Do[ ca=tdc[[ia]];If[Length[ca]==1 ,coeH=Join[coeH,{ca}]  ]          ,{ia,1,Length[tdc]}];
coeH2=DeleteDuplicates[coeH];  (*Print["coeH2=",coeH2];*)

kg=msym[[i,2]];
(*ma=msym\[LeftDoubleBracket]i,1\[RightDoubleBracket].Ham[kg].Inverse[msym\[LeftDoubleBracket]i,1\[RightDoubleBracket]]-Ham[k];*)
ma=msym[[i,1]] . Ham[k] . Inverse[msym[[i,1]]]-Ham[kg];

coea=Simplify@Flatten[CoefficientList[ma,k]];(*Print["coea=",coea]; *)
isov=Solve[coea==Table[0,{i,1,Length[coea]}],coeH2];    coed=FullSimplify[tdc/.isov[[1]]];   (*Print["isov=",isov];*)  (*Abort[];*)

Do[dc[ia]=coed[[ia]],{ia,1,Length[tdc]}];

,{i,1,Length[msym]}];

(*time reversal symmetry----------------------------------------------*)
If[iTsym!=0,msym=MRmsymc[[2]];

Do[
Clear[ma,coea,isov,coed];
coeH={};
Do[ ca=tdc[[ia]];If[Length[ca]==1 ,coeH=Join[coeH,{ca}]  ]          ,{ia,1,Length[tdc]}];
coeH2=DeleteDuplicates[coeH];  (*Print["coeH2=",coeH2];*)

kg=msym[[i,2]];
(*ma=msym\[LeftDoubleBracket]i,1\[RightDoubleBracket].HamConj[kg].Inverse[msym\[LeftDoubleBracket]i,1\[RightDoubleBracket]]-Ham[k];*)
ma=msym[[i,1]] . HamConj[k] . Inverse[msym[[i,1]]]-Ham[kg];

coea=Simplify@Flatten[CoefficientList[ma,k]];
isov=Solve[coea==Table[0,{i,1,Length[coea]}],coeH2];    coed=FullSimplify[tdc/.isov[[1]]];    (*Print["isov=",isov]; Abort[];*)

Do[dc[ia]=coed[[ia]],{ia,1,Length[tdc]}];

,{i,1,Length[msym]}];
];

Clear[ma,coea,isov,coed(*,c*)];
coeH={};
Do[ ca=tdc[[ia]];If[Length[ca]==1 ,coeH=Join[coeH,{ca}]  ]          ,{ia,1,Length[tdc]}];
coeH2=DeleteDuplicates[coeH];(*Print["Princ Coefficiens=",coeH2];*)
(*isov3=Solve[coeH2==Table[Subscript[c, i],{i,1,Length[coeH2]}],coeH2];*)
isov3=Solve[coeH2==Table[C[i],{i,1,Length[coeH2]}],coeH2];
(*Print["ham is over",MatrixForm[Ham[k]/.isov3\[LeftDoubleBracket]1\[RightDoubleBracket]]]; Abort[];*)

Return[Association[{"ham"->ToRadicals[(*FullSimplify[*)Ham[k]/.isov3[[1](*]*)]]},"korder"->ikOrder,"dim"->id,"NumberOfParameters"->Length[coeH2]]]
];



(* ::Subsection:: *)
(*Direct Product Decomposition Method*)


(* ::Text:: *)
(* D. Gresch, Identifying Topological Semimetals, Ph.D. thesis (ETH Zurich) (2018) *)


DirectProductDecomposition[ikOrder_,input_] := Module[
{
 rotK,kOrder,
fkg,Fg,tmp,ndim,X,basedim,Gg,FrobeniusInnerProduct,Gf,Gfc,
getIntersection,nullspace,toHam,
symham,symHgk,gHMinusHg,klist,toCartesianCoordinates,
check,bases,GM,
MRmsymc,
matop,np
},
(*klist=kz^#[[1]]ky^#[[2]]kx^#[[3]]&/@SortBy[Values[Solve[{a+b+c<=ikOrder},{a,b,c},NonNegativeIntegers]],{Total[#]}&];*)
If[IntegerQ[ikOrder],kOrder=Range[0,ikOrder],kOrder=ikOrder];


MRmsymc=Values/@(If[MissingQ[#1],Association[],#1]&/@{input["Unitary"],input["Anitunitary"],input["AntisymmetryUnitaryTest"],input["AntisymmetryAnitunitaryTest"]});
(*Print[MRmsymc];*)
ndim = Length[MRmsymc[[1,1,1]]];
(*grepmats=N@grepmats;*)
(*tmpre[i_,j_]=RandomInteger[{-10,10}];
tmpim[i_,j_]=RandomInteger[{-10,10}];*)

GM={{{1,0,0},{0,1,0},{0,0,1}},{{0,-I,0},{I,0,0},{0,0,0}},{{0,0,-I},{0,0,0},{I,0,0}},{{0,0,0},{0,0,-I},{0,I,0}},{{0,1,0},{1,0,0},{0,0,0}},{{1,0,0},{0,-1,0},{0,0,0}},{{0,0,1},{0,0,0},{1,0,0}},{{0,0,0},{0,0,1},{0,1,0}},1/Sqrt[3] {{1,0,0},{0,1,0},{0,0,-2}}};

(*X=Flatten[Join[Table[tmp=Table[0,ndim,ndim];tmp[[i,j]]=1;tmp[[j,i]]=1;tmp,{i,ndim},{j,i}],
Table[tmp=Table[0,ndim,ndim];tmp[[i,j]]=-I;tmp[[j,i]]=I;tmp,{i,2,ndim},{j,i-1}]],1];*)
bases[n_]:=Block[{},
Which[n==1,{{1}},
n==2,Table[PauliMatrix[i],{i,0,3}],
n==3,GM,
n==4,Flatten[Table[KroneckerProduct[PauliMatrix[i],PauliMatrix[j]],{i,0,3},{j,0,3}],1],
n==6,Flatten[Table[KroneckerProduct[PauliMatrix[i],gm],{i,0,3},{gm,GM}],1],
n==8,Flatten[Table[KroneckerProduct[PauliMatrix[i],PauliMatrix[j],PauliMatrix[k]],{i,0,3},{j,0,3},{k,0,3}],2],
True,{IdentityMatrix[n]}~Join~Table[ SparseArray[Table[{j,j}->1,{j,1,l}]~Join~{{l+1,l+1}->-l},{n,n}],{l,1,n-1}]~Join~Flatten[Table[SparseArray[{{j,k}->1,{k,j}->1},{n,n}],{k,2,n},{j,1,k-1}],1]~Join~Flatten[Table[SparseArray[{{j,k}->-I,{k,j}->+I},{n,n}],{k,2,n},{j,1,k-1}],1];
]

];
X=bases[ndim];

FrobeniusInnerProduct=FullSimplify@Tr[ConjugateTranspose[#1] . #2]&;
getIntersection[l1_,l2_]:=Module[{n,ker,coeffs,res},
If[l1=={}||l2=={},Return[{}]];
ker=Simplify@NullSpace[Transpose[Join[l1,l2]]];
n=Length[l1];
coeffs=Map[Function[v,v[[1;;n]]],ker];
(*Print[coeffs];*)
res=Simplify@Map[Function[v,v . l1],coeffs](*;
#/First@Cases[#,_?(#\[NotEqual]0&),1,1]&/@res*)
];

toHam[nullspace_,order_]:=Block[{c},Total@Table[Subscript[C,order,line] Total@(Flatten[X # &/@ klist,1]nullspace[[line]]),{line,Length[nullspace]}]];
(*Print[{ndim,basedim}];*)

Gf[g_,x_]:=g . x . Inverse[g];
Gfc[g_,x_]:=g . Conjugate[x] . Inverse[g];
symham=0 IdentityMatrix[ndim];
np=0;
Do[
(*klist=kz^#[[1]]ky^#[[2]]kx^#[[3]]&/@SortBy[Values[Solve[Or@@Table[a+b+c==i,{i,#}],{a,b,c},NonNegativeIntegers]&[{Order}]],{Total[#]}&];*)
klist=kz^#[[1]]ky^#[[2]]kx^#[[3]]&/@Sort@Flatten[Permutations/@IntegerPartitions[Order,{3},Range[0,Order]],1];
basedim=Length[klist];
nullspace=IdentityMatrix[ndim^2 basedim];
Do[
Fg=fd[{Order,iop[[2]]}];
If[MissingQ[Fg],
rotK = Thread[{kx, ky, kz} -> iop[[2]]];
fkg=klist/.rotK;
Fg=Table[tmp[i,j],{i,basedim},{j,basedim}]/.First@SolveAlways[klist ==fkg . Table[tmp[i,j],{i,basedim},{j,basedim}],{kx,ky,kz}];
];

Gg=Chop@Table[FrobeniusInnerProduct[Gf[iop[[1]],i],j]/FrobeniusInnerProduct[j,j],{j,X},{i,X}];

(*Print[Dimensions@(Simplify@(KroneckerProduct[Fg,Gg]-IdentityMatrix[ndim^2 basedim]))];*)
nullspace=getIntersection[Simplify@NullSpace[Simplify@(KroneckerProduct[Fg,Gg]-IdentityMatrix[ndim^2 basedim])],nullspace];

nullspace=FullSimplify[nullspace];
(*Print[Dimensions[nullspace]];*)
(*Print[nullspace];*)
,{iop, First[MRmsymc]}];

Do[
(*matop=Normal[CoefficientArrays[iop[[2]],{kx,ky,kz}][[2]]];
(*Print[gmlg[[iop]]];*)
rotK = Thread[{kx, ky, kz} -> Inverse[matop] . {kx, ky, kz}];
(*Print[gmlg[[iop]]];*)
(*rotK = Thread[{kx, ky, kz} -> iop[[2]]];*)
fkg=klist/.rotK;
(*Print[Transpose@#[[2]] . {kx, ky, kz} &@(gmlg[[iop]])];*)
Fg=Table[tmp[i,j],{i,basedim},{j,basedim}]/.First@SolveAlways[klist . Table[tmp[i,j],{i,basedim},{j,basedim}]==fkg,{kx,ky,kz}];*)
Fg=fd[{Order,iop[[2]]}];
If[MissingQ[Fg],
rotK = Thread[{kx, ky, kz} -> iop[[2]]];
fkg=klist/.rotK;
Fg=Table[tmp[i,j],{i,basedim},{j,basedim}]/.First@SolveAlways[klist ==fkg . Table[tmp[i,j],{i,basedim},{j,basedim}],{kx,ky,kz}];
];


Gg=Chop@Table[FrobeniusInnerProduct[Gfc[iop[[1]],i],j]/FrobeniusInnerProduct[j,j],{j,X},{i,X}];
(*Print[Dimensions@(Simplify@(KroneckerProduct[Fg,Gg]-IdentityMatrix[ndim^2 basedim]))];*)
nullspace=getIntersection[Simplify@NullSpace[Simplify@(KroneckerProduct[Fg,Gg]-IdentityMatrix[ndim^2 basedim])],nullspace];
nullspace=FullSimplify[nullspace];
(*Print[Dimensions@nullspace];*)
,{iop, MRmsymc[[2]]}];


Do[
Fg=fd[{Order,iop[[2]]}];
If[MissingQ[Fg],
rotK = Thread[{kx, ky, kz} -> iop[[2]]];
fkg=klist/.rotK;
Fg=Table[tmp[i,j],{i,basedim},{j,basedim}]/.First@SolveAlways[klist ==fkg . Table[tmp[i,j],{i,basedim},{j,basedim}],{kx,ky,kz}];
];


(*Print[Dimensions@(Simplify@(KroneckerProduct[Fg,Gg]-IdentityMatrix[ndim^2 basedim]))];*)
nullspace=getIntersection[Simplify@NullSpace[Simplify@(KroneckerProduct[Fg,Gg]+IdentityMatrix[ndim^2 basedim])],nullspace];
nullspace=FullSimplify[nullspace];
(*Print[Dimensions@nullspace];*)
,{iop, MRmsymc[[3]]}];


Do[
Fg=fd[{Order,iop[[2]]}];
If[MissingQ[Fg],
rotK = Thread[{kx, ky, kz} -> iop[[2]]];
fkg=klist/.rotK;
Fg=Table[tmp[i,j],{i,basedim},{j,basedim}]/.First@SolveAlways[klist ==fkg . Table[tmp[i,j],{i,basedim},{j,basedim}],{kx,ky,kz}];
];


Gg=Chop@Table[FrobeniusInnerProduct[Gfc[iop[[1]],i],j]/FrobeniusInnerProduct[j,j],{j,X},{i,X}];
(*Print[Dimensions@(Simplify@(KroneckerProduct[Fg,Gg]-IdentityMatrix[ndim^2 basedim]))];*)
nullspace=getIntersection[Simplify@NullSpace[Simplify@(KroneckerProduct[Fg,Gg]+IdentityMatrix[ndim^2 basedim])],nullspace];
nullspace=FullSimplify[nullspace];
(*Print[Dimensions@nullspace];*)
,{iop, MRmsymc[[4]]}];


(*nullspace=(#/First@Cases[#,_?(#!=0&),1,1]&/@nullspace);
nullspace=SortBy[#,FirstPosition[#,_?(#!=0&)]&]&[nullspace];*)
(*Print[Dimensions[nullspace],MatrixForm@nullspace];*)
(*symham=Chop[toHam[nullspace]];*)
symham+=toHam[nullspace,Order];
np+=Length[nullspace];
,{Order,kOrder}];

(*Print[MatrixForm@symham];*)
Return[Association[{"ham"->symham,"korder"->kOrder,"dim"-> ndim,"NumberOfParameters"->np(*,"NullSpace"->nullspace,"kbase"->klist,"base"->X*)}]]
];



(* ::Subsection:: *)
(*kp Ham*)


Options[kpHam] = {"Method" -> "IterativeSimplification"};
kpHam[ikOrder_,input_,OptionsPattern[]]:=Which[
OptionValue["Method"]=="IterativeSimplification2016"&&NumberQ[ikOrder]&&(MissingQ[input["AntisymmetryUnitaryTest"]])&&MissingQ[input["AntisymmetryAnitunitaryTest"]],
GETkpHam[ikOrder,input],
OptionValue["Method"]=="IterativeSimplification",IterativelySimplify[ikOrder,input],
OptionValue["Method"]=="DirectProductDecomposition",DirectProductDecomposition[ikOrder,input],
True,Print["Plz chk your inputs"];Abort[]
]


(* ::Subsection::Closed:: *)
(*End*)


End[]


EndPackage[]
