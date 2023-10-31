(* ::Package:: *)

BeginPackage["MagneticKP`"]
Clear["MagneticKP`*"]


(* ::Subsubsection:: *)
(*Contrcut kp*)


GETkpHam::usage = "Old version of iteratively simplify algorithm"
IterativelySimplify::usage = "Iteratively simplify algorithm"
DirectProductDecomposition::usage = "Direct product decomposition algorithm"
kpHam::usage = "Get the kp Hamiltonian"


(* ::Subsubsection:: *)
(*IO*)


(*TakeSymmetryInfoSG::usage = "Init the program"
TakeSymmetryInfoMSG::usage = "Init the program"
interfacekp::usage = "Init the program"
TakeSymmetryInfo::usage = "Init the program"*)
interfaceRep::usage = "Interface for SpaceGroupIrep and MSGCorep"


(* ::Subsubsection:: *)
(*Utilities*)


(*getMLGSimpleRepGenerator::usage = "Init the program"
getMLGGenerator::usage = "Init the program"*)
(*GenerateGroup::usage = "Init the program"*)
getGenerator::usage = "Use greedy algorithm to get the generator of a group, may not give the minima generator set"
kx::usage = "kx"
ky::usage = "ky"
kz::usage = "kz"
directSum::usage = "directSum"
(*texOutput::usage = "Init the program"*)
(*c::usage = "Init the program"*)
(*Min*)
bandManipulate::usage = "directSum"

bandplot::usage = "bandplot"



(* ::Subsection:: *)
(*Fitting*)


vaspEig::usage = "Read the VASP EIGENVAL file"
bandManipulateWithVASP::usage = ""
fittingKP::usage = ""


(* ::Subsection:: *)
(*End*)


EndPackage[]
