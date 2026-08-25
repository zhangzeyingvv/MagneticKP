(* ::Package:: *)

BeginPackage["MagneticKP`"]
Clear["MagneticKP`*"]


(* ::Subsubsection:: *)
(*Contrcut kp*)


GETkpHam::usage = "GETkpHam[n, input] constructs only the homogeneous order-n part of a k.p Hamiltonian using the legacy calling convention."
IterativelySimplify::usage = "IterativelySimplify[korder, input] constructs a symmetry-allowed k.p Hamiltonian through the default solver."
DirectProductDecomposition::usage = "DirectProductDecomposition[korder, input] constructs a symmetry-allowed k.p Hamiltonian through the direct-product solver."
kpHam::usage = "kpHam[korder, input] constructs a symmetry-allowed k.p Hamiltonian and returns an Association containing the Hamiltonian, orders, dimension, and parameter count."

MagneticKP`e::usage = "e is the coefficient family used for momentum-independent terms in a MagneticKP Hamiltonian."
MagneticKP`t::usage = "t is the coefficient family used for first-order momentum terms in a MagneticKP Hamiltonian."
MagneticKP`r::usage = "r is the coefficient family used for second-order momentum terms in a MagneticKP Hamiltonian."
MagneticKP`s::usage = "s is the coefficient family used for third-order momentum terms in a MagneticKP Hamiltonian."


(* ::Subsubsection:: *)
(*IO*)


(*TakeSymmetryInfoSG::usage = "Init the program"
TakeSymmetryInfoMSG::usage = "Init the program"
interfacekp::usage = "Init the program"
TakeSymmetryInfo::usage = "Init the program"*)
interfaceRep::usage = "interfaceRep[group, k, rep] converts representation data from SpaceGroupIrep or MSGCorep into input accepted by kpHam."


(* ::Subsubsection:: *)
(*Utilities*)


(*getMLGSimpleRepGenerator::usage = "Init the program"
getMLGGenerator::usage = "Init the program"*)
(*GenerateGroup::usage = "Init the program"*)
getGenerator::usage = "getGenerator[elements, identity, multiply] returns a generating subset for a finite group described by its elements and multiplication function."
kx::usage = "kx is the x component of momentum used in MagneticKP Hamiltonians."
ky::usage = "ky is the y component of momentum used in MagneticKP Hamiltonians."
kz::usage = "kz is the z component of momentum used in MagneticKP Hamiltonians."
directSum::usage = "directSum[a, b] returns the block-diagonal direct sum of matrices a and b."
(*texOutput::usage = "Init the program"*)
(*c::usage = "Init the program"*)
(*Min*)
bandManipulate::usage = "bandManipulate[path, n, ham] creates an interactive band plot with controls for the free parameters in ham."

bandplot::usage = "bandplot[path, n, ham, rules] draws the numerical eigenvalues of ham along a labeled reciprocal-fractional path after applying parameter rules."
highSymmetryKPath::usage = "highSymmetryKPath[input] derives symmetry-inequivalent local momentum lines from the complete little magnetic co-group in interfaceRep input and returns a path without plotting it."



(* ::Subsection:: *)
(*Fitting*)


vaspEig::usage = "vaspEig[file, efermi, spin, startband, endband] imports selected bands from a VASP EIGENVAL file and shifts energies by efermi."
bandManipulateWithVASP::usage = "bandManipulateWithVASP[ham, eigdata] interactively compares a k.p Hamiltonian with imported VASP eigenvalues."
fittingKP::usage = "fittingKP[ham, eigdata, krange, initialRules] fits the free parameters of ham to selected reference bands and returns fitted rules, error diagnostics, and a comparison plot."


(* ::Subsection:: *)
(*End*)


EndPackage[]
