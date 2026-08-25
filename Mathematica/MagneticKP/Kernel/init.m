(* ::Package:: *)

(* Mathematica Init File *)
With[
  {
    packageRoot = DirectoryName[$InputFileName]
  },
  Get[FileNameJoin[{packageRoot, "Usage.wl"}]];
  Get[FileNameJoin[{packageRoot, "Utilities.wl"}]];
  Get[FileNameJoin[{packageRoot, "CoreAlgorithms.wl"}]];
  Get[FileNameJoin[{packageRoot, "magkp.wl"}]];
  Get[FileNameJoin[{packageRoot, "IO.m"}]];
  Get[FileNameJoin[{packageRoot, "HighSymmetryKPath.wl"}]];
  Get[FileNameJoin[{packageRoot, "Fitting.wl"}]]
]
(* Get[FileNameJoin[{packageRoot, "plot.wl"}]] *)

