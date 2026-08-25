(* Verify generated documentation coverage, links, built pages, and core help
   examples without installing the paclet into the user's repository. *)

$HistoryLength = 0;
ClearAll["Global`*"];

scriptDirectory = DirectoryName[ExpandFileName[$InputFileName]];
testingRoot = DirectoryName[DirectoryName[scriptDirectory]];
sourcePacletRoot = FileNameJoin[{testingRoot, "MagneticKP"}];
builtPacletRoot = FileNameJoin[{testingRoot, "build", "MagneticKP"}];
sourcePaclet = Get[FileNameJoin[{sourcePacletRoot, "PacletInfo.wl"}]];
pacletVersion = sourcePaclet["Version"];
archivePath = FileNameJoin[{
  testingRoot,
  "build",
  "MagneticKP-" <> pacletVersion <> ".paclet"
}];
manualPath = FileNameJoin[{
  testingRoot,
  "Documentation",
  "MagneticKP" <> FromCharacterCode[{20013, 25991, 25163, 20876}] <> ".tex"
}];

expectedSymbols = Sort@Map[
  StringReplace[FileNameTake[#], ".source.wl" ~~ EndOfString -> ""] &,
  FileNames[
    "*.source.wl",
    FileNameJoin[{scriptDirectory, "Sources", "ReferencePages", "Symbols"}]
  ]
];
expectedNotebookCount = 2 (Length[expectedSymbols] + 2);

fail[message_, value_: Missing["NotSupplied"]] := (
  Print[message, If[MissingQ[value], "", InputForm[value]]];
  Exit[1]
);

If[!DirectoryQ[builtPacletRoot] || !FileExistsQ[archivePath],
  fail["Build the MagneticKP paclet before verifying documentation."]
];

(* Only the release-summary entries are ordered: newest version first. This
   check deliberately ignores the positions of every other manual section. *)
If[!FileExistsQ[manualPath], fail["The Chinese manual source is missing."]];
manualText = Import[manualPath, "Text", CharacterEncoding -> "UTF-8"];
manualReleaseVersions = StringCases[
  manualText,
  RegularExpression[
    "\\\\(?:section|subsection)\\{" <>
      FromCharacterCode[{29256, 26412}] <>
      "([0-9]+(?:\\.[0-9]+)*)" <>
      FromCharacterCode[{26356, 26032, 25688, 35201}] <>
      "\\}"
  ] -> "$1"
];
versionKey[value_String] := ToExpression /@ StringSplit[value, "."];
expectedReleaseOrder = Reverse@SortBy[manualReleaseVersions, versionKey];
If[manualReleaseVersions =!= expectedReleaseOrder,
  fail[
    "Chinese manual release summaries must be newest first: ",
    manualReleaseVersions
  ]
];

referenceNames[language_String] := Sort[
  FileBaseName /@ FileNames[
    "*.nb",
    FileNameJoin[{
      builtPacletRoot,
      "Documentation",
      language,
      "ReferencePages",
      "Symbols"
    }]
  ]
];

If[referenceNames["English"] =!= expectedSymbols,
  fail["English symbol page coverage mismatch: ", referenceNames["English"]]
];
If[referenceNames["ChineseSimplified"] =!= expectedSymbols,
  fail[
    "ChineseSimplified symbol page coverage mismatch: ",
    referenceNames["ChineseSimplified"]
  ]
];

authoringNotebooks = FileNames[
  "*.nb",
  FileNameJoin[{sourcePacletRoot, "Documentation"}],
  Infinity
];
builtNotebooks = FileNames[
  "*.nb",
  FileNameJoin[{builtPacletRoot, "Documentation"}],
  Infinity
];
If[
  Length[authoringNotebooks] =!= expectedNotebookCount ||
    Length[builtNotebooks] =!= expectedNotebookCount,
  fail[
    "Unexpected authoring or built notebook count: ",
    {Length[authoringNotebooks], Length[builtNotebooks]}
  ]
];
If[!And @@ (MatchQ[Quiet@Check[Get[#], $Failed], _Notebook] & /@ builtNotebooks),
  fail["At least one built documentation notebook could not be parsed."]
];

(* Guides are navigation pages. Every executable reference page and tutorial
   must carry at least one real saved result in both authoring and built form. *)
runnableNotebooks = Select[
  Join[authoringNotebooks, builtNotebooks],
  !StringContainsQ[#, $PathnameSeparator <> "Guides" <> $PathnameSeparator] &
];
notebookOutputCount[path_String] := Length@Cases[
  Quiet@Check[Get[path], $Failed],
  Cell[_, "Output", ___],
  Infinity
];
outputlessNotebooks = Select[runnableNotebooks, notebookOutputCount[#] < 1 &];
If[outputlessNotebooks =!= {},
  fail["Executable help notebooks without saved Output cells: ", outputlessNotebooks]
];

interfaceRepNotebooks = Select[
  runnableNotebooks,
  FileNameTake[#] === "interfaceRep.nb" &
];

notebookUsesSymbol[path_String, symbol_String] := Module[{notebook},
  notebook = Quiet@Check[Get[path], $Failed];
  !FreeQ[
    notebook,
    value_String /; StringContainsQ[value, symbol <> "["],
    Infinity
  ] || !FreeQ[
    notebook,
    HoldPattern[RowBox[{symbol, "[", ___, "]"}]],
    Infinity
  ]
];

highSymmetryKPathNotebooks = Select[
  runnableNotebooks,
  FileNameTake[#] === "highSymmetryKPath.nb" &
];
If[
  Length[highSymmetryKPathNotebooks] =!= 4 ||
    !And @@ (notebookOutputCount[#] === 2 & /@ highSymmetryKPathNotebooks) ||
    !And @@ (
      notebookUsesSymbol[#, "interfaceRep"] &&
        notebookUsesSymbol[#, "highSymmetryKPath"] &&
        notebookUsesSymbol[#, "bandplot"] & /@
        highSymmetryKPathNotebooks
    ),
  fail[
    "Each authoring and built highSymmetryKPath page must contain only the single K-point example, with the saved path and plot outputs: ",
    highSymmetryKPathNotebooks
  ]
];
If[
  Length[interfaceRepNotebooks] =!= 4 ||
    !And @@ (notebookOutputCount[#] === 3 & /@ interfaceRepNotebooks),
  fail[
    "Each authoring and built interfaceRep page must contain three saved outputs: ",
    notebookOutputCount /@ interfaceRepNotebooks
  ]
];

dualInitializationNames = {
  "interfaceRep.nb",
  "highSymmetryKPath.nb",
  "bandManipulateWithVASP.nb",
  "fittingKP.nb"
};
dualInitializationNotebooks = Select[
  builtNotebooks,
  MemberQ[dualInitializationNames, FileNameTake[#]] &
];
notebookHasDualInitialization[path_String] := Module[
  {notebook, prologs, strings},
  notebook = Quiet@Check[Get[path], $Failed];
  prologs = Cases[
    notebook,
    HoldPattern[CellProlog :> expression_] :> HoldComplete[expression],
    Infinity
  ];
  strings = ToString[#, InputForm] & /@ prologs;
  prologs =!= {} && And @@ Map[
    StringContainsQ[#, "Needs[\"MagneticKP`\"]"] &&
      StringContainsQ[#, "Needs[\"MSGCorep`\"]"] &,
    strings
  ]
];
If[
  Length[dualInitializationNotebooks] =!= 8 ||
    !And @@ (notebookHasDualInitialization /@ dualInitializationNotebooks),
  fail[
    "Built database-backed help inputs must initialize both MagneticKP and MSGCorep: ",
    dualInitializationNotebooks
  ]
];

bandManipulateNotebooks = Select[
  runnableNotebooks,
  FileNameTake[#] === "bandManipulate.nb" &
];
notebookUsesKpHam[path_String] := Module[{notebook, sourceInputStrings},
  notebook = Quiet@Check[Get[path], $Failed];
  sourceInputStrings = Cases[
    notebook,
    Cell[BoxData[code_String], "Input", ___] :> code,
    Infinity
  ];
  AnyTrue[
    sourceInputStrings,
    StringContainsQ[#, "ham = kpHam["] &
  ] || !FreeQ[
    notebook,
    HoldPattern[RowBox[{"kpHam", "[", ___, "]"}]],
    Infinity
  ]
];
If[
  Length[bandManipulateNotebooks] =!= 4 ||
    !And @@ (notebookUsesKpHam /@ bandManipulateNotebooks),
  fail[
    "Each authoring and built bandManipulate page must construct ham with kpHam: ",
    bandManipulateNotebooks
  ]
];

bandManipulateWithVASPNotebooks = Select[
  runnableNotebooks,
  FileNameTake[#] === "bandManipulateWithVASP.nb" &
];
If[
  Length[bandManipulateWithVASPNotebooks] =!= 4 ||
    !And @@ (notebookOutputCount[#] === 4 & /@ bandManipulateWithVASPNotebooks) ||
    !And @@ (
      notebookUsesKpHam[#] &&
        notebookUsesSymbol[#, "interfaceRep"] &&
        notebookUsesSymbol[#, "vaspEig"] &&
        notebookUsesSymbol[#, "bandManipulateWithVASP"] & /@
        bandManipulateWithVASPNotebooks
    ),
  fail[
    "Each authoring and built bandManipulateWithVASP page must import EIGENVAL data and use a kpHam model: ",
    bandManipulateWithVASPNotebooks
  ]
];

bandplotNotebooks = Select[
  runnableNotebooks,
  FileNameTake[#] === "bandplot.nb" &
];
notebookUsesOption[path_String, option_String] := Module[
  {notebook, sourceInputStrings},
  notebook = Quiet@Check[Get[path], $Failed];
  sourceInputStrings = Cases[
    notebook,
    Cell[BoxData[code_String], "Input", ___] :> code,
    Infinity
  ];
  AnyTrue[
    sourceInputStrings,
    StringContainsQ[#, option <> " ->"] &
  ] || Cases[
    notebook,
    RowBox[items_List] /;
      MemberQ[items, option] && MemberQ[items, "->"],
    Infinity
  ] =!= {}
];
bandplotOptionNames = {
  "\"PlotRange\"",
  "\"YTicks\"",
  "FontSize",
  "FontFamily",
  "ImageSize"
};
If[
  Length[bandplotNotebooks] =!= 4 ||
    !And @@ (notebookOutputCount[#] === 6 & /@ bandplotNotebooks) ||
    !And @@ (Function[notebook,
      notebookUsesKpHam[notebook] &&
        notebookUsesSymbol[notebook, "bandplot"] &&
        And @@ (
          notebookUsesOption[notebook, #] & /@ bandplotOptionNames
        )
    ] /@
        bandplotNotebooks
    ),
  fail[
    "Each authoring and built bandplot page must demonstrate every option and contain six saved outputs: ",
    bandplotNotebooks
  ]
];

fittingKPNotebooks = Select[
  runnableNotebooks,
  FileNameTake[#] === "fittingKP.nb" &
];
If[
  Length[fittingKPNotebooks] =!= 4 ||
    !And @@ (notebookOutputCount[#] === 5 & /@ fittingKPNotebooks) ||
    !And @@ (
      notebookUsesKpHam[#] &&
        notebookUsesSymbol[#, "showMLGCorep"] &&
        notebookUsesSymbol[#, "interfaceRep"] &&
        notebookUsesSymbol[#, "vaspEig"] &&
        notebookUsesSymbol[#, "fittingKP"] & /@
        fittingKPNotebooks
    ),
  fail[
    "Each authoring and built fittingKP page must run the original VASP/k.p workflow and contain five saved outputs: ",
    fittingKPNotebooks
  ]
];

privateContextLeaks = Select[
  runnableNotebooks,
  !FreeQ[
    Quiet@Check[Get[#], $Failed],
    value_String /; StringContainsQ[value, "MagneticKP`Private`"],
    Infinity
  ] &
];
If[privateContextLeaks =!= {},
  fail[
    "Saved documentation output exposes the MagneticKP private context: ",
    privateContextLeaks
  ]
];

eigenvalSource = FileNameJoin[{testingRoot, "Examples", "EIGENVAL"}];
eigenvalCopies = Flatten@Table[
  FileNameJoin[{
    root,
    "Documentation",
    language,
    "ReferencePages",
    "Symbols",
    "EIGENVAL"
  }],
  {root, {sourcePacletRoot, builtPacletRoot}},
  {language, {"English", "ChineseSimplified"}}
];
If[
  !FileExistsQ[eigenvalSource] ||
    !And @@ (FileExistsQ[#] &&
      FileHash[#, "SHA256"] === FileHash[eigenvalSource, "SHA256"] & /@ eigenvalCopies),
  fail["The bundled EIGENVAL documentation samples are missing or inconsistent."]
];

sourcePages = Get /@ authoringNotebooks;
pageLinks = DeleteDuplicates@Cases[
  sourcePages,
  HoldPattern[ButtonData -> value_String] :> value,
  Infinity
];

uriTarget[uri_String, language_String] := Module[{relative},
  relative = StringReplace[uri, "paclet:MagneticKP/" -> ""];
  Which[
    StringStartsQ[relative, "ref/"],
      FileNameJoin[{
        builtPacletRoot,
        "Documentation",
        language,
        "ReferencePages",
        "Symbols",
        StringDrop[relative, 4] <> ".nb"
      }],
    StringStartsQ[relative, "guide/"],
      FileNameJoin[{
        builtPacletRoot,
        "Documentation",
        language,
        "Guides",
        StringDrop[relative, 6] <> ".nb"
      }],
    StringStartsQ[relative, "tutorial/"],
      FileNameJoin[{
        builtPacletRoot,
        "Documentation",
        language,
        "Tutorials",
        StringDrop[relative, 9] <> ".nb"
      }],
    True,
      Missing["UnknownURI", uri]
  ]
];

unresolvedLinks = Select[
  Flatten@Table[uriTarget[uri, language], {uri, pageLinks}, {language, {"English", "ChineseSimplified"}}],
  MissingQ[#] || !FileExistsQ[#] &
];
If[unresolvedLinks =!= {}, fail["Unresolved documentation links: ", unresolvedLinks]];

Needs["PacletTools`"];
PacletDirectoryLoad[builtPacletRoot];
builtPaclet = SelectFirst[
  PacletFind["MagneticKP"],
  ExpandFileName[#["Location"]] === ExpandFileName[builtPacletRoot] &,
  Missing["NotFound"]
];
If[MissingQ[builtPaclet], fail["The built MagneticKP paclet was not registered."]];
If[PacletExtensions[builtPaclet, "Documentation"] === {},
  fail["The built paclet has no Documentation extension."]
];

Needs["MagneticKP`"];
actualPublicSymbols = Sort@StringDelete[
  Select[Names["MagneticKP`*"], !StringContainsQ[#, "$"] &],
  StartOfString ~~ "MagneticKP`"
];
intentionallyUndocumentedSymbols = {
  "kx",
  "ky",
  "kz",
  "e",
  "t",
  "r",
  "s",
  "DirectProductDecomposition"
};
expectedPublicSymbols = Sort@Join[
  expectedSymbols,
  intentionallyUndocumentedSymbols
];
If[actualPublicSymbols =!= expectedPublicSymbols,
  fail["Public symbols and reference pages are out of sync: ", actualPublicSymbols]
];
manualInput = <|
  "Unitary" -> <|
    "C2z" -> {{{0, 1}, {-1, 0}}, {-MagneticKP`kx, -MagneticKP`ky, MagneticKP`kz}},
    "Inversion" -> {{{1, 0}, {0, -1}}, {-MagneticKP`kx, -MagneticKP`ky, -MagneticKP`kz}}
  |>
|>;
model = MagneticKP`kpHam[1, manualInput];
If[
  !AssociationQ[model] ||
    !And @@ (KeyExistsQ[model, #] & /@ {"ham", "korder", "dim", "NumberOfParameters"}) ||
    Dimensions[model["ham"]] =!= {2, 2},
  fail["The documented kpHam example did not return the expected result.", model]
];
If[Dimensions[MagneticKP`directSum[IdentityMatrix[2], {{1}}]] =!= {3, 3},
  fail["The documented directSum example failed."]
];
If[
  MagneticKP`getGenerator[Range[0, 3], 0, Mod[#1 + #2, 4] &] === {},
  fail["The documented getGenerator example failed."]
];

Print[InputForm@<|
  "Status" -> "Success",
  "PublicSymbolPagesPerLanguage" -> Length[expectedSymbols],
  "AuthoringNotebooks" -> Length[authoringNotebooks],
  "BuiltNotebooks" -> Length[builtNotebooks],
  "NotebooksWithSavedOutputs" -> Length[runnableNotebooks],
  "InterfaceRepOutputsPerPage" -> 3,
  "HighSymmetryKPathOutputsPerPage" -> 2,
  "BandManipulateUsesKpHam" -> True,
  "BandManipulateWithVASPUsesRawDataAndKpHam" -> True,
  "BandplotUsesKpHam" -> True,
  "BandplotOptionExamplesPerPage" -> Length[bandplotOptionNames],
  "FittingKPUsesOriginalVASPWorkflow" -> True,
  "ManualReleaseOrder" -> manualReleaseVersions,
  "CheckedLinks" -> Length[pageLinks],
  "Archive" -> archivePath
|>];
Exit[0];
