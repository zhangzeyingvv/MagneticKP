(* Generate simplified-Chinese authoring sources from the canonical English
   MagneticKP sources. Code strings, symbols, contexts, and URIs are kept. *)

ClearAll["Global`*"];

scriptDirectory = DirectoryName[ExpandFileName[$InputFileName]];
projectRoot = DirectoryName[DirectoryName[scriptDirectory]];
documentationRoot = FileNameJoin[{
  projectRoot,
  "MagneticKP",
  "Documentation",
  "ChineseSimplified"
}];
englishSourceRoot = FileNameJoin[{scriptDirectory, "Sources"}];
englishExpectedRoot = FileNameJoin[{scriptDirectory, "Expected"}];
sourceRoot = FileNameJoin[{scriptDirectory, "SourcesChineseSimplified"}];
expectedRoot = FileNameJoin[{scriptDirectory, "ExpectedChineseSimplified"}];

translations = Get[FileNameJoin[{scriptDirectory, "ChineseTranslations.wl"}]];
If[!AssociationQ[translations],
  Print["ChineseTranslations.wl must evaluate to an Association."];
  Exit[2]
];

undocumentedHelpSymbols = {
  "kx",
  "ky",
  "kz",
  "e",
  "t",
  "r",
  "s",
  "DirectProductDecomposition"
};
Scan[
  Function[name,
    Scan[
      If[FileExistsQ[#], DeleteFile[#]] &,
      {
        FileNameJoin[{
          sourceRoot,
          "ReferencePages",
          "Symbols",
          name <> ".source.wl"
        }],
        FileNameJoin[{
          expectedRoot,
          "ReferencePages",
          "Symbols",
          name <> ".expected.txt"
        }],
        FileNameJoin[{
          documentationRoot,
          "ReferencePages",
          "Symbols",
          name <> ".nb"
        }]
      }
    ]
  ],
  undocumentedHelpSymbols
];

(* The translated vaspEig page uses the same real sample as the English page.
   Keep it beside the generated notebook so NotebookDirectory[] is portable. *)
eigenvalSource = FileNameJoin[{projectRoot, "Examples", "EIGENVAL"}];
eigenvalTarget = FileNameJoin[{
  documentationRoot,
  "ReferencePages",
  "Symbols",
  "EIGENVAL"
}];
If[!FileExistsQ[eigenvalSource],
  Print["Missing documentation EIGENVAL sample: ", eigenvalSource];
  Exit[3]
];
If[!DirectoryQ[DirectoryName[eigenvalTarget]],
  CreateDirectory[DirectoryName[eigenvalTarget], CreateIntermediateDirectories -> True]
];
If[
  !FileExistsQ[eigenvalTarget] ||
    FileHash[eigenvalSource, "SHA256"] =!= FileHash[eigenvalTarget, "SHA256"],
  CopyFile[eigenvalSource, eigenvalTarget, OverwriteTarget -> True]
];

writeUTF8TextIfChanged[path_String, text_String] := Module[
  {existing, normalized},
  normalized = StringReplace[
    text,
    RegularExpression["(?m)[ \\t]+$"] -> ""
  ];
  If[!DirectoryQ[DirectoryName[path]],
    CreateDirectory[DirectoryName[path], CreateIntermediateDirectories -> True]
  ];
  existing = If[
    FileExistsQ[path],
    FromCharacterCode[BinaryReadList[path], "UTF8"],
    Missing["NotFound"]
  ];
  If[existing =!= normalized,
    Export[path, normalized, "Text", CharacterEncoding -> "UTF-8"]
  ];
];

translateString[value_String] := Lookup[translations, value, value];

(* Translate visible authoring content only. Cell styles, CellLabel metadata,
   ButtonData URIs, option names/defaults, and executable BoxData strings must
   remain byte-for-byte compatible with the English authoring page. *)
ClearAll[
  translateInline,
  translateEmbeddedCells,
  translateCellContent,
  translateCell,
  translateNotebookOption,
  translateNotebook
];

translateInline[value_String] := translateString[value];
translateInline[value_List] := translateInline /@ value;
translateInline[ButtonBox[label_, options___]] := ButtonBox[
  If[StringQ[label], translateString[label], translateInline[label]],
  options
];
translateInline[value_Cell] := translateCell[value];
translateInline[value_] := value;

translateEmbeddedCells[value_] := value /. cell_Cell :> translateCell[cell];

translateCellContent[value_String] := translateString[value];
translateCellContent[TextData[value_]] := TextData[translateInline[value]];
translateCellContent[BoxData[value_]] := BoxData[translateEmbeddedCells[value]];
translateCellContent[value_] := value;

translateCell[Cell[CellGroupData[cells_List, state_], options___]] := Cell[
  CellGroupData[translateCell /@ cells, state],
  options
];
translateCell[Cell[content_, style_String, options___]] := Cell[
  translateCellContent[content],
  style,
  options
];
translateCell[value_] := value;

translateNotebookOption[WindowTitle -> value_String] :=
  WindowTitle -> translateString[value];
translateNotebookOption[value_] := value;

translateNotebook[Notebook[cells_List, options___]] := Notebook[
  translateCell /@ cells,
  Sequence @@ (translateNotebookOption /@ {options})
];

englishSources = Sort@FileNames["*.source.wl", englishSourceRoot, Infinity];
If[englishSources === {},
  Print["No English documentation sources were found."];
  Exit[4]
];

Scan[
  Function[englishPath,
    relativePath = FileNameDrop[englishPath, FileNameDepth[englishSourceRoot]];
    relativeBase = StringReplace[relativePath, ".source.wl" ~~ EndOfString -> ""];
    notebook = Get[englishPath];
    If[Head[notebook] =!= Notebook,
      Print["Invalid English Notebook source: ", englishPath];
      Exit[5]
    ];
    translated = translateNotebook[notebook];
    writeUTF8TextIfChanged[
      FileNameJoin[{sourceRoot, relativeBase <> ".source.wl"}],
      ToString[translated, InputForm, PageWidth -> 120]
    ];

    englishExpectedPath = FileNameJoin[{englishExpectedRoot, relativeBase <> ".expected.txt"}];
    expectedLines = Select[
      StringSplit[Import[englishExpectedPath, "Text", CharacterEncoding -> "UTF-8"], "\n"],
      StringLength[#] > 0 &
    ];
    translatedExpected = translateString /@ expectedLines;
    writeUTF8TextIfChanged[
      FileNameJoin[{expectedRoot, relativeBase <> ".expected.txt"}],
      StringRiffle[translatedExpected, "\n"] <> "\n"
    ];
    Print[FileNameJoin[{sourceRoot, relativeBase <> ".source.wl"}]]
  ],
  englishSources
];

Print[
  "Generated simplified-Chinese MagneticKP documentation sources: ",
  Length[englishSources],
  " pages."
];
Exit[0];
