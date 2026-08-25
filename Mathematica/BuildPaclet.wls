(* Build a bilingual MagneticKP paclet with the official PacletTools pipeline.

   Mathematica 14 builds one authoring language tree at a time. English is
   built by PacletBuild; ChineseSimplified is built in an isolated
   documentation-only stage and then installed into the final paclet tree.
*)

$HistoryLength = 0;
ClearAll["Global`*"];

Needs["PacletTools`"];

testingRoot = ExpandFileName[DirectoryName[$InputFileName]];
sourceDirectory = FileNameJoin[{testingRoot, "MagneticKP"}];
documentationDirectory = FileNameJoin[{sourceDirectory, "Documentation"}];
englishDocumentation = FileNameJoin[{documentationDirectory, "English"}];
chineseDocumentation = FileNameJoin[{documentationDirectory, "ChineseSimplified"}];
buildDirectory = FileNameJoin[{testingRoot, "build"}];
stageDirectory = FileNameJoin[{buildDirectory, "PacletSource"}];
chineseStageDirectory = FileNameJoin[{buildDirectory, "ChineseDocumentationSource"}];
chineseBuildDirectory = FileNameJoin[{buildDirectory, "ChineseDocumentationBuild"}];
builtPacletDirectory = FileNameJoin[{buildDirectory, "MagneticKP"}];

requiredPaths = {
  FileNameJoin[{sourceDirectory, "PacletInfo.wl"}],
  FileNameJoin[{sourceDirectory, "Kernel"}],
  englishDocumentation,
  chineseDocumentation
};
If[!And @@ (FileExistsQ /@ requiredPaths),
  Print["Missing MagneticKP build input: ", Select[requiredPaths, Not@*FileExistsQ]];
  Exit[2]
];

pacletObject = Get[FileNameJoin[{sourceDirectory, "PacletInfo.wl"}]];
If[!MatchQ[pacletObject, _PacletObject],
  Print["MagneticKP/PacletInfo.wl did not evaluate to a PacletObject."];
  Exit[3]
];
pacletVersion = pacletObject["Version"];

englishPageCount = Length@FileNames["*.nb", englishDocumentation, Infinity];
chinesePageCount = Length@FileNames["*.nb", chineseDocumentation, Infinity];
If[englishPageCount === 0 || chinesePageCount =!= englishPageCount,
  Print[
    "The bilingual authoring trees are incomplete or out of sync: ",
    {englishPageCount, chinesePageCount}
  ];
  Exit[4]
];

If[!DirectoryQ[buildDirectory],
  CreateDirectory[buildDirectory, CreateIntermediateDirectories -> True]
];
Scan[
  If[DirectoryQ[#], DeleteDirectory[#, DeleteContents -> True]] &,
  {
    stageDirectory,
    chineseStageDirectory,
    chineseBuildDirectory,
    builtPacletDirectory
  }
];
archiveCandidate = FileNameJoin[{
  buildDirectory,
  "MagneticKP-" <> pacletVersion <> ".paclet"
}];
If[FileExistsQ[archiveCandidate], DeleteFile[archiveCandidate]];

CopyDirectory[sourceDirectory, stageDirectory];
DeleteDirectory[
  FileNameJoin[{stageDirectory, "Documentation", "ChineseSimplified"}],
  DeleteContents -> True
];
Scan[DeleteFile, FileNames[".DS_Store", stageDirectory, Infinity]];

(* Avoid an irrelevant WolframLanguageData lookup for third-party symbols
   while retaining the official documentation skeletonizer and search index. *)
Needs["DocumentationSearch`Skeletonizer`"];
SetOptions[
  DocumentationSearch`Skeletonizer`Skeletonize,
  "UseWolframLanguageData" -> False
];
SetOptions[
  DocumentationSearch`Skeletonizer`Private`skeletonizeFiles,
  "UseWolframLanguageData" -> False
];

Print["[1/4] Building package and English documentation with PacletBuild..."];
{buildResult, buildMessages} = Block[{$MessageList = {}},
  {
    PacletTools`PacletBuild[stageDirectory, buildDirectory],
    $MessageList
  }
];
If[buildMessages =!= {},
  Print["Unexpected messages during PacletBuild: ", buildMessages];
  Exit[5]
];
If[!MatchQ[buildResult, _Success],
  Print[InputForm[buildResult]];
  Exit[6]
];

Print["[2/4] Building ChineseSimplified documentation..."];
CreateDirectory[
  FileNameJoin[{chineseStageDirectory, "Documentation"}],
  CreateIntermediateDirectories -> True
];
CopyDirectory[
  chineseDocumentation,
  FileNameJoin[{chineseStageDirectory, "Documentation", "English"}]
];
Put[
  PacletObject[<|
    "Name" -> "MagneticKP",
    "Version" -> pacletVersion,
    "WolframVersion" -> "12.1+",
    "Extensions" -> {
      {
        "Documentation",
        "LinkBase" -> "MagneticKP",
        "MainPage" -> "Guides/MagneticKP"
      }
    }
  |>],
  FileNameJoin[{chineseStageDirectory, "PacletInfo.wl"}]
];
{chineseBuildResult, chineseBuildMessages} = Block[{$MessageList = {}},
  {
    PacletTools`PacletDocumentationBuild[
      chineseStageDirectory,
      chineseBuildDirectory
    ],
    $MessageList
  }
];
If[chineseBuildMessages =!= {},
  Print[
    "Unexpected messages during Chinese documentation build: ",
    chineseBuildMessages
  ];
  Exit[7]
];
If[!MatchQ[chineseBuildResult, _Success],
  Print[InputForm[chineseBuildResult]];
  Exit[8]
];

builtChineseEnglish = FileNameJoin[{
  chineseBuildDirectory,
  "MagneticKP",
  "Documentation",
  "English"
}];
builtChineseTarget = FileNameJoin[{
  builtPacletDirectory,
  "Documentation",
  "ChineseSimplified"
}];
If[!DirectoryQ[builtChineseEnglish],
  Print["The Chinese documentation build directory was not created."];
  Exit[9]
];
If[DirectoryQ[builtChineseTarget],
  DeleteDirectory[builtChineseTarget, DeleteContents -> True]
];
CopyDirectory[builtChineseEnglish, builtChineseTarget];

(* PacletDocumentationBuild copies generated pages and indices but omits
   arbitrary files stored beside authoring notebooks. Restore the real VASP
   sample explicitly so the documented NotebookDirectory[] example remains
   runnable after installation. *)
documentationResourcePath = {
  "ReferencePages",
  "Symbols",
  "EIGENVAL"
};
Scan[
  Function[language,
    resourceSource = FileNameJoin[Join[
      {documentationDirectory, language},
      documentationResourcePath
    ]];
    resourceTarget = FileNameJoin[Join[
      {builtPacletDirectory, "Documentation", language},
      documentationResourcePath
    ]];
    If[!FileExistsQ[resourceSource],
      Print["Missing documentation resource: ", resourceSource];
      Exit[10]
    ];
    CopyFile[resourceSource, resourceTarget, OverwriteTarget -> True]
  ],
  {"English", "ChineseSimplified"}
];

builtEnglishCount = Length@FileNames[
  "*.nb",
  FileNameJoin[{builtPacletDirectory, "Documentation", "English"}],
  Infinity
];
builtChineseCount = Length@FileNames["*.nb", builtChineseTarget, Infinity];
If[
  builtEnglishCount =!= englishPageCount ||
    builtChineseCount =!= chinesePageCount,
  Print[
    "Unexpected built documentation page counts: ",
    {builtEnglishCount, builtChineseCount}
  ];
  Exit[11]
];

Print["[3/4] Rebuilding the manifest for the final bilingual tree..."];
manifestExcludedFiles = {"PacletInfo.wl", "PacletManifest.wl"};
manifestRelativeFiles = Sort@Map[
  FileNameDrop[#, FileNameDepth[builtPacletDirectory]] &,
  Select[
    FileNames["*", builtPacletDirectory, Infinity],
    FileType[#] === File &
  ]
];
manifestRelativeFiles = Complement[
  manifestRelativeFiles,
  manifestExcludedFiles
];
completeManifest = <|
  "Files" -> Map[
    Function[path,
      <|
        "File" -> File[path],
        "Hash" -> FileHash[
          FileNameJoin[{builtPacletDirectory, path}],
          "SHA256",
          "HexString"
        ]
      |>
    ],
    manifestRelativeFiles
  ]
|>;
Put[
  completeManifest,
  FileNameJoin[{builtPacletDirectory, "PacletManifest.wl"}]
];

Print["[4/4] Creating the final bilingual paclet archive..."];
archivePath = buildResult["PacletArchive"];
If[FileExistsQ[archivePath], DeleteFile[archivePath]];
archivePath = CreatePacletArchive[builtPacletDirectory, buildDirectory];
If[archivePath === $Failed || !FileExistsQ[archivePath],
  Print["Unable to create the MagneticKP paclet archive."];
  Exit[12]
];

Print[InputForm@<|
  "Status" -> "Success",
  "PacletArchive" -> archivePath,
  "EnglishDocumentationPages" -> builtEnglishCount,
  "ChineseSimplifiedDocumentationPages" -> builtChineseCount,
  "ManifestFiles" -> Length[manifestRelativeFiles]
|>];
Exit[0];
