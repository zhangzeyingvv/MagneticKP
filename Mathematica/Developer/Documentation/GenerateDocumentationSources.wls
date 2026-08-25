(* Generate the editable English authoring sources for MagneticKP help pages.

   The generated .source.wl files are the readable source boundary. Run
   BuildNotebooks.sh to canonicalize and validate the corresponding .nb files
   under MagneticKP/Documentation/English. Do not edit generated .nb files.
*)

ClearAll["Global`*"];

scriptDirectory = DirectoryName[ExpandFileName[$InputFileName]];
projectRoot = DirectoryName[DirectoryName[scriptDirectory]];
documentationRoot = FileNameJoin[{
  projectRoot,
  "MagneticKP",
  "Documentation",
  "English"
}];
sourceRoot = FileNameJoin[{scriptDirectory, "Sources"}];
expectedRoot = FileNameJoin[{scriptDirectory, "Expected"}];

Scan[
  If[!DirectoryQ[#], CreateDirectory[#, CreateIntermediateDirectories -> True]] &,
  {
    FileNameJoin[{sourceRoot, "Guides"}],
    FileNameJoin[{sourceRoot, "ReferencePages", "Symbols"}],
    FileNameJoin[{sourceRoot, "Tutorials"}],
    FileNameJoin[{expectedRoot, "Guides"}],
    FileNameJoin[{expectedRoot, "ReferencePages", "Symbols"}],
    FileNameJoin[{expectedRoot, "Tutorials"}],
    FileNameJoin[{documentationRoot, "ReferencePages", "Symbols"}]
  }
];

(* These remain public computational symbols, but do not have standalone help
   pages. Remove any previously generated authoring artifacts so they cannot
   survive a rebuild as stale pages. *)
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

(* Keep the real VASP sample beside the vaspEig page. NotebookDirectory[]
   therefore works in both the authoring tree and an installed paclet. *)
eigenvalSource = FileNameJoin[{projectRoot, "Examples", "EIGENVAL"}];
eigenvalTarget = FileNameJoin[{
  documentationRoot,
  "ReferencePages",
  "Symbols",
  "EIGENVAL"
}];
If[!FileExistsQ[eigenvalSource],
  Print["Missing documentation EIGENVAL sample: ", eigenvalSource];
  Exit[2]
];
If[
  !FileExistsQ[eigenvalTarget] ||
    FileHash[eigenvalSource, "SHA256"] =!= FileHash[eigenvalTarget, "SHA256"],
  CopyFile[eigenvalSource, eigenvalTarget, OverwriteTarget -> True]
];

packageEntryPoint = FileNameJoin[{
  projectRoot,
  "MagneticKP",
  "Kernel",
  "init.m"
}];
If[Quiet@Check[Get[packageEntryPoint]; True, False] =!= True,
  Print["Unable to load MagneticKP while generating documentation outputs."];
  Exit[3]
];

writeUTF8TextIfChanged[path_String, text_String] := Module[
  {existing, normalized},
  normalized = StringReplace[
    text,
    RegularExpression["(?m)[ \\t]+$"] -> ""
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

documentationLink[label_String, uri_String] := ButtonBox[
  label,
  BaseStyle -> "Link",
  ButtonData -> "paclet:" <> uri
];

linkListCell[heading_String, style_String, links_List] := Cell[CellGroupData[{
  Cell[heading, style <> "Section"],
  Cell[
    TextData@Riffle[
      documentationLink[#[[1]], #[[2]]] & /@ links,
      " \[FilledSmallSquare] "
    ],
    style
  ]
}, Open]];

historyCell[] := Cell[CellGroupData[{
  Cell["History", "HistorySection"],
  Cell[
    TextData[{
      "New in: ", Cell["1.1.0", "HistoryData", CellTags -> "New"],
      " | Modified in: ", Cell["1.2.0", "HistoryData", CellTags -> "Modified"],
      " | Obsolete in: ", Cell[" ", "HistoryData", CellTags -> "Obsolete"]
    }],
    "History"
  ]
}, Closed]];

categorizationCell[type_String, uri_String] := Cell[CellGroupData[{
  Cell["Categorization", "CategorizationSection"],
  Cell[type, "Categorization", CellLabel -> "Entity Type"],
  Cell["MagneticKP", "Categorization", CellLabel -> "Paclet Name"],
  Cell["MagneticKP`", "Categorization", CellLabel -> "Context"],
  Cell[uri, "Categorization", CellLabel -> "URI"]
}, Closed]];

keywordsCell[keywords_List] := Cell[CellGroupData[
  Prepend[Cell[#, "Keywords"] & /@ keywords, Cell["Keywords", "KeywordsSection"]],
  Closed
]];

usageCell[syntax_String, description_String] := Cell[
  TextData[{
    Cell["   ", "ModInfo"],
    Cell[BoxData[syntax], "InlineFormula"],
    " \[LineSeparator]",
    description
  }],
  "Usage"
];

optionsCell[options_List] := Cell[
  BoxData@GridBox[
    ({
      Cell["      ", "TableRowIcon"],
      #[[1]],
      #[[2]],
      Cell[#[[3]], "TableText"]
    } &) /@ options,
    GridBoxAlignment -> {"Columns" -> {{Left}}}
  ],
  "3ColumnTableMod"
];

outputBoxes[value_] := ToBoxes[value, StandardForm];

(* Store the exact boxes returned by Manipulate behind a compact payload.
   This preserves a working interactive result without tying the notebook to
   the kernel session that generated it. *)
portableDynamicOutput[panel_] := Module[
  {payload = Compress[outputBoxes[panel]]},
  RawBoxes@With[
    {encoded = payload},
    DynamicBox[
      Uncompress[encoded],
      SingleEvaluation -> True,
      TrackedSymbols :> {}
    ]
  ]
];

savedOutputCell[value_] := Module[{boxes},
  boxes = outputBoxes[value];
  Cell[
    BoxData@If[
      FreeQ[boxes, _DynamicModuleBox, Infinity],
      boxes,
      outputBoxes[portableDynamicOutput[value]]
    ],
    "Output"
  ]
];

evaluateExample[spec_Association] := Module[{code, value},
  code = Lookup[spec, "EvaluationCode", spec["ExampleCode"]];
  Print["Evaluating documentation example: ", spec["Name"]];
  value = Block[
    {Print = (Null &)},
    Quiet@Check[ToExpression[code, InputForm], $Failed]
  ];
  If[value === $Failed,
    Print["Documentation example failed: ", spec["Name"]];
    Exit[4]
  ];
  value
];

referencePage[spec_Association] := Module[
  {
    name = spec["Name"],
    usage = spec["Usage"],
    notes = spec["Notes"],
    options = Lookup[spec, "Options", {}],
    related = Lookup[spec, "Related", {}],
    initialization = Lookup[
      spec,
      "Initialization",
      "Needs[\"MagneticKP`\"]"
    ],
    examples,
    exampleCells
  },
  examples = Lookup[
    spec,
    "Examples",
    {KeyTake[spec, {"ExampleText", "ExampleCode", "EvaluationCode"}]}
  ];
  exampleCells = Flatten@MapIndexed[
    Function[{example, index},
      Module[{evaluationSpec, exampleOutput},
        evaluationSpec = Join[
          <|"Name" -> name <> "/" <> ToString[First[index]]|>,
          example
        ];
        exampleOutput = evaluateExample[evaluationSpec];
        {
          Cell[example["ExampleText"], "ExampleText"],
          Cell[BoxData[example["ExampleCode"]], "Input"],
          savedOutputCell[exampleOutput]
        }
      ]
    ],
    examples
  ];
  Notebook[
    Join[
      {
        Cell[CellGroupData[
          Join[
            {Cell[name, "ObjectName"]},
            (usageCell[#[[1]], #[[2]]] & /@ usage),
            (Cell[#, "Notes"] & /@ notes),
            If[options === {}, {}, {optionsCell[options]}]
          ],
          Open
        ]],
        Cell[CellGroupData[{
          Cell["Examples Initialization", "ExamplesInitializationSection"],
          Cell[BoxData[initialization], "ExampleInitialization"]
        }, Open]],
        Cell[CellGroupData[
          Prepend[exampleCells, Cell["Basic Examples", "PrimaryExamplesSection"]],
          Open
        ]]
      },
      If[
        related === {},
        {},
        {linkListCell["Related Functions", "SeeAlso", ({#, "MagneticKP/ref/" <> #} &) /@ related]}
      ],
      {
        linkListCell[
          "Related Guides",
          "MoreAbout",
          {{"MagneticKP", "MagneticKP/guide/MagneticKP"}}
        ],
        linkListCell[
          "Related Tutorials",
          "Tutorials",
          {{"Getting Started with MagneticKP", "MagneticKP/tutorial/GettingStarted"}}
        ],
        historyCell[],
        categorizationCell["Symbol", "MagneticKP/ref/" <> name],
        keywordsCell[spec["Keywords"]]
      }
    ],
    TaggingRules -> <|"Paclet" -> "MagneticKP"|>,
    WindowTitle -> name,
    StyleDefinitions -> FrontEnd`FileName[
      {"Wolfram"},
      "FunctionPageStylesExt.nb",
      CharacterEncoding -> "UTF-8"
    ]
  ]
];

commonKPOptions = {
  {
    "\"GroupInput\"",
    "Automatic",
    "specifies whether the supplied symmetry data contains generators or all elements; Automatic uses metadata in input when available."
  },
  {
    "\"MaxGroupOrder\"",
    "512",
    "sets the largest finite group order accepted while completing generator input."
  },
  {
    "\"ReturnDiagnostics\"",
    "False",
    "specifies whether the returned Association also includes per-order diagnostic information."
  }
};

manualInputCode = StringRiffle[{
  "input = <|",
  "  \"Unitary\" -> <|",
  "    \"C2z\" -> {",
  "      {{0, 1}, {-1, 0}},",
  "      {-kx, -ky, kz}",
  "    },",
  "    \"Inversion\" -> {",
  "      {{1, 0}, {0, -1}},",
  "      {-kx, -ky, -kz}",
  "    }",
  "  |>",
  "|>;"
}, "\n"];

referenceSpecs = {
  <|
    "Name" -> "kpHam",
    "Usage" -> {
      {
        "kpHam[korder, input]",
        "constructs the symmetry-allowed k.p Hamiltonian through korder and returns its data in an Association."
      },
      {
        "kpHam[{n1, n2, ...}, input]",
        "constructs only the listed homogeneous orders."
      }
    },
    "Notes" -> {
      "input is an Association whose \"Unitary\" and optional \"Antiunitary\" entries associate symmetry labels with representation matrices and transformed momentum vectors.",
      "The result contains \"ham\", \"korder\", \"dim\", and \"NumberOfParameters\". Use result[\"ham\"] to obtain the Hamiltonian matrix.",
      "The historical input key \"Anitunitary\" remains accepted, but new input should use \"Antiunitary\"."
    },
    "Options" -> Join[
      {{
        "\"Method\"",
        "\"IterativeSimplification\"",
        "selects \"IterativeSimplification\", \"DirectProductDecomposition\", \"Reynolds\", \"Automatic\", or the compatibility method \"IterativeSimplification2016\"."
      }},
      commonKPOptions
    ],
    "ExampleText" -> "Construct a first-order two-band model from explicitly supplied symmetry data:",
    "ExampleCode" -> manualInputCode <> "\nresult = kpHam[1, input];\nMatrixForm[result[\"ham\"]]",
    "Related" -> {"interfaceRep", "IterativelySimplify", "GETkpHam"},
    "Keywords" -> {"k.p Hamiltonian", "symmetry", "effective model", "band theory"}
  |>,
  <|
    "Name" -> "interfaceRep",
    "Usage" -> {{
      "interfaceRep[group, k, rep]",
      "converts SpaceGroupIrep or MSGCorep representation data into input accepted by kpHam."
    }},
    "Notes" -> {
      "Use an integer group number with SpaceGroupIrep and a BNS pair {spaceGroup, serial} with MSGCorep.",
      "k may be a high-symmetry-point label or exact coordinates. rep may select one representation or a list whose direct sum is required.",
      "For MSGCorep, each integer in rep selects the row with that Index in showMLGCorep; here {2} selects B_2.",
      "The returned Association records the selected k-point label and coordinate system so highSymmetryKPath can label a local path automatically.",
      "Load SpaceGroupIrep or MSGCorep before calling interfaceRep."
    },
    "Options" -> {
      {"\"CartesianCoordinates\"", "True", "specifies whether symmetry operations are converted to Cartesian momentum coordinates."},
      {"\"CalculateGenerators\"", "True", "specifies whether the returned input contains a generating subset instead of every group element."}
    },
    "Initialization" -> "CompoundExpression[Needs[\"MagneticKP`\"], Needs[\"MSGCorep`\"]]",
    "Examples" -> {
      <|
        "ExampleText" -> "List the corepresentations at the B point of magnetic space group 13.65:",
        "ExampleCode" -> "showMLGCorep[{13, 65}, \"B\"]"
      |>,
      <|
        "ExampleText" -> "The second data row is labeled B_2, so select the second corepresentation with {2}. interfaceRep returns the complete input Association:",
        "ExampleCode" -> "input = interfaceRep[{13, 65}, \"B\", {2}]"
      |>,
      <|
        "ExampleText" -> "Pass the returned Association directly to kpHam to construct the first-order Hamiltonian:",
        "ExampleCode" -> "ham = kpHam[1, input][\"ham\"];\nMatrixForm[ham]"
      |>
    },
    "Related" -> {"kpHam", "highSymmetryKPath", "directSum"},
    "Keywords" -> {"MSGCorep", "SpaceGroupIrep", "representation", "magnetic space group"}
  |>,
  <|
    "Name" -> "highSymmetryKPath",
    "Usage" -> {{
      "highSymmetryKPath[input]",
      "returns a local labeled momentum path obtained from the complete little magnetic co-group contained in interfaceRep input."
    }},
    "Notes" -> {
      "input is normally the Association returned by interfaceRep. A manually prepared input with three-dimensional momentum transformations is also accepted.",
      "The default result is an ordinary path list. highSymmetryKPath does not call bandplot; pass the returned path to bandplot or bandManipulate in a separate step.",
      "The path is centered at the local k.p expansion origin {0, 0, 0}. Each retained line has the complete largest magnetic subgroup that fixes a generic point on that line, and symmetry-equivalent lines are represented once by default.",
      "Unitary and antiunitary little-co-group elements remain distinct even when they act on momentum in the same way. In returned diagnostics, SubgroupOrder is the full magnetic subgroup order, while MomentumActionOrder counts distinct momentum-action matrices.",
      "highSymmetryKPath identifies symmetry lines; whether bands stay degenerate along a line depends on the selected representations and Hamiltonian.",
      "Use \"ReturnData\" -> True to inspect the directions, subgroup orders, little co-group order, and the generated path together."
    },
    "Options" -> {
      {"\"Radius\"", "1/10", "sets the extent of each path arm in the coordinate system used by input."},
      {"\"Center\"", "{0, 0, 0}", "sets the local center of the returned path."},
      {"\"CenterLabel\"", "Automatic", "uses the k-point label recorded by interfaceRep, or sets an explicit label string."},
      {"\"EquivalentDirections\"", "\"Representatives\"", "returns one representative of each symmetry-equivalent line class; \"All\" keeps every distinct line."},
      {"\"IncludeNegativeDirections\"", "True", "includes both arms of every line through the local center."},
      {"\"ReturnData\"", "False", "returns only the path, or an Association containing the path and symmetry diagnostics."},
      {"\"MaxGroupOrder\"", "512", "sets the largest little co-group order accepted while completing the supplied transformations."},
      {"\"MaxSubgroups\"", "4096", "sets the largest number of subgroups enumerated while finding fixed lines."}
    },
    "Initialization" -> "CompoundExpression[Needs[\"MagneticKP`\"], Needs[\"MSGCorep`\"]]",
    "Examples" -> {
      <|
        "ExampleText" -> "Generate and display the path for K5 and K6 at K of magnetic space group 191.234:",
        "ExampleCode" -> "input = interfaceRep[\n  {191, 234}, \"K\", {5, 6},\n  \"CartesianCoordinates\" -> False\n];\npath = highSymmetryKPath[input]"
      |>,
      <|
        "ExampleText" -> "Construct the first-order Hamiltonian and pass the displayed path directly to bandplot:",
        "ExampleCode" -> "ham = kpHam[1, input][\"ham\"];\nparams = Complement[Variables[ham], {kx, ky, kz}];\nrules = Thread[params -> N[Range[Length[params]]/10]];\nbandplot[\n  path, 30, ham, rules,\n  \"PlotRange\" -> All, FontSize -> 16, ImageSize -> 500\n]"
      |>
    },
    "Related" -> {"interfaceRep", "kpHam", "bandplot", "bandManipulate"},
    "Keywords" -> {"little co-group", "high-symmetry line", "momentum path", "subgroup"}
  |>,
  <|
    "Name" -> "IterativelySimplify",
    "Usage" -> {{
      "IterativelySimplify[korder, input]",
      "constructs a symmetry-allowed k.p Hamiltonian with the default solver."
    }},
    "Notes" -> {
      "An integer korder includes every order from zero through korder; a list selects only the listed homogeneous orders.",
      "The return format is the same Association used by kpHam."
    },
    "Options" -> commonKPOptions,
    "ExampleText" -> "Construct all terms through first order:",
    "ExampleCode" -> manualInputCode <> "\nIterativelySimplify[1, input][\"ham\"] // MatrixForm",
    "Related" -> {"kpHam", "GETkpHam"},
    "Keywords" -> {"k.p Hamiltonian", "symmetry constraints", "effective model"}
  |>,
  <|
    "Name" -> "GETkpHam",
    "Usage" -> {{
      "GETkpHam[n, input]",
      "constructs only the homogeneous order-n part of a k.p Hamiltonian using the legacy calling convention."
    }},
    "Notes" -> {
      "GETkpHam is retained for notebooks written for earlier MagneticKP releases.",
      "Use kpHam for new work, especially when multiple orders or diagnostic output are required."
    },
    "ExampleText" -> "Construct only the first-order contribution:",
    "ExampleCode" -> manualInputCode <> "\nGETkpHam[1, input][\"ham\"] // MatrixForm",
    "Related" -> {"kpHam", "IterativelySimplify"},
    "Keywords" -> {"legacy interface", "homogeneous order", "k.p Hamiltonian"}
  |>,
  <|
    "Name" -> "bandManipulate",
    "Usage" -> {{
      "bandManipulate[path, n, ham]",
      "creates an interactive band plot with controls for the free parameters in ham."
    }},
    "Notes" -> {
      "Each path segment has the form {{k1, k2}, {label1, label2}} and n gives the sampling count per segment.",
      "The momentum coordinates in path are multiplied by 2 Pi before substitution into ham."
    },
    "Options" -> {
      {"\"PlotRange\"", "All", "sets the displayed energy range."},
      {"\"ManipulateRange\"", "1", "sets a symmetric slider range, or may be given explicitly as {min, max}."}
    },
    "ExampleText" -> "Construct a two-band Hamiltonian with kpHam, then explore it with interactive parameter controls:",
    "ExampleCode" -> manualInputCode <> "\nham = kpHam[1, input][\"ham\"];\npath = {{{{0, 0, 0}, {1/3, 0, 0}}, {\"G\", \"X\"}}};\nbandManipulate[path, 20, ham, \"PlotRange\" -> {-2, 2}]",
    "Related" -> {"bandplot", "bandManipulateWithVASP", "kpHam"},
    "Keywords" -> {"band plot", "interactive", "parameters", "momentum path"}
  |>,
  <|
    "Name" -> "bandplot",
    "Usage" -> {{
      "bandplot[path, n, ham, rules]",
      "draws the numerical eigenvalues of ham along a labeled reciprocal-fractional path after applying parameter rules."
    }},
    "Notes" -> {
      "Each path segment has the form {{k1, k2}, {label1, label2}}. The two momentum points contain three reciprocal-fractional coordinates.",
      "n is the number of subdivisions per segment, so each segment contains n + 1 sampled points. Momentum coordinates are multiplied by 2 Pi before substitution into ham.",
      "rules must supply numerical values for every free parameter in ham.",
      "At a discontinuity between two consecutive segments, the boundary tick combines the two labels as end|start.",
      "The result is a framed Graphics object. Invalid paths, sampling counts, rules, tick specifications, font settings, or nonnumeric Hamiltonians return $Failed with a message."
    },
    "Options" -> {
      {"\"PlotRange\"", "All", "sets the displayed energy range."},
      {"\"YTicks\"", "Automatic", "sets the vertical frame ticks to Automatic, None, numeric positions, or standard tick specifications."},
      {"FontSize", "24", "sets the size of path labels and frame ticks."},
      {"FontFamily", "\"Times\"", "sets the font family of path labels and frame ticks."},
      {"ImageSize", "Automatic", "sets the final Graphics image size."}
    },
    "Examples" -> {
      <|
        "ExampleText" -> "Construct a first-order two-band Hamiltonian with kpHam and assign a numerical value to every model parameter:",
        "ExampleCode" -> manualInputCode <> "\nham = kpHam[1, input][\"ham\"];\nparams = Complement[Variables[ham], {kx, ky, kz}];\nrules = Thread[params -> N[Range[Length[params]]/5]];\npath = {\n  {{{0, 0, 0}, {1/3, 0, 0}}, {\"G\", \"X\"}},\n  {{{1/3, 0, 0}, {1/3, 1/3, 0}}, {\"X\", \"M\"}}\n};\n{Dimensions[ham], rules}"
      |>,
      <|
        "ExampleText" -> "\"PlotRange\" - restrict the displayed energy interval:",
        "ExampleCode" -> "bandplot[path, 30, ham, rules, \"PlotRange\" -> {-2, 2}]"
      |>,
      <|
        "ExampleText" -> "\"YTicks\" - place explicit labeled ticks on the energy axis:",
        "ExampleCode" -> "bandplot[\n  path, 30, ham, rules,\n  \"YTicks\" -> {{-2, \"-2 eV\"}, {0, \"0\"}, {2, \"2 eV\"}}\n]"
      |>,
      <|
        "ExampleText" -> "FontSize - enlarge the path labels and frame ticks:",
        "ExampleCode" -> "bandplot[path, 30, ham, rules, FontSize -> 32]"
      |>,
      <|
        "ExampleText" -> "FontFamily - use Helvetica for the path labels and frame ticks:",
        "ExampleCode" -> "bandplot[path, 30, ham, rules, FontFamily -> \"Helvetica\"]"
      |>,
      <|
        "ExampleText" -> "ImageSize - set the final plot width:",
        "ExampleCode" -> "bandplot[path, 30, ham, rules, ImageSize -> 320]"
      |>
    },
    "Related" -> {"highSymmetryKPath", "bandManipulate", "kpHam"},
    "Keywords" -> {"band plot", "dispersion", "momentum path"}
  |>,
  <|
    "Name" -> "vaspEig",
    "Usage" -> {{
      "vaspEig[file, efermi, spin, startband, endband]",
      "imports selected bands from a VASP EIGENVAL file and shifts their energies by efermi."
    }},
    "Notes" -> {
      "spin selects the energy column to import, and startband through endband select an inclusive band range.",
      "The result is a list of {k, energies} pairs suitable for bandManipulateWithVASP and fittingKP."
    },
    "ExampleText" -> "Import the first two bands from the EIGENVAL sample stored beside this help page:",
    "ExampleCode" -> "file = FileNameJoin[{NotebookDirectory[], \"EIGENVAL\"}];\neigdata = vaspEig[file, 0, 1, 1, 2];\n{Length[eigdata], First[eigdata]}",
    "EvaluationCode" -> "file = " <> ToString[eigenvalSource, InputForm] <> ";\neigdata = vaspEig[file, 0, 1, 1, 2];\n{Length[eigdata], First[eigdata]}",
    "Related" -> {"bandManipulateWithVASP", "fittingKP"},
    "Keywords" -> {"VASP", "EIGENVAL", "band data", "import"}
  |>,
  <|
    "Name" -> "bandManipulateWithVASP",
    "Usage" -> {{
      "bandManipulateWithVASP[ham, eigdata]",
      "interactively compares the eigenvalues of ham with imported VASP eigenvalues."
    }},
    "Notes" -> {
      "eigdata should be the result returned by vaspEig or data in the same {k, energies} form.",
      "The model bands are drawn in black and imported bands in red.",
      "When kpHam is constructed in fractional coordinates about a nonzero expansion point, shift its momentum variables to the absolute fractional coordinates stored in EIGENVAL before comparing the bands."
    },
    "Options" -> {
      {"\"PlotRange\"", "All", "sets the displayed energy range."},
      {"\"ManipulateRange\"", "2 Pi", "sets a symmetric slider range, or may be given explicitly as {min, max}."}
    },
    "Initialization" -> "CompoundExpression[Needs[\"MagneticKP`\"], Needs[\"MSGCorep`\"]]",
    "Examples" -> {
      <|
        "ExampleText" -> "Inspect the corepresentations at K used by the original MagneticKP fitting example:",
        "ExampleCode" -> "showMLGCorep[{156, 50}, {1/3, 1/3, 0}]"
      |>,
      <|
        "ExampleText" -> "Select corepresentations 2 and 3 in fractional coordinates and construct the second-order k.p Hamiltonian:",
        "ExampleCode" -> "input = interfaceRep[{156, 50}, \"K\", {2, 3}, \"CartesianCoordinates\" -> False];\nham = kpHam[2, input][\"ham\"];\nMatrixForm[ham]"
      |>,
      <|
        "ExampleText" -> "Read bands 21 and 22 from the original VASP EIGENVAL file using the Fermi energy from the original example:",
        "ExampleCode" -> "file = FileNameJoin[{NotebookDirectory[], \"EIGENVAL\"}];\neigdata = vaspEig[file, -0.0072, 1, 21, 22];\n{Length[eigdata], First[eigdata]}",
        "EvaluationCode" -> "file = " <> ToString[eigenvalSource, InputForm] <> ";\neigdata = vaspEig[file, -0.0072, 1, 21, 22];\n{Length[eigdata], First[eigdata]}"
      |>,
      <|
        "ExampleText" -> "Shift the local k.p expansion from K to the absolute fractional coordinates in EIGENVAL, then tune the model against the VASP bands:",
        "ExampleCode" -> "hamshift = ham /. {kx -> kx - 1/3, ky -> ky - 1/3};\nbandManipulateWithVASP[hamshift, eigdata, \"PlotRange\" -> {-2, 3}, \"ManipulateRange\" -> {-4, 5}]"
      |>
    },
    "Related" -> {"vaspEig", "fittingKP", "bandManipulate"},
    "Keywords" -> {"VASP", "band comparison", "interactive", "fitting"}
  |>,
  <|
    "Name" -> "fittingKP",
    "Usage" -> {{
      "fittingKP[ham, eigdata, krange, initialRules]",
      "fits the free parameters of ham to selected reference bands."
    }},
    "Notes" -> {
      "eigdata uses the {{kx, ky, kz}, {energy1, energy2, ...}} records returned by vaspEig. krange is a nonempty list of 1-based record indices.",
      "initialRules may be a rule list or Association and must give a real numeric starting value for every free parameter in ham.",
      "fittingKP substitutes the coordinates in eigdata directly into ham; it does not multiply them by 2 Pi. Shift a local k.p Hamiltonian to the absolute fractional coordinates stored in EIGENVAL before fitting.",
      "At each selected k point, sorted model eigenvalue j is compared with reference energy j. Inspect \"ComparisonPlot\" after fitting, especially near band crossings.",
      "Use result[\"FittedParams\"] as the parameter rules for bandplot. \"LSQForInitParams\" and \"LSQForFittedParams\" report the sum of squared band residuals before and after fitting.",
      "The result also contains \"ComparisonPlot\", \"Converged\", \"Iterations\", the effective data selection, and the number of residuals used."
    },
    "Options" -> {
      {"MaxIterations", "100", "sets the maximum number of fitting iterations."},
      {"TimeConstraint", "60", "sets the fitting time limit in seconds, or may be Infinity."},
      {"\"EnergyWindow\"", "All", "uses every selected reference energy, or only energies in the inclusive interval {emin, emax}."},
      {"\"KPointNeighborhood\"", "All", "uses every index in krange, or selects a Radius or componentwise Range around a fractional-coordinate Center."},
      {"\"BandSelection\"", "All", "uses every reference band position, or selects 1-based positions with an integer, increasing list, or positive-step Span."},
      {"\"ResidualWeights\"", "Automatic", "uses equal weights, or a nonnegative record-by-band matrix; a zero entry excludes that residual."},
      {"\"FiniteDifferenceStep\"", "10^-5", "sets the relative numerical-difference step used during fitting."},
      {"\"InitialDamping\"", "10^-3", "sets the initial damping used by the fitter."},
      {"\"FitTolerance\"", "10^-8", "sets the stopping tolerance for the fit."}
    },
    "Initialization" -> "CompoundExpression[Needs[\"MagneticKP`\"], Needs[\"MSGCorep`\"]]",
    "Examples" -> {
      <|
        "ExampleText" -> "Inspect the K-point corepresentations used by the original MagneticKP fitting example:",
        "ExampleCode" -> "showMLGCorep[{156, 50}, {1/3, 1/3, 0}]"
      |>,
      <|
        "ExampleText" -> "Construct the second-order model in fractional coordinates, shift its local K-point momenta to the absolute coordinates in EIGENVAL, and prepare starting values for every generated parameter:",
        "ExampleCode" -> "input = interfaceRep[{156, 50}, \"K\", {2, 3}, \"CartesianCoordinates\" -> False];\nham = kpHam[2, input][\"ham\"];\nhamshift = ham /. {kx -> kx - 1/3, ky -> ky - 1/3};\nparams = Complement[Variables[hamshift], {kx, ky, kz}];\ninitialRules = Thread[params -> {0., 0.19, 4.33, 4.33, 0., 0., 0., 0., 0.}];\nMatrixForm[hamshift]"
      |>,
      <|
        "ExampleText" -> "Read VASP bands 21 and 22 from the supplied original EIGENVAL file:",
        "ExampleCode" -> "file = FileNameJoin[{NotebookDirectory[], \"EIGENVAL\"}];\neigdata = vaspEig[file, -0.0072, 1, 21, 22];\n{Length[eigdata], First[eigdata]}",
        "EvaluationCode" -> "file = " <> ToString[eigenvalSource, InputForm] <> ";\neigdata = vaspEig[file, -0.0072, 1, 21, 22];\n{Length[eigdata], First[eigdata]}"
      |>,
      <|
        "ExampleText" -> "Fit five points near K and inspect the fitted rules, error reduction, and convergence status:",
        "ExampleCode" -> "fit = fittingKP[\n  hamshift, eigdata, Range[50, 70, 5], initialRules\n];\nKeyTake[\n  fit,\n  {\"FittedParams\", \"LSQForInitParams\", \"LSQForFittedParams\", \"Converged\", \"Iterations\"}\n]"
      |>,
      <|
        "ExampleText" -> "Compare the selected VASP points with the initial and fitted k.p bands:",
        "ExampleCode" -> "fit[\"ComparisonPlot\"]"
      |>
    },
    "Related" -> {"vaspEig", "bandManipulateWithVASP", "bandplot"},
    "Keywords" -> {"parameter fitting", "VASP", "least squares", "k.p Hamiltonian"}
  |>,
  <|
    "Name" -> "getGenerator",
    "Usage" -> {{
      "getGenerator[elements, identity, multiply]",
      "returns a generating subset for a finite group described by its elements and multiplication function."
    }},
    "Notes" -> {
      "elements should contain the complete finite group and identity should be its identity element.",
      "The returned subset generates the group, but it is not guaranteed to have the minimum possible size."
    },
    "ExampleText" -> "Find a generator for integers modulo four under addition:",
    "ExampleCode" -> "getGenerator[Range[0, 3], 0, Mod[#1 + #2, 4] &]",
    "Related" -> {"interfaceRep"},
    "Keywords" -> {"finite group", "generator", "group elements"}
  |>,
  <|
    "Name" -> "directSum",
    "Usage" -> {{
      "directSum[a, b]",
      "returns the block-diagonal direct sum of matrices a and b."
    }},
    "Notes" -> {
      "directSum is useful when combining representation matrices for several irreducible representations."
    },
    "ExampleText" -> "Form the direct sum of a two-dimensional and a one-dimensional matrix:",
    "ExampleCode" -> "directSum[{{1, 0}, {0, -1}}, {{I}}] // MatrixForm",
    "Related" -> {"interfaceRep", "kpHam"},
    "Keywords" -> {"direct sum", "block diagonal", "representation matrix"}
  |>
};

(* Evaluate the tutorial fixtures once and save their real results below. *)
tutorialResult = evaluateExample@<|
  "Name" -> "GettingStarted/model",
  "ExampleCode" -> manualInputCode <> "\nkpHam[1, input]"
|>;
tutorialKeys = Keys[tutorialResult];
tutorialHamiltonian = tutorialResult["ham"];
tutorialPath = {
  {{{0, 0, 0}, {1/3, 0, 0}}, {"G", "X"}},
  {{{0, 0, 0}, {0, 1/3, 0}}, {"G", "Y"}},
  {{{0, 0, 0}, {0, 0, 1/3}}, {"G", "Z"}}
};
tutorialPanel = Block[
  {Print = (Null &)},
  bandManipulate[
    tutorialPath,
    40,
    tutorialHamiltonian,
    "PlotRange" -> {-2, 2}
  ]
];
tutorialInterfaceSummary = evaluateExample@<|
  "Name" -> "GettingStarted/interfaceRep",
  "ExampleCode" ->
    "Needs[\"MSGCorep`\"];\ninput = interfaceRep[{13, 65}, \"B\", {2}];\n{input[\"lable\"], input[\"GroupInput\"]}"
|>;

guideEntry[label_String, description_String] := Cell[
  TextData[{
    documentationLink[label, "MagneticKP/ref/" <> label],
    " \[LongDash] " <> description
  }],
  "GuideText"
];

guideSection[title_String, entries_List] := Cell[CellGroupData[
  Prepend[entries, Cell[title, "GuideFunctionsSubsection"]],
  Open
]];

guideNotebook = Notebook[{
  Cell["MagneticKP", "GuideTitle"],
  Cell[
    "MagneticKP constructs symmetry-allowed k.p Hamiltonians, connects to representation packages, plots model bands, and fits model parameters to VASP data.",
    "GuideAbstract"
  ],
  Cell[CellGroupData[{
    Cell["", "GuideFunctionsSection"],
    guideSection["Constructing k.p models", {
      guideEntry["kpHam", "construct a symmetry-allowed k.p Hamiltonian."],
      guideEntry["interfaceRep", "prepare kpHam input from SpaceGroupIrep or MSGCorep."],
      guideEntry["IterativelySimplify", "call the default construction method directly."],
      guideEntry["GETkpHam", "use the legacy homogeneous-order interface."]
    }],
    Cell["\t", "GuideDelimiter"],
    guideSection["Bands and parameter fitting", {
      guideEntry["highSymmetryKPath", "generate a local momentum path from little-co-group subgroups."],
      guideEntry["bandManipulate", "explore a model with interactive parameter controls."],
      guideEntry["bandplot", "plot a numerical model along a labeled momentum path."],
      guideEntry["vaspEig", "import selected bands from a VASP EIGENVAL file."],
      guideEntry["bandManipulateWithVASP", "compare model bands with imported VASP data."],
      guideEntry["fittingKP", "fit model parameters to selected band data."]
    }],
    Cell["\t", "GuideDelimiter"],
    guideSection["Utilities", {
      guideEntry["getGenerator", "find a generating subset of a finite group."],
      guideEntry["directSum", "form a block-diagonal direct sum."]
    }]
  }, Open]],
  Cell[CellGroupData[{
    Cell["Tutorials", "GuideTutorialsSection"],
    Cell[
      TextData[{
        documentationLink[
          "Getting Started with MagneticKP",
          "MagneticKP/tutorial/GettingStarted"
        ]
      }],
      "GuideTutorial"
    ]
  }, Open]],
  historyCell[],
  categorizationCell["Guide", "MagneticKP/guide/MagneticKP"],
  keywordsCell[{"k.p Hamiltonian", "symmetry", "band structure", "VASP"}]
},
  TaggingRules -> <|"Paclet" -> "MagneticKP"|>,
  WindowTitle -> "MagneticKP",
  StyleDefinitions -> FrontEnd`FileName[
    {"Wolfram"},
    "GuidePageStylesExt.nb",
    CharacterEncoding -> "UTF-8"
  ]
];

tutorialNotebook = Notebook[{
  Cell["Getting Started with MagneticKP", "Title"],
  Cell[
    "This tutorial shows the shortest complete workflow: load MagneticKP, enter symmetry data, construct a Hamiltonian, and inspect its bands.",
    "Text"
  ],
  Cell["Load MagneticKP", "Section"],
  Cell["Load the package after installing the MagneticKP paclet:", "Text"],
  Cell[BoxData["Needs[\"MagneticKP`\"]"], "Input"],
  Cell["Enter symmetry data", "Section"],
  Cell[
    "For each symmetry operation, supply its representation matrix and its action on momentum. Put unitary operations under \"Unitary\" and add \"Antiunitary\" only when needed:",
    "Text"
  ],
  Cell[BoxData[manualInputCode], "Input"],
  Cell["Construct the Hamiltonian", "Section"],
  Cell[
    "Use an integer order to include every order from zero through that value. The returned Association keeps the Hamiltonian and basic summary data together:",
    "Text"
  ],
  Cell[BoxData["result = kpHam[1, input];\nKeys[result]"], "Input"],
  savedOutputCell[tutorialKeys],
  Cell["Extract the matrix from the \"ham\" entry:", "Text"],
  Cell[BoxData["ham = result[\"ham\"];\nMatrixForm[ham]"], "Input"],
  savedOutputCell[MatrixForm[tutorialHamiltonian]],
  Cell["Inspect model bands", "Section"],
  Cell[
    "Define a labeled momentum path and use bandManipulate while choosing parameter values:",
    "Text"
  ],
  Cell[BoxData[StringRiffle[{
    "path = {",
    "  {{{0, 0, 0}, {1/3, 0, 0}}, {\"G\", \"X\"}},",
    "  {{{0, 0, 0}, {0, 1/3, 0}}, {\"G\", \"Y\"}},",
    "  {{{0, 0, 0}, {0, 0, 1/3}}, {\"G\", \"Z\"}}",
    "};",
    "bandManipulate[path, 40, ham, \"PlotRange\" -> {-2, 2}]"
  }, "\n"]], "Input"],
  savedOutputCell[tutorialPanel],
  Cell["Use representation databases", "Section"],
  Cell[
    "When representation matrices are available from MSGCorep or SpaceGroupIrep, interfaceRep prepares the same input format automatically:",
    "Text"
  ],
  Cell[BoxData["Needs[\"MSGCorep`\"];\ninput = interfaceRep[{13, 65}, \"B\", {2}];\n{input[\"lable\"], input[\"GroupInput\"]}"], "Input"],
  savedOutputCell[tutorialInterfaceSummary],
  Cell["Next steps", "Section"],
  Cell[
    TextData[{
      "See ", documentationLink["kpHam", "MagneticKP/ref/kpHam"],
      " for all construction options, ", documentationLink["bandplot", "MagneticKP/ref/bandplot"],
      " for a fixed publication-style band plot, and ", documentationLink["vaspEig", "MagneticKP/ref/vaspEig"],
      " with ", documentationLink["fittingKP", "MagneticKP/ref/fittingKP"],
      " for fitting to VASP data."
    }],
    "Text"
  ],
  historyCell[],
  categorizationCell["Tech Note", "MagneticKP/tutorial/GettingStarted"],
  keywordsCell[{"getting started", "k.p Hamiltonian", "symmetry input", "band plot"}]
},
  TaggingRules -> <|"Paclet" -> "MagneticKP"|>,
  WindowTitle -> "Getting Started with MagneticKP",
  StyleDefinitions -> FrontEnd`FileName[
    {"Wolfram"},
    "TechNotePageStylesExt.nb",
    CharacterEncoding -> "UTF-8"
  ]
];

writePage[relativePath_String, notebook_Notebook, expected_List] := Module[
  {sourcePath, expectedPath},
  sourcePath = FileNameJoin[{sourceRoot, relativePath <> ".source.wl"}];
  expectedPath = FileNameJoin[{expectedRoot, relativePath <> ".expected.txt"}];
  writeUTF8TextIfChanged[
    sourcePath,
    ToString[notebook, InputForm, PageWidth -> 120]
  ];
  writeUTF8TextIfChanged[
    expectedPath,
    StringRiffle[DeleteDuplicates[expected], "\n"] <> "\n"
  ];
  Print[sourcePath];
];

writePage[
  FileNameJoin[{"Guides", "MagneticKP"}],
  guideNotebook,
  {
    "MagneticKP constructs symmetry-allowed k.p Hamiltonians",
    "Constructing k.p models",
    "Bands and parameter fitting",
    "Getting Started with MagneticKP"
  }
];

Scan[
  Function[spec,
    writePage[
      FileNameJoin[{"ReferencePages", "Symbols", spec["Name"]}],
      referencePage[spec],
      Join[
        {
          spec["Name"],
          spec["Usage"][[1, 2]],
          spec["Name"] <> "["
        },
          If[
            spec["Name"] === "interfaceRep",
          {
            "showMLGCorep[{13, 65}",
            "input = interfaceRep[{13, 65}",
            "ham = kpHam[1, input]"
          },
          If[
            spec["Name"] === "highSymmetryKPath",
            {
              "input = interfaceRep[",
              "{191, 234}",
              "path = highSymmetryKPath[input]",
              "ham = kpHam[1, input]",
              "bandplot["
            },
            If[
              spec["Name"] === "bandManipulate",
              {"ham = kpHam[1, input]"},
              If[
                spec["Name"] === "bandplot",
                {
                  "ham = kpHam[1, input]",
                  "rules = Thread[params",
                  "bandplot[path, 30, ham, rules",
                  "YTicks"
                },
                If[
                  spec["Name"] === "bandManipulateWithVASP",
                  {
                    "showMLGCorep[{156, 50}",
                    "input = interfaceRep[{156, 50}",
                    "ham = kpHam[2, input]",
                    "eigdata = vaspEig[file, -0.0072, 1, 21, 22]",
                    "hamshift = ham /.",
                    "bandManipulateWithVASP[hamshift, eigdata"
                  },
                  {}
                ]
              ]
            ]
          ]
        ]
      ]
    ]
  ],
  referenceSpecs
];

writePage[
  FileNameJoin[{"Tutorials", "GettingStarted"}],
  tutorialNotebook,
  {
    "Getting Started with MagneticKP",
    "This tutorial shows the shortest complete workflow",
    "result = kpHam[1, input]",
    "bandManipulate[path, 40, ham"
  }
];

Print[
  "Generated English MagneticKP documentation sources: ",
  Length[referenceSpecs] + 2,
  " pages."
];
Exit[0];
