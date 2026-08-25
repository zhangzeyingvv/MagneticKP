# MagneticKP for Wolfram Language

[English](README.md) | [简体中文](README.zh.md)

This directory contains the stable Mathematica distribution of MagneticKP.
The ready-to-install package is
[`build/MagneticKP-1.2.0.paclet`](build/MagneticKP-1.2.0.paclet).

## Install the paclet

Open Mathematica, evaluate the following code, and select the downloaded
`MagneticKP-1.2.0.paclet` file:

```wl
pacletFile = SystemDialogInput["FileOpen"];
PacletInstall[pacletFile, ForceVersionInstall -> True]
```

`ForceVersionInstall -> True` replaces an existing installation of the same
version, which is useful when reinstalling the release package. Start a new
kernel after installation, then load MagneticKP:

```wl
Needs["MagneticKP`"]
```

Confirm the installed paclet and open its main help page with:

```wl
PacletFind["MagneticKP"]
Documentation`HelpLookup["paclet:MagneticKP/guide/MagneticKP"]
```

The installed help contains English and simplified-Chinese pages. Function
pages and tutorials include saved outputs, so users can see the expected result
before evaluating an example themselves. No manual copy into
`$UserBaseDirectory/Applications` is required.

## Development and testing

The remaining directories contain the source, examples, documentation, and validation tools:

* `MagneticKP/`: paclet source;
* `Tests/`: unit, regression, integration, and exhaustive tests;
* `Benchmarks/`: benchmark and prototype scripts;
* `Examples/`: notebooks and sample data;
* `Documentation/`: PDF manual sources and rendered manual;
* `Developer/Documentation/`: generated-help specifications, translations,
  expected-text fixtures, and notebook build scripts.

## Generate the help notebooks

The files under `MagneticKP/Documentation/` are generated authoring notebooks.
Do not edit them directly. Regenerate the English specification, the simplified
Chinese translation, and the canonical notebooks from the maintenance workspace
root with:

```console
env OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /Applications/Mathematica.app/Contents/MacOS/WolframKernel \
  -noinit -script MagneticKP/Mathematica/Developer/Documentation/GenerateDocumentationSources.wls
env OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /Applications/Mathematica.app/Contents/MacOS/WolframKernel \
  -noinit -script MagneticKP/Mathematica/Developer/Documentation/GenerateChineseDocumentationSources.wls
MagneticKP/Mathematica/Developer/Documentation/BuildNotebooks.sh
```

Edit `GenerateDocumentationSources.wls` for English content and page structure,
and `ChineseTranslations.wl` for simplified-Chinese wording. `BuildNotebooks.sh`
canonicalizes and validates every generated notebook.

## Build the paclet

Build the distributable paclet from the maintenance workspace root with:

```console
wolframscript -file MagneticKP/Mathematica/BuildPaclet.wls
```

The archive is written to `MagneticKP/Mathematica/build/`.

After building, verify bilingual page coverage, help links, built notebook
parsing, paclet registration, and the core documented examples with:

```console
env OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /Applications/Mathematica.app/Contents/MacOS/WolframKernel \
  -noinit -script MagneticKP/Mathematica/Developer/Documentation/VerifyDocumentation.wls
```

From the maintenance workspace root, run the regular test suite with:

```console
wolframscript -file MagneticKP/Mathematica/Tests/RunTests.wls
```
