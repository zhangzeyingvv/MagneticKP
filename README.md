# MagneticKP

[简体中文](README.zh-CN.md)

MagneticKP constructs symmetry-constrained k·p Hamiltonians from symmetry operations and (co)representation matrices. It provides Mathematica, Python, and Web interfaces.

## Installation

### Mathematica

Download `MagneticKP-1.2.0.paclet` from the project [Releases page](https://github.com/zhangzeyingvv/MagneticKP/releases), then run in Mathematica:

```wl
PacletInstall[SystemDialogInput["FileOpen"]]
Needs["MagneticKP`"]
```

Select the downloaded `.paclet` file in the file dialog. Mathematica 12.1 or later is required.

### Python/Web

Enter the `Python-Web` directory. It is recommended to create and use a virtual environment in this directory:

```bash
python3 -m venv env
source env/bin/activate
python -m pip install --upgrade pip
python -m pip install --find-links dist dist/magnetickp-1.2.0-py3-none-any.whl
magnetickp-web
```

The installation command automatically installs `msgcorep` from the same `dist/` directory; no separate installation is required. Then open `http://127.0.0.1:8000` in a browser.

When finished, leave the virtual environment with:

```bash
deactivate
```

To use MagneticKP again later, enter the `Python-Web` directory and run `source env/bin/activate` before starting Python or `magnetickp-web`.

## MagneticKP capabilities

* Construct k·p Hamiltonians to arbitrary order;
* Support magnetic, non-magnetic, and spin-orbit-coupled systems;
* Support the internal continuous spin rotation `CInfinity`;
* Select magnetic space groups, k points, and corepresentations through MSGCorep;
* Automatically calculate local high-symmetry k paths from the little co-group at the selected k point;
* Plot bands and fit parameters to first-principles band structures;
* Provide Mathematica, Python, and Web interfaces.

## How to open the manuals

After installing the Paclet, open the built-in Mathematica help with:

```wl
Documentation`HelpLookup["paclet:MagneticKP/guide/MagneticKP"]
```

The built-in help includes English and simplified-Chinese pages.

For Python/Web, run `magnetickp-web`, then open `http://127.0.0.1:8000/help/` in a browser.

For more details, see [Mathematica/README.md](Mathematica/README.md) and [Python-Web/README.md](Python-Web/README.md).

## Simple examples

### Mathematica

```wl
input = <|
  "Unitary" -> <|
    "C2z" -> {PauliMatrix[3], {-kx, -ky, kz}}
  |>
|>;

kpHam[2, input]
```

### Python

```python
from sympy import Matrix
from magnetickp import kpHam

input_data = {
    "Unitary": {
        "C2z": (
            Matrix([[1, 0], [0, -1]]),
            Matrix([[-1, 0, 0], [0, -1, 0], [0, 0, 1]]),
        )
    }
}

result = kpHam(input_data).getkpHam(2)
print(result["ham"])
```

## Release notes

### v1.2.0 — 2026/8/20

* Added the Python/Web interface for selecting magnetic space groups and corepresentations;
* Added automatic local high-symmetry k-path calculation from the little co-group at the selected k point;
* Added extensive help documentation for Mathematica and Python/Web;

### v1.1.0 — 2025/10/31

* Accelerated the momentum transformation matrix calculation;
* Added the Reynolds solver;
* Added support for the continuous spin rotation `CInfinity`;

### v1.04 — 2023/10/31

* Added `vaspEig`, `bandManipulateWithVASP`, and `fittingKP` for fitting k·p Hamiltonians.

### v1.03 — 2023/03/26

* Added the Python version.

### v1.02 — 2022/10/15

* Added `bandManipulate` and `bandplot` for plotting band structures;
* Added support for MSGCorep v1.0.0.

### v1.01 — 2022/06/30

* Added support for one- and two-dimensional k·p Hamiltonians.
