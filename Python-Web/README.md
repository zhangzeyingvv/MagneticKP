# MagneticKP for Python/Web

[简体中文](README.zh-CN.md) | English

This directory contains the stable Python/Web distribution:

* `magnetickp/`: Python package and bilingual Web interface;
* `tests/`: Python and browser-side tests;
* `setup.py` and `MANIFEST.in`: packaging configuration.

MagneticKP and the bundled Python MSGCorep distribution are pure Python
packages and include their browser assets and databases. Their dependencies
include platform-specific NumPy, SciPy, and spglib wheels, but the installed
application does not require Mathematica.

The `dist/` directory contains the matching MagneticKP and MSGCorep wheels.
Create and activate a virtual environment before installation so that the
package and its dependencies do not modify the system Python environment:

```console
cd /path/to/MagneticKP1.2/MagneticKP/Python-Web
python3 -m venv env
source env/bin/activate
python -m pip install --upgrade pip
python -m pip install --find-links dist \
  dist/magnetickp-1.2.0-py3-none-any.whl
magnetickp-web
```

The installation command installs the Python core, MSGCorep catalog, Web
application, static plotting, and VASP fitting together. Run `deactivate`
when you want to leave the virtual environment. To use MagneticKP again,
return to this directory and run `source env/bin/activate` first.

The application listens on `http://127.0.0.1:8000` by default. Set
`MAGNETICKP_WEB_HOST` or `MAGNETICKP_WEB_PORT` before starting it to override
the address.

After constructing `ham` and `path`, plotting is one call:

```python
from magnetickp import bandplot
bandplot(ham, path, parameter_values)
```

To build redistributable source and wheel packages:

```console
python3 -m pip install build
python3 -m build
```

Place `msgcorep-0.1.0a0-py3-none-any.whl` beside the generated MagneticKP
wheel in `dist/` before distributing the complete local installation bundle.

For development, install the editable package and test dependencies from the
same local wheel directory with:

```console
python3 -m pip install --find-links dist -e '.[test]'
python3 -m pytest -q
```
