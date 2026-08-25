.. MagneticKP documentation master file, created by
   sphinx-quickstart on Fri Mar 24 17:57:44 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to MagneticKP's documentation!
======================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

Introduction
============
MagneticKP is an open source package for constructing symmetry-constrained
:math:`k \cdot p` Hamiltonians from (co)representation matrices.  Version
1.2.0 constructs the momentum and adjoint representations dynamically, no
longer loads the precomputed ``fd.pkl`` table, and includes the MSGCorep
catalog workflow used by the Web interface.


Installation
============

Before using MagneticKP, first unzip the "MagneticKP-main.zip" file and then go into the ``Python-Web``
directory and run

.. code-block:: console

   python setup.py install

That's all for installation.

If you do not want to change the SymPy version in the base environment, create
a conda environment:

.. code-block:: console

   conda create -n magkp python=3
   conda activate magkp
   python setup.py install

After you finished using MagneticKP, exit the environment by:

.. code-block:: console

   conda deactivate magkp

and run 

.. code-block:: console

   conda activate magkp

to activate MagneticKP.

Web interface
=============

Install MagneticKP and start the bilingual local interface:

.. code-block:: console

   python -m pip install -e .
   python -m magnetickp.web

Then open ``http://127.0.0.1:8000``.  The interface follows the
HilbertBalance visual language and provides matrix editors, verified examples,
``CInfinity`` input, LaTeX/JSON export, solver diagnostics and an interactive
band plot.  Band eigenvalues and parameter sliders are evaluated locally in
the browser, so the plotting workflow does not add a separate
``bandManipulate`` function to the Python API.  High-symmetry paths follow the
segment-list semantics in ``Mathematica/Examples/GeneralExamples.nb``.  For example,
``Γ:-2/3,1/3,0 -> K:0,0,0; K:0,0,0 -> M:0,1/2,0`` represents two segments.
The browser accepts exact fractional coordinates and keeps disconnected
segments visually disconnected.

The default catalog mode uses the installed MSGCorep package and lets users
select a space group, its BNS magnetic group, a
high-symmetry point or line, and one or more little-group corepresentations.
Both Cartesian momentum coordinates (the Mathematica default) and reciprocal
fractional coordinates are supported.  MSGCorep remains a separate dependency
behind a narrow adapter, so the crystallographic database and MagneticKP can
be maintained independently.  Manual matrix input remains available.

The adapter is regression-tested against the Mathematica ``interfaceRep`` and
``kpHam`` pipeline.  The comparison
checks exact coefficient-subspace equality rather than requiring the two
implementations to choose the same free-parameter basis.

Web expressions use a strict allowlist grammar and symbolic calculations run
in disposable worker processes with a hard timeout.  Set
``MAGNETICKP_WEB_TIMEOUT`` to change the default 30-second limit and
``MAGNETICKP_BASE_PATH`` when serving below a reverse-proxy subpath.


Capability of MagneticKP
========================
.. automodule:: magnetickp
    :members:

Solver backends
===============

The default call is backward compatible and uses the iterative generator
solver:

.. code-block:: python

    result = model.getkpHam(2)

The symbolic coefficients are grouped by momentum order and numbered from
one: ``e_i`` for order 0, ``t_i`` for order 1, ``r_i`` for order 2, and
``s_i`` for order 3. Orders 4 and above retain the previous ``C`` naming
convention.

The solver can be selected explicitly:

.. code-block:: python

    result = model.getkpHam(
        2,
        method="reynolds",
        group_input="generators",
        max_group_order=512,
        return_diagnostics=True,
    )

The accepted backends are ``"iterative"``, ``"direct_product"``,
``"reynolds"`` and ``"automatic"``.  The Mathematica-style spelling
``"DirectProductDecomposition"`` is also accepted.

For Reynolds projection, set ``group_input="generators"`` when the input
contains generators.  MagneticKP closes the exact finite group in paired
``(F, J, eta)`` form before averaging.  Set ``group_input="elements"`` only
when the input already contains all group elements.  The average is verified;
if projective representatives induce a finite central extension, MagneticKP
completes that paired-factor closure and sets
``element_closure_completed=True``.  If exact closure cannot be established
within ``max_group_order``, or if the completed average still fails the
projector checks, the implementation falls back to the iterative solver and
reports the reason in ``Diagnostics``.

``method="automatic"`` currently uses the conservative iterative backend.
Reynolds performance is problem dependent: constructing and row-reducing a
group average may cost more than imposing a small generator set iteratively.
Select ``method="reynolds"`` explicitly when a benchmark for the target
representation justifies it.

On the release test machine, the new cold-cache ``F`` constructor was about
4.8x, 26.4x and 289.8x faster than the legacy linear solve at homogeneous
orders 2, 4 and 8.  The batched ``J`` constructor was about 1.9x, 13.1x and
41.3x faster at state dimensions 2, 4 and 8.  In contrast, explicit Reynolds
projection was slower than generator iteration in the tested 48-element
``O_h`` cases, which motivates the conservative automatic policy.

Input and output conventions
============================

Each symmetry item is a pair ``(D, R)``.  ``D`` is the representation or
corepresentation matrix.  ``R`` may be a momentum matrix or the corresponding
list of linear momentum images.  The corrected key ``"Antiunitary"`` and the
legacy key ``"Anitunitary"`` are both accepted.

One internal continuous spin rotation may instead be supplied as the
three-item record ``(D_theta, (kx, ky, kz), theta)`` under the key
``"CInfinity"`` (the alias ``"Cinfty"`` is also accepted):

.. code-block:: python

    from sympy import I, Matrix, exp, symbols
    from magnetickp.magkp import kx, ky, kz

    theta = symbols("theta")
    D_theta = Matrix.diag(
        exp(-I * theta / 2),
        exp(I * theta / 2),
    )
    model = magnetickp.kpHam({
        "Unitary": {
            "CInfinity": (D_theta, (kx, ky, kz), theta),
            "SpinFlip": (Matrix([[0, 1], [1, 0]]), (kx, ky, kz)),
        }
    })
    result = model.getkpHam(1, return_diagnostics=True)

The third item declares the only continuous group parameter.  MagneticKP
extracts the Hermitian infinitesimal generator, restricts the operator basis
to its commutant, and then applies the ordinary two-item discrete records.
This first implementation supports exactly one continuous unitary entry,
requires trivial momentum action and no additional free symbols in
``D_theta``, and verifies that every discrete operation maps the generator to
``+Q`` or ``-Q``.  Reynolds projection, when selected, is applied only to the
remaining finite discrete factor group.

For an ``n``-dimensional representation, the default operator basis contains
exactly :math:`n^2` Hermitian matrices.  The historical order is retained for
``n = 1, 2, 3, 4, 6, 8``; all other positive dimensions use diagonal matrix
units followed by symmetric and imaginary-antisymmetric off-diagonal pairs.

Examples
========
Here we provide two examples for constructing the :math:`k \cdot p` model.
First is timing and constucting the :math:`k \cdot p` model for 226.123 magnetic space group L4L4 corep

.. code-block:: python

    import time
    from sympy.matrices import  eye,zeros
    from sympy import sqrt, I
    from sympy.matrices import Matrix
    from sympy.parsing.mathematica import parse_mathematica
    import magnetickp
    # Matrix representation for symmetry operators
    C31=(Matrix(parse_mathematica("{{-1, 0, 0, 0}, {0, -1, 0, 0}, {0, 0, -1, 0}, {0, 0, 0, -1}}")))
    C2b=(Matrix(parse_mathematica("{{I, 0, 0, 0}, {0, -I, 0, 0}, {0, 0, I, 0}, {0, 0, 0, -I}}")))
    Inv=(Matrix(parse_mathematica("{{0, 1, 0, 0}, {1, 0, 0, 0}, {0, 0, 0, 1}, {0, 0, 1, 0}}")))
    T=(Matrix(parse_mathematica("{{0, 0, 0, -1}, {0, 0, -1, 0}, {0, 1, 0, 0}, {1, 0, 0, 0}}")))

    #  initialization
    a=magnetickp.kpHam({'Unitary':
                {'C31':(C31,Matrix([[0,0,1],[1,0,0],[0,1,0]])),
                 'C2b':(C2b,Matrix([[0,-1,0],[-1,0,0],[0,0,-1]])),
                 'Inv':(Inv,Matrix([[-1,0,0],[0,-1,0],[0,0,-1]])),
                 },
                'Anitunitary':
                {'T':(T,Matrix([[-1,0,0],[0,-1,0],[0,0,-1]]))
                    }
                })
    start_time = time.perf_counter()

    # Construct the 2nd order model for 226.123 magnetic space group L4L4 corep
    kp=(a.getkpHam(2))
    
    end_time = time.perf_counter()
    runtime = end_time - start_time
    
    print(f"The code took {runtime:.6f} seconds to execute.")
    
    print(kp)
    
    
Second is timing and constucting the :math:`k \cdot p` model for 218.82 magnetic space group R5R5 corep

.. code-block:: python
 
    import time
    from sympy.matrices import  eye,zeros
    from sympy import sqrt, I
    from sympy.matrices import Matrix
    from sympy.parsing.mathematica import parse_mathematica
    import magnetickp
    
    C2x=(Matrix(parse_mathematica("{{-1, 0, 0, 0, 0, 0}, {0, -1, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0}, {0, 0, 0, -1, 0, 0}, {0, 0, 0, 0, -1, 0}, {0, 0, 0, 0, 0, 1}}")))
    C2y=(Matrix(parse_mathematica("{{-1, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 0, 0}, {0, 0, -1, 0, 0, 0}, {0, 0, 0, -1, 0, 0}, {0, 0, 0, 0, 1, 0}, {0, 0, 0, 0, 0, -1}}")))
    C31=(Matrix(parse_mathematica("{{0, 1, 0, 0, 0, 0}, {0, 0, -1, 0, 0, 0}, {-1, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 1, 0}, {0, 0, 0, 0, 0, -1}, {0, 0, 0, -1, 0, 0}}")))
    S4x=(Matrix(parse_mathematica("{{0,-I,0,0,0,0},{I,0,0,0,0,0},{0,0,-I,0,0,0},{0,0,0,0,I,0},{0,0,0,-I,0,0},{0,0,0,0,0,I}}")))
    T=(Matrix(parse_mathematica("{{0, 0, 0, 1, 0, 0}, {0, 0, 0, 0, 1, 0}, {0, 0, 0, 0, 0, 1}, {1, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0}}")))
    
    a=magnetickp.kpHam({'Unitary':
                {'C2x':(C2x,Matrix([[1,0,0],[0,-1,0],[0,0,-1]])),
                 'C2y':(C2y,Matrix([[-1,0,0],[0,1,0],[0,0,-1]])),
                 'C31':(C31,Matrix([[0,0,1],[1,0,0],[0,1,0]])),
                 'S4x':(S4x,Matrix([[-1,0,0],[0,0,1],[0,-1,0]]))
                 },
                'Anitunitary':
                {'T':(T,Matrix([[-1,0,0],[0,-1,0],[0,0,-1]]))
                    }
                })
    start_time = time.perf_counter()
    kp=(a.getkpHam(2))
    
    end_time = time.perf_counter()
    runtime = end_time - start_time
    
    print(f"The code took {runtime:.6f} seconds to execute.")
    print(kp)


Feedback
========

Please send comments or suggestions for improvement to <zhangzeyingvv@gmail.com> or open an issue at <https://github.com/zhangzeyingvv/MagneticKP>.
