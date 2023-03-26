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
MagneticKP is an open source package for constructing :math:`k \cdot p` Hamiltonian for arbitrary magnetic space group.


Installition
============

Before using MagneticKP, first unzip the "MagneticKP-main.zip" file and then go into the MagneticKP-python
directory and run

.. code-block:: console

   python setup.py install

That's all for Installition.

If you don't want to change the version of sympy in the base environment, you can creates a conda environment:

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

to active MagneticKP.


Capability of MagneticKP
========================
.. automodule:: magkp
    :members:

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
