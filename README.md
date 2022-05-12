# MagneticKP

A package using iterative simplification algorithm for quickly constructing k·p models.


## Installation

The steps of installing MagneticKP is exactly the same as installing MagneticTB:
Unzip the "MagneticKP-main.zip" file and copy the MagneticKP directory to any directory in $Path. e.g.,
copy to ```FileNameJoin[{$UserBaseDirectory, "Applications"}]```.


Then one can use the package after running ```Needs["MagneticKP`"]```.
The version of Mathematica should higher or equal to 11.3.

## Capabilities of MagneticKP

* Construct the k·p model for given (co)representation matrices
* Both iterative simplification algorithm and direct-product decomposition algorithm are implemented in MagneticKP
* By Interfaceing with [```SpaceGroupIrep```](https://github.com/goodluck1982/SpaceGroupIrep) or ```MSGCorep``` packages, it can directly output the k · p Hamiltonian around arbitrary momentum and expanded to arbitrary order in k.


See [Comput. Phys. Commun. 270, 108153 (2022)](https://www.sciencedirect.com/science/article/abs/pii/S0010465521002654) [(arXiv:2012.08871)](https://arxiv.org/abs/2105.09504) for detail (please cite this paper if you use our code for your research).
