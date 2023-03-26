# MagneticKP

A package use iterative simplification algorithm for quickly constructing k·p models.


## Installation

The steps of installing MagneticKP is exactly the same as installing [MagneticTB](https://github.com/zhangzeyingvv/MagneticTB):
Unzip the "MagneticKP-main.zip" file and copy the MagneticKP directory to any directory in $Path. e.g.,
copy to ```FileNameJoin[{$UserBaseDirectory, "Applications"}]```.


Then one can use the package after running ```Needs["MagneticKP`"]```.
The version of Mathematica should higher or equal to 11.3.

## Capabilities of MagneticKP

* Construct the k·p model for given (co)representation matrices.
* Both iterative simplification algorithm and direct-product decomposition algorithm are implemented in MagneticKP.
* By Interfaceing with [SpaceGroupIrep](https://github.com/goodluck1982/SpaceGroupIrep) or [MSGCorep](https://github.com/goodluck1982/MSGCorep) packages, it can directly output the k·p Hamiltonian around arbitrary momentum, expanded to arbitrary order in k.

See [arXiv:2205.05830](https://arxiv.org/abs/2205.05830) for detail (please cite this preprint if you use our code for your research).

## Examples

See examples.nb.

## Release Notes

v1.01   2022/06/30

* Support low-dimensional (one- and two-dimensional ) kp Hamiltonian ([Phys. Rev. B 107, 075405 (2023)](https://link.aps.org/doi/10.1103/PhysRevB.107.075405) ([arXiv:2210.11080](https://arxiv.org/abs/2210.11080))). See example.nb file for examples.

v1.02   2022/10/15

* Add ```bandManipulate``` and ```bandplot``` functions for plot the band structure of kp Hamiltonian. See example.nb file for example.
* Support [MSGCorep](https://github.com/goodluck1982/MSGCorep) v1.0.0

v1.03   2023/03/26

* Add a python version of MagneticKP
* Please refer to MagneticKP-python\doc\html\index.html for detailed information on the installation, capabilities, and examples of the Python version of MagneticKP.


