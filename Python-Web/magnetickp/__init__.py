from .magkp import kpHam
from .kpath import HighSymmetryPathError, format_band_path, high_symmetry_k_path
from .plotting import BandPlotError, bandplot
from .vasp import VaspEigenvalError, read_vasp_eigenval
from .fitting import KpFittingError, fit_kp_to_vasp, fittingKP

__author__ = 'Zhang Zeying'
__email__ = 'zhangzeyingvv@gmail.com'
__version__ = '1.2.0'

__all__ = [
    "BandPlotError",
    "HighSymmetryPathError",
    "KpFittingError",
    "VaspEigenvalError",
    "bandplot",
    "format_band_path",
    "fit_kp_to_vasp",
    "fittingKP",
    "high_symmetry_k_path",
    "kpHam",
    "read_vasp_eigenval",
]


if __name__ == '__main__':
    print('test')
