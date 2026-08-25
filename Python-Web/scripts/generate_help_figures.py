"""Generate the computed figures embedded in the Python/Web help page."""

from pathlib import Path
import os
import sys
import tempfile


os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(tempfile.gettempdir()) / "magnetickp-matplotlib"),
)

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from sympy import Matrix, Rational


PROJECT_ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = PROJECT_ROOT.parent
OUTPUT_DIRECTORY = PROJECT_ROOT / "magnetickp/web/static/help/assets"
EIGENVAL_PATH = REPOSITORY_ROOT / "Mathematica-Testing/Examples/EIGENVAL"

sys.path.insert(0, str(PROJECT_ROOT))

from magnetickp import bandplot, fittingKP, kpHam, read_vasp_eigenval  # noqa: E402
from magnetickp.catalog import build_catalog_core_input  # noqa: E402
from magnetickp.transforms import kx, ky, kz  # noqa: E402


plt.rcParams.update(
    {
        "axes.edgecolor": "#9aa6bc",
        "axes.labelcolor": "#26334f",
        "axes.titlecolor": "#1f2a44",
        "figure.facecolor": "white",
        "font.family": "DejaVu Sans",
        "savefig.facecolor": "white",
        "svg.fonttype": "path",
        "svg.hashsalt": "magnetickp-help",
        "xtick.color": "#46536d",
        "ytick.color": "#46536d",
    }
)


def save_svg(figure, filename):
    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    figure.savefig(
        OUTPUT_DIRECTORY / filename,
        format="svg",
        bbox_inches="tight",
        metadata={"Date": None},
    )
    plt.close(figure)


def bandplot_example():
    inputdict = {
        "Unitary": {
            "C2z": (
                Matrix([[0, 1], [-1, 0]]),
                (-kx, -ky, kz),
            ),
            "Inversion": (
                Matrix([[1, 0], [0, -1]]),
                (-kx, -ky, -kz),
            ),
        }
    }
    ham = kpHam(inputdict).getkpHam(1)["ham"]
    parameters = sorted(ham.free_symbols - {kx, ky, kz}, key=str)
    rules = dict(zip(parameters, (0.2, 0.4, 0.6, 0.8)))
    path = [
        ((0, 0, 0), (1 / 3, 0, 0), ("G", "X")),
        ((1 / 3, 0, 0), (1 / 3, 1 / 3, 0), ("X", "M")),
    ]
    return ham, path, rules


def generate_bandplot_figure():
    ham, path, rules = bandplot_example()
    figure, axes = plt.subplots(2, 3, figsize=(12, 7), constrained_layout=True)
    cases = [
        ("PlotRange", {"ylim": (-2, 2)}),
        ("YTicks", {"yticks": ([-2, 0, 2], ["-2 eV", "0", "2 eV"])}),
        ("FontSize", {"fontsize": 16}),
        ("FontFamily", {"fontfamily": "Helvetica"}),
        ("ImageSize: 320 x 240", {}),
    ]
    for axis, (title, options) in zip(axes.flat, cases):
        bandplot(
            ham,
            path,
            rules,
            subdivisions=30,
            ax=axis,
            show=False,
            **options,
        )
        axis.set_title(title, fontsize=12, pad=8)
    axes.flat[-1].axis("off")
    figure.suptitle("bandplot output examples", fontsize=16, fontweight="bold")
    save_svg(figure, "bandplot-examples.svg")


def calculate_fitting_data():
    core_input, _, group_input = build_catalog_core_input(
        space_group=156,
        magnetic_group="156.50",
        k_name="K",
        representation_family="double",
        representations=[{"index": 2}, {"index": 3}],
        coordinate_basis="fractional",
    )
    ham = kpHam(core_input).getkpHam(2, group_input=group_input)["ham"]
    center = (Rational(1, 3), Rational(1, 3), 0)
    eigdata = read_vasp_eigenval(
        EIGENVAL_PATH,
        efermi=-0.0072,
        spin=1,
        start_band=21,
        end_band=22,
    )
    indices = np.arange(50, 71, 5)
    selected = [eigdata[index - 1] for index in indices]
    fit = fittingKP(
        ham,
        EIGENVAL_PATH,
        indices,
        efermi=-0.0072,
        spin=1,
        start_band=21,
        kpoint_center=center,
    )

    def predictions(parameter_rules):
        result = []
        for kpoint, _ in selected:
            local_point = tuple(
                coordinate - float(center_coordinate)
                for coordinate, center_coordinate in zip(kpoint, center)
            )
            matrix = ham.subs(parameter_rules).subs(
                {kx: local_point[0], ky: local_point[1], kz: local_point[2]}
            )
            result.extend(
                np.linalg.eigvalsh(np.asarray(matrix.evalf(), np.complex128))
            )
        return np.asarray(result)

    target = np.asarray(
        [energy for _, energies in selected for energy in energies]
    )
    return (
        indices,
        target.reshape(-1, 2),
        predictions(fit["initial_parameters"]).reshape(-1, 2),
        predictions(fit["fitted_parameters"]).reshape(-1, 2),
        fit["fitted_least_squares"],
    )


def generate_fitting_figure():
    x, reference, initial, fitted, least_squares_error = calculate_fitting_data()
    figure, axis = plt.subplots(figsize=(8, 4.8), constrained_layout=True)
    for band_index in range(2):
        axis.scatter(
            x,
            reference[:, band_index],
            color="#5269e8",
            edgecolor="white",
            linewidth=0.8,
            s=58,
            zorder=3,
            label="VASP" if band_index == 0 else None,
        )
        axis.plot(
            x,
            initial[:, band_index],
            "--",
            color="#8d97aa",
            linewidth=1.7,
            label="Initial k.p" if band_index == 0 else None,
        )
        axis.plot(
            x,
            fitted[:, band_index],
            color="#e25562",
            linewidth=2.1,
            label="Fitted k.p" if band_index == 0 else None,
        )
    axis.set_title("fittingKP comparison", fontsize=15, fontweight="bold")
    axis.set_xlabel("k-point index")
    axis.set_ylabel("Energy (eV)")
    axis.set_xticks(x)
    axis.grid(axis="y", color="#e4e8f1", linewidth=0.8)
    axis.legend(frameon=False, ncol=3, loc="upper center")
    axis.text(
        0.98,
        0.04,
        "least squares = {:.6f}".format(least_squares_error),
        transform=axis.transAxes,
        ha="right",
        color="#5d6980",
        fontsize=9,
    )
    save_svg(figure, "fittingkp-comparison.svg")


def main():
    if not EIGENVAL_PATH.is_file():
        raise FileNotFoundError(EIGENVAL_PATH)
    generate_bandplot_figure()
    generate_fitting_figure()
    print(OUTPUT_DIRECTORY / "bandplot-examples.svg")
    print(OUTPUT_DIRECTORY / "fittingkp-comparison.svg")


if __name__ == "__main__":
    main()
