"""Figures for the Iberian pines manuscript.

Inputs: ../results/plot_level/*.csv (written by export_plot_level.py) and
data/spain_provinces.geojson (Natural Earth 1:10m admin-1, public domain,
mainland Spain and Balearic Islands, reprojected to ETRS89 / UTM 30N).
Outputs: fig*.pdf (vector, 18 cm wide) and fig*.png (300 dpi) in this folder.
"""

from pathlib import Path

import geopandas as gpd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "results" / "plot_level"
W = 18 / 2.54  # journal minimum width, inches

INK, INK2, MUTED, AXIS, HAIR = "#0b0b0b", "#52514e", "#898781", "#c3c2b7", "#e1e0d9"
BLUE, ORANGE, NEUTRAL = "#2a78d6", "#eb6834", "#f0efec"
mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Liberation Sans", "Arial", "DejaVu Sans"],
    "font.size": 8, "axes.titlesize": 8.5, "axes.labelsize": 8, "xtick.labelsize": 7,
    "ytick.labelsize": 7, "legend.fontsize": 7.5, "axes.edgecolor": AXIS, "axes.linewidth": 0.6,
    "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK2, "ytick.labelcolor": INK2,
    "axes.labelcolor": INK2, "text.color": INK, "axes.spines.top": False,
    "axes.spines.right": False, "axes.grid": False, "pdf.fonttype": 42, "savefig.dpi": 300,
    "mathtext.fontset": "custom", "mathtext.it": "Liberation Sans:italic",
    "mathtext.rm": "Liberation Sans",
})

ISO_TO_INE = {"AB": "02", "AV": "05", "BA": "06", "B": "08", "CC": "10", "CA": "11", "CR": "13",
              "CO": "14", "C": "15", "CU": "16", "GI": "17", "GU": "19", "J": "23", "L": "25",
              "LU": "27", "M": "28", "MA": "29", "OR": "32", "PO": "36", "SA": "37", "T": "43",
              "TO": "45", "ZA": "49"}
NAME_TO_INE = {"Ourense": "32", "Pontevedra": "36", "Ciudad Real": "13", "Toledo": "45",
               "Guadalajara": "19", "Cáceres": "10", "Badajoz": "06", "Albacete": "02"}
CASES = ["Pinus pinaster (Atlantic)", "Pinus pinaster (Mediterranean)", "Pinus sylvestris",
         "Pinus pinea", "Pinus halepensis", "Pinus nigra"]
LABEL = {"Pinus pinaster (Atlantic)": "P. pinaster, Atlantic",
         "Pinus pinaster (Mediterranean)": "P. pinaster, Mediterranean",
         "Pinus sylvestris": "P. sylvestris", "Pinus pinea": "P. pinea",
         "Pinus halepensis": "P. halepensis", "Pinus nigra": "P. nigra"}
FIT_AGE = {CASES[0]: ("15", "27"), CASES[1]: ("28", "37"), CASES[2]: ("27",),
           CASES[3]: ("06",), CASES[4]: ("02", "16"), CASES[5]: ("25",)}
FIT_TREE = {CASES[1]: ("28", "37"), CASES[2]: ("25", "28"), CASES[3]: ("06", "28"),
            CASES[4]: ("02", "08"), CASES[5]: ("16", "25")}


def code(value) -> str:
    text = str(value)
    return NAME_TO_INE.get(text, text.split(".")[0].zfill(2))


def load():
    spain = gpd.read_file(HERE / "data" / "spain_provinces.geojson")
    spain["ine"] = spain.iso_3166_2.str[3:].map(ISO_TO_INE)
    age = pd.read_csv(DATA / "age_only_forecasts.csv")
    tree = pd.read_csv(DATA / "measured_tree_forecasts.csv")
    age["ine"], tree["ine"] = age.province_code.map(code), tree.province_code.map(code)
    return spain, age, tree


def save(fig, name):
    fig.savefig(HERE / f"{name}.pdf")
    fig.savefig(HERE / f"{name}.png")
    plt.close(fig)


def italic(label: str) -> str:
    species, _, rest = label.partition(",")
    spaced = species.replace(" ", "\\ ")
    return "$\\mathit{" + spaced + "}$" + (f",{rest}" if rest else "")


def two_line(label: str) -> str:
    species, _, rest = label.partition(", ")
    return italic(species) + "\n" + rest


def rmse(error) -> float:
    return float(np.sqrt(np.mean(np.asarray(error) ** 2)))


def base_map(ax, spain):
    spain.plot(ax=ax, color="white", edgecolor=AXIS, linewidth=0.3)
    ax.set_axis_off()
    ax.set_aspect("equal")


def figure_design(spain, age, tree):
    fig, axes = plt.subplots(2, 6, figsize=(W, 2.95))
    held = spain[spain.ine.notna()]
    for row, (frame, fit) in enumerate(((age, FIT_AGE), (tree, FIT_TREE))):
        for col, case in enumerate(CASES):
            ax = axes[row, col]
            base_map(ax, spain)
            if row == 0:
                ax.set_title(two_line(LABEL[case]), fontsize=7.5, pad=2, linespacing=1.15)
            if case not in fit:
                ax.text(0.5, 0.5, "not\nanalysed", transform=ax.transAxes, ha="center",
                        va="center", color=INK2, fontsize=7.5)
                continue
            held.plot(ax=ax, color=HAIR, edgecolor="white", linewidth=0.3)
            test = set(frame.loc[frame.case == case, "ine"]) - set(fit[case])
            held[held.ine.isin(test)].plot(ax=ax, color=BLUE, edgecolor="white", linewidth=0.3)
            held[held.ine.isin(fit[case])].plot(ax=ax, color=ORANGE, edgecolor="white", linewidth=0.3)
    for row, text in enumerate(("Age-only task", "Measured-tree task")):
        axes[row, 0].text(-0.08, 0.5, text, transform=axes[row, 0].transAxes, rotation=90,
                          ha="center", va="center", fontsize=8, color=INK)
    fig.legend(handles=[Patch(color=ORANGE, label="Fitting provinces"),
                        Patch(color=BLUE, label="Testing provinces"),
                        Patch(color=HAIR, label="Inventory data held, no eligible plots"),
                        Patch(facecolor="white", edgecolor=AXIS, label="Not included")],
               loc="lower center", ncol=4, frameon=False, handlelength=1.2, columnspacing=1.6)
    fig.subplots_adjust(left=0.035, right=0.995, top=0.9, bottom=0.1, wspace=0.0, hspace=0.0)
    save(fig, "fig1_design")


def panel_grid(n):
    fig, axes = plt.subplots(2, 3, figsize=(W, 4.9))
    axes = axes.ravel()
    for ax in axes[n:]:
        ax.set_axis_off()
    return fig, axes


def figure_age_only(age):
    fig, axes = panel_grid(6)
    for ax, case in zip(axes, CASES):
        g = age[age.case == case]
        lo, hi = 0.5, 2000
        ax.plot([lo, hi], [lo, hi], color=AXIS, linewidth=0.8, zorder=1)
        ax.scatter(g.observed, g["default"], s=11, color=ORANGE, edgecolor="white",
                   linewidth=0.3, zorder=3)
        ax.scatter(g.observed, g.measured_inputs, s=11, color=BLUE, edgecolor="white",
                   linewidth=0.3, zorder=2)
        ax.set(xscale="log", yscale="log", xlim=(lo, hi), ylim=(lo, hi))
        ax.set_aspect("equal")
        ax.set_title(italic(LABEL[case]), loc="left", pad=4)
        ax.text(0.97, 0.05, f"n = {len(g)}\nRMSE {rmse(g['default'] - g.observed):.0f} | "
                f"{rmse(g.measured_inputs - g.observed):.0f}", transform=ax.transAxes,
                ha="right", va="bottom", fontsize=7, color=INK2, linespacing=1.3)
        ax.set_xticks([1, 10, 100, 1000], ["1", "10", "100", "1000"])
        ax.set_yticks([1, 10, 100, 1000], ["1", "10", "100", "1000"])
        ax.minorticks_off()
    for ax in axes[3:]:
        ax.set_xlabel("Observed biomass (Mg ha$^{-1}$)")
    for ax in axes[::3]:
        ax.set_ylabel("Forecast biomass (Mg ha$^{-1}$)")
    fig.legend(handles=[
        Line2D([], [], marker="o", linestyle="", color=ORANGE, markersize=4.5,
               label="Default site and density"),
        Line2D([], [], marker="o", linestyle="", color=BLUE, markersize=4.5,
               label="Measured height and stem number"),
        Line2D([], [], color=AXIS, linewidth=0.8, label="Forecast = observed")],
        loc="lower center", ncol=3, frameon=False, columnspacing=1.8)
    fig.subplots_adjust(left=0.07, right=0.99, top=0.95, bottom=0.14, wspace=0.18, hspace=0.32)
    save(fig, "fig2_age_only")


def figure_measured_tree(tree):
    fig, axes = panel_grid(5)
    for ax, case in zip(axes, CASES[1:]):
        g = tree[tree.case == case]
        hi = 400
        ax.plot([0, hi], [0, hi], color=AXIS, linewidth=0.8, zorder=1)
        ax.scatter(g.agb4, g.agb3, s=7, color=MUTED, alpha=0.45, linewidth=0, zorder=2)
        ax.scatter(g.agb4, g.pred_total, s=7, color=BLUE, alpha=0.7, linewidth=0, zorder=3)
        ax.set(xlim=(0, hi), ylim=(0, hi), xticks=[0, 100, 200, 300, 400],
               yticks=[0, 100, 200, 300, 400])
        ax.set_aspect("equal")
        ax.set_title(italic(LABEL[case]), loc="left", pad=4)
        ax.text(0.97, 0.05, f"n = {len(g)}\nRMSE {rmse(g.pred_total - g.agb4):.1f}\n"
                f"no change {rmse(g.agb3 - g.agb4):.1f}", transform=ax.transAxes,
                ha="right", va="bottom", fontsize=7, color=INK2, linespacing=1.3)
        ax.set_xlabel("Observed biomass at IFN4 (Mg ha$^{-1}$)")
    for ax in axes[:5:3]:
        ax.set_ylabel("Forecast biomass at IFN4 (Mg ha$^{-1}$)")
    axes[5].legend(handles=[
        Line2D([], [], marker="o", linestyle="", color=BLUE, markersize=4.5,
               label="Forecast from the IFN3 tree list"),
        Line2D([], [], marker="o", linestyle="", color=MUTED, markersize=4.5,
               label="No change (biomass at IFN3)"),
        Line2D([], [], color=AXIS, linewidth=0.8, label="Forecast = observed")],
        loc="center left", frameon=False, labelspacing=0.9)
    fig.subplots_adjust(left=0.07, right=0.99, top=0.95, bottom=0.09, wspace=0.18, hspace=0.42)
    save(fig, "fig3_measured_tree")


def figure_bias_map(spain, tree):
    fig, axes = plt.subplots(1, 5, figsize=(W, 2.1))
    cmap = mpl.colors.LinearSegmentedColormap.from_list("bias", [ORANGE, NEUTRAL, BLUE])
    norm = mpl.colors.Normalize(-30, 30)
    for ax, case in zip(axes, CASES[1:]):
        base_map(ax, spain)
        g = tree[tree.case == case].assign(error=lambda d: d.pred_total - d.agb4)
        by = g.groupby("ine").error.agg(["size", "mean"])
        by = by[by["size"] >= 10]
        shown = spain.merge(by, left_on="ine", right_index=True)
        shown.plot(ax=ax, color=[cmap(norm(v)) for v in shown["mean"]], edgecolor=MUTED,
                   linewidth=0.35)
        ax.set_title(two_line(LABEL[case]), fontsize=7.5, pad=2, linespacing=1.15)
    fig.subplots_adjust(left=0.005, right=0.995, top=0.86, bottom=0.27, wspace=0.0)
    cax = fig.add_axes([0.3, 0.2, 0.4, 0.035])
    bar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal",
                       ticks=[-30, -15, 0, 15, 30])
    bar.outline.set_visible(False)
    bar.set_label("Mean forecast error by province (Mg ha$^{-1}$): "
                  "orange, forecast too low; blue, forecast too high")
    save(fig, "fig4_bias_by_province")


if __name__ == "__main__":
    provinces, age_only, measured_tree = load()
    figure_design(provinces, age_only, measured_tree)
    figure_age_only(age_only)
    figure_measured_tree(measured_tree)
    figure_bias_map(provinces, measured_tree)
