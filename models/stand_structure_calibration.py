"""
Calibrates the mean-tree / dominant-tree biomass ratio used by spain_stand_agb.py.

Why: every growth curve we have describes the DOMINANT tree (pine and cork oak: dominant
height/diameter; holm oak: tallest-tree juvenile anchor + IFN adult rate). Multiplying a
dominant tree's biomass by all stems/ha overestimates stand AGB. The ratio

    mean aboveground biomass per tree / aboveground biomass of the dominant trees

is measured here on IFN3 plots (Spanish National Forest Inventory, ~2000), using each tree's
measured d and h in the Ruiz-Peinado equations and the IFN concentric-plot expansion factors.
"Dominant" = Assmann's top height convention: the 100 thickest trees per hectare.

Plot selection: the target species holds >= 80% of the plot's basal area and >= 100 of its
trees/ha (below that every tree is "dominant" and the ratio is 1 by construction).
Regions: Galicia (15, 27, 36) for pine, Extremadura (06, 10) for the oaks, i.e. where
Sylvestris works.

CAVEAT: IFN counts only trees >= 7.5 cm dbh, and its stands are mostly adult and often
uneven-aged, while a plantation is even-aged. Young plantations are probably more uniform
(ratio closer to 1), so applying these ratios to them is conservative (underestimates AGB).

Run: python stand_structure_calibration.py  (prints the table hardcoded in spain_stand_agb.py)
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ifn_remeasurement_validation import _find, _mdb, _num, trees_per_ha  # noqa: E402
from spain_biomass_equations import (  # noqa: E402
    pinus_pinaster_biomass,
    quercus_ilex_biomass,
    quercus_suber_biomass,
)

TARGETS = {
    "pinus_pinaster": (26, ("15", "27", "36"), pinus_pinaster_biomass),
    "quercus_ilex": (45, ("06", "10"), quercus_ilex_biomass),
    "quercus_suber": (46, ("06", "10"), quercus_suber_biomass),
}
MIN_SPECIES_BA_SHARE = 0.8
MIN_SPECIES_TREES_PER_HA = 100.0
DOMINANT_TREES_PER_HA = 100.0


def load_ifn3_trees(code: str) -> pd.DataFrame:
    t = _mdb(_find((os.environ.get("IFN_DATA_ROOT", os.path.expanduser("~/ifn-data")) + f"/ifn3/p{code}/Ifn3p{code}.accdb")), "PCMayores")
    t["plot"] = code + "_" + t.Estadillo.str.strip() + t.Cla.str.strip() + t.Subclase.fillna("").str.strip()
    t["species_code"] = _num(t.Especie)
    t["d_cm"] = (_num(t.Dn1) + _num(t.Dn2)) / 20.0
    t["h_m"] = _num(t.Ht)
    t = t[(t.d_cm >= 7.5) & (t.h_m > 1.3)].copy()
    t["w"] = trees_per_ha(t.d_cm.to_numpy())
    t["ba_m2_ha"] = t.w * np.pi * (t.d_cm / 200.0) ** 2
    return t[["plot", "species_code", "d_cm", "h_m", "w", "ba_m2_ha"]]


def dominant_mean(values: np.ndarray, d: np.ndarray, w: np.ndarray, n_dom: float = DOMINANT_TREES_PER_HA) -> float:
    """Weighted mean of `values` over the n_dom thickest trees/ha (partial weight for the last one)."""
    order = np.argsort(-d)
    v, w = values[order], w[order]
    take = np.minimum(w, np.maximum(n_dom - np.concatenate([[0.0], np.cumsum(w)[:-1]]), 0.0))
    return float((v * take).sum() / take.sum())


def plot_ratios(species: str) -> pd.DataFrame:
    code, provinces, biomass = TARGETS[species]
    trees = pd.concat([load_ifn3_trees(p) for p in provinces], ignore_index=True)
    ba_share = trees[trees.species_code == code].groupby("plot").ba_m2_ha.sum() / trees.groupby("plot").ba_m2_ha.sum()
    sp = trees[trees.species_code == code].copy()
    sp["agb_kg"] = [biomass(d, h).aboveground_kg for d, h in zip(sp.d_cm, sp.h_m)]
    rows = []
    for plot, g in sp.groupby("plot"):
        n = g.w.sum()
        if n < MIN_SPECIES_TREES_PER_HA or ba_share.get(plot, 0.0) < MIN_SPECIES_BA_SHARE:
            continue
        mean_agb = float((g.agb_kg * g.w).sum() / n)
        dom_agb = dominant_mean(g.agb_kg.to_numpy(), g.d_cm.to_numpy(), g.w.to_numpy())
        rows.append({"plot": plot, "trees_per_ha": n, "ratio": mean_agb / dom_agb})
    return pd.DataFrame(rows)


if __name__ == "__main__":
    for species in TARGETS:
        r = plot_ratios(species)
        q = r.ratio.quantile([0.25, 0.5, 0.75])
        print(f"{species:15s} plots={len(r):5d}  ratio median={q[0.5]:.3f}  IQR={q[0.25]:.3f}-{q[0.75]:.3f}")
        bins = pd.cut(r.trees_per_ha, [100, 200, 400, 800, 1600, 1e6])
        print("   by trees/ha:", r.groupby(bins, observed=True).ratio.agg(["median", "size"]).round(3).to_dict("index"))
