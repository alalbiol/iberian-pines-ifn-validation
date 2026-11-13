"""
Recalibration of the pine stand chain on IFN4 Galicia, after ifn4_age_validation.py showed the all-defaults
model at ~0.2x measured AGB on known-age plots. Three pieces, each fitted on plots NOT used to test it:

  1. Dominant height -> dominant diameter. The nationwide h-d inverse gives the diameter of an AVERAGE tree of
     height H0; dominant trees are thicker (measured D0 / inverse = 1/0.78 on the aged plots). Fit the same M4
     form (h = 1.3 + a ln(1+d)^b) on plot-level (D0, H0) pairs, D0/H0 = mean d/h of the 100 thickest trees/ha.
     Fit set: all pine-dominated IFN4 plots in the fit provinces.
  2. Mean/dominant tree AGB ratio for EVEN-AGED stands (FPMasa 1/2), by trees/ha, fit provinces. The IFN3 table
     (all stand structures) runs 0.1-0.15 lower at 300-1000 trees/ha.
  3. Default site: median site index (dominant height at 20 y on the interior curve's shape) of known-age plots
     aged >= 10 y (younger ones carry residual older trees).

Split by province: ALL three are fitted on A Coruña + Lugo only; Ourense + Pontevedra are held out for the
test in ifn4_age_validation.py (TEST_REGIONS there). A split by plot (known-age vs not) was tried first but
leaves only 30 even-aged plots for the ratio, since almost every even-aged plot carries an age.

Run: python ifn4_pine_calibration.py  (prints the constants hardcoded in spain_stand_agb.py)
"""

from __future__ import annotations

import math
import os
import sys

import numpy as np
import pandas as pd
from scipy.optimize import brentq, curve_fit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ifn4_age_validation import DATABASES, IFN4_DIR, _plot_key, load_aged_plots  # noqa: E402
from ifn_remeasurement_validation import _mdb, _num, trees_per_ha  # noqa: E402
from spain_biomass_equations import pinus_pinaster_biomass  # noqa: E402
from spain_growth_curves import pinus_pinaster_project  # noqa: E402
from spain_stand_agb import PLANTING_LAG_YEARS  # noqa: E402
from stand_structure_calibration import dominant_mean  # noqa: E402

GALICIA = ["A Coruña", "Lugo", "Ourense", "Pontevedra"]
FIT_REGIONS = ["A Coruña", "Lugo"]
MIN_BA_SHARE = 0.8
N_BINS = [100, 200, 300, 450, 600, 800, 1000, 1400, 3000]


def pine_plots(region: str, database: str | None = None) -> pd.DataFrame:
    """Per pine-dominated plot: D0, H0, trees/ha, mean/dominant AGB ratio, even-aged flag."""
    db = database or DATABASES[region]
    t = _mdb(os.path.join(IFN4_DIR, f"{db}.accdb"), "PCMayores")
    t["plot_id"] = _plot_key(t)
    t["species"] = _num(t.Especie)
    t["d"] = (_num(t.Dn1) + _num(t.Dn2)) / 20.0
    t["h"] = _num(t.Ht)
    t = t[(t.d >= 7.5) & (t.h > 1.3)].copy()
    t["w"] = trees_per_ha(t.d.to_numpy())
    t["ba"] = t.w * math.pi * (t.d / 200.0) ** 2
    e = _mdb(os.path.join(IFN4_DIR, f"{db}.accdb"), "PCEspParc")
    e["plot_id"] = _plot_key(e)
    e = e[_num(e.PosEsp) == 1]
    even = set(e.plot_id[_num(e.FPMasa).isin([1, 2])])
    share = t[t.species == 26].groupby("plot_id").ba.sum() / t.groupby("plot_id").ba.sum()
    rows = []
    for pid, g in t[t.species == 26].groupby("plot_id"):
        if share.get(pid, 0) < MIN_BA_SHARE:
            continue
        n = g.w.sum()
        agb = np.array([pinus_pinaster_biomass(d, h).aboveground_kg for d, h in zip(g.d, g.h)])
        d0 = dominant_mean(g.d.to_numpy(), g.d.to_numpy(), g.w.to_numpy())
        h0 = dominant_mean(g.h.to_numpy(), g.d.to_numpy(), g.w.to_numpy())
        ratio = (agb * g.w).sum() / n / dominant_mean(agb, g.d.to_numpy(), g.w.to_numpy()) if n >= 100 else np.nan
        rows.append({"plot_id": pid, "region": region, "d0": d0, "h0": h0, "trees_per_ha": n, "ratio": ratio,
                     "even_aged": pid in even})
    return pd.DataFrame(rows)


def hd_model(d, a, b):
    return 1.3 + a * np.log1p(d) ** b


def site_index_20(h: float, total_age: float, coastal: bool = False) -> float:
    """Dominant height at 20 y of the ecoregion curve passing through (total_age, h).

    `total_age` is on the PAPER's clock (age from germination), which is how we read IFN4 `Edad` —
    see the AGE CONVENTION block in spain_stand_agb.py. Before 2026-09-26 this was called with
    `Edad - 1`, one year off the convention the model itself now uses."""
    return pinus_pinaster_project(h, total_age, 20.0, ecoregion_coastal=coastal)


if __name__ == "__main__":
    aged = pd.concat([load_aged_plots(r, DATABASES[r]) for r in FIT_REGIONS])
    plots = pd.concat([pine_plots(r) for r in FIT_REGIONS], ignore_index=True)
    fit = plots
    print(f"Pine-dominated IFN4 plots in {FIT_REGIONS}: {len(fit)}")

    # 1. dominant h-d
    (a, b), _ = curve_fit(hd_model, fit.d0, fit.h0, p0=(1.5, 1.5))
    resid = fit.h0 - hd_model(fit.d0, a, b)
    print(f"\n1. Dominant h-d (D0, H0), n={len(fit)}: a={a:.4f}, b={b:.4f}, RMSE {np.sqrt((resid**2).mean()):.2f} m")
    for h in (8, 12, 16, 20):
        d_dom = brentq(lambda d: hd_model(d, a, b) - h, 0.1, 200)
        print(f"   H0 {h} m -> D0 {d_dom:.1f} cm")

    # 2. even-aged mean/dominant ratio
    ev = fit[fit.even_aged & fit.ratio.notna()]
    g = ev.groupby(pd.cut(ev.trees_per_ha, N_BINS), observed=True)
    table = pd.DataFrame({"n": g.size(), "N_median": g.trees_per_ha.median(), "q25": g.ratio.quantile(0.25),
                          "median": g.ratio.median(), "q75": g.ratio.quantile(0.75)})
    print(f"\n2. Even-aged mean/dominant AGB ratio, n={len(ev)} plots:")
    print(table.round(2).to_string())

    # 3. default site from known-age plots >= 10 y, A Coruña + Lugo only
    ap = aged[aged.Especie == 26].merge(plots[["plot_id", "h0"]], on="plot_id")
    ap["age"] = ap.Edad - PLANTING_LAG_YEARS
    ap = ap[ap.age >= 10]
    si = np.array([site_index_20(h, t) for h, t in zip(ap.h0, ap.Edad)])
    print(f"\n3. Site index (H0 at 20 y, interior curve shape), known-age plots >= 10 y in {FIT_REGIONS}: "
          f"n={len(si)}, median {np.median(si):.2f} m (IQR {np.percentile(si, 25):.2f}-{np.percentile(si, 75):.2f}); "
          f"the paper's mean stand gives {site_index_20(10.64, 19.3):.2f} m")
