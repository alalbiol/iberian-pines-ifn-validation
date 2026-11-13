"""Exploratory external transfer check of Bravo et al. (2008) pine ingrowth.

IFN4 field guide annex 39 identifies OrdenIf3=000 as a newly sampled tree.
We score its IFN4 basal area per five elapsed years against the published
five-year equation. This interval normalization is approximate: new trees
may have grown after entering the sample, and events need not be linear in
time. The result can reject simple transfer but cannot validate an exact
15-year ingrowth trajectory.
"""

from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "use_case_1_simulator/agb_forecast_analysis"))

from fit_mediterranean_density import PROVINCES, load_remeasurements  # noqa: E402
from evaluate_mediterranean_pinaster import FIT_PROVINCES  # noqa: E402
from ifn4_age_validation import IFN4_DIR, _plot_key  # noqa: E402
from ifn_remeasurement_validation import _mdb, _num, trees_per_ha  # noqa: E402
from mediterranean_pinaster_ingrowth import (  # noqa: E402
    SOURCE_BA_RANGE_M2_HA, SOURCE_INTERVAL_YEARS, SOURCE_QMD_RANGE_CM,
    expected_basal_area_m2_ha, occurrence_probability, threshold_basal_area_m2_ha,
)
from spain_biomass_equations import pinus_pinaster_biomass  # noqa: E402


def load_new_trees() -> pd.DataFrame:
    frames = []
    for code, (_, database) in PROVINCES.items():
        path = Path(IFN4_DIR) / f"{database}.accdb"
        trees = _mdb(str(path), "PCMayores")
        trees["plot_id"] = _plot_key(trees, code)
        trees = trees[trees.plot_id.str.startswith(f"{code}_")].copy()
        trees["species"] = _num(trees.Especie)
        trees["ord3"] = _num(trees.OrdenIf3)
        trees["ord4"] = _num(trees.OrdenIf4)
        trees["d_cm"] = (_num(trees.Dn1) + _num(trees.Dn2)) / 20.0
        trees["h_m"] = _num(trees.Ht)
        trees = trees[(trees.species == 26) & (trees.ord3 == 0) &
                      (trees.ord4 > 0) & ~trees.ord4.isin((444, 666, 777, 888, 999)) &
                      (trees.d_cm >= 7.5)].copy()
        trees["w"] = trees_per_ha(trees.d_cm.to_numpy())
        trees["new_ba"] = trees.w * np.pi * (trees.d_cm / 200.0) ** 2
        trees["new_agb"] = [w * pinus_pinaster_biomass(d, h).aboveground_kg / 1000.0
                            if h > 1.3 else 0.0
                            for w, d, h in zip(trees.w, trees.d_cm, trees.h_m)]
        frames.append(trees.groupby("plot_id").agg(new_n_ha=("w", "sum"),
                                                     new_ba_m2_ha=("new_ba", "sum"),
                                                     new_agb_mg_ha=("new_agb", "sum")))
    return pd.concat(frames).groupby(level=0).sum().reset_index()


def metrics(df: pd.DataFrame, prediction: np.ndarray) -> dict:
    observed = df.new_ba_per_five_years.to_numpy()
    error = prediction - observed
    return {"n": len(df), "bias_m2_ha": round(float(error.mean()), 3),
            "mae_m2_ha": round(float(np.abs(error).mean()), 3),
            "rmse_m2_ha": round(float(np.sqrt(np.mean(error**2))), 3)}


def stem_metrics(df: pd.DataFrame, prediction: np.ndarray) -> dict:
    observed = df.new_n_ha.to_numpy()
    error = prediction - observed
    return {"n": len(df), "bias_stems_ha": round(float(error.mean()), 3),
            "mae_stems_ha": round(float(np.abs(error).mean()), 3),
            "rmse_stems_ha": round(float(np.sqrt(np.mean(error**2))), 3)}


def evaluate() -> dict:
    plots = load_remeasurements().merge(load_new_trees(), on="plot_id", how="left", validate="one_to_one")
    plots[["new_n_ha", "new_ba_m2_ha"]] = plots[["new_n_ha", "new_ba_m2_ha"]].fillna(0.0)
    plots = plots[plots.ba3.between(*SOURCE_BA_RANGE_M2_HA) &
                  plots.qmd_cm3.between(*SOURCE_QMD_RANGE_CM)].copy()
    plots["new_ba_per_five_years"] = plots.new_ba_m2_ha * SOURCE_INTERVAL_YEARS / plots.years
    fit = plots[plots.province.isin(FIT_PROVINCES)]
    held = plots[~plots.province.isin(FIT_PROVINCES)]
    if len(fit) < 20 or len(held) < 20:
        raise RuntimeError(f"Insufficient source-domain ingrowth plots: {len(fit)}/{len(held)}")
    published_expected = np.array([expected_basal_area_m2_ha(ba, qmd)
                                   for ba, qmd in zip(held.ba3, held.qmd_cm3)])
    published_threshold = np.array([threshold_basal_area_m2_ha(ba, qmd)
                                    for ba, qmd in zip(held.ba3, held.qmd_cm3)])
    fit_mean = float(fit.new_ba_per_five_years.mean())
    fit_stems_per_new_ba = float(fit.new_n_ha.sum() / fit.new_ba_m2_ha.sum())
    fit_annual_new_stems = float((fit.new_n_ha / fit.years).mean())
    published_new_stems = (published_expected * held.years.to_numpy() /
                           SOURCE_INTERVAL_YEARS * fit_stems_per_new_ba)
    probability_five_year = np.array([occurrence_probability(ba, qmd)
                                      for ba, qmd in zip(held.ba3, held.qmd_cm3)])
    probability_interval = 1.0 - (1.0 - probability_five_year) ** (held.years.to_numpy() / SOURCE_INTERVAL_YEARS)
    observed_occurrence = (held.new_ba_m2_ha.to_numpy() > 0).astype(float)
    fit_occurrence = float((fit.new_ba_m2_ha > 0).mean())
    return {
        "design": "Published central-Spain five-year ingrowth equation versus new IFN4 pine trees identified by OrdenIf3=000; IFN3 BA/QMD source domain; Madrid+Salamanca baseline fit, other provinces test",
        "fit_n": len(fit), "heldout_n": len(held),
        "heldout_by_province": held.province.value_counts().to_dict(),
        "heldout_interval_years": [round(float(held.years.min()), 2), round(float(held.years.max()), 2)],
        "observed_heldout_new_tree_occurrence": round(float(observed_occurrence.mean()), 3),
        "observed_heldout_median_new_stems_ha": round(float(held.new_n_ha.median()), 3),
        "observed_heldout_mean_new_ba_per_five_years_m2_ha": round(float(held.new_ba_per_five_years.mean()), 3),
        "fit_mean_new_ba_per_five_years_m2_ha": round(fit_mean, 3),
        "fit_stems_per_new_basal_area_m2": round(fit_stems_per_new_ba, 3),
        "heldout_results": {
            "published_probability_weighted": metrics(held, published_expected),
            "published_threshold_rule": metrics(held, published_threshold),
            "fit_mean_baseline": metrics(held, np.full(len(held), fit_mean)),
            "zero_ingrowth": metrics(held, np.zeros(len(held))),
        },
        "heldout_new_stem_results": {
            "published_basal_area_with_fit_conversion": stem_metrics(held, published_new_stems),
            "fit_constant_annual_ingrowth": stem_metrics(held, fit_annual_new_stems * held.years.to_numpy()),
            "zero_ingrowth": stem_metrics(held, np.zeros(len(held))),
        },
        "occurrence_brier": {
            "published_independent_five_year_hazards": round(float(np.mean((probability_interval - observed_occurrence)**2)), 3),
            "fit_constant_interval_probability": round(float(np.mean((fit_occurrence - observed_occurrence)**2)), 3),
        },
    }


if __name__ == "__main__":
    import json
    print(json.dumps(evaluate(), indent=2, ensure_ascii=False))
