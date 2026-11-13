"""Exploratory IFN4 test of the mature Mediterranean Pinaster AGB chain.

The same Madrid/Salamanca fit and other-province split as the height study
is used. The held-out residuals informed the later density term, so this
set is no longer untouched for the full chain. Test A uses fit-province
default site and standing
density; test B supplies measured H0 and large-stem count to isolate the
height→diameter→allometry→aggregation chain. IFN AGB uses the same
Ruiz-Peinado allometry, so this does not validate allometry independently.
"""

from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "models"))

from evaluate_mediterranean_pinaster import FIT_PROVINCES, load_eligible  # noqa: E402
from mediterranean_pinaster import site_index_70  # noqa: E402
from mediterranean_pinaster_stand import (  # noqa: E402
    DEFAULT_LARGE_STEMS_HA, DEFAULT_SITE_INDEX_70_M, estimate_stand_agb,
)


def metrics(df: pd.DataFrame, column: str) -> dict:
    error = df[column] - df.agb_mg_ha
    ratio = df[column] / df.agb_mg_ha
    return {
        "n": len(df), "bias_mg_ha": round(float(error.mean()), 3),
        "mae_mg_ha": round(float(error.abs().mean()), 3),
        "rmse_mg_ha": round(float(np.sqrt(np.mean(error**2))), 3),
        "median_model_over_ifn": round(float(ratio.median()), 3),
        "ratio_iqr": [round(float(ratio.quantile(.25)), 3), round(float(ratio.quantile(.75)), 3)],
    }


def evaluate() -> dict:
    df = load_eligible()
    fit = df[df.region.isin(FIT_PROVINCES)]
    held = df[~df.region.isin(FIT_PROVINCES)].copy()
    fit_si = float(np.median([site_index_70(h, age) for h, age in zip(fit.h0_m, fit.Edad)]))
    fit_n = float(fit["trees_per_ha_ge7.5"].median())
    if not np.isclose(DEFAULT_SITE_INDEX_70_M, fit_si, atol=.001):
        raise RuntimeError("hard-coded site default differs from calibration")
    if not np.isclose(DEFAULT_LARGE_STEMS_HA, fit_n, atol=.001):
        raise RuntimeError("hard-coded density default differs from calibration")
    for name, use_h, use_n in (("defaults", False, False),
                               ("measured_height", True, False),
                               ("measured_density", False, True),
                               ("measured_height_and_density", True, True)):
        held[name] = [estimate_stand_agb(age,
                      dominant_height_m=h if use_h else None,
                      large_stems_per_ha=n if use_n else None).aboveground_biomass_mg_ha
                      for age, h, n in zip(held.Edad, held.h0_m, held["trees_per_ha_ge7.5"])]
    return {
        "design": "Mature pooled 2007 height + Mediterranean IFN4 structure + shared allometry; fit Madrid/Salamanca, test other provinces",
        "fit_n": len(fit), "heldout_n": len(held),
        "fit_site_index_70_m": round(fit_si, 3),
        "fit_large_stems_per_ha": round(fit_n, 3),
        "heldout_by_province": held.region.value_counts().to_dict(),
        "heldout_results": {name: metrics(held, name) for name in
                            ("defaults", "measured_height", "measured_density", "measured_height_and_density")},
        "heldout_province_results": {province: {"defaults": metrics(g, "defaults"),
                                                  "measured_height_and_density": metrics(g, "measured_height_and_density")}
                                     for province, g in held.groupby("region")},
    }


if __name__ == "__main__":
    import json
    print(json.dumps(evaluate(), indent=2, ensure_ascii=False))
