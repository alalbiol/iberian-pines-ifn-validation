"""Provisional age-based AGB stand chains for four other Spanish pines.

Fit only site-index default, H0→D0, mean/dominant AGB ratio, and density on
designated provinces. Published height shape and Ruiz-Peinado allometry stay
fixed. Held-out provinces assess a cross-sectional inventory snapshot, not
a planting-to-harvest forecast or an independent allometry validation.
"""

from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "models"))

from evaluate_other_pine_height_remeasurement import PINEA_REGION  # noqa: E402
from other_pine_dominant_height import REFERENCE_AGES, project_height  # noqa: E402
from other_pine_stand import (  # noqa: E402
    LARGE_STEMS_QUARTILES,
    PARAMETERS,
    SITE_INDEX_QUARTILES,
    estimate_stand,
)
from other_pine_tree_growth import BIOMASS  # noqa: E402

FIT_CODES = {"pinus_sylvestris": ("27",), "pinus_pinea": ("06",),
             "pinus_halepensis": ("02", "16"), "pinus_nigra": ("25",)}
SOURCE_MIN_AGE = {"pinus_sylvestris": 10, "pinus_pinea": 20,
                  "pinus_halepensis": 10, "pinus_nigra": 30}


def metrics(frame: pd.DataFrame, predicted: str, observed: str) -> dict:
    error = frame[predicted]-frame[observed]
    result = {"plots": len(frame), "bias": round(float(error.mean()), 3),
              "mae": round(float(error.abs().mean()), 3),
              "rmse": round(float(np.sqrt(np.mean(error**2))), 3)}
    if observed == "agb_mg_ha":
        result["median_model_over_observed"] = round(float((frame[predicted]/frame[observed]).median()), 3)
    return result


def evaluate() -> dict:
    data = pd.read_csv(Path(__file__).parent / "results/other_pines_ifn4_stands.csv",
                       dtype={"plot_id": str, "province_code": str})
    data.province_code = data.province_code.str.zfill(2)
    data = data[data.direct_age & data.pine_dominated_measured & data.h0_m.gt(1.3) &
                data.d0_cm.ge(7.5) & data.agb_mg_ha.gt(0)].copy()
    output = {"design": "IFN4 direct-age pine-dominated stands; published height shape; fit SI/D(H)/AGB ratio/density in stated province(s), test all other source-region plots; reference AGB uses same allometry",
              "species": {}}
    for species in FIT_CODES:
        frame = data[data.species == species.replace("_", " ").capitalize()].copy()
        frame = frame[frame.Edad >= SOURCE_MIN_AGE[species]]
        if species == "pinus_sylvestris":
            frame = frame[frame.province_code.isin(("15", "27", "32", "36"))]
        if species == "pinus_pinea":
            frame = frame[frame.province_code.isin(PINEA_REGION)]
        frame["pinea_region"] = frame.province_code.map(PINEA_REGION)
        fit = frame[frame.province_code.isin(FIT_CODES[species])].copy()
        held = frame[~frame.province_code.isin(FIT_CODES[species])].copy()
        if len(fit) < 10 or len(held) < 5:
            raise RuntimeError(f"Insufficient age-based fit/test stands for {species}: {len(fit)}/{len(held)}")
        fit["si"] = [project_height(species, h, age, REFERENCE_AGES[species],
                                     pinea_region=region)
                     for h, age, region in zip(fit.h0_m, fit.Edad, fit.pinea_region)]
        default_si = float(fit.si.median())
        default_n = float(fit.pine_stems_ha.median())
        x = np.log(fit.h0_m.to_numpy()/10.0)
        y = np.log(fit.d0_cm.to_numpy())
        fit_hd = least_squares(lambda b: b[0]+b[1]*x-y,
                               x0=[np.log(float(fit.d0_cm.median())), 1.0],
                               bounds=([0, -2], [6, 4]), loss="soft_l1", f_scale=.15)
        if not fit_hd.success:
            raise RuntimeError(f"H→D fit failed for {species}")
        dominant_bio = np.array([BIOMASS[species](d,h).aboveground_kg
                                 for d,h in zip(fit.d0_cm, fit.h0_m)])
        ratio = float(np.median(fit.agb_mg_ha.to_numpy()*1000.0 /
                                (fit.pine_stems_ha.to_numpy()*dominant_bio)))
        if not 0 < ratio < 1.5:
            raise RuntimeError(f"Implausible biomass ratio for {species}: {ratio}")
        frozen = PARAMETERS[species]
        if not np.allclose((default_si, default_n, *fit_hd.x, ratio),
                           (frozen.default_site_index_m, frozen.default_large_stems_ha,
                            frozen.log_diameter_at_10m, frozen.log_height_exponent,
                            frozen.mean_to_dominant_agb_ratio), atol=1e-7):
            raise RuntimeError(f"Frozen stand parameters differ from fit for {species}")
        site_quartiles = tuple(float(fit.si.quantile(q)) for q in (.25, .5, .75))
        density_quartiles = tuple(float(fit.pine_stems_ha.quantile(q)) for q in (.25, .5, .75))
        if not (np.allclose(site_quartiles, SITE_INDEX_QUARTILES[species], atol=1e-7) and
                np.allclose(density_quartiles, LARGE_STEMS_QUARTILES[species], atol=1e-7)):
            raise RuntimeError(f"Frozen fit-stand quartiles differ from fit for {species}")
        defaults = [estimate_stand(species, age, pinea_region=region)
                    for age,region in zip(held.Edad,held.pinea_region)]
        held["default_h0_m"] = [p.dominant_height_m for p in defaults]
        held["default_agb"] = [p.aboveground_biomass_mg_ha for p in defaults]
        held["measured_hn_agb"] = [estimate_stand(
            species, age, pinea_region=region, dominant_height_m=h, large_stems_ha=n
        ).aboveground_biomass_mg_ha for age,region,h,n in zip(
            held.Edad,held.pinea_region,held.h0_m,held.pine_stems_ha)]
        held["fit_median_agb"] = float(fit.agb_mg_ha.median())
        output["species"][species] = {
            "fit_plots": len(fit), "heldout_plots": len(held),
            "fit_provinces": list(FIT_CODES[species]),
            "heldout_by_province": held.province_code.value_counts().to_dict(),
            "heldout_age_range": [float(held.Edad.min()), float(held.Edad.max())],
            "default_site_index_m": default_si,
            "reference_age_years": REFERENCE_AGES[species],
            "default_large_stems_ha": default_n,
            "fit_site_index_quartiles_m": list(site_quartiles),
            "fit_large_stems_quartiles_ha": list(density_quartiles),
            "height_to_diameter_log_parameters": [float(v) for v in fit_hd.x],
            "fit_mean_to_dominant_agb_ratio": ratio,
            "heldout_height": metrics(held, "default_h0_m", "h0_m"),
            "heldout_agb": {name: metrics(held, name, "agb_mg_ha")
                            for name in ("default_agb", "measured_hn_agb", "fit_median_agb")},
        }
    return output


if __name__ == "__main__":
    import json
    print(json.dumps(evaluate(), indent=2, ensure_ascii=False))
