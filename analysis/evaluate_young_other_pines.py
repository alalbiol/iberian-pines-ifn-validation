"""Operational transfer check of published young-afforestation pine models.

Documented planted IFN4 stands with direct ages are screened to the source
age and density ranges. IFN height uses trees >=7.5 cm, whereas the source
uses plantation mean height including smaller seedlings; biomass references
also use different tree allometry. This is not an isolated source-equation
test. IFN total age minus one nursery year is an explicit assumption.
"""

from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "use_case_1_simulator/agb_forecast_analysis"))

from evaluate_other_pine_stands import FIT_CODES  # noqa: E402
from young_other_pines import PARAMETERS, aboveground_biomass_kg_ha, site_index_30  # noqa: E402


def metrics(frame: pd.DataFrame, prediction: str) -> dict:
    error = frame[prediction]-frame.agb_mg_ha
    return {"plots": len(frame), "bias_mg_ha": round(float(error.mean()),3),
            "mae_mg_ha": round(float(error.abs().mean()),3),
            "rmse_mg_ha": round(float(np.sqrt(np.mean(error**2))),3),
            "median_model_over_observed": round(float((frame[prediction]/frame.agb_mg_ha).median()),3)}


def evaluate() -> dict:
    data = pd.read_csv(Path(__file__).parent / "results/other_pines_ifn4_stands.csv",
                       dtype={"plot_id": str, "province_code": str})
    data.province_code = data.province_code.str.zfill(2)
    data["planting_age"] = data.Edad-1.0
    output = {"design": "Documented planted IFN4, direct age, source age/density range; age=IFN Edad-1 assumption; site index from >=7.5cm tree mean height or fit-province median; IFN reference AGB uses different allometry",
              "species": {}}
    for species,p in PARAMETERS.items():
        frame = data[(data.species == species.replace("_"," ").capitalize()) &
                     data.direct_age & data.planted & data.pine_dominated_measured &
                     data.planting_age.between(*p.age_range) &
                     data.pine_stems_ha.between(*p.density_range) &
                     data.mean_large_h_m.gt(1.3) & data.agb_mg_ha.gt(0)].copy()
        if frame.empty:
            continue
        frame["measured_si30"] = [site_index_30(species,h,age) for h,age in
                                  zip(frame.mean_large_h_m,frame.planting_age)]
        fit = frame[frame.province_code.isin(FIT_CODES[species])].copy()
        held = frame[~frame.province_code.isin(FIT_CODES[species])].copy()
        if len(fit) < 5 or len(held) < 5:
            raise RuntimeError(f"Insufficient young fit/test stands for {species}: {len(fit)}/{len(held)}")
        fit_si = float(fit.measured_si30.median())
        held["default_site_agb"] = [aboveground_biomass_kg_ha(
            species,age,n,fit_si)/1000.0 for age,n in
            zip(held.planting_age,held.pine_stems_ha)]
        held["measured_site_agb"] = [aboveground_biomass_kg_ha(
            species,age,n,si)/1000.0 for age,n,si in
            zip(held.planting_age,held.pine_stems_ha,held.measured_si30)]
        held["fit_median_agb"] = float(fit.agb_mg_ha.median())
        output["species"][species] = {
            "fit_plots": len(fit), "heldout_plots": len(held),
            "fit_provinces": list(FIT_CODES[species]),
            "heldout_by_province": held.province_code.value_counts().to_dict(),
            "heldout_planting_age_range": [float(held.planting_age.min()),
                                            float(held.planting_age.max())],
            "fit_median_si30_from_large_trees_m": fit_si,
            "heldout_results": {name: metrics(held,name) for name in
                                ("default_site_agb","measured_site_agb","fit_median_agb")},
        }
    return output


if __name__ == "__main__":
    import json
    print(json.dumps(evaluate(),indent=2,ensure_ascii=False))
