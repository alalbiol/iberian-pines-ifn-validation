"""Published height-curve transfer on strict same-tree dominant cohorts.

Uses direct IFN4 stand ages and an IFN3→IFN4 linked top-100/ha survivor
cohort. This is a practical trajectory check, not repeated whole-stand H0:
the future survivors were selected by linkage. P. pinea region labels are
geographic approximations for provinces outside the four source regions.
"""

from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "use_case_1_simulator/agb_forecast_analysis"))

from fit_other_pine_tree_growth import SPECIES, load_pairs  # noqa: E402
from ifn_remeasurement_validation import trees_per_ha  # noqa: E402
from other_pine_dominant_height import project_height  # noqa: E402
from stand_structure_calibration import dominant_mean  # noqa: E402
from young_other_pines import PARAMETERS as YOUNG_PARAMETERS, project_mean_height  # noqa: E402

PINEA_REGION = {
    "06": "west_andalusia", "11": "west_andalusia", "14": "west_andalusia",
    "29": "west_andalusia", "05": "central_range", "13": "central_range",
    "16": "central_range", "28": "central_range", "45": "central_range",
    "37": "northern_plateau", "49": "northern_plateau",
    "08": "catalonia", "17": "catalonia", "25": "catalonia", "43": "catalonia",
}


def metrics(frame: pd.DataFrame, prediction: str) -> dict:
    error = frame[prediction]-frame.h4_m
    return {"plots": len(frame), "bias_m": round(float(error.mean()), 3),
            "mae_m": round(float(error.abs().mean()), 3),
            "rmse_m": round(float(np.sqrt(np.mean(error**2))), 3)}


def evaluate() -> dict:
    stands = pd.read_csv(Path(__file__).parent / "results/other_pines_ifn4_stands.csv",
                         dtype={"plot_id": str, "province_code": str})
    stands.province_code = stands.province_code.str.zfill(2)
    stands = stands[stands.direct_age & stands.pine_dominated_measured].copy()
    pairs = load_pairs()
    pairs = pairs.merge(stands[["plot_id", "Edad", "species", "province_code", "planted"]],
                        on="plot_id", suffixes=("", "_stand"), validate="many_to_one")
    pairs = pairs[pairs.species_code.map(SPECIES) == pairs.species]
    pairs["w"] = trees_per_ha(pairs.d3_cm.to_numpy())
    rows = []
    young_rows = []
    for plot_id, group in pairs.groupby("plot_id"):
        if len(group) < 2 or group.w.sum() < 100:
            continue
        row = group.iloc[0]
        age4 = float(row.Edad)
        years = float(group.years.median())
        age3 = age4-years
        if age3 <= 0:
            continue
        species = row.species.lower().replace(" ", "_")
        province = row.province_code
        young = YOUNG_PARAMETERS[species]
        if (row.planted and young.age_range[0] <= age3-1 <= young.age_range[1]
                and young.age_range[0] <= age4-1 <= young.age_range[1]):
            mean3 = float(np.average(group.h3_m,weights=group.w))
            mean4 = float(np.average(group.h4_m,weights=group.w))
            pred = project_mean_height(species,mean3,age3-1,age4-1)
            young_rows.append((species,province,plot_id,mean3,mean4,pred))
        if species == "pinus_sylvestris" and province not in ("15", "27", "32", "36"):
            continue  # Galician source model only
        if species == "pinus_pinea" and (age3 < 20 or province not in PINEA_REGION):
            continue
        if species == "pinus_nigra" and age3 < 30:
            continue
        h3 = dominant_mean(group.h3_m.to_numpy(), group.d3_cm.to_numpy(), group.w.to_numpy())
        h4 = dominant_mean(group.h4_m.to_numpy(), group.d3_cm.to_numpy(), group.w.to_numpy())
        h_pred = project_height(species, h3, age3, age4,
                                pinea_region=PINEA_REGION.get(province))
        rows.append((species, province, plot_id, age3, age4, bool(row.planted), h3, h4, h_pred))
    result = pd.DataFrame(rows, columns=("species", "province", "plot_id", "age3", "age4",
                                         "planted", "h3_m", "h4_m", "published_h4_m"))
    if result.empty:
        raise RuntimeError("No eligible direct-age remeasured cohorts")
    young_result = pd.DataFrame(young_rows,columns=("species","province","plot_id",
                                                   "mean3","mean4","young_pred"))
    def young_metrics(group: pd.DataFrame, column: str) -> dict:
        error = group[column]-group.mean4
        return {"plots": len(group), "bias_m": round(float(error.mean()),3),
                "mae_m": round(float(error.abs().mean()),3),
                "rmse_m": round(float(np.sqrt(np.mean(error**2))),3)}
    return {"design": "Strict linked top-100/ha survivor cohort; direct IFN4 age minus actual interval; source-region/age screens; no fit parameters",
            "species": {name: {
                "plots": len(group), "planted": int(group.planted.sum()),
                "by_province": group.province.value_counts().to_dict(),
                "age3_range": [round(float(group.age3.min()), 2), round(float(group.age3.max()), 2)],
                "published_curve": metrics(group, "published_h4_m"),
                "no_growth": metrics(group.assign(no_growth=group.h3_m), "no_growth"),
            } for name, group in result.groupby("species")},
            "young_afforestation_survivor_cohorts": {
                name: {"plots": len(group),"by_province": group.province.value_counts().to_dict(),
                       "published": young_metrics(group,"young_pred"),
                       "no_growth": young_metrics(group.assign(no_growth=group.mean3),"no_growth")}
                for name,group in young_result.groupby("species")}}


if __name__ == "__main__":
    import json
    print(json.dumps(evaluate(), indent=2, ensure_ascii=False))
