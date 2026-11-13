"""Province-held-out survivor growth tests for four Spanish pine species.

Uses strict IFN3→IFN4 same-tree links. Fit provinces were chosen before
examining outcome errors from the two largest usable geographic groups for
each species; every other province is held out. These are conditional-on-
survival tree models, not stand mortality, recruitment, or planting models.
"""

import os  # Configurable external inventory paths.

from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "use_case_1_simulator/agb_forecast_analysis"))
from other_pine_tree_growth import PARAMETERS  # noqa: E402

SPECIES = {21: "Pinus sylvestris", 23: "Pinus pinea",
           24: "Pinus halepensis", 25: "Pinus nigra"}
FIT_CODES = {21: ("25", "28"), 23: ("06", "28"),
             24: ("02", "08"), 25: ("16", "25")}
PAIRS_DIR = Path((os.environ.get("IFN_DATA_ROOT", os.path.expanduser("~/ifn-data")) + "/ifn_validation"))


def load_pairs() -> pd.DataFrame:
    frames = []
    for path in sorted(PAIRS_DIR.glob("pairs34_p*.csv")):
        code = path.stem[-2:]
        d = pd.read_csv(path, dtype={"province": str, "plot_key": str})
        d = d[d.species_code.isin(SPECIES)].copy()
        d["province_code"] = code
        d["plot_id"] = str(int(code)) + "_" + d.plot_key.str.replace("_", "", regex=False)
        frames.append(d)
    df = pd.concat(frames, ignore_index=True)
    df = df[df.d3_cm.between(7.5, 80) & df.d4_cm.between(7.5, 100) &
            df.h3_m.between(1.31, 35) & df.h4_m.between(1.31, 45) &
            df.years.between(5, 25) & df.density_per_ha.between(20, 3000)].copy()
    df["annual_dbh_cm"] = (df.d4_cm - df.d3_cm)/df.years
    df["annual_height_m"] = (df.h4_m - df.h3_m)/df.years
    return df


def predict_rate(df: pd.DataFrame, coefficients: np.ndarray, kind: str) -> np.ndarray:
    x = [np.log(df.d3_cm.to_numpy()/20.0),
         np.log(df.density_per_ha.to_numpy()/500.0)]
    if kind == "height":
        x.insert(0, np.log(df.h3_m.to_numpy()/10.0))
    return np.exp(coefficients[0] + sum(b*v for b, v in zip(coefficients[1:], x)))


def metrics(frame: pd.DataFrame, prediction: np.ndarray, target: str) -> dict:
    error = prediction - frame[target].to_numpy()
    return {"trees": len(frame), "plots": int(frame.plot_id.nunique()),
            "bias": round(float(error.mean()), 3),
            "mae": round(float(np.abs(error).mean()), 3),
            "rmse": round(float(np.sqrt(np.mean(error**2))), 3)}


def fit_component(fit: pd.DataFrame, held: pd.DataFrame, kind: str,
                  species: str) -> dict:
    annual = "annual_dbh_cm" if kind == "diameter" else "annual_height_m"
    initial = "d3_cm" if kind == "diameter" else "h3_m"
    target = "d4_cm" if kind == "diameter" else "h4_m"
    flat = float(np.median(fit.groupby("plot_id")[annual].median()))
    weights = 1.0/np.sqrt(fit.groupby("plot_id").plot_id.transform("size").to_numpy())
    n_parameters = 3 if kind == "diameter" else 4
    result = least_squares(
        lambda b: (predict_rate(fit, b, kind)-fit[annual].to_numpy())*weights,
        x0=np.array([np.log(max(flat, .01))]+[0.0]*(n_parameters-1)),
        bounds=([-5.0]+[-3.0]*(n_parameters-1), [1.0]+[3.0]*(n_parameters-1)),
        loss="soft_l1", f_scale=.1 if kind == "diameter" else .15)
    if not result.success:
        raise RuntimeError(result.message)
    frozen = PARAMETERS[species]
    target_parameters = ((frozen.log_dbh_rate, frozen.dbh_size_exponent,
                          frozen.dbh_density_exponent) if kind == "diameter" else
                         (frozen.log_height_rate, frozen.height_size_exponent,
                          frozen.height_dbh_exponent, frozen.height_density_exponent))
    if not np.allclose(result.x, target_parameters, atol=1e-7):
        raise RuntimeError(f"Frozen {species} {kind} coefficients differ from refit")
    if kind == "height" and not np.isclose(flat, frozen.flat_height_rate_m_year, atol=1e-7):
        raise RuntimeError(f"Frozen {species} flat height increment differs from refit")
    actual = held[target].to_numpy()
    years = held.years.to_numpy()
    initial_value = held[initial].to_numpy()
    preds = {"size_density_model": initial_value + years*predict_rate(held, result.x, kind),
             "fit_constant_rate": initial_value + years*flat,
             "no_growth": initial_value}
    by_province = {}
    for code, group in held.groupby("province_code"):
        pred = group[initial].to_numpy()+group.years.to_numpy()*predict_rate(group, result.x, kind)
        by_province[code] = metrics(group, pred, target)
    return {"fit_trees": len(fit), "fit_plots": int(fit.plot_id.nunique()),
            "heldout_trees": len(held), "heldout_plots": int(held.plot_id.nunique()),
            "fit_constant_annual_rate": flat,
            "size_density_parameters": [float(x) for x in result.x],
            "heldout_results": {name: metrics(held, pred, target) for name, pred in preds.items()},
            "heldout_by_province": by_province,
            "negative_observed_increment_fraction": round(float((fit[annual] < 0).mean()), 4)}


def evaluate() -> dict:
    df = load_pairs()
    output = {"design": "Strict surviving IFN3→IFN4 trees; plot-equal robust fit on prespecified provinces, other provinces held out; DBH cm and height m",
              "species": {}}
    for species_code, species in SPECIES.items():
        group = df[df.species_code == species_code]
        fit = group[group.province_code.isin(FIT_CODES[species_code])].reset_index(drop=True)
        held = group[~group.province_code.isin(FIT_CODES[species_code])].reset_index(drop=True)
        if len(fit) < 100 or len(held) < 100:
            raise RuntimeError(f"Insufficient fit/test pairs for {species}")
        output["species"][species] = {
            "fit_provinces": list(FIT_CODES[species_code]),
            "diameter": fit_component(fit, held, "diameter", species.lower().replace(" ", "_")),
            "height": fit_component(fit, held, "height", species.lower().replace(" ", "_")),
        }
    return output


if __name__ == "__main__":
    import json
    print(json.dumps(evaluate(), indent=2, ensure_ascii=False))
