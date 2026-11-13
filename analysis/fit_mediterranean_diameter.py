"""Fit/test a provisional Mediterranean pine survivor-tree DBH model.

Strict IFN3→IFN4 linked pairs are supplied by ifn3_ifn4_pairs.py. Madrid and
Salamanca are the fit provinces, as in the height analysis; all other locally
available central provinces are held out. The fit weights plots equally and
uses a robust loss. It does not use IFN4 DBH as an input at prediction time.
"""

import os  # Configurable external inventory paths.

from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "models"))

from mediterranean_pinaster_diameter import DiameterGrowthParameters, MEDITERRANEAN_PARAMS  # noqa: E402

PAIRS_DIR = Path((os.environ.get("IFN_DATA_ROOT", os.path.expanduser("~/ifn-data")) + "/ifn_validation"))
FIT_CODES = ("28", "37")
HELD_CODES = ("02", "06", "10", "13", "16", "19", "45")
PROVINCES = {"02": "Albacete", "06": "Badajoz", "10": "Cáceres", "13": "Ciudad Real",
             "16": "Cuenca", "19": "Guadalajara", "28": "Madrid", "37": "Salamanca", "45": "Toledo"}


def load_pairs() -> pd.DataFrame:
    frames = []
    for code in (*FIT_CODES, *HELD_CODES):
        path = PAIRS_DIR / f"pairs34_p{code}.csv"
        if not path.exists():
            raise FileNotFoundError(f"missing strict-pair cache {path}")
        df = pd.read_csv(path, dtype={"province": str, "plot_key": str})
        df = df[df.species_code == 26].copy()
        df["province_code"] = code
        df["plot_id"] = code + "_" + df.plot_key
        frames.append(df)
    df = pd.concat(frames, ignore_index=True)
    # Prespecified input/interval domain; no exclusion based on observed growth.
    df = df[df.d3_cm.between(7.5, 80) & df.d4_cm.between(7.5, 100) &
            df.years.between(5, 25) & df.density_per_ha.between(20, 3000)].copy()
    df["annual_observed_cm"] = (df.d4_cm - df.d3_cm) / df.years
    return df


def predict_rate(df: pd.DataFrame, coefficients: np.ndarray) -> np.ndarray:
    x_d = np.log(df.d3_cm.to_numpy() / 20.0)
    x_n = np.log(df.density_per_ha.to_numpy() / 500.0)
    return np.exp(coefficients[0] + coefficients[1] * x_d + coefficients[2] * x_n)


def metrics(df: pd.DataFrame, prediction: np.ndarray) -> dict:
    error = prediction - df.d4_cm.to_numpy()
    annual_error = error / df.years.to_numpy()
    plot_bias = pd.Series(error).groupby(df.plot_id.reset_index(drop=True)).mean()
    return {
        "trees": len(df), "plots": df.plot_id.nunique(),
        "dbh_bias_cm": round(float(error.mean()), 3),
        "dbh_mae_cm": round(float(np.abs(error).mean()), 3),
        "dbh_rmse_cm": round(float(np.sqrt(np.mean(error**2))), 3),
        "annual_mae_cm_per_year": round(float(np.abs(annual_error).mean()), 3),
        "plot_mean_bias_rmse_cm": round(float(np.sqrt(np.mean(plot_bias**2))), 3),
    }


def evaluate() -> dict:
    df = load_pairs()
    fit = df[df.province_code.isin(FIT_CODES)].copy().reset_index(drop=True)
    held = df[df.province_code.isin(HELD_CODES)].copy().reset_index(drop=True)
    plot_sizes = fit.groupby("plot_id").plot_id.transform("size").to_numpy()
    weights = 1.0 / np.sqrt(plot_sizes)
    observed = fit.annual_observed_cm.to_numpy()
    flat_rate = float(np.median(fit.groupby("plot_id").annual_observed_cm.median()))
    result = least_squares(
        lambda b: (predict_rate(fit, b) - observed) * weights,
        x0=np.array([np.log(max(flat_rate, .01)), 0.0, 0.0]),
        bounds=([-5.0, -3.0, -3.0], [1.0, 3.0, 3.0]),
        loss="soft_l1", f_scale=.1,
    )
    if not result.success:
        raise RuntimeError(result.message)
    p = DiameterGrowthParameters(*map(float, result.x))
    if not np.allclose(result.x, [MEDITERRANEAN_PARAMS.log_rate_at_20cm_500ha,
                                  MEDITERRANEAN_PARAMS.log_diameter_exponent,
                                  MEDITERRANEAN_PARAMS.log_density_exponent], atol=1e-6):
        raise RuntimeError("Frozen Mediterranean diameter parameters differ from refit")
    pred_model = held.d3_cm.to_numpy() + held.years.to_numpy() * predict_rate(held, result.x)
    pred_flat = held.d3_cm.to_numpy() + held.years.to_numpy() * flat_rate
    pred_none = held.d3_cm.to_numpy()
    by_province = {}
    for code, g in held.groupby("province_code"):
        by_province[PROVINCES[code]] = {
            "model": metrics(g.reset_index(drop=True),
                             g.d3_cm.to_numpy() + g.years.to_numpy() * predict_rate(g, result.x)),
            "median_annual_observed_cm": round(float(g.annual_observed_cm.median()), 3),
        }
    return {
        "design": "Strict surviving-tree IFN3→IFN4 pairs; Madrid+Salamanca fit, seven other provinces held out; plot-equal robust fit",
        "fit_trees": len(fit), "fit_plots": fit.plot_id.nunique(),
        "heldout_trees": len(held), "heldout_plots": held.plot_id.nunique(),
        "negative_observed_growth_fraction": round(float((df.annual_observed_cm < 0).mean()), 4),
        "parameters": {"log_rate_at_20cm_500ha": p.log_rate_at_20cm_500ha,
                       "log_diameter_exponent": p.log_diameter_exponent,
                       "log_density_exponent": p.log_density_exponent},
        "constant_rate_baseline_cm_per_year": round(flat_rate, 4),
        "heldout_results": {"diameter_density_model": metrics(held, pred_model),
                            "constant_rate": metrics(held, pred_flat),
                            "no_growth": metrics(held, pred_none)},
        "heldout_by_province": by_province,
    }


if __name__ == "__main__":
    import json
    print(json.dumps(evaluate(), indent=2, ensure_ascii=False))
