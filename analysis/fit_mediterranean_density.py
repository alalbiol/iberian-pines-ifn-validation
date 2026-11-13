"""Explore net Pinus pinaster large-stem density change in IFN3→IFN4 plots.

Whole-plot stems >=7.5 cm are expanded from PCMayores, with pine basal-area
share >=80% in both cycles. A plot must contain at least two strict same-tree
pine links in the existing pair cache. The density change includes ingrowth,
mortality and management; it is not a pure survival estimate. Parameters are
fit in Madrid/Salamanca and scored without refitting in other provinces.
"""

import os  # Configurable external inventory paths.

from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "use_case_1_simulator/agb_forecast_analysis"))

from evaluate_mediterranean_pinaster import FIT_PROVINCES  # noqa: E402
from ifn4_age_validation import IFN4_DIR, _plot_key  # noqa: E402
from ifn_remeasurement_validation import _mdb, _num, trees_per_ha  # noqa: E402

IFN3_DIR = Path((os.environ.get("IFN_DATA_ROOT", os.path.expanduser("~/ifn-data")) + "/ifn3"))
PAIRS_DIR = Path((os.environ.get("IFN_DATA_ROOT", os.path.expanduser("~/ifn-data")) + "/ifn_validation"))
PROVINCES = {
    2: ("Albacete", "Ifn4_Albacete"),
    6: ("Badajoz", "Ifn4_Extremadura"),
    10: ("Cáceres", "Ifn4_Extremadura"),
    13: ("Ciudad Real", "Ifn4_CiudadReal"),
    16: ("Cuenca", "Ifn4p_Cuenca"),
    19: ("Guadalajara", "Ifn4_Guadalajara"),
    28: ("Madrid", "Ifn4p28-madrid"),
    37: ("Salamanca", "Ifn4p37-salamanca"),
    45: ("Toledo", "Ifn4_Toledo"),
}


def plot_stems(path: Path, code: int) -> pd.DataFrame:
    trees = _mdb(str(path), "PCMayores")
    trees["plot_id"] = _plot_key(trees, code)
    trees["species"] = _num(trees.Especie)
    trees["d_cm"] = (_num(trees.Dn1) + _num(trees.Dn2)) / 20.0
    trees["h_m"] = _num(trees.Ht)
    current_order = _num(trees.OrdenIf4 if "OrdenIf4" in trees else trees.OrdenIf3)
    # PCMayores also retains records of trees cut or dead since the prior
    # inventory. Those records may carry old diameters; a current measured
    # total height is required for the live standing-stock tally.
    trees = trees[(trees.d_cm >= 7.5) & (trees.h_m > 1.3) &
                  (current_order > 0) & ~current_order.isin((444, 666, 777, 888, 999))].copy()
    trees["w"] = trees_per_ha(trees.d_cm.to_numpy())
    trees["ba"] = trees.w * np.pi * (trees.d_cm / 200.0) ** 2
    total_ba = trees.groupby("plot_id").ba.sum()
    pine = trees[trees.species == 26].groupby("plot_id").agg(
        n=("w", "sum"), ba=("ba", "sum"))
    pine["share"] = pine.ba / total_ba
    pine["qmd_cm"] = 200.0 * np.sqrt(pine.ba / (np.pi * pine.n))
    return pine.reset_index()


def load_remeasurements() -> pd.DataFrame:
    frames = []
    for code, (region, ifn4_db) in PROVINCES.items():
        p3 = IFN3_DIR / f"p{code:02d}" / f"Ifn3p{code:02d}.accdb"
        p4 = Path(IFN4_DIR) / f"{ifn4_db}.accdb"
        cache = PAIRS_DIR / f"pairs34_p{code:02d}.csv"
        if not (p3.exists() and p4.exists() and cache.exists()):
            continue
        p3_stems = plot_stems(p3, code)
        p4_stems = plot_stems(p4, code)
        matched = p3_stems.merge(p4_stems, on="plot_id", suffixes=("3", "4"), validate="one_to_one")
        pairs = pd.read_csv(cache, usecols=["plot_key", "species_code", "years"])
        pairs = pairs[(pairs.species_code == 26) & pairs.years.between(5, 25)].copy()
        pairs["plot_id"] = str(code) + "_" + pairs.plot_key.str.replace("_", "", regex=False)
        links = pairs.groupby("plot_id").agg(strict_pine_links=("years", "size"), years=("years", "median"))
        matched = matched.merge(links, on="plot_id", validate="one_to_one")
        matched = matched[(matched.strict_pine_links >= 2) &
                          (matched.share3 >= .8) & (matched.share4 >= .8) &
                          (matched.n3 >= 100) & (matched.n4 > 0)].copy()
        matched["province"] = region
        frames.append(matched)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def predict_density(n3: np.ndarray, years: np.ndarray, a: float, b: float) -> np.ndarray:
    annual_rate = a + b * np.log(n3 / 500.0)
    return n3 * np.exp(years * annual_rate)


def metrics(df: pd.DataFrame, prediction: np.ndarray) -> dict:
    observed = df.n4.to_numpy()
    error = prediction - observed
    return {
        "n": len(df),
        "bias_stems_ha": round(float(error.mean()), 3),
        "mae_stems_ha": round(float(np.abs(error).mean()), 3),
        "rmse_stems_ha": round(float(np.sqrt(np.mean(error**2))), 3),
        "median_prediction_over_observed": round(float(np.median(prediction / observed)), 3),
    }


def evaluate() -> dict:
    df = load_remeasurements()
    if df.empty:
        raise RuntimeError("No eligible IFN3→IFN4 pine plots")
    fit = df[df.province.isin(FIT_PROVINCES)].copy()
    held = df[~df.province.isin(FIT_PROVINCES)].copy()
    if len(fit) < 20 or len(held) < 20:
        raise RuntimeError(f"Insufficient fit/test plots: {len(fit)}/{len(held)}")
    fit_rate = np.log(fit.n4.to_numpy() / fit.n3.to_numpy()) / fit.years.to_numpy()
    constant_rate = float(np.median(fit_rate))
    x = np.log(fit.n3.to_numpy() / 500.0)
    model = least_squares(lambda p: p[0] + p[1] * x - fit_rate,
                          x0=[constant_rate, -.01], loss="soft_l1", f_scale=.03)
    if not model.success:
        raise RuntimeError(model.message)
    predictions = {
        "persistence": held.n3.to_numpy(),
        "constant_annual_rate": predict_density(held.n3.to_numpy(), held.years.to_numpy(), constant_rate, 0.0),
        "density_dependent_rate": predict_density(held.n3.to_numpy(), held.years.to_numpy(), *model.x),
    }
    return {
        "design": "Same-plot IFN3 to IFN4 pine >=80% basal area both cycles, >=2 strict pine links; Madrid/Salamanca fit; other provinces test",
        "fit_n": len(fit), "heldout_n": len(held),
        "fit_by_province": fit.province.value_counts().to_dict(),
        "heldout_by_province": held.province.value_counts().to_dict(),
        "heldout_interval_years": [round(float(held.years.min()), 2), round(float(held.years.max()), 2)],
        "fit_median_annual_log_rate": constant_rate,
        "fit_density_rate_parameters": {"annual_log_rate_at_500_stems_ha": float(model.x[0]),
                                        "slope_per_log_density": float(model.x[1])},
        "heldout_results": {name: metrics(held, prediction) for name, prediction in predictions.items()},
        "heldout_province_results": {
            province: {name: metrics(g, prediction[held.province == province])
                       for name, prediction in predictions.items()}
            for province, g in held.groupby("province")
        },
        "heldout_observed_change": {
            "median_n3_stems_ha": round(float(held.n3.median()), 3),
            "median_n4_stems_ha": round(float(held.n4.median()), 3),
            "median_n4_over_n3": round(float((held.n4 / held.n3).median()), 3),
            "share_with_net_increase": round(float((held.n4 > held.n3).mean()), 3),
        },
    }


if __name__ == "__main__":
    import json
    print(json.dumps(evaluate(), indent=2, ensure_ascii=False))
