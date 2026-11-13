"""Forecast the concentric-plot expansion change of surviving pine trees.

Trees crossing a DBH threshold move to a larger IFN subplot and receive a
smaller stems/ha expansion weight. This analysis conditions on survivors
identified at IFN4 to isolate the weight-change term; it is not a full
prospective density forecast because future survivors are unknown.
"""

from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "use_case_1_simulator/agb_forecast_analysis"))

from audit_mediterranean_turnover import load_turnover  # noqa: E402
from evaluate_mediterranean_pinaster import FIT_PROVINCES  # noqa: E402
from fit_mediterranean_diameter import FIT_CODES, load_pairs  # noqa: E402
from ifn_remeasurement_validation import trees_per_ha  # noqa: E402
from mediterranean_pinaster_diameter import MEDITERRANEAN_PARAMS  # noqa: E402


def metrics(df: pd.DataFrame, prediction: np.ndarray) -> dict:
    error = prediction - df.survivor_expansion_delta.to_numpy()
    return {"n": len(df), "bias_stems_ha": round(float(error.mean()), 3),
            "mae_stems_ha": round(float(np.abs(error).mean()), 3),
            "rmse_stems_ha": round(float(np.sqrt(np.mean(error**2))), 3)}


def expected_weight(d3: np.ndarray, predicted_d4: np.ndarray,
                    residual_mean: float, residual_sd: float) -> np.ndarray:
    """Expected IFN expansion weight with fitted DBH error at circle thresholds."""
    weight = np.full(len(d3), float(trees_per_ha(np.array([7.5]))[0]))
    for threshold, next_d in ((12.5, 12.5), (22.5, 22.5), (42.5, 42.5)):
        before = float(trees_per_ha(np.array([threshold - .01]))[0])
        after = float(trees_per_ha(np.array([next_d]))[0])
        probability = np.where(d3 >= threshold, 1.0,
                               norm.sf((threshold - predicted_d4 - residual_mean) / residual_sd))
        weight += (after - before) * probability
    return weight


def evaluate() -> dict:
    plots, trees = load_turnover()
    survivors = trees[trees.category == "survived"].merge(
        plots[["plot_id", "n3", "years"]], on="plot_id", validate="many_to_one")
    params = MEDITERRANEAN_PARAMS
    rate = np.exp(params.log_rate_at_20cm_500ha
                  + params.log_diameter_exponent * np.log(survivors.d3_cm.to_numpy() / 20.0)
                  + params.log_density_exponent * np.log(survivors.n3.to_numpy() / 500.0))
    fit_pairs = load_pairs()
    fit_pairs = fit_pairs[fit_pairs.province_code.isin(FIT_CODES)]
    flat_rate = float(np.median(fit_pairs.groupby("plot_id").annual_observed_cm.median()))
    fit_pred = fit_pairs.d3_cm.to_numpy() + fit_pairs.years.to_numpy() * np.exp(
        params.log_rate_at_20cm_500ha
        + params.log_diameter_exponent * np.log(fit_pairs.d3_cm.to_numpy() / 20.0)
        + params.log_density_exponent * np.log(fit_pairs.density_per_ha.to_numpy() / 500.0))
    residual = fit_pairs.d4_cm.to_numpy() - fit_pred
    residual_mean, residual_sd = float(residual.mean()), float(residual.std())
    predicted_d4 = survivors.d3_cm.to_numpy() + survivors.years.to_numpy() * rate
    survivors["model_shift"] = (trees_per_ha(predicted_d4)
                                - survivors.w3.to_numpy())
    survivors["expected_model_shift"] = (expected_weight(survivors.d3_cm.to_numpy(), predicted_d4,
                                                          residual_mean, residual_sd)
                                         - survivors.w3.to_numpy())
    survivors["flat_shift"] = (trees_per_ha(survivors.d3_cm.to_numpy() + survivors.years.to_numpy() * flat_rate)
                               - survivors.w3.to_numpy())
    shifts = survivors.groupby("plot_id").agg(model_shift=("model_shift", "sum"),
                                               expected_model_shift=("expected_model_shift", "sum"),
                                               flat_shift=("flat_shift", "sum"))
    plots = plots.merge(shifts, on="plot_id", validate="one_to_one")
    held = plots[~plots.province.isin(FIT_PROVINCES)]
    return {
        "design": "IFN3→IFN4 pine-dominated plots; observed IFN4 survivor set supplied; Madrid/Salamanca DBH model, other provinces test",
        "survivor_trees": len(survivors), "heldout_n": len(held),
        "heldout_by_province": held.province.value_counts().to_dict(),
        "fit_flat_diameter_rate_cm_per_year": round(flat_rate, 4),
        "fit_diameter_residual_mean_cm": round(residual_mean, 3),
        "fit_diameter_residual_sd_cm": round(residual_sd, 3),
        "heldout_results": {
            "error_integrated_diameter_model": metrics(held, held.expected_model_shift.to_numpy()),
            "size_density_diameter_model": metrics(held, held.model_shift.to_numpy()),
            "flat_diameter_rate": metrics(held, held.flat_shift.to_numpy()),
            "no_expansion_shift": metrics(held, np.zeros(len(held))),
        },
        "heldout_mean_observed_shift_stems_ha": round(float(held.survivor_expansion_delta.mean()), 3),
    }


if __name__ == "__main__":
    import json
    print(json.dumps(evaluate(), indent=2, ensure_ascii=False))
