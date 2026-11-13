"""First independent IFN4 height test for Mediterranean P. pinaster.

Calibration: Madrid + Salamanca only. Held out: other central provinces.
Only pine-dominated stands with measured H0, FInfor 5/8 and total age
within the published 2007 fitting range are used. A single default site
index is calibrated on fit provinces; published growth parameters are not
refit. This is cross-sectional height validation, not longitudinal growth
validation. Run from the repository root with the sim_case1 environment.
"""

from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "use_case_1_simulator/agb_forecast_analysis"))

from audit_pine_ifn4 import CENTRAL_DATABASES, PINASTER_CODE  # noqa: E402
from ifn4_age_validation import load_aged_plots, measured_stands  # noqa: E402
from mediterranean_pinaster import (  # noqa: E402
    REFERENCE_AGE, SOURCE_STAND_AGE_RANGE, project_height, site_index_70,
)
from spain_growth_curves import pinus_pinaster_project  # noqa: E402

FIT_PROVINCES = ("Madrid", "Salamanca")


def load_eligible() -> pd.DataFrame:
    frames = []
    for province, database in CENTRAL_DATABASES.items():
        aged = load_aged_plots(province, database)
        aged = aged[(aged.Especie == PINASTER_CODE) & aged.FInfor.isin((5, 8)) &
                    aged.Edad.between(*SOURCE_STAND_AGE_RANGE)].copy()
        if aged.empty:
            continue
        measured = measured_stands(database, aged)
        if measured.empty:
            continue
        measured = measured.merge(aged[["plot_id", "Edad", "Fiabil"]], on="plot_id", validate="one_to_one")
        measured = measured[(measured.h0_m > 1.3) & np.isfinite(measured.h0_m)].copy()
        frames.append(measured)
    return pd.concat(frames, ignore_index=True)


def metrics(df: pd.DataFrame, column: str) -> dict:
    error = df[column] - df.h0_m
    return {
        "n": len(df),
        "bias_m": round(float(error.mean()), 3),
        "mae_m": round(float(error.abs().mean()), 3),
        "rmse_m": round(float(np.sqrt(np.mean(error**2))), 3),
        "median_predicted_over_observed": round(float((df[column] / df.h0_m).median()), 3),
    }


def evaluate() -> dict:
    df = load_eligible()
    fit = df[df.region.isin(FIT_PROVINCES)].copy()
    held = df[~df.region.isin(FIT_PROVINCES)].copy()
    if len(fit) < 10 or len(held) < 10:
        raise RuntimeError(f"Insufficient calibration/test stands: {len(fit)}/{len(held)}")

    # Match each model's default site to fit provinces only. The Atlantic
    # Galician interior curve is shown solely as a regional-transfer baseline.
    med_si = float(np.median([site_index_70(h, age) for h, age in zip(fit.h0_m, fit.Edad)]))
    gal_si = float(np.median([pinus_pinaster_project(h, age, 20.0, False)
                              for h, age in zip(fit.h0_m, fit.Edad)]))
    held["med_h0"] = [project_height(med_si, REFERENCE_AGE, age) for age in held.Edad]
    held["gal_h0"] = [pinus_pinaster_project(gal_si, 20.0, age, False) for age in held.Edad]
    held["constant_h0"] = float(fit.h0_m.median())
    return {
        "design": "Madrid and Salamanca fit; all other eligible provinces held out; strict age source 5/8; total age 25-180",
        "fit_n": len(fit),
        "heldout_n": len(held),
        "fit_mediterranean_si70_m": round(med_si, 3),
        "fit_galician_interior_si20_m": round(gal_si, 3),
        "heldout_by_province": held.region.value_counts().to_dict(),
        "heldout_total_age_range": [float(held.Edad.min()), float(held.Edad.max())],
        "heldout_results": {
            "mediterranean_2007_pooled": metrics(held, "med_h0"),
            "galicia_2005_interior_transfer": metrics(held, "gal_h0"),
            "fit_median_height_no_age": metrics(held, "constant_h0"),
        },
        "heldout_province_results": {
            province: {"mediterranean": metrics(g, "med_h0"), "galicia_transfer": metrics(g, "gal_h0")}
            for province, g in held.groupby("region")
        },
    }


if __name__ == "__main__":
    import json
    print(json.dumps(evaluate(), indent=2, ensure_ascii=False))
