"""External first-pass test of Calama et al. (2024) young pine biomass.

IFN4 stands are screened for documented planting, direct age source, the
source's 4-41-year age range and 323-2281 stems/ha density range. The source
uses plantation mean height and years since planting. Its site index is
estimated from IFN mean height, or set to the fit-province median to emulate
an unmeasured site. IFN allometry differs from the source's young-tree
allometry, so this is an operational transfer comparison, not an isolated
validation of the source biomass equation.
"""

from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "models"))

from audit_pine_ifn4 import CENTRAL_DATABASES, PINASTER_CODE  # noqa: E402
from evaluate_mediterranean_pinaster import FIT_PROVINCES  # noqa: E402
from ifn4_age_validation import load_aged_plots, measured_stands  # noqa: E402
from young_pinaster_afforestation import (  # noqa: E402
    BIOMASS_AGE_RANGE, BIOMASS_DENSITY_RANGE, SITE_REFERENCE_AGE,
    aboveground_biomass_kg_ha, site_index_30,
)


def load_eligible() -> pd.DataFrame:
    frames = []
    for province, database in CENTRAL_DATABASES.items():
        aged = load_aged_plots(province, database)
        aged = aged[(aged.Especie == PINASTER_CODE) & aged.FInfor.isin((5, 8))].copy()
        if aged.empty:
            continue
        measured = measured_stands(database, aged)
        if measured.empty:
            continue
        measured = measured[(measured.planted) &
                            measured.age.between(*BIOMASS_AGE_RANGE) &
                            measured.trees_per_ha.between(*BIOMASS_DENSITY_RANGE) &
                            (measured.mean_h_m > 1.3) &
                            np.isfinite(measured.mean_h_m) &
                            (measured.agb_mg_ha > 0)].copy()
        if len(measured):
            frames.append(measured)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def metrics(df: pd.DataFrame, column: str) -> dict:
    error = df[column] - df.agb_mg_ha
    ratio = df[column] / df.agb_mg_ha
    return {
        "n": len(df),
        "bias_mg_ha": round(float(error.mean()), 3),
        "mae_mg_ha": round(float(error.abs().mean()), 3),
        "rmse_mg_ha": round(float(np.sqrt(np.mean(error**2))), 3),
        "median_model_over_ifn": round(float(ratio.median()), 3),
        "ratio_iqr": [round(float(ratio.quantile(.25)), 3), round(float(ratio.quantile(.75)), 3)],
    }


def evaluate() -> dict:
    df = load_eligible()
    if df.empty:
        raise RuntimeError("No eligible documented planted stands")
    df["si30"] = [site_index_30(h, age) for h, age in zip(df.mean_h_m, df.age)]
    fit = df[df.region.isin(FIT_PROVINCES)].copy()
    held = df[~df.region.isin(FIT_PROVINCES)].copy()
    if fit.empty or held.empty:
        raise RuntimeError(f"Need both fit and held-out planted stands, got {len(fit)}/{len(held)}")
    fit_si = float(fit.si30.median())
    for frame in (df, held):
        frame["published_with_measured_site"] = [
            aboveground_biomass_kg_ha(age, n, si) / 1000.0
            for age, n, si in zip(frame.age, frame.trees_per_ha, frame.si30)
        ]
        frame["published_with_fit_site"] = [
            aboveground_biomass_kg_ha(age, n, fit_si) / 1000.0
            for age, n in zip(frame.age, frame.trees_per_ha)
        ]
    return {
        "design": "Documented planting, direct age, 4-41 years, 323-2281 trees/ha; Madrid+Salamanca median SI30 fit; other provinces held out",
        "eligible_n": len(df), "fit_n": len(fit), "heldout_n": len(held),
        "fit_si30_mean_height_m": round(fit_si, 3),
        "heldout_by_province": held.region.value_counts().to_dict(),
        "heldout_age_range": [float(held.age.min()), float(held.age.max())],
        "heldout_results": {
            "published_with_fit_site": metrics(held, "published_with_fit_site"),
            "published_with_measured_site": metrics(held, "published_with_measured_site"),
        },
    }


if __name__ == "__main__":
    import json
    print(json.dumps(evaluate(), indent=2, ensure_ascii=False))
