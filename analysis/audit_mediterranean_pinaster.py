"""Audit Mediterranean maritime-pine IFN4 stands before calibration.

Uses the same age and measured-stand filters as the Galician analysis.
Reports age-source and equation-domain attrition, with no parameter fitting.
Run with the sim_case1 environment from the repository root.
"""

from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "models"))

from ifn4_age_validation import PLANTING_LAG_YEARS, load_aged_plots, measured_stands  # noqa: E402
from mediterranean_pinaster import SOURCE_STAND_AGE_RANGE, site_index_70  # noqa: E402
from audit_pine_ifn4 import CENTRAL_DATABASES, PINASTER_CODE  # noqa: E402


def province_audit(region: str, database: str) -> dict:
    aged = load_aged_plots(region, database)
    aged = aged[(aged.Especie == PINASTER_CODE) &
                (aged.Edad - PLANTING_LAG_YEARS >= 10)].copy()
    stands = measured_stands(database, aged) if len(aged) else pd.DataFrame()
    if len(stands):
        stands = stands.merge(aged[["plot_id", "Edad", "Fiabil"]], on="plot_id", validate="one_to_one")
        stands = stands[(stands.h0_m > 1.3) & np.isfinite(stands.h0_m)].copy()
        stands["si70"] = [site_index_70(h, age) for h, age in zip(stands.h0_m, stands.Edad)]
    else:
        stands = pd.DataFrame(columns=["age_source", "planted", "Edad", "Fiabil", "h0_m", "si70"])
    source_counts = {str(int(k)): int(v) for k, v in stands.age_source.value_counts().items()}
    strict = stands[stands.age_source.isin((5, 8))]
    in_range = stands[stands.Edad.between(*SOURCE_STAND_AGE_RANGE)]
    return {
        "province": region,
        "screened_age_10plus": len(aged),
        "pine_dominated_measured_with_h0": len(stands),
        "recorded_planted": int(stands.planted.sum()),
        "age_source_counts": source_counts,
        "whorl_or_borer_age": len(strict),
        "whorl_or_borer_in_source_range": int(strict.Edad.between(*SOURCE_STAND_AGE_RANGE).sum()),
        "all_in_source_age_range": len(in_range),
        "median_total_age": float(stands.Edad.median()) if len(stands) else None,
        "median_h0_m": float(stands.h0_m.median()) if len(stands) else None,
        "median_si70_m": float(stands.si70.median()) if len(stands) else None,
        "si70_iqr_m": [float(stands.si70.quantile(.25)), float(stands.si70.quantile(.75))] if len(stands) else None,
    }


if __name__ == "__main__":
    import json
    results = [province_audit(region, database) for region, database in CENTRAL_DATABASES.items()]
    print(json.dumps(results, indent=2, ensure_ascii=False))
