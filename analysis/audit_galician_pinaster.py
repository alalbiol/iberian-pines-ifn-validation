"""Reproduce and audit the held-out Galician Pinus pinaster IFN4 stand result.

The age shifts below are sensitivity scenarios, not inferred corrections to
individual IFN age records. Stem number and observed AGB refer to trees >= 7.5 cm
only (large_tree_view), as in the other pine cases and the stand-chain calibration.
Run with the sim_case1 Python environment.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "models"))

from ifn4_age_validation import (  # noqa: E402
    DATABASES, IFN4_DIR, TEST_REGIONS, _plot_key, load_aged_plots,
    large_tree_view, measured_stands, run_models,
)
from ifn_remeasurement_validation import _mdb, _num  # noqa: E402

OUT = Path(__file__).resolve().parent / "results" / "galician_pinaster_audit.json"


def ratio_summary(frame: pd.DataFrame, field: str) -> dict:
    values = frame[field] / frame.agb_mg_ha
    return {
        "plots": int(len(values)),
        "median": round(float(values.median()), 3),
        "q25": round(float(values.quantile(0.25)), 3),
        "q75": round(float(values.quantile(0.75)), 3),
    }


def evaluate() -> dict:
    frames = []
    for region in TEST_REGIONS:
        db = DATABASES[region]
        aged = load_aged_plots(region, db)
        frames.append(measured_stands(db, aged))
    stands = large_tree_view(pd.concat(frames, ignore_index=True))
    stands = stands[(stands.species == "pinus_pinaster") &
                    (stands.age >= 10) & (stands.agb_mg_ha > 1)].copy()
    predicted = run_models(stands)
    direct = predicted[predicted.age_source.isin((5, 8))]
    historical_rows = 0
    historical_plots: set[str] = set()
    for region in TEST_REGIONS:
        db = DATABASES[region]
        keys = set(predicted.loc[predicted.region == region, "plot_id"])
        trees = _mdb(str(Path(IFN4_DIR) / f"{db}.accdb"), "PCMayores")
        trees["plot_id"] = _plot_key(trees)
        trees = trees[trees.plot_id.isin(keys)].copy()
        trees["d_cm"] = (_num(trees.Dn1) + _num(trees.Dn2)) / 20.0
        trees["h_m"] = _num(trees.Ht)
        trees = trees[(trees.d_cm >= 7.5) & (trees.h_m > 1.3)]
        order = _num(trees.OrdenIf4)
        bad = trees[(order.isna()) | (order <= 0) | order.isin((444, 888))]
        historical_rows += len(bad)
        historical_plots.update(bad.plot_id)

    result = {
        "design": "Ourense and Pontevedra held out; trees >= 7.5 cm only for stem number and AGB; direct age=IFN4 FInfor 5 or 8; source 7 retained only with Fiabil 3/4; age shifts are scenarios",
        "heldout_plots": len(predicted),
        "by_province": {str(k): int(v) for k, v in predicted.region.value_counts().items()},
        "age_source_counts": {str(int(k)): int(v) for k, v in predicted.age_source.value_counts().items()},
        "planted_plots": int(predicted.planted.sum()),
        "eligible_historical_tree_rows": historical_rows,
        "plots_with_eligible_historical_rows": len(historical_plots),
        "all_age_sources": {
            "default_agb_ratio": ratio_summary(predicted, "modelA_agb"),
            "measured_height_density_agb_ratio": ratio_summary(predicted, "modelB_agb"),
        },
        "direct_age_only": {
            "default_agb_ratio": ratio_summary(direct, "modelA_agb"),
            "measured_height_density_agb_ratio": ratio_summary(direct, "modelB_agb"),
        },
        "age_shift_sensitivity": {},
    }
    for shift in (-2, -1, 0, 1, 2):
        shifted = stands.copy()
        shifted["age"] += shift
        scored = run_models(shifted)
        result["age_shift_sensitivity"][str(shift)] = {
            "default_agb_ratio": ratio_summary(scored, "modelA_agb"),
            "measured_height_density_agb_ratio": ratio_summary(scored, "modelB_agb"),
        }
    return result


if __name__ == "__main__":
    result = evaluate()
    OUT.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
