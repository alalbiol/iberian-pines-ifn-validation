"""Export plot-level forecasts behind the manuscript tables and figures.

Reruns the existing evaluations unchanged and captures the held-out plot
frames they score, so the figures and the extra reference rows use exactly
the plots and forecasts already reported in results/*.json. Writes CSV files
to results/plot_level/. Run from the repository root with `sim_case1`.
"""

from pathlib import Path
import sys

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "models"))
sys.path.insert(0, str(HERE))
OUT = HERE / "results" / "plot_level"
OUT.mkdir(parents=True, exist_ok=True)


def measured_tree_other_pines() -> pd.DataFrame:
    import evaluate_other_pine_inventory_projection as mod
    original, frames = mod.metrics, []

    def capture(frame, predicted, observed):
        if predicted == "total_pred" and observed == "agb4":
            frames.append(frame.copy())
        return original(frame, predicted, observed)

    mod.metrics = capture
    parts = []
    for code, name in mod.SPECIES.items():
        frames.clear()
        mod.evaluate_species(code, name)
        held = frames[-1]
        parts.append(pd.DataFrame({
            "case": name, "plot_id": held.plot_id, "province_code": held.province_code,
            "years": held.years, "agb3": held.initial_agb, "pred_no_entry": held.old_pred,
            "pred_total": held.total_pred, "agb4": held.agb4}))
    mod.metrics = original
    return pd.concat(parts, ignore_index=True)


def measured_tree_mediterranean() -> pd.DataFrame:
    import evaluate_mediterranean_tree_biomass as mod
    from fit_mediterranean_diameter import PROVINCES
    original, frames = mod.metrics, []

    def capture(frame, predicted, observed):
        if predicted == "with_published_ingrowth" and observed == "agb4_mg_ha":
            frames.append(frame.copy())
        return original(frame, predicted, observed)

    mod.metrics = capture
    mod.evaluate()
    mod.metrics = original
    held = max(frames, key=len)
    code = {name: c for c, name in PROVINCES.items()}
    return pd.DataFrame({
        "case": "Pinus pinaster (Mediterranean)", "plot_id": held.plot_id,
        "province_code": held.province.map(code), "years": held.years, "agb3": held.agb3,
        "pred_no_entry": held.agb4_cohort_pred, "pred_total": held.with_published_ingrowth,
        "agb4": held.agb4_mg_ha})


def age_only_other_pines() -> pd.DataFrame:
    import evaluate_other_pine_stands as mod
    original, frames = mod.metrics, {}

    def capture(frame, predicted, observed):
        if predicted == "default_agb":
            frames[frame.species.iloc[0]] = frame.copy()
        return original(frame, predicted, observed)

    mod.metrics = capture
    mod.evaluate()
    mod.metrics = original
    return pd.concat([pd.DataFrame({
        "case": name, "plot_id": f.plot_id, "province_code": f.province_code, "age": f.Edad,
        "observed": f.agb_mg_ha, "default": f.default_agb, "measured_inputs": f.measured_hn_agb,
        "reference": f.fit_median_agb}) for name, f in frames.items()], ignore_index=True)


def age_only_galicia() -> pd.DataFrame:
    from ifn4_age_validation import (DATABASES, TEST_REGIONS, large_tree_view, load_aged_plots,
                                     measured_stands, run_models)
    from ifn4_pine_calibration import FIT_REGIONS

    def stands(regions):
        s = large_tree_view(pd.concat([measured_stands(DATABASES[r], load_aged_plots(r, DATABASES[r]))
                                       for r in regions], ignore_index=True))
        return s[(s.species == "pinus_pinaster") & (s.age >= 10) & (s.agb_mg_ha > 1)].copy()

    reference = float(stands(FIT_REGIONS).agb_mg_ha.median())
    p = run_models(stands(TEST_REGIONS))
    return pd.DataFrame({
        "case": "Pinus pinaster (Atlantic)", "plot_id": p.plot_id, "province_code": p.region,
        "age": p.age, "observed": p.agb_mg_ha, "default": p.modelA_agb,
        "measured_inputs": p.modelB_agb, "reference": reference,
        "direct_age": p.age_source.isin((5, 8))})


def age_only_mediterranean() -> pd.DataFrame:
    import evaluate_mediterranean_ifn3 as mod
    from evaluate_mediterranean_pinaster import FIT_PROVINCES, load_eligible
    captured = {}
    original = mod.metrics

    def capture(frame, name):
        if name == "defaults":
            captured.setdefault("frames", []).append(frame.copy())
        return original(frame, name)

    mod.metrics = capture
    mod.evaluate()
    mod.metrics = original
    fresh = max(captured["frames"], key=len)
    from ifn4_age_validation import large_tree_view
    eligible = large_tree_view(load_eligible())
    reference = float(eligible[eligible.region.isin(FIT_PROVINCES)].agb_mg_ha.median())
    return pd.DataFrame({
        "case": "Pinus pinaster (Mediterranean)", "plot_id": fresh.plot_id,
        "province_code": fresh.region, "age": fresh.Edad, "observed": fresh.agb_mg_ha,
        "default": fresh.defaults, "measured_inputs": fresh.measured_height_and_density,
        "reference": reference})


def young_mediterranean_reference() -> dict:
    import evaluate_young_afforestation as mod
    from evaluate_mediterranean_pinaster import FIT_PROVINCES
    df = mod.load_eligible()
    fit, held = df[df.region.isin(FIT_PROVINCES)], df[~df.region.isin(FIT_PROVINCES)]
    reference = float(fit.agb_mg_ha.median())
    error = reference - held.agb_mg_ha
    return {"fit_plots": len(fit), "heldout_plots": len(held), "fit_median_agb_mg_ha": reference,
            "reference_rmse_mg_ha": float(np.sqrt(np.mean(error**2)))}


if __name__ == "__main__":
    import json
    measured = pd.concat([measured_tree_mediterranean(), measured_tree_other_pines()],
                         ignore_index=True)
    measured.to_csv(OUT / "measured_tree_forecasts.csv", index=False)
    age = pd.concat([age_only_galicia(), age_only_mediterranean(), age_only_other_pines()],
                    ignore_index=True)
    age.to_csv(OUT / "age_only_forecasts.csv", index=False)
    young = young_mediterranean_reference()
    (OUT / "young_mediterranean_reference.json").write_text(json.dumps(young, indent=2) + "\n")
    print(measured.groupby("case").size(), age.groupby("case").size(), young, sep="\n")
