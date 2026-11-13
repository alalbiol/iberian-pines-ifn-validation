"""Prospective IFN3→IFN4 pine tree-list biomass forecast on uncut plots.

Uses only initial tree d/h, density and interval at prediction time. The
components are Madrid/Salamanca-fitted DBH growth, constant height increment
and natural mortality, plus published ingrowth basal area converted to AGB
with fit-province data. Other provinces provide the test. This is a
15-year interval experiment, not a validated planting-to-harvest simulator.
"""

from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "models"))

from audit_mediterranean_turnover import load_turnover  # noqa: E402
from evaluate_mediterranean_expansion import expected_weight  # noqa: E402
from evaluate_mediterranean_pinaster import FIT_PROVINCES  # noqa: E402
from fit_mediterranean_density import PROVINCES  # noqa: E402
from fit_mediterranean_diameter import FIT_CODES, load_pairs  # noqa: E402
from ifn4_age_validation import IFN4_DIR, _plot_key  # noqa: E402
from ifn_remeasurement_validation import _mdb, _num, trees_per_ha  # noqa: E402
from mediterranean_pinaster_diameter import MEDITERRANEAN_PARAMS  # noqa: E402
from mediterranean_pinaster_ingrowth import SOURCE_BA_RANGE_M2_HA, SOURCE_QMD_RANGE_CM  # noqa: E402
from mediterranean_pinaster_tree_forecast import InventoryTree, project_inventory  # noqa: E402
from spain_biomass_equations import pinus_pinaster_biomass  # noqa: E402


def pine_biomass_kg(diameter_cm: np.ndarray, height_m: np.ndarray) -> np.ndarray:
    return np.array([pinus_pinaster_biomass(float(d), float(h)).aboveground_kg
                     for d, h in zip(diameter_cm, height_m)])


def ifn4_pine_agb(eligible_ids: set[str]) -> pd.DataFrame:
    frames = []
    for code, (_, database) in PROVINCES.items():
        ids = {pid for pid in eligible_ids if pid.startswith(f"{code}_")}
        if not ids:
            continue
        path = Path(IFN4_DIR) / f"{database}.accdb"
        trees = _mdb(str(path), "PCMayores")
        trees["plot_id"] = _plot_key(trees, code)
        trees = trees[trees.plot_id.isin(ids)].copy()
        trees["species"] = _num(trees.Especie)
        trees["ord4"] = _num(trees.OrdenIf4)
        trees["d_cm"] = (_num(trees.Dn1) + _num(trees.Dn2)) / 20.0
        trees["h_m"] = _num(trees.Ht)
        trees = trees[(trees.species == 26) & trees.d_cm.ge(7.5) & trees.h_m.gt(1.3) &
                      (trees.ord4 > 0) & ~trees.ord4.isin((444, 666, 777, 888, 999))].copy()
        trees["w"] = trees_per_ha(trees.d_cm.to_numpy())
        trees["agb_mg_ha"] = trees.w.to_numpy() * pine_biomass_kg(
            trees.d_cm.to_numpy(), trees.h_m.to_numpy()) / 1000.0
        frames.append(trees.groupby("plot_id").agb_mg_ha.sum())
    return pd.concat(frames).rename("agb4_mg_ha").reset_index()


def metrics(df: pd.DataFrame, predicted: str, observed: str) -> dict:
    error = df[predicted] - df[observed]
    return {"n": len(df), "bias_mg_ha": round(float(error.mean()), 3),
            "mae_mg_ha": round(float(error.abs().mean()), 3),
            "rmse_mg_ha": round(float(np.sqrt(np.mean(error**2))), 3),
            "median_model_over_observed": round(float((df[predicted] / df[observed]).median()), 3)}


def evaluate() -> dict:
    plots, trees = load_turnover()
    plots = plots[(plots.cut == 0) & (plots.other == 0) & (plots.unrecorded == 0)].copy()
    trees = trees[trees.plot_id.isin(set(plots.plot_id))].merge(
        plots[["plot_id", "province", "n3", "years"]], on="plot_id", validate="many_to_one")
    trees["lost"] = trees.category.isin(("dead", "not_found")).astype(float)
    fit_trees = trees[trees.province.isin(FIT_PROVINCES)]
    fit_weight = fit_trees.w3.to_numpy() / fit_trees.n3.to_numpy()
    fit_target = fit_trees.lost.to_numpy()

    def survival_loss(log_rate: float) -> float:
        p = np.clip(-np.expm1(-fit_trees.years.to_numpy() * np.exp(log_rate)), 1e-9, 1-1e-9)
        return float(-np.sum(fit_weight * (fit_target*np.log(p)+(1-fit_target)*np.log1p(-p))))

    survival_fit = minimize_scalar(survival_loss, bounds=(-8, -1), method="bounded")
    if not survival_fit.success:
        raise RuntimeError("Survival fit failed")
    annual_hazard = float(np.exp(survival_fit.x))
    fit_pairs = load_pairs()
    fit_pairs = fit_pairs[fit_pairs.province_code.isin(FIT_CODES)]
    flat_height_rate = float(np.median(
        fit_pairs.assign(annual_height=(fit_pairs.h4_m-fit_pairs.h3_m)/fit_pairs.years)
        .groupby("plot_id").annual_height.median()))
    p = MEDITERRANEAN_PARAMS
    rate_fit = np.exp(p.log_rate_at_20cm_500ha
                      + p.log_diameter_exponent*np.log(fit_pairs.d3_cm.to_numpy()/20.0)
                      + p.log_density_exponent*np.log(fit_pairs.density_per_ha.to_numpy()/500.0))
    residual = fit_pairs.d4_cm.to_numpy() - (
        fit_pairs.d3_cm.to_numpy() + fit_pairs.years.to_numpy()*rate_fit)
    residual_mean, residual_sd = float(residual.mean()), float(residual.std())

    d3, h3 = trees.d3_cm.to_numpy(), trees.h3_m.to_numpy()
    years, n3 = trees.years.to_numpy(), trees.n3.to_numpy()
    dbh_rate = np.exp(p.log_rate_at_20cm_500ha
                      + p.log_diameter_exponent*np.log(d3/20.0)
                      + p.log_density_exponent*np.log(n3/500.0))
    d4_pred = d3 + years*dbh_rate
    h4_pred = h3 + years*flat_height_rate
    w4_pred = expected_weight(d3, d4_pred, residual_mean, residual_sd)
    survival_probability = np.exp(-annual_hazard*years)
    trees["agb3"] = trees.w3.to_numpy()*pine_biomass_kg(d3,h3)/1000.0
    trees["agb4_cohort_pred_no_loss"] = w4_pred*pine_biomass_kg(d4_pred,h4_pred)/1000.0
    trees["agb4_cohort_pred"] = survival_probability*trees.agb4_cohort_pred_no_loss.to_numpy()
    trees["agb4_cohort_pred_oracle_survival"] = np.where(
        trees.category == "survived", trees.agb4_cohort_pred_no_loss, 0.0)
    current_live = ((trees.category == "survived") & (trees.w4 > 0) &
                    trees.d4_cm.ge(7.5) & trees.h4_m.gt(1.3))
    d4_safe = np.where(current_live, trees.d4_cm.to_numpy(), d3)
    h4_safe = np.where(current_live, trees.h4_m.to_numpy(), h3)
    trees["agb4_cohort_observed"] = np.where(
        current_live, trees.w4.to_numpy()*pine_biomass_kg(d4_safe, h4_safe)/1000.0,
        0.0)
    sums = trees.groupby("plot_id")[["agb3", "agb4_cohort_pred_no_loss", "agb4_cohort_pred",
                                      "agb4_cohort_pred_oracle_survival", "agb4_cohort_observed"]].sum()
    plots = plots.merge(sums, on="plot_id", validate="one_to_one")
    plots = plots.merge(ifn4_pine_agb(set(plots.plot_id)), on="plot_id", validate="one_to_one")
    plots["cohort_plus_observed_entry"] = plots.agb4_cohort_pred + plots.new_agb_mg_ha
    source_domain = plots.ba3.between(*SOURCE_BA_RANGE_M2_HA) & plots.qmd_cm3.between(*SOURCE_QMD_RANGE_CM)
    source = plots[source_domain].copy()
    fit_source = source[source.province.isin(FIT_PROVINCES)]
    held_source = source[~source.province.isin(FIT_PROVINCES)].copy()
    if len(fit_source) < 20 or len(held_source) < 20:
        raise RuntimeError("Insufficient uncut source-domain biomass plots")
    fit_agb_per_new_ba = float(fit_source.new_agb_mg_ha.sum()/fit_source.new_ba_m2_ha.sum())
    fit_annual_new_agb = float((fit_source.new_agb_mg_ha/fit_source.years).mean())
    # Exercise the reusable prospective API on the untouched provinces.
    projections = {}
    held_source_ids = set(held_source.plot_id)
    for plot_id, group in trees.groupby("plot_id"):
        if plot_id not in held_source_ids:
            continue
        row = held_source.loc[held_source.plot_id == plot_id].iloc[0]
        initial = [InventoryTree(float(t.d3_cm), float(t.h3_m)) for t in group.itertuples()]
        projections[plot_id] = project_inventory(
            initial, float(row.years), annual_natural_hazard=annual_hazard,
            annual_height_increment_m=flat_height_rate,
            dbh_residual_mean_cm=residual_mean, dbh_residual_sd_cm=residual_sd,
            basal_area_m2_ha=float(row.ba3), quadratic_mean_diameter_cm=float(row.qmd_cm3),
            new_agb_per_new_ba_mg_m2=fit_agb_per_new_ba)
    held_source["with_published_ingrowth"] = held_source.plot_id.map(
        lambda pid: projections[pid].total_agb_mg_ha)
    parity_error = max(abs(projections[row.plot_id].surviving_cohort_agb_mg_ha - row.agb4_cohort_pred)
                       for row in held_source.itertuples())
    if parity_error > 1e-8:
        raise RuntimeError(f"Reusable tree-list forecast differs from evaluation: {parity_error}")
    held_source["with_fit_constant_ingrowth"] = (held_source.agb4_cohort_pred
        + fit_annual_new_agb*held_source.years.to_numpy())
    held = plots[~plots.province.isin(FIT_PROVINCES)]
    province_results = {province: {
        "cohort": metrics(group, "agb4_cohort_pred", "agb4_cohort_observed"),
        "total": metrics(held_source[held_source.province == province],
                         "with_published_ingrowth", "agb4_mg_ha")
        if (held_source.province == province).any() else None,
    } for province, group in held.groupby("province")}
    # Paired plot bootstrap: uncertainty in the RMSE difference, not a new fit.
    rng = np.random.default_rng(20261022)
    truth = held_source.agb4_mg_ha.to_numpy()
    persistence = held_source.agb3.to_numpy() - truth
    cohort_only = held_source.agb4_cohort_pred.to_numpy() - truth
    forecast = held_source.with_published_ingrowth.to_numpy() - truth
    samples = rng.integers(0, len(truth), size=(2000, len(truth)))
    rmse_gain = np.sqrt(np.mean(persistence[samples]**2, axis=1)) - np.sqrt(
        np.mean(forecast[samples]**2, axis=1))
    ingrowth_gain = np.sqrt(np.mean(cohort_only[samples]**2, axis=1)) - np.sqrt(
        np.mean(forecast[samples]**2, axis=1))
    return {
        "design": "Uncut IFN3 pine tree-list to IFN4 pine AGB over 15-18 years; Madrid/Salamanca component fit, other provinces test; IFN4 pine trees >=7.5cm reference allometry",
        "fit_plots": int(plots.province.isin(FIT_PROVINCES).sum()),
        "heldout_plots": len(held), "heldout_source_domain_plots": len(held_source),
        "fit_annual_natural_hazard": annual_hazard,
        "fit_flat_height_increment_m_per_year": flat_height_rate,
        "fit_new_agb_per_new_ba_mg_m2": fit_agb_per_new_ba,
        "fit_diameter_residual_mean_cm": residual_mean,
        "fit_diameter_residual_sd_cm": residual_sd,
        "heldout_cohort_results": {
            name: metrics(held, name, "agb4_cohort_observed") for name in
            ("agb3", "agb4_cohort_pred_no_loss", "agb4_cohort_pred",
             "agb4_cohort_pred_oracle_survival")},
        "heldout_total_results": {
            name: metrics(held_source, name, "agb4_mg_ha") for name in
            ("agb3", "agb4_cohort_pred", "cohort_plus_observed_entry",
             "with_published_ingrowth", "with_fit_constant_ingrowth")},
        "heldout_accounting": metrics(held_source.assign(accounted=lambda x: x.agb4_cohort_observed+x.new_agb_mg_ha),
                                      "accounted", "agb4_mg_ha"),
        "heldout_by_province": province_results,
        "paired_plot_bootstrap_rmse_gain_over_initial_mg_ha_95pct":
            [round(float(x), 3) for x in np.quantile(rmse_gain, [.025, .975])],
        "paired_plot_bootstrap_rmse_gain_from_ingrowth_mg_ha_95pct":
            [round(float(x), 3) for x in np.quantile(ingrowth_gain, [.025, .975])],
    }


if __name__ == "__main__":
    import json
    print(json.dumps(evaluate(), indent=2, ensure_ascii=False))
