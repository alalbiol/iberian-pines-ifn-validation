"""Prospective IFN3→IFN4 AGB projection without recorded target-pine cuts.

The input is the complete IFN3 measured tree list for each pine-dominated
plot, not future survivor identities. Fit constant natural-death hazard and
mean annual new-tree AGB on designated provinces only. Other provinces test
the tree-growth, mortality, sampling-weight and entry chain. Future cutting
is excluded here; published species-specific ingrowth remains to be tested.
"""

import os  # Configurable external inventory paths.

from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "models"))

from audit_other_pines import normalized_plot_key  # noqa: E402
from fit_other_pine_tree_growth import FIT_CODES, SPECIES  # noqa: E402
from ifn_remeasurement_validation import _mdb, _num, trees_per_ha  # noqa: E402
from other_pine_inventory_forecast import (InventoryTree, PARAMETERS,
                                           project_inventory)  # noqa: E402
from other_pine_tree_growth import BIOMASS  # noqa: E402

IFN3_DIR = Path((os.environ.get("IFN_DATA_ROOT", os.path.expanduser("~/ifn-data")) + "/ifn3"))
IFN4_DIR = Path((os.environ.get("IFN_DATA_ROOT", os.path.expanduser("~/ifn-data")) + "/ifn4"))
PAIRS_DIR = Path((os.environ.get("IFN_DATA_ROOT", os.path.expanduser("~/ifn-data")) + "/ifn_validation"))
IFN4_FILE = {
    "02": "Ifn4_Albacete", "05": "Ifn4p05-avila", "06": "Ifn4_Extremadura",
    "08": "Ifn4_Cataluna", "10": "Ifn4_Extremadura", "11": "Ifn4p11",
    "13": "Ifn4_CiudadReal", "14": "Ifn4p14", "15": "Ifn4p15-coruna",
    "16": "Ifn4p_Cuenca", "17": "Ifn4_Cataluna", "19": "Ifn4_Guadalajara",
    "23": "ifn4p23", "25": "Ifn4_Cataluna", "27": "Ifn4p27-lugo",
    "28": "Ifn4p28-madrid", "29": "Ifn4p29", "32": "Ifn4p32-ourense",
    "36": "Ifn4p36-pontevedra", "37": "Ifn4p37-salamanca",
    "43": "Ifn4_Cataluna", "45": "Ifn4_Toledo", "49": "Ifn4p49-zamora",
}
BAD_ORDER = (444, 666, 777, 888, 999)


def trees_from_db(path: Path, code: str, cycle: int) -> pd.DataFrame:
    t = _mdb(str(path), "PCMayores")
    t["plot_id"] = normalized_plot_key(t, int(code))
    t = t[t.plot_id.str.startswith(f"{int(code)}_")].copy()
    t["species"] = _num(t.Especie)
    t["d_cm"] = (_num(t.Dn1)+_num(t.Dn2))/20.0
    t["h_m"] = _num(t.Ht)
    t["ord3"] = _num(t.OrdenIf3)
    if cycle == 4:
        t["ord4"] = _num(t.OrdenIf4)
        order = t.ord4
    else:
        order = t.ord3
    t["live"] = (t.d_cm >= 7.5) & (t.h_m > 1.3) & (order > 0) & ~order.isin(BAD_ORDER)
    t["w"] = np.where(t.live, trees_per_ha(t.d_cm.to_numpy()), 0.0)
    t["ba"] = t.w*np.pi*(t.d_cm/200.0)**2
    return t


def plot_stocks(t: pd.DataFrame, species_code: int, suffix: str) -> pd.DataFrame:
    live = t[t.live]
    total_ba = live.groupby("plot_id").ba.sum()
    all_n = live.groupby("plot_id").w.sum()
    pine = live[live.species == species_code].groupby("plot_id").agg(
        **{f"n{suffix}": ("w", "sum"), f"ba{suffix}": ("ba", "sum")})
    pine[f"share{suffix}"] = pine[f"ba{suffix}"]/total_ba
    pine[f"all_n{suffix}"] = all_n
    return pine.reset_index()


def metrics(frame: pd.DataFrame, predicted: str, observed: str) -> dict:
    error = frame[predicted]-frame[observed]
    return {"plots": len(frame), "bias_mg_ha": round(float(error.mean()), 3),
            "mae_mg_ha": round(float(error.abs().mean()), 3),
            "rmse_mg_ha": round(float(np.sqrt(np.mean(error**2))), 3)}


def evaluate_species(species_code: int, species_name: str) -> dict:
    species = species_name.lower().replace(" ", "_")
    plots_parts, tree_parts = [], []
    for code, ifn4_name in IFN4_FILE.items():
        p3 = IFN3_DIR / f"p{code}" / f"Ifn3p{code}.accdb"
        p4 = IFN4_DIR / f"{ifn4_name}.accdb"
        cache = PAIRS_DIR / f"pairs34_p{code}.csv"
        if not (p3.exists() and p4.exists() and cache.exists()):
            continue
        pairs = pd.read_csv(cache, usecols=["plot_key", "species_code", "years"])
        pairs = pairs[(pairs.species_code == species_code) & pairs.years.between(5,25)]
        if pairs.empty:
            continue
        pairs["plot_id"] = str(int(code)) + "_" + pairs.plot_key.str.replace("_", "", regex=False)
        links = pairs.groupby("plot_id").agg(links=("years", "size"), years=("years", "median"))
        t3 = trees_from_db(p3, code, 3)
        t4 = trees_from_db(p4, code, 4)
        stocks = plot_stocks(t3, species_code, "3").merge(
            plot_stocks(t4, species_code, "4"), on="plot_id", validate="one_to_one")
        stocks = stocks.merge(links, on="plot_id", validate="one_to_one")
        stocks = stocks[(stocks.links >= 2) & (stocks.share3 >= .8) &
                        (stocks.share4 >= .8) & (stocks.n3 >= 100) & (stocks.n4 > 0)].copy()
        if stocks.empty:
            continue
        ids = set(stocks.plot_id)
        old = t3[t3.plot_id.isin(ids) & t3.live & (t3.species == species_code)].copy()
        future = t4[t4.plot_id.isin(ids)].copy()
        former = future[(future.ord3 > 0) & ~future.ord3.isin(BAD_ORDER)]
        ambiguous = set(old.loc[old.duplicated(["plot_id", "ord3"], keep=False), "plot_id"])
        ambiguous.update(former.loc[former.duplicated(["plot_id", "ord3"], keep=False), "plot_id"])
        if ambiguous:
            stocks = stocks[~stocks.plot_id.isin(ambiguous)]
            old = old[~old.plot_id.isin(ambiguous)]
            future = future[~future.plot_id.isin(ambiguous)]
            former = former[~former.plot_id.isin(ambiguous)]
        if stocks.empty:
            continue
        status = former[["plot_id", "ord3", "ord4", "species", "live"]].rename(
            columns={"species": "species4", "live": "live4"})
        old = old.merge(status, on=["plot_id", "ord3"], how="left", validate="one_to_one")
        old["category"] = np.select(
            [old.ord4 == 0, old.ord4 == 888, old.ord4 == 444,
             old.ord4.isin((666,777,999)),
             (old.ord4 > 0) & old.live4.astype("boolean").fillna(False) & (old.species4 == species_code)],
            ["cut", "dead", "not_found", "other", "survived"], default="unrecorded")
        new = future[future.live & (future.species == species_code) & (future.ord3 == 0)].copy()
        new["new_agb"] = [w*BIOMASS[species](d,h).aboveground_kg/1000.0
                          for d,h,w in zip(new.d_cm,new.h_m,new.w)]
        new_agb = new.groupby("plot_id").new_agb.sum().rename("new_agb")
        all4 = future[future.live & (future.species == species_code)].copy()
        all4["agb4"] = [w*BIOMASS[species](d,h).aboveground_kg/1000.0
                        for d,h,w in zip(all4.d_cm,all4.h_m,all4.w)]
        agb4 = all4.groupby("plot_id").agb4.sum().rename("agb4")
        category = old.pivot_table(index="plot_id", columns="category", values="w",
                                   aggfunc="sum", fill_value=0.0)
        for name in ("cut","dead","not_found","other","unrecorded"):
            if name not in category:
                category[name] = 0.0
        stocks = stocks.join(category, on="plot_id")
        stocks = stocks.join(new_agb, on="plot_id").join(agb4, on="plot_id")
        stocks["new_agb"] = stocks.new_agb.fillna(0.0)
        stocks["province_code"] = code
        old["province_code"] = code
        plots_parts.append(stocks)
        tree_parts.append(old)
    if not plots_parts:
        raise RuntimeError(f"No repeated pine-dominated plots for {species}")
    plots = pd.concat(plots_parts, ignore_index=True)
    trees = pd.concat(tree_parts, ignore_index=True)
    all_plots, all_trees = plots.copy(), trees.copy()
    plots = plots[(plots.cut == 0) & (plots.other == 0) & (plots.unrecorded == 0)].copy()
    trees = trees[trees.plot_id.isin(set(plots.plot_id))].merge(
        plots[["plot_id", "years", "n3"]], on="plot_id", validate="many_to_one")
    fit = plots[plots.province_code.isin(FIT_CODES[species_code])]
    held = plots[~plots.province_code.isin(FIT_CODES[species_code])]
    if len(fit) < 15 or len(held) < 15:
        raise RuntimeError(f"Insufficient uncut fit/test plots for {species}: {len(fit)}/{len(held)}")
    fit_trees = trees[trees.province_code.isin(FIT_CODES[species_code])]
    lost = fit_trees.category.isin(("dead", "not_found")).to_numpy().astype(float)
    weights = fit_trees.w.to_numpy()/fit_trees.n3.to_numpy()
    years = fit_trees.years.to_numpy()
    def neg_log_likelihood(log_hazard: float) -> float:
        p = np.clip(-np.expm1(-years*np.exp(log_hazard)), 1e-9, 1-1e-9)
        return float(-np.sum(weights*(lost*np.log(p)+(1-lost)*np.log1p(-p))))
    hazard_result = minimize_scalar(neg_log_likelihood, bounds=(-9,-1), method="bounded")
    if not hazard_result.success:
        raise RuntimeError(f"Hazard fit failed for {species}")
    hazard = float(np.exp(hazard_result.x))
    annual_new_agb = float((fit.new_agb/fit.years).mean())
    # Residual distribution uses only designated fit-province strict pairs.
    # The independent survivor-pair calibration is saved by the separate
    # cohort evaluation; its residual moments are fit-province only.
    import json
    cohort_file = Path(__file__).parent / "results/other_pines_cohort_agb.json"
    cohort = json.loads(cohort_file.read_text())["species"][species_name]
    residual_mean = cohort["fit_dbh_residual_mean_cm"]
    residual_sd = cohort["fit_dbh_residual_sd_cm"]
    frozen = PARAMETERS[species]
    if not np.allclose((hazard,annual_new_agb,residual_mean,residual_sd),
                       (frozen.annual_natural_hazard,frozen.annual_new_tree_agb_mg_ha,
                        frozen.dbh_residual_mean_cm,frozen.dbh_residual_sd_cm),atol=1e-7):
        raise RuntimeError(f"Frozen inventory parameters differ from fit for {species}")
    forecast_rows = []
    for plot_id, group in trees.groupby("plot_id"):
        row = plots.loc[plots.plot_id == plot_id].iloc[0]
        initial = [InventoryTree(float(t.d_cm),float(t.h_m)) for t in group.itertuples()]
        forecast = project_inventory(species,initial,float(row.all_n3),float(row.years))
        forecast_rows.append((plot_id,forecast.initial_agb_mg_ha,
                              forecast.surviving_cohort_agb_mg_ha))
    forecast = pd.DataFrame(forecast_rows,columns=("plot_id","initial_agb","old_pred"))
    plots = plots.merge(forecast,on="plot_id",validate="one_to_one")
    plots["total_pred"] = plots.old_pred + annual_new_agb*plots.years
    plots["old_observed"] = plots.agb4-plots.new_agb
    held = plots[~plots.province_code.isin(FIT_CODES[species_code])]
    rng = np.random.default_rng(20261023+species_code)
    sample = rng.integers(0,len(held),size=(2000,len(held)))
    initial_error = held.initial_agb.to_numpy()-held.agb4.to_numpy()
    total_error = held.total_pred.to_numpy()-held.agb4.to_numpy()
    cohort_error = held.old_pred.to_numpy()-held.agb4.to_numpy()
    gain_initial = np.sqrt(np.mean(initial_error[sample]**2,axis=1)) - np.sqrt(
        np.mean(total_error[sample]**2,axis=1))
    gain_entry = np.sqrt(np.mean(cohort_error[sample]**2,axis=1)) - np.sqrt(
        np.mean(total_error[sample]**2,axis=1))
    cut_plots = all_plots[(all_plots.cut > 0) & (all_plots.other == 0) &
                          (all_plots.unrecorded == 0) &
                          ~all_plots.province_code.isin(FIT_CODES[species_code])].copy()
    cut_rows = []
    for row in cut_plots.itertuples():
        group = all_trees[all_trees.plot_id == row.plot_id]
        initial = [InventoryTree(float(t.d_cm),float(t.h_m),
                                 1.0 if t.category == "cut" else 0.0)
                   for t in group.itertuples()]
        forecast = project_inventory(species,initial,float(row.all_n3),float(row.years))
        cut_trees = group[group.category == "cut"]
        cut_initial_agb = float(sum(
            t.w*BIOMASS[species](t.d_cm,t.h_m).aboveground_kg/1000.0
            for t in cut_trees.itertuples()))
        cut_rows.append((row.plot_id,row.agb4-row.new_agb,
                         forecast.surviving_cohort_agb_mg_ha,
                         forecast.initial_agb_mg_ha-cut_initial_agb))
    cut_frame = pd.DataFrame(cut_rows,columns=("plot_id","old_observed","model","initial_less_cut"))
    return {"fit_plots": len(fit), "heldout_plots": len(held),
            "heldout_interval_years": [round(float(held.years.min()),2),
                                       round(float(held.years.max()),2)],
            "fit_annual_natural_hazard": hazard,
            "fit_annual_new_tree_agb_mg_ha": annual_new_agb,
            "heldout_by_province": held.province_code.value_counts().to_dict(),
            "heldout_mean_new_tree_agb_mg_ha": round(float(held.new_agb.mean()),3),
            "heldout_total_rmse_gain_over_initial_95pct":
                [round(float(x),3) for x in np.quantile(gain_initial,[.025,.975])],
            "heldout_ingrowth_rmse_gain_over_cohort_95pct":
                [round(float(x),3) for x in np.quantile(gain_entry,[.025,.975])],
            "heldout_old_cohort": {
                "forecast": metrics(held,"old_pred","old_observed"),
                "initial_agb": metrics(held,"initial_agb","old_observed")},
            "heldout_total": {
                "forecast_with_fit_mean_entry": metrics(held,"total_pred","agb4"),
                "cohort_only": metrics(held,"old_pred","agb4"),
                "initial_agb": metrics(held,"initial_agb","agb4")},
            "heldout_supplied_cut_scenario": {
                "plots": len(cut_frame),
                "forecast_old_cohort": metrics(cut_frame,"model","old_observed")
                if len(cut_frame) else None,
                "initial_less_cut": metrics(cut_frame,"initial_less_cut","old_observed")
                if len(cut_frame) else None}}


def evaluate() -> dict:
    return {"design": "Pine-dominated IFN3→IFN4 plots without recorded target-pine cuts; growth, constant natural hazard, IFN sampling, fit-province mean new-tree AGB; other provinces held out; no age input",
            "species": {name: evaluate_species(code, name) for code,name in SPECIES.items()}}


if __name__ == "__main__":
    import json
    print(json.dumps(evaluate(), indent=2, ensure_ascii=False))
