"""Separate recorded cuts, deaths, disappearances and new pine trees.

IFN4 field guide annex 39 defines OrdenIf4 000=cut, 444=not found with no
clear cutting evidence, 888=dead tree present, and OrdenIf3 000=new sample
tree. This audit is retrospective, since future cuts are management inputs,
not values a growth model can know in advance.
"""

from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "use_case_1_simulator/agb_forecast_analysis"))

from evaluate_mediterranean_ingrowth import load_new_trees  # noqa: E402
from fit_mediterranean_density import IFN3_DIR, PROVINCES, load_remeasurements  # noqa: E402
from ifn4_age_validation import IFN4_DIR, _plot_key  # noqa: E402
from ifn_remeasurement_validation import _mdb, _num, trees_per_ha  # noqa: E402

CATEGORIES = ("cut", "dead", "not_found", "survived", "other", "unrecorded")


def prior_tree_status(code: int, database: str, eligible_ids: set[str]) -> pd.DataFrame:
    p3 = IFN3_DIR / f"p{code:02d}" / f"Ifn3p{code:02d}.accdb"
    p4 = Path(IFN4_DIR) / f"{database}.accdb"
    old = _mdb(str(p3), "PCMayores")
    old["plot_id"] = _plot_key(old, code)
    old = old[old.plot_id.isin(eligible_ids)].copy()
    old["species"] = _num(old.Especie)
    old["d_cm"] = (_num(old.Dn1) + _num(old.Dn2)) / 20.0
    old["h_m"] = _num(old.Ht)
    old["ord3"] = _num(old.OrdenIf3)
    old = old[(old.species == 26) & (old.d_cm >= 7.5) & (old.h_m > 1.3) &
              (old.ord3 > 0) & ~old.ord3.isin((444, 666, 777, 888, 999))].copy()
    old["w3"] = trees_per_ha(old.d_cm.to_numpy())
    current = _mdb(str(p4), "PCMayores")
    current["plot_id"] = _plot_key(current, code)
    current = current[current.plot_id.isin(eligible_ids)].copy()
    current["ord3"] = _num(current.OrdenIf3)
    current["ord4"] = _num(current.OrdenIf4)
    current["species4"] = _num(current.Especie)
    current["d4_cm"] = (_num(current.Dn1) + _num(current.Dn2)) / 20.0
    current["h4_m"] = _num(current.Ht)
    current = current[(current.ord3 > 0) &
                      ~current.ord3.isin((444, 666, 777, 888, 999))].copy()
    current_live = (current.ord4 > 0) & ~current.ord4.isin((444, 666, 777, 888, 999)) & (
        (current.species4 == 26) & (current.d4_cm >= 7.5) & (current.h4_m > 1.3))
    current["w4"] = np.where(current_live, trees_per_ha(current.d4_cm.to_numpy()), 0.0)
    if old.duplicated(["plot_id", "ord3"]).any() or current.duplicated(["plot_id", "ord3"]).any():
        raise RuntimeError(f"Ambiguous previous-tree order in province {code}")
    joined = old[["plot_id", "ord3", "w3", "d_cm", "h_m"]].rename(
        columns={"d_cm": "d3_cm", "h_m": "h3_m"}).merge(
        current[["plot_id", "ord3", "ord4", "w4", "d4_cm", "h4_m"]], on=["plot_id", "ord3"],
        how="left", validate="one_to_one")
    order = joined.ord4
    joined["category"] = np.select(
        [order == 0, order == 888, order == 444,
         (order > 0) & ~order.isin((444, 666, 777, 888, 999)),
         order.isin((666, 777, 999))],
        ["cut", "dead", "not_found", "survived", "other"],
        default="unrecorded",
    )
    joined["survivor_expansion_delta"] = np.where(
        joined.category == "survived", joined.w4.fillna(0.0) - joined.w3, 0.0)
    return joined


def load_turnover() -> tuple[pd.DataFrame, pd.DataFrame]:
    plots = load_remeasurements()
    new_trees = load_new_trees()
    plots = plots.merge(new_trees, on="plot_id", how="left", validate="one_to_one")
    plots[["new_n_ha", "new_ba_m2_ha", "new_agb_mg_ha"]] = plots[
        ["new_n_ha", "new_ba_m2_ha", "new_agb_mg_ha"]].fillna(0.0)
    parts = []
    for code, (_, database) in PROVINCES.items():
        ids = set(plots.plot_id[plots.plot_id.str.startswith(f"{code}_")])
        if ids:
            parts.append(prior_tree_status(code, database, ids))
    trees = pd.concat(parts, ignore_index=True)
    component = trees.pivot_table(index="plot_id", columns="category", values="w3", aggfunc="sum", fill_value=0.0)
    for category in CATEGORIES:
        if category not in component:
            component[category] = 0.0
    plots = plots.merge(component.reset_index(), on="plot_id", validate="one_to_one")
    shift = trees.groupby("plot_id").survivor_expansion_delta.sum().rename("survivor_expansion_delta")
    plots = plots.merge(shift, on="plot_id", validate="one_to_one")
    # Check that the category decomposition covers the exact initial trees
    # used in the density analysis.
    total = plots[list(CATEGORIES)].sum(axis=1)
    if not np.allclose(total, plots.n3, rtol=1e-8, atol=1e-6):
        raise RuntimeError("IFN3 turnover categories do not sum to initial pine density")
    plots["accounted_n4"] = plots.survived + plots.survivor_expansion_delta + plots.new_n_ha
    plots["accounting_residual"] = plots.accounted_n4 - plots.n4
    return plots, trees


def evaluate() -> dict:
    plots, trees = load_turnover()
    out = {
        "design": "IFN3→IFN4 pine-dominated strict-link plots; IFN4 annex-39 tree outcome codes; weights at IFN3 diameter",
        "plots": len(plots), "prior_pine_trees": len(trees),
        "prior_tree_status_counts": trees.category.value_counts().to_dict(),
        "weighted_pine_stems_ha_median": {c: round(float(plots[c].median()), 3) for c in CATEGORIES},
        "weighted_pine_stems_ha_mean": {c: round(float(plots[c].mean()), 3) for c in CATEGORIES},
        "new_stems_ha_median": round(float(plots.new_n_ha.median()), 3),
        "mean_survivor_expansion_delta_stems_ha": round(float(plots.survivor_expansion_delta.mean()), 3),
        "median_survivor_expansion_delta_stems_ha": round(float(plots.survivor_expansion_delta.median()), 3),
        "accounting_residual_rmse_stems_ha": round(float(np.sqrt(np.mean(plots.accounting_residual**2))), 3),
        "accounting_residual_mean_stems_ha": round(float(plots.accounting_residual.mean()), 3),
        "share_plots_with_recorded_cut": round(float((plots.cut > 0).mean()), 3),
        "share_plots_with_recorded_death": round(float((plots.dead > 0).mean()), 3),
        "share_plots_with_unrecorded_prior_pine": round(float((plots.unrecorded > 0).mean()), 3),
        "by_province": {},
    }
    for province, group in plots.groupby("province"):
        out["by_province"][province] = {
            "plots": len(group),
            "mean_cut_stems_ha": round(float(group.cut.mean()), 3),
            "mean_dead_stems_ha": round(float(group.dead.mean()), 3),
            "mean_not_found_stems_ha": round(float(group.not_found.mean()), 3),
            "mean_new_stems_ha": round(float(group.new_n_ha.mean()), 3),
            "share_with_cut": round(float((group.cut > 0).mean()), 3),
        }
    return out


if __name__ == "__main__":
    import json
    print(json.dumps(evaluate(), indent=2, ensure_ascii=False))
