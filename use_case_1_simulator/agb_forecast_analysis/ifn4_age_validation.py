"""
Stand AGB block (spain_stand_agb.py) against IFN4 plots of KNOWN AGE — the validation ESA CCI cannot give
(spain_stand_cci_comparison.py: CCI has a ~50 Mg/ha floor on cleared stands).

Data: IFN4 field databases (MITECO), Galicia (A Coruña, Lugo, Ourense, Pontevedra; surveyed 2008-09) and
Extremadura, in /media/HD/Data/ifn4/. Stand age comes from PCEspParc (main species only, PosEsp = 1):
`Edad` is recorded only for even-aged stands (FPMasa 1 coetánea / 2 regular). Kept when the age is trustworthy:
not "by appearance" (FInfor 1), and either bored (FInfor 8), whorl-counted (5, valid for young conifers) or
rated medium/high reliability (Fiabil 3/4). See the IFN4 field manual (docs in /media/HD/Data/ifn4/).

Measured stand per plot (main species only): trees >= 7.5 cm from PCMayores (concentric-plot expansion
factors, measured d and h) plus "pies menores" 2.5-7.5 cm from PCRegenera (category 4, counted in the 5 m
subplot, mean height recorded; d taken as 5 cm). AGB from the same Ruiz-Peinado equations the model uses, so
the test is about growth, density and stand structure, not the biomass equations. Plots kept if the species
holds >= 80% of the basal area. H0 = mean height of the 100 thickest trees/ha (Assmann).

Age since planting = Edad - PLANTING_LAG_YEARS, the SHARED constant in spain_stand_agb.py (1-year nursery
stock), which the model then adds back to reach the papers' total age. See the AGE CONVENTION block there:
the IFN4 manual never says what `Edad` counts from, and the origin depends on how it was obtained
(whorl count ~ total age, boring ~ age at the boring height, consultation ~ planting year). We treat it as
a total age; that is an assumption, and the spread between conventions (0-4 y) exceeds the nursery lag.

Tests:
  A. All defaults at the plot's age vs measured (AGB, trees/ha): the tool as a user without data would run it.
  B. Pine only: model given the plot's measured H0 (site anchor) and trees/ha (no survival loss, no thinning):
     isolates the height -> diameter -> biomass -> mean/dominant-ratio chain.
  C. Mean/dominant tree AGB ratio measured on these even-aged plots vs the IFN3 table the model uses.
"""

from __future__ import annotations

import os  # Configurable external inventory paths.

import math
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ifn_remeasurement_validation import _mdb, _num, trees_per_ha  # noqa: E402
from spain_biomass_equations import pinus_pinaster_biomass, quercus_ilex_biomass, quercus_suber_biomass  # noqa: E402
from spain_growth_curves import pine_ecoregion_from_province  # noqa: E402
from spain_height_diameter import hd_params, height_from_diameter  # noqa: E402
from spain_stand_agb import (  # noqa: E402
    PLANTING_LAG_YEARS,
    StandInputs,
    mean_to_dominant_ratio,
    project_stand_agb,
)
from stand_structure_calibration import dominant_mean  # noqa: E402

IFN4_DIR = (os.environ.get("IFN_DATA_ROOT", os.path.expanduser("~/ifn-data")) + "/ifn4")
OUT_DIR = (os.environ.get("IFN_DATA_ROOT", os.path.expanduser("~/ifn-data")) + "/ifn_validation")
DATABASES = {
    "A Coruña": "Ifn4p15-coruna", "Lugo": "Ifn4p27-lugo", "Ourense": "Ifn4p32-ourense",
    "Pontevedra": "Ifn4p36-pontevedra", "Extremadura": "Ifn4_Extremadura",
}
# INE province code per Galician region, for the pine ecoregion (A3; Vega et al. 1993 name Pontevedra
# and A Coruña as the coastal improvement zone). Extremadura is outside the curve's range entirely.
PROVINCE_CODE = {"A Coruña": 15, "Lugo": 27, "Ourense": 32, "Pontevedra": 36}
SPECIES = {26: ("pinus_pinaster", pinus_pinaster_biomass), 45: ("quercus_ilex", quercus_ilex_biomass),
           46: ("quercus_suber", quercus_suber_biomass)}
MIN_BA_SHARE = 0.8
SMALL_TREE_D_CM = 5.0
SMALL_TREE_PLOT_RADIUS_M = 5.0
AGE_BINS = [0, 5, 10, 15, 20, 30, 45, 80]
# Pine chain recalibrated on A Coruña + Lugo (ifn4_pine_calibration.py); only these provinces are a clean test.
TEST_REGIONS = ["Ourense", "Pontevedra"]


def _plot_key(df: pd.DataFrame, province_code: int | None = None) -> pd.Series:
    if "Provincia" in df:
        province = df.Provincia.astype(str).str.strip()
    elif province_code is not None:
        province = pd.Series(str(province_code), index=df.index)
    else:
        raise ValueError("province_code is required for inventories without Provincia")
    return (province + "_" + df.Estadillo.str.strip() + df.Cla.str.strip()
            + df.Subclase.fillna("").astype(str).str.strip())


def load_aged_plots(region: str, db: str, *, database_path: str | None = None,
                    province_code: int | None = None) -> pd.DataFrame:
    path = database_path or os.path.join(IFN4_DIR, f"{db}.accdb")
    e = _mdb(path, "PCEspParc")
    e["plot_id"] = _plot_key(e, province_code)
    for c in ("PosEsp", "Especie", "FPMasa", "Edad", "FInfor", "Fiabil", "OrgMasa1", "OrgMasa2"):
        e[c] = _num(e[c])
    e = e[(e.PosEsp == 1) & e.Especie.isin(list(SPECIES)) & e.FPMasa.isin([1, 2]) & e.Edad.notna()]
    reliable = (e.FInfor != 1) & (e.FInfor.isin([5, 8]) | e.Fiabil.isin([3, 4]))
    return e[reliable][["plot_id", "Especie", "Edad", "FInfor", "Fiabil", "OrgMasa1", "OrgMasa2"]].assign(region=region)


def measured_stands(db: str, aged: pd.DataFrame, *, database_path: str | None = None,
                    province_code: int | None = None) -> pd.DataFrame:
    path = database_path or os.path.join(IFN4_DIR, f"{db}.accdb")
    t = _mdb(path, "PCMayores")
    t["plot_id"] = _plot_key(t, province_code)
    t = t[t.plot_id.isin(set(aged.plot_id))].copy()
    t["species"] = _num(t.Especie)
    t["d_cm"] = (_num(t.Dn1) + _num(t.Dn2)) / 20.0
    t["h_m"] = _num(t.Ht)
    t = t[(t.d_cm >= 7.5) & (t.h_m > 1.3)]
    t["w"] = trees_per_ha(t.d_cm.to_numpy())
    t["ba"] = t.w * math.pi * (t.d_cm / 200.0) ** 2
    r = _mdb(path, "PCRegenera")
    r["plot_id"] = _plot_key(r, province_code)
    r = r[r.plot_id.isin(set(aged.plot_id)) & (_num(r.CatDes) == 4)].copy()
    r["species"], r["n"], r["h_m"] = _num(r.Especie), _num(r.NumPies), _num(r.Hm) / 10.0
    small_w = 10000.0 / (math.pi * SMALL_TREE_PLOT_RADIUS_M**2)

    rows = []
    for _, a in aged.iterrows():
        code = int(a.Especie)
        sp, bio = SPECIES[code]
        g = t[t.plot_id == a.plot_id]
        tot_ba = g.ba.sum() + 0.0
        m = g[g.species == code]
        small = r[(r.plot_id == a.plot_id) & (r.species == code)]
        n_small = float(small.n.fillna(0).sum()) * small_w
        if tot_ba > 0 and m.ba.sum() / tot_ba < MIN_BA_SHARE:
            continue
        if m.empty and n_small == 0:
            continue
        agb_trees = np.array([bio(d, h).aboveground_kg for d, h in zip(m.d_cm, m.h_m)])
        agb_large = float((agb_trees * m.w).sum()) / 1000.0
        agb_small = 0.0
        h_small = float("nan")
        if n_small > 0:
            h_small = small.h_m.dropna().mean()
            if not h_small > 1.3:
                h_small = height_from_diameter(SMALL_TREE_D_CM, hd_params("pinus_pinaster", "Atlantic", plantation=True))
            agb_small = n_small * bio(SMALL_TREE_D_CM, h_small).aboveground_kg / 1000.0
        n_large = float(m.w.sum())
        mean_h = ((m.h_m * m.w).sum() + (n_small * h_small if n_small > 0 else 0.0)) / (n_large + n_small)
        h0 = dominant_mean(m.h_m.to_numpy(), m.d_cm.to_numpy(), m.w.to_numpy()) if n_large > 0 else float("nan")
        ratio = float("nan")
        if n_large >= 100:
            ratio = (agb_large * 1000.0 / n_large) / dominant_mean(agb_trees, m.d_cm.to_numpy(), m.w.to_numpy())
        rows.append({
            "plot_id": a.plot_id, "region": a.region, "species": sp, "age": a.Edad - PLANTING_LAG_YEARS,
            # IFN4 manual annex 10: OrgMasa1=2 is artificial origin, whereas
            # OrgMasa2=2 specifically records planting (as opposed to sowing).
            "planted": a.OrgMasa1 == 2 and a.OrgMasa2 == 2, "age_source": int(a.FInfor), "trees_per_ha": n_large + n_small,
            "trees_per_ha_ge7.5": n_large, "agb_mg_ha": agb_large + agb_small,
            "agb_ge7.5_mg_ha": agb_large,
            "h0_m": h0, "mean_h_m": mean_h, "measured_ratio": ratio,
        })
    return pd.DataFrame(rows)


def large_tree_view(stands: pd.DataFrame) -> pd.DataFrame:
    """Restrict measured stands to trees >= 7.5 cm DBH, the IFN measured-tree set.

    `measured_stands` adds counted small trees (2.5-7.5 cm, PCRegenera) to stem
    number and AGB. The pine stand chain is calibrated on trees >= 7.5 cm only
    (ifn4_pine_calibration.py), so feeding it total stems mixes two populations:
    a plot with a sparse overstorey over thousands of saplings gets the overstorey
    tree's biomass multiplied by every sapling. This view makes stem number and
    observed AGB both refer to trees >= 7.5 cm and drops plots without such trees.
    """
    if stands.empty:
        return stands.copy()
    out = stands[stands["trees_per_ha_ge7.5"] > 0].copy()
    out["trees_per_ha"] = out["trees_per_ha_ge7.5"]
    out["agb_mg_ha"] = out["agb_ge7.5_mg_ha"]
    return out


def run_models(df: pd.DataFrame) -> pd.DataFrame:
    """Test A (all defaults), test B (measured H0 + trees/ha) and their per-plot-ecoregion variants (A3).

    A3 asks whether assigning each plot the ecoregion Vega et al. (1993) give its province beats the
    model's single hard-coded interior curve. "eco" columns are the same two tests with
    pine_ecoregion set per plot."""
    out = df.copy()
    a_agb, a_n, b_agb, model_ratio = [], [], [], []
    a_eco, b_eco, ecoregion, a_obs = [], [], [], []
    for _, p in out.iterrows():
        age = max(p.age, 1.0)
        res = project_stand_agb(StandInputs(p.species), [age])
        a_agb.append(res.agb_mg_ha[0])
        a_n.append(res.trees_per_ha[0])
        model_ratio.append(mean_to_dominant_ratio(p["trees_per_ha_ge7.5"]) if p["trees_per_ha_ge7.5"] > 0 else float("nan"))
        eco = pine_ecoregion_from_province(PROVINCE_CODE.get(p.region, -1)) if p.species == "pinus_pinaster" else None
        ecoregion.append(eco)
        if eco is not None:
            a_eco.append(project_stand_agb(StandInputs(p.species, pine_ecoregion=eco), [age]).agb_mg_ha[0])
            a_obs.append(project_stand_agb(
                StandInputs(p.species, pine_ecoregion=eco, management="galicia_observed"), [age]).agb_mg_ha[0])
        else:
            a_eco.append(float("nan"))
            a_obs.append(float("nan"))
        if p.species == "pinus_pinaster" and p.h0_m > 1.3 and p.trees_per_ha > 0:
            kw = dict(planting_density_per_ha=p.trees_per_ha, survival_checks={0.5: 1.0},
                      annual_mortality_after_last_check=0.0, thinnings={}, pine_site_anchor=(p.h0_m, age))
            b_agb.append(project_stand_agb(StandInputs("pinus_pinaster", **kw), [age]).agb_mg_ha[0])
            b_eco.append(project_stand_agb(StandInputs("pinus_pinaster", pine_ecoregion=eco or "interior", **kw),
                                           [age]).agb_mg_ha[0] if eco else float("nan"))
        else:
            b_agb.append(float("nan"))
            b_eco.append(float("nan"))
    out["modelA_agb"], out["modelA_trees"], out["modelB_agb"], out["model_ratio_ifn3"] = a_agb, a_n, b_agb, model_ratio
    out["modelA_eco_agb"], out["modelB_eco_agb"], out["ecoregion"] = a_eco, b_eco, ecoregion
    out["modelA_eco_obs_agb"] = a_obs
    return out


def _group(row) -> str:
    if row.species != "pinus_pinaster":
        return row.species
    if row.region in TEST_REGIONS:
        return "pinus_pinaster HELD-OUT (Ourense, Pontevedra)"
    return "pinus_pinaster Extremadura" if row.region == "Extremadura" else "pinus_pinaster fit provinces (A Coruña, Lugo)"


def summarize(df: pd.DataFrame) -> None:
    df = df.assign(group=df.apply(_group, axis=1))
    for sp, g in df.groupby("group"):
        g = g.assign(age_class=pd.cut(g.age, AGE_BINS))
        print(f"\n=== {sp}: {len(g)} plots ({g.region.value_counts().to_dict()}), planted {int(g.planted.sum())}")
        t = g.groupby("age_class", observed=True).agg(
            n=("plot_id", "size"), trees_meas=("trees_per_ha", "median"), trees_modelA=("modelA_trees", "median"),
            agb_meas=("agb_mg_ha", "median"), agb_modelA=("modelA_agb", "median"), agb_modelB=("modelB_agb", "median"),
            agb_modelA_obs=("modelA_eco_obs_agb", "median"),
        )
        print(t.round(1).to_string())
        for col, name in (("modelA_agb", "A (defaults)"), ("modelB_agb", "B (measured H0 + trees/ha)"),
                         ("modelA_eco_agb", "A + per-plot ecoregion"), ("modelB_eco_agb", "B + per-plot ecoregion"),
                         ("modelA_eco_obs_agb", "A + ecoregion + galicia_observed thinning")):
            ok = g[col].notna() & (g.agb_mg_ha > 1) & (g.age >= 10)
            if ok.sum():
                r = (g[col][ok] / g.agb_mg_ha[ok])
                print(f"  ages >= 10: model {name} / measured AGB: median {r.median():.2f} (IQR {r.quantile(.25):.2f}-{r.quantile(.75):.2f}), n={ok.sum()}")
        c = g.dropna(subset=["measured_ratio"])
        if len(c):
            c = c.assign(n_class=pd.cut(c["trees_per_ha_ge7.5"], [100, 300, 600, 1000, 3000]))
            print("  test C, mean/dominant AGB ratio (even-aged IFN4 vs IFN3 table the model uses):")
            print(c.groupby("n_class", observed=True).agg(n=("plot_id", "size"), measured=("measured_ratio", "median"),
                                                         ifn3_table=("model_ratio_ifn3", "median")).round(2).to_string())


if __name__ == "__main__":
    frames = []
    for region, db in DATABASES.items():
        aged = load_aged_plots(region, db)
        frames.append(measured_stands(db, aged))
    df = run_models(pd.concat(frames, ignore_index=True))
    os.makedirs(OUT_DIR, exist_ok=True)
    df.to_csv(os.path.join(OUT_DIR, "ifn4_aged_plots_vs_model.csv"), index=False)
    summarize(df)
