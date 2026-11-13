"""
Validation of the Sylvestris growth model against Spanish National Forest Inventory
REMEASUREMENTS: the same trees measured in IFN2 (Caceres 1991-93) and IFN3 (2001),
~9-10 years apart. This is independent data — none of the growth equations in
spain_growth_curves.py were fit on it.

Linking: IFN3 trees in class A1 plots (plots re-found from IFN2) carry OrdenIf2, the
tree's IFN2 order number in the same plot. A link is kept only if it is the same
physical tree by three independent criteria: same species code, bearing within 10
grads and distance within 0.5 m of the IFN2 record ("strict pairs"). Caceres links
are good (86% species / 84% bearing agreement before filtering); Pontevedra's are
poor (57%) and are not used.

Why age-free tests only: IFN2/IFN3 record no tree age (IFN2 has a stand-age field on
~380 of 3,492 Caceres plots, mostly young stands), so the age-dependent curves (E1,
pine and cork oak site curves) cannot be tested tree by tree. What CAN be tested:
  1. Q. ilex age-independent diameter model T1 (Gea Izquierdo et al. 2007): needs only
     current diameter, stand density and elapsed years.
  2. (planned) spatial climate response of observed growth vs Gea-Izquierdo et al. 2013.

Data: /media/HD/Data/ifn2/<prov>/ (MITECO IFN2 "estadillos" DBF) and
/media/HD/Data/ifn3/<prov>/ (IFN3 Access DB, exported with mdbtools).
"""

from __future__ import annotations

import os  # Configurable external inventory paths.

import glob
import math
import os
import subprocess
import sys
from io import StringIO

import numpy as np
import pandas as pd
from dbfread import DBF

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from spain_growth_curves import (  # noqa: E402
    CORK_THICKNESS_ASSUMED_CM,
    pinus_pinaster_project,
    quercus_ilex_age_independent_project,
    quercus_ilex_diameter_from_planting,
    quercus_suber_diameter_default,
    quercus_suber_from_site_class,
)
from spain_height_diameter import diameter_from_height, hd_params  # noqa: E402
from spain_climate_modifier import quercus_ilex_growth_index  # noqa: E402

PLOT_CLIMATE = (os.environ.get("IFN_DATA_ROOT", os.path.expanduser("~/ifn-data")) + "/ifn_validation/plot_climate_caceres.csv")  # dataset/fetch_ifn_plot_climate.py

MDB_EXPORT = os.environ.get("MDB_EXPORT", "mdb-export")
# Province code -> data dirs. All 50 provinces downloaded 2026-09-25 from MITECO (IFN2 "estadillos"
# pages ifn2_parcelas_1_25/26_50, IFN3 pages ifn3_base_datos_1_25/26_50) into
# /media/HD/Data/ifn2/pXX/ and /media/HD/Data/ifn3/pXX/.
PROVINCE_NAMES = {"caceres": "10", "pontevedra": "36"}
PAIRS_CACHE = (os.environ.get("IFN_DATA_ROOT", os.path.expanduser("~/ifn-data")) + "/ifn_validation/pairs_p{code}.csv")


def _find(pattern):
    hits = glob.glob(pattern) + glob.glob(pattern.upper()) + glob.glob(pattern.lower())
    if not hits:
        raise FileNotFoundError(pattern)
    return hits[0]
SPECIES = {26: "pinus_pinaster", 45: "quercus_ilex", 46: "quercus_suber"}


def _num(s):
    # IFN2 DBFs in some provinces use a decimal comma ("02,3"); IFN3 files use a point.
    if s.dtype == object:
        s = s.astype(str).str.strip().str.replace(",", ".", regex=False)
    return pd.to_numeric(s, errors="coerce")


def _mdb(path, table):
    out = subprocess.run([MDB_EXPORT, path, table], capture_output=True, text=True, check=True).stdout
    return pd.read_csv(StringIO(out), dtype=str)


def plot_radius_m(d_cm):
    """IFN concentric subplots: trees >=7.5/12.5/22.5/42.5 cm are measured within 5/10/15/25 m."""
    return np.select([d_cm >= 42.5, d_cm >= 22.5, d_cm >= 12.5], [25.0, 15.0, 10.0], 5.0)


def trees_per_ha(d_cm):
    return 10000.0 / (math.pi * plot_radius_m(d_cm) ** 2)


def load_pairs(prov: str, use_cache: bool = True) -> pd.DataFrame:
    """Strict IFN2-IFN3 tree pairs for a province (code "06" or a name in PROVINCE_NAMES)."""
    code = PROVINCE_NAMES.get(prov, prov)
    cache = PAIRS_CACHE.format(code=code)
    if use_cache and os.path.exists(cache):
        return pd.read_csv(cache, parse_dates=["date2", "date3"])
    t2 = pd.DataFrame(iter(DBF(_find((os.environ.get("IFN_DATA_ROOT", os.path.expanduser("~/ifn-data")) + f"/ifn2/p{code}/PIESMA{code}.DBF")), encoding="latin1")))
    p2 = pd.DataFrame(iter(DBF(_find((os.environ.get("IFN_DATA_ROOT", os.path.expanduser("~/ifn-data")) + f"/ifn2/p{code}/DATEST{code}.DBF")), encoding="latin1")))
    accdb = _find((os.environ.get("IFN_DATA_ROOT", os.path.expanduser("~/ifn-data")) + f"/ifn3/p{code}/Ifn3p{code}.accdb"))
    t3 = _mdb(accdb, "PCMayores").rename(columns={"Distancia": "Distanci"})  # name differs by province
    p3 = _mdb(accdb, "PCParcelas")

    t2["est"], t2["ord2"] = _num(t2.ESTADILLO), _num(t2.NUMORDEN)
    t2["d2_cm"] = (_num(t2.DIAMETRO1) + _num(t2.DIAMETRO2)) / 20.0  # mm -> cm, mean of two
    # IFN2 plot density at the start of the interval (all species, trees >= 7.5 cm).
    dens = t2.assign(n=trees_per_ha(t2.d2_cm)).groupby("est").n.sum().rename("density_per_ha")

    p2["est"] = _num(p2.ESTADILLO)
    fech = _num(p2.FECHGR).astype("Int64").astype(str).str.zfill(6)  # some provinces store it as float
    p2["date2"] = pd.to_datetime(fech, format="%y%m%d", errors="coerce")
    p3 = p3[(p3.Cla == "A") & (p3.Subclase == "1")].copy()
    p3["est"] = _num(p3.Estadillo)
    p3["date3"] = pd.to_datetime(p3.FechaIni, format="%m/%d/%y %H:%M:%S", errors="coerce")

    t3 = t3[(t3.Cla == "A") & (t3.Subclase == "1")].copy()
    t3["est"], t3["ord2"] = _num(t3.Estadillo), _num(t3.OrdenIf2)
    t3 = t3[t3.ord2 > 0]
    t3["d3_cm"] = (_num(t3.Dn1) + _num(t3.Dn2)) / 20.0

    m = t3.merge(t2, on=["est", "ord2"], suffixes=("_3", "_2"))
    link_stats = {
        "province": code,
        "linked": len(m),
        "species_agree": float((_num(m.Especie) == _num(m.ESPECIE)).mean()) if len(m) else float("nan"),
    }
    bearing_diff = ((_num(m.Rumbo) - _num(m.RUMBO) + 200) % 400 - 200).abs()  # centesimal degrees
    strict = (
        (_num(m.Especie) == _num(m.ESPECIE))
        & (bearing_diff <= 10)
        & ((_num(m.Distanci) - _num(m.DISTANCI)).abs() <= 0.5)
    )
    m = m[strict].copy()
    m["species_code"] = _num(m.Especie).astype(int)
    m["h2_m"], m["h3_m"] = _num(m.ALTURA), _num(m.Ht)
    m["cork_param_ifn3"] = m.ParEsp.fillna("").str.strip()
    m = m.merge(p2[["est", "date2"]], on="est").merge(p3[["est", "date3"]], on="est").join(dens, on="est")
    m["years"] = (m.date3 - m.date2).dt.days / 365.25
    cols = ["est", "ord2", "species_code", "d2_cm", "d3_cm", "h2_m", "h3_m", "date2", "date3", "years", "density_per_ha", "cork_param_ifn3"]
    out = m[cols].dropna(subset=["d2_cm", "d3_cm", "years"]).reset_index(drop=True)
    out.attrs["link_stats"] = {**link_stats, "strict_pairs": len(out)}
    os.makedirs(os.path.dirname(cache), exist_ok=True)
    out.to_csv(cache, index=False)
    return out


def test_quercus_ilex_t1(pairs: pd.DataFrame) -> pd.DataFrame:
    """T1 (Gea Izquierdo et al. 2007) vs observed diameter, against two baselines it must beat:
    no growth (persistence) and a constant mean annual increment (the population average)."""
    q = pairs[(pairs.species_code == 45) & (pairs.d2_cm >= 10)].copy()
    n_all = len(q)
    # T1's domain: its fitted density range (Gea Izquierdo et al. 2007, Table 2: 39.5-210 trees/ha) and
    # trees below its density-dependent asymptote 1/a = N / ad1 (above it the age inversion is undefined).
    q = q[q.density_per_ha.between(39.5, 210.4) & (q.d2_cm < q.density_per_ha / 0.505583)]
    print(f"Q. ilex T1 domain: kept {len(q)} of {n_all} trees (rest outside fitted density range or above asymptote)")
    q["pred_t1"] = [
        quercus_ilex_age_independent_project(d, n, y) for d, n, y in zip(q.d2_cm, q.density_per_ha, q.years)
    ]
    mean_rate = ((q.d3_cm - q.d2_cm) / q.years).mean()
    q["pred_const"] = q.d2_cm + mean_rate * q.years
    q["obs_inc"], q["t1_inc"] = q.d3_cm - q.d2_cm, q.pred_t1 - q.d2_cm
    return q


def report_t1(q: pd.DataFrame) -> None:
    print(f"Q. ilex T1 test: {len(q)} trees, {q.est.nunique()} plots, interval {q.years.median():.1f} y, "
          f"density median {q.density_per_ha.median():.0f}/ha (T1 fit range 40-210)")
    for name, col in [("persistence (no growth)", "d2_cm"), ("constant mean increment", "pred_const"), ("T1", "pred_t1")]:
        err = q[col] - q.d3_cm
        print(f"  {name:26s} bias {err.mean():+.2f} cm  RMSE {np.sqrt((err**2).mean()):.2f} cm")
    r = np.corrcoef(q.obs_inc, q.t1_inc)[0, 1]
    print(f"  increment: observed mean {q.obs_inc.mean():.2f} cm, T1 mean {q.t1_inc.mean():.2f} cm, r(obs, T1) = {r:.2f}")
    q = q.assign(dclass=pd.cut(q.d2_cm, [10, 20, 30, 40, 60, 200]))
    print(q.groupby("dclass", observed=True).agg(n=("obs_inc", "size"), obs=("obs_inc", "median"), t1=("t1_inc", "median")).round(2).to_string())


# --- Default-curve prediction test -----------------------------------------------
# The deployed model, when a parcel's site is unknown, runs each species' default curve.
# Test exactly that: infer each tree's "equivalent age" on the default curve from its IFN2
# size, advance it by the observed interval, compare with IFN3. Dominant-tree models, so
# only each plot's largest trees of the species are used (top 2 by IFN2 diameter).

def _pine_height(t):  # interior ecoregion (Caceres), anchor = Alvarez-Gonzalez et al. 2005 Table I mean
    return pinus_pinaster_project(10.64, 19.3, t, ecoregion_coastal=False)


_PINE_HD = hd_params("pinus_pinaster", "Mediterranean")
CURVES = {
    # name: (species_code, size column pair, curve t -> size, note)
    "Q. ilex diameter (bridge + E1)": (45, ("d2_cm", "d3_cm"), quercus_ilex_diameter_from_planting, ""),
    "Q. suber diameter (2005 default site)": (
        46, ("d2_cm", "d3_cm"),
        lambda t: quercus_suber_diameter_default(t) + 2 * CORK_THICKNESS_ASSUMED_CM,
        "never-debarked trees only; over-cork = under-cork + 2 x assumed cork",
    ),
    "Q. suber height (2010 class III)": (46, ("h2_m", "h3_m"), lambda t: quercus_suber_from_site_class("III", t), ""),
    "P. pinaster height (interior, Table I site)": (26, ("h2_m", "h3_m"), _pine_height, ""),
    "P. pinaster diameter (height curve + h-d inverse)": (
        26, ("d2_cm", "d3_cm"), lambda t: diameter_from_height(max(_pine_height(t), 1.31), _PINE_HD), "",
    ),
}


def _invert_age(curve, size, lo=1.0, hi=400.0):
    """Age at which a monotone curve reaches `size` (bisection); NaN if outside the curve's range."""
    if not (curve(lo) <= size <= curve(hi)):
        return float("nan")
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if curve(mid) < size else (lo, mid)
    return 0.5 * (lo + hi)


# Gross-error screen. IFN height pairs have a heavy non-biological tail (Caceres pine: 16.5% of
# pairs lose height, 10% lose > 7.5 m: broken tops, a different tree, recording errors; IFN2
# heights are rounded to 0.5 m). Means and RMSE are dominated by these, so the test uses
# medians and median absolute error, after dropping losses beyond measurement noise.
MAX_LOSS = {"h2_m": 1.0, "d2_cm": 1.0}


def test_default_curves(pairs: pd.DataFrame, prov: str, curves: dict) -> None:
    print(f"\nDefault-curve prediction test, {prov} (largest 2 trees per plot and species; robust stats):")
    for name, (code, (c1, c2), curve, note) in curves.items():
        q = pairs[pairs.species_code == code].dropna(subset=[c1, c2])
        q = q[q[c1] > 1.3]
        if code == 46 and c1 == "d2_cm":
            q = q[q.cork_param_ifn3 == "1000"]  # "con corcho bornizo todo el arbol": never debarked
        q = q.sort_values("d2_cm", ascending=False).groupby("est").head(2)
        n0 = len(q)
        q = q[(q[c2] - q[c1]) >= -MAX_LOSS[c1]].copy()
        ages = np.array([_invert_age(curve, v) for v in q[c1]])
        ok = ~np.isnan(ages)
        q, ages = q[ok], ages[ok]
        if len(q) < 30:
            print(f"  {name}: n={len(q)}, too few")
            continue
        q["pred"] = [curve(a + y) for a, y in zip(ages, q.years)]
        obs_inc, pred_inc = q[c2] - q[c1], q.pred - q[c1]
        base = q[c1] + (obs_inc / q.years).median() * q.years  # median annual rate x interval
        mae = lambda e: np.median(np.abs(e))
        print(f"  {name}: n={len(q)} (dropped {n0 - len(q) - 0} of {n0}: gross loss or outside curve range){'; ' + note if note else ''}")
        print(f"    median increment obs {obs_inc.median():.2f} vs model {pred_inc.median():.2f} (ratio {pred_inc.median() / obs_inc.median():.2f}); "
              f"median abs error model {mae(q.pred - q[c2]):.2f} vs median-rate baseline {mae(base - q[c2]):.2f}; "
              f"Spearman(obs, model inc) {pd.Series(obs_inc.values).corr(pd.Series(pred_inc.values), method='spearman'):.2f}")
        cls = pd.qcut(q[c1], 4, duplicates="drop")
        tab = pd.DataFrame({"obs": obs_inc, "model": pred_inc, "cls": cls}).groupby("cls", observed=True).median()
        print("    by size quartile (median increment obs / model):",
              "; ".join(f"{k.left:.0f}-{k.right:.0f}: {r.obs:.2f} / {r.model:.2f}" for k, r in tab.iterrows()))


# --- Climate-response test (holm oak, Gea-Izquierdo et al. 2013) --------------------------
# Gea-Izquierdo's model was fitted ACROSS sites (4 sites along a temperature gradient), so a
# spatial test across Caceres plots is a fair one: size-adjusted basal-area increment (BAI)
# should scale with g1(P) * g2(Tmin), i.e. peak near Tmin ~9 C and fall on both sides. Caceres
# plots straddle that optimum (1990-2001 mean daily minimum 4.7-11.5 C, median 10.1 C).
# Model: log(BAI) = a + b log(d) + c log(density) + gamma log(g); gamma = 1 if the published response
# holds in magnitude, gamma = 0 if climate has no effect. Plot-level means, plot-bootstrap CIs.

def _plot_climate_index(clim_row: pd.Series, date2, date3) -> tuple[float, float, float]:
    """Mean hydro-year precip, mean daily Tmin and mean Gea index over the growth years between
    surveys: hydrological years Y (Nov Y-1..Oct Y) with Nov Y-1 after IFN2 and Oct Y before IFN3."""
    first = date2.year + 1 if date2.month < 11 else date2.year + 2
    last = date3.year if date3.month >= 10 else date3.year - 1
    years = [y for y in range(first, last + 1) if f"hydro_tmin_c_{y}" in clim_row.index]
    p = [clim_row[f"hydro_precip_mm_{y}"] for y in years]
    t = [clim_row[f"hydro_tmin_c_{y}"] for y in years]
    g = [quercus_ilex_growth_index(pp, tt) for pp, tt in zip(p, t)]
    return float(np.mean(p)), float(np.mean(t)), float(np.mean(g))


def test_climate_response(pairs: pd.DataFrame, n_boot: int = 2000, seed: int = 0) -> None:
    clim = pd.read_csv(PLOT_CLIMATE).set_index("est")
    q = pairs[(pairs.species_code == 45) & (pairs.d2_cm >= 10) & (pairs.d3_cm > pairs.d2_cm)].copy()
    n_nonpos = ((pairs.species_code == 45) & (pairs.d2_cm >= 10) & (pairs.d3_cm <= pairs.d2_cm)).sum()
    q["bai"] = math.pi / 4 * (q.d3_cm**2 - q.d2_cm**2) / q.years  # cm2/yr, over bark
    plots = q.groupby("est").agg(
        bai=("bai", "median"), d=("d2_cm", "median"), n=("bai", "size"),
        density=("density_per_ha", "first"), date2=("date2", "first"), date3=("date3", "first"),
    )
    plots = plots[plots.index.isin(clim.index)]
    ci = [_plot_climate_index(clim.loc[e], r.date2, r.date3) for e, r in plots.iterrows()]
    plots["precip"], plots["tmin"], plots["g"] = zip(*ci)
    print(f"\nClimate-response test, Q. ilex (Gea-Izquierdo et al. 2013): {len(q)} trees with positive growth "
          f"({n_nonpos} with zero/negative growth excluded), {len(plots)} plots")
    print(f"  plot climate 1990s: precip {plots.precip.min():.0f}-{plots.precip.max():.0f} mm, "
          f"mean daily Tmin {plots.tmin.min():.1f}-{plots.tmin.max():.1f} C; corr(precip, Tmin) = {plots.precip.corr(plots.tmin):.2f}")

    X = np.column_stack([np.ones(len(plots)), np.log(plots.d), np.log(plots.density), np.log(plots.g)])
    y = np.log(plots.bai.values)
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    rng = np.random.default_rng(seed)
    boots = [np.linalg.lstsq(X[i], y[i], rcond=None)[0][3] for i in (rng.integers(0, len(y), len(y)) for _ in range(n_boot))]
    lo, hi = np.percentile(boots, [2.5, 97.5])
    print(f"  log BAI = a + {beta[1]:.2f} log d + {beta[2]:.2f} log density + gamma log g:  "
          f"gamma = {beta[3]:.2f}  (95% plot-bootstrap CI {lo:.2f} .. {hi:.2f}; 1 = published response, 0 = no climate effect)")

    # Nonparametric view: size/density-adjusted growth by Tmin band. Does it peak near 9 C?
    Xs = X[:, :3]
    resid = y - Xs @ np.linalg.lstsq(Xs, y, rcond=None)[0]
    plots["adj"] = np.exp(resid)
    band = pd.cut(plots.tmin, [4, 7, 8, 9, 10, 10.5, 11, 12])
    tab = plots.groupby(band, observed=True).agg(plots=("adj", "size"), adj_growth=("adj", "median"), model_g=("g", "median"))
    tab["model_rel"] = tab.model_g / tab.model_g.max()
    tab["obs_rel"] = tab.adj_growth / tab.adj_growth.max()
    print("  size/density-adjusted growth by mean daily Tmin band (relative to the band maximum):")
    print(tab[["plots", "obs_rel", "model_rel"]].round(2).to_string())


PLOT_CLIMATE_BY_PROVINCE = (os.environ.get("IFN_DATA_ROOT", os.path.expanduser("~/ifn-data")) + "/ifn_validation/plot_climate_p{code}.csv")  # fetch_ifn_plot_climate.py


def national_plot_table(species_code: int = 45, min_d: float = 10.0, clim_template: str = PLOT_CLIMATE_BY_PROVINCE) -> pd.DataFrame:
    """Plot-level table across all provinces with pairs + climate: plot-median BAI of the species'
    positive-growth trees, median dbh, IFN2 density, interval climate (P, Tmin, Gea index), elevation."""
    rows = []
    for path in sorted(glob.glob(clim_template.format(code="*"))):
        code = path[-6:-4]
        clim = pd.read_csv(path).set_index("est")
        q = load_pairs(code)
        q = q[(q.species_code == species_code) & (q.d2_cm >= min_d) & (q.d3_cm > q.d2_cm)].copy()
        if q.empty:
            continue
        q["bai"] = math.pi / 4 * (q.d3_cm**2 - q.d2_cm**2) / q.years
        pl = q.groupby("est").agg(
            bai=("bai", "median"), d=("d2_cm", "median"), n=("bai", "size"),
            density=("density_per_ha", "first"), date2=("date2", "first"), date3=("date3", "first"),
        )
        pl = pl[pl.index.isin(clim.index)]
        if pl.empty:
            continue
        ci = [_plot_climate_index(clim.loc[e], r.date2, r.date3) for e, r in pl.iterrows()]
        pl["precip"], pl["tmin"], pl["g"] = zip(*ci)
        pl["elev"] = clim.loc[pl.index, "elev_m"].values
        pl["lat"] = clim.loc[pl.index, "lat"].values
        pl["province"] = code
        rows.append(pl.reset_index())
    return pd.concat(rows, ignore_index=True)


def _ols_boot(X, y, cols, groups=None, n_boot=1000, seed=0):
    """OLS coefficients for the last `cols` columns, with bootstrap CIs resampling plots."""
    b = np.linalg.lstsq(X, y, rcond=None)[0]
    rng = np.random.default_rng(seed)
    boots = np.array([np.linalg.lstsq(X[i], y[i], rcond=None)[0][-cols:] for i in (rng.integers(0, len(y), len(y)) for _ in range(n_boot))])
    lo, hi = np.percentile(boots, [2.5, 97.5], axis=0)
    return b[-cols:], lo, hi


def test_climate_response_national(species_code: int = 45, clim_template: str = PLOT_CLIMATE_BY_PROVINCE) -> pd.DataFrame:
    from spain_climate_modifier import ILEX_P1_MM, ILEX_P2, ILEX_TOPT_C, ILEX_TSCALE_C

    pl = national_plot_table(species_code, clim_template=clim_template)
    pl = pl[(pl.density > 0) & np.isfinite(pl.tmin)].copy()
    print(f"\nNATIONAL climate-response test, species {species_code}, climate {os.path.basename(clim_template)}: {len(pl)} plots, {pl.province.nunique()} provinces, "
          f"{int(pl.n.sum())} trees; mean daily Tmin {pl.tmin.min():.1f}-{pl.tmin.max():.1f} C, precip {pl.precip.min():.0f}-{pl.precip.max():.0f} mm")
    y = np.log(pl.bai.values)
    prov = pd.get_dummies(pl.province, drop_first=True).values.astype(float)
    base = [np.ones(len(pl)), np.log(pl.d), np.log(pl.density)]
    g1 = np.log(1 / (1 + (pl.precip / ILEX_P1_MM) ** ILEX_P2))
    g2 = -np.abs(pl.tmin - ILEX_TOPT_C) / ILEX_TSCALE_C
    specs = [
        ("gamma on log(g1*g2)", [np.log(pl.g)], ["gamma"]),
        ("Gea terms separately", [g1, g2], ["g1 precip", "g2 temp"]),
        ("Gea terms + elevation/100", [pl.elev / 100, g1, g2], ["elev", "g1 precip", "g2 temp"]),
        ("linear Tmin, precip/100", [pl.tmin, pl.precip / 100], ["Tmin per C", "precip per 100mm"]),
    ]
    for fe_name, fe in [("no province FE", []), ("province fixed effects", [prov])]:
        for name, cols, labels in specs:
            X = np.column_stack(base + cols) if not fe else np.column_stack(base + fe + cols)
            b, lo, hi = _ols_boot(X, y, len(cols))
            print(f"  [{fe_name}] {name}: " + ", ".join(f"{l} {bb:+.3f} ({a:+.3f}..{h:+.3f})" for l, bb, a, h in zip(labels, b, lo, hi)))
    # Nonparametric: growth adjusted for size, density and province, by Tmin band, vs the model's g2 shape.
    Xs = np.column_stack(base + [prov])
    resid = y - Xs @ np.linalg.lstsq(Xs, y, rcond=None)[0]
    pl["adj"] = np.exp(resid)
    band = pd.cut(pl.tmin, [-5, 2, 4, 6, 7, 8, 9, 10, 11, 12, 16])
    tab = pl.groupby(band, observed=True).agg(plots=("adj", "size"), obs=("adj", "median"),
                                               model=("tmin", lambda t: float(np.median(np.exp(-np.abs(t - ILEX_TOPT_C) / ILEX_TSCALE_C)))))
    tab["obs_rel"], tab["model_rel"] = tab.obs / tab.obs.max(), tab.model / tab.model.max()
    print("  growth adjusted for size, density, province, by mean daily Tmin band (relative to band max):")
    print(tab[["plots", "obs_rel", "model_rel"]].round(2).to_string())
    return pl


def _pine_height_coastal(t):  # coastal ecoregion; anchor Ho = 14 m at 20 y (a mid coastal site, illustrative)
    return pinus_pinaster_project(14.0, 20.0, t, ecoregion_coastal=True)


PONTEVEDRA_CURVES = {
    "P. pinaster height (interior, Table I site)": (26, ("h2_m", "h3_m"), _pine_height, ""),
    "P. pinaster height (coastal, Ho 14 m at 20 y)": (26, ("h2_m", "h3_m"), _pine_height_coastal, ""),
    "P. pinaster diameter (interior height curve + h-d inverse, Atlantic)": (
        26, ("d2_cm", "d3_cm"),
        lambda t: diameter_from_height(max(_pine_height(t), 1.31), hd_params("pinus_pinaster", "Atlantic")), "",
    ),
}


if __name__ == "__main__":
    pairs = load_pairs("caceres")
    print(f"Caceres strict pairs: {len(pairs)}; by species: "
          f"{ {SPECIES.get(k, k): v for k, v in pairs.species_code.value_counts().head(6).items()} }")
    report_t1(test_quercus_ilex_t1(pairs))
    test_default_curves(pairs, "Caceres", CURVES)
    test_climate_response(pairs)
    # Pontevedra: IFN2 (1989) -> IFN3 (1998) links are weaker; strict pairs only. Pine is the
    # Sylvestris species there, and the pine curves were fit in Galicia.
    test_default_curves(load_pairs("pontevedra"), "Pontevedra", PONTEVEDRA_CURVES)
