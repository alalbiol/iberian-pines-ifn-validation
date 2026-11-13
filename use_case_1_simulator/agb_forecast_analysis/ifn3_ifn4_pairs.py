"""
IFN3 -> IFN4 remeasured tree pairs, the second growth interval (2026-09-27).

Why. `ifn_remeasurement_validation.py` builds IFN2 -> IFN3 pairs (457,825 strict pairs, ~10 y apart,
1986-96 -> 1997-2007). That interval is now 20-30 years old and its trees were measured with 1980s
protocols. IFN4 (2008-2024) re-measures the same permanently-marked plots, and IFN4's `PCMayores`
carries `OrdenIf3`, the tree's IFN3 number, exactly as IFN3 carried `OrdenIf2`. So a second, more
recent interval is available for every province whose IFN4 field database we hold.

What it is for, in order of value tonight:
  * CORK OAK: the stand-level test is stuck at four known-age IFN4 plots and cannot grow (the age field
    exists only for even-aged stands, and oak woodland is almost never classed that way). The GROWTH
    curve, though, is testable without any age at all, and this doubles the data behind it.
  * PINE: an independent check on the large-tree growth that the IFN2 -> IFN3 pairs left ambiguous -
    measured dominant increment is flat with size, but the big trees in that sample are plausibly
    just good sites. A second interval with different plots and a different decade helps.
  * HOLM OAK: more small-tree increments, the B3 question.

Strict pairs only, the same rule as IFN2 -> IFN3: same species, bearing within 10 grads and distance
within 0.5 m of the earlier record. Loose links are dropped rather than repaired.

CAUTION, and it differs from IFN2 -> IFN3: the IFN3 -> IFN4 interval is NOT a uniform ten years. IFN4
field work ran from 2008 (Navarra) to 2024 (Andalucia) and was slowed by the 2010s crisis, so for some
provinces the gap exceeds 15 years while for Galicia it is closer to 10. Every pair therefore carries
its own `years` from the two plot dates, and nothing here assumes a nominal interval.

Run: python ifn3_ifn4_pairs.py [--provinces 06 10 ...]
Writes /media/HD/Data/ifn_validation/pairs34_p<XX>.csv per province plus a link-quality summary.
"""

from __future__ import annotations

import os  # Configurable external inventory paths.

import argparse
import glob
import math
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ifn_remeasurement_validation import _mdb, _num, trees_per_ha  # noqa: E402

IFN3_DIR = (os.environ.get("IFN_DATA_ROOT", os.path.expanduser("~/ifn-data")) + "/ifn3")
IFN4_DIR = (os.environ.get("IFN_DATA_ROOT", os.path.expanduser("~/ifn-data")) + "/ifn4")
OUT_DIR = (os.environ.get("IFN_DATA_ROOT", os.path.expanduser("~/ifn-data")) + "/ifn_validation")
CACHE = os.path.join(OUT_DIR, "pairs34_p{code}.csv")
SPECIES_OF_INTEREST = {26: "pinus_pinaster", 45: "quercus_ilex", 46: "quercus_suber"}


def _find(pattern: str) -> str | None:
    hits = glob.glob(pattern) + glob.glob(pattern.upper()) + glob.glob(pattern.lower())
    return hits[0] if hits else None


def _plot_key(df: pd.DataFrame) -> pd.Series:
    """Estadillo + Cla + Subclase identifies a plot within a province in both cycles."""
    return (df.Estadillo.astype(str).str.strip() + "_" + df.Cla.astype(str).str.strip()
            + df.Subclase.fillna("").astype(str).str.strip())


def _dates(path: str, table: str = "PCParcelas") -> pd.DataFrame | None:
    """plot_key -> survey date. Column names and formats differ between cycles and provinces."""
    try:
        p = _mdb(path, table)
    except Exception:
        return None
    col = next((c for c in ("FechaIni", "FechaFin", "Fecha") if c in p.columns), None)
    if col is None:
        return None
    p = p.copy()
    p["plot_key"] = _plot_key(p)
    d = pd.to_datetime(p[col], format="%m/%d/%y %H:%M:%S", errors="coerce")
    if d.isna().all():
        d = pd.to_datetime(p[col], errors="coerce", dayfirst=True)
    return pd.DataFrame({"plot_key": p.plot_key, "date": d}).dropna().drop_duplicates("plot_key")


def ifn4_databases() -> dict[str, list[str]]:
    """province code -> the IFN4 database files that contain it (regional files hold several)."""
    out: dict[str, list[str]] = {}
    for path in sorted(glob.glob(os.path.join(IFN4_DIR, "*.accdb"))):
        try:
            t = _mdb(path, "PCMayores")
        except Exception as e:
            print(f"  ! cannot read {os.path.basename(path)}: {e}")
            continue
        if "Provincia" not in t.columns:
            print(f"  ! {os.path.basename(path)} has no Provincia column, skipped")
            continue
        for code in sorted(t.Provincia.dropna().astype(str).str.strip().str.zfill(2).unique()):
            out.setdefault(code, []).append(path)
    return out


def build_pairs(code: str, ifn4_paths: list[str]) -> pd.DataFrame:
    """Strict IFN3 -> IFN4 tree pairs for one province."""
    ifn3_path = _find(os.path.join(IFN3_DIR, f"p{code}", f"Ifn3p{code}.accdb"))
    if ifn3_path is None:
        raise FileNotFoundError(f"no IFN3 database for province {code}")
    t3 = _mdb(ifn3_path, "PCMayores")
    t3 = t3.rename(columns={"Distancia": "Distanci"})
    t3["plot_key"] = _plot_key(t3)
    t3["ord3"] = _num(t3.OrdenIf3 if "OrdenIf3" in t3.columns else t3.nArbol)
    t3["d3_cm"] = (_num(t3.Dn1) + _num(t3.Dn2)) / 20.0
    t3["h3_m"] = _num(t3.Ht)
    t3["sp3"] = _num(t3.Especie)
    t3["rumbo3"], t3["dist3"] = _num(t3.Rumbo), _num(t3.Distanci)
    # Cork oak debarking state. A debarked tree's diameter growth is depressed for years afterwards,
    # and a NEW plantation is never-debarked, so the two groups must not be pooled.
    t3["cork3"] = t3.ParEsp.fillna("").astype(str).str.strip() if "ParEsp" in t3.columns else ""
    t3 = t3[t3.ord3 > 0]
    # Stand density at the START of the interval, all species, trees >= 7.5 cm - the control variable.
    dens = (t3[t3.d3_cm >= 7.5].assign(n=lambda d: trees_per_ha(d.d3_cm.to_numpy()))
            .groupby("plot_key").n.sum().rename("density_per_ha"))
    dates3 = _dates(ifn3_path)

    frames = []
    for path in ifn4_paths:
        t4 = _mdb(path, "PCMayores").rename(columns={"Distancia": "Distanci"})
        t4 = t4[t4.Provincia.astype(str).str.strip().str.zfill(2) == code].copy()
        if t4.empty:
            continue
        t4["plot_key"] = _plot_key(t4)
        t4["ord3"] = _num(t4.OrdenIf3)
        t4["d4_cm"] = (_num(t4.Dn1) + _num(t4.Dn2)) / 20.0
        t4["h4_m"] = _num(t4.Ht)
        t4["sp4"] = _num(t4.Especie)
        t4["rumbo4"], t4["dist4"] = _num(t4.Rumbo), _num(t4.Distanci)
        t4["cork4"] = t4.ParEsp.fillna("").astype(str).str.strip() if "ParEsp" in t4.columns else ""
        t4 = t4[t4.ord3 > 0]
        dates4 = _dates(path)

        m = t4.merge(t3, on=["plot_key", "ord3"], suffixes=("_4", "_3"))
        if m.empty:
            continue
        linked = len(m)
        bearing = ((m.rumbo4 - m.rumbo3 + 200) % 400 - 200).abs()  # centesimal degrees
        strict = (m.sp4 == m.sp3) & (bearing <= 10) & ((m.dist4 - m.dist3).abs() <= 0.5)
        m = m[strict].copy()
        m["province"] = code
        m["species_code"] = m.sp3.astype(int)
        if dates3 is not None and dates4 is not None:
            m = (m.merge(dates3.rename(columns={"date": "date3"}), on="plot_key", how="left")
                   .merge(dates4.rename(columns={"date": "date4"}), on="plot_key", how="left"))
            m["years"] = (m.date4 - m.date3).dt.days / 365.25
        else:
            m["date3"] = m["date4"] = pd.NaT
            m["years"] = np.nan
        m = m.join(dens, on="plot_key")
        m.attrs["linked"] = linked
        frames.append(m)

    if not frames:
        return pd.DataFrame()
    out = pd.concat(frames, ignore_index=True)
    cols = ["province", "plot_key", "ord3", "species_code", "d3_cm", "d4_cm", "h3_m", "h4_m",
            "date3", "date4", "years", "density_per_ha", "cork3", "cork4"]
    out = out[cols].dropna(subset=["d3_cm", "d4_cm"]).reset_index(drop=True)
    out.attrs["linked"] = sum(f.attrs.get("linked", 0) for f in frames)
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--provinces", nargs="*", default=None, help="INE codes; default = every province with an IFN4 file")
    args = ap.parse_args()

    os.makedirs(OUT_DIR, exist_ok=True)
    by_prov = ifn4_databases()
    codes = args.provinces or sorted(by_prov)
    print(f"IFN4 databases cover {len(by_prov)} provinces; building {len(codes)}\n")

    summary = []
    for code in codes:
        if code not in by_prov:
            print(f"{code}: no IFN4 data, skipped")
            continue
        try:
            df = build_pairs(code, by_prov[code])
        except Exception as e:
            print(f"{code}: FAILED ({type(e).__name__}: {e})")
            summary.append({"province": code, "strict_pairs": 0, "note": f"failed: {e}"})
            continue
        if df.empty:
            print(f"{code}: no links")
            summary.append({"province": code, "strict_pairs": 0, "note": "no links"})
            continue
        df.to_csv(CACHE.format(code=code), index=False)
        yrs = df.years.dropna()
        counts = {SPECIES_OF_INTEREST[s]: int((df.species_code == s).sum())
                  for s in SPECIES_OF_INTEREST if (df.species_code == s).any()}
        print(f"{code}: {len(df):6d} strict pairs of {df.attrs.get('linked', 0):6d} linked"
              f" | interval {yrs.median():.1f} y (P10-P90 {yrs.quantile(.1):.1f}-{yrs.quantile(.9):.1f})"
              if len(yrs) else f"{code}: {len(df):6d} strict pairs (no usable dates)")
        if counts:
            print(f"      of interest: {counts}")
        summary.append({"province": code, "strict_pairs": len(df), "linked": df.attrs.get("linked", 0),
                        "median_years": float(yrs.median()) if len(yrs) else np.nan, **counts})

    s = pd.DataFrame(summary)
    s.to_csv(os.path.join(OUT_DIR, "link_quality_ifn3_ifn4.csv"), index=False)
    print("\n=== totals ===")
    for sp in SPECIES_OF_INTEREST.values():
        if sp in s.columns:
            print(f"  {sp:16s} {int(s[sp].fillna(0).sum()):7d} strict pairs")
    print(f"  ALL SPECIES      {int(s.strict_pairs.sum()):7d}")
    print(f"\nwrote {os.path.join(OUT_DIR, 'link_quality_ifn3_ifn4.csv')}")
