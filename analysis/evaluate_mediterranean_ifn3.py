"""Exploratory cross-cycle test of the IFN4-fitted Mediterranean pine stand chain.

IFN3 provinces are read through the same tree expansion and allometry used
for IFN4. Exclude Madrid/Salamanca (fit provinces) and plot identifiers seen
in any available IFN4 age record before scoring. IFN3 `Edad` time origin
remains ambiguous; source documentation confirms `Ht` in metres and
`FInfor=5` as whorl count.
"""

import os  # Configurable external inventory paths.

from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "models"))

from audit_pine_ifn4 import CENTRAL_DATABASES, PINASTER_CODE  # noqa: E402
from evaluate_mediterranean_pinaster import FIT_PROVINCES, metrics as height_metrics  # noqa: E402
from evaluate_mediterranean_stand import metrics  # noqa: E402
from ifn4_age_validation import IFN4_DIR, _plot_key, large_tree_view, load_aged_plots, measured_stands  # noqa: E402
from ifn_remeasurement_validation import _mdb  # noqa: E402
from mediterranean_pinaster import REFERENCE_AGE, SOURCE_STAND_AGE_RANGE, project_height  # noqa: E402
from mediterranean_pinaster_stand import DEFAULT_LARGE_STEMS_HA, DEFAULT_SITE_INDEX_70_M, estimate_stand_agb  # noqa: E402

IFN3_DIR = Path((os.environ.get("IFN_DATA_ROOT", os.path.expanduser("~/ifn-data")) + "/ifn3"))
PROVINCES = {
    2: "Albacete", 6: "Badajoz", 10: "Cáceres", 13: "Ciudad Real",
    16: "Cuenca", 19: "Guadalajara", 28: "Madrid", 37: "Salamanca",
    45: "Toledo",
}


def load_ifn3() -> pd.DataFrame:
    frames = []
    for code, region in PROVINCES.items():
        if region in FIT_PROVINCES:
            continue
        path = IFN3_DIR / f"p{code:02d}" / f"Ifn3p{code:02d}.accdb"
        if not path.exists():
            continue
        aged = load_aged_plots(region, path.stem, database_path=str(path), province_code=code)
        # IFN3 defines 5 as whorl count; unlike IFN4, its field document
        # does not define age-source code 8 as boring.
        aged = aged[(aged.Especie == PINASTER_CODE) & (aged.FInfor == 5) &
                    aged.Edad.between(*SOURCE_STAND_AGE_RANGE)].copy()
        if aged.empty:
            continue
        # Observed AGB from trees >= 7.5 cm only, matching the input density and the other cases.
        stands = large_tree_view(measured_stands(path.stem, aged, database_path=str(path), province_code=code))
        if stands.empty:
            continue
        stands = stands.merge(aged[["plot_id", "Edad", "Fiabil"]], on="plot_id", validate="one_to_one")
        frames.append(stands[(stands.h0_m > 1.3) & np.isfinite(stands.h0_m)].copy())
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def evaluate() -> dict:
    df = load_ifn3()
    if df.empty:
        raise RuntimeError("No eligible IFN3 stands")
    # Exclude every IFN4 plot key, including plots without pine dominance or
    # reliable age, so this cross-cycle test is spatially disjoint from all
    # development analyses in these provinces.
    seen = set()
    for region, database in CENTRAL_DATABASES.items():
        inventory = _mdb(str(Path(IFN4_DIR) / f"{database}.accdb"), "PCEspParc")
        seen.update(_plot_key(inventory))
    overlap = df.plot_id.isin(seen)
    fresh = df[~overlap].copy()
    if fresh.empty:
        raise RuntimeError("All eligible IFN3 plots overlap the explored IFN4 plots")
    fresh["default_height"] = [project_height(DEFAULT_SITE_INDEX_70_M, REFERENCE_AGE, age)
                               for age in fresh.Edad]
    for name, use_h, use_n in (("defaults", False, False),
                               ("measured_height_and_density", True, True)):
        fresh[name] = [estimate_stand_agb(age,
                       dominant_height_m=h if use_h else None,
                       large_stems_per_ha=n if use_n else None).aboveground_biomass_mg_ha
                       for age, h, n in zip(fresh.Edad, fresh.h0_m, fresh["trees_per_ha_ge7.5"])]
    return {
        "design": "IFN3 whorl-count age-source 5 pine-dominated stands outside fit provinces and with plot keys absent from all available IFN4 age records",
        "ifn3_eligible_n": len(df), "excluded_ifn4_plot_key_overlap_n": int(overlap.sum()),
        "fresh_n": len(fresh), "fresh_by_province": fresh.region.value_counts().to_dict(),
        "fresh_total_age_range": [float(fresh.Edad.min()), float(fresh.Edad.max())],
        "fresh_height": height_metrics(fresh, "default_height"),
        "default_large_stems_per_ha": DEFAULT_LARGE_STEMS_HA,
        "fresh_large_stems_per_ha_median": round(float(fresh["trees_per_ha_ge7.5"].median()), 3),
        "fresh_results": {name: metrics(fresh, name) for name in ("defaults", "measured_height_and_density")},
        "fresh_province_results": {province: {name: metrics(g, name) for name in
                                   ("defaults", "measured_height_and_density")}
                                   for province, g in fresh.groupby("region")},
    }


if __name__ == "__main__":
    import json
    print(json.dumps(evaluate(), indent=2, ensure_ascii=False))
