"""IFN4 age, species dominance, and live-tree screen for four other pines.

No model fitting. Uses the same main-species, even-aged, age-source, and
80% basal-area conditions as the Pinaster development work. Historical
cut/dead records retained in PCMayores are excluded from current trees.
"""

from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "models"))

from ifn4_age_validation import IFN4_DIR, _plot_key  # noqa: E402
from ifn_remeasurement_validation import _mdb, _num, trees_per_ha  # noqa: E402
from other_pine_tree_growth import BIOMASS  # noqa: E402
from stand_structure_calibration import dominant_mean  # noqa: E402

SPECIES = {21: "Pinus sylvestris", 23: "Pinus pinea",
           24: "Pinus halepensis", 25: "Pinus nigra"}
BAD_CURRENT_ORDER = (444, 666, 777, 888, 999)


def normalized_plot_key(frame: pd.DataFrame, province_code: int | None = None) -> pd.Series:
    raw = _plot_key(frame, province_code)
    parts = raw.str.extract(r"^0*([0-9]+)_(.*)$")
    if parts.isna().any().any():
        raise RuntimeError("Unparseable IFN plot key")
    return parts[0] + "_" + parts[1]


def audit_database(path: Path) -> pd.DataFrame:
    ages = _mdb(str(path), "PCEspParc")
    ages["plot_id"] = normalized_plot_key(ages)
    for field in ("PosEsp", "Especie", "FPMasa", "Edad", "FInfor", "Fiabil",
                  "OrgMasa1", "OrgMasa2"):
        ages[field] = _num(ages[field])
    ages = ages[(ages.PosEsp == 1) & ages.Especie.isin(SPECIES) &
                ages.FPMasa.isin((1, 2)) & ages.Edad.notna() & (ages.Edad > 0)].copy()
    if ages.empty:
        return pd.DataFrame()
    ages["direct_age"] = ages.FInfor.isin((5, 8))
    ages["age_screen"] = (ages.FInfor != 1) & (ages.direct_age | ages.Fiabil.isin((3, 4)))
    ages = ages[ages.age_screen].copy()
    if ages.empty:
        return pd.DataFrame()
    if ages.duplicated("plot_id").any():
        raise RuntimeError(f"Duplicate main-species age rows: {path}")
    trees = _mdb(str(path), "PCMayores")
    trees["plot_id"] = normalized_plot_key(trees)
    trees = trees[trees.plot_id.isin(set(ages.plot_id))].copy()
    trees["species"] = _num(trees.Especie)
    trees["d_cm"] = (_num(trees.Dn1) + _num(trees.Dn2)) / 20.0
    trees["h_m"] = _num(trees.Ht)
    trees["ord4"] = _num(trees.OrdenIf4)
    trees = trees[(trees.d_cm >= 7.5) & (trees.h_m > 1.3) &
                  (trees.ord4 > 0) & ~trees.ord4.isin(BAD_CURRENT_ORDER)].copy()
    if trees.empty:
        return pd.DataFrame()
    trees["w"] = trees_per_ha(trees.d_cm.to_numpy())
    trees["ba"] = trees.w * np.pi * (trees.d_cm / 200.0)**2
    total = trees.groupby("plot_id").agg(total_ba=("ba", "sum"),
                                          all_trees=("ba", "size"))
    main = trees.merge(ages[["plot_id", "Especie"]].rename(columns={"Especie": "main_species"}),
                       on="plot_id", validate="many_to_one")
    main = main[main.species == main.main_species]
    main["agb_mg_ha"] = [w*BIOMASS[SPECIES[int(sp)].lower().replace(" ", "_")](d, h).aboveground_kg/1000.0
                          for sp, d, h, w in zip(main.species, main.d_cm, main.h_m, main.w)]
    h0 = main.groupby("plot_id").apply(
        lambda g: dominant_mean(g.h_m.to_numpy(), g.d_cm.to_numpy(), g.w.to_numpy()),
        include_groups=False).rename("h0_m")
    mean_h = main.groupby("plot_id").apply(
        lambda g: float(np.average(g.h_m.to_numpy(), weights=g.w.to_numpy())),
        include_groups=False).rename("mean_large_h_m")
    d0 = main.groupby("plot_id").apply(
        lambda g: dominant_mean(g.d_cm.to_numpy(), g.d_cm.to_numpy(), g.w.to_numpy()),
        include_groups=False).rename("d0_cm")
    main = main.groupby("plot_id").agg(main_ba=("ba", "sum"),
                                        pine_trees=("ba", "size"),
                                        pine_stems_ha=("w", "sum"),
                                        agb_mg_ha=("agb_mg_ha", "sum"))
    main = main.join(h0).join(d0).join(mean_h)
    out = ages.merge(total, on="plot_id", how="left", validate="one_to_one")
    out = out.merge(main, on="plot_id", how="left", validate="one_to_one")
    out["ba_share"] = out.main_ba/out.total_ba
    out["pine_dominated_measured"] = out.ba_share.ge(0.8) & out.pine_trees.ge(1)
    out["species"] = out.Especie.map(SPECIES)
    out["database"] = path.stem
    out["province_code"] = out.plot_id.str.split("_").str[0]
    out["planted"] = (out.OrgMasa1 == 2) & (out.OrgMasa2 == 2)
    return out[["plot_id", "province_code", "database", "species", "Edad",
                "FInfor", "Fiabil", "direct_age", "planted", "total_ba",
                "main_ba", "ba_share", "pine_trees", "pine_stems_ha",
                "h0_m", "mean_large_h_m", "d0_cm", "agb_mg_ha",
                "pine_dominated_measured"]]


def load_data() -> pd.DataFrame:
    frames = [audit_database(path) for path in sorted(Path(IFN4_DIR).glob("*.accdb"))]
    data = pd.concat([frame for frame in frames if not frame.empty], ignore_index=True)
    if data.duplicated("plot_id").any():
        raise RuntimeError("IFN4 plot keys occur in multiple database files")
    return data


def audit(data: pd.DataFrame) -> dict:
    result = {"design": "IFN4 main pine species, even-aged, direct or reliable age, live >=7.5cm trees, >=80% pine basal area; all locally available databases",
              "databases": int(data.database.nunique()), "species": {}}
    for species, group in data.groupby("species"):
        eligible = group[group.pine_dominated_measured].copy()
        direct = eligible[eligible.direct_age]
        result["species"][species] = {
            "age_screen": len(group),
            "pine_dominated_measured": len(eligible),
            "direct_age": len(direct),
            "direct_age_planted": int(direct.planted.sum()),
            "direct_age_under_25": int((direct.Edad < 25).sum()),
            "direct_age_25plus": int((direct.Edad >= 25).sum()),
            "direct_age_age_range": ([float(direct.Edad.min()), float(direct.Edad.max())]
                                     if len(direct) else None),
            "direct_age_by_province": direct.province_code.value_counts().to_dict(),
            "direct_age_source_counts": {str(int(k)): int(v) for k, v in
                                         direct.FInfor.value_counts().items()},
        }
    return result


if __name__ == "__main__":
    import json
    data = load_data()
    data.to_csv(Path(__file__).parent / "results/other_pines_ifn4_stands.csv", index=False)
    print(json.dumps(audit(data), indent=2, ensure_ascii=False))
