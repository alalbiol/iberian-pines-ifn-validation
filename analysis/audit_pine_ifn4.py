"""Inventory-age and stand-eligibility screen for the next pine case (no fitting).

Run with the sim_case1 Python environment. Reads the local MITECO IFN4 Access
databases and applies the filters from ifn4_age_validation.py. The first table
is an age screen; the second also requires pine dominance and measured trees.
"""

import os  # Configurable external inventory paths.

from io import StringIO
from pathlib import Path
import subprocess
import sys

import pandas as pd


IFN4_DIR = Path((os.environ.get("IFN_DATA_ROOT", os.path.expanduser("~/ifn-data")) + "/ifn4"))
MDB_EXPORT = Path(os.environ.get("MDB_EXPORT", "mdb-export"))
PINASTER_CODE = 26
MODEL_DIR = Path(__file__).resolve().parents[1] / "models"
sys.path.insert(0, str(MODEL_DIR))
from ifn4_age_validation import (  # noqa: E402
    PLANTING_LAG_YEARS, load_aged_plots, measured_stands,
)

CENTRAL_DATABASES = {
    "Albacete": "Ifn4_Albacete", "Ciudad Real": "Ifn4_CiudadReal",
    "Extremadura": "Ifn4_Extremadura", "Guadalajara": "Ifn4_Guadalajara",
    "Toledo": "Ifn4_Toledo", "Ávila": "Ifn4p05-avila",
    "Madrid": "Ifn4p28-madrid", "Salamanca": "Ifn4p37-salamanca",
    "Zamora": "Ifn4p49-zamora", "Cuenca": "Ifn4p_Cuenca",
}


def age_counts(path: Path) -> dict[str, int | str]:
    result = subprocess.run(
        [str(MDB_EXPORT), str(path), "PCEspParc"],
        capture_output=True, text=True, check=True,
    )
    rows = pd.read_csv(StringIO(result.stdout), low_memory=False)
    for field in ("PosEsp", "Especie", "FPMasa", "Edad", "FInfor", "Fiabil", "OrgMasa1", "OrgMasa2"):
        rows[field] = pd.to_numeric(rows[field], errors="coerce")
    pine = rows[(rows.PosEsp == 1) & (rows.Especie == PINASTER_CODE)]
    aged = pine[pine.FPMasa.isin((1, 2)) & pine.Edad.notna()]
    reliable = aged[
        (aged.FInfor != 1) & (aged.FInfor.isin((5, 8)) | aged.Fiabil.isin((3, 4)))
    ]
    older = reliable[reliable.Edad - PLANTING_LAG_YEARS >= 10]
    return {
        "database": path.stem,
        "main_species_plots": len(pine),
        "known_age_even_aged": len(aged),
        "age_screen_pass": len(reliable),
        "age_screen_pass_10plus": len(older),
        "age_screen_pass_10plus_planted": int(((older.OrgMasa1 == 2) & (older.OrgMasa2 == 2)).sum()),
    }


def stand_counts(region: str, database: str) -> dict[str, int | str]:
    aged = load_aged_plots(region, database)
    aged = aged[(aged.Especie == PINASTER_CODE) &
                (aged.Edad - PLANTING_LAG_YEARS >= 10)]
    stands = measured_stands(database, aged) if len(aged) else pd.DataFrame()
    return {
        "region": region,
        "age_screen_pass_10plus": len(aged),
        "pine_dominated_measured_10plus": len(stands),
        "of_which_recorded_planted": int(stands.planted.sum()) if len(stands) else 0,
    }


if __name__ == "__main__":
    summary = pd.DataFrame(age_counts(path) for path in sorted(IFN4_DIR.glob("*.accdb")))
    print("AGE SCREEN")
    print(summary.to_csv(index=False), end="")
    print("\nCENTRAL-SPAIN STAND SCREEN")
    stands = pd.DataFrame(stand_counts(region, database)
                          for region, database in CENTRAL_DATABASES.items())
    print(stands.to_csv(index=False), end="")
