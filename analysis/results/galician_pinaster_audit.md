# Galician maritime pine held-out result: reproduction and age audit

The current code produces the 144-plot Ourense/Pontevedra IFN4 result
(corrected; see below). Run from the repository root with the `sim_case1` environment:

```bash
python use_case_1_simulator/agb_forecast_analysis/ifn4_pine_calibration.py
python use_case_1_simulator/agb_forecast_analysis/ifn4_age_validation.py
python analysis/audit_galician_pinaster.py
```

The final command writes `galician_pinaster_audit.json`. The fit provinces
are A Coruña and Lugo; neither held-out province supplies calibration
parameters. The held-out result retains plots aged at least 10 years with
pine-dominated measured trees and IFN-derived AGB above 1 Mg/ha.

| Held-out subset | Plots | Default model / IFN AGB, median (IQR) | Measured-height-and-density model / IFN AGB, median (IQR) |
|---|---:|---:|---:|
| All accepted ages | 144 | 0.836 (0.491–1.464) | 0.974 (0.711–1.386) |
| Direct age source only | 109 | 0.835 (0.486–1.501) | 0.971 (0.707–1.377) |

**Correction, 11 November 2026.** Stem number and observed AGB now count
trees ≥7.5 cm only (`large_tree_view` in `ifn4_age_validation.py`). The
earlier version added counted small trees (2.5–7.5 cm, `PCRegenera`) to both.
The stand chain is calibrated on trees ≥7.5 cm, so a plot with a sparse
overstorey over thousands of saplings received the overstorey tree's biomass
multiplied by every sapling (one plot: 159 trees/ha ≥7.5 cm plus ~4,300
saplings, forecast 1,332 Mg/ha against 60 observed). The other pine cases
already used trees ≥7.5 cm. With the correction, measured-input RMSE falls
from 154.1 to 65.3 Mg/ha; default RMSE is 45.6 Mg/ha either way. Two plots
without trees ≥7.5 cm drop out (146 → 144). Previous figures were
0.833 (0.493–1.448) and 1.082 (0.770–1.556) on 146 plots.

The 144 plots comprise 86 in Ourense and 58 in Pontevedra; 82 are recorded
as planted. IFN4 `FInfor` is 8 (boring) on 95, 5 (whorl count) on 14,
and 7 (other sources) on 35. Source 7 enters only with `Fiabil` 3 or 4
(medium or high reported reliability). Definitions are in annexes 8–9 of
the local MITECO IFN4 field database guide. Restricting to the 109 direct
sources barely changes the aggregate result. No historical/cut/dead rows
with eligible diameter and height were found in the held-out pine plots or
the two fit-province tree tables.

The guide defines `Edad` as even-aged stand age but does **not** specify
whether it starts at germination, planting or the tree-coring height. The
model treats it as total age and uses a one-year nursery lag to translate
the planting-age input. As a sensitivity scenario, shifting all ages by
−2 or +2 years changes the default median model/IFN ratio to 0.739 or
0.909, respectively. The measured-height-and-density ratio stays 0.974
because the curve is anchored at each plot's observed height. These shifts
are scenarios, not corrections supported by individual plot records.

The remaining limitation is therefore **age-origin uncertainty**, not a
missing computation or a missing age-source count. A definitive origin
would require additional plot-level field documentation or an independent
planting-date record. The paper can report the result now with this stated
limitation and its sensitivity analysis.
