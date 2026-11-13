# Mediterranean Pinus pinaster first pass

Run from the repository root with the `sim_case1` environment and locally
available MITECO IFN4 databases and IFN3→IFN4 strict-pair caches:

```bash
python analysis/audit_mediterranean_pinaster.py
python analysis/evaluate_mediterranean_pinaster.py
python analysis/evaluate_mediterranean_remeasurement.py
python analysis/fit_mediterranean_diameter.py
python analysis/fit_mediterranean_structure.py
python analysis/evaluate_mediterranean_stand.py
python analysis/evaluate_mediterranean_ifn3.py
python analysis/evaluate_young_afforestation.py
python analysis/fit_mediterranean_density.py
python analysis/fit_mediterranean_height_increment.py
python analysis/evaluate_mediterranean_ingrowth.py
python analysis/audit_mediterranean_turnover.py
python analysis/fit_mediterranean_natural_loss.py
python analysis/evaluate_mediterranean_expansion.py
python analysis/evaluate_mediterranean_uncut_density.py
python analysis/fit_mediterranean_tree_survival.py
python analysis/evaluate_mediterranean_tree_biomass.py
python analysis/evaluate_mediterranean_thinning_scenario.py
```

The corresponding JSON files here are the captured outputs. Inputs are
`/media/HD/Data/ifn4/*.accdb`, `/media/HD/Data/ifn3/pXX/Ifn3pXX.accdb` and
`/media/HD/Data/ifn_validation/pairs34_p*.csv`. The latter were built by
`models/ifn3_ifn4_pairs.py` using
species, bearing and distance checks. Source age/method/origin definitions
are in the local MITECO IFN4 field database manual, annexes 8–10.

The height model is Bravo-Oviedo, del Río and Montero (2007),
*Forest Ecology and Management* 247:107–119,
[doi:10.1016/j.foreco.2007.04.034](https://doi.org/10.1016/j.foreco.2007.04.034),
equations 16–17 and Table 3. The Galician comparator is Álvarez-González
et al. (2005), and is used only to quantify regional transfer. See
`../../paper_plan.md` for the interpretation and limitations.

The mature stand prototype is
`models/mediterranean_pinaster_stand.py`.
It combines the published height projection with IFN4-fitted dominant
diameter, mean/dominant biomass ratio, and Ruiz-Peinado tree allometry.
The IFN4 fit used Madrid and Salamanca. On 16 other-province known-age
plots, default inputs give AGB RMSE 57.97 Mg/ha; measured dominant height
and density give 21.01 Mg/ha. The density term was adjusted after viewing
these IFN4 residuals, making this a development check.

`evaluate_mediterranean_ifn3.py` excludes all IFN3 plot keys present in the
available IFN4 databases, as well as the Madrid/Salamanca fit provinces.
Its 85 remaining known-age, pine-dominated stands give default AGB RMSE
39.04 Mg/ha and measured-input RMSE 17.58 Mg/ha. They are ages 25–40 y
and mostly in Ciudad Real and Toledo, so this does not establish full
25–180 y transfer or a planting-to-harvest forecast. The
[MITECO IFN3 field database guide](https://portal-miteco-stage.adobecqms.net/content/dam/miteco/es/biodiversidad/servicios/banco-datos-naturaleza/documentador_bdcampo_ifn3_tcm30-282240.pdf)
confirms `FInfor=5` means whorl count, tree height is in metres, and
diameters are in millimetres. It does not resolve the exact origin of
recorded stand age. Both inventories'
reference AGB uses the same allometry as the model.

A separate young-afforestation candidate from Calama et al. (2024),
*Forests* 15:827, [doi:10.3390/f15050827](https://doi.org/10.3390/f15050827),
is implemented in `young_pinaster_afforestation.py`. On five held-out
documented planted remeasurement cohorts, its mean-height forecast has
RMSE 6.02 m and +5.92 m bias, worse than the no-growth baseline (2.48 m).
Its biomass equation has 54.77 Mg/ha RMSE on 11 held-out planted stands
when using the fit-province site index. These transfer failures prevent
using it as a validated young-stand trajectory.

Two further repeated-inventory experiments use Madrid/Salamanca for fitting
and other central provinces for testing. Whole-plot density change is
measured in 640 pine-dominated plots with at least two strict pine tree
links and both cycles available. The extractor requires current live-tree
height and order codes, excluding historical cut/dead records retained in
`PCMayores`. On 534 test plots, keeping density fixed gives 331 stems/ha
RMSE, a fit-province constant annual loss rate gives 256 stems/ha, and a
density-dependent annual rate gives 263 stems/ha.
About 27% of test plots gain large pine stems. This observational net
change includes ingrowth, mortality and management, so neither fitted rate
is accepted as a reliable stand-density forecast. The strict-link and
pine-dominance screens also exclude stands that lost every pine tree.
On 11,689 strictly linked test pine trees, an exploratory height-size-density increment model
has 2.21 m final-height RMSE; a constant annual increment has 1.89 m.
These survivor-tree results do not validate a young-stand height curve.

The [published five-year ingrowth equation](https://doi.org/10.5424/srf/2008173-01039)
is implemented in `mediterranean_pinaster_ingrowth.py`. On 382 other-province
plots within its basal-area and quadratic-mean-diameter range, comparing
new-tree IFN4 basal area divided by elapsed years and multiplied by five
gives 1.40 m²/ha RMSE for its probability-weighted prediction versus 2.01
m²/ha for zero ingrowth. This interval scaling is approximate and its event
probability transfers poorly (Brier 0.255 versus 0.052 for a fit-province
constant probability). Converting basal area to new stems/ha with a fitted
factor gives 152 stems/ha RMSE and −54 stems/ha bias. Thus the source
equation is a useful ingrowth amount candidate, not an accepted stem-count
trajectory.

IFN4 annex 39 records `OrdenIf3=000` for new sample trees and separate
`OrdenIf4` codes for cut (000), not found without cutting evidence (444),
and dead (888) prior trees. The retrospective turnover audit covers 640
plots: 47% have a recorded cut; mean expanded pine loss is 136 stems/ha
from cuts and 41 from recorded deaths. On 276 uncut test plots, a fitted
natural-loss rate reduces loss RMSE from 143 to 110 stems/ha. Combining
that loss rate with predicted ingrowth does not improve final density.
Surviving trees also change concentric-subplot expansion weights as their
DBH crosses thresholds: their mean shift is −88 stems/ha, and the full
observed-outcome balance reconstructs IFN4 density within 4.6 stems/ha
RMSE. A diameter model with fitted residual uncertainty predicts that
shift at 76 stems/ha RMSE on 532 other-province plots, conditional on
knowing which trees survived. This sampling effect must be represented in
any future individual-tree validation; future cuts and survival remain
unknown inputs or processes.

The prospective inventory-anchored tree-list API is
`models/mediterranean_pinaster_tree_forecast.py`.
It takes measured pine DBH and height, a projection interval, and optional
per-tree thinning fractions supplied by the user. Optional initial basal
area and quadratic mean diameter activate the published ingrowth component
within its source range. It reports expected IFN-expanded surviving-cohort
stems and pine AGB. Thinning timing and its effect on growth are unestimated.

```python
from mediterranean_pinaster_tree_forecast import InventoryTree, project_inventory

trees = [InventoryTree(20.0, 12.0), InventoryTree(30.0, 18.0)]
forecast = project_inventory(trees, years=15.0)
thinned = project_inventory([InventoryTree(20.0, 12.0, 0.25),
                             InventoryTree(30.0, 18.0, 0.25)], years=15.0)
```

Each input tree is an observed IFN sample tree, not one tree per hectare.
The optional third field is its supplied removal fraction. The current
prototype uses initial IFN density in both scenarios when computing DBH
growth, pending a tested post-thinning response.

`fit_mediterranean_tree_survival.py` fits 0.00809 annual natural-death hazard
on 44 uncut Madrid/Salamanca plots. On 276 uncut other-province plots, it
predicts old-tree loss with 110 stems/ha RMSE, versus 143 for zero loss.
`evaluate_mediterranean_tree_biomass.py` combines this with frozen DBH
increment, a fitted constant height increment, expected IFN expansion
weights, and published ingrowth BA. On the same 276 held-out plots,
surviving prior-tree AGB RMSE is 11.97 Mg/ha versus 18.89 for keeping initial
AGB fixed. On 179 plots within the ingrowth paper's BA/QMD range, total
pine AGB RMSE is 24.80 Mg/ha versus 42.62 for keeping initial AGB fixed.
The paired-plot bootstrap 95% interval for that RMSE reduction is
14.20–21.27 Mg/ha. Adding predicted ingrowth to the projected old cohort
reduces RMSE from 27.08 to 24.80 Mg/ha, but its bootstrap interval for the
gain includes zero (−0.17–4.51 Mg/ha). Province biases vary: +14.99 Mg/ha
in Cuenca and −15.25 in Cáceres. Allometry is shared by predictions and
reference AGB, so these errors test tree growth and aggregation, not direct
biomass measurement. Fifteen-to-eighteen-year interval scaling of the
five-year ingrowth source remains approximate. This does not validate a
planting-to-harvest or managed-stand projection.

A further retrospective management check supplies the IFN4 identities of
cut trees as exact per-tree removal fractions. On 232 other-province plots
with recorded cuts, projected surviving prior-tree AGB has 12.11 Mg/ha RMSE,
versus 19.17 Mg/ha for initial AGB after subtracting cut trees at their
initial size. This tests the forecast conditional on a perfectly known
cut-tree scenario. It does not predict which trees will be cut, when the
cut occurs, or how remaining trees respond afterward.

For young planted stands, seven Madrid/Salamanca linked cohorts give a
median mean-height increment of 0.1955 m/year. That simple measured-stand
bridge has 0.98 m RMSE on five held-out planted cohorts, versus 2.48 m
for no growth and 6.02 m for the published Calama et al. (2024) curve.
The held-out cohorts start roughly 10–24 years after planting; this small
test does not establish seedling-to-tree growth or a full young-stand path.
Only one of the 19 documented-planted, direct-age IFN4 stands in the
published young-biomass source range is younger than ten years; two are
younger than fifteen. The local inventory subset is therefore especially
thin at establishment ages.
