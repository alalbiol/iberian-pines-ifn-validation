# Other Spanish pines: development status

The [planting-to-harvest validation audit](planting_to_harvest_validation.md)
records which independent tests are complete and why the four additional
species do not yet have validated planting-to-harvest forecasts.

Species: *Pinus sylvestris* (IFN 21), *P. pinea* (23), *P. halepensis*
(24), and *P. nigra* (25). Codes follow the [MITECO IFN species list](https://www.miteco.gob.es/es/biodiversidad/temas/inventarios-nacionales/mapa-forestal-espana/especies_arboreas.html).
The input data are the local IFN3/IFN4 Access databases and strict
tree-pair caches in `/media/HD/Data/ifn_validation/`.

Run in order from the repository root with the `sim_case1` environment:

```bash
python analysis/audit_other_pines.py
python analysis/fit_other_pine_tree_growth.py
python analysis/evaluate_other_pine_cohort_agb.py
python analysis/evaluate_other_pine_height_remeasurement.py
python analysis/evaluate_other_pine_stands.py
python analysis/evaluate_other_pine_inventory_projection.py
python analysis/evaluate_young_other_pines.py
```

The stand screen requires an even-aged main species, directly measured
whorl or bore age for the stricter tests, at least 80% target-pine basal
area, and live measured trees at least 7.5 cm DBH. There are 207 direct-age
*P. sylvestris* stands, 62 *P. halepensis*, 42 *P. pinea*, and 26
*P. nigra* across the available provinces. Current-tree order and height
filters exclude historical cut/dead records retained in `PCMayores`.
Province prefixes are normalized before matching IFN3 and IFN4. Plots
with duplicated old-tree order numbers are excluded from the turnover
analysis.

The height shapes come from published regional work:
[Galician Scots pine](https://doi.org/10.1051/forest:2005006),
[regional stone pine](https://doi.org/10.1051/forest:2003017),
[Aleppo pine plantations](https://doi.org/10.5424/720), and
[Spanish black pine](https://doi.org/10.1016/j.foreco.2008.06.024).
Four corresponding aboveground tree allometries are taken from
[Ruiz-Peinado et al. (2011), Table 3](https://doi.org/10.5424/fs/2011201-11643).
IFN reference AGB uses these same allometries, so reported errors do not
validate the allometric equations against harvested biomass.

## Held-out results

The strict surviving-tree increment fits use two designated provinces for
each species and all other available provinces for testing. They condition
on a tree surviving both cycles. The reusable code is
`other_pine_tree_growth.py`. RMSE units are cm for DBH and m for height.

| Species | Held-out trees | DBH model / no growth | Height model / no growth |
|---|---:|---:|---:|
| *P. sylvestris* | 14,132 | 2.78 / 4.67 | 2.01 / 2.91 |
| *P. pinea* | 6,106 | 3.64 / 6.10 | 1.45 / 2.25¹ |
| *P. halepensis* | 9,365 | 2.99 / 4.90 | 1.59 / 2.25 |
| *P. nigra* | 8,867 | 2.48 / 4.11 | 1.67 / 2.47 |

¹ The fitted constant height increment transferred better for *P. pinea*
than the larger size/density model (1.45 versus 1.56 m RMSE). That choice
was made after inspecting these test errors, so the selected API needs a
fresh confirmatory test. The fitted *P. nigra* DBH model barely improves
on its constant-rate baseline (2.482 versus 2.485 cm).

On known surviving tree cohorts, projected IFN-expanded AGB reduces
plot-level RMSE compared with holding initial cohort AGB fixed: 7.56 versus
14.00 Mg/ha for *P. sylvestris*, 7.92 versus 11.54 for *P. pinea*,
6.79 versus 12.72 for *P. halepensis*, and 7.09 versus 13.63 for
*P. nigra*. This calculation knows the future survivor set and therefore
is not a whole-stand forecast.

The prospective inventory-anchored projection instead starts with the
complete IFN3 target-pine tree list, forecasts natural survival with a
fit-province constant hazard, includes expected changes in IFN concentric
plot expansion, and adds the fit-province mean annual new-tree AGB.
It tests only plots without recorded target-pine cuts. The entry term is
an empirical constant, not a published species-specific ingrowth model.
The reusable API is `other_pine_inventory_forecast.py`, which requires
the measured target-pine tree list, **all-species** initial plot density,
and elapsed years. Its per-tree removal fractions are supplied management
inputs. The tested intervals span about 5–21 years by species; the API
rejects extrapolation beyond 25 years.

| Species | Fit / held-out plots | Total AGB forecast / initial AGB RMSE (Mg/ha) | Paired-plot 95% interval for RMSE gain |
|---|---:|---:|---:|
| *P. sylvestris* | 182 / 327 | 24.31 / 40.30 | 14.48–17.41 |
| *P. pinea* | 71 / 167 | 23.95 / 36.71 | 7.84–17.66 |
| *P. halepensis* | 283 / 433 | 18.61 / 29.35 | 8.67–12.65 |
| *P. nigra* | 283 / 287 | 23.78 / 40.93 | 14.59–19.57 |

The paired plot bootstrap resamples held-out plots 2,000 times. The
predicted entry term gives a clear RMSE gain over the projected old cohort
for *P. sylvestris*, *P. halepensis*, and *P. nigra*. For *P. pinea*, the
95% gain interval includes zero (−1.72–3.60 Mg/ha). All four forecasts
remain imprecise at individual plots. They apply to measured existing
trees over an IFN remeasurement interval, not a planting-to-harvest path.

In a retrospective cut-tree check, actual IFN4 cut identities are supplied
as perfect removal fractions. Prior-tree AGB RMSE then improves over
initial AGB after subtracting those trees at their original size:
*P. sylvestris* 16.18 versus 21.09 Mg/ha (200 held-out plots),
*P. pinea* 10.19 versus 14.90 (114), *P. halepensis* 13.34 versus
17.91 (96), and *P. nigra* 18.14 versus 26.11 (55). This conditional
check does not predict cut identities, timing, or post-cut growth response.

The published dominant-height equations are implemented in
`other_pine_dominant_height.py`. A strict linked top-100/ha survivor-cohort
check finds 59 Galician *P. sylvestris* plots (1.87 m height RMSE), seven
*P. halepensis* (1.06 m), three *P. nigra* (1.39 m), and five
*P. pinea* (1.69 m). The latter three samples are too small to establish
reliable trajectory transfer. The cohort height is not a repeated full
stand dominant height.

The age-based `other_pine_stand.py` is a **research prototype** using
fit-province median site index and density plus fitted H→D and
mean/dominant biomass ratio. Its default AGB RMSE is 53.02 Mg/ha on
86 held-out Galician *P. sylvestris* stands, 33.92 on nine *P. pinea*,
27.31 on 26 *P. halepensis*, and 105.14 on 14 *P. nigra*. The *P. pinea*
and *P. nigra* defaults lose to a fit-province median-AGB baseline;
the *P. nigra* median predicted/measured AGB ratio is 2.69. Supplying
measured height and density reduces *P. nigra* RMSE to 24.76 Mg/ha.
These defaults must not be used as reliable planting-to-harvest forecasts.

`pine_stand_forecast.py` now exposes all five pines through one forecast
interface. The four additional species use published adult height curves,
their frozen IFN-fitted stand parameters, Ruiz-Peinado aboveground/root
allometry, and species carbon fractions. It accepts planting age, coordinates,
area, current stems ≥7.5 cm/ha, dominant height/diameter, and the *P. pinea*
source region. For these species, current large-tree density is fixed except
for optional Scots-pine PS1/PS2 density caps; planting density is rejected because establishment and
entry into the ≥7.5 cm class are not modeled. *P. pinea* and *P. nigra*
require measured density and site height/index. Coordinates do not supply
a site index. Results carry `conditional_stand_snapshots` status and cannot
be interpreted as validated planting-to-harvest or net-sequestration
trajectories. The live `use_case_1_spain/demo/` Simulate tab calls this API.

Remaining work: validate establishment survival and growth below the
measured-tree threshold, species-specific ingrowth, actual thinning timing and
post-thinning growth response, age-source conventions, and stronger site/density defaults.
Galicia publishes [Scots-pine PS1 and PS2 silvicultural plans](https://mediorural.xunta.gal/es/temas/forestal/ordenacion/adhesion-a-ms-y-rbfp-orientativos).
The demo now applies their quantified density steps as optional conditional
scenarios; the [planting-to-harvest audit](planting_to_harvest_validation.md)
explains the omitted steps and why this is not a validated complete forecast.
The source range for the black-pine height paper begins at age 30 and the
stone-pine paper at age 20; younger-age projections are excluded here.

## Published young-afforestation transfer

The species-specific mean-height and age/density/SI30 biomass equations
from [Calama et al. (2024)](https://doi.org/10.3390/f15050827) are now in
`young_other_pines.py`. Its source uses years since planting, mean height
including small trees, and biomass from root-collar-based allometry. The
IFN4 test assumes total inventory age minus one nursery year and uses
mean height of recorded trees at least 7.5 cm. It is an operational
transfer check with a target mismatch, especially at young ages.

| Species | Fit / held-out planted stands | Default-SI30 biomass RMSE | Measured-SI30 biomass RMSE | Fit-median-AGB baseline |
|---|---:|---:|---:|---:|
| *P. sylvestris* | 19 / 55 | 44.91 | 35.98 | 49.20 |
| *P. pinea* | 7 / 11 | 20.68 | 11.31 | 22.47 |
| *P. halepensis* | 20 / 11 | 29.90 | 30.77 | 30.49 |
| *P. nigra* | 10 / 6 | 74.98 | 47.85 | 57.19 |

All RMSE values are Mg/ha. The Aleppo-pine equation overpredicts by
17.78 Mg/ha on its 11 held-out plots. On strict linked planted survivor
cohorts within the source age range, published mean-height RMSE is
2.37 m for *P. sylvestris* (33 plots), 1.37 m for *P. halepensis*
(three plots), and 3.45 m for *P. nigra* (seven plots). The last loses
to no growth (3.02 m). No eligible *P. pinea* young remeasurement
cohort was found. Thus these equations are implemented, but they do not
yet validate a continuous planting-to-harvest trajectory for all four.
