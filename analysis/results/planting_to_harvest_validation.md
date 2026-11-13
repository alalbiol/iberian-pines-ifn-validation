# Planting-to-harvest validation audit — four additional Iberian pines

Checked 28 October 2026 against the local IFN3/IFN4 analyses and the frozen
forecast interface. This audit concerns *P. sylvestris*, *P. pinea*,
*P. halepensis*, and *P. nigra*. The Galician *P. pinaster* age chain has its
own separate validation and limitations.

## Decision

**A planting-to-harvest forecast is not validated for any of these four
species.** The available tests support narrower claims: age-based biomass
checks on independently held-out stands, young-afforestation equation transfer
checks, survivor-tree height and diameter increments, and inventory-anchored
IFN3→IFN4 biomass forecasts. None observes and predicts the full path from
planted seedlings through establishment, stand development, management, and
harvest using only information available at planting.

The shared `pine_stand_forecast.py` interface therefore returns
`conditional_stand_snapshots` for these species. It requires the **current**
number of trees at least 7.5 cm DBH, or uses a weak fit-province default for
Scots and Aleppo pine. Density stays fixed unless a user selects the
illustrative Scots-pine PS1/PS2 caps described below. It does not forecast
establishment, mortality, ingrowth, thinning decisions or post-thinning
growth. A clearfell input removes standing stock at the supplied
age; it does not predict harvest timing or harvested wood.

**Scots-pine management source found 10 November 2026.** Galicia's official
[PS1 and PS2 silvicultural models](https://mediorural.xunta.gal/es/temas/forestal/ordenacion/adhesion-a-ms-y-rbfp-orientativos)
do prescribe thinning for *P. sylvestris*. In particular,
[PS2](https://mediorural.xunta.gal/sites/default/files/temas/forestal/adhesions-ms/modelos/PS2_cas.pdf)
specifies removing 550–600 stems/ha at age 25 and 300 stems/ha at age 35,
with a 45–50-year final cut on suitable high-quality sites. Its stated
objective is short-rotation technical timber; it is not a universal Scots-pine
schedule. This source corrects the earlier suggestion that no published
Scots-pine plan was available. The demo now offers None, PS1 and PS2 as
**conditional management scenarios**. PS1 caps future density at 825 trees/ha
at age 25, representing the midpoint of the published early 750–900 range;
later repeated thinnings lack numerical tree-count targets and are not drawn.
PS2 caps future density at 800 at age 25 and 500 at age 35. For a stand with
1,375 trees/ha before the first cut, those caps correspond to the source's
approximate removals of 575 then 300 trees/ha. A sparser stand loses less or
none, and a measured current stand does not have past cuts applied again.
These options supply no establishment, post-thinning growth response, or
validation. A forester should still choose which plan and stand are suitable.

## Independent evidence available

All biomass RMSE values below are Mg/ha. The first comparison is the adult
age-based model against an IFN4 stand at a known age. The second is the
published young-afforestation equation transferred to planted IFN4 stands;
its measurement target and allometry differ from the source study. The
inventory-anchored test starts with an observed IFN3 tree list and therefore
cannot validate a forecast made at planting.

| Species | Adult age-based held-out plots; default / fit-median RMSE | Young planted held-out plots; default / fit-median RMSE | Young linked height cohorts | Uncut inventory-anchored held-out plots; forecast / no-change RMSE |
|---|---:|---:|---:|---:|
| *P. sylvestris* | 86; 53.02 / 58.50 | 55; 44.91 / 49.20 | 33 | 327; 24.31 / 40.30 |
| *P. pinea* | 9; 33.92 / 23.24 | 11; 20.68 / 22.47 | 0 | 167; 23.95 / 36.71 |
| *P. halepensis* | 26; 27.31 / 36.24 | 11; 29.90 / 30.49 | 3 | 433; 18.61 / 29.35 |
| *P. nigra* | 14; 105.14 / 82.44 | 6; 74.98 / 57.19 | 7 | 287; 23.78 / 40.93 |

The age-based defaults for stone and black pine lose to the simple
fit-province median. The young default for black pine also loses to that
baseline; Aleppo pine has a 17.78 Mg/ha positive bias on its 11 held-out
young stands. The published young-height projection loses to no growth for
black pine (3.45 versus 3.02 m RMSE); stone pine has no eligible linked young
cohort. These outcomes rule out treating the young and adult components as a
validated continuous chain.

The local IFN4 stand extract contains 200, 39, 53, and 26 directly aged,
planted, pine-dominated, measured stands for Scots, stone, Aleppo, and black
pine, respectively. Their observed total-age ranges are 6–60, 9–50, 14–50,
and 28–57 years. These are cross-sectional observations of trees that have
already entered the IFN measured class. They do not observe planting stock,
early survival, or a harvest outcome for the same stand.

## Requirements for a full validation

An independent planted-stand dataset must record planting date and initial
stems, site and provenance, repeated counts and sizes including trees below
7.5 cm DBH, management dates and removals, and a measured pre-harvest or
harvest biomass outcome. The model must predict survival, entry into the
measured class, growth, thinning response, and harvest from inputs known at
planting. Freeze that model before evaluating separate stands or regions;
report errors by species, age, and management with sample sizes and uncertainty.
Until those inputs and components exist, the correct scientific claim is the
narrower inventory-anchored and age-specific performance reported above.

## Reproduction

Using `~/miniconda3/envs/sim_case1/bin/python` from the repository root,
`evaluate_young_other_pines.py` and `evaluate_other_pine_stands.py` both
reproduced their saved JSON files exactly on 28 October 2026. The other
figures are from `other_pines_inventory_projection.json` and
`other_pines_height_remeasurement.json`; those full database evaluations were
not rerun for this audit. The source data and scripts are described in
`other_pines_README.md`.
