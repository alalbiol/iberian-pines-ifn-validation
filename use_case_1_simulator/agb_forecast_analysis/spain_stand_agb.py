"""
Stand-level AGB block (Mg/ha): turns the per-tree chain (growth curve -> d, h -> Ruiz-Peinado
biomass, see sylvestris_species_growth_plan.md T9) into stand AGB by multiplying by the number
of trees standing at each age:

    AGB(t) [Mg/ha] = trees_per_ha(t) x dominant-tree AGB(t) [kg] x mean/dominant ratio(N) / 1000

The trees/ha trajectory is built from inputs Sylvestris can supply (planting density, survival
checks, thinning plan, parcel area, site). Every input they do not supply falls back to a
literature default, and the result carries a warning saying what was assumed and what it does
to the number. When defaults are used, a low/high range is computed by re-running the model with
the low/high ends of the ASSUMED inputs only (inputs Sylvestris gave are held fixed).

DEFAULTS AND SOURCES (PDFs in use_case_1_simulator/docs/ where available):
  Pinus pinaster (Galicia) — Xunta de Galicia silvicultural model PP2 (Orden 19/05/2014, annex I
    as amended by Orden 09/02/2021, DOG 39 of 26/02/2021): planting 3.5 x 2.5 to 3.5 x 3 m =
    1,140-952 plants/ha; first thinning at 15-18 y to 650-700 trees/ha; second at 20-25 y to
    400 trees/ha; clearfell at 35-40 y. Early survival: NO sourced value found; 0.90 at 5 y is our
    assumption (it only matters before the first thinning, which resets the count).
  Quercus ilex / Q. suber (Extremadura) — planting density: mean planted density of pure stands in
    Extremadura under EU Reg. 2080/92, 434 (holm oak) and 421 (cork oak) trees/ha (Campos, Ovando
    & Montero, Rev. Esp. Estud. Agrosoc. Pesq. 214, Cuadro 3; mapa_2007_forestaciones_encina_
    alcornoque_2080.pdf). Range 320-625: Montero & Cañellas (1999, Manual de reforestación y cultivo
    de alcornoque; montero_1999_manual_reforestacion_alcornoque.pdf) recommend >= 625 (4 x 4 m),
    or 320-360 on flat land with tree shelters and rigorous replacement of failures.
    Survival: Montero & Cañellas (1999) — losses of 15-25% in the first 8-10 y; in cork oak a
    further 20-25% decline between 10 and 20-25 y (animals, natural mortality and selection
    thinning); then natural mortality "not usually below 1% per year". Holm oak: García-Fayos et
    al. (2024, New Forests) measured 80% survival 25 y after planting in 153 plantations, i.e. no
    second decline, so holm oak keeps 0.80 from 10 to 25 y. The 1%/yr after 25 y is the cork oak
    figure applied to holm oak too (no holm oak value found).
    Cork oak density cap: Montero & Cañellas (1999) Table 4, target trees/ha by stem
    circumference for cork-producing stands (upper end of each range), applied as a ceiling.
  Mean/dominant tree ratio — measured on IFN3 plots, stand_structure_calibration.py (pine in
    Galicia + both oaks in Extremadura pooled: the two curves agree within ~0.05 in every bin with
    data). It is NOT a default Sylvestris could replace; it is always applied, with an info warning.

CARBON: belowground biomass from the same Ruiz-Peinado root equations (per tree, scaled like AGB),
then CO2 = (AGB + BGB) x carbon fraction x 44/12. Carbon fractions are species-specific, from
Montero, Ruiz-Peinado & Muñoz (2005), "Producción de biomasa y fijación de CO2 por los bosques
españoles", INIA Monografías Serie Forestal 13, Table 3 (from CREAF, Ibáñez et al.):
P. pinaster 51.1%, Q. ilex 47.5%, Q. suber 47.2%. PDF: docs/montero_2005_biomasa_co2_bosques_espanoles.pdf.
Only tree biomass: no soil, litter, dead wood or understory, and cork itself is outside the equations.

CLIMATE: `growth_modifier` multiplies the growth curves' increments from `climate_from_age` on
(spain_stand_scenarios.py turns a climate scenario into this number). Survival/mortality do not
respond to climate.
"""

from __future__ import annotations

import math
import os
import sys
from dataclasses import dataclass, field, replace

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from spain_biomass_equations import (  # noqa: E402
    pinus_pinaster_biomass,
    quercus_ilex_biomass,
    quercus_suber_biomass,
)
from spain_growth_curves import (  # noqa: E402
    CORK_OAK_DEFAULT_ANCHOR,
    CORK_THICKNESS_ASSUMED_CM,
    CORK_THICKNESS_PHASE_IN,
    PINASTER_ECOREGION_SOURCE,
    PINASTER_REFERENCE_AGE,
    pine_ecoregion_from_climate,
    PINASTER_FIT_STAND_AGE_MIN,
    PINASTER_FIT_WELL_SAMPLED_AGE_MAX,
    pine_ecoregion_from_location,
    pinus_pinaster_project,
    quercus_ilex_diameter_from_planting,
    quercus_suber_cork_thickness_cm,
    quercus_suber_diameter_default,
    quercus_suber_diameter_over_cork,
    quercus_suber_diameter_project,
    quercus_suber_site_index,
    quercus_suber_from_site_class,
)
from spain_height_diameter import (  # noqa: E402
    HDParams,
    diameter_from_height,
    hd_params_dehesa,
    height_from_diameter,
)

SPECIES = ("pinus_pinaster", "quercus_ilex", "quercus_suber")

CARBON_FRACTION = {"pinus_pinaster": 0.511, "quercus_ilex": 0.475, "quercus_suber": 0.472}  # Montero et al. 2005, Table 3
CO2_PER_C = 44.0 / 12.0

# --- Mean/dominant tree AGB ratio vs trees/ha (IFN3, stand_structure_calibration.py, 2026-09-25) ---
# Columns: median trees/ha of the bin, ratio Q25, median, Q75. 1,674 plots pooled. At <= 100 trees/ha
# every tree is among the 100 thickest, so the ratio is 1 by definition.
# Likely CONSERVATIVE for young even-aged plantations: IFN stands are adult and often uneven-aged, so
# their small-tree tail is heavier than a plantation's.
_RATIO_TABLE = np.array([
    [100.0, 1.00, 1.00, 1.00],
    [124.0, 0.85, 0.91, 0.97],
    [171.0, 0.66, 0.72, 0.79],
    [255.0, 0.51, 0.59, 0.67],
    [350.0, 0.43, 0.51, 0.58],
    [483.0, 0.36, 0.43, 0.51],
    [705.0, 0.29, 0.36, 0.44],
    [951.0, 0.23, 0.33, 0.41],
    [1414.0, 0.28, 0.37, 0.43],
])

# Pine, EVEN-AGED stands only (IFN4 A Coruña + Lugo, FPMasa 1/2, 135 plots; ifn4_pine_calibration.py,
# 2026-09-25). Replaces the IFN3 table for pine: on known-age IFN4 plots the IFN3 table ran 0.1-0.15 low
# (0.35 vs 0.49 measured at 600-1000 trees/ha). The 800-1000 and 1000-1400 bins are merged-in as measured.
# Oaks keep the IFN3 table: too few even-aged oak plots in IFN4 Extremadura to fit one.
_RATIO_TABLE_PINE_EVEN_AGED = np.array([
    [100.0, 1.00, 1.00, 1.00],
    [127.0, 0.79, 0.92, 1.00],
    [269.0, 0.56, 0.63, 0.72],
    [378.0, 0.53, 0.58, 0.66],
    [529.0, 0.42, 0.47, 0.55],
    [711.0, 0.46, 0.50, 0.59],
    [888.0, 0.34, 0.46, 0.55],
    [1105.0, 0.39, 0.41, 0.44],
    [1639.0, 0.36, 0.48, 0.52],
])


def mean_to_dominant_ratio(trees_per_ha: float, quantile: str = "median", species: str | None = None,
                           table: str = "auto") -> float:
    """Mean-tree / dominant-tree AGB ratio.

    table="auto" (default) picks the even-aged IFN4 table for pine and the IFN3 all-structures table for
    anything else, which is the behaviour from 2026-09-25. "even_aged" and "ifn3" force one of them.

    B4c (2026-09-26): the nine known-age oak PLANTATION plots found in B1 measure ratios of 0.68-1.00
    (median ~0.80) at 130-520 trees/ha, against 0.51-0.57 from the IFN3 table the oaks use and 0.47-0.63
    from the even-aged pine table. So the IFN3 table does err low for plantations, as suspected, and even
    the even-aged table is low — but nine plots (three of which hold so few stems that the ratio is near 1
    by construction) cannot support an oak table of its own. The IFN3 table stays the default for oaks and
    "even_aged" is offered as the better-supported alternative."""
    if table == "even_aged":
        table_arr = _RATIO_TABLE_PINE_EVEN_AGED
    elif table == "ifn3":
        table_arr = _RATIO_TABLE
    elif table == "auto":
        table_arr = _RATIO_TABLE_PINE_EVEN_AGED if species == "pinus_pinaster" else _RATIO_TABLE
    else:
        raise ValueError(f'table must be "auto", "ifn3" or "even_aged", got {table!r}')
    table = table_arr
    col = {"low": 1, "median": 2, "high": 3}[quantile]
    if trees_per_ha <= table[0, 0]:
        return 1.0
    return float(np.interp(math.log(trees_per_ha), np.log(table[:, 0]), table[:, col]))


# --- Cork oak density ceiling (Montero & Cañellas 1999, Table 4; upper end of each range) ---
# Stem circumference at 1.3 m (cm, over cork) -> trees/ha. Below 60 cm: no ceiling.
_SUBER_CAP_CIRC_CM = np.array([60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200])
_SUBER_CAP_TREES = np.array([500, 450, 350, 300, 250, 200, 175, 150, 125, 100, 95, 90, 85, 80, 75])


def quercus_suber_density_cap(d_over_cork_cm: float) -> float:
    circ = math.pi * d_over_cork_cm
    if circ < _SUBER_CAP_CIRC_CM[0]:
        return math.inf
    return float(np.interp(circ, _SUBER_CAP_CIRC_CM, _SUBER_CAP_TREES))


# --- Defaults: central value and the (low, high) ends used for the AGB range ---
@dataclass(frozen=True)
class SpeciesDefaults:
    planting_density: tuple[float, float, float]  # low, central, high (trees/ha)
    survival: dict[float, tuple[float, float, float]]  # age -> (low, central, high) fraction of planted alive
    annual_mortality_after_last_check: float
    thinnings: dict[float, float]  # age -> residual trees/ha
    sources: dict[str, str]


DEFAULTS = {
    "pinus_pinaster": SpeciesDefaults(
        planting_density=(952.0, 1046.0, 1140.0),
        survival={5.0: (0.85, 0.90, 0.95)},
        annual_mortality_after_last_check=0.0,
        thinnings={16.5: 675.0, 22.5: 400.0},
        sources={
            "planting_density": "Xunta de Galicia model PP2 (DOG 39, 26/02/2021): 952-1,140 plants/ha",
            "survival": ("OUR ASSUMPTION, still unsourced after a literature search on 2026-09-26: 90% at 5 y. "
                         "For context, known-age IFN4 Galician pine plots hold a median 712 trees/ha at 5-10 y "
                         "(IQR 382-1,114) against a default planting of 1,046, but stems below 2.5 cm are not "
                         "counted and those plots' own planting densities are unknown, so that is a lower bound, "
                         "not a measurement. Awaiting Sylvestris's own establishment data"),
            "management": "Xunta PP2: thin to 650-700 trees/ha at 15-18 y and to 400 at 20-25 y",
        },
    ),
    "quercus_ilex": SpeciesDefaults(
        planting_density=(320.0, 434.0, 625.0),
        survival={10.0: (0.75, 0.80, 0.85), 25.0: (0.75, 0.80, 0.85)},
        annual_mortality_after_last_check=0.01,
        thinnings={},
        sources={
            "planting_density": "mean holm oak planting density in Extremadura under EU Reg. 2080/92 (Campos et al., MAPA); range 320-625 from Montero & Cañellas 1999",
            "survival": "15-25% losses in the first 10 y (Montero & Cañellas 1999); 80% alive at 25 y (García-Fayos et al. 2024); then 1%/yr (cork oak figure, Montero & Cañellas 1999)",
            "management": "no thinning assumed",
        },
    ),
    "quercus_suber": SpeciesDefaults(
        planting_density=(320.0, 421.0, 625.0),
        survival={10.0: (0.75, 0.80, 0.85), 25.0: (0.75 * 0.75, 0.80 * 0.775, 0.85 * 0.80)},
        annual_mortality_after_last_check=0.01,
        thinnings={},
        sources={
            "planting_density": "mean cork oak planting density in Extremadura under EU Reg. 2080/92 (Campos et al., MAPA); range 320-625 from Montero & Cañellas 1999",
            "survival": "15-25% losses in the first 10 y, a further 20-25% by 25 y, then >= 1%/yr (Montero & Cañellas 1999)",
            "management": "density kept at or below Montero & Cañellas (1999) Table 4 target for tree size",
        },
    ),
}

PINE_PP2_CLEARFELL_AGE = 35.0

# --- Named management plans (A5, 2026-09-26) -----------------------------------------------------
# name -> (thinnings {age: residual trees/ha}, source string). "custom" is not listed: it means the
# caller passed `thinnings` explicitly. The default per species is DEFAULT_MANAGEMENT below, and it
# reproduces exactly what the model did before this was introduced.
#
# `galicia_observed` is what IFN4 actually finds standing in Galician pine stands, not a prescription:
# medians of the measured stand (trees >= 2.5 cm) on the known-age pine plots, expressed as residual
# densities at the start of each band. Fitted on the FIT provinces only (A Coruña + Lugo), like every
# other piece of the pine calibration, so the held-out test stays honest: 782 trees/ha at 10-20 y
# (n=41), 707 at 20-30 y (n=25), 423 at 30-45 y (n=23). The held-out provinces come out at 870/764/440,
# i.e. the preset travels. The PP2 model thins far harder (400 trees/ha from 22.5 y), which is the main
# reason the all-defaults model falls short at 15-30 y.
# DECISION PENDING (Alberto): which of the two should be the default? PP2 (the legal model, today's
# default) understates standing carbon at 15-30 y; galicia_observed matches the inventory but describes
# stands that were often simply not thinned on schedule, and says nothing about what Sylvestris will do.
MANAGEMENT_PRESETS: dict[str, dict[str, tuple[dict[float, float], str]]] = {
    "quercus_ilex": {
        "literature": ({}, "no thinning: holm oak afforestations under EU Reg. 2080/92 are not thinned on a "
                           "published schedule (Montero & Cañellas 1999; Campos et al.)"),
        "none": ({}, "no thinning: survival and adult mortality only"),
    },
    "quercus_suber": {
        "literature": ({}, "no thinning schedule; density is instead held at or below the Montero & Cañellas "
                           "(1999) Table 4 target trees/ha for the tree's size (the cork oak density ceiling)"),
        "none": ({}, "no thinning and NO density ceiling: survival and adult mortality only"),
    },
    "pinus_pinaster": {
        "pp2": ({16.5: 675.0, 22.5: 400.0},
                "Xunta de Galicia model PP2 (DOG 39, 26/02/2021): thin to 650-700 trees/ha at 15-18 y and to 400 at 20-25 y"),
        "galicia_observed": ({10.0: 782.0, 20.0: 707.0, 30.0: 423.0},
                             "median measured density of known-age IFN4 pine plots in A Coruña + Lugo (782 trees/ha at "
                             "10-20 y, 707 at 20-30 y, 423 at 30-45 y; n=89) — observed practice, not a prescription"),
        "none": ({}, "no thinning: survival and adult mortality only"),
    },
}
DEFAULT_MANAGEMENT = {"pinus_pinaster": "pp2", "quercus_ilex": "literature", "quercus_suber": "literature"}
# Clearfell ages the sources give, used only for the existing warning, never applied by default.
CLEARFELL_AGE_SOURCE = {"pinus_pinaster": "Xunta PP2 clearfells at 35-40 y"}


def management_plan(species: str, name: str) -> tuple[dict[float, float], str]:
    """(thinnings, source) for a named plan. Raises on an unknown name."""
    presets = MANAGEMENT_PRESETS.get(species, {})
    if name not in presets:
        known = tuple(presets) + ("custom",)
        raise ValueError(f"management for {species} must be one of {known}, got {name!r}")
    thin, source = presets[name]
    return dict(thin), source


# --- Inputs, warnings, result ---
@dataclass
class StandInputs:
    """What Sylvestris can supply. Leave a field as None to use the default (with a warning)."""

    species: str
    planting_density_per_ha: float | None = None
    survival_checks: dict[float, float] | None = None  # age -> fraction of planted trees alive
    annual_mortality_after_last_check: float | None = None
    thinnings: dict[float, float] | None = None  # age -> residual trees/ha; {} means "never thinned"
    # Named management plan (A5): see MANAGEMENT_PRESETS. None uses DEFAULT_MANAGEMENT for the species;
    # passing `thinnings` instead is the "custom" plan. Setting both is an error.
    management: str | None = None
    # Optional final felling: standing biomass drops to 0 from this age on. None = not felled (today's
    # behaviour), which for pine keeps the existing no_final_felling warning.
    clearfell_age: float | None = None
    # Density effect on dominant DIAMETER growth (A6), pine only. None = OFF (today's behaviour).
    # Set to PINE_DENSITY_EXPONENTS["remeasured"] or ["cross_sectional"], or your own exponent.
    density_growth_exponent: float | None = None
    density_reference_trees_per_ha: float | None = None  # None uses PINE_DENSITY_REFERENCE_TREES_PER_HA
    # Mean/dominant AGB ratio table (B4c): "auto" (even-aged for pine, IFN3 for the oaks), "ifn3" or
    # "even_aged". Measured oak PLANTATION ratios sit above both, so "even_aged" is the closer choice
    # for an oak plantation - left opt-in because nine plots cannot settle it.
    mean_dominant_table: str = "auto"
    area_ha: float | None = None
    # Site: pine dominant height (m) at a known age; cork oak site class I-V.
    pine_site_anchor: tuple[float, float] | None = None  # (height_m, age)
    pine_site_preset: str | None = None  # "poor" | "median" | "good"; see PINE_SITE_PRESETS (A4)
    cork_oak_site_class: str | None = None
    # A MEASURED (dominant dbh over cork, age since planting) pair. The 2005 diameter model is
    # base-age invariant, so a real measurement can anchor it directly instead of using the paper's
    # own sample mean at 65 y. This is the most useful thing Sylvestris can give for cork oak.
    cork_oak_diameter_anchor: tuple[float, float] | None = None
    # Holm oak has NO site parameter: the juvenile bridge is one fixed rate and the adult phase is
    # another, neither indexed to site quality. So a measurement cannot set the site the way it can
    # for the other two - all it can do is reset the clock, placing the tree at the age the curve
    # says it reaches that diameter. Honest, and much weaker than a site index.
    holm_oak_diameter_anchor: tuple[float, float] | None = None
    # Age of the planting stock, bridging "years since planting" to the papers' total age (see the
    # AGE CONVENTION block). None uses PLANTING_LAG_YEARS; 0.0 recovers the pre-2026-09-26 behaviour.
    planting_lag_years: float | None = None
    # Pine ecoregion (A3): "coastal" or "interior" (Álvarez-González et al. 2005's dummy I). None
    # derives it from lon_lat/elevation_m if given, else falls back to PINE_DEFAULT_ECOREGION.
    pine_ecoregion: str | None = None
    # (mean annual precipitation mm, mean daily minimum °C) for the parcel. Preferred over lon_lat:
    # it applies Vega et al.'s own climatic criterion instead of our geometric proxy for it.
    # spain_stand_scenarios fills this in automatically when a parcel's climate is loaded.
    ecoregion_climate: tuple[float, float] | None = None
    lon_lat: tuple[float, float] | None = None  # WGS84 (lon, lat) of the parcel
    elevation_m: float | None = None


@dataclass
class AgbWarning:
    code: str
    severity: str  # "high" | "medium" | "info"
    message: str


@dataclass
class StandAgbResult:
    species: str
    ages: list[float]
    trees_per_ha: list[float]
    dominant_tree_agb_kg: list[float]
    mean_to_dominant_ratio: list[float]
    agb_mg_ha: list[float]
    agb_mg_ha_low: list[float]
    agb_mg_ha_high: list[float]
    bgb_mg_ha: list[float]
    co2_t_ha: list[float]  # (AGB + BGB) x carbon fraction x 44/12
    co2_t_ha_low: list[float]
    co2_t_ha_high: list[float]
    total_agb_mg: list[float] | None
    total_co2_t: list[float] | None
    growth_modifier: float
    assumed_inputs: list[str]
    warnings: list[AgbWarning] = field(default_factory=list)


# --- Per-tree dominant trajectory (same defaults as the per-tree chain in the plan doc) ---
# Pine default site: dominant height at 20 y on the interior-ecoregion curve shape, the median of known-age
# IFN4 plots >= 10 y in A Coruña + Lugo (ifn4_pine_calibration.py). The paper's mean stand
# (Álvarez-González et al. 2005 Table I, 10.64 m at 19.3 y -> 11.0 m at 20 y) was the original default;
# measured IFN4 Galician H0 ran 1.32x that curve.
# 2026-09-26 (A2): re-derived as 12.41 m, from 12.76 m. Nothing about the plots changed — the site index is
# now computed at the plot's `Edad` rather than `Edad - 1`, matching the total-age convention the model itself
# uses (AGE CONVENTION above). The 20 y here is the paper's reference age, i.e. a TOTAL age.
PINE_DEFAULT_ANCHOR = (12.41, 20.0)
# Pine DOMINANT height -> DOMINANT diameter (2026-09-25): M4 form fitted on plot-level (D0, H0) of 654
# pine-dominated IFN4 plots in A Coruña + Lugo, RMSE 2.9 m (ifn4_pine_calibration.py). Replaces the nationwide
# mean-tree h-d inverse (Atlantic, plantation), which gave 0.78x the measured D0: it is the diameter of an
# average tree of that height, and dominant trees are thicker. Used both ways (h -> d and back).
PINE_DOMINANT_HD = HDParams(0.2072, 3.4069)
# B6 (2026-09-26): the default cork oak site class moves from III to I for PLANTATIONS.
# Evidence: all four known-age cork oak plantation plots in IFN4 (Extremadura 3, Ciudad Real 1, 17-22 y)
# measure a dominant height of 5.6-7.4 m, and every one of them matches site class I (9 m at 40 y), the
# BEST of the five published classes; class III predicts 2.2 m at 19 y against 6.1 m measured, and that
# single factor (0.36) is essentially the whole cork oak gap (ifn4_oak_diagnosis.py). In fact even class I
# is slightly low (4.9 m at 19 y). Plantations on prepared ground out-grow the natural cork oak forest the
# site classes were built from, and the four plots sit in exactly Sylvestris's region.
# Counter-signal, recorded honestly: IFN2 -> IFN3 remeasurements of ADULT Cáceres dehesa cork oaks showed
# the height model already running 2.1x too fast (0.5 m/decade measured). Those are old, pruned, debarked,
# open-grown trees near their asymptote, so both can be true - but it means class I should not be trusted
# far beyond the ages tested here.
# DECISION PENDING (Alberto): keep class I as the plantation default, or revert to III?
# CORK_OAK_NATURAL_SITE_CLASS restores the old behaviour.
CORK_OAK_DEFAULT_SITE_CLASS = "I"
CORK_OAK_NATURAL_SITE_CLASS = "III"


# --- AGE CONVENTION (A2, 2026-09-26) -------------------------------------------------------------
# Three different "ages" meet in this model. They are NOT the same quantity, and until now the
# difference was buried in one constant inside one validation script.
#
# (1) THE PAPERS' AGE = TOTAL AGE FROM GERMINATION.
#     * Pine, Álvarez-González et al. 2005 (F5014.pdf §2.1, Fig. 5): dominant trees sectioned at the
#       stump, "the age at each section height determined in the laboratory", site index plotted
#       against "total age".
#     * Cork oak, Sánchez-González et al. 2005 (f5074.pdf Table II footnote): "Age: number of rings
#       at stump height (years)" — total age too.
#     * Holm oak juvenile bridge is the EXCEPTION: García-Fayos et al. (2024) measure the tallest
#       tree "25 years after planting", so the bridge is already in years-since-planting and must
#       NOT be shifted. The adult phase is a rate (cm/yr), so no age origin enters it at all.
# (2) IFN4 `Edad` — the field manual (ANEXO 8) only says "edad en las formas coetáneas y/o regulares"
#     and lists how it was obtained; it never states the origin. Worse, the origin DEPENDS on the
#     method (`FInfor`), and we keep three of them:
#       FInfor 5, whorl count  -> whorls are counted up the stem from the base, so ~total age;
#       FInfor 8, bored        -> a core gives the rings at the boring height. The manual does not say
#                                 at what height the age core is taken. If it is breast height (1.3 m,
#                                 as for every diameter in the IFN) this is breast-height age, which
#                                 for P. pinaster in Galicia is roughly 3-4 y less than total age;
#       Fiabil 3/4 by consultation (FInfor 2-4) -> the informant gives the PLANTING year, i.e.
#                                 years since planting.
#     So IFN4 `Edad` is a mixture, and the spread between the conventions (0 to ~4 y) is larger than
#     the nursery lag below. ASSUMPTION, unresolvable from the manual: we treat `Edad` as a total age.
# (3) SYLVESTRIS'S AGE = YEARS SINCE PLANTING. This is the tool's input, because it is the only age a
#     reforestation operator actually knows.
#
# Bridging (2)/(3) to (1) needs the age of the planting stock: total age = years since planting +
# PLANTING_LAG_YEARS. Galician P. pinaster and Spanish oak reforestation both use 1-year nursery
# stock as the norm, so 1 y is the central value. It is OUR assumption, not a sourced figure, and it
# is overridable per run (`StandInputs.planting_lag_years`).
#
# Before this change the model fed years-since-planting straight into curves whose t is a total age,
# i.e. it implicitly used a lag of 0, while `ifn4_age_validation.py` subtracted 1 y from `Edad` — so
# the model and the validation disagreed by exactly one year. Both now share this constant.
PLANTING_LAG_YEARS = 1.0
# Species whose published curve is parameterised on total age, so the lag applies. The holm oak
# bridge is in years since planting already (see above).
_TOTAL_AGE_SPECIES = ("pinus_pinaster", "quercus_suber")


# Pine ecoregion fallback (A3). "interior" keeps the pre-2026-09-26 behaviour, when `dominant_tree`
# hard-coded ecoregion_coastal=False. It is now only a LAST resort: a parcel with its own climate gets
# Vega et al.'s actual climatic criterion (pine_ecoregion_from_climate), and one with coordinates gets
# the geometric stand-in for it. Measured on the Maceira pilot tile (2026-09-27), all 777 parcels come
# out COASTAL on the climatic test, and not marginally — median 1,490 mm/yr against Vega's ~1,500 and a
# mean daily minimum of 10.75 °C against his 9.4-10.8 °C. So for Sylvestris's own site the fallback here
# is simply wrong, and the climate route answers the question this constant used to leave open.
# It still matters what a caller with NO parcel data gets, which is what this constant sets.
# DECIDED (Alberto, 2026-09-27): the no-parcel fallback STAYS "interior" — the conservative reading, about
# 30% below coastal at 40 y. It is not a claim about any real site: a Sylvestris parcel always has its own
# climate and so takes the coastal/interior decision from that, never from this constant.
PINE_DEFAULT_ECOREGION = "interior"
PINE_ECOREGIONS = ("interior", "coastal")

# --- Pine site presets (A4, 2026-09-26) ----------------------------------------------------------
# Dominant height (m) at the paper's reference age of 20 y — i.e. the site index itself. Values are the
# quartiles of the site index of known-age IFN4 plots >= 10 y in the FIT provinces (A Coruña + Lugo,
# n=94; ifn4_pine_calibration.py), so "median" is exactly PINE_DEFAULT_ANCHOR and the presets bracket
# the real spread of Galician plantations rather than a guess.
# NOT split per ecoregion: fitting one anchor per ecoregion was tried (A3) and does not generalise —
# the coastal anchor fitted on A Coruña (11.5 m) is 4.7 m below held-out coastal Pontevedra (16.2 m),
# because province-level site differences swamp the ecoregion label.
# A site index quoted AT the reference age is shape-independent, so these presets are valid for both
# ecoregion curves; only the trajectory away from 20 y differs.
PINE_SITE_PRESETS = {"poor": 10.79, "median": 12.41, "good": 15.71}
PINE_SITE_PRESET_SOURCE = (
    "quartiles of the site index of 94 known-age IFN4 plots >= 10 y in A Coruña + Lugo "
    "(P25 10.79, median 12.41, P75 15.71 m dominant height at 20 y)"
)
# Below this age a measured dominant height SELECTS a preset instead of setting the curve directly.
# Design rule agreed with Alberto (2026-09-26): young measurements are too noisy to fix a site index.
# A1 shows why — measured H0 at 5-8 y runs 1.7x the curve and the implied site index falls monotonically
# with stand age, so a single young measurement projected to age 20 is badly biased upward.
PINE_YOUNG_MEASUREMENT_MAX_AGE = 10.0


@dataclass(frozen=True)
class PineSiteChoice:
    """What a measured (dominant height, age) pair implies about the site."""

    site_index_20_m: float  # continuous site index: dominant height at 20 y on this ecoregion's curve
    preset: str  # nearest named preset
    recommended: str  # "preset" (young, noisy measurement) or "anchor" (old enough to use directly)
    why: str

    def to_inputs(self, species: str = "pinus_pinaster", **kw) -> StandInputs:
        """A StandInputs carrying whichever of the two the rule recommends."""
        if self.recommended == "preset":
            return StandInputs(species, pine_site_preset=self.preset, **kw)
        return StandInputs(species, pine_site_anchor=(self.site_index_20_m, PINASTER_REFERENCE_AGE), **kw)


def pine_site_from_measurement(
    dominant_height_m: float,
    age_since_planting: float,
    ecoregion: str = PINE_DEFAULT_ECOREGION,
    planting_lag_years: float = PLANTING_LAG_YEARS,
) -> PineSiteChoice:
    """Turn a measured dominant height at a known age into a site index and the nearest preset.

    The age is years since planting, as everywhere else in this module; it is converted to the paper's
    total age before projecting to the reference age. Measurements younger than
    PINE_YOUNG_MEASUREMENT_MAX_AGE come back recommending the PRESET, not the continuous index."""
    if dominant_height_m <= 1.3:
        raise ValueError("dominant height must exceed breast height (1.3 m)")
    total = curve_age("pinus_pinaster", age_since_planting, planting_lag_years)
    si = pinus_pinaster_project(dominant_height_m, total, PINASTER_REFERENCE_AGE, ecoregion == "coastal")
    preset = min(PINE_SITE_PRESETS, key=lambda k: abs(PINE_SITE_PRESETS[k] - si))
    if age_since_planting < PINE_YOUNG_MEASUREMENT_MAX_AGE:
        return PineSiteChoice(si, preset, "preset", (
            f"measured at {age_since_planting:.0f} y, below {PINE_YOUNG_MEASUREMENT_MAX_AGE:.0f} y: the implied "
            f"site index ({si:.1f} m) is too noisy to use directly — on IFN4 plots the implied index falls from "
            f"16.3 m at 5-8 y to 13.7 m at 10-15 y — so it only selects the '{preset}' preset "
            f"({PINE_SITE_PRESETS[preset]:.2f} m)"))
    return PineSiteChoice(si, preset, "anchor", (
        f"measured at {age_since_planting:.0f} y: site index {si:.1f} m used directly "
        f"(nearest preset would be '{preset}', {PINE_SITE_PRESETS[preset]:.2f} m)"))


def resolve_pine_ecoregion(inputs: StandInputs) -> tuple[str, str, bool]:
    """(ecoregion, why, is_assumed) for a pine stand: explicit input, else coordinates, else default."""
    if inputs.pine_ecoregion is not None:
        if inputs.pine_ecoregion not in PINE_ECOREGIONS:
            raise ValueError(f"pine_ecoregion must be one of {PINE_ECOREGIONS}, got {inputs.pine_ecoregion!r}")
        return inputs.pine_ecoregion, "supplied", False
    if inputs.ecoregion_climate is not None:
        eco, why = pine_ecoregion_from_climate(*inputs.ecoregion_climate)
        return eco, why, True
    if inputs.lon_lat is not None:
        eco, why = pine_ecoregion_from_location(*inputs.lon_lat, inputs.elevation_m)
        return eco, why, True
    return PINE_DEFAULT_ECOREGION, "no ecoregion and no coordinates given", True


def curve_age(species: str, age_since_planting: float, planting_lag_years: float = PLANTING_LAG_YEARS) -> float:
    """Age to feed the published growth curve, from years since planting (see AGE CONVENTION above)."""
    if species in _TOTAL_AGE_SPECIES:
        return age_since_planting + planting_lag_years
    return age_since_planting


# --- Density effect on dominant DIAMETER growth (A6, 2026-09-26) ----------------------------------
# Until now per-tree growth ignored stand density, so thinning never sped up the remaining trees and
# heavier thinning always showed less carbon. The effect belongs on the DIAMETER, not the height:
# dominant height is conventionally density-independent over a wide range of stockings, which is why
# site index is built from it (Álvarez-González et al. 2005 §1, citing Ortega & Montero), and Castedo
# Dorado et al. (hal-00883904.pdf §2.2.3, §3) tested a thinning term in the h-d relationship and found
# "no benefits", while stating that thinning "has a very positive influence on residual tree diameter
# increment and fewer (if any) effects on tree height growth". So: g(N) multiplies the dominant
# diameter's yearly increment, the dominant height stays on the site curve, and a thinning from below
# causes no jump in dominant diameter — only faster growth afterwards. That is why it must be applied
# path-dependently, year by year, rather than in closed form.
#
# TWO ESTIMATES, BOTH NEGATIVE AND TIGHT, BUT A FACTOR 3.6 APART (a6 in the plan doc):
#  (a) remeasurements: IFN2 -> IFN3 dominant P. pinaster in the four Galician provinces (1,151 trees,
#      627 plots), log(annual diameter increment) ~ log(d) + log(N), plot-bootstrap CIs:
#      d log(increment) / d log(N) = -0.175 (95% CI -0.219 to -0.132), and -0.178 with province fixed
#      effects, so the estimate is not a province artefact. This is the elasticity of the INCREMENT.
#  (b) cross-sectional: a density term in the dominant h-d fit on the same 654 IFN4 plots that produced
#      PINE_DOMINANT_HD, D0 = f(H0) x (N/N_ref)^gamma: gamma = -0.048 (95% CI -0.057 to -0.039), and it
#      lowers the D0 RMSE from 5.69 to 5.29 cm. This is the elasticity of the LEVEL.
# They are not the same quantity — a stand's current N is a poor proxy for the density it grew under, so
# a cumulative level effect is diluted relative to the instantaneous one — but they cannot both be right
# for a model that holds N fixed. THE FEATURE IS THEREFORE OFF BY DEFAULT (see PINE_DENSITY_EXPONENTS).
# DECISION PENDING (Alberto): which exponent, and should it be on?
PINE_DENSITY_EXPONENTS = {"remeasured": -0.175, "cross_sectional": -0.048}
# The density the current calibration implicitly represents. The diameter path is H0 -> D0 through
# PINE_DOMINANT_HD, so it carries the density of the sample that relation was fitted on: the median
# trees/ha of those 654 IFN4 pine plots in A Coruña + Lugo.
PINE_DENSITY_REFERENCE_TREES_PER_HA = 387.0


def pine_density_factor(trees_per_ha: float, exponent: float,
                        reference: float = PINE_DENSITY_REFERENCE_TREES_PER_HA) -> float:
    """g(N)/g(N_ref): how much faster (>1) or slower (<1) the dominant diameter grows at this density."""
    if trees_per_ha <= 0:
        return 1.0
    return (trees_per_ha / reference) ** exponent


def _climate_scaled(base, t: float, modifier: float, from_age: float) -> float:
    """Size at age t when every increment after `from_age` is multiplied by `modifier`."""
    if modifier == 1.0 or t <= from_age:
        return base(t)
    start = base(from_age) if from_age > 0 else 0.0
    return start + modifier * (base(t) - start)


def quercus_ilex_equivalent_age(d_cm: float) -> float:
    """Years since planting at which the holm oak curve reaches this diameter (bisection)."""
    lo, hi = 0.0, 400.0
    for _ in range(120):
        mid = (lo + hi) / 2
        if quercus_ilex_diameter_from_planting(mid) < d_cm:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def _under_cork(d_over_cm: float) -> float:
    """Invert the size-dependent cork thickness: dbh over cork -> dbh under cork."""
    x = d_over_cm
    for _ in range(12):
        x = d_over_cm - 2.0 * quercus_suber_cork_thickness_cm(max(x, 0.0))
    return max(x, 0.1)


def _pine_site_curve(inputs: StandInputs, lag: float):
    """(d_base, hd, coastal): dominant diameter as a function of the PAPER's age, for this stand's site."""
    if inputs.pine_site_anchor is not None:
        h_anchor, t_anchor = inputs.pine_site_anchor
        t_anchor = curve_age("pinus_pinaster", t_anchor, lag)
    elif inputs.pine_site_preset is not None:
        h_anchor, t_anchor = PINE_SITE_PRESETS[inputs.pine_site_preset], PINASTER_REFERENCE_AGE
    else:
        h_anchor, t_anchor = PINE_DEFAULT_ANCHOR
    hd = PINE_DOMINANT_HD
    coastal = resolve_pine_ecoregion(inputs)[0] == "coastal"

    def d_base(t: float) -> float:
        h = pinus_pinaster_project(h_anchor, t_anchor, t, ecoregion_coastal=coastal) if t > 0 else 0.0
        return diameter_from_height(h, hd) if h > 1.3 else 0.0

    return d_base, hd, coastal


def pine_dominant_path(
    ages: list[float],
    inputs: StandInputs,
    growth_modifier: float = 1.0,
    climate_from_age: float = 0.0,
    trees_at_age=None,
) -> list[tuple[float, float, float, float]]:
    """(d_cm, h_m, aboveground kg, root kg) of the dominant pine at each age, stepping year by year.

    Two diameters are carried: one scaled by climate only, which the HEIGHT comes from (dominant height
    is density-independent, so the density factor must never feed back into h), and one scaled by both,
    which the biomass comes from. With the density effect off this reproduces the closed-form
    `_climate_scaled` result exactly, so nothing changes.

    `trees_at_age(age)` supplies N(t) for the density factor; it is ignored when
    `inputs.density_growth_exponent` is None (the default)."""
    lag = PLANTING_LAG_YEARS if inputs.planting_lag_years is None else inputs.planting_lag_years
    d_base, hd, _ = _pine_site_curve(inputs, lag)
    exponent = inputs.density_growth_exponent
    n_ref = inputs.density_reference_trees_per_ha or PINE_DENSITY_REFERENCE_TREES_PER_HA
    m = growth_modifier
    tc0 = curve_age("pinus_pinaster", climate_from_age, lag) if climate_from_age > 0 else 0.0

    # Step on the paper's clock. The grid holds every whole year, every requested age and the age the
    # climate scenario starts at, so the climate result is identical to the closed form.
    wanted = {curve_age("pinus_pinaster", a, lag) for a in ages if a > 0}
    top = max(wanted) if wanted else 0.0
    grid = sorted({0.0, tc0} | wanted | {float(y) for y in range(1, int(math.ceil(top)) + 1) if y <= top})
    d_clim = d_dens = 0.0
    out: dict[float, tuple[float, float]] = {}
    for u, v in zip(grid, grid[1:]):
        inc = d_base(v) - d_base(u)
        f_clim = m if u >= tc0 else 1.0
        d_clim += inc * f_clim
        f_dens = 1.0
        if exponent is not None and trees_at_age is not None:
            # N at the midpoint of the step, on the planting clock (a thinning at age a acts from a on).
            f_dens = pine_density_factor(trees_at_age(max((u + v) / 2.0 - lag, 0.0)), exponent, n_ref)
        d_dens += inc * f_clim * f_dens
        out[v] = (d_dens, d_clim)

    result = []
    for a in ages:
        if a <= 0:
            result.append((0.0, 0.0, 0.0, 0.0))
            continue
        d, dc = out[curve_age("pinus_pinaster", a, lag)]
        if d <= 0:
            result.append((0.0, 0.0, 0.0, 0.0))
            continue
        h = height_from_diameter(dc, hd)
        b = pinus_pinaster_biomass(d, h)
        result.append((d, h, b.aboveground_kg, b.roots_kg))
    return result


def dominant_tree(
    species: str, age: float, inputs: StandInputs, growth_modifier: float = 1.0, climate_from_age: float = 0.0
) -> tuple[float, float, float, float]:
    """(d_cm, h_m, aboveground kg, root kg) of the dominant tree at `age` years after planting.

    The climate modifier scales the DIAMETER increment for all three species, because every source
    it comes from measures radial growth (ring width or basal-area increment), not height. Pine:
    diameter from the height curve via the h-d inverse, scaled, and height back from the forward h-d
    model (so h stays on the h-d curve; identical to the base curve when the modifier is 1). Scaling
    height instead was tried first and is wrong: the h-d inverse is steep, so +20% height gave +100%
    AGB. Cork oak: diameter under cork scaled; height stays on its site curve (height growth
    response not quantified)."""
    if age <= 0:
        return 0.0, 0.0, 0.0, 0.0
    m = growth_modifier
    lag = PLANTING_LAG_YEARS if inputs.planting_lag_years is None else inputs.planting_lag_years
    # Everything below works on the PAPER's clock (see AGE CONVENTION). climate_from_age is a
    # years-since-planting age too, except that 0 means "from planting", which must stay at size 0.
    tc = curve_age(species, age, lag)
    tc0 = curve_age(species, climate_from_age, lag) if climate_from_age > 0 else 0.0
    if species == "pinus_pinaster":
        if inputs.pine_site_anchor is not None:
            # A measured (dominant height, age) pair from Sylvestris: its age is years since planting.
            h_anchor, t_anchor = inputs.pine_site_anchor
            t_anchor = curve_age(species, t_anchor, lag)
        elif inputs.pine_site_preset is not None:
            h_anchor, t_anchor = PINE_SITE_PRESETS[inputs.pine_site_preset], PINASTER_REFERENCE_AGE
        else:
            h_anchor, t_anchor = PINE_DEFAULT_ANCHOR  # already on the paper's reference age (total age)
        hd = PINE_DOMINANT_HD
        coastal = resolve_pine_ecoregion(inputs)[0] == "coastal"

        def d_base(t: float) -> float:
            h = pinus_pinaster_project(h_anchor, t_anchor, t, ecoregion_coastal=coastal) if t > 0 else 0.0
            return diameter_from_height(h, hd) if h > 1.3 else 0.0

        d = _climate_scaled(d_base, tc, m, tc0)
        if d <= 0:
            return 0.0, 0.0, 0.0, 0.0
        h = height_from_diameter(d, hd)
        b = pinus_pinaster_biomass(d, h)
    elif species == "quercus_ilex":
        # No lag: the juvenile bridge's anchor is itself "25 years after planting" (García-Fayos et
        # al. 2024) and the adult phase is a rate, so both are already on the planting clock.
        shift = 0.0
        if inputs.holm_oak_diameter_anchor is not None:
            d_obs, t_obs = inputs.holm_oak_diameter_anchor
            shift = quercus_ilex_equivalent_age(d_obs) - t_obs  # + if the tree is ahead of the curve
        base = (lambda t: quercus_ilex_diameter_from_planting(max(t + shift, 0.0))) if shift else \
            quercus_ilex_diameter_from_planting
        d = _climate_scaled(base, tc, m, tc0)
        h = height_from_diameter(d, hd_params_dehesa("quercus_ilex"))
        b = quercus_ilex_biomass(d, h)
    elif species == "quercus_suber":
        h = quercus_suber_from_site_class(inputs.cork_oak_site_class or CORK_OAK_DEFAULT_SITE_CLASS, tc)
        if inputs.cork_oak_diameter_anchor is not None:
            d_over, t_obs = inputs.cork_oak_diameter_anchor
            d_under = _under_cork(d_over)
            si = quercus_suber_site_index(
                quercus_suber_from_site_class(inputs.cork_oak_site_class or CORK_OAK_DEFAULT_SITE_CLASS,
                                              curve_age(species, t_obs, lag)),
                curve_age(species, t_obs, lag))
            hd_ratio = CORK_OAK_DEFAULT_ANCHOR["hd_ratio"]
            base = lambda t: quercus_suber_diameter_project(  # noqa: E731
                d_under, curve_age(species, t_obs, lag), max(t, 0.5), si, hd_ratio)
        else:
            base = quercus_suber_diameter_default
        d = quercus_suber_diameter_over_cork(_climate_scaled(base, tc, m, tc0))
        b = quercus_suber_biomass(d, h)
    else:
        raise ValueError(f"species must be one of {SPECIES}, got {species!r}")
    return d, h, b.aboveground_kg, b.roots_kg


# --- Trees/ha trajectory ---
def _survival(age: float, checks: dict[float, float], annual_mortality: float) -> float:
    """Piecewise linear from 1.0 at planting through the checks, then compound annual mortality."""
    pts = sorted(checks.items())
    ages = [0.0] + [a for a, _ in pts]
    fracs = [1.0] + [f for _, f in pts]
    if age <= ages[-1]:
        return float(np.interp(age, ages, fracs))
    return fracs[-1] * (1.0 - annual_mortality) ** (age - ages[-1])


def trees_per_ha_trajectory(
    ages: list[float],
    planting_density: float,
    survival_checks: dict[float, float],
    annual_mortality: float,
    thinnings: dict[float, float],
    density_cap=None,
    clearfell_age: float | None = None,
) -> list[float]:
    """Trees/ha at each age. A thinning at age a scales the stand down to its residual density
    (if it is above it); survival keeps acting on what is left. `density_cap(age)` is an optional
    ceiling (cork oak: Montero & Cañellas Table 4 by tree size). From `clearfell_age` on the stand is
    felled, so the count is 0 (standing biomass only: the harvested wood is not tracked)."""

    def n_before_cap(t: float) -> float:
        factor = 1.0  # product of the thinnings already applied
        for a, residual in sorted(thinnings.items()):
            if t < a:
                break
            n_at_a = planting_density * _survival(a, survival_checks, annual_mortality) * factor
            if n_at_a > residual:
                factor *= residual / n_at_a
        return planting_density * _survival(t, survival_checks, annual_mortality) * factor

    out = []
    for t in ages:
        if clearfell_age is not None and t >= clearfell_age:
            out.append(0.0)
            continue
        n = n_before_cap(t)
        if density_cap is not None:
            n = min(n, density_cap(t))
        out.append(n)
    return out


# --- Main entry point ---
def project_stand_agb(
    inputs: StandInputs, ages: list[float], growth_modifier: float = 1.0, climate_from_age: float = 0.0
) -> StandAgbResult:
    sp = inputs.species
    if sp not in SPECIES:
        raise ValueError(f"species must be one of {SPECIES}, got {sp!r}")
    dflt = DEFAULTS[sp]
    warnings: list[AgbWarning] = []
    assumed: list[str] = []

    # Resolve each input: user value (fixed across the range) or default (central + low/high ends).
    if inputs.planting_density_per_ha is not None:
        dens = (inputs.planting_density_per_ha,) * 3
    else:
        dens = (dflt.planting_density[0], dflt.planting_density[1], dflt.planting_density[2])
        assumed.append("planting_density")
        warnings.append(AgbWarning(
            "assumed_planting_density", "high",
            f"Planting density not provided; assumed {dens[1]:.0f} trees/ha ({dflt.sources['planting_density']}). "
            f"Mg/ha scales in proportion to it: the range shown uses {dens[0]:.0f}-{dens[2]:.0f} trees/ha.",
        ))

    if inputs.survival_checks is not None:
        surv = (dict(inputs.survival_checks),) * 3
    else:
        surv = tuple({a: v[i] for a, v in dflt.survival.items()} for i in range(3))
        assumed.append("survival")
        first_age = min(dflt.survival)
        warnings.append(AgbWarning(
            "assumed_survival", "high",
            f"Survival not provided; assumed {surv[1][first_age]:.0%} of planted trees alive at {first_age:.0f} y "
            f"({dflt.sources['survival']}). If losses were higher and dead trees were not replaced, AGB is overestimated."
            + (" This is the one default in the pine chain with NO published source behind it." if sp == "pinus_pinaster" else ""),
        ))

    if inputs.annual_mortality_after_last_check is not None:
        mort = inputs.annual_mortality_after_last_check
    else:
        mort = dflt.annual_mortality_after_last_check
        if mort > 0:
            assumed.append("adult_mortality")
            warnings.append(AgbWarning(
                "assumed_adult_mortality", "medium",
                f"Adult mortality not provided; assumed {mort:.0%}/yr after the last survival check "
                "(Montero & Cañellas 1999 give this as a minimum), so old-stand AGB may be overestimated.",
            ))

    # Management (A5): an explicit `thinnings` dict is the "custom" plan; otherwise a named preset, and
    # if neither is given, the species default — which is the same plan the model used before A5.
    if inputs.thinnings is not None and inputs.management not in (None, "custom"):
        raise ValueError("give either `thinnings` (custom plan) or `management` (named plan), not both")
    if inputs.thinnings is not None:
        thin, plan_name = dict(inputs.thinnings), "custom"
        warnings.append(AgbWarning(
            "management_custom", "info",
            f"Custom management plan: {thin if thin else 'no thinning'}. Beyond ~20 y the standing stock depends "
            "mostly on thinning decisions.",
        ))
    elif inputs.management is not None:
        plan_name = inputs.management
        thin, plan_source = management_plan(sp, plan_name)
        warnings.append(AgbWarning(
            "management_preset", "info" if plan_name != "none" else "medium",
            f"Management plan '{plan_name}': {plan_source}.",
        ))
    else:
        plan_name = DEFAULT_MANAGEMENT[sp]
        thin = dict(dflt.thinnings)
        assumed.append("management")
        warnings.append(AgbWarning(
            "assumed_management", "high" if sp == "pinus_pinaster" else "medium",
            f"No management plan provided; assumed '{plan_name}': {dflt.sources['management']}. "
            "Beyond ~20 y the standing stock depends mostly on thinning decisions."
            + (" On known-age IFN4 plots real Galician stands hold about twice this at 15-30 y; the "
               "'galicia_observed' plan reproduces them." if sp == "pinus_pinaster" else ""),
        ))

    # B7: the cork oak density ceiling belongs to the "literature" plan; "none" and "custom" drop it.
    cap = None
    if sp == "quercus_suber" and inputs.thinnings is None and plan_name == "literature":
        cap = lambda t: quercus_suber_density_cap(dominant_tree(sp, t, inputs, growth_modifier, climate_from_age)[0])  # noqa: E731

    if inputs.clearfell_age is not None:
        warnings.append(AgbWarning(
            "clearfell", "info",
            f"Final felling at {inputs.clearfell_age:.0f} y: standing biomass is 0 from then on. Only STANDING "
            "carbon is modelled, so the harvested wood and any carbon kept in products are not counted.",
        ))
    elif sp == "pinus_pinaster" and max(ages) > PINE_PP2_CLEARFELL_AGE:
        warnings.append(AgbWarning(
            "no_final_felling", "medium",
            f"The Xunta PP2 model clearfells at 35-40 y; ages above {PINE_PP2_CLEARFELL_AGE:.0f} y assume the "
            "stand is NOT felled. If it is harvested, standing AGB drops to ~0 at that age: set clearfell_age.",
        ))

    if sp == "pinus_pinaster":
        eco, why, eco_assumed = resolve_pine_ecoregion(inputs)
        if eco_assumed:
            assumed.append("ecoregion")
            warnings.append(AgbWarning(
                "assumed_ecoregion", "medium",
                f"Pine ecoregion not given; derived as {eco} ({why}). "
                + ("Derived from the parcel's own climate, which is Vega et al.'s actual criterion."
                   if inputs.ecoregion_climate is not None else
                   f"{PINASTER_ECOREGION_SOURCE}. Supply the parcel's climate (or use project_parcel_scenarios "
                   "with a parcel) to apply Vega's climatic criterion instead of this geometric stand-in.")
                + " The choice only acts through the default site — given a measured dominant height it "
                "makes no difference — but on defaults it is worth about 30% of AGB at 40 y.",
            ))

    if sp == "pinus_pinaster" and inputs.pine_site_anchor is None and inputs.pine_site_preset is None:
        assumed.append("site")
        warnings.append(AgbWarning(
            "assumed_site", "medium",
            f"Site quality not provided; pine uses the '{'median'}' preset, the median Galician site of known-age "
            f"IFN4 plots ({PINE_SITE_PRESETS['median']:.2f} m dominant height at 20 y; poor {PINE_SITE_PRESETS['poor']:.2f}, "
            f"good {PINE_SITE_PRESETS['good']:.2f}). A measured dominant height at a known age would replace it: "
            "pine_site_from_measurement() turns one into a site index or the nearest preset.",
        ))
    if sp == "pinus_pinaster" and inputs.pine_site_anchor is not None:
        anchor_age = inputs.pine_site_anchor[1]
        if anchor_age < PINE_YOUNG_MEASUREMENT_MAX_AGE:
            warnings.append(AgbWarning(
                "young_site_anchor", "medium",
                f"The site anchor was measured at {anchor_age:.0f} y, below {PINE_YOUNG_MEASUREMENT_MAX_AGE:.0f} y. "
                "Young dominant heights are too noisy to fix a site index (on IFN4 plots the implied index falls "
                "from 16.3 m at 5-8 y to 13.7 m at 10-15 y, and young plots still carry residual older trees). "
                "Prefer pine_site_from_measurement(), which turns a young measurement into a preset instead.",
            ))
    if sp == "pinus_pinaster" and inputs.pine_site_preset is not None:
        if inputs.pine_site_preset not in PINE_SITE_PRESETS:
            raise ValueError(f"pine_site_preset must be one of {tuple(PINE_SITE_PRESETS)}, got {inputs.pine_site_preset!r}")
        warnings.append(AgbWarning(
            "site_preset", "info",
            f"Site from the '{inputs.pine_site_preset}' preset ({PINE_SITE_PRESETS[inputs.pine_site_preset]:.2f} m "
            f"dominant height at 20 y): {PINE_SITE_PRESET_SOURCE}.",
        ))
    if inputs.cork_oak_site_class is None and sp == "quercus_suber":
        assumed.append("site")
        warnings.append(AgbWarning(
            "assumed_site", "medium",
            f"Site class not provided; cork oak uses class {CORK_OAK_DEFAULT_SITE_CLASS} (5 m at 40 y) and the "
            f"2005 diameter model's average site, with an assumed {CORK_THICKNESS_ASSUMED_CM:.0f} cm cork thickness.",
        ))

    if inputs.planting_density_per_ha is not None:
        lo, hi = dflt.planting_density[0], dflt.planting_density[2]
        if not 0.5 * lo <= inputs.planting_density_per_ha <= 2.0 * hi:
            warnings.append(AgbWarning(
                "density_outside_typical", "medium",
                f"Planting density {inputs.planting_density_per_ha:.0f} trees/ha is far from typical ({lo:.0f}-{hi:.0f}). "
                "Per-tree growth is not adjusted for density, so AGB is less reliable.",
            ))

    if sp == "quercus_suber":
        warnings.append(AgbWarning(
            "cork_thickness", "info",
            f"Cork thickness is phased in with tree size (0 at planting, ~1 cm at {CORK_THICKNESS_PHASE_IN[1][0]:.0f} cm "
            f"under cork, {CORK_THICKNESS_ASSUMED_CM:.0f} cm from {CORK_THICKNESS_PHASE_IN[2][0]:.1f} cm, the size of the "
            "first debarking) rather than the flat 2 cm added to every tree before 2026-09-26, which inflated young "
            "trees' diameter and roots. The two anchors are Montero & Cañellas (1999); the interpolation is ours.",
        ))
    if sp in ("quercus_ilex", "quercus_suber") and inputs.mean_dominant_table == "auto":
        warnings.append(AgbWarning(
            "oak_ratio_table", "medium",
            "The mean/dominant AGB ratio for oaks comes from IFN3 plots of ALL stand structures. On the nine "
            "known-age oak PLANTATION plots found in IFN4 the measured ratio is 0.68-1.00 (median ~0.80) at "
            "130-520 trees/ha against 0.51-0.57 from that table, so stand AGB here is probably 1.4-1.6x too low "
            'for a plantation. Set mean_dominant_table="even_aged" for the closer (still low) alternative.',
        ))
    warnings.append(AgbWarning(
        "mean_tree_correction", "info",
        "Growth curves describe dominant trees; stand AGB applies a mean/dominant tree ratio measured on even-aged "
        "IFN4 Galician pine plots." if sp == "pinus_pinaster" else
        "Growth curves describe dominant trees; stand AGB applies a mean/dominant tree ratio measured on IFN3 plots "
        "of all stand structures. Even-aged plantations are more uniform (for pine the even-aged ratio ran 0.1-0.15 "
        "higher), so this likely errs low.",
    ))

    if sp in _TOTAL_AGE_SPECIES and inputs.planting_lag_years is None:
        assumed.append("planting_lag")
        warnings.append(AgbWarning(
            "assumed_planting_lag", "info",
            f"Ages are read as years since planting, but the {sp.replace('_', ' ')} growth curve is "
            f"parameterised on total age from germination, so {PLANTING_LAG_YEARS:.0f} y of nursery stock is "
            "added (our assumption; standard Spanish reforestation stock is 1 y). Set planting_lag_years to "
            "change it, or 0 to treat the input as a total age. Worth about one year of growth: at 40 y pine "
            "gains ~4 Mg/ha by it.",
        ))

    if sp == "pinus_pinaster":
        # A1 (2026-09-26): the height curve's own fitting data. Álvarez-González et al. (2005) sampled
        # even-aged STANDS of 8-50 y and state that above 30 y "the sample is not representative (less
        # than 30 trees)". Both ends matter here: the tool is routinely run from planting to 100 y.
        lag = PLANTING_LAG_YEARS if inputs.planting_lag_years is None else inputs.planting_lag_years
        curve_ages = [curve_age(sp, t, lag) for t in ages if t > 0]
        if curve_ages and min(curve_ages) < PINASTER_FIT_STAND_AGE_MIN:
            warnings.append(AgbWarning(
                "pine_curve_young_extrapolation", "medium",
                f"Ages below {PINASTER_FIT_STAND_AGE_MIN:.0f} y are outside the height curve's stand-age range "
                f"(Álvarez-González et al. 2005 sampled stands of 8-50 y). Its stem-analysis pairs do reach "
                "younger tree ages, so the shape is informed, but no stand that young was measured. On known-age "
                "IFN4 Galician plots the measured dominant height at 5-8 y ran 1.7x this curve (7.7 vs 4.2 m), "
                "part of which is residual older trees, so young AGB here is uncertain in both directions.",
            ))
        if curve_ages and max(curve_ages) > PINASTER_FIT_WELL_SAMPLED_AGE_MAX:
            warnings.append(AgbWarning(
                "pine_curve_old_extrapolation", "medium",
                f"Ages above {PINASTER_FIT_WELL_SAMPLED_AGE_MAX:.0f} y rest on fewer than 30 sampled trees "
                "(Álvarez-González et al. 2005, §3.2) and the fitting data stop at 50 y. IFN4 known-age plots "
                "older than 30 y also imply a lower site index (11.4 m at 30-45 y, 8.8 m at 45-80 y against "
                "13-15 m at 10-30 y), so the default site is optimistic for old stands.",
            ))

    # The dominant-tree trajectory. For pine it is stepped year by year so the density effect (A6) can
    # act on each year's increment; with the effect off it equals the closed form exactly. It depends on
    # N(t) when the effect is on, so it has to be rebuilt inside each range run.
    def tree_path(planting_density: float, survival: dict[float, float]):
        if sp != "pinus_pinaster":
            return [dominant_tree(sp, t, inputs, growth_modifier, climate_from_age) for t in ages]
        n_at = None
        if inputs.density_growth_exponent is not None:
            def n_at(t: float, _d=planting_density, _s=survival) -> float:
                return trees_per_ha_trajectory([t], _d, _s, mort, thin, cap, inputs.clearfell_age)[0]
        return pine_dominant_path(ages, inputs, growth_modifier, climate_from_age, n_at)

    tree = tree_path(dens[1], surv[1])
    if inputs.density_growth_exponent is not None and sp == "pinus_pinaster":
        warnings.append(AgbWarning(
            "density_growth_effect", "medium",
            f"Dominant diameter growth is scaled by (N/{inputs.density_reference_trees_per_ha or PINE_DENSITY_REFERENCE_TREES_PER_HA:.0f})"
            f"^{inputs.density_growth_exponent:+.3f} year by year, so thinning speeds up the trees left standing. Two "
            "estimates of that exponent disagree by a factor 3.6 (-0.175 from IFN2->IFN3 dominant-tree increments, "
            "-0.048 from a density term in the dominant h-d fit), so treat the magnitude as uncertain.",
        ))
    elif sp == "pinus_pinaster" and (inputs.management is not None or inputs.thinnings is not None):
        warnings.append(AgbWarning(
            "density_growth_not_modelled", "medium",
            "Per-tree growth does NOT respond to stand density, so thinning never speeds up the trees left "
            "standing and differences between management plans are overstated: a heavier thinning always shows "
            "less standing carbon here. Set density_growth_exponent to switch the effect on (A6).",
        ))
    if any(0 < d < 7.5 for d, _, _, _ in tree):
        warnings.append(AgbWarning(
            "small_tree_extrapolation", "info",
            "At young ages trees are below 7.5 cm dbh, the smallest size in the biomass equations' data; "
            "AGB there is extrapolated (small in absolute terms).",
        ))
    warnings.append(AgbWarning(
        "carbon_scope", "info",
        f"CO2 counts tree biomass above and below ground only (roots scaled like AGB; carbon fraction "
        f"{CARBON_FRACTION[sp]:.1%}, Montero et al. 2005). Soil, litter, dead wood and understory are not included.",
    ))

    co2_factor = CARBON_FRACTION[sp] * CO2_PER_C

    def run(i: int, ratio_q: str):
        n = trees_per_ha_trajectory(ages, dens[i], surv[i], mort, thin, cap, inputs.clearfell_age)
        tr = tree if i == 1 else tree_path(dens[i], surv[i])
        r = [mean_to_dominant_ratio(x, ratio_q, sp, inputs.mean_dominant_table) for x in n]
        agb = [ni * kg * ri / 1000.0 for ni, (_, _, kg, _), ri in zip(n, tr, r)]
        bgb = [ni * root * ri / 1000.0 for ni, (_, _, _, root), ri in zip(n, tr, r)]
        co2 = [(a + b) * co2_factor for a, b in zip(agb, bgb)]
        return n, r, agb, bgb, co2

    n, ratio, agb, bgb, co2 = run(1, "median")
    _, _, agb_lo, _, co2_lo = run(0, "low")
    _, _, agb_hi, _, co2_hi = run(2, "high")
    if inputs.area_ha is not None:
        total_agb = [a * inputs.area_ha for a in agb]
        total_co2 = [c * inputs.area_ha for c in co2]
    else:
        total_agb = total_co2 = None
        warnings.append(AgbWarning("no_area", "info", "Parcel area not provided; results are per hectare only."))

    order = {"high": 0, "medium": 1, "info": 2}
    warnings.sort(key=lambda w: order[w.severity])
    return StandAgbResult(
        species=sp, ages=list(ages), trees_per_ha=n, dominant_tree_agb_kg=[kg for _, _, kg, _ in tree],
        mean_to_dominant_ratio=ratio, agb_mg_ha=agb, agb_mg_ha_low=agb_lo, agb_mg_ha_high=agb_hi,
        bgb_mg_ha=bgb, co2_t_ha=co2, co2_t_ha_low=co2_lo, co2_t_ha_high=co2_hi,
        total_agb_mg=total_agb, total_co2_t=total_co2, growth_modifier=growth_modifier,
        assumed_inputs=assumed, warnings=warnings,
    )


if __name__ == "__main__":
    checks = []
    ages = [5, 10, 20, 30, 40, 60, 80, 100]

    for sp in SPECIES:
        res = project_stand_agb(StandInputs(sp), ages)
        print(f"\n=== {sp} (all defaults) ===")
        print(" age  trees/ha  dom.kg  ratio   AGB Mg/ha (low-high)   BGB Mg/ha   tCO2/ha (low-high)")
        for row in zip(ages, res.trees_per_ha, res.dominant_tree_agb_kg, res.mean_to_dominant_ratio,
                       res.agb_mg_ha, res.agb_mg_ha_low, res.agb_mg_ha_high, res.bgb_mg_ha,
                       res.co2_t_ha, res.co2_t_ha_low, res.co2_t_ha_high):
            print(" {:3d}  {:8.0f}  {:6.0f}  {:5.2f}  {:6.1f} ({:.1f}-{:.1f})   {:9.1f}   {:7.1f} ({:.1f}-{:.1f})".format(*row))
        print(f"  implied tCO2 per Mg AGB at 40 y: {res.co2_t_ha[4] / res.agb_mg_ha[4]:.2f}")
        for w in res.warnings:
            print(f"  [{w.severity}] {w.message}")
        checks.append((f"{sp}: low <= central <= high at every age",
                       all(lo <= c <= hi for lo, c, hi in zip(res.agb_mg_ha_low, res.agb_mg_ha, res.agb_mg_ha_high)), 0.0))
        checks.append((f"{sp}: every assumed input has a warning", len(res.assumed_inputs) <= len(res.warnings), 0.0))
        co2_ok = all(abs(c - (a + b) * CARBON_FRACTION[sp] * CO2_PER_C) < 1e-9
                     for a, b, c in zip(res.agb_mg_ha, res.bgb_mg_ha, res.co2_t_ha))
        checks.append((f"{sp}: CO2 = (AGB + BGB) x C fraction x 44/12", co2_ok, 0.0))

        # Climate multiplier: 1.0 is neutral, >1 grows more, and nothing changes before the scenario starts.
        wet = project_stand_agb(StandInputs(sp), ages, growth_modifier=1.2)
        checks.append((f"{sp}: modifier 1.2 raises AGB at every age", all(w > b for w, b in zip(wet.agb_mg_ha, res.agb_mg_ha)), wet.agb_mg_ha[4]))
        late = project_stand_agb(StandInputs(sp), ages, growth_modifier=0.8, climate_from_age=30)
        same_before = all(abs(x - y) < 1e-9 for t, x, y in zip(ages, late.agb_mg_ha, res.agb_mg_ha) if t <= 30)
        lower_after = all(x < y for t, x, y in zip(ages, late.agb_mg_ha, res.agb_mg_ha) if t > 30)
        checks.append((f"{sp}: modifier from age 30 leaves ages <= 30 unchanged, lowers later ones", same_before and lower_after, 0.0))

    # Supplying everything removes the high/medium warnings; the range then reflects only the
    # mean/dominant ratio spread (IFN3 Q25-Q75), which no Sylvestris input replaces.
    full = StandInputs("quercus_ilex", planting_density_per_ha=500, survival_checks={10: 0.9},
                       annual_mortality_after_last_check=0.005, thinnings={}, area_ha=12.0)
    res = project_stand_agb(full, ages)
    # Supplying every stand input must clear every warning about a MISSING input. What can remain is a
    # caveat about the model's own calibration, which no Sylvestris input replaces - here the oak
    # mean/dominant ratio table (B4c). Those are listed so a new one cannot slip in unnoticed.
    MODEL_CAVEATS = {"oak_ratio_table"}
    left = [w for w in res.warnings if w.severity != "info" and w.code not in MODEL_CAVEATS]
    checks.append(("all inputs given: no high/medium warnings about missing inputs", not left,
                   float(len(left))))
    ratio_only = all(
        abs(lo - c * mean_to_dominant_ratio(n, "low", "quercus_ilex") / r) < 1e-9
        and abs(hi - c * mean_to_dominant_ratio(n, "high", "quercus_ilex") / r) < 1e-9
        for lo, c, hi, n, r in zip(res.agb_mg_ha_low, res.agb_mg_ha, res.agb_mg_ha_high, res.trees_per_ha, res.mean_to_dominant_ratio)
    )
    checks.append(("all inputs given: range comes from the ratio spread only", ratio_only, 0.0))
    checks.append(("total = Mg/ha x area", abs(res.total_agb_mg[-1] - 12.0 * res.agb_mg_ha[-1]) < 1e-9, res.total_agb_mg[-1]))

    # Age convention (A2): the lag adds one year of growth; 0.0 recovers the pre-2026-09-26 behaviour.
    lagged = project_stand_agb(StandInputs("pinus_pinaster"), [40])
    unlagged = project_stand_agb(StandInputs("pinus_pinaster", planting_lag_years=0.0), [40])
    checks.append(("pine: planting lag raises 40 y AGB by 2-6%", 1.02 < lagged.agb_mg_ha[0] / unlagged.agb_mg_ha[0] < 1.06,
                   lagged.agb_mg_ha[0] / unlagged.agb_mg_ha[0]))
    checks.append(("holm oak: no planting lag (bridge is already since-planting)",
                   abs(project_stand_agb(StandInputs("quercus_ilex"), [40]).agb_mg_ha[0]
                       - project_stand_agb(StandInputs("quercus_ilex", planting_lag_years=0.0), [40]).agb_mg_ha[0]) < 1e-9, 0.0))

    # Ecoregion (A3): an input, a coordinate default, and interior when nothing is known.
    checks.append(("pine: default ecoregion is interior (unchanged behaviour)",
                   resolve_pine_ecoregion(StandInputs("pinus_pinaster"))[0] == "interior", 0.0))
    checks.append(("pine: explicit ecoregion is not flagged as assumed",
                   resolve_pine_ecoregion(StandInputs("pinus_pinaster", pine_ecoregion="coastal")) == ("coastal", "supplied", False), 0.0))
    # Maceira (Pontevedra, Sylvestris's own site) is in Vega et al.'s coastal zone; Ourense city is not.
    maceira = resolve_pine_ecoregion(StandInputs("pinus_pinaster", lon_lat=(-8.60, 42.47), elevation_m=250.0))
    checks.append(("pine: Maceira/Pontevedra coordinates give coastal", maceira[0] == "coastal", 0.0))
    ourense = resolve_pine_ecoregion(StandInputs("pinus_pinaster", lon_lat=(-7.86, 42.34), elevation_m=130.0))
    checks.append(("pine: Ourense coordinates give interior", ourense[0] == "interior", 0.0))
    high = resolve_pine_ecoregion(StandInputs("pinus_pinaster", lon_lat=(-8.60, 42.47), elevation_m=900.0))
    checks.append(("pine: above 600 m near the coast is interior", high[0] == "interior", 0.0))
    # The climatic rule (2026-09-27), which applies Vega et al.'s own criterion. Numbers below are the
    # medians actually measured on the tiles: Pontevedra/Maceira 1,490 mm and 10.75 °C; Cáceres 592 mm
    # and 11.47 °C (warm but far too dry); upland Galicia 1,649 mm and 8.18 °C (wet but too cold).
    maceira_clim = resolve_pine_ecoregion(StandInputs("pinus_pinaster", ecoregion_climate=(1490.0, 10.75)))
    checks.append(("pine climate: Maceira's own climate gives coastal", maceira_clim[0] == "coastal", 0.0))
    caceres = resolve_pine_ecoregion(StandInputs("pinus_pinaster", ecoregion_climate=(592.0, 11.47)))
    checks.append(("pine climate: warm but dry (Cáceres) gives interior", caceres[0] == "interior", 0.0))
    checks.append(("pine climate: a one-sided miss is flagged BORDERLINE", "BORDERLINE" in caceres[1], 0.0))
    upland = resolve_pine_ecoregion(StandInputs("pinus_pinaster", ecoregion_climate=(1649.0, 8.18)))
    checks.append(("pine climate: wet but cold (upland Galicia) gives interior", upland[0] == "interior", 0.0))
    checks.append(("pine climate: neither criterion met is not called borderline",
                   "BORDERLINE" not in resolve_pine_ecoregion(StandInputs("pinus_pinaster", ecoregion_climate=(500.0, 7.0)))[1], 0.0))
    # Climate must outrank the geometric proxy, and an explicit ecoregion must outrank both.
    both = StandInputs("pinus_pinaster", ecoregion_climate=(500.0, 7.0), lon_lat=(-8.60, 42.47), elevation_m=250.0)
    checks.append(("pine: climate outranks the coordinate proxy", resolve_pine_ecoregion(both)[0] == "interior", 0.0))
    checks.append(("pine: an explicit ecoregion outranks the climate",
                   resolve_pine_ecoregion(replace(both, pine_ecoregion="coastal")) == ("coastal", "supplied", False), 0.0))

    coastal_agb = project_stand_agb(StandInputs("pinus_pinaster", pine_ecoregion="coastal"), [40]).agb_mg_ha[0]
    checks.append(("pine: coastal curve differs from interior at 40 y", abs(coastal_agb - lagged.agb_mg_ha[0]) > 1.0, coastal_agb))
    try:
        resolve_pine_ecoregion(StandInputs("pinus_pinaster", pine_ecoregion="atlantic"))
        bad_eco = False
    except ValueError:
        bad_eco = True
    checks.append(("pine: an unknown ecoregion raises", bad_eco, 0.0))

    # Site presets (A4): the median preset IS the default, and they order correctly.
    preset_agb = {k: project_stand_agb(StandInputs("pinus_pinaster", pine_site_preset=k), [40]).agb_mg_ha[0]
                  for k in PINE_SITE_PRESETS}
    checks.append(("pine: 'median' preset reproduces the default exactly",
                   abs(preset_agb["median"] - lagged.agb_mg_ha[0]) < 1e-9, preset_agb["median"]))
    checks.append(("pine: poor < median < good at 40 y",
                   preset_agb["poor"] < preset_agb["median"] < preset_agb["good"], preset_agb["good"]))
    # A measured height at the reference age reads back as itself, whatever the ecoregion.
    for eco in PINE_ECOREGIONS:
        c = pine_site_from_measurement(15.0, PINASTER_REFERENCE_AGE - PLANTING_LAG_YEARS, eco)
        checks.append((f"pine site: 15 m at the reference age reads back as SI 15 ({eco})",
                       abs(c.site_index_20_m - 15.0) < 1e-6 and c.recommended == "anchor", c.site_index_20_m))
    young = pine_site_from_measurement(6.0, 6.0)
    checks.append(("pine site: a 6 y measurement selects a preset, not the curve", young.recommended == "preset", young.site_index_20_m))
    checks.append(("pine site: the young choice builds StandInputs with a preset and no anchor",
                   young.to_inputs().pine_site_preset == young.preset and young.to_inputs().pine_site_anchor is None, 0.0))
    old = pine_site_from_measurement(18.0, 25.0)
    checks.append(("pine site: a 25 y measurement anchors directly", old.recommended == "anchor", old.site_index_20_m))
    checks.append(("pine site: the old choice builds StandInputs anchored at the reference age",
                   old.to_inputs().pine_site_anchor == (old.site_index_20_m, PINASTER_REFERENCE_AGE), 0.0))
    checks.append(("pine site: a poor young stand picks 'poor'", pine_site_from_measurement(2.2, 6.0).preset == "poor", 0.0))
    warned = project_stand_agb(StandInputs("pinus_pinaster", pine_site_anchor=(6.0, 6.0)), [40])
    checks.append(("pine site: a young direct anchor is warned about",
                   any(w.code == "young_site_anchor" for w in warned.warnings), 0.0))
    try:
        project_stand_agb(StandInputs("pinus_pinaster", pine_site_preset="excellent"), [40])
        bad_preset = False
    except ValueError:
        bad_preset = True
    checks.append(("pine site: an unknown preset raises", bad_preset, 0.0))

    # Pine default thinning: 400 trees/ha after the second thinning (survival is 1 after it, mortality 0).
    res = project_stand_agb(StandInputs("pinus_pinaster"), [16, 17, 23, 30])
    checks.append(("pine: count drops at the 16.5 y thinning", res.trees_per_ha[1] < res.trees_per_ha[0], res.trees_per_ha[1]))
    checks.append(("pine: 400 trees/ha after the 22.5 y thinning", abs(res.trees_per_ha[2] - 400) < 1e-6, res.trees_per_ha[2]))

    # Management presets (A5): the default is pp2 and identical to the pre-A5 behaviour; the others order
    # as expected, and heavier thinning always lowers AGB because growth ignores density (see A6).
    ages_m = [15, 20, 25, 30, 40, 60]
    pp2 = project_stand_agb(StandInputs("pinus_pinaster", management="pp2"), ages_m)
    obs = project_stand_agb(StandInputs("pinus_pinaster", management="galicia_observed"), ages_m)
    non = project_stand_agb(StandInputs("pinus_pinaster", management="none"), ages_m)
    dflt_m = project_stand_agb(StandInputs("pinus_pinaster"), ages_m)
    checks.append(("pine: the default plan is pp2, bit for bit",
                   all(abs(a - b) < 1e-9 for a, b in zip(dflt_m.agb_mg_ha, pp2.agb_mg_ha)), 0.0))
    checks.append(("pine: galicia_observed holds ~707 trees/ha at 20-25 y",
                   abs(obs.trees_per_ha[1] - 707) < 1e-6 and abs(obs.trees_per_ha[2] - 707) < 1e-6, obs.trees_per_ha[2]))
    checks.append(("pine: 'none' keeps the post-survival density at every age",
                   len(set(round(x, 6) for x in non.trees_per_ha)) == 1, non.trees_per_ha[-1]))
    checks.append(("pine: none >= galicia_observed >= pp2 at 60 y",
                   non.agb_mg_ha[-1] >= obs.agb_mg_ha[-1] >= pp2.agb_mg_ha[-1], non.agb_mg_ha[-1]))
    checks.append(("pine: galicia_observed beats pp2 at 25 y (the known shortfall)",
                   obs.agb_mg_ha[2] > pp2.agb_mg_ha[2] * 1.2, obs.agb_mg_ha[2] / pp2.agb_mg_ha[2]))
    for plan in MANAGEMENT_PRESETS["pinus_pinaster"]:
        res_m = project_stand_agb(StandInputs("pinus_pinaster", management=plan), [30])
        checks.append((f"pine: plan '{plan}' names its source in a warning",
                       any(w.code.startswith("management") for w in res_m.warnings), 0.0))
    try:
        project_stand_agb(StandInputs("pinus_pinaster", management="heavy"), [30])
        bad_plan = False
    except ValueError:
        bad_plan = True
    checks.append(("pine: an unknown management plan raises", bad_plan, 0.0))
    try:
        project_stand_agb(StandInputs("pinus_pinaster", management="pp2", thinnings={20: 500}), [30])
        both = False
    except ValueError:
        both = True
    checks.append(("giving both thinnings and management raises", both, 0.0))

    # Clearfell (A5): optional, never default, and it zeroes the stand from that age on.
    felled = project_stand_agb(StandInputs("pinus_pinaster", clearfell_age=38.0), [30, 38, 40, 60])
    checks.append(("pine: clearfell at 38 y leaves 30 y untouched", abs(felled.agb_mg_ha[0] - pp2.agb_mg_ha[3]) < 1e-9, felled.agb_mg_ha[0]))
    checks.append(("pine: clearfell zeroes AGB from its age on", felled.agb_mg_ha[1:] == [0.0, 0.0, 0.0], 0.0))
    checks.append(("pine: no clearfell by default (the warning stays)",
                   any(w.code == "no_final_felling" for w in dflt_m.warnings)
                   and not any(w.code == "clearfell" for w in dflt_m.warnings), 0.0))

    # Density effect on dominant diameter growth (A6). OFF by default; when on it must behave correctly.
    checks.append(("density: g(N_ref) == 1", abs(pine_density_factor(PINE_DENSITY_REFERENCE_TREES_PER_HA, -0.175) - 1.0) < 1e-12, 0.0))
    checks.append(("density: denser grows slower, sparser faster",
                   pine_density_factor(1500, -0.175) < 1.0 < pine_density_factor(200, -0.175), 0.0))
    ages_d = [1, 5, 10, 16, 17, 20, 23, 30, 40, 60]
    off = project_stand_agb(StandInputs("pinus_pinaster"), ages_d)
    # The stepper with the effect off must reproduce the closed form bit for bit, with and without climate.
    for mod, from_age in ((1.0, 0.0), (1.2, 0.0), (0.8, 30.0)):
        stepped = pine_dominant_path(ages_d, StandInputs("pinus_pinaster"), mod, from_age)
        closed = [dominant_tree("pinus_pinaster", t, StandInputs("pinus_pinaster"), mod, from_age) for t in ages_d]
        same = all(abs(a[i] - b[i]) < 1e-9 for a, b in zip(stepped, closed) for i in range(4))
        checks.append((f"density off: the yearly stepper equals the closed form (modifier {mod}, from {from_age:.0f} y)", same, 0.0))
    on = project_stand_agb(StandInputs("pinus_pinaster", density_growth_exponent=-0.175), ages_d)
    checks.append(("density on: it changes the result", any(abs(a - b) > 0.1 for a, b in zip(on.agb_mg_ha, off.agb_mg_ha)), on.agb_mg_ha[-1]))
    checks.append(("density on: a warning names the uncertainty",
                   any(w.code == "density_growth_effect" for w in on.warnings), 0.0))
    checks.append(("density off with a management plan: the overstatement is warned",
                   any(w.code == "density_growth_not_modelled"
                       for w in project_stand_agb(StandInputs("pinus_pinaster", management="none"), ages_d).warnings), 0.0))
    # A thinning from below keeps the dominant trees, so dominant diameter must not jump at 16.5 y -
    # it may only grow faster afterwards.
    n_pp2 = lambda t: trees_per_ha_trajectory([t], 1046.0, {5.0: 0.9}, 0.0, {16.5: 675.0, 22.5: 400.0})[0]  # noqa: E731
    thinned = StandInputs("pinus_pinaster", density_growth_exponent=-0.175)
    around = pine_dominant_path([16.45, 16.55, 22.45, 22.55], thinned, trees_at_age=n_pp2)
    step_1 = around[1][0] - around[0][0]
    step_2 = around[3][0] - around[2][0]
    # A tenth of a year of growth is ~0.1 cm; a jump would be whole centimetres.
    checks.append(("density on: dominant diameter is continuous across both thinnings",
                   0 < step_1 < 0.2 and 0 < step_2 < 0.2, max(step_1, step_2)))
    # ... and growth AFTERWARDS is faster than it would have been without the thinning.
    after_thin = pine_dominant_path([23, 30], thinned, trees_at_age=n_pp2)
    after_none = pine_dominant_path([23, 30], thinned, trees_at_age=lambda t: 941.4)
    checks.append(("density on: the trees left standing grow faster after a thinning",
                   (after_thin[1][0] - after_thin[0][0]) > (after_none[1][0] - after_none[0][0]),
                   (after_thin[1][0] - after_thin[0][0]) - (after_none[1][0] - after_none[0][0])))
    slow = pine_dominant_path([30], StandInputs("pinus_pinaster", density_growth_exponent=-0.175),
                              trees_at_age=lambda t: 1500.0)[0][0]
    fast = pine_dominant_path([30], StandInputs("pinus_pinaster", density_growth_exponent=-0.175),
                              trees_at_age=lambda t: 300.0)[0][0]
    checks.append(("density on: a stand held at 1,500 trees/ha ends thinner than one at 300", slow < fast, fast - slow))
    # Dominant HEIGHT must stay on the site curve whatever the density does (that is the whole point).
    h_dense = pine_dominant_path([30], StandInputs("pinus_pinaster", density_growth_exponent=-0.175),
                                 trees_at_age=lambda t: 1500.0)[0][1]
    h_sparse = pine_dominant_path([30], StandInputs("pinus_pinaster", density_growth_exponent=-0.175),
                                  trees_at_age=lambda t: 300.0)[0][1]
    checks.append(("density on: dominant height is unaffected by density", abs(h_dense - h_sparse) < 1e-9, h_dense))
    # Climate and density compose multiplicatively on the increment.
    both = pine_dominant_path([40], StandInputs("pinus_pinaster", density_growth_exponent=-0.175), 1.2,
                              trees_at_age=lambda t: 800.0)[0][0]
    clim_only = pine_dominant_path([40], StandInputs("pinus_pinaster"), 1.2)[0][0]
    dens_only = pine_dominant_path([40], StandInputs("pinus_pinaster", density_growth_exponent=-0.175), 1.0,
                                   trees_at_age=lambda t: 800.0)[0][0]
    plain = pine_dominant_path([40], StandInputs("pinus_pinaster"))[0][0]
    checks.append(("density x climate compose multiplicatively", abs(both * plain - clim_only * dens_only) < 1e-6,
                   both * plain - clim_only * dens_only))

    # B4 (oaks): cork thickness phased in with size, the 10 cm holm oak kink smoothed, and an
    # even-aged ratio table offered for oak plantations.
    from spain_growth_curves import quercus_ilex_diameter_from_planting, quercus_suber_cork_thickness_cm
    checks.append(("cork: no cork on a seedling", quercus_suber_cork_thickness_cm(0.0) == 0.0, 0.0))
    checks.append(("cork: ~1 cm at the 5 cm stem Montero & Cañellas photograph",
                   abs(quercus_suber_cork_thickness_cm(5.0) - 1.0) < 1e-9, quercus_suber_cork_thickness_cm(5.0)))
    checks.append((f"cork: the full {CORK_THICKNESS_ASSUMED_CM:.0f} cm from the first-debarking size on",
                   abs(quercus_suber_cork_thickness_cm(CORK_THICKNESS_PHASE_IN[2][0]) - CORK_THICKNESS_ASSUMED_CM) < 1e-9
                   and quercus_suber_cork_thickness_cm(40.0) == CORK_THICKNESS_ASSUMED_CM, 0.0))
    checks.append(("cork: thickness never decreases with size",
                   all(quercus_suber_cork_thickness_cm(d) <= quercus_suber_cork_thickness_cm(d + 0.5) + 1e-12
                       for d in np.arange(0, 30, 0.5)), 0.0))
    checks.append(("cork: phase_in=False restores the flat 2 cm",
                   quercus_suber_cork_thickness_cm(1.0, phase_in=False) == CORK_THICKNESS_ASSUMED_CM, 0.0))
    from spain_growth_curves import quercus_suber_diameter_default, quercus_suber_diameter_over_cork
    d_big = CORK_THICKNESS_PHASE_IN[2][0] + 5.0
    checks.append(("cork: above the first-debarking size the phase-in changes nothing",
                   abs(quercus_suber_diameter_over_cork(d_big)
                       - quercus_suber_diameter_over_cork(d_big, CORK_THICKNESS_ASSUMED_CM)) < 1e-12, 0.0))
    d_young = quercus_suber_diameter_default(curve_age("quercus_suber", 5.0, PLANTING_LAG_YEARS))
    thin_young = quercus_suber_diameter_over_cork(d_young)
    fat_young = quercus_suber_diameter_over_cork(d_young, CORK_THICKNESS_ASSUMED_CM)
    checks.append(("cork: at 5 y the phase-in gives a much thinner stem than the flat 2 cm",
                   thin_young < 0.7 * fat_young, thin_young / fat_young))
    # Holm oak: the rate must fall smoothly, with no step anywhere.
    rates = [quercus_ilex_diameter_from_planting(t + 0.25) - quercus_ilex_diameter_from_planting(t - 0.25)
             for t in np.arange(1.0, 80.0, 0.5)]
    biggest_step = max(abs(b - a) for a, b in zip(rates, rates[1:]))
    checks.append(("holm oak: the diameter growth rate has no step (smoothed kink)", biggest_step < 0.002, biggest_step))
    step_rates = [quercus_ilex_diameter_from_planting(t + 0.25, smooth=False)
                  - quercus_ilex_diameter_from_planting(t - 0.25, smooth=False) for t in np.arange(1.0, 80.0, 0.5)]
    checks.append(("holm oak: smooth=False still shows the old step",
                   max(abs(b - a) for a, b in zip(step_rates, step_rates[1:])) > 0.02, 0.0))
    checks.append(("holm oak: the García-Fayos anchor survives the smoothing",
                   abs(quercus_ilex_diameter_from_planting(25.0) - quercus_ilex_diameter_from_planting(25.0, smooth=False)) < 1e-9, 0.0))
    checks.append(("holm oak: the smoothed curve reaches the adult rate by 60 y",
                   abs((quercus_ilex_diameter_from_planting(60.5) - quercus_ilex_diameter_from_planting(59.5)) - 0.197) < 1e-3, 0.0))
    # B4c ratio table
    ev = project_stand_agb(StandInputs("quercus_ilex", mean_dominant_table="even_aged"), [40])
    au = project_stand_agb(StandInputs("quercus_ilex"), [40])
    checks.append(("oaks: the even-aged ratio table raises plantation AGB", ev.agb_mg_ha[0] > au.agb_mg_ha[0],
                   ev.agb_mg_ha[0] / au.agb_mg_ha[0]))
    checks.append(("oaks: the IFN3 default is warned about",
                   any(w.code == "oak_ratio_table" for w in au.warnings)
                   and not any(w.code == "oak_ratio_table" for w in ev.warnings), 0.0))
    try:
        mean_to_dominant_ratio(500, "median", "quercus_ilex", "made_up")
        bad_table = False
    except ValueError:
        bad_table = True
    checks.append(("oaks: an unknown ratio table raises", bad_table, 0.0))

    # B6/B7 (oaks): the cork oak plantation site class and the oak management plans.
    suber_ages = [20, 40, 100]
    new_cls = project_stand_agb(StandInputs("quercus_suber"), suber_ages)
    old_cls = project_stand_agb(StandInputs("quercus_suber", cork_oak_site_class=CORK_OAK_NATURAL_SITE_CLASS), suber_ages)
    checks.append((f"cork oak: the default site class is now {CORK_OAK_DEFAULT_SITE_CLASS} (plantations)",
                   CORK_OAK_DEFAULT_SITE_CLASS == "I", 0.0))
    checks.append(("cork oak: class I gives more than the old class III at every age",
                   all(a > b for a, b in zip(new_cls.agb_mg_ha, old_cls.agb_mg_ha)),
                   new_cls.agb_mg_ha[1] / old_cls.agb_mg_ha[1]))
    checks.append(("cork oak: the old class III is still reachable", old_cls.agb_mg_ha[1] > 0, old_cls.agb_mg_ha[1]))
    for sp_oak in ("quercus_ilex", "quercus_suber"):
        for plan in MANAGEMENT_PRESETS[sp_oak]:
            res_o = project_stand_agb(StandInputs(sp_oak, management=plan), [40])
            checks.append((f"{sp_oak}: plan '{plan}' runs and names its source",
                           res_o.agb_mg_ha[0] > 0 and any(w.code.startswith("management") for w in res_o.warnings), 0.0))
    checks.append(("oaks: the default plan is 'literature' and equals the pre-B7 behaviour",
                   all(abs(a - b) < 1e-9 for a, b in zip(
                       project_stand_agb(StandInputs("quercus_ilex"), suber_ages).agb_mg_ha,
                       project_stand_agb(StandInputs("quercus_ilex", management="literature"), suber_ages).agb_mg_ha)), 0.0))
    # The cork oak density ceiling belongs to "literature" and is dropped by "none". It only bites at
    # densities well above the default planting, so the test forces one.
    dense = dict(planting_density_per_ha=1200.0)
    capped = project_stand_agb(StandInputs("quercus_suber", management="literature", **dense), [100])
    uncapped = project_stand_agb(StandInputs("quercus_suber", management="none", **dense), [100])
    checks.append(("cork oak: the density ceiling binds at 1,200 planted and 'none' drops it",
                   capped.trees_per_ha[0] < uncapped.trees_per_ha[0], uncapped.trees_per_ha[0] - capped.trees_per_ha[0]))

    # Measurement anchors (2026-09-28): what a measured (diameter, age) pair can do for each oak.
    # Cork oak: it anchors the 2005 curve, so the curve passes through the measurement exactly.
    anch = StandInputs("quercus_suber", cork_oak_diameter_anchor=(20.0, 20.0))
    checks.append(("cork oak: an anchored curve passes through the measurement",
                   abs(dominant_tree("quercus_suber", 20.0, anch)[0] - 20.0) < 0.05,
                   dominant_tree("quercus_suber", 20.0, anch)[0]))
    checks.append(("cork oak: anchoring on a thicker tree raises later diameter",
                   dominant_tree("quercus_suber", 40.0, anch)[0]
                   > dominant_tree("quercus_suber", 40.0, StandInputs("quercus_suber"))[0], 0.0))
    # Holm oak has no site parameter, so a measurement can only shift the clock - but it must still
    # reproduce the measurement exactly at the age it was taken.
    hanch = StandInputs("quercus_ilex", holm_oak_diameter_anchor=(12.0, 20.0))
    checks.append(("holm oak: a shifted curve passes through the measurement",
                   abs(dominant_tree("quercus_ilex", 20.0, hanch)[0] - 12.0) < 0.05,
                   dominant_tree("quercus_ilex", 20.0, hanch)[0]))
    checks.append(("holm oak: equivalent age inverts the curve",
                   abs(quercus_ilex_equivalent_age(10.0)
                       - next(t for t in np.arange(0, 200, 0.01)
                              if quercus_ilex_diameter_from_planting(t) >= 10.0)) < 0.05, 0.0))
    checks.append(("oak anchors: leaving them unset changes nothing",
                   dominant_tree("quercus_suber", 40.0, StandInputs("quercus_suber"))[0]
                   == dominant_tree("quercus_suber", 40.0, StandInputs("quercus_suber", cork_oak_diameter_anchor=None))[0], 0.0))
    checks.append(("cork oak: under-cork inversion round-trips",
                   abs(_under_cork(quercus_suber_diameter_over_cork(11.0)) - 11.0) < 0.01, 0.0))

    # Two thinnings compose correctly with ongoing mortality.
    n = trees_per_ha_trajectory([0, 9, 11, 19, 21], 1000, {5: 0.9}, 0.02, {10: 600, 20: 300})
    checks.append(("custom thinning: 600 just after first thinning (minus 1 y mortality)", abs(n[2] - 600 * 0.98) < 1e-6, n[2]))
    checks.append(("custom thinning: 300 just after second thinning (minus 1 y mortality)", abs(n[4] - 300 * 0.98) < 1e-6, n[4]))

    # Unit sanity: mature stands land in the tens-to-low-hundreds Mg/ha, as ESA CCI reports for Spain.
    for sp in SPECIES:
        agb80 = project_stand_agb(StandInputs(sp), [80]).agb_mg_ha[0]
        checks.append((f"{sp}: 80 y default AGB within 5-300 Mg/ha", 5 <= agb80 <= 300, agb80))

    print()
    for name, ok, value in checks:
        print(f"{'OK  ' if ok else 'FAIL'} {name}: {value:.1f}")
    if not all(ok for _, ok, _ in checks):
        raise SystemExit(1)
