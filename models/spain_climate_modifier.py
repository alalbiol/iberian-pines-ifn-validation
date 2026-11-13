"""
Climate-sensitivity modifier for the Sylvestris species growth curves (T4 in
sylvestris_species_growth_plan.md), to be layered on top of the T3 age/diameter
growth curves (none of which take a climate input themselves — confirmed).

Design, informed by the T4 literature pass (direction solid, magnitude NOT yet
confirmed — see the plan doc's "T4 results" section for caveats and sources):

  - A water-balance drought index (precipitation minus potential evaporation,
    both directly available as ERA5-Land monthly bands — no need to estimate
    PET from temperature). Standardized against each parcel's own long-run
    mean/SD over a rolling window, giving a SPEI-LIKE index (explicitly not
    the accredited SPEI methodology, which fits a log-logistic distribution
    per calendar month — this is a documented simplification, flag before
    treating the number as a real SPEI value).
  - Season-specific precipitation windows: previous-Sep-to-May for the oaks
    (antecedent moisture), current May-Jun for the pine (growing-season
    timing) — per species, not one national window.
  - Summer (Jun-Aug) mean max temperature for the pine specifically.

The modifier is a multiplicative factor on the growth curve's INCREMENT (not
the absolute size), roughly monotonic in the drought index, clipped to a
plausible range so a single bad year cannot blow up the curve.

Status as of 2026-09-23, after reading both previously-blocked papers in full
(PDFs supplied by the user, now in the original project reference library/):

- **Quercus suber**: magnitude now REAL, not a placeholder. Camarero, J.J.,
  Sánchez-Miranda, A., Colangelo, M., Matías, L. (2024). "Climatic drivers of
  cork growth depend on site aridity." Science of the Total Environment 912,
  169574. https://doi.org/10.1016/j.scitotenv.2023.169574 — regresses cork
  ring-width index directly on a standardized SPEI drought index across 4
  sites of varying aridity (NE Spain to SW Morocco): slope 0.35 (95% CI
  0.24-0.50) at the driest/most continental site (14-month August SPEI),
  0.28 (CI 0.18-0.39) at a wetter site (1-4 month SPEI). None of the 4
  studied sites is Sylvestris's actual Galicia/Pontevedra region, so this is
  an analog transfer across the aridity gradient, not a site-specific
  number — still real literature, not a guess. See quercus_suber_modifier
  for how this is used.
- **Quercus ilex**: magnitude now REAL (2026-09-24). Gea-Izquierdo, G.,
  Fernández-de-Uña, L., Cañellas, I. (2013). "Growth projections reveal local
  vulnerability of Mediterranean oaks with rising temperatures." Forest
  Ecology and Management 305, 282-293.
  https://doi.org/10.1016/j.foreco.2013.05.058 — a fitted multiplicative
  model of basal-area increment on hydrological-year (Nov-Oct) precipitation
  and mean minimum temperature, from 91 dominant holm oaks at 4 sites along a
  temperature gradient in central-west Spain (Table 3, model #6; validated on
  independent 1925-1950 data). See quercus_ilex_modifier. Natalini et al.
  (2016, Dendrochronologia 38:51-60) still supports the direction and the
  antecedent-precipitation window, but has no slope of its own.
- **Pinus pinaster**: magnitude now DERIVED (2026-09-24), not read off a
  table — see PINE_SPEI_SLOPE for the derivation and its caveats. Sources:
  Rozas, V., Lamas, S., García-González, I. (2009). "Differential tree-growth
  responses to local and large-scale climatic variation in two Pinus and two
  Quercus species in northwest Spain." Écoscience 16(3), 299-310
  (https://doi.org/10.2980/16-3-3212) — the Galicia match; and Bogino, S.,
  Bravo, F. (2008). "Growth response of Pinus pinaster Ait. to climatic
  variables in central Spanish forests." Annals of Forest Science 65, 506
  (https://doi.org/10.1051/forest:2008025) — the range check. Neither paper
  gives a per-unit slope: both report response functions (standardized
  coefficients per month, shown only as figures), so the slope is derived from
  the variance they explain and the chronology's variability.

TEMPERATURE VARIABLE WARNING (found 2026-09-24): fetch_spain_climate_monthly.py
reads ERA5-Land MONTHLY_AGGR temperature_2m_min/_max, which are the month's
EXTREME hourly values, not means of daily minima/maxima. Evidence: Cáceres
tiles give an annual mean of temp_min_c ~5 °C and a July temp_max_c ~38.5 °C,
against ~10-11 °C mean daily minimum and ~34-35 °C mean July daily maximum
for the region. The dendro papers use mean daily minima/maxima (CRU TS tmn/
tmx). So temp_min_c / temp_max_c must NOT be fed to a modifier fitted on
those. The holm oak modifier needs mean daily minimum; the pipeline does not
provide it yet (needs an ERA5-Land DAILY_AGGR fetch averaged per month).
"""

from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass
class MonthlyClimate:
    """One month of climate at a parcel. Matches fetch_spain_climate_monthly.py's output columns.

    temp_max_c / temp_min_c are the month's extreme hourly values (ERA5-Land
    MONTHLY_AGGR), not means of daily maxima/minima — see module docstring.
    """

    year: int
    month: int
    temp_mean_c: float
    temp_max_c: float
    temp_min_c: float
    precip_mm: float
    pet_mm: float

    @property
    def water_balance_mm(self) -> float:
        return self.precip_mm - self.pet_mm


def rolling_water_balance(months: list[MonthlyClimate], end_index: int, window: int) -> float:
    """Sum of water balance over `window` months ending at end_index (inclusive)."""
    start = end_index - window + 1
    if start < 0:
        raise ValueError("not enough history for this window")
    return sum(m.water_balance_mm for m in months[start : end_index + 1])


def spei_like(months: list[MonthlyClimate], end_index: int, window: int) -> float:
    """
    Standardized rolling water balance ("SPEI-like", see module docstring caveat).
    Standardizes against the mean/SD of the SAME window computed at every valid
    end_index in the series (a full climatology would standardize per calendar
    month across many years — this pools all months, a simplification for v1).
    """
    series = [rolling_water_balance(months, i, window) for i in range(window - 1, len(months))]
    mean = sum(series) / len(series)
    sd = (sum((x - mean) ** 2 for x in series) / len(series)) ** 0.5
    if sd == 0:
        return 0.0
    value = rolling_water_balance(months, end_index, window)
    return (value - mean) / sd


def season_precip_mm(months: list[MonthlyClimate], year: int, start_month: tuple[int, int], end_month: tuple[int, int]) -> float:
    """
    Sum precip over a season window given as (year_offset, month), e.g. the
    oak antecedent window is start=(-1, 9) [previous September] end=(0, 5)
    [current May], relative to `year`.
    """
    by_key = {(m.year, m.month): m for m in months}
    (so, sm), (eo, em) = start_month, end_month
    total = 0.0
    y, mo = year + so, sm
    while (y, mo) <= (year + eo, em):
        if (y, mo) in by_key:
            total += by_key[(y, mo)].precip_mm
        mo += 1
        if mo > 12:
            mo = 1
            y += 1
    return total


def summer_max_temp_c(months: list[MonthlyClimate], year: int) -> float:
    vals = [m.temp_max_c for m in months if m.year == year and m.month in (6, 7, 8)]
    if not vals:
        raise ValueError(f"no Jun-Aug data for {year}")
    return sum(vals) / len(vals)


# --- Species-specific modifiers -------------------------------------------

# Derived slope of ring-width index per SPEI-like unit for P. pinaster.
# A residual chronology with near-zero autocorrelation has index SD ~ MS / 1.128
# (MS = mean sensitivity; for iid normal noise E|x_t - x_t-1| = 1.128 SD), and
# the index change per 1 SD of the climate signal it explains is ~sqrt(r2) * SD.
#   Rozas et al. 2009, Galicia, total ring: MS 0.253, AC 0.016, r2 0.664
#     -> SD ~0.224 -> sqrt(0.664) * 0.224 = 0.18.
#   Bogino & Bravo 2008, 9 central-Spain sites: SD 0.165-0.311, R2 0.43-0.65
#     -> 0.11-0.25.
# Rozas's dominant driver for P. pinaster is water balance (P - PET) from winter
# to summer of the growth year, the same quantity spei_like() standardizes, so the
# slope is applied to a 9-month SPEI-like index. Caveats: response-function r2
# comes from ~20 monthly predictors on ~50 years, so it is inflated and this
# slope leans high; it lumps temperature into the water-balance term, so no
# separate temperature term is added (it would double-count, as for cork oak).
PINE_SPEI_SLOPE = 0.18


def pinus_pinaster_modifier(spei9: float) -> float:
    """Derived magnitude (see PINE_SPEI_SLOPE): +/-0.18 growth per SPEI-like unit, capped."""
    return max(0.5, min(1.5, 1.0 + PINE_SPEI_SLOPE * max(min(spei9, 2.5), -2.5)))


# Ring-derived slopes for the four other Iberian pines (added 2026-11-11), replacing the transfer of
# PINE_SPEI_SLOPE to these species. Each is the change in ring-width index (relative annual radial
# growth) per +1 SD of the site's own hydrological-year water balance (TerraClimate P - PET),
# pooled over NOAA ITRDB raw ring-width sites with a random-effects model that treats co-located
# sites as one unit. Full analysis, robustness grid and site audit:
# analysis/ring_climate/results/RING_CLIMATE_FINDINGS.md
#   halepensis 0.20 (95% 0.16-0.24): 26 sites (2 Spanish, rest N. Africa/Mediterranean); sensitivity
#              falls with wetter climate, value is the aridity meta-regression at the median IFN
#              Aleppo-plot aridity (P/PET 0.41).
#   pinea      0.12 (95% 0.08-0.16): the 5 Spanish sites. The all-site aridity regression at the IFN
#              median aridity gives 0.107, so the value does not hinge on excluding the wetter
#              Italian coast sites.
#   nigra      0.04 (95% 0.02-0.06): 30 Spanish sites. Archive sites sit above most IFN stands and
#              lower sites respond more; at the IFN median elevation ~0.05-0.07 (extrapolated).
#   sylvestris 0.02 (95% 0.005-0.04): 32 Spanish sites; marginal (interval includes zero with
#              ERA5-Land rain). Lower sites respond more; ~0.04-0.056 at the IFN median elevation.
# These are year-to-year sensitivities: a shift held for decades, mortality and drought damage
# are outside what they measure. Scale the shift by the species' own water-balance SD
# (other_pine_climate_reference.json), not by the Pontevedra reference used for maritime pine.
OTHER_PINE_WB_SLOPE = {
    "pinus_halepensis": 0.20,
    "pinus_pinea": 0.12,
    "pinus_nigra": 0.04,
    "pinus_sylvestris": 0.02,
}
OTHER_PINE_WB_SLOPE_CI = {
    "pinus_halepensis": (0.16, 0.24),
    "pinus_pinea": (0.08, 0.16),
    "pinus_nigra": (0.02, 0.06),
    "pinus_sylvestris": (0.005, 0.04),
}


def other_pine_modifier(species: str, wb_shift_sd: float) -> float:
    """Ring-derived multiplier on diameter increment for a water-balance shift in local SD units.

    Same form and caps as pinus_pinaster_modifier: linear in the shift, shift capped at +/-2.5 SD,
    multiplier clipped to 0.5-1.5.
    """
    if species not in OTHER_PINE_WB_SLOPE:
        raise ValueError(f"No ring-derived climate slope for {species!r}")
    slope = OTHER_PINE_WB_SLOPE[species]
    return max(0.5, min(1.5, 1.0 + slope * max(min(wb_shift_sd, 2.5), -2.5)))


# Gea-Izquierdo et al. 2013, Table 3, Q. ilex (QUIL), maximum-likelihood estimates.
# The text extraction drops minus signs; p2 is negative by the paper's own
# definition ("negative to make the function monotonically increasing"), and its
# support interval (-1.1012, -0.6684) brackets -0.8373.
ILEX_P1_MM = 333.0670  # half-saturation precipitation (hydrological year)
ILEX_P2 = -0.8373
ILEX_TOPT_C = 8.9677  # optimum mean minimum temperature (t01)
ILEX_TSCALE_C = 4.6216  # Laplace scale (t02)
# Calibration range of the fit (Table 1): Ppt 210-1066 mm, Tmin 5.4-12.0 °C.
ILEX_PPT_RANGE_MM = (210.1, 1066.0)
ILEX_TMIN_RANGE_C = (5.4, 12.0)
HYDRO_YEAR_WINDOW = ((-1, 11), (0, 10))  # previous November -> current October

# Temperature term OFF by default (decision 2026-09-25, from IFN2->IFN3 remeasurement tests; see
# ifn_remeasurement_validation.py and sylvestris_species_growth_plan.md). Across Caceres (952 plots) and
# nationally (thousands of plots, province fixed effects), size- and density-adjusted holm oak growth
# shows the paper's PRECIPITATION response (fitted strength ~0.7-2 vs 1 published) but NO temperature
# optimum: growth is flat from ~4 to ~12 C mean daily minimum, where the Laplace term predicts a peak
# at 9 C and a ~45% drop on either side. The term is kept behind this switch, not deleted.
ILEX_USE_TEMPERATURE_TERM = False


def quercus_ilex_growth_index(hydro_precip_mm: float, hydro_mean_tmin_c: float, use_temperature: bool = True) -> float:
    """Unnormalised g1(Ppt) * g2(Tmin) from Gea-Izquierdo et al. 2013 (MG dropped).

    g1 is a logistic saturation in precipitation; g2 is a Laplace peak in mean
    minimum temperature (optimum ~9 °C, growth falls on either side). With
    use_temperature=False, g2 is dropped (see ILEX_USE_TEMPERATURE_TERM).
    """
    g1 = 1.0 / (1.0 + (max(hydro_precip_mm, 1e-6) / ILEX_P1_MM) ** ILEX_P2)
    if not use_temperature:
        return g1
    g2 = math.exp(-abs(hydro_mean_tmin_c - ILEX_TOPT_C) / ILEX_TSCALE_C)
    return g1 * g2


def quercus_ilex_modifier(
    hydro_precip_mm: float,
    hydro_mean_tmin_c: float,
    baseline_precip_mm: float,
    baseline_mean_tmin_c: float,
    use_temperature: bool = ILEX_USE_TEMPERATURE_TERM,
) -> float:
    """
    Real, cited magnitude (Gea-Izquierdo et al. 2013): growth relative to the
    parcel's own baseline climate, as the ratio of the paper's climate terms.
    The paper models basal-area increment; applied here to the T3 diameter
    increment as a multiplicative factor, which is exact for proportional BAI
    changes at fixed tree size.

    hydro_mean_tmin_c must be the MEAN DAILY MINIMUM over Nov-Oct (CRU tmn
    equivalent), NOT MonthlyClimate.temp_min_c (a monthly extreme, ~5 °C lower
    here). This is direction-critical, not just scale: at a parcel above the
    ~9 °C optimum (most of Extremadura) warming reduces growth, but feeding the
    colder extreme would place it below the optimum and reverse the sign.

    Wider clip than the other species ([0.25, 2.0] vs [0.5, 1.5]): the paper's
    Laplace term is steep (x0.65 per +2 °C above the optimum) and the tighter
    clip would silently flatten a published scenario response.
    """
    ratio = quercus_ilex_growth_index(hydro_precip_mm, hydro_mean_tmin_c, use_temperature) / quercus_ilex_growth_index(
        baseline_precip_mm, baseline_mean_tmin_c, use_temperature
    )
    return max(0.25, min(2.0, ratio))


CORK_SPEI_SLOPE = 0.3  # Camarero et al. 2024: 0.28-0.35 across sites, see module docstring.


def quercus_suber_modifier(spei_value: float) -> float:
    """
    Real, cited magnitude (Camarero et al. 2024, see module docstring): cork
    ring-width index regresses on a standardized SPEI-like drought index with
    slope ~0.3 (site range 0.28-0.35). `spei_value` should come from this
    module's own spei_like() — the paper's mechanism (warm spring/early-summer
    -> higher PET -> lower water balance -> less cork) is a temperature effect
    already embedded IN the SPEI water-balance calculation (precip minus PET),
    so no separate temperature term is added here: doing so on top of SPEI
    would double-count warming's effect.
    """
    return max(0.5, min(1.5, 1.0 + CORK_SPEI_SLOPE * max(min(spei_value, 2.5), -2.5)))


OAK_ANTECEDENT_WINDOW = ((-1, 9), (0, 5))  # previous September -> current May
PINE_GROWING_WINDOW = ((0, 5), (0, 6))  # current May -> current June
