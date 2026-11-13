"""
Runnable implementations of the T3 growth-curve equations for the 3 Sylvestris
species, transcribed directly from the papers now in
use_case_1_simulator/docs/ (see sylvestris_species_growth_plan.md, "T3
results" for the full citations and caveats). Each function is checked
against a number quoted directly in its source paper (see `if __name__ ==
"__main__"` at the bottom) — this catches transcription typos, it does NOT
mean these curves are validated for our actual use case (see the plan doc's
open questions: young-stand coverage for Q. ilex, real climate-modifier
magnitude, anchoring approach).

All three are "dynamic" (base-age-invariant) equations: given a KNOWN
(age, size) pair, they project to a size at another age. Per T6's simplified
design (age is now a user-supplied input, not estimated), two ways to use
these:
  (a) the user also supplies a current observed height/diameter -> use
      `project_from_observation` to solve the implied site index and predict
      forward/backward to any other age.
  (b) the user supplies only species + age, no observation -> assume a
      default (or user-chosen) site-quality class and evaluate directly from
      age 0 via `from_site_class`.
"""

from __future__ import annotations

import math
from dataclasses import dataclass


# --- Pinus pinaster (Álvarez-González et al. 2005, Ann. For. Sci. 62(2):115-127) ---
# Eq. 11, Goelz-Burk-derived Bertalanffy-Richards form. Reference age 20 years.
# I = 0 interior ecoregion, I = 1 coastal ecoregion (Galicia). Which applies to
# Sylvestris's Maceira/Pontevedra site is NOT confirmed — see plan doc.
#
# WHICH AGE (checked 2026-09-26, F5014.pdf; see PINASTER_AGE_DEFINITION below): the paper's t is
# TOTAL age counted from germination. Its data are 204 destructively sampled dominant trees,
# sectioned at the stump, at breast height, at 2.0 m and then at 1 m intervals, with "the age at
# each section height determined in the laboratory" (§2.1) and Figure 5 plotting site index against
# "total age". So t is neither breast-height age nor years since planting.
#
# RANGE OF THE FITTING DATA (F5014.pdf Table I, and §3.2 on the reference age):
#   * STAND ages 8-50 y (coast 8-39, interior 9-50; mean 19.3 y), 102 even-aged plots.
#   * Above 30 y "the sample is not representative (less than 30 trees)" — the paper says this while
#     choosing the reference age, so 30 y is the top of the well-sampled range, not 50.
#   * Within each sampled tree, stem analysis yields height-age pairs at 2-year increments from the
#     stump upwards, so the 22,554 fitted (H1, t1, H2, t2) pairs DO reach ages below 8 y; what is
#     missing below 8 y is any whole stand of that age. Using the curve below 8 y is therefore an
#     extrapolation in stand terms, not in tree terms.
# The stand block turns these into warnings (spain_stand_agb.py: pine_curve_young_extrapolation,
# pine_curve_old_extrapolation).
PINASTER_FIT_STAND_AGE_MIN = 8.0
PINASTER_FIT_STAND_AGE_MAX = 50.0
PINASTER_FIT_WELL_SAMPLED_AGE_MAX = 30.0
PINASTER_REFERENCE_AGE = 20.0
PINASTER_AGE_DEFINITION = "total age from germination (stem analysis, stump section; F5014.pdf §2.1, Fig. 5)"

# --- WHICH ECOREGION? (A3, 2026-09-26) -----------------------------------------------------------
# Neither F5014.pdf nor hal-00883904.pdf defines the two ecoregions; both defer to Vega et al. (1993),
# "Mejora del Pinus pinaster Ait. en Galicia", I Congreso Forestal Español vol. 2, 129-134, now in
# docs/vega_1993_mejora_pinus_pinaster_galicia.pdf. What that paper actually says:
#   * Galician climate studies agree on three zones — coast, an intermediate Mediterranean zone and
#     mountain — and P. pinaster occupies the first two. Site-quality work (Bará & Toval 1983,
#     confirmed by CIF Lourizán) found growth in the coastal and intermediate zones follows different
#     patterns, and the two P. pinaster improvement zones were drawn on that split. So the "interior"
#     ecoregion is the intermediate Mediterranean zone, not the mountains.
#   * The COASTAL zone is seed-classification region 0 of Toval & Vega (1982), characterised by:
#     maximum temperatures 25-36 °C, MINIMUM temperatures 9.4-10.8 °C, mean annual precipitation
#     ~1,500 mm (up to 3,000 mm in the Eume basin and southern Pontevedra), 150-175 rain days a year
#     and mean annual relative humidity 72-79%.
#   * Geographically the coastal improvement zone is "la franja costera de las Rías Bajas, comprendida
#     entre los 0 y 600 mts. de altitud (Provincias de Pontevedra y La Coruña)".
# So the published criterion is CLIMATIC, with a geographic envelope: the Atlantic strip of Pontevedra
# and A Coruña below 600 m. It is not a clean rule on coordinates alone, and we have no digitised
# version of the Toval & Vega map. `pine_ecoregion_from_location` below implements the envelope as an
# explicit, labelled approximation; everything about the numbers in it is ours, not Vega's.
PINASTER_COASTAL_MAX_ELEVATION_M = 600.0
PINASTER_COASTAL_PROVINCES = (15, 36)  # A Coruña, Pontevedra (INE codes)
PINASTER_INTERIOR_PROVINCES = (27, 32)  # Lugo, Ourense
PINASTER_ECOREGION_SOURCE = (
    "Vega et al. 1993 (I Congreso Forestal Español): coastal = Toval & Vega 1982 seed region 0, the "
    "Rías Baixas coastal strip of Pontevedra and A Coruña below 600 m (Tmin 9.4-10.8 °C, P ~1,500 mm)"
)
# Our approximation of the Galician Atlantic/Cantabrian coastline, as (lon, lat) in WGS84: enough
# vertices to put a parcel on the right side of the 600 m contour question, not a real coastline.
_GALICIA_COAST_POINTS = (
    (-8.87, 43.76), (-8.30, 43.58), (-7.70, 43.75), (-7.05, 43.55),  # Rías Altas / Cantabrian edge
    (-9.18, 43.39), (-9.29, 43.01), (-9.03, 42.75), (-8.87, 42.54),  # Costa da Morte / Muros-Noia
    (-8.82, 42.28), (-8.74, 42.05), (-8.87, 41.89),                  # Rías Baixas to the Miño mouth
)
# Below this distance to the coastline our proxy calls a parcel coastal, provided it is also below
# PINASTER_COASTAL_MAX_ELEVATION_M. OUR NUMBER: the Rías Baixas strip Vega describes is roughly this
# wide before the land rises over 600 m. Not in any paper.
PINASTER_COASTAL_MAX_DISTANCE_KM = 40.0

# The CLIMATIC test, which is what Vega et al. actually specify (2026-09-27). Preferred over the
# geometric proxy above whenever a parcel's own climate is available, because it applies the published
# criterion instead of our invented coastal strip.
#   * Mean daily MINIMUM temperature 9.4-10.8 °C. Use a true mean of daily minima (TerraClimate tmmn, or
#     ERA5-Land daily minima averaged -- NOT ERA5-Land's monthly `temperature_2m_min`, which is the
#     month's single coldest hour and reads ~5 °C where the real figure is 10-11 °C).
#   * Mean annual precipitation about 1,500 mm, up to 3,000 mm in the Eume basin and southern Pontevedra.
# Vega quotes a central value, not a threshold, so the cut-off below is OURS: 1,200 mm, low enough to
# admit the drier edge of the coastal strip and still well above the intermediate Mediterranean zone.
# The upper Tmin bound is NOT applied as an exclusion -- warmer than 10.8 °C is still coastal-like, and
# Vega's figures are 1982 station climatology while these data are modern, so present-day minima run
# warmer than his for the same place. That asymmetry is deliberate and biases nothing towards interior.
PINASTER_COASTAL_TMIN_RANGE_C = (9.4, 10.8)
PINASTER_COASTAL_TYPICAL_PRECIP_MM = 1500.0
PINASTER_COASTAL_MIN_PRECIP_MM = 1200.0  # OUR cut-off, not Vega's


def pine_ecoregion_from_climate(annual_precip_mm: float, mean_daily_tmin_c: float) -> tuple[str, str]:
    """("coastal" | "interior", why) from the parcel's own climate, applying Vega et al.'s criteria.

    Both criteria must hold for "coastal". When they disagree the reason says so, because a parcel that
    passes one and fails the other is genuinely borderline and worth a human look."""
    warm = mean_daily_tmin_c >= PINASTER_COASTAL_TMIN_RANGE_C[0]
    wet = annual_precip_mm >= PINASTER_COASTAL_MIN_PRECIP_MM
    obs = f"mean daily minimum {mean_daily_tmin_c:.1f} °C, {annual_precip_mm:.0f} mm/yr"
    crit = (f"Vega et al. 1993 coastal zone = Toval & Vega 1982 seed region 0: minima "
            f"{PINASTER_COASTAL_TMIN_RANGE_C[0]}-{PINASTER_COASTAL_TMIN_RANGE_C[1]} °C, "
            f"~{PINASTER_COASTAL_TYPICAL_PRECIP_MM:.0f} mm/yr (our rainfall cut-off {PINASTER_COASTAL_MIN_PRECIP_MM:.0f} mm)")
    if warm and wet:
        return "coastal", f"{obs}: meets both coastal criteria ({crit})"
    if warm or wet:
        which = "the temperature criterion but not the rainfall one" if warm else "the rainfall criterion but not the temperature one"
        return "interior", f"{obs}: BORDERLINE — meets {which}, so interior is the cautious reading ({crit})"
    return "interior", f"{obs}: meets neither coastal criterion ({crit})"


def pine_ecoregion_from_location(lon: float, lat: float, elevation_m: float | None = None) -> tuple[str, str]:
    """("coastal" | "interior", why) for a Galician parcel, approximating Vega et al.'s (1993) zones.

    PREFER `pine_ecoregion_from_climate` when the parcel's climate is available: it applies Vega's own
    criterion rather than this geometric stand-in for it. This function is for parcels with coordinates
    but no climate.

    Coastal if within PINASTER_COASTAL_MAX_DISTANCE_KM of our coastline polyline AND (when an elevation
    is given) below PINASTER_COASTAL_MAX_ELEVATION_M. The published criterion is climatic; this is a
    geographic stand-in for it, and the returned reason says so. An elevation of None is treated as
    "unknown, assume it passes", which biases towards coastal near the sea."""
    d_km = min(
        math.hypot((lon - clon) * 111.32 * math.cos(math.radians(lat)), (lat - clat) * 110.57)
        for clon, clat in _GALICIA_COAST_POINTS
    )
    high = elevation_m is not None and elevation_m > PINASTER_COASTAL_MAX_ELEVATION_M
    if d_km <= PINASTER_COASTAL_MAX_DISTANCE_KM and not high:
        return "coastal", (
            f"{d_km:.0f} km from the coast"
            + (f" at {elevation_m:.0f} m" if elevation_m is not None else ", elevation unknown")
            + f"; approximates {PINASTER_ECOREGION_SOURCE}"
        )
    why = f"{d_km:.0f} km from the coast" if d_km > PINASTER_COASTAL_MAX_DISTANCE_KM else f"above {PINASTER_COASTAL_MAX_ELEVATION_M:.0f} m"
    return "interior", f"{why}; approximates {PINASTER_ECOREGION_SOURCE}"


def pine_ecoregion_from_province(province_code: int) -> str | None:
    """("coastal" | "interior" | None) from an INE province code, for the IFN4 validation.

    Vega et al. name Pontevedra and A Coruña as the coastal improvement zone and leave Lugo and
    Ourense to the interior one. Crude — inland A Coruña and Lugo are not coastal — but it is the
    split the paper itself states, so it is what the validation tests."""
    if province_code in PINASTER_COASTAL_PROVINCES:
        return "coastal"
    if province_code in PINASTER_INTERIOR_PROVINCES:
        return "interior"
    return None


def _pinaster_z(h1: float, t1: float, ecoregion_coastal: bool) -> float:
    i = 1.0 if ecoregion_coastal else 0.0
    return (0.1352 + 0.0276 * i) * (h1 / t1) ** (0.9831 + 0.0940 * i) * t1 ** (-0.2108 - 0.0929 * i)


def pinus_pinaster_project(h1: float, t1: float, t2: float, ecoregion_coastal: bool) -> float:
    """Dominant height (m) at age t2, given height h1 (m) observed at age t1 (years)."""
    if t1 <= 0 or t2 <= 0:
        raise ValueError("ages must be > 0 (model undefined at age 0)")
    i = 1.0 if ecoregion_coastal else 0.0
    z = _pinaster_z(h1, t1, ecoregion_coastal)
    exponent = 1.4202 + 0.0801 * i
    # exponent applies to the WHOLE ratio (Goelz-Burk M9 form), not just the denominator —
    # an earlier version of this function had that wrong; caught by the self-check below.
    ratio = (1 - math.exp(-z * t2)) / (1 - math.exp(-z * t1))
    return h1 * ratio**exponent


# --- Quercus suber (Sánchez-González, Stiti, Chaar & Cañellas 2010, Forest Systems 19(3):285-298) ---
# Eq. 8, Hossfeld-derived, base-age invariant. Reference age 40 years. Same
# model validated for both Spain and Tunisia (no country difference found).
CORK_OAK_ASYMPTOTE = 23.49
CORK_OAK_SHAPE = 1.35
CORK_OAK_SITE_CLASSES_AT_40Y = {"I": 9.0, "II": 7.0, "III": 5.0, "IV": 3.5, "V": 2.0}


def quercus_suber_project(h0: float, t0: float, t: float) -> float:
    """Dominant height (m) at age t (years), given height h0 (m) observed at age t0 (years)."""
    if t0 <= 0 or t <= 0:
        raise ValueError("ages must be > 0")
    a, b = CORK_OAK_ASYMPTOTE, CORK_OAK_SHAPE
    return a / (1 - (1 - a / h0) * (t0 / t) ** b)


def quercus_suber_from_site_class(site_class: str, t: float) -> float:
    """Dominant height (m) at age t, from one of the 5 standard site-quality classes (at reference age 40)."""
    h40 = CORK_OAK_SITE_CLASSES_AT_40Y[site_class]
    return quercus_suber_project(h40, 40.0, t)


# --- Quercus suber DIAMETER (Sánchez-González, Tomé & Montero 2005, Ann. For. Sci. 62:633-643) ---
# Closes T9's cork oak gap: the height->diameter inversion does not work in dehesa
# (see spain_height_diameter.py), so diameter comes from this paper's own dominant-tree
# diameter-growth model instead. PDF: use_case_1_simulator/docs/f5074.pdf. Stem analysis of
# 85 dominant trees (Los Alcornocales + Catalonia), one model for both regions.
#   Height, eq. 12 (McDill-Amateis): h2 = 20.7216 / (1 - (1 - 20.7216/h1) * (t1/t2)**1.4486).
#     Used ONLY to compute this paper's site index (dominant height at 80 years; classes
#     I-V = 14/12/10/8/6 m). Keeps SI consistent with the one the diameter model was fit on;
#     heights themselves still come from the 2010 curve above.
#   Diameter, eq. 14 (Richards, n free): d2 = A * (d1/A)**r, r = ln(1-e^-k t2) / ln(1-e^-k t1),
#     A = 83.20 + 5.28*SI - 1.53*h/d, k = 0.0063 (Table V). d is dbh UNDER cork (cm); h/d is
#     height (cm) / dbh (cm), a stand-density proxy (sample mean ~34.5, range 20-61).
CORK_OAK_2005_HEIGHT_ASYMPTOTE = 20.7216
CORK_OAK_2005_HEIGHT_SHAPE = 1.4486
CORK_OAK_SI_BASE_AGE = 80.0
CORK_OAK_DIAM_K = 0.0063
# Paper's Table II sample mean, Los Alcornocales (dominant trees): used as the default
# diameter anchor when no measured (d, age) pair is available.
# An anchor fixes the site: projecting the SAME (d, age) point with a different SI makes the
# curves cross (better sites come out thinner when young: 3.1 vs 7.0 cm at 10 y for SI 14 vs 6).
# So the default anchor is only used with its own SI, from the same trees' mean height
# (8.3 m at 65.2 y -> SI ~9.8 m). Other site classes need a measured (d, age) pair.
CORK_OAK_DEFAULT_ANCHOR = {"d_under_cork_cm": 26.7, "age": 65.2, "hd_ratio": 34.3, "height_m": 8.3}
# ASSUMPTION, not from these papers: cork thickness added on each side to convert dbh under
# cork to dbh over cork (what the Ruiz-Peinado biomass equations and the IFN measure).
# Reproduction cork reaches roughly 3 cm at a 9-year harvest, so ~2 cm is a mid-rotation value;
# virgin cork on never-debarked trees differs. Replace with a measurement if Sylvestris has one.
CORK_THICKNESS_ASSUMED_CM = 2.0
# B4a (2026-09-26): the flat 2 cm was added to every tree including seedlings, which inflated young
# cork oak diameters (6 cm of "tree" on a 2 cm stem) and so their biomass and especially their roots.
# Cork thickness is phased in with tree size instead, through two figures from Montero & Cañellas
# (1999, montero_1999_manual_reforestacion_alcornoque.pdf):
#   * a 5-year-old plantation whose stem bases measure 4-7 cm already carries "corcho bornizo de
#     aproximadamente 1 cm" (p. 64-65 photo caption; the text's rabbit-damage note says the bornizo
#     there was over 1 cm thick);
#   * the first debarking is due when the tree reaches 60 cm circumference OVER cork at 1.3 m, i.e.
#     19.1 cm over cork ~ 15.1 cm under cork, reached at 20-25 y on good sites and 35-40 y on poor
#     ones. CORK_THICKNESS_ASSUMED_CM is a mid-rotation value, so it belongs at that size.
# Between and below those points the thickness is interpolated linearly from 0 at planting. The
# 0 -> 1 cm and 1 -> 2 cm anchors are sourced; the straight lines between them are ours.
CORK_THICKNESS_PHASE_IN = ((0.0, 0.0), (5.0, 1.0), (15.1, CORK_THICKNESS_ASSUMED_CM))
CORK_FIRST_DEBARK_CIRCUMFERENCE_CM = 60.0


def quercus_suber_cork_thickness_cm(d_under_cork_cm: float, phase_in: bool = True) -> float:
    """Cork thickness (cm) on one side, at a given dbh under cork. See CORK_THICKNESS_PHASE_IN."""
    if not phase_in:
        return CORK_THICKNESS_ASSUMED_CM
    xs = [d for d, _ in CORK_THICKNESS_PHASE_IN]
    ys = [t for _, t in CORK_THICKNESS_PHASE_IN]
    if d_under_cork_cm >= xs[-1]:
        return ys[-1]
    for (x0, y0), (x1, y1) in zip(CORK_THICKNESS_PHASE_IN, CORK_THICKNESS_PHASE_IN[1:]):
        if d_under_cork_cm <= x1:
            return y0 + (y1 - y0) * (d_under_cork_cm - x0) / (x1 - x0)
    return ys[-1]


def quercus_suber_site_index(h: float, t: float) -> float:
    """Site index (dominant height at 80 years, m) from a dominant height h (m) at age t, via 2005 eq. 12."""
    if t <= 0:
        raise ValueError("age must be > 0")
    a, b = CORK_OAK_2005_HEIGHT_ASYMPTOTE, CORK_OAK_2005_HEIGHT_SHAPE
    return a / (1 - (1 - a / h) * (t / CORK_OAK_SI_BASE_AGE) ** b)


def quercus_suber_diameter_asymptote(site_index_m: float, hd_ratio: float) -> float:
    return 83.20 + 5.28 * site_index_m - 1.53 * hd_ratio


def quercus_suber_diameter_project(d1_cm: float, t1: float, t2: float, site_index_m: float, hd_ratio: float) -> float:
    """Dominant dbh under cork (cm) at age t2, given d1_cm under cork at age t1 (2005 eq. 14)."""
    if t1 <= 0 or t2 <= 0:
        raise ValueError("ages must be > 0")
    a = quercus_suber_diameter_asymptote(site_index_m, hd_ratio)
    if not 0 < d1_cm < a:
        raise ValueError(f"d1_cm={d1_cm} must be in (0, asymptote={a:.1f})")
    k = CORK_OAK_DIAM_K
    r = math.log(1 - math.exp(-k * t2)) / math.log(1 - math.exp(-k * t1))
    return a * (d1_cm / a) ** r


def quercus_suber_diameter_default(t: float) -> float:
    """Dominant dbh under cork (cm) at age t for the paper's average sampled site (see CORK_OAK_DEFAULT_ANCHOR)."""
    a = CORK_OAK_DEFAULT_ANCHOR
    si = quercus_suber_site_index(a["height_m"], a["age"])
    return quercus_suber_diameter_project(a["d_under_cork_cm"], a["age"], t, si, a["hd_ratio"])


def quercus_suber_diameter_over_cork(
    d_under_cm: float, cork_thickness_cm: float | None = None, phase_in: bool = True
) -> float:
    """Dbh over cork (cm), what the IFN measures and the biomass equations expect.

    `cork_thickness_cm` forces one thickness on both sides (the pre-2026-09-26 behaviour is
    `cork_thickness_cm=CORK_THICKNESS_ASSUMED_CM`); by default it is phased in with size (B4a)."""
    t = quercus_suber_cork_thickness_cm(d_under_cm, phase_in) if cork_thickness_cm is None else cork_thickness_cm
    return d_under_cm + 2.0 * t


# --- Quercus ilex (Gea Izquierdo, Cañellas, Garriga & Montero 2007, Cuad. SECF 23:175-180) ---
# Diameter (cm), not height. Two variants: E1 (age-dependent, GADA Hossfeld IV,
# NOT VALID below age 25 or diameter 10 cm — fit on mature dehesa trees) and T1
# (age-independent, needs only current diameter + stand density, no age at all).

QILEX_E1_B1 = 273.2923
QILEX_E1_C = 1.248261
QILEX_MIN_VALID_AGE = 25.0
QILEX_MIN_VALID_DIAMETER_CM = 10.0

QILEX_T1_AD1 = 0.505583
QILEX_T1_BD1 = 0.025352
QILEX_T1_C = 1.030452


def quercus_ilex_age_dependent_project(y1_cm: float, t1: float, t2: float) -> float:
    """
    Diameter (cm) at age t2, given diameter y1_cm observed at age t1 — E1, GADA
    Hossfeld IV, not density-expanded. NOT VALID for t1 or t2 < 25 years, or
    y1_cm < 10 cm — the paper explicitly restricts this range (mature dehesa
    trees only); raises ValueError outside it rather than silently extrapolating.
    """
    if t1 < QILEX_MIN_VALID_AGE or t2 < QILEX_MIN_VALID_AGE:
        raise ValueError(f"E1 is not valid below age {QILEX_MIN_VALID_AGE} (paper's stated range)")
    if y1_cm < QILEX_MIN_VALID_DIAMETER_CM:
        raise ValueError(f"E1 is not valid below diameter {QILEX_MIN_VALID_DIAMETER_CM} cm")
    c, b1 = QILEX_E1_C, QILEX_E1_B1
    x0 = (t1**c / y1_cm) / (b1 + t1**c)
    return t2**c / (x0 * (b1 + t2**c))


# --- Quercus ilex JUVENILE BRIDGE (ages 0-25+, where E1 is not valid) ---------------------
# Our construction, anchored on one published measurement — flag it as such to reviewers.
# García-Fayos et al. (2024), New Forests, doi:10.1007/s11056-024-10031-3: in 153 CAP holm oak
# plantations on former cropland, the tallest tree 25 years after planting averages 4.2 m
# (95% CI 4.04-4.37). Converted to diameter with the Caceres-dehesa h-d fit
# (spain_height_diameter.hd_params_dehesa): 4.2 m -> ~8.3 cm. The bridge is a straight line
# from 0 cm at planting through that point, continued until it reaches E1's 10 cm floor
# (~30 years), where E1 takes over, so the handover sits inside E1's valid range.
# Check against the two other dominant heights García-Fayos cites from other plantations
# (~2.4 m at 10 y, ~3.45 m at 15 y): the bridge + dehesa h-d gives ~2.9 m and ~3.4 m.
# Caveats: planted cropland (deep soils) may grow faster than forest plantings; a straight
# line ignores the slow establishment years; one anchor only.
QILEX_JUVENILE_ANCHOR_AGE = 25.0
QILEX_JUVENILE_ANCHOR_HEIGHT_M = 4.2

# ADULT PHASE, calibrated on inventory remeasurements (2026-09-25; ifn_remeasurement_validation.py).
# IFN2 -> IFN3 strict tree pairs show holm oak diameter growth does NOT slow with size: in Extremadura
# (Badajoz + Caceres, 7,506 trees >= 10 cm, ~10-year interval) the median is 0.197 cm/yr (bootstrap
# CI 0.191-0.199), flat from 10 cm to > 60 cm (class medians 1.8-2.2 cm/decade; the plot's 2 largest
# trees give the same 2.0). Nationally (45,678 trees, 38 provinces) 1.6-2.2 cm/decade, slightly rising
# with size. E1 instead slows to ~0.4 cm/decade by 35-43 cm and cannot even reach the sizes of half the
# dominant Caceres trees. So by default the adult phase grows at this measured rate; "e1" keeps the
# published curve for comparison.
# CAVEAT: flat-with-size in a cross-section does not prove flat-with-age for one tree. Large trees
# are partly old and partly on good sites / fast growers; that mixing can flatten the cross-sectional
# curve even if each tree slows with age. What the IFN data establish firmly is the regional LEVEL
# (~0.2 cm/yr); the trajectory SHAPE is a pragmatic choice. The two options differ by 10-15% in
# diameter at 40-100 y (ifn 12.0/15.9/23.8 cm vs e1 13.2/18.6/26.3 cm at 40/60/100 y).
QILEX_IFN_RATE_CM_PER_YR = 0.197


def _qilex_juvenile_rate_cm_per_yr() -> float:
    from spain_height_diameter import diameter_from_height, hd_params_dehesa

    d25 = diameter_from_height(QILEX_JUVENILE_ANCHOR_HEIGHT_M, hd_params_dehesa("quercus_ilex"))
    return d25 / QILEX_JUVENILE_ANCHOR_AGE


# B4b (2026-09-26): the bridge (0.334 cm/yr) used to meet the adult rate (0.197) as a STEP at 10 cm.
# IFN2 -> IFN3 remeasurements in Extremadura show the decline is gradual, not a step: median annual
# diameter increment by IFN2 size class is 0.249 cm/yr below 10 cm, 0.214 at 10-15, 0.200 at 15-20,
# 0.193 at 20-30 and 0.192 above 30 cm (7,506 trees; b3 in the plan doc). So the rate now falls
# linearly IN DIAMETER from the bridge rate to the adult rate between the juvenile anchor's own
# diameter and QILEX_SMOOTH_TO_DIAMETER_CM, which removes the kink and keeps the curve C1.
# The lower end is the anchor diameter itself, so García-Fayos's measured 4.2 m at 25 y is untouched;
# the upper end is where the measured class medians have flattened. Set smooth=False for the old step.
QILEX_SMOOTH_TO_DIAMETER_CM = 15.0


def _qilex_rate_at(d_cm: float, r0: float, r1: float, d0: float, d1: float) -> float:
    """Annual diameter increment at diameter d: r0 below d0, r1 above d1, linear between."""
    if d_cm <= d0:
        return r0
    if d_cm >= d1:
        return r1
    return r0 + (r1 - r0) * (d_cm - d0) / (d1 - d0)


def quercus_ilex_diameter_from_planting(t: float, adult: str = "ifn", smooth: bool = True) -> float:
    """Dominant diameter (cm) at t years after planting.

    Juvenile bridge, then the adult phase: the IFN-measured Extremadura rate (adult="ifn", default) or
    the published E1 curve (adult="e1"). With smooth=True (default, adult="ifn" only) the juvenile rate
    eases into the adult rate between the anchor diameter and QILEX_SMOOTH_TO_DIAMETER_CM instead of
    stepping down at 10 cm; smooth=False restores the pre-2026-09-26 step."""
    if t < 0:
        raise ValueError("age must be >= 0")
    r0 = _qilex_juvenile_rate_cm_per_yr()
    if adult == "ifn" and smooth:
        r1 = QILEX_IFN_RATE_CM_PER_YR
        d0 = r0 * QILEX_JUVENILE_ANCHOR_AGE  # the anchor's own diameter, so the anchor is preserved
        d1 = QILEX_SMOOTH_TO_DIAMETER_CM
        t0 = QILEX_JUVENILE_ANCHOR_AGE  # time to reach d0, by construction
        if t <= t0:
            return r0 * t
        # Between d0 and d1 the rate is r0 + k (d - d0), so t - t0 = ln(1 + k (d - d0) / r0) / k.
        k = (r1 - r0) / (d1 - d0)
        t1 = t0 + math.log(1 + k * (d1 - d0) / r0) / k
        if t <= t1:
            return d0 + r0 * (math.exp(k * (t - t0)) - 1.0) / k
        return d1 + r1 * (t - t1)
    t_handover = QILEX_MIN_VALID_DIAMETER_CM / r0
    if t <= t_handover:
        return r0 * t
    if adult == "ifn":
        return QILEX_MIN_VALID_DIAMETER_CM + QILEX_IFN_RATE_CM_PER_YR * (t - t_handover)
    if adult == "e1":
        return quercus_ilex_age_dependent_project(QILEX_MIN_VALID_DIAMETER_CM, t_handover, t)
    raise ValueError(f"adult must be 'ifn' or 'e1', got {adult!r}")


def quercus_ilex_age_independent_project(y1_cm: float, density_per_ha: float, years_elapsed: float) -> float:
    """
    Diameter (cm) after `years_elapsed` more years, given current diameter
    y1_cm and stand density (trees/ha) — T1, no age needed at all. Less
    accurate than the age-dependent form (paper's cross-validation RMSE 4.59
    vs 1.72 cm) but usable when age genuinely isn't known.
    """
    a = QILEX_T1_AD1 / density_per_ha
    b = QILEX_T1_BD1 * density_per_ha
    c = QILEX_T1_C
    t1 = (y1_cm * b / (1 - a * y1_cm)) ** (1.0 / c)
    t2 = t1 + years_elapsed
    return t2**c / (b + a * t2**c)


if __name__ == "__main__":
    # Checks against numbers quoted directly in each source paper.
    checks = []

    # Pine: paper's own worked example is the ecoregion parameter tables, not a single
    # numeric worked height — sanity-check instead that the curve is monotonic and
    # roughly matches the paper's Table I mean dominant height/age for the interior
    # ecoregion (age~19.3, Ho~10.64m -> project forward to age 20 should stay close).
    h20 = pinus_pinaster_project(h1=10.64, t1=19.3, t2=20.0, ecoregion_coastal=False)
    checks.append(("pine ~20yr from ~19.3yr obs stays close to 10.64m", abs(h20 - 10.64) < 1.0, h20))
    identity = pinus_pinaster_project(h1=14.0, t1=20.0, t2=20.0, ecoregion_coastal=True)
    checks.append(("pine base-age-invariant identity: same age in/out returns same height", abs(identity - 14.0) < 1e-6, identity))
    identity_suber = quercus_suber_project(h0=6.3, t0=27.0, t=27.0)
    checks.append(("cork oak identity: same age in/out returns same height", abs(identity_suber - 6.3) < 1e-6, identity_suber))

    # Cork oak: paper's site index definition IS the number to check against —
    # site class I should read back exactly 9m at the reference age of 40.
    h40 = quercus_suber_from_site_class("I", 40.0)
    checks.append(("cork oak class I at age 40 == 9.0m", abs(h40 - 9.0) < 1e-6, h40))

    # Q. ilex age-dependent: paper's Table 2 mean sample is DAP~30.8cm at age~77 —
    # project that same tree to a nearby age and check it stays in a plausible range.
    d80 = quercus_ilex_age_dependent_project(y1_cm=30.81, t1=77.0, t2=80.0)
    checks.append(("Q. ilex E1: DAP grows (not shrinks) from age 77->80", d80 > 30.81, d80))

    # Cork oak diameter (2005): site index definition reads back exactly at 80 years.
    si = quercus_suber_site_index(10.0, 80.0)
    checks.append(("cork oak 2005 SI: h=10 m at 80 y -> SI 10", abs(si - 10.0) < 1e-9, si))
    a = CORK_OAK_DEFAULT_ANCHOR
    ident = quercus_suber_diameter_project(a["d_under_cork_cm"], a["age"], a["age"], 10.0, a["hd_ratio"])
    checks.append(("cork oak diameter identity at anchor age", abs(ident - a["d_under_cork_cm"]) < 1e-9, ident))
    # Reciprocity (paper: pairs can be swapped): 65->100->65 returns the start.
    fwd = quercus_suber_diameter_project(26.7, 65.2, 100.0, 10.0, 34.3)
    back = quercus_suber_diameter_project(fwd, 100.0, 65.2, 10.0, 34.3)
    checks.append(("cork oak diameter reciprocity 65->100->65", abs(back - 26.7) < 1e-6, back))
    # Paper's sample range at 30-158 y is 10.8-52.5 cm: from the anchor, age 30 and 150 stay inside it.
    d30 = quercus_suber_diameter_project(26.7, 65.2, 30.0, 10.0, 34.3)
    d150 = quercus_suber_diameter_project(26.7, 65.2, 150.0, 10.0, 34.3)
    checks.append(("cork oak diameter at 30 y within sampled 10.8-52.5 cm", 10.8 <= d30 <= 52.5, d30))
    checks.append(("cork oak diameter at 150 y within sampled 10.8-52.5 cm", 10.8 <= d150 <= 52.5, d150))

    # Q. ilex juvenile bridge: passes through the anchor, and is continuous at the E1 handover.
    from spain_height_diameter import height_from_diameter, hd_params_dehesa

    p = hd_params_dehesa("quercus_ilex")
    h25 = height_from_diameter(quercus_ilex_diameter_from_planting(25.0), p)
    checks.append(("Q. ilex bridge: height at 25 y == García-Fayos 4.2 m", abs(h25 - 4.2) < 1e-6, h25))
    h15 = height_from_diameter(quercus_ilex_diameter_from_planting(15.0), p)
    checks.append(("Q. ilex bridge: height at 15 y within 0.5 m of cited ~3.45 m", abs(h15 - 3.45) < 0.5, h15))
    t_h = QILEX_MIN_VALID_DIAMETER_CM / _qilex_juvenile_rate_cm_per_yr()
    for adult in ("ifn", "e1"):
        jump = quercus_ilex_diameter_from_planting(t_h + 1e-6, adult) - quercus_ilex_diameter_from_planting(t_h - 1e-6, adult)
        checks.append((f"Q. ilex bridge continuous at handover age {t_h:.1f} y ({adult})", abs(jump) < 1e-3, jump))
    rate_100 = quercus_ilex_diameter_from_planting(101.0) - quercus_ilex_diameter_from_planting(100.0)
    checks.append(("Q. ilex adult (ifn): 0.197 cm/yr at age 100, no slowdown", abs(rate_100 - QILEX_IFN_RATE_CM_PER_YR) < 1e-9, rate_100))

    for name, ok, value in checks:
        print(f"{'OK  ' if ok else 'FAIL'} {name}: {value:.3f}")
    if not all(ok for _, ok, _ in checks):
        raise SystemExit(1)
