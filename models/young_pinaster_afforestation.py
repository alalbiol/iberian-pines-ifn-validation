"""Published young-afforestation models for Spanish Pinus pinaster.

Calama et al. (2024), Forests 15:827, doi:10.3390/f15050827.
Height: model M1 in Tables 3 and 5, calibrated on 51 stem analyses aged
4-39 years. It predicts *mean plantation height*, not dominant height.
Biomass: equation 2 and Table 7, fitted to 58 P. pinaster plot records
(plantation ages 4-41 years, densities 323-2281 stems/ha). It predicts
aboveground biomass in kg/ha from current density, age since planting and
site index (mean height at 30 years). The published fit is for young planted
afforestations, not established natural or thinned stands.
"""

from __future__ import annotations

import math


SOURCE_DOI = "10.3390/f15050827"
HEIGHT_AGE_RANGE = (4.0, 39.0)
BIOMASS_AGE_RANGE = (4.0, 41.0)
BIOMASS_DENSITY_RANGE = (323.0, 2281.0)
SITE_REFERENCE_AGE = 30.0
HEIGHT_A1 = -0.0111
HEIGHT_A2_M = 65.4742
BIOMASS_A_KG_HA = 111.0
BIOMASS_B0 = 1.34620
BIOMASS_B1_PER_STEM_HA = 0.00023
BIOMASS_B2_PER_M = 0.03543


def project_mean_height(height_m: float, age1_years: float, age2_years: float) -> float:
    """M1 Chapman-Richards ADA projection from one planted-age/mean-height pair."""
    if not all(math.isfinite(x) and x > 0 for x in (height_m, age1_years, age2_years)):
        raise ValueError("height and planting ages must be finite and positive")
    if height_m >= HEIGHT_A2_M:
        raise ValueError("observed mean height exceeds the published curve asymptote")
    log_growth_ratio = math.log(-math.expm1(HEIGHT_A1 * age2_years)) / math.log(
        -math.expm1(HEIGHT_A1 * age1_years))
    return HEIGHT_A2_M * (height_m / HEIGHT_A2_M) ** log_growth_ratio


def site_index_30(mean_height_m: float, age_since_planting: float) -> float:
    return project_mean_height(mean_height_m, age_since_planting, SITE_REFERENCE_AGE)


def aboveground_biomass_kg_ha(age_since_planting: float, trees_per_ha: float,
                               site_index_30_m: float) -> float:
    """Published Eq. 2/Table 7, using current stocking and SI30 of *mean* height."""
    if not all(math.isfinite(x) and x > 0 for x in
               (age_since_planting, trees_per_ha, site_index_30_m)):
        raise ValueError("age, density and site index must be finite and positive")
    exponent = (BIOMASS_B0 + BIOMASS_B1_PER_STEM_HA * trees_per_ha
                + BIOMASS_B2_PER_M * site_index_30_m)
    return BIOMASS_A_KG_HA * age_since_planting ** exponent


if __name__ == "__main__":
    for h in (3.0, 6.0, 12.0, 20.0):
        for t1 in (4.0, 10.0, 20.0, 39.0):
            assert math.isclose(project_mean_height(h, t1, t1), h, rel_tol=1e-12)
            for t2 in (4.0, 10.0, 20.0, 39.0):
                h2 = project_mean_height(h, t1, t2)
                assert math.isclose(project_mean_height(h2, t2, t1), h, rel_tol=1e-12)
                assert (h2 >= h) == (t2 >= t1)
    print("Young Pinaster M1 anchor, inverse and monotonic checks passed")
