"""Published young-afforestation height and biomass models for four pines.

Calama et al. (2024), Forests 15:827, doi:10.3390/f15050827,
Tables 3, 5 and 7. Height is plantation *mean* height; age is years since
planting. Biomass is kg dry aboveground matter/ha. The original plots were
largely unthinned afforestations and use root-collar-based tree allometry.
"""

from __future__ import annotations

from dataclasses import dataclass
import math


@dataclass(frozen=True)
class YoungParameters:
    height_model: str
    height_a1: float
    height_a2: float
    biomass_a: float
    biomass_b: float
    biomass_density_b: float
    biomass_si_b: float
    age_range: tuple[float, float]
    density_range: tuple[float, float]


PARAMETERS = {
    "pinus_sylvestris": YoungParameters("M1", -0.0070, 75.0863,
        1156.3, 0.78544, 0.00006, 0.03286, (5,38), (504,2900)),
    "pinus_pinea": YoungParameters("M6", 0.0022, -3.6977,
        250.9, 0.77303, 0.00031, 0.08417, (1,35), (206,2375)),
    "pinus_halepensis": YoungParameters("M2", 0.0759, 0.6520,
        52.2, 1.39838, 0.00025, 0.05534, (2,35), (178,2375)),
    "pinus_nigra": YoungParameters("M1", -0.0120, 52.2884,
        233.8, 1.01399, 0.00021, 0.05138, (3,44), (380,2083)),
}
REFERENCE_AGE_YEARS = 30.0

# Calama et al. (2024), Table 1: mean age, plantation mean height, and living
# stems/ha across sampled young afforestations. Illustrative defaults only;
# they are not measurements of the parcel being simulated.
SAMPLE_MEANS = {
    "pinus_sylvestris": (19.0, 6.5, 1443.0),
    "pinus_pinea": (15.0, 3.2, 884.0),
    "pinus_halepensis": (15.0, 3.5, 862.0),
    "pinus_nigra": (20.0, 6.1, 1208.0),
}


def project_mean_height(species: str, height_m: float,
                        age1_years: float, age2_years: float) -> float:
    """Project an observed plantation mean-height anchor within source ages."""
    if species not in PARAMETERS:
        raise ValueError(f"Unsupported species {species!r}")
    p = PARAMETERS[species]
    if not (math.isfinite(height_m) and 0 < height_m and
            math.isfinite(age1_years) and math.isfinite(age2_years) and
            age1_years > 0 and age2_years > 0):
        raise ValueError("Height and ages must be finite and positive")
    if not (p.age_range[0] <= age1_years <= p.age_range[1] and
            p.age_range[0] <= age2_years <= p.age_range[1]):
        raise ValueError(f"Ages must be within the published {p.age_range} year range")
    if p.height_model == "M1":
        if height_m >= p.height_a2:
            raise ValueError("Observed height exceeds source curve asymptote")
        log_ratio = math.log(-math.expm1(p.height_a1*age2_years)) / math.log(
            -math.expm1(p.height_a1*age1_years))
        return p.height_a2*(height_m/p.height_a2)**log_ratio
    if p.height_model == "M2":
        shape2 = -math.expm1(-p.height_a1*age2_years)
        shape1 = -math.expm1(-p.height_a1*age1_years)
        return height_m*(shape2/shape1)**p.height_a2
    if p.height_model == "M6":
        f1 = math.log(-math.expm1(-p.height_a1*age1_years))
        x0 = (math.log(height_m)-p.height_a2*f1)/(1.0+f1)
        base2 = -math.expm1(-p.height_a1*age2_years)
        return math.exp(x0)*base2**(p.height_a2+x0)
    raise AssertionError("Unknown published height model")


def site_index_30(species: str, mean_height_m: float, age_since_planting: float) -> float:
    return project_mean_height(species,mean_height_m,age_since_planting,REFERENCE_AGE_YEARS)


def aboveground_biomass_kg_ha(species: str, age_since_planting: float,
                              stems_ha: float, site_index_30_m: float) -> float:
    """Source Eq. 2: a × age^(b + b1×density + b2×SI30)."""
    if species not in PARAMETERS:
        raise ValueError(f"Unsupported species {species!r}")
    if not all(math.isfinite(x) and x > 0 for x in
               (age_since_planting,stems_ha,site_index_30_m)):
        raise ValueError("Age, density and site index must be finite and positive")
    p = PARAMETERS[species]
    if not p.age_range[0] <= age_since_planting <= p.age_range[1]:
        raise ValueError(f"Age must be within the published {p.age_range} year range")
    if not p.density_range[0] <= stems_ha <= p.density_range[1]:
        raise ValueError(f"Density must be within the published {p.density_range} stems/ha range")
    exponent = p.biomass_b+p.biomass_density_b*stems_ha+p.biomass_si_b*site_index_30_m
    return p.biomass_a*age_since_planting**exponent
