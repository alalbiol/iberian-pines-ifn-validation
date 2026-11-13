"""Research-only age-based stand snapshots for four Spanish pine species.

Combines published dominant-height shape, fit-province site-index and stand
structure defaults, and Ruiz-Peinado tree allometry. Density is fixed at a
snapshot value, not forecast from planting, mortality or management. IFN4
held-out tests in evaluate_other_pine_stands.py show that default AGB is
not yet reliable for P. pinea and P. nigra.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

from other_pine_dominant_height import REFERENCE_AGES, project_height
from other_pine_tree_growth import BIOMASS


@dataclass(frozen=True)
class StandParameters:
    default_site_index_m: float
    default_large_stems_ha: float
    log_diameter_at_10m: float
    log_height_exponent: float
    mean_to_dominant_agb_ratio: float
    minimum_source_age_years: float


PARAMETERS = {
    "pinus_sylvestris": StandParameters(
        15.618229356763319, 668.2385543951713, 3.1909166016988904,
        0.6258782691935689, 0.5572888603184454, 10.0),
    "pinus_pinea": StandParameters(
        17.87054444388174, 381.9718634205488, 3.437236761872603,
        0.8664405265164256, 0.688963027792526, 20.0),
    "pinus_halepensis": StandParameters(
        14.433618761453637, 286.4788975654116, 3.048335455882719,
        0.8799265533216261, 0.8362161670789326, 10.0),
    "pinus_nigra": StandParameters(
        19.60719525778724, 962.0032115776784, 3.098583845763465,
        0.7012340221040557, 0.5254403601363108, 30.0),
}

# Lower quartile, median and upper quartile of the same fit-province stands that
# produced PARAMETERS above. The medians are the defaults by construction, and
# evaluate_other_pine_stands.py checks all three against a fresh fit. They give
# the demo a site and density range taken from measured stands instead of an
# invented spread; this is a sample spread, not a confidence interval.
SITE_INDEX_QUARTILES = {
    "pinus_sylvestris": (13.503088992869488, 15.618229356763319, 19.293980038528133),
    "pinus_pinea": (16.33198643243524, 17.87054444388174, 18.9173755028416),
    "pinus_halepensis": (9.906806541131019, 14.433618761453637, 16.564798505194183),
    "pinus_nigra": (18.611766878894787, 19.60719525778724, 23.37505958840828),
}
LARGE_STEMS_QUARTILES = {
    "pinus_sylvestris": (477.7124036316068, 668.2385543951713, 1021.9162057104942),
    "pinus_pinea": (172.31175172082536, 381.9718634205488, 666.6823727293838),
    "pinus_halepensis": (151.1971959373006, 286.4788975654116, 534.0532534861377),
    "pinus_nigra": (767.4805033542509, 962.0032115776784, 1133.53687246561),
}


@dataclass(frozen=True)
class StandSnapshot:
    dominant_height_m: float
    dominant_diameter_cm: float
    large_stems_ha: float
    aboveground_biomass_mg_ha: float


def estimate_stand(species: str, age_years: float, *,
                   pinea_region: str | None = None,
                   site_index_m: float | None = None,
                   dominant_height_m: float | None = None,
                   large_stems_ha: float | None = None) -> StandSnapshot:
    """Estimate an even-aged stand snapshot within the tested age domain."""
    if species not in PARAMETERS:
        raise ValueError(f"Unsupported species {species!r}")
    p = PARAMETERS[species]
    if not math.isfinite(age_years) or age_years < p.minimum_source_age_years:
        raise ValueError(f"{species} requires age >= {p.minimum_source_age_years:g} years")
    h = dominant_height_m
    if h is None:
        si = p.default_site_index_m if site_index_m is None else site_index_m
        if not math.isfinite(si) or si <= 0:
            raise ValueError("Site index must be finite and positive")
        h = project_height(species, si, REFERENCE_AGES[species], age_years,
                           pinea_region=pinea_region)
    n = p.default_large_stems_ha if large_stems_ha is None else large_stems_ha
    if not (math.isfinite(h) and h > 1.3 and math.isfinite(n) and n > 0):
        raise ValueError("Dominant height and large-stem density must be positive")
    d = math.exp(p.log_diameter_at_10m + p.log_height_exponent*math.log(h/10.0))
    agb = n*p.mean_to_dominant_agb_ratio*BIOMASS[species](d,h).aboveground_kg/1000.0
    return StandSnapshot(h,d,n,agb)
