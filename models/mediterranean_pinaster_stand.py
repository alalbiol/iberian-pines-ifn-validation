"""Provisional mature Mediterranean Pinus pinaster stand AGB chain.

Uses Bravo-Oviedo et al. (2007) dominant height, IFN4 Mediterranean
height-to-diameter and mean/dominant ratio fits, and Ruiz-Peinado pine
allometry. The default site index and standing density are medians of 22
known-age Madrid/Salamanca fit plots. Constant standing density is a
snapshot assumption, not a mortality or management model. Use only for
the published mature-stand age domain (25-180 total years); this is not a
validated planting-to-harvest trajectory.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

from mediterranean_pinaster import REFERENCE_AGE, SOURCE_STAND_AGE_RANGE, project_height
from mediterranean_pinaster_structure import dominant_diameter_cm, mean_to_dominant_agb_ratio
from spain_biomass_equations import pinus_pinaster_biomass


DEFAULT_SITE_INDEX_70_M = 18.711
DEFAULT_LARGE_STEMS_HA = 381.547449


@dataclass(frozen=True)
class StandEstimate:
    total_age_years: float
    dominant_height_m: float
    dominant_diameter_cm: float
    large_stems_per_ha: float
    mean_to_dominant_ratio: float
    aboveground_biomass_mg_ha: float


def estimate_stand_agb(total_age_years: float, *, dominant_height_m: float | None = None,
                       large_stems_per_ha: float | None = None,
                       site_index_70_m: float = DEFAULT_SITE_INDEX_70_M) -> StandEstimate:
    """Estimate AGB at one mature stand age, optionally anchored to field data."""
    if not math.isfinite(total_age_years) or not (SOURCE_STAND_AGE_RANGE[0] <= total_age_years <= SOURCE_STAND_AGE_RANGE[1]):
        raise ValueError(f"total age must be within published stand range {SOURCE_STAND_AGE_RANGE}")
    h = (project_height(site_index_70_m, REFERENCE_AGE, total_age_years)
         if dominant_height_m is None else dominant_height_m)
    n = DEFAULT_LARGE_STEMS_HA if large_stems_per_ha is None else large_stems_per_ha
    if not math.isfinite(n) or n <= 0:
        raise ValueError("large stems per hectare must be finite and positive")
    d = dominant_diameter_cm(h, n)
    ratio = mean_to_dominant_agb_ratio(n)
    agb = n * ratio * pinus_pinaster_biomass(d, h).aboveground_kg / 1000.0
    return StandEstimate(total_age_years, h, d, n, ratio, agb)
