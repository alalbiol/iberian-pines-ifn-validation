"""Exploratory age-independent diameter increment for Mediterranean P. pinaster.

Fit on strict IFN3→IFN4 linked trees from Madrid and Salamanca only; see
fit_mediterranean_diameter.py for the reproducible calibration and external
province test. This estimates survivor-tree DBH change, not stand mortality
or a complete stand biomass trajectory.
"""

from __future__ import annotations

from dataclasses import dataclass
import math


@dataclass(frozen=True)
class DiameterGrowthParameters:
    log_rate_at_20cm_500ha: float
    log_diameter_exponent: float
    log_density_exponent: float


# Frozen Madrid/Salamanca exploratory fit from fit_mediterranean_diameter.py.
MEDITERRANEAN_PARAMS = DiameterGrowthParameters(
    -0.8060244090159404, -0.2968745445571137, -0.09069442208159328)


def annual_increment_cm(diameter_cm: float, trees_per_ha: float,
                        params: DiameterGrowthParameters) -> float:
    if not (math.isfinite(diameter_cm) and math.isfinite(trees_per_ha)):
        raise ValueError("diameter and density must be finite")
    if diameter_cm <= 0 or trees_per_ha <= 0:
        raise ValueError("diameter and density must be positive")
    return math.exp(params.log_rate_at_20cm_500ha
                    + params.log_diameter_exponent * math.log(diameter_cm / 20.0)
                    + params.log_density_exponent * math.log(trees_per_ha / 500.0))


def project_diameter(diameter_cm: float, trees_per_ha: float, years: float,
                     params: DiameterGrowthParameters) -> float:
    if not math.isfinite(years) or years < 0:
        raise ValueError("years must be finite and nonnegative")
    return diameter_cm + years * annual_increment_cm(diameter_cm, trees_per_ha, params)
