"""Published five-year basal-area ingrowth for central-Spain Pinus pinaster.

Bravo et al. (2008), Forest Systems 17:250-260,
doi:10.5424/srf/2008173-01039, Table 2. The source uses backdated
five-year changes in 88 plots, and defines ingrowth in the concentric plot
sample, including diameter-threshold crossings. Output is new basal area
in m2/ha, not stem count or a full density trajectory. The source model has
not been externally validated for IFN3→IFN4 intervals.
"""

from __future__ import annotations

import math


SOURCE_BA_RANGE_M2_HA = (13.25, 55.59)
SOURCE_QMD_RANGE_CM = (9.64, 43.49)
SOURCE_INTERVAL_YEARS = 5.0
OCCURRENCE_THRESHOLD = 0.38


def occurrence_probability(basal_area_m2_ha: float, quadratic_mean_diameter_cm: float) -> float:
    """Probability of any ingrowth over five years, source Table 2 logit."""
    if not all(math.isfinite(x) and x > 0 for x in (basal_area_m2_ha, quadratic_mean_diameter_cm)):
        raise ValueError("basal area and quadratic mean diameter must be finite and positive")
    z = 12.3424 + 0.1108 * basal_area_m2_ha - 0.6154 * quadratic_mean_diameter_cm
    return 1.0 / (1.0 + math.exp(-z)) if z >= 0 else math.exp(z) / (1.0 + math.exp(z))


def conditional_basal_area_m2_ha(quadratic_mean_diameter_cm: float) -> float:
    """New basal area, conditional on ingrowth, from source Table 2."""
    if not math.isfinite(quadratic_mean_diameter_cm) or quadratic_mean_diameter_cm <= 0:
        raise ValueError("quadratic mean diameter must be finite and positive")
    return max(0.0, 6.7389 - 0.2235 * quadratic_mean_diameter_cm)


def expected_basal_area_m2_ha(basal_area_m2_ha: float, quadratic_mean_diameter_cm: float) -> float:
    """Continuous expectation of source two-stage five-year model."""
    return (occurrence_probability(basal_area_m2_ha, quadratic_mean_diameter_cm)
            * conditional_basal_area_m2_ha(quadratic_mean_diameter_cm))


def threshold_basal_area_m2_ha(basal_area_m2_ha: float, quadratic_mean_diameter_cm: float) -> float:
    """Paper's deterministic two-stage output using its 0.38 cutoff."""
    p = occurrence_probability(basal_area_m2_ha, quadratic_mean_diameter_cm)
    return conditional_basal_area_m2_ha(quadratic_mean_diameter_cm) if p >= OCCURRENCE_THRESHOLD else 0.0


if __name__ == "__main__":
    assert math.isclose(conditional_basal_area_m2_ha(26.28), 0.86532)
    assert 0 < occurrence_probability(36.14, 26.28) < 1
    assert expected_basal_area_m2_ha(36.14, 26.28) > 0
    print("Published Pinaster five-year ingrowth equations passed basic checks")
