"""Provisional Mediterranean P. pinaster stand-structure calibration.

Calibrated on pine-dominated IFN4 plots from Madrid and Salamanca by
fit_mediterranean_structure.py. Held-out central provinces show dominant
height-to-diameter RMSE 5.10 cm (949 plots) and even-aged mean/dominant
tree biomass-ratio RMSE 0.081 (225 plots). These are empirical relations,
not published physiological laws. They apply to trees >=7.5 cm and should
not be extrapolated to seedlings below breast height.
"""

from __future__ import annotations

import math


HD_SCALE = 6.932591731737465
HD_EXPONENT = 0.631322281996406
HD_DENSITY_EXPONENT = -0.11555670416574702
RATIO_LOGIT_500 = 0.4181079317972758
RATIO_LOG_DENSITY_SLOPE = -0.8883217291048581


def dominant_diameter_cm(dominant_height_m: float, trees_per_ha: float) -> float:
    if not math.isfinite(dominant_height_m) or dominant_height_m <= 1.3:
        raise ValueError("dominant height must be finite and above breast height")
    if not math.isfinite(trees_per_ha) or trees_per_ha <= 0:
        raise ValueError("trees per hectare must be finite and positive")
    return (HD_SCALE * (dominant_height_m - 1.3) ** HD_EXPONENT
            * (trees_per_ha / 500.0) ** HD_DENSITY_EXPONENT)


def mean_to_dominant_agb_ratio(trees_per_ha: float) -> float:
    if not math.isfinite(trees_per_ha) or trees_per_ha <= 0:
        raise ValueError("trees per hectare must be finite and positive")
    z = RATIO_LOGIT_500 + RATIO_LOG_DENSITY_SLOPE * math.log(trees_per_ha / 500.0)
    return 1.0 / (1.0 + math.exp(-z))
