"""Measured-tree growth components for four Spanish pine species.

These coefficients were fitted on strict IFN3→IFN4 surviving-tree pairs with
province-held-out tests in fit_other_pine_tree_growth.py. They require a
measured initial tree and plot density; they do not supply planting growth,
mortality, recruitment, or management. Diameter is cm, height is m, and
density is all-species live trees/ha at the initial inventory.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

from spain_biomass_equations import (pinus_halepensis_biomass, pinus_nigra_biomass,
                                     pinus_pinea_biomass, pinus_sylvestris_biomass)


@dataclass(frozen=True)
class PineGrowthParameters:
    log_dbh_rate: float
    dbh_size_exponent: float
    dbh_density_exponent: float
    log_height_rate: float
    height_size_exponent: float
    height_dbh_exponent: float
    height_density_exponent: float
    flat_height_rate_m_year: float
    use_flat_height: bool = False


PARAMETERS = {
    "pinus_sylvestris": PineGrowthParameters(
        -1.4182876662817079, -0.03624704379507389, -0.1710713850316293,
        -2.1531796530652496, -0.1322694616859551, -0.22913629332625024,
        0.05130826840236677, 0.1009407043048469),
    "pinus_pinea": PineGrowthParameters(
        -0.981684486913155, -0.05133662241657413, -0.03540891754368781,
        -1.881898680027767, -0.016127655835421643, -0.2578094726548059,
        0.14961384272780695, 0.12855260582930056, True),
    "pinus_halepensis": PineGrowthParameters(
        -1.3288614663739589, 0.11103420757867741, -0.12032793811474538,
        -2.243302512177598, 0.373685375207522, -0.5358422022468806,
        -0.09724066370905196, 0.09266944545125047),
    "pinus_nigra": PineGrowthParameters(
        -1.6846908355090597, -0.16823265415061028, -0.13226239593037942,
        -2.408835100664844, 0.2115540147655431, -0.6070239985463675,
        -0.011225024936819287, 0.076751802782552),
}

BIOMASS = {
    "pinus_sylvestris": pinus_sylvestris_biomass,
    "pinus_pinea": pinus_pinea_biomass,
    "pinus_halepensis": pinus_halepensis_biomass,
    "pinus_nigra": pinus_nigra_biomass,
}


@dataclass(frozen=True)
class TreeProjection:
    diameter_cm: float
    height_m: float
    aboveground_biomass_kg: float


def project_tree(species: str, diameter_cm: float, height_m: float,
                 initial_density_ha: float, years: float) -> TreeProjection:
    """Project a surviving measured tree over a nonnegative interval."""
    if species not in PARAMETERS:
        raise ValueError(f"Unknown species {species!r}; choose {tuple(PARAMETERS)}")
    if not all(math.isfinite(v) for v in (diameter_cm, height_m, initial_density_ha, years)):
        raise ValueError("All inputs must be finite")
    if diameter_cm < 7.5 or height_m <= 1.3 or initial_density_ha <= 0 or years < 0:
        raise ValueError("DBH >=7.5 cm, height >1.3 m, density >0, years >=0 required")
    p = PARAMETERS[species]
    dbh_rate = math.exp(p.log_dbh_rate + p.dbh_size_exponent*math.log(diameter_cm/20.0)
                        + p.dbh_density_exponent*math.log(initial_density_ha/500.0))
    height_rate = (p.flat_height_rate_m_year if p.use_flat_height else
                   math.exp(p.log_height_rate + p.height_size_exponent*math.log(height_m/10.0)
                            + p.height_dbh_exponent*math.log(diameter_cm/20.0)
                            + p.height_density_exponent*math.log(initial_density_ha/500.0)))
    d4 = diameter_cm + years*dbh_rate
    h4 = height_m + years*height_rate
    return TreeProjection(d4, h4, BIOMASS[species](d4, h4).aboveground_kg)
