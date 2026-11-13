"""Inventory-anchored AGB forecasts for four Spanish pines.

Requires measured IFN sample trees >=7.5 cm and initial all-species plot
density. The model projects expected future IFN sample expansion and pine
AGB over a 5–25 year interval. The future cut fraction of each tree can be
supplied as a management scenario. New-tree AGB is a fitted constant annual
amount and does not respond to cutting; it is not a regeneration model.
This is not a planting-to-harvest simulator.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Sequence

import numpy as np

from ifn_remeasurement_validation import trees_per_ha
from mediterranean_pinaster_tree_forecast import expected_ifn_weight
from other_pine_tree_growth import BIOMASS, project_tree


@dataclass(frozen=True)
class InventoryTree:
    diameter_cm: float
    height_m: float
    removal_fraction: float = 0.0


@dataclass(frozen=True)
class InventoryParameters:
    annual_natural_hazard: float
    annual_new_tree_agb_mg_ha: float
    dbh_residual_mean_cm: float
    dbh_residual_sd_cm: float


PARAMETERS = {
    "pinus_sylvestris": InventoryParameters(0.006076947240821352, 1.4345956035984921,
                                            0.112, 2.2201),
    "pinus_pinea": InventoryParameters(0.004730952759216419, 1.1871946061958825,
                                       -0.4276, 2.787),
    "pinus_halepensis": InventoryParameters(0.005456223561553413, 0.7007679214458313,
                                            -0.1565, 2.7597),
    "pinus_nigra": InventoryParameters(0.002187010945228589, 0.8580391608084003,
                                       -0.2381, 1.9577),
}


@dataclass(frozen=True)
class InventoryProjection:
    initial_agb_mg_ha: float
    surviving_cohort_agb_mg_ha: float
    surviving_cohort_stems_ha: float
    expected_new_tree_agb_mg_ha: float
    total_agb_mg_ha: float


def project_inventory(species: str, trees: Sequence[InventoryTree],
                      initial_all_species_density_ha: float, years: float,
                      *, annual_new_tree_agb_mg_ha: float | None = None) -> InventoryProjection:
    """Forecast AGB conditional on an initial inventory and management inputs."""
    if species not in PARAMETERS:
        raise ValueError(f"Unsupported species {species!r}")
    if not trees or not (math.isfinite(initial_all_species_density_ha) and
                         initial_all_species_density_ha > 0):
        raise ValueError("Provide trees and positive finite all-species density")
    if not math.isfinite(years) or not 0 <= years <= 25:
        raise ValueError("Years must be finite and between 0 and 25")
    p = PARAMETERS[species]
    new_rate = p.annual_new_tree_agb_mg_ha if annual_new_tree_agb_mg_ha is None else annual_new_tree_agb_mg_ha
    if not math.isfinite(new_rate) or new_rate < 0:
        raise ValueError("Annual new-tree AGB must be finite and nonnegative")
    survival = math.exp(-p.annual_natural_hazard*years)
    initial_agb = cohort_agb = cohort_stems = 0.0
    for tree in trees:
        if not (math.isfinite(tree.diameter_cm) and tree.diameter_cm >= 7.5 and
                math.isfinite(tree.height_m) and tree.height_m > 1.3 and
                math.isfinite(tree.removal_fraction) and 0 <= tree.removal_fraction <= 1):
            raise ValueError("Each tree needs DBH >=7.5 cm, height >1.3 m, removal fraction in [0,1]")
        w3 = float(trees_per_ha(np.array([tree.diameter_cm]))[0])
        initial_agb += w3*BIOMASS[species](tree.diameter_cm,tree.height_m).aboveground_kg/1000.0
        future = project_tree(species,tree.diameter_cm,tree.height_m,
                              initial_all_species_density_ha,years)
        w4 = (w3 if years == 0 else expected_ifn_weight(
            tree.diameter_cm,future.diameter_cm,p.dbh_residual_mean_cm,p.dbh_residual_sd_cm))
        retained = survival*(1.0-tree.removal_fraction)
        cohort_agb += retained*w4*future.aboveground_biomass_kg/1000.0
        cohort_stems += retained*w4
    entry = years*new_rate
    return InventoryProjection(initial_agb,cohort_agb,cohort_stems,entry,cohort_agb+entry)
