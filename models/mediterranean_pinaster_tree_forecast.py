"""Inventory-anchored Mediterranean Pinus pinaster tree-list projection.

Inputs are measured IFN pine trees >=7.5 cm, their plot density, and a future
interval. Outputs use IFN concentric-plot expansion, so projected stems/ha are
an expected future inventory estimate, not a physical tree census. The
calibration applies to roughly 15–18 year remeasurements of existing stands.
Thinning fractions are supplied management scenarios; their timing and effects
on subsequent growth have not been estimated.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Sequence

import numpy as np
from scipy.stats import norm

from ifn_remeasurement_validation import trees_per_ha
from mediterranean_pinaster_diameter import MEDITERRANEAN_PARAMS, project_diameter
from mediterranean_pinaster_ingrowth import (SOURCE_BA_RANGE_M2_HA, SOURCE_INTERVAL_YEARS,
                                            SOURCE_QMD_RANGE_CM, expected_basal_area_m2_ha)
from spain_biomass_equations import pinus_pinaster_biomass


@dataclass(frozen=True)
class InventoryTree:
    diameter_cm: float
    height_m: float
    removal_fraction: float = 0.0


@dataclass(frozen=True)
class InventoryProjection:
    initial_agb_mg_ha: float
    surviving_cohort_agb_mg_ha: float
    surviving_cohort_stems_ha: float
    new_tree_agb_mg_ha: float
    total_agb_mg_ha: float


def expected_ifn_weight(initial_diameter_cm: float, projected_diameter_cm: float,
                        residual_mean_cm: float, residual_sd_cm: float) -> float:
    """Expected IFN expansion weight for a survivor crossing DBH thresholds."""
    if not (math.isfinite(residual_mean_cm) and math.isfinite(residual_sd_cm)
            and residual_sd_cm > 0):
        raise ValueError("DBH residual mean must be finite and SD positive")
    weight = float(trees_per_ha(np.array([7.5]))[0])
    for threshold in (12.5, 22.5, 42.5):
        before = float(trees_per_ha(np.array([threshold - .01]))[0])
        after = float(trees_per_ha(np.array([threshold]))[0])
        probability = (1.0 if initial_diameter_cm >= threshold else
                       float(norm.sf((threshold - projected_diameter_cm - residual_mean_cm) /
                                     residual_sd_cm)))
        weight += (after - before) * probability
    return weight


def project_inventory(trees: Sequence[InventoryTree], years: float, *,
                      annual_natural_hazard: float = 0.0080866,
                      annual_height_increment_m: float = 0.169389,
                      dbh_residual_mean_cm: float = -0.245,
                      dbh_residual_sd_cm: float = 4.014,
                      basal_area_m2_ha: float | None = None,
                      quadratic_mean_diameter_cm: float | None = None,
                      new_agb_per_new_ba_mg_m2: float = 3.535956) -> InventoryProjection:
    """Forecast expected pine AGB; optional ingrowth requires both BA and QMD.

    ``removal_fraction`` is the fraction of each initial tree's represented
    cohort removed by a supplied thinning scenario. No growth response to
    thinning is assumed. Ingrowth uses linear scaling of a published five-year
    expectation, an unvalidated approximation for longer intervals.
    """
    if not trees or not math.isfinite(years) or years < 0:
        raise ValueError("Provide at least one tree and finite nonnegative years")
    if not all(math.isfinite(x) and x >= 0 for x in
               (annual_natural_hazard, annual_height_increment_m, new_agb_per_new_ba_mg_m2)):
        raise ValueError("Rates must be finite and nonnegative")
    if (basal_area_m2_ha is None) != (quadratic_mean_diameter_cm is None):
        raise ValueError("Provide both basal area and quadratic mean diameter for ingrowth")
    if basal_area_m2_ha is not None and not (
        SOURCE_BA_RANGE_M2_HA[0] <= basal_area_m2_ha <= SOURCE_BA_RANGE_M2_HA[1]
        and SOURCE_QMD_RANGE_CM[0] <= quadratic_mean_diameter_cm <= SOURCE_QMD_RANGE_CM[1]):
        raise ValueError("Ingrowth inputs are outside the published BA/QMD range")
    if not math.isfinite(dbh_residual_mean_cm) or not math.isfinite(dbh_residual_sd_cm) or dbh_residual_sd_cm <= 0:
        raise ValueError("DBH residual mean must be finite and SD positive")
    for tree in trees:
        if not (math.isfinite(tree.diameter_cm) and tree.diameter_cm >= 7.5 and
                math.isfinite(tree.height_m) and tree.height_m > 1.3 and
                math.isfinite(tree.removal_fraction) and 0 <= tree.removal_fraction <= 1):
            raise ValueError("Trees require DBH >=7.5 cm, height >1.3 m, and removal fraction in [0,1]")
    initial_density = sum(float(trees_per_ha(np.array([t.diameter_cm]))[0]) for t in trees)
    initial_agb = cohort_agb = cohort_stems = 0.0
    survival = math.exp(-annual_natural_hazard * years)
    for tree in trees:
        initial_weight = float(trees_per_ha(np.array([tree.diameter_cm]))[0])
        initial_agb += initial_weight * pinus_pinaster_biomass(
            tree.diameter_cm, tree.height_m).aboveground_kg / 1000.0
        d4 = project_diameter(tree.diameter_cm, initial_density, years, MEDITERRANEAN_PARAMS)
        h4 = tree.height_m + annual_height_increment_m * years
        weight4 = (initial_weight if years == 0 else
                   expected_ifn_weight(tree.diameter_cm, d4, dbh_residual_mean_cm, dbh_residual_sd_cm))
        retained = survival * (1.0 - tree.removal_fraction)
        cohort_stems += retained * weight4
        cohort_agb += retained * weight4 * pinus_pinaster_biomass(d4, h4).aboveground_kg / 1000.0
    new_agb = 0.0
    if basal_area_m2_ha is not None:
        new_agb = (expected_basal_area_m2_ha(basal_area_m2_ha, quadratic_mean_diameter_cm)
                   * years / SOURCE_INTERVAL_YEARS * new_agb_per_new_ba_mg_m2)
    return InventoryProjection(initial_agb, cohort_agb, cohort_stems, new_agb, cohort_agb + new_agb)
