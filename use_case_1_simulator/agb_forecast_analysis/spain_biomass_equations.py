"""
Per-tree biomass equations linking the T3 growth curves (which output height
or diameter) to AGB (which is what our ESA CCI satellite data and the
project's carbon-credit tooling ultimately need).

Source: Ruiz-Peinado, del Rio & Montero — the standard INIA/CIFOR reference
for Spanish national-forest-inventory-compatible biomass equations, same
research group as the T3 growth papers:
  - Pine: Ruiz-Peinado, del Rio & Montero (2011), "New models for estimating
    the carbon sink capacity of Spanish softwood species," Forest Systems
    20(1):176-188. PDF: use_case_1_simulator/docs/ruiz_peinado_2011_softwood_biomass.pdf
  - Oaks: Ruiz-Peinado, Montero & del Rio (2012), "Biomass models to estimate
    carbon stocks for hardwood tree species," Forest Systems 21(1):42-52.
    PDF: use_case_1_simulator/docs/ruiz_peinado_2012_hardwood_biomass.pdf

Both papers fit weight (kg, dry) as a function of diameter at breast height
(d, cm) and/or total height (h, m), per species, per biomass component (stem,
thick/medium/thin branches, roots), using the exact same predictor variables
Spain's National Forest Inventory records — deliberately, so these equations
apply directly to NFI trees, and to ours.

IMPORTANT ASYMMETRY, resolved 2026-09-24 by spain_height_diameter.py
(see sylvestris_species_growth_plan.md, T9):
  - Quercus ilex's equations run almost entirely on diameter alone — a close
    match to our T3 Q. ilex curve, which already outputs diameter directly.
    Height for the thick-branch term: spain_height_diameter.hd_params_dehesa.
  - Pinus pinaster: T3 gives height only; diameter now comes from inverting
    the nationwide h-d model (spain_height_diameter), which fits IFN3 well.
  - Quercus suber: diameter from the dominant-tree diameter-growth model
    (spain_growth_curves.quercus_suber_diameter_project, Sanchez-Gonzalez et al.
    2005). It gives dbh UNDER cork; convert with quercus_suber_diameter_over_cork
    (assumed cork thickness) before calling quercus_suber_biomass.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass
class BiomassComponents:
    stem_kg: float
    thick_branches_kg: float
    medium_branches_kg: float
    thin_branches_foliage_kg: float
    roots_kg: float

    @property
    def aboveground_kg(self) -> float:
        return self.stem_kg + self.thick_branches_kg + self.medium_branches_kg + self.thin_branches_foliage_kg

    @property
    def total_kg(self) -> float:
        return self.aboveground_kg + self.roots_kg


# --- Pinus pinaster (Ruiz-Peinado et al. 2011, Table 3) --------------------
# NEEDS BOTH d (cm) and h (m) — see module docstring's open gap.

def pinus_pinaster_biomass(d_cm: float, h_m: float) -> BiomassComponents:
    stem = 0.0278 * d_cm**2.115 * h_m**0.618
    thick_and_medium = 0.000381 * d_cm**3.141  # paper reports thick+medium as one combined term
    thin_and_foliage = 0.0129 * d_cm**2.320
    roots = 0.00444 * d_cm**2.804
    return BiomassComponents(stem, thick_and_medium, 0.0, thin_and_foliage, roots)


PINUS_PINASTER_ROOT_SHOOT_RATIO = 0.285  # alternative to the roots_kg term above, if only aboveground is known


# --- Other Spanish pines: Ruiz-Peinado et al. (2011), Table 3 -----------
# d in cm, total height h in m; all output components in kg dry matter.
# DOI: 10.5424/fs/2011201-11643. The source sampling DBH minima vary by
# species (6.2–10 cm), so 7.5 cm IFN trees can be below some source minima.

def pinus_halepensis_biomass(d_cm: float, h_m: float) -> BiomassComponents:
    stem = 0.0139 * d_cm**2 * h_m
    thick = 3.926 * (d_cm - 27.5) if d_cm > 27.5 else 0.0
    medium = 4.257 + 0.00506 * d_cm**2 * h_m - 0.0722 * d_cm * h_m
    thin = 6.197 + 0.00932 * d_cm**2 * h_m - 0.0686 * d_cm * h_m
    roots = 0.0785 * d_cm**2
    return BiomassComponents(stem, thick, medium, thin, roots)


def pinus_nigra_biomass(d_cm: float, h_m: float) -> BiomassComponents:
    stem = 0.0403 * d_cm**1.838 * h_m**0.945
    thick = 0.228 * (d_cm - 32.5)**2 if d_cm > 32.5 else 0.0
    medium = 0.0521 * d_cm**2
    thin = 0.0720 * d_cm**2
    roots = 0.0189 * d_cm**2.445
    return BiomassComponents(stem, thick, medium, thin, roots)


def pinus_pinea_biomass(d_cm: float, h_m: float) -> BiomassComponents:
    stem = 0.0224 * d_cm**1.923 * h_m**1.0193
    thick = 0.247 * (d_cm - 22.5)**2 if d_cm > 22.5 else 0.0
    medium = 0.0525 * d_cm**2
    thin = 21.927 + 0.0707 * d_cm**2 - 2.827 * h_m
    roots = 0.117 * d_cm**2
    return BiomassComponents(stem, thick, medium, thin, roots)


def pinus_sylvestris_biomass(d_cm: float, h_m: float) -> BiomassComponents:
    stem = 0.0154 * d_cm**2 * h_m
    thick = ((0.540 - 0.0119 * h_m) * (d_cm - 37.5)**2
             if d_cm > 37.5 else 0.0)
    medium = 0.0295 * d_cm**2.742 * h_m**-0.899
    thin = 0.530 * d_cm**2.199 * h_m**-1.153
    roots = 0.130 * d_cm**2
    return BiomassComponents(stem, thick, medium, thin, roots)


# --- Quercus ilex (Ruiz-Peinado et al. 2012, Table 2) ----------------------
# Runs almost entirely on d (cm) alone — matches our T3 Q. ilex curve output
# directly. Height only matters for the thick-branch term, and only once
# d > 12.5 cm.

def quercus_ilex_biomass(d_cm: float, h_m: float | None = None) -> BiomassComponents:
    stem = 0.143 * d_cm**2
    if d_cm > 12.5 and h_m is not None:
        thick = 0.0684 * (d_cm - 12.5) ** 2 * h_m
    elif d_cm > 12.5 and h_m is None:
        raise ValueError("d > 12.5 cm needs h_m for the thick-branch term (Ruiz-Peinado et al. 2012)")
    else:
        thick = 0.0
    medium = 0.0898 * d_cm**2
    thin_and_foliage = 0.0824 * d_cm**2
    roots = 0.254 * d_cm**2
    return BiomassComponents(stem, thick, medium, thin_and_foliage, roots)


# --- Quercus suber (Ruiz-Peinado et al. 2012, Table 2) ---------------------
# NEEDS BOTH d (cm) and h (m) — see module docstring's open gap.

def quercus_suber_biomass(d_cm: float, h_m: float) -> BiomassComponents:
    stem = 0.00525 * d_cm**2 * h_m + 0.278 * d_cm * h_m
    thick = 0.0135 * d_cm**2 * h_m
    medium = 0.127 * d_cm * h_m
    thin_and_foliage = 0.0463 * d_cm * h_m
    roots = 0.0829 * d_cm**2
    return BiomassComponents(stem, thick, medium, thin_and_foliage, roots)


def agb_per_hectare_mg(biomass_kg_per_tree: float, trees_per_ha: float) -> float:
    """kg/tree x trees/ha -> Mg/ha, the unit ESA CCI AGB reports in — lets our
    growth-model output be sanity-checked directly against the satellite data
    already collected in spain_forecast_plan.md's S1."""
    return biomass_kg_per_tree * trees_per_ha / 1000.0


if __name__ == "__main__":
    # Sanity checks against the paper's own reported sample ranges (Table 1 /
    # sample description), not a real validation — just catching transcription
    # errors and unit mistakes, same spirit as spain_growth_curves.py's checks.
    checks = []

    # Pine sample: dbh 6.1-64.0 cm, height 5-20.7 m, aboveground biomass 5-1,521 kg.
    b = pinus_pinaster_biomass(d_cm=19.7, h_m=11.6)  # paper's own reported MEAN d, h
    checks.append(("pine at mean sampled size stays within sampled aboveground range (5-1521 kg)", 5 <= b.aboveground_kg <= 1521, b.aboveground_kg))

    # Q. ilex sample: dbh 7.8-85.9 cm, aboveground biomass max 2,506 kg, mean d=40.9.
    b = quercus_ilex_biomass(d_cm=40.9, h_m=7.5)
    checks.append(("Q. ilex at mean sampled diameter stays under sampled max aboveground (2506 kg)", b.aboveground_kg < 2506, b.aboveground_kg))

    # Q. suber sample: dbh 10.5-83.0 cm, height 3.5-14.0 m, aboveground biomass max 2,412 kg, mean d=36.5, h=8.7.
    b = quercus_suber_biomass(d_cm=36.5, h_m=8.7)
    checks.append(("Q. suber at mean sampled size stays under sampled max aboveground (2412 kg)", b.aboveground_kg < 2412, b.aboveground_kg))

    # Unit sanity: a mature pine stand, ~400 trees/ha, ~150 kg/tree aboveground -> should land in a
    # plausible forest AGB range (tens to low hundreds of Mg/ha), not off by orders of magnitude.
    per_ha = agb_per_hectare_mg(150.0, 400.0)
    checks.append(("plausible stand-level AGB order of magnitude (10-300 Mg/ha)", 10 <= per_ha <= 300, per_ha))

    for name, ok, value in checks:
        print(f"{'OK  ' if ok else 'FAIL'} {name}: {value:.1f}")
    if not all(ok for _, ok, _ in checks):
        raise SystemExit(1)
