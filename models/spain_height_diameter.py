"""
Height-diameter (h-d) link for the Sylvestris growth model (T9 in
sylvestris_species_growth_plan.md): closes the gap between the T3 growth
curves for Pinus pinaster and Quercus suber (which output dominant HEIGHT
only) and the Ruiz-Peinado biomass equations in spain_biomass_equations.py
(which need BOTH diameter and height for those two species).

Source: Vazquez-Veloso, Yang, Bullock & Bravo (2025), "One model to rule them
all: A nationwide height-diameter model for 91 Spanish forest species,"
Forest Ecology and Management 595, 122981.
https://doi.org/10.1016/j.foreco.2025.122981 (open access, CC BY).
PDF: the original project reference library/1-s2.0-S037811272500489X-main.pdf
Fit on 1.51M h-d pairs from the Spanish NFI (IFN2-4), so it uses the same
dbh-over-bark convention as the Ruiz-Peinado biomass equations. Debarked cork
oaks and resin-tapped pines were excluded from the fit.

Model (their M4, El Mamoun et al. 2013 base form):
    h = 1.3 + a * ln(1 + d) ** b,   d in cm, h in m
with a, b = species value + fixed-effect adjustments for biogeographic region
(Mediterranean is the reference), stand origin (natural is the reference) and
mixture (pure is the reference).

Coefficients are NOT printed in the paper. Table 4 there gives bootstrap
medians for the population-level terms, while the per-species values sit in
the supplementary data. The numbers below are copied from the authors' own
repository (github.com/aitorvv/height-diameter_models_Spain, MIT licence,
tools/data/ab_coefs.csv and fe_coefs.csv), which the published calculator
uses. Their per-species a, b are FULL values, not deviations, and the
region/origin/mixture adjustments are added on top, exactly as in the
authors' tools/one_model_to_rule_them_all.py.

Why the nationwide model and not the two species-specific ones supplied at the
same time: Sanchez-Gonzalez et al. (2007, cork oak) and Castedo-Dorado et al.
(2005, Galician P. pinaster) are both "generalized" h-d models constrained to
pass through (D0, H0), i.e. they need dominant DIAMETER as an input. They
predict the height of other trees in a stand whose D0 is already known, so they
cannot supply D0 from H0, which is the thing we are missing.

HOW WE USE IT, AND THE CAVEATS THAT COME WITH IT:
  - We invert the model: d = exp(((h - 1.3) / a) ** (1 / b)) - 1. That turns the
    T3 dominant height H0 into an estimate of the dominant tree's diameter D0.
  - Inverting a mean-height-given-diameter regression is not the same as a
    mean-diameter-given-height regression (regression dilution): it is the
    diameter at which an average tree reaches height h. Treat it as a central
    estimate, not an unbiased predictor of D0.
  - The model was fit on all tree sizes, not just dominant trees. Dominant trees
    in a closed stand are taller for their diameter than open-grown ones, so
    this inversion probably OVERestimates D0 for dense young stands (d for a
    given h is larger than in a closed stand) and is closer to right for open
    dehesa-type stands. Direction stated, size untested.
  - D0 describes the largest trees. Multiplying per-tree biomass at D0 by the
    full stems/ha overestimates stand AGB; mean-tree biomass needs a separate
    stand-structure assumption (see spain_biomass_equations.agb_per_hectare_mg).
  - The paper reports the model is weakest at extreme (over-mature) sizes.

STATUS BY SPECIES (2026-09-24, checked against IFN3 Caceres/Pontevedra):
  - P. pinaster: nationwide model fits IFN3 well (measured/model height
    1.01-1.22 by diameter class). Inverting H0 -> D0 is usable.
  - Q. ilex: never needs the inverse (T3 curve already gives diameter). For the
    forward model (h for Ruiz-Peinado's thick-branch term) use
    hd_params_dehesa() in dehesa; the nationwide values give ~2x the height.
  - Q. suber: do NOT use this module's inverse. The nationwide inverse is built
    on straight, never-debarked trees, so in dehesa it underestimates D0 by ~2x,
    and the dehesa calibration cannot be inverted (see _DEHESA_AB). Diameter
    comes from the published dominant-tree diameter-growth model instead
    (Sanchez-Gonzalez, Tome & Montero 2005; spain_growth_curves.
    quercus_suber_diameter_project, added 2026-09-24).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

# Population-level fixed-effect adjustments (fe_coefs.csv). Added to the
# species a/b below. Reference levels: Mediterranean, natural, pure.
_FE_A = {
    "mixed": 0.0179742918743663,
    "plantation": -0.0119812862492191,
    "Atlantic": 0.201486366723887,
    "Alpine": -0.0323944527070033,
    "Macaronesian": 0.183767256321656,
}
_FE_B = {
    "mixed": -0.0390335215450275,
    "plantation": 0.0382250816531315,
    "Atlantic": -0.0690993469623317,
    "Alpine": 0.034801026756618,
    "Macaronesian": -0.199457488760704,
}

# Species values (ab_coefs.csv), keyed like species_factors.py / spain_growth_curves.py.
# NFI species codes in comments.
_SPECIES_AB = {
    "pinus_pinaster": (0.41751165848068, 2.67387058532806),  # 26
    "quercus_ilex": (0.761853664036472, 2.16848083211578),  # 45, Q. ilex ssp. ballota
    "quercus_suber": (0.856920106185962, 1.93647716842393),  # 46
}

REGIONS = ("Mediterranean", "Atlantic", "Alpine", "Macaronesian")

# DEHESA calibration — our own fit, not from the paper (2026-09-24). Same M4 form,
# plain nonlinear least squares on IFN3 Caceres (PCMayores, species 045/046, all
# forms, h > 1.3 m, d >= 7.5 cm): Q. ilex n=4058 trees / 1056 plots, Q. suber
# n=1763 / 359 plots. Why it exists: at the same diameter, Caceres oaks are only
# 0.50-0.61x (Q. ilex) and 0.66-0.78x (Q. suber) the nationwide model's height, per
# diameter class. The nationwide fit excluded non-straight trees and debarked
# cork oaks; ~99% of Caceres Q. ilex are IFN form 4/5 and most Caceres Q. suber
# are debarked (ParEsp 11xx). Pine is fine nationally (IFN3/model 1.01-1.22).
# Standard errors ignore plot clustering, so they are optimistic; RMSE 1.3 m
# (Q. ilex) and 1.9 m (Q. suber).
# DO NOT INVERT these for Q. suber (height -> diameter): dehesa height saturates
# (9.9 m at d=60), so the inverse explodes (h=11 m -> d~84 cm, h=13 m -> ~148 cm).
_DEHESA_AB = {
    "quercus_ilex": (0.9275, 1.4179),
    "quercus_suber": (0.9366, 1.5687),
}


def hd_params_dehesa(species: str) -> HDParams:
    """Caceres-dehesa calibration (see _DEHESA_AB). Use for the forward model only."""
    if species not in _DEHESA_AB:
        raise KeyError(f"no dehesa calibration for {species!r}; have {sorted(_DEHESA_AB)}")
    return HDParams(*_DEHESA_AB[species])


@dataclass(frozen=True)
class HDParams:
    a: float
    b: float


def hd_params(species: str, region: str = "Mediterranean", plantation: bool = False, mixed: bool = False) -> HDParams:
    """Species + stand-condition parameters, as in the authors' select_pars()."""
    if species not in _SPECIES_AB:
        raise KeyError(f"no h-d coefficients for {species!r}; have {sorted(_SPECIES_AB)}")
    if region not in REGIONS:
        raise ValueError(f"region must be one of {REGIONS}, got {region!r}")
    a, b = _SPECIES_AB[species]
    for key, on in ((region, region != "Mediterranean"), ("plantation", plantation), ("mixed", mixed)):
        if on:
            a += _FE_A[key]
            b += _FE_B[key]
    return HDParams(a, b)


def height_from_diameter(d_cm: float, params: HDParams) -> float:
    """Forward model: mean total height (m) of a tree with dbh d_cm."""
    if d_cm < 0:
        raise ValueError("d_cm must be >= 0")
    return 1.3 + params.a * math.log(1.0 + d_cm) ** params.b


def diameter_from_height(h_m: float, params: HDParams) -> float:
    """Inverse model: dbh (cm) at which the mean tree reaches height h_m.

    See the module docstring for why this is a central estimate of D0, not an
    unbiased one. Undefined at or below breast height (1.3 m), where dbh does
    not exist, so that raises instead of returning 0.
    """
    if h_m <= 1.3:
        raise ValueError(f"h_m={h_m} m is at or below breast height (1.3 m); dbh is undefined")
    return math.exp(((h_m - 1.3) / params.a) ** (1.0 / params.b)) - 1.0


if __name__ == "__main__":
    checks = []

    # 1. Reproduce the authors' own calculator on their own test data (tools/data/test_hd_df.csv):
    #    P. pinaster, d=20 cm, pure natural stand, four regions. Expected values computed by running
    #    their one_model_to_rule_them_all.py (2026-09-24), not taken from the paper.
    expected = {"Alpine": 9.157558, "Atlantic": 12.549834, "Macaronesian": 10.751531, "Mediterranean": 9.494761}
    for region, h_ref in expected.items():
        h = height_from_diameter(20.0, hd_params("pinus_pinaster", region))
        checks.append((f"matches authors' code, P. pinaster d=20 {region} (ref {h_ref} m)", abs(h - h_ref) < 1e-3, h))

    # 2. The inverse is an exact inverse.
    for sp in _SPECIES_AB:
        p = hd_params(sp, "Atlantic", plantation=True)
        d = diameter_from_height(height_from_diameter(30.0, p), p)
        checks.append((f"{sp}: inverse(forward(30 cm)) == 30 cm", abs(d - 30.0) < 1e-9, d))

    # 3. Paper's reported direction for P. pinaster: Atlantic trees are the tallest for a given dbh.
    hs = {r: height_from_diameter(30.0, hd_params("pinus_pinaster", r)) for r in REGIONS}
    checks.append(("P. pinaster: Atlantic tallest region at d=30 (paper Fig. 3c)", max(hs, key=hs.get) == "Atlantic", hs["Atlantic"]))

    # 4. Independent size check against IFN3 (/media/HD/Data/ifn3, PCMayores): median measured height
    #    in the 20-30 cm diameter class vs the model at the median diameter of that class. Nationwide
    #    model for pine (both provinces) and Pontevedra cork oak; dehesa calibration for Caceres oaks
    #    (the nationwide model fails there by ~2x — see _DEHESA_AB). Tolerance 3.1 m = paper's RMSE.
    ifn3_20_30 = [
        ("pinus_pinaster", hd_params("pinus_pinaster", "Mediterranean"), 12.7),  # Caceres
        ("pinus_pinaster", hd_params("pinus_pinaster", "Atlantic"), 17.0),  # Pontevedra
        ("quercus_suber", hd_params("quercus_suber", "Atlantic"), 11.0),  # Pontevedra
        ("quercus_ilex", hd_params_dehesa("quercus_ilex"), 6.0),  # Caceres
        ("quercus_suber", hd_params_dehesa("quercus_suber"), 7.0),  # Caceres
    ]
    for sp, p, h_ifn in ifn3_20_30:
        h = height_from_diameter(25.0, p)
        checks.append((f"{sp} a={p.a:.3f}: h(25 cm) within 3.1 m of IFN3 20-30 cm median h={h_ifn}", abs(h - h_ifn) <= 3.1, h))

    for name, ok, value in checks:
        print(f"{'OK  ' if ok else 'FAIL'} {name}: {value:.3f}")
    if not all(ok for _, ok, _ in checks):
        raise SystemExit(1)
