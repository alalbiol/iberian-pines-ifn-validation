"""Published Spanish dominant-height trajectories for four pine species.

These project an observed stand height from age t1 to t2. Source contexts:
P. sylvestris Galicia (Diéguez-Aranda et al. 2005; 10–55 y);
P. pinea four Spanish regions (Calama et al. 2003; 20–214 y);
P. halepensis Spanish plantations (Montero et al. 2001; age-80 SI);
P. nigra three Spanish regions (Martín-Benito et al. 2008; 30–183 y).
Age units are years. IFN stand-age origin can vary by source method.
"""

from __future__ import annotations

import math


PINEA_REGIONS = {
    "west_andalusia": (3.8290, -0.5199),
    "northern_plateau": (5.2613, -0.2576),
    "central_range": (4.0410, -0.3844),
    "catalonia": (4.0352, -0.4138),
}
REFERENCE_AGES = {"pinus_sylvestris": 40.0, "pinus_pinea": 100.0,
                  "pinus_halepensis": 80.0, "pinus_nigra": 80.0}


def project_height(species: str, height_m: float, age1_years: float,
                   age2_years: float, *, pinea_region: str | None = None) -> float:
    """Project stand dominant height from an observed age-height anchor."""
    if not all(math.isfinite(x) and x > 0 for x in (height_m, age1_years, age2_years)):
        raise ValueError("Height and ages must be finite and positive")
    if species == "pinus_sylvestris":
        a, power = 51.39, 1.277
        denominator = 1.0 - (1.0-a/height_m)*(age1_years/age2_years)**power
        if denominator <= 0:
            raise ValueError("Height/age anchor is outside the Sylvestris equation domain")
        return a/denominator
    if species == "pinus_pinea":
        if pinea_region not in PINEA_REGIONS:
            raise ValueError(f"P. pinea requires a region from {tuple(PINEA_REGIONS)}")
        a, power = PINEA_REGIONS[pinea_region]
        return math.exp(a + (math.log(height_m)-a)*(age2_years/age1_years)**power)
    if species == "pinus_halepensis":
        # Montero's 0.203954 coefficient uses age in decades; 0.0203954
        # is the same rate with age in years. The latter is confirmed in
        # later Spanish plantation work using the original site curve.
        def shape(age: float) -> float:
            return (-math.expm1(-0.0203954*age))**(1.0/1.046295)
        return height_m*shape(age2_years)/shape(age1_years)
    if species == "pinus_nigra":
        b1, exponent = 43.735, 0.380
        log_h1 = math.log(height_m)
        x0 = 0.5*(log_h1 + math.sqrt(log_h1**2 + 4*b1*age1_years**(-exponent)))
        return math.exp(x0-(b1/x0)*age2_years**(-exponent))
    raise ValueError(f"Unsupported pine species {species!r}")


def site_index(species: str, height_m: float, age_years: float,
               *, pinea_region: str | None = None) -> float:
    return project_height(species, height_m, age_years, REFERENCE_AGES[species],
                          pinea_region=pinea_region)
