"""Published Mediterranean Pinus pinaster dominant-height growth.

Bravo-Oviedo, del Río & Montero (2007), Forest Ecology and Management
247:107-119, doi:10.1016/j.foreco.2007.04.034, equations 16-17 and
Table 3. The 2007 paper explicitly supersedes the biased 2004 site curves.

The default parameters pool the four Mediterranean regions. Regional
parameters are available only when the stand's source region is documented;
an IFN province by itself does not establish genetic provenance.
The source stand ages are 25-180 years (Table 1); younger projections are
extrapolations. Age is total age from stem analysis, not years since planting.
"""

from __future__ import annotations

from dataclasses import dataclass
import math


SOURCE_DOI = "10.1016/j.foreco.2007.04.034"
SOURCE_STAND_AGE_RANGE = (25.0, 180.0)
REFERENCE_AGE = 70.0


@dataclass(frozen=True)
class HeightParameters:
    a: float
    b: float
    c: float


# Table 3, Hossfeld model. NPS = Northern Plateau, CMS = Central Mountain
# Range, IMS = Iberian Mountain Range, SES = Segura-Alcaraz Mountains.
PARAMETERS = {
    "pooled": HeightParameters(14.45, 8.418, 1.341),
    "NPS": HeightParameters(21.54, 6.42, 1.271),
    "CMS": HeightParameters(17.99, 8.62, 1.400),
    "IMS": HeightParameters(14.46, 8.60, 1.313),
    "SES": HeightParameters(17.63, 7.86, 1.421),
}


def project_height(h0_m: float, age0_years: float, age_years: float,
                   region: str = "pooled") -> float:
    """Project dominant height from an observed total-age/height pair.

    Implements source equations 16-17 exactly. ``region='pooled'`` uses the
    paper's combined Mediterranean fit, with no claim about seed provenance.
    """
    if not all(math.isfinite(v) and v > 0 for v in (h0_m, age0_years, age_years)):
        raise ValueError("height and both total ages must be finite and positive")
    if region not in PARAMETERS:
        raise ValueError(f"unknown region {region!r}; use {tuple(PARAMETERS)}")
    p = PARAMETERS[region]
    tc0 = age0_years ** p.c
    tc = age_years ** p.c
    exp_b = math.exp(p.b)
    r0 = h0_m - p.a + math.sqrt((h0_m - p.a) ** 2 + 2 * h0_m * exp_b / tc0)
    return h0_m * tc * (tc0 * r0 + exp_b) / (tc0 * (tc * r0 + exp_b))


def site_index_70(h_m: float, total_age_years: float, region: str = "pooled") -> float:
    return project_height(h_m, total_age_years, REFERENCE_AGE, region)


if __name__ == "__main__":
    for region in PARAMETERS:
        for si in (10.0, 15.0, 20.0):
            assert math.isclose(project_height(si, REFERENCE_AGE, REFERENCE_AGE, region), si,
                                rel_tol=1e-12)
            heights = []
            for age in (25.0, 40.0, 100.0, 180.0):
                h = project_height(si, REFERENCE_AGE, age, region)
                heights.append(h)
                assert math.isclose(project_height(h, age, REFERENCE_AGE, region), si,
                                    rel_tol=1e-12)
                assert h > 0
            assert heights == sorted(heights)
    print("Mediterranean pinaster height equation: anchor and inverse checks passed")
