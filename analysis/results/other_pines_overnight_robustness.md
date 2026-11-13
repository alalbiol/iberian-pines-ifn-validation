# Overnight climate-window robustness audit

Species-specific climate terms were tested on two inventory intervals at the same plot. 
Each province was held out in turn; training included an inventory-period intercept and 
estimated the climate slope from differences among plots within provinces. 
A positive RMSE gain means improvement over the period-only baseline. 
Province bootstrap intervals reflect geographic uncertainty. These are tree-growth diagnostics, 
not validated stand-carbon factors.

| Species | Climate feature | Plots | Province folds | RMSE gain % | 95% province bootstrap | Improved provinces |
|---|---|---:|---:|---:|---:|---:|
| Pinus sylvestris | `hyd_pr` | 1482 | 9 | -0.36 | -1.53 to +0.82 | 5/9 |
| Pinus sylvestris | `sepjun_wb` | 1482 | 9 | -2.29 | -6.97 to +2.00 | 4/9 |
| Pinus sylvestris | `spring_pr` | 1482 | 9 | +0.61 | -0.14 to +1.65 | 4/9 |
| Pinus sylvestris | `hyd_wb` | 1482 | 9 | -0.72 | -2.74 to +1.66 | 5/9 |
| Pinus pinea | `sepjun_wb` | 932 | 11 | -0.24 | -3.04 to +3.50 | 5/11 |
| Pinus pinea | `spring_pr` | 932 | 11 | +0.10 | -0.21 to +0.50 | 6/11 |
| Pinus pinea | `sepjun_pr` | 932 | 11 | +0.32 | -1.79 to +2.71 | 5/11 |
| Pinus pinea | `hyd_pr` | 932 | 11 | -1.38 | -6.29 to +3.71 | 5/11 |
| Pinus halepensis | `sepjun_wb` | 1636 | 9 | +6.15 | +0.79 to +10.68 | 7/9 |
| Pinus halepensis | `sepjun_pr` | 1636 | 9 | +7.10 | +1.64 to +11.60 | 6/9 |
| Pinus halepensis | `hyd_wb` | 1636 | 9 | +5.23 | -0.94 to +10.13 | 7/9 |
| Pinus halepensis | `spring_pr` | 1636 | 9 | +1.11 | -0.43 to +2.08 | 5/9 |
| Pinus nigra | `spring_pr` | 1752 | 8 | +0.33 | -0.00 to +0.73 | 4/8 |
| Pinus nigra | `sepjun_wb` | 1752 | 8 | +2.61 | +1.49 to +3.62 | 6/8 |
| Pinus nigra | `hyd_pr` | 1752 | 8 | +1.21 | +0.74 to +1.72 | 7/8 |
| Pinus nigra | `hyd_wb` | 1752 | 8 | +2.83 | +0.92 to +4.32 | 6/8 |

The stable-density sensitivity repeats the analysis with plots whose median inferred density 
changed by no more than a factor of two between inventory periods. Full fold-level 
and stable-density results are in the companion JSON. Selecting the best row after seeing 
these held-out results would overfit the test set; use a fresh data period or independent 
annual ring series before adopting a factor. Survival and stand biomass are not modeled.
