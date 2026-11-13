# Overnight climate-window robustness audit

Species-specific climate terms were tested on two inventory intervals at the same plot. 
Each province was held out in turn; training included an inventory-period intercept and 
estimated the climate slope from differences among plots within provinces. 
A positive RMSE gain means improvement over the period-only baseline. 
Province bootstrap intervals reflect geographic uncertainty. These are tree-growth diagnostics, 
not validated stand-carbon factors.

| Species | Climate feature | Plots | Province folds | RMSE gain % | 95% province bootstrap | Improved provinces |
|---|---|---:|---:|---:|---:|---:|
| Pinus sylvestris | `hyd_pr` | 1482 | 9 | +0.24 | -0.07 to +0.76 | 5/9 |
| Pinus sylvestris | `sepjun_wb` | 1482 | 9 | +0.29 | -0.57 to +1.24 | 5/9 |
| Pinus sylvestris | `spring_pr` | 1482 | 9 | +0.39 | -0.02 to +1.00 | 7/9 |
| Pinus sylvestris | `hyd_wb` | 1482 | 9 | +0.57 | -0.25 to +1.71 | 5/9 |
| Pinus pinea | `sepjun_wb` | 932 | 11 | -0.00 | -0.01 to +0.00 | 6/11 |
| Pinus pinea | `spring_pr` | 932 | 11 | -0.11 | -0.28 to +0.00 | 6/11 |
| Pinus pinea | `sepjun_pr` | 932 | 11 | -0.00 | -0.00 to +0.00 | 5/11 |
| Pinus pinea | `hyd_pr` | 932 | 11 | -0.07 | -0.69 to +0.89 | 6/11 |
| Pinus halepensis | `sepjun_wb` | 1636 | 9 | -0.75 | -1.57 to +0.18 | 3/9 |
| Pinus halepensis | `sepjun_pr` | 1636 | 9 | -1.32 | -2.60 to +0.26 | 2/9 |
| Pinus halepensis | `hyd_wb` | 1636 | 9 | +0.45 | -1.89 to +1.99 | 6/9 |
| Pinus halepensis | `spring_pr` | 1636 | 9 | -0.06 | -0.21 to +0.01 | 2/9 |
| Pinus nigra | `spring_pr` | 1752 | 8 | +0.28 | -0.28 to +0.79 | 4/8 |
| Pinus nigra | `sepjun_wb` | 1752 | 8 | +0.08 | -1.71 to +1.29 | 3/8 |
| Pinus nigra | `hyd_pr` | 1752 | 8 | +0.10 | -1.32 to +0.99 | 4/8 |
| Pinus nigra | `hyd_wb` | 1752 | 8 | +0.29 | -3.73 to +2.53 | 5/8 |

The stable-density sensitivity repeats the analysis with plots whose median inferred density 
changed by no more than a factor of two between inventory periods. Full fold-level 
and stable-density results are in the companion JSON. Selecting the best row after seeing 
these held-out results would overfit the test set; use a fresh data period or independent 
annual ring series before adopting a factor. Survival and stand biomass are not modeled.
