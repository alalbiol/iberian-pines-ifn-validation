"""Recompute descriptive held-out metrics from the distributed plot forecasts."""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'analysis/results/plot_level'
OUT = ROOT / 'outputs'


def summarize(frame, observed, predictions):
    rows = []
    for case, group in frame.groupby('case', sort=False):
        if group.plot_id.isna().any() or group.plot_id.duplicated().any():
            raise ValueError(f'Missing or duplicate plot IDs for {case}')
        for name, column in predictions.items():
            truth = group[observed].to_numpy(dtype=float)
            prediction = group[column].to_numpy(dtype=float)
            if not np.isfinite(truth).all() or not np.isfinite(prediction).all():
                raise ValueError(f'Non-finite values for {case}: {name}')
            error = prediction - truth
            rows.append(dict(case=case, model=name, n=len(group),
                             rmse_mg_ha=np.sqrt(np.mean(error**2)),
                             bias_mg_ha=error.mean(), mae_mg_ha=np.abs(error).mean()))
    return pd.DataFrame(rows)


def main():
    OUT.mkdir(exist_ok=True)
    jobs = [('age_only', 'observed', {'default': 'default', 'measured_inputs': 'measured_inputs',
                                     'fit_median_reference': 'reference'}),
            ('measured_tree', 'agb4', {'no_change': 'agb3', 'without_entry': 'pred_no_entry',
                                      'with_entry': 'pred_total'})]
    for name, observed, predictions in jobs:
        frame = pd.read_csv(DATA / f'{name}_forecasts.csv', dtype={'plot_id': str, 'province_code': str})
        result = summarize(frame, observed, predictions)
        result.to_csv(OUT / f'{name}_metrics.csv', index=False)
        print(result.to_string(index=False, float_format=lambda v: f'{v:.3f}'))
        print()


if __name__ == '__main__':
    main()
