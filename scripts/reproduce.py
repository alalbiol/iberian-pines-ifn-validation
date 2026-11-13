"""Reproduce and check the distributed metrics, then regenerate the figures."""
import argparse
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import subprocess
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
FIGURES = ROOT / 'analysis/figures'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--metrics-only', action='store_true')
    args = parser.parse_args()
    report = ROOT / 'outputs/reproduction_check.json'
    report.unlink(missing_ok=True)
    subprocess.run([sys.executable, str(ROOT / 'scripts/summarize_results.py')], check=True)
    for name in ['age_only_metrics.csv', 'measured_tree_metrics.csv']:
        actual = pd.read_csv(ROOT / 'outputs' / name)
        expected = pd.read_csv(ROOT / 'expected' / name)
        pd.testing.assert_frame_equal(actual, expected, rtol=1e-9, atol=1e-9)
    if not args.metrics_only:
        env = dict(os.environ, MPLBACKEND='Agg')
        subprocess.run([sys.executable, str(FIGURES / 'make_figures.py')], env=env, check=True)
        for name in ['fig1_design', 'fig2_age_only', 'fig3_measured_tree', 'fig4_bias_by_province']:
            for suffix in ['.pdf', '.png']:
                if not (FIGURES / (name + suffix)).is_file():
                    raise RuntimeError(f'Missing generated figure: {name}{suffix}')
    report.write_text(json.dumps({
        'status': 'passed', 'python': platform.python_version(),
        'metrics_match_reference': True, 'rtol': 1e-9, 'atol': 1e-9,
        'figures_regenerated': not args.metrics_only,
        'dependencies': {name: importlib.metadata.version(name) for name in
                         ['numpy', 'pandas', 'scipy', 'matplotlib', 'geopandas', 'dbfread']},
    }, indent=2) + '\n')
    print(f'PASS: reference metrics reproduced. Report: {report}')


if __name__ == '__main__':
    main()
