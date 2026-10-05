"""Numerical checks of docs/equivalence-derivation.md.

Run from the repository root:  python3 checks/check_equivalence.py
Needs numpy, scipy, pandas, statsmodels (pip install pandas statsmodels).
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

from equivalence import min_passing_delta, PAIRS

rng = np.random.default_rng(0)


def brute_force_power(N, Delta, true_dev, n_sims, sd_scan=.05, sd_noise=.10,
                      alpha=.05, age_min=25, age_max=65):
    """Simulate data with the true slope differences built in and apply the
    15 TOSTs directly (fitting each difference-score regression by OLS)."""
    tq = stats.t.ppf(1 - alpha, N - 2)
    passed = 0
    for _ in range(n_sims):
        age = rng.uniform(age_min, age_max, N); x = age - age.mean()
        y = (rng.normal(0, sd_scan, (N, 1, 3)) + rng.normal(0, sd_noise, (N, 5, 3))
             + x[:, None, None] * true_dev * Delta)
        ok = True
        for r in range(5):
            for j, k in PAIRS:
                D = y[:, r, j] - y[:, r, k]
                b = (x * D).sum() / (x ** 2).sum()
                res = D - D.mean() - b * x
                se = np.sqrt((res ** 2).sum() / (N - 2) / (x ** 2).sum())
                ok &= (b - tq * se > -Delta) and (b + tq * se < Delta)
        passed += ok
    return passed / n_sims


# 1. Closed-form smallest passing Delta (with true differences) vs brute force.
dev = np.zeros((5, 3)); dev[:, 2] = 0.25                 # default scenario
mixed = np.array([[0, .3, -.2], [.1, -.4, 0], [0, 0, 0],    # mixed signs and sizes
                  [-.2, .2, .5], [.3, 0, -.3]])
for label, d, N, Delta in [("default", dev, 10, 0.022), ("default", dev, 20, 0.015),
                           ("mixed-sign", mixed, 15, 0.025)]:
    closed = (min_passing_delta(N, true_dev=d, n_sims=40000, rng=1) < Delta).mean()
    brute = brute_force_power(N, Delta, d, 3000)
    print(f"1. {label} scenario, N={N}, Delta={Delta}: closed form {closed:.3f}, "
          f"brute force {brute:.3f} (MC SE {np.sqrt(brute * (1 - brute) / 3000):.3f})")

# 2. Per-site mixed-model slopes equal per-region OLS slopes (GLS = OLS).
N = 15
age = rng.uniform(25, 65, N); x = age - age.mean()
y = (rng.normal(0, .3, (N, 1)) + rng.normal(0, .2, (N, 5))
     + np.linspace(-.01, .01, 5) * x[:, None])
d = pd.DataFrame(dict(y=y.ravel(), region=np.tile(np.arange(5), N),
                      age_c=np.repeat(x, 5), subject=np.repeat(np.arange(N), 5)))
fit = smf.mixedlm("y ~ 0 + C(region) + C(region):age_c", d, groups="subject").fit(reml=True)
mixed = np.array([fit.params[f"C(region)[{r}]:age_c"] for r in range(5)])
ols = np.array([np.polyfit(x, y[:, r], 1)[0] for r in range(5)])
print(f"2. max |mixed-model slope - OLS slope| = {np.abs(mixed - ols).max():.1e}")

# 3. Scaling: Delta* proportional to sqrt(sd_scan^2 + sd_noise^2) / age range.
base = np.quantile(min_passing_delta(5, true_dev=dev, rng=2, n_sims=40000), .8)
scaled = np.quantile(min_passing_delta(5, sd_scan=.15, sd_noise=.30, age_min=25,
                                       age_max=45, true_dev=dev, rng=2, n_sims=40000), .8)
print(f"3. noise x3 and age range /2 multiplies Delta* by {scaled / base:.3f} (expect 6)")
