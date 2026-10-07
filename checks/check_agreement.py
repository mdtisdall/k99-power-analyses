"""Numerical checks of docs/agreement-derivation.md.

Run from the repository root:  python3 checks/check_agreement.py
Needs numpy, scipy, pandas, statsmodels (pip install pandas statsmodels).
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.anova import anova_lm

from agreement import (smallest_margins, exact_pooled_power, scenario, ratio, mean_squares,
                       ORTHO, PAIR, N_SITES, N_REGIONS)

rng = np.random.default_rng(0)

# 1. Mean squares match a standard two-way ANOVA (subject, site, subject:site).
N, n = 8, 2
y = (rng.normal(0, 1, (N, N_SITES, n)) + rng.normal(0, 1, (N, 1, 1))
     + np.array([0, .5, -.3])[None, :, None])                  # site offsets included
d = pd.DataFrame(dict(y=y.ravel(), subject=np.repeat(np.arange(N), N_SITES * n),
                      site=np.tile(np.repeat(np.arange(N_SITES), n), N)))
tab = anova_lm(smf.ols("y ~ C(subject) * C(site)", d).fit())
ss_between = tab.loc["C(site)", "sum_sq"] + tab.loc["C(subject):C(site)", "sum_sq"]
ms_site, ms_rep = mean_squares(y[None, ..., None], n)
print(f"1. between-site MS {ms_site[0, 0]:.6f} vs ANOVA {ss_between / (N * (N_SITES - 1)):.6f}; "
      f"repeat MS {ms_rep[0, 0]:.6f} vs ANOVA {tab.loc['Residual', 'mean_sq']:.6f}")

# 2. Pooled sum of squares over orthonormal contrasts = (sum over the 10 pairs) / 5.
Y = rng.normal(0, 1, (1, N, N_SITES, n, N_REGIONS))
a = mean_squares(Y @ ORTHO.T, n)[0].sum(); b = mean_squares(Y @ PAIR.T, n)[0].sum()
print(f"2. pooled orthonormal MS sum {a:.6f} vs pairwise sum / 5 {b / 5:.6f}")

# 3. Simulated pooled power (no offsets) vs the exact F formula.
for N, R, R0 in [(10, 1.1, 1.5), (20, 1.25, 1.5), (8, 1.25, 2.0)]:
    pool = smallest_margins(N, 2.1, np.sqrt(R - 1) * 2.1, None, n_sims=40000, rng=N)[0]
    print(f"3. N={N}, R={R}, R0={R0}: simulated {(pool < R0).mean():.3f}, "
          f"exact {exact_pooled_power(N, R, R0):.3f}")

# 4. Level at the margin (true R = R0), whatever the split between offsets and
#    scatter and wherever the offsets are. Rejection rates must be <= 0.05.
spread = np.zeros((N_SITES, N_REGIONS)); spread[1] = [1, -1, 0, .5, 0]; spread[2] = [0, 1, -1, 0, .5]
one = np.zeros((N_SITES, N_REGIONS)); one[2, 0] = 1
for label, pat in [("one region at one site", one), ("spread over regions", spread)]:
    out = []
    for share in (0, .5, 1):
        sd_site, off = scenario(1.5, share, pat, 2.1)
        assert abs(ratio(2.1, sd_site, off).mean() - 1.5) < 1e-12
        rates = [(smallest_margins(N, 2.1, sd_site, off, n_sims=20000, rng=N)[0] < 1.5).mean()
                 for N in (6, 20, 50)]
        out.append(f"share {share}: " + "/".join(f"{r:.3f}" for r in rates))
    print(f"4. pooled level (N = 6/20/50), offsets {label}: " + "; ".join(out))
# Every-pair test: its level is set by the worst pair, so put R0 at the worst pair's R.
sd_site, off = scenario(1.5, .5, one, 2.1)
worst = ratio(2.1, sd_site, off, PAIR).max()
rates = [(smallest_margins(N, 2.1, sd_site, off, n_sims=20000, rng=N)[1] < worst).mean()
         for N in (6, 20, 50)]
print(f"4. every-pair level at the worst pair's R ({worst:.2f}), N = 6/20/50: "
      + "/".join(f"{r:.3f}" for r in rates))
