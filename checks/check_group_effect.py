"""Numerical checks of docs/group-effect-derivation.md.

Run from the repository root:  python3 checks/check_group_effect.py [n_sims]
Needs numpy, scipy, pandas, statsmodels (pip install pandas statsmodels).
The mixed-model simulation takes about 0.4 s per study; the derivation quotes
1500 studies (about 10 minutes). The default is 200.
"""
import os, sys, warnings
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

from group_effect import make_design, group_se, dummies

warnings.filterwarnings("ignore")
rng = np.random.default_rng(0)

# 1. Closed-form SE vs brute-force GLS (X' V^-1 X)^-1 with known variances.
for S, n, K, R, rho in [(20, 10, 2, 10, 0.5), (20, 3, 2, 10, 0.9), (4, 7, 3, 6, 0.8),
                        (3, 5, 2, 1, 0.3)]:
    site, group, x = make_design(S, n, K, 25, 65, rng)
    V = rho * np.ones((R, R)) + (1 - rho) * np.eye(R)
    Vi = np.linalg.inv(V)
    XtVX = 0
    for i in range(site.size):
        g = dummies(group[i:i + 1], K)[0]
        s = dummies(site[i:i + 1], S)[0]
        Xi = np.hstack([np.eye(R), np.eye(R) * x[i], np.kron(g, np.eye(R)),
                        np.outer(np.ones(R), s)])
        XtVX = XtVX + Xi.T @ Vi @ Xi
    gls = np.sqrt(np.diag(np.linalg.inv(XtVX))[2 * R:2 * R + R * (K - 1)]).reshape(K - 1, R)
    se, df = group_se(site, group, x, S, K, R, 1.0, rho)
    print(f"1. S={S} n={n} K={K} R={R} region_corr={rho}: "
          f"max |GLS - closed form| = {np.abs(gls - se[:, None]).max():.1e}")

# 2. Power and false positives from real mixed-model fits at the default design.
S, n, K, R, rho, effect, a = 20, 10, 2, 10, 0.5, 0.5193, 0.005
n_sims = int(sys.argv[1]) if len(sys.argv) > 1 else 200
hit = false_pos = 0
for _ in range(n_sims):
    site, group, x = make_design(S, n, K, 25, 65, rng)
    N = site.size
    y = (np.sqrt(rho) * rng.normal(size=(N, 1)) + np.sqrt(1 - rho) * rng.normal(size=(N, R))
         + rng.normal(size=R)[None]                       # region means
         + np.linspace(-.02, .02, R)[None] * x[:, None]   # region age slopes
         + 0.4 * rng.normal(size=S)[site][:, None]        # site offsets
         + effect * (np.arange(R) == 0)[None] * group[:, None])  # effect in region 0
    d = pd.DataFrame(dict(y=y.ravel(), region=np.tile(np.arange(R), N),
                          age_c=np.repeat(x, R), site=np.repeat(site, R),
                          group=np.repeat(group, R), subject=np.repeat(np.arange(N), R)))
    p = smf.mixedlm("y ~ 0 + C(region) + C(region):age_c + C(site) + C(region):group",
                    d, groups="subject").fit(reml=True).pvalues
    hit += p["C(region)[0]:group"] < a
    false_pos += p["C(region)[1]:group"] < a
print(f"2. {n_sims} mixed-model fits: power {hit / n_sims:.3f} (expect 0.80, MC SE "
      f"{np.sqrt(.16 / n_sims):.3f}); false positives {false_pos / n_sims:.4f} (expect {a})")
