"""Numerical checks of docs/group-effect-derivation.md.
See also check_group_effect_lmer.R, which checks repeated visits with lmerTest.

Run from the repository root:
    python3 checks/check_group_effect.py [n_sims] [subjects_per_site]
Needs numpy, scipy, pandas, statsmodels (pip install pandas statsmodels).
The mixed-model simulation takes about 0.4 s per study; the derivation quotes
1500 studies (about 10 minutes). The default is 200.
"""
import os, sys, warnings
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

from scipy import stats

from group_effect import (make_design, group_se, dummies, power_curve_inputs,
                          detectable_effect, variance_components)

warnings.filterwarnings("ignore")
rng = np.random.default_rng(0)

# 1. Closed-form SE vs brute-force GLS (X' V^-1 X)^-1 with known variances,
#    for one visit and for repeated visits (visit_years, retest_corr,
#    visit_region_corr).
for S, n, K, R, rho, years, rt, vrc in [
        (20, 10, 2, 10, 0.5, [0], 1.0, 0.0), (20, 3, 2, 10, 0.9, [0], 1.0, 0.0),
        (4, 7, 3, 6, 0.8, [0], 1.0, 0.0), (3, 5, 2, 1, 0.3, [0], 1.0, 0.0),
        (14, 6, 2, 10, 0.5, [0, 1, 2], 0.8, 0.5), (20, 3, 2, 10, 0.9, [0, 1, 2], 0.95, 0.2),
        (4, 7, 3, 6, 0.8, [0, 1, 3, 7], 0.6, 0.9), (3, 5, 2, 1, 0.3, [0, 2], 0.5, 0.5)]:
    site, group, x = make_design(S, n, K, 25, 65, rng)
    T = len(years)
    vc = variance_components(1.0, rho, rt, vrc)
    JT, JR, IT, IR = np.ones((T, T)), np.ones((R, R)), np.eye(T), np.eye(R)
    V = (vc["subject"] * np.kron(JT, JR) + vc["subject_region"] * np.kron(JT, IR)
         + vc["visit"] * np.kron(IT, JR) + vc["visit_region"] * np.kron(IT, IR))
    Vi = np.linalg.inv(V)
    dt = np.array(years) - np.mean(years)
    XtVX = 0
    for i in range(site.size):
        g = dummies(group[i:i + 1], K)[0]
        s = dummies(site[i:i + 1], S)[0]
        Xi = np.vstack([np.hstack([np.eye(R), np.eye(R) * (x[i] + dt[t]),
                                   np.kron(g, np.eye(R)), np.outer(np.ones(R), s)])
                        for t in range(T)])
        XtVX = XtVX + Xi.T @ Vi @ Xi
    gls = np.sqrt(np.diag(np.linalg.inv(XtVX))[2 * R:2 * R + R * (K - 1)]).reshape(K - 1, R)
    se, df = group_se(site, group, x, S, K, R, 1.0, rho, years, rt, vrc)
    print(f"1. S={S} n={n} K={K} R={R} region_corr={rho} visit_years={years}: "
          f"max |GLS - closed form| = {np.abs(gls - se[:, None]).max():.1e}")

# 2. Power and false positives from real mixed-model fits, one scan per subject
#    (check_group_effect_lmer.R checks repeated visits with lmerTest). Each fit's Wald
#    statistic (estimate / SE) is compared with the t critical value at the
#    closed-form Satterthwaite df, as lmerTest would; MixedLM's own p-values are
#    z-based and would be too liberal at small df.
n_sims = int(sys.argv[1]) if len(sys.argv) > 1 else 200
n = int(sys.argv[2]) if len(sys.argv) > 2 else 10           # subjects per site
S, K, R, rho, a = 20, 2, 10, 0.5, 0.005
se_d, df_d = power_curve_inputs(S, n, K, R, 1.0, rho, 25, 65, 500, np.random.default_rng(1))
effect = detectable_effect(se_d[:, 0], df_d[:, 0], a, 0.8)
hit = false_pos = 0
for _ in range(n_sims):
    site, group, x = make_design(S, n, K, 25, 65, rng)
    N = site.size
    tc = stats.t.ppf(1 - a / 2, group_se(site, group, x, S, K, R, 1.0, rho)[1][0])
    y = (np.sqrt(rho) * rng.normal(size=(N, 1)) + np.sqrt(1 - rho) * rng.normal(size=(N, R))
         + rng.normal(size=R)[None]                       # region means
         + np.linspace(-.02, .02, R)[None] * x[:, None]   # region age slopes
         + 0.4 * rng.normal(size=S)[site][:, None]        # site offsets
         + effect * (np.arange(R) == 0)[None] * group[:, None])  # effect in region 0
    d = pd.DataFrame(dict(y=y.ravel(), region=np.tile(np.arange(R), N),
                          age_c=np.repeat(x, R), site=np.repeat(site, R),
                          group=np.repeat(group, R), subject=np.repeat(np.arange(N), R)))
    fit = smf.mixedlm("y ~ 0 + C(region) + C(region):age_c + C(site) + C(region):group",
                      d, groups="subject").fit(reml=True)
    t = fit.params / fit.bse
    hit += abs(t["C(region)[0]:group"]) > tc
    false_pos += abs(t["C(region)[1]:group"]) > tc
print(f"2. {S} sites x {n} subjects, effect d = {effect:.3f} in one region "
      f"(closed-form power 0.80), {n_sims} mixed-model fits: power {hit / n_sims:.3f} "
      f"(MC SE {np.sqrt(.16 / n_sims):.3f}); false positives {false_pos / n_sims:.4f} "
      f"(expect {a})")
