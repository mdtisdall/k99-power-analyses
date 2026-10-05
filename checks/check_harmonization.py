"""Numerical checks of docs/harmonization-derivation.md.

Run from the repository root:  python3 checks/check_harmonization.py
Needs numpy, scipy, pandas, statsmodels (pip install pandas statsmodels).
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

from harmonization import z_stats, disagreement, PAIRS

rng = np.random.default_rng(0)
R, S = 5, 3
ONE = np.ones((R, S))


def spread(p):
    g = ONE.copy(); g[:, 1], g[:, 2] = 1 - p, 1 + p
    return g


def brute_force(N, gain_A, gain_B, n_sims, slope=1.0, sd_between=20., sd_scan=1.,
                sd_noise=2.1, sd_site=2.9, n_rep=2, noise_corr=.9, alpha=.05):
    """Simulate every scan explicitly (n_rep per site), average them, fit each
    site's slope separately, form the slope differences, and apply the pooled
    test. Returns the rejection rate, and the mean of the pooled estimate with its
    Monte Carlo SE."""
    Kpair = np.zeros((3, 2)); Kpair[0, 0] = Kpair[1, 1] = 1; Kpair[2] = [-1, 1]
    rej, ests = 0, []
    for _ in range(n_sims):
        x = rng.uniform(25, 65, N); x = x - x.mean()
        chi = slope * x[:, None] + rng.normal(0, sd_between, (N, R))
        site = rng.normal(0, sd_site, (2, N, R, S))                 # persists across scans
        scans = []
        for _ in range(n_rep):
            w = rng.normal(0, sd_scan, (2, N, 1, S))
            e = rng.normal(0, sd_noise, (2, N, R, S))
            scans.append(w + e)
        scan_noise = np.mean(scans, axis=0)
        noise_B = site[0] + scan_noise[0]
        noise_A = noise_corr * noise_B + np.sqrt(1 - noise_corr ** 2) * (site[1] + scan_noise[1])
        y = {"A": gain_A * chi[..., None] + noise_A, "B": gain_B * chi[..., None] + noise_B}
        z, resid = [], []
        for r in range(R):
            for m in "AB":
                b = np.array([np.polyfit(x, y[m][:, r, s], 1)[0] for s in range(S)])
                for j, k in [(0, 1), (0, 2)]:
                    D = y[m][:, r, j] - y[m][:, r, k]
                    g = (x * D).sum() / (x ** 2).sum()
                    assert abs(g - (b[j] - b[k])) < 1e-9          # Step 2 identity
                    z.append(g); resid.append(D - D.mean() - g * x)
        z, resid = np.array(z), np.array(resid)
        cov = resid @ resid.T / (N - 2) / (x ** 2).sum()
        K = np.array([[2., -1.], [-1., 2.]])
        Q = np.kron(np.eye(R), np.block([[-K, np.zeros((2, 2))], [np.zeros((2, 2)), K]]))
        est = z @ Q @ z - np.trace(Q @ cov)
        var = 2 * np.trace(Q @ cov @ Q @ cov) + 4 * z @ Q @ cov @ Q @ z
        ests.append(est)
        rej += est / np.sqrt(var) > stats.norm.isf(alpha)
    return rej / n_sims, np.mean(ests), np.std(ests) / np.sqrt(n_sims)


crit = stats.norm.isf(.05)

# 1. Vectorized script vs explicit brute force (repeats simulated scan by scan).
for N, gA, gB in [(20, ONE, spread(.2)), (40, spread(.05), spread(.1))]:
    vec = (z_stats(N, gain_A=gA, gain_B=gB, n_sims=20000, rng=1)[1] > crit).mean()
    brute, mean_est, se_est = brute_force(N, gA, gB, 2000)
    truth = (disagreement(gB) - disagreement(gA)).sum()
    print(f"1. N={N}: power script {vec:.3f}, brute force {brute:.3f} "
          f"(MC SE {np.sqrt(brute * (1 - brute) / 2000):.3f}); mean estimate "
          f"{mean_est:.4f} (MC SE {se_est:.4f}) vs true difference {truth:.4f}")

# 2. Size at the null boundary (equal disagreement): rejection rate must be <= 0.05.
mirror = ONE.copy(); mirror[:, 1], mirror[:, 2] = 1.2, 0.8
for label, gA, gB, kw in [("both perfect", ONE, ONE, {}),
                          ("same +-20% gains", spread(.2), spread(.2), {}),
                          ("mirrored gains, noise_corr 0.5", mirror, spread(.2),
                           dict(noise_corr=.5))]:
    rates = []
    for N in (10, 20, 50, 100):
        z_reg, z_all = z_stats(N, gain_A=gA, gain_B=gB, n_sims=20000, rng=N, **kw)
        rates.append(f"N={N}: pooled {(z_all > crit).mean():.3f}, "
                     f"region {(z_reg > crit).mean():.3f}")
    print(f"2. size, {label}: " + "; ".join(rates))

# 3. Per-site mixed-model slopes equal per-region OLS slopes (GLS = OLS).
N = 15
age = rng.uniform(25, 65, N); x = age - age.mean()
y = (rng.normal(0, .3, (N, 1)) + rng.normal(0, .2, (N, 5))
     + np.linspace(-.01, .01, 5) * x[:, None])
d = pd.DataFrame(dict(y=y.ravel(), region=np.tile(np.arange(5), N),
                      age_c=np.repeat(x, 5), subject=np.repeat(np.arange(N), 5)))
fit = smf.mixedlm("y ~ 0 + C(region) + C(region):age_c", d, groups="subject").fit(reml=True)
mixed = np.array([fit.params[f"C(region)[{r}]:age_c"] for r in range(5)])
ols = np.array([np.polyfit(x, y[:, r], 1)[0] for r in range(5)])
print(f"3. max |mixed-model slope - OLS slope| = {np.abs(mixed - ols).max():.1e}")
