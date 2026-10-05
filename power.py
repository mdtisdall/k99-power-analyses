"""Power for cross-site equivalence of region-specific age effects.

See README.md for the design, method, and inputs. Each of the
5 regions x 3 site pairs = 15 comparisons regresses the per-subject difference
score D on centered age; sites are "the same" only if every 90% CI lies inside
(-Delta, +Delta). Power = P(all 15 pass).

All parameter values below are placeholders. Replace them with estimates from
prior (ideally traveling-subject / test-retest) data.
"""

import numpy as np
from scipy import stats

# ---- Parameters -------------------------------------------------------------

Delta = 0.003                     # equivalence margin, outcome units per year
s_ss = 0.05                       # SD of subject:site offset (shared across regions in a scan)
sigma = np.array([.10, .10, .10])  # region-level measurement-noise SD at each site
age_min, age_max = 50, 80         # age distribution: uniform(age_min, age_max)
dslope = np.zeros((5, 3))         # true slope deviation, regions (rows) x sites (cols)
alpha = 0.05                      # one-sided level for each TOST (gives a 90% CI)
N_grid = range(20, 201, 20)
n_sims = 2000
seed = 1

PAIRS = [(0, 1), (0, 2), (1, 2)]

# ---- Simulation -------------------------------------------------------------


def power(N, Delta, s_ss=.05, sigma=(.10, .10, .10), dslope=np.zeros((5, 3)),
          alpha=.05, age_min=50, age_max=80, n_sims=2000, rng=None):
    """P(all 15 TOSTs pass), estimated from n_sims simulated studies at once."""
    rng = np.random.default_rng(rng)
    sigma, dslope = np.asarray(sigma), np.asarray(dslope)

    # Simulate only the terms that do not cancel in between-site differences.
    age = rng.uniform(age_min, age_max, (n_sims, N))
    age_c = age - age.mean(axis=1, keepdims=True)
    e = (rng.normal(0, s_ss, (n_sims, N, 1, 3))                 # subject:site offset
         + rng.normal(0, 1, (n_sims, N, 5, 3)) * sigma          # measurement noise
         + age_c[:, :, None, None] * dslope)                    # true site slope deviation

    D = np.stack([e[..., j] - e[..., k] for j, k in PAIRS], axis=-1)  # (sims, N, 5, 3)

    # OLS of D on age_c (centered, so the slope needs no intercept adjustment).
    x = age_c[:, :, None, None]
    sxx = (x ** 2).sum(axis=1)
    b = (x * D).sum(axis=1) / sxx
    resid = D - D.mean(axis=1, keepdims=True) - b[:, None] * x
    se = np.sqrt((resid ** 2).sum(axis=1) / (N - 2) / sxx)

    tq = stats.t.ppf(1 - alpha, N - 2)
    ok = (b - tq * se > -Delta) & (b + tq * se < Delta)          # (sims, 5, 3)
    return ok.all(axis=(1, 2)).mean()


# ---- Power curve ------------------------------------------------------------

if __name__ == "__main__":
    rng = np.random.default_rng(seed)
    rows = [(N, power(N, Delta, s_ss, sigma, dslope, alpha, age_min, age_max,
                      n_sims, rng)) for N in N_grid]
    print("   N  power")
    for N, p in rows:
        print(f"{N:4d}  {p:.4f}")
    np.savetxt("power_curve.csv", rows, delimiter=",", header="N,power",
               comments="", fmt=["%d", "%.4f"])
