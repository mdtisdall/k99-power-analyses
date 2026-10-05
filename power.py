"""Smallest equivalence margin (Delta) detectable for cross-site age effects.

See README.md for the design, method, and inputs. Each of the
5 regions x 3 site pairs = 15 comparisons regresses the per-subject difference
score D on centered age; sites are "the same" only if every 90% CI lies inside
(-Delta, +Delta). For a given N, this script reports the smallest Delta for
which P(all 15 pass) reaches the target power.

All parameter values below are placeholders. Replace them with estimates from
prior (ideally traveling-subject / test-retest) data.
"""

import numpy as np
from scipy import stats

# ---- Parameters -------------------------------------------------------------

N_grid = [5]                       # number of subjects (each scanned at all 3 sites)
target_power = 0.80                # required P(all 15 TOSTs pass)
s_ss = 0.05                        # SD of subject:site offset (shared across regions in a scan)
sigma = 0.10                       # region-level measurement-noise SD (all sites and regions)
age_min, age_max = 25, 65          # age distribution: uniform(age_min, age_max)
dslope = np.zeros((5, 3))          # true slope deviation, regions (rows) x sites (cols)
alpha = 0.05                       # one-sided level for each TOST (gives a 90% CI)
n_sims = 10000
seed = 1

PAIRS = [(0, 1), (0, 2), (1, 2)]

# ---- Simulation -------------------------------------------------------------


def min_passing_delta(N, s_ss=.05, sigma=.10, dslope=np.zeros((5, 3)),
                      alpha=.05, age_min=25, age_max=65, n_sims=10000, rng=None):
    """For each of n_sims simulated studies, the smallest Delta at which all 15
    TOSTs would pass: max over comparisons of |slope diff| + t * SE."""
    rng = np.random.default_rng(rng)
    dslope = np.asarray(dslope)

    # Simulate only the terms that do not cancel in between-site differences.
    age = rng.uniform(age_min, age_max, (n_sims, N))
    age_c = age - age.mean(axis=1, keepdims=True)
    e = (rng.normal(0, s_ss, (n_sims, N, 1, 3))                 # subject:site offset
         + rng.normal(0, sigma, (n_sims, N, 5, 3))              # measurement noise
         + age_c[:, :, None, None] * dslope)                    # true site slope deviation

    D = np.stack([e[..., j] - e[..., k] for j, k in PAIRS], axis=-1)  # (sims, N, 5, 3)

    # OLS of D on age_c (centered, so the slope needs no intercept adjustment).
    x = age_c[:, :, None, None]
    sxx = (x ** 2).sum(axis=1)
    b = (x * D).sum(axis=1) / sxx
    resid = D - D.mean(axis=1, keepdims=True) - b[:, None] * x
    se = np.sqrt((resid ** 2).sum(axis=1) / (N - 2) / sxx)

    # A comparison passes iff its 90% CI is inside (-Delta, +Delta),
    # i.e. iff |b| + tq * se < Delta.
    tq = stats.t.ppf(1 - alpha, N - 2)
    return (np.abs(b) + tq * se).max(axis=(1, 2))


def detectable_delta(N, target_power=.80, **kw):
    """Smallest Delta with P(all 15 TOSTs pass) >= target_power."""
    return np.quantile(min_passing_delta(N, **kw), target_power)


def power(N, Delta, **kw):
    """P(all 15 TOSTs pass) at a given Delta."""
    return (min_passing_delta(N, **kw) < Delta).mean()


# ---- Detectable Delta by N --------------------------------------------------

if __name__ == "__main__":
    rng = np.random.default_rng(seed)
    kw = dict(s_ss=s_ss, sigma=sigma, dslope=dslope, alpha=alpha,
              age_min=age_min, age_max=age_max, n_sims=n_sims)
    rows = [(N, detectable_delta(N, target_power, rng=rng, **kw)) for N in N_grid]
    print(f"Smallest Delta (outcome units/year) with {target_power:.0%} power:")
    print("   N  Delta")
    for N, d in rows:
        print(f"{N:4d}  {d:.4g}")
    np.savetxt("detectable_delta.csv", rows, delimiter=",", header="N,Delta",
               comments="", fmt=["%d", "%.6g"])
