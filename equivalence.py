"""Power to show cross-site equivalence of regional age effects.

See docs/equivalence-derivation.md for the derivation. Every subject is scanned
at all 3 sites. Each of the 5 regions x 3 site pairs = 15 comparisons regresses
the per-subject difference between two sites on centered age; the sites are
"the same" only if every 90% CI for that slope lies inside (-Delta, +Delta).

For each N, this script reports the power to show equivalence at margin Delta,
and the smallest Delta for which P(all 15 pass) reaches the target power. Both
assume the planning scenario in true_dev (by default, site 3's age slope
differs from the others by 0.25 * Delta in every region). It saves the table to
equivalence_curve.csv and plots both curves in equivalence_curve.png.

All parameter values below are placeholders. Replace them with estimates from
traveling-subject data (the same subjects scanned at each site).
"""

import numpy as np
from scipy import stats

from plotting import plot_curves

N_REGIONS, N_SITES = 5, 3
PAIRS = [(0, 1), (0, 2), (1, 2)]

# ---- Parameters -------------------------------------------------------------

N_grid = [5, 10, 15, 20, 25, 30, 40, 50]  # numbers of subjects (each scanned at all 3 sites)
target_power = 0.80                # required P(all 15 equivalence tests pass)
Delta = 0.015                      # margin (outcome units/year) for the power curve
sd_scan = 0.05                     # SD of a whole-scan offset (shared by all regions in one scan)
sd_noise = 0.10                    # SD of region-level measurement noise (all sites and regions)
age_min, age_max = 25, 65          # age distribution: uniform(age_min, age_max)
true_dev = np.zeros((N_REGIONS, N_SITES))  # planning scenario: each site's true age-slope
true_dev[:, 2] = 0.25              # deviation, as a fraction of Delta (regions x sites)
alpha = 0.05                       # one-sided level for each test (gives a 90% CI)
n_sims = 10000
seed = 1

# ---- Simulation -------------------------------------------------------------


def _check(N, true_dev):
    if N < 3:
        raise ValueError(f"N must be at least 3 (got {N}); each regression has N - 2 df.")
    if np.shape(true_dev) != (N_REGIONS, N_SITES):
        raise ValueError(f"true_dev must have shape {(N_REGIONS, N_SITES)} "
                         f"(regions x sites), got {np.shape(true_dev)}.")
    g = np.array([true_dev[:, j] - true_dev[:, k] for j, k in PAIRS]).T
    if np.abs(g).max() >= 1:
        raise ValueError("Every true site difference must be smaller than Delta "
                         "(|true_dev[:, j] - true_dev[:, k]| < 1), or equivalence "
                         "can never be shown.")
    return g                                                   # (regions, pairs)


def min_passing_delta(N, sd_scan=.05, sd_noise=.10, true_dev=None, alpha=.05,
                      age_min=25, age_max=65, n_sims=10000, rng=None, chunk=2000):
    """For each of n_sims simulated studies, the smallest Delta at which all 15
    tests would pass, given true site differences of true_dev * Delta."""
    if true_dev is None:
        true_dev = np.zeros((N_REGIONS, N_SITES))
    g = _check(N, true_dev)
    rng = np.random.default_rng(rng)
    tq = stats.t.ppf(1 - alpha, N - 2)
    out = []
    for start in range(0, n_sims, chunk):              # chunks bound memory use
        n = min(chunk, n_sims - start)

        # Simulate only the terms that do not cancel in between-site differences.
        age = rng.uniform(age_min, age_max, (n, N))
        age_c = age - age.mean(axis=1, keepdims=True)
        e = (rng.normal(0, sd_scan, (n, N, 1, N_SITES))           # whole-scan offset
             + rng.normal(0, sd_noise, (n, N, N_REGIONS, N_SITES)))  # measurement noise
        D = np.stack([e[..., j] - e[..., k] for j, k in PAIRS], axis=-1)

        # OLS of D on age_c (centered, so the slope needs no intercept adjustment).
        x = age_c[:, :, None, None]
        sxx = (x ** 2).sum(axis=1)
        b = (x * D).sum(axis=1) / sxx                   # noise part of each slope
        resid = D - D.mean(axis=1, keepdims=True) - b[:, None] * x
        se = np.sqrt((resid ** 2).sum(axis=1) / (N - 2) / sxx)

        # A true difference g * Delta shifts the slope estimate by exactly that
        # and leaves the SE unchanged. The 90% CI is inside (-Delta, +Delta) iff
        # (1 - g) Delta > b + tq se  and  (1 + g) Delta > tq se - b.
        need = np.maximum((b + tq * se) / (1 - g), (tq * se - b) / (1 + g))
        out.append(need.max(axis=(1, 2)))
    return np.concatenate(out)


def detectable_delta(N, target_power=.80, **kw):
    """Smallest Delta with P(all 15 tests pass) >= target_power."""
    return np.quantile(min_passing_delta(N, **kw), target_power)


def power(N, Delta, **kw):
    """P(all 15 tests pass) at margin Delta."""
    return (min_passing_delta(N, **kw) < Delta).mean()


# ---- Power curve over N -----------------------------------------------------

if __name__ == "__main__":
    kw = dict(sd_scan=sd_scan, sd_noise=sd_noise, true_dev=true_dev, alpha=alpha,
              age_min=age_min, age_max=age_max, n_sims=n_sims)
    rows = []
    for N in N_grid:
        # Seeded per N, so a row does not change when N_grid is edited.
        dmin = min_passing_delta(N, rng=[seed, N], **kw)
        rows.append((N, (dmin < Delta).mean(), np.quantile(dmin, target_power)))

    ci = f"{100 * (1 - 2 * alpha):g}%"
    worst = np.abs(_check(3, true_dev)).max()
    scenario = (f"true site differences up to {worst:g}\u0394" if worst else
                "no true site differences")
    print(f"Planning scenario: {scenario}".replace("\u0394", " x Delta"))
    print(f"Power to show equivalence at Delta = {Delta:g}, and smallest Delta "
          f"(outcome units/year) with {target_power:.0%} power:")
    print("   N  Power at Delta  Smallest Delta")
    for N, p, d in rows:
        print(f"{N:4d}  {p:14.3f}  {d:14.4g}")
    np.savetxt("equivalence_curve.csv", rows, delimiter=",",
               header="N,power_at_Delta,smallest_Delta", comments="",
               fmt=["%d", "%.4f", "%.6g"])

    plot_curves(
        [r[0] for r in rows], "Subjects (each scanned at all 3 sites)",
        [dict(y=[r[1] for r in rows],
              title=f"Power to show equivalence at Δ = {Delta:g}",
              ylim=(0, 1), ref=(target_power, f"{target_power:.0%} target")),
         dict(y=[r[2] for r in rows],
              title=f"Smallest Δ at {target_power:.0%} power (units/year)")],
        f"3 sites, 5 regions; all 15 {ci} CIs must fall inside ±Δ; {scenario}",
        "equivalence_curve.png")
    print("\nSaved equivalence_curve.csv and equivalence_curve.png")
