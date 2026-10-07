"""Power to show that scanners preserve each subject's between-region differences
about as well as repeat scans on one scanner do.

See docs/agreement-derivation.md for the derivation. Every subject is scanned
n_rep times at each of 3 sites. The outcome is every pairwise difference between
the 5 regions within one scan (10 differences). For each, the target ratio is

    R = E[(difference between two sites, single scans)^2]
        / E[(difference between two repeats at one site)^2],

which counts both fixed scanner offsets and subject-specific scanner scatter as
disagreement. R is at least 1. The claim is R < R0 (a margin), tested with an F
ratio of the between-site mean square within subjects (site offsets not
removed) to the repeat mean square.

- Primary test: pooled over all region pairs (equivalently, over 4 orthonormal
  contrasts), one-sided, at margin R0_pooled.
- Secondary test: every one of the 10 pairs must pass at margin R0_pair
  (intersection-union, no multiplicity correction).

For each planning scenario and N, the script reports the power of both tests and
the smallest margin each test can show with the target power. It saves the table
to agreement_curve.csv and plots power in agreement_curve.png.
"""

from itertools import combinations

import numpy as np
from scipy import stats

from plotting import plot_curves

N_REGIONS, N_SITES = 5, 3
REGION_PAIRS = list(combinations(range(N_REGIONS), 2))   # the 10 pairwise differences
SITE_PAIRS = list(combinations(range(N_SITES), 2))

# ---- Parameters -------------------------------------------------------------

N_grid = [6, 8, 10, 15, 20, 30, 40, 50]  # numbers of subjects (each scanned at all 3 sites)
n_rep = 2                          # scans per subject at each site
target_power = 0.80
alpha = 0.05                       # one-sided level of each test
R0_pooled = 1.5                    # margin for the pooled test
R0_pair = 2.0                      # margin for every-pair test
sd_noise = 2.1                     # SD (ppb) of region-level scan-to-scan noise (Naji et al. 2022);
                                   # only used to express scenarios in ppb
scenarios = {"R = 1.1": 1.1,       # planning scenarios: true ratio R
             "R = 1.25": 1.25}     # (at most two, for the plot)
offset_share = 0.5                 # fraction of R - 1 due to fixed scanner offsets
offset_pattern = np.zeros((N_SITES, N_REGIONS))
offset_pattern[2, 0] = 1.0         # where the offsets are (scaled to fit offset_share):
                                   # site 3 shifts region 1 relative to the others
n_sims = 20000
seed = 1

# ---- Model ------------------------------------------------------------------

PAIR = np.zeros((len(REGION_PAIRS), N_REGIONS))
for c, (a, b) in enumerate(REGION_PAIRS):
    PAIR[c, a], PAIR[c, b] = 1, -1
# 4 orthonormal contrasts spanning the same space as the 10 pairwise differences.
ORTHO = np.linalg.qr(np.column_stack([np.ones(N_REGIONS), np.eye(N_REGIONS)[:, :4]]))[0][:, 1:].T


def offset_excess(offsets, M):
    """Mean over site pairs of the squared offset difference, for each contrast
    in the rows of M (ppb^2)."""
    oc = np.asarray(offsets, float) @ M.T
    return np.mean([(oc[j] - oc[k]) ** 2 for j, k in SITE_PAIRS], axis=0)


def ratio(sd_noise, sd_site, offsets, M=ORTHO):
    """True R for each contrast m in M:
    1 + (sd_site^2 + offset excess / (2 |m|^2)) / sd_noise^2.
    The pooled R is the mean over the orthonormal contrasts."""
    norm2 = (M ** 2).sum(axis=1)
    return 1 + (sd_site ** 2 + offset_excess(offsets, M) / (2 * norm2)) / sd_noise ** 2


def scenario(R, offset_share, offset_pattern, sd_noise):
    """sd_site and offsets (ppb) giving pooled ratio R, with offset_share of
    R - 1 from fixed offsets in the shape of offset_pattern."""
    unit = offset_excess(offset_pattern, ORTHO).mean()      # per orthonormal contrast
    if offset_share and not unit:
        raise ValueError("offset_pattern must shift some region relative to the others.")
    scale = np.sqrt(2 * offset_share * (R - 1) * sd_noise ** 2 / unit) if offset_share else 0.0
    sd_site = np.sqrt((1 - offset_share) * (R - 1)) * sd_noise
    return sd_site, scale * np.asarray(offset_pattern, float)


# ---- Simulation -------------------------------------------------------------


def _check(N, n_rep):
    if int(N) != N or N < 2:
        raise ValueError(f"N must be a whole number of at least 2 (got {N}).")
    if int(n_rep) != n_rep or n_rep < 2:
        raise ValueError(f"n_rep must be at least 2 (got {n_rep}): repeats estimate "
                         f"the within-scanner noise.")


def mean_squares(Y, n_rep):
    """Y: (sims, N, sites, n_rep, k). Between-site mean square within subjects
    (site offsets not removed; df N (S - 1)) and repeat mean square
    (df N S (n_rep - 1)), for each of the k contrasts."""
    N, S = Y.shape[1], Y.shape[2]
    cell = Y.mean(axis=3)
    ms_site = n_rep * ((cell - cell.mean(axis=2, keepdims=True)) ** 2).sum(axis=(1, 2))
    ms_rep = ((Y - cell[:, :, :, None]) ** 2).sum(axis=(1, 2, 3))
    return ms_site / (N * (S - 1)), ms_rep / (N * S * (n_rep - 1))


def smallest_margins(N, sd_noise=2.1, sd_site=0.0, offsets=None, n_rep=2, alpha=.05,
                     n_sims=20000, rng=None, chunk=5000):
    """For each simulated study, the smallest margin R0 the pooled test and the
    every-pair test would pass. Returns two arrays of length n_sims."""
    _check(N, n_rep)
    N, n_rep = int(N), int(n_rep)
    offsets = np.zeros((N_SITES, N_REGIONS)) if offsets is None else np.asarray(offsets, float)
    rng = np.random.default_rng(rng)
    d1, d2 = N * (N_SITES - 1), N * N_SITES * (n_rep - 1)
    q_pair = stats.f.ppf(alpha, d1, d2)
    q_pool = stats.f.ppf(alpha, ORTHO.shape[0] * d1, ORTHO.shape[0] * d2)
    pooled, pairs = [], []
    for start in range(0, n_sims, chunk):
        n = min(chunk, n_sims - start)
        # Subject patterns and whole-scan offsets cancel in within-scan differences
        # and in the mean squares, so only these terms are simulated.
        y = (rng.normal(0, sd_site, (n, N, N_SITES, 1, N_REGIONS))
             + rng.normal(0, sd_noise, (n, N, N_SITES, n_rep, N_REGIONS))
             + offsets[None, None, :, None, :])
        # Test passes at margin R0 iff F < (1 + n_rep (R0 - 1)) q, so the smallest
        # margin passed is 1 + (F / q - 1) / n_rep.
        s, r = mean_squares(y @ ORTHO.T, n_rep)
        pooled.append(1 + (s.sum(1) / r.sum(1) / q_pool - 1) / n_rep)
        s, r = mean_squares(y @ PAIR.T, n_rep)
        pairs.append(1 + ((s / r).max(axis=1) / q_pair - 1) / n_rep)
    return np.concatenate(pooled), np.concatenate(pairs)


def exact_pooled_power(N, R, R0, n_rep=2, alpha=.05):
    """Pooled power with no fixed offsets (exact): F = (1 + n_rep (R - 1)) F(d1, d2)."""
    d1, d2 = 4 * N * (N_SITES - 1), 4 * N * N_SITES * (n_rep - 1)
    tau0, tau = 1 + n_rep * (R0 - 1), 1 + n_rep * (R - 1)
    return stats.f.cdf(tau0 / tau * stats.f.ppf(alpha, d1, d2), d1, d2)


# ---- Power curves over N ----------------------------------------------------

if __name__ == "__main__":
    rows, curves = [], {}
    for name, R in scenarios.items():
        sd_site, offsets = scenario(R, offset_share, offset_pattern, sd_noise)
        rms = np.sqrt(offset_excess(offsets, PAIR).mean())
        print(f"{name}: sd_site = {sd_site:.2f} ppb; fixed offsets give an RMS "
              f"between-site difference of {rms:.2f} ppb in the regional differences "
              f"(largest {np.sqrt(offset_excess(offsets, PAIR).max()):.2f} ppb); "
              f"per-pair R from {ratio(sd_noise, sd_site, offsets, PAIR).min():.2f} "
              f"to {ratio(sd_noise, sd_site, offsets, PAIR).max():.2f}")
        curves[name] = ([], [])
        for N in N_grid:
            # Seeded per scenario and N, so a row does not change when N_grid is edited.
            pool, pair = smallest_margins(N, sd_noise, sd_site, offsets, n_rep, alpha,
                                          n_sims, rng=[seed, N, int(1000 * R)])
            row = (R, N, (pool < R0_pooled).mean(), np.quantile(pool, target_power),
                   (pair < R0_pair).mean(), np.quantile(pair, target_power))
            rows.append(row)
            curves[name][0].append(row[2]); curves[name][1].append(row[4])

    print(f"\nPower at the margins (pooled R0 = {R0_pooled:g}, every pair R0 = {R0_pair:g}), "
          f"and the smallest margin shown with {target_power:.0%} power:")
    print("     R    N  Pooled power  Pooled margin  Every-pair power  Every-pair margin")
    for R, N, pp, pm, ap, am in rows:
        print(f"{R:6.2f} {N:4d}  {pp:12.3f}  {pm:13.3f}  {ap:16.3f}  {am:17.3f}")
    np.savetxt("agreement_curve.csv", rows, delimiter=",",
               header="R_true,N,power_pooled,smallest_margin_pooled,power_every_pair,"
                      "smallest_margin_every_pair",
               comments="", fmt=["%.4g", "%d", "%.4f", "%.4f", "%.4f", "%.4f"])

    scans = f"{n_rep} scans at each of 3 sites"
    plot_curves(
        N_grid, f"Subjects ({scans})",
        [dict(series=[(f"True {k}", v[0]) for k, v in curves.items()],
              title=f"Pooled over region pairs, margin R0 = {R0_pooled:g}",
              ylim=(0, 1), ref=(target_power, f"{target_power:.0%} target")),
         dict(series=[(f"True {k}", v[1]) for k, v in curves.items()],
              title=f"Every region pair, margin R0 = {R0_pair:g}",
              ylim=(0, 1), ref=(target_power, f"{target_power:.0%} target"))],
        f"Power to show between-scanner / repeat disagreement in regional differences "
        f"< R0 ({offset_share:.0%} of the excess from fixed offsets)",
        "agreement_curve.png")
    print("\nSaved agreement_curve.csv and agreement_curve.png")
