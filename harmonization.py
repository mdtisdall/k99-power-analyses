"""Power to show that harmonization A leaves less between-scanner disagreement in
regional age slopes than harmonization B.

See docs/harmonization-derivation.md for the derivation. Every subject is scanned
n_rep times at each of 3 sites; each site's value is the average of those scans.
Both approaches are applied to the same scans. For each approach, region, and
pair of sites, the age-slope difference between the two sites is the slope of the
per-subject between-site difference regressed on centered age. An approach's
disagreement in a region is the sum of squared slope differences over the 3 site
pairs (3 times the variance of that region's age slope across sites).

Primary test: A's disagreement summed over the 5 regions is smaller than B's
(one-sided, level alpha). Secondary tests: the same comparison in each region,
Holm-corrected across regions. Each test removes the bias that sampling noise adds
to squared slope differences, and uses a plug-in standard error.

For each N, this script reports the power of the primary test in the planning
scenario (gain_A, gain_B), the chance that each region's secondary test succeeds,
and the smallest gain mismatch in B that the primary test detects at the target
power. It saves the table to harmonization_curve.csv and plots the curves in
harmonization_curve.png.

The outcome is mean magnetic susceptibility (ppb) from QSM. age_slope and
sd_between are from Li et al. (2023), NeuroImage 269:119923; the noise values are
from Naji et al. (2022), NMR Biomed 35:e4788 (see the README).
"""

import numpy as np
from scipy import stats

from plotting import plot_curves

N_REGIONS, N_SITES = 5, 3
PAIRS = [(0, 1), (0, 2), (1, 2)]

# ---- Parameters -------------------------------------------------------------

N_grid = [10, 20, 30, 40, 50, 75, 100]  # numbers of subjects (each scanned at all 3 sites)
target_power = 0.80                # required power of the primary test
alpha = 0.05                       # one-sided level of each test
age_slope = 1.0                    # true age slope in every region (ppb/year), Li et al. 2023
sd_between = 20.0                  # SD (ppb) of true susceptibility across same-age subjects
age_min, age_max = 25, 65          # age distribution: uniform(age_min, age_max)
n_rep = 2                          # scans per subject at each site, averaged before analysis
sd_scan = 1.0                      # SD of a whole-scan offset (ppb; shared by all regions in one scan)
sd_noise = 2.1                     # SD of region-level scan-to-scan noise (ppb)
sd_site = 2.9                      # SD of region-level subject-by-site deviation that repeats
                                   # in every scan at that site (ppb); n_rep does not reduce it
noise_corr = 0.9                   # correlation between A's and B's noise in the same scans
gain_A = np.ones((N_REGIONS, N_SITES))  # planning scenario: each approach's remaining gain
gain_B = np.ones((N_REGIONS, N_SITES))  # (scale) at each site, by region (rows) and site
gain_B[:, 1], gain_B[:, 2] = 0.80, 1.20  # (columns); 1 everywhere = perfectly harmonized
n_sims = 10000
seed = 1

# ---- Simulation -------------------------------------------------------------

# Disagreement in one region is g' K g, where g = (gamma_12, gamma_13) are the slope
# differences for site pairs (1,2) and (1,3); gamma_23 = gamma_13 - gamma_12.
K = np.array([[2., -1.], [-1., 2.]])
Q_REGION = np.zeros((4, 4))
Q_REGION[:2, :2], Q_REGION[2:, 2:] = -K, K     # z = (A12, A13, B12, B13): B minus A
Q_ALL = np.kron(np.eye(N_REGIONS), Q_REGION)


def disagreement(gain, slope=1.0):
    """True disagreement in each region: sum over site pairs of squared slope
    differences, (gain_j - gain_k)^2 * slope^2."""
    gain = np.asarray(gain, dtype=float)
    return slope ** 2 * sum((gain[:, j] - gain[:, k]) ** 2 for j, k in PAIRS)


def largest_gap(gain, slope=1.0):
    """Largest between-site slope difference over regions and site pairs."""
    gain = np.asarray(gain, dtype=float)
    return slope * max(np.abs(gain[:, j] - gain[:, k]).max() for j, k in PAIRS)


def _check(N, gain_A, gain_B, n_rep, noise_corr, age_min, age_max):
    if int(N) != N or N < 3:
        raise ValueError(f"N must be a whole number of at least 3 (got {N}); "
                         f"each regression has N - 2 df.")
    if int(n_rep) != n_rep or n_rep < 1:
        raise ValueError(f"n_rep must be a whole number of at least 1 (got {n_rep}).")
    if not -1 <= noise_corr <= 1:
        raise ValueError("noise_corr must be between -1 and 1.")
    if not age_max > age_min:
        raise ValueError("age_max must be greater than age_min")
    for name, g in (("gain_A", gain_A), ("gain_B", gain_B)):
        if np.shape(g) != (N_REGIONS, N_SITES):
            raise ValueError(f"{name} must have shape {(N_REGIONS, N_SITES)} "
                             f"(regions x sites), got {np.shape(g)}.")


def z_stats(N, gain_A=gain_A, gain_B=gain_B, slope=1.0, sd_between=20.0, sd_scan=1.0,
            sd_noise=2.1, sd_site=2.9, n_rep=2, noise_corr=0.9, age_min=25, age_max=65,
            n_sims=10000, rng=None, chunk=1000):
    """Test statistics for n_sims simulated studies: (n_sims, 5) for the regions
    and (n_sims,) for the pooled test. Large values favor 'A disagrees less'."""
    _check(N, gain_A, gain_B, n_rep, noise_corr, age_min, age_max)
    N, n_rep = int(N), int(n_rep)
    gain_A, gain_B = np.asarray(gain_A, float), np.asarray(gain_B, float)
    rng = np.random.default_rng(rng)
    # Averaging n_rep scans shrinks the scan-to-scan terms by sqrt(n_rep) but not
    # the subject-by-site term. Sums of normals are normal, so simulating the
    # averages directly is exact.
    sd_w = sd_scan / np.sqrt(n_rep)
    sd_e = np.sqrt(sd_noise ** 2 / n_rep + sd_site ** 2)

    def noise(n):
        return (rng.normal(0, sd_w, (n, N, 1, N_SITES))             # whole-scan offsets
                + rng.normal(0, sd_e, (n, N, N_REGIONS, N_SITES)))  # region-level noise

    z_reg, z_all = [], []
    for start in range(0, n_sims, chunk):              # chunks bound memory use
        n = min(chunk, n_sims - start)
        age = rng.uniform(age_min, age_max, (n, N))
        x = age - age.mean(axis=1, keepdims=True)
        chi = slope * x[..., None] + rng.normal(0, sd_between, (n, N, N_REGIONS))
        e_B = noise(n)
        e_A = noise_corr * e_B + np.sqrt(1 - noise_corr ** 2) * noise(n)
        # Constant site offsets cancel in the differences, so they are left out.
        y_A = gain_A * chi[..., None] + e_A
        y_B = gain_B * chi[..., None] + e_B
        D = np.stack([y_A[..., 0] - y_A[..., 1], y_A[..., 0] - y_A[..., 2],
                      y_B[..., 0] - y_B[..., 1], y_B[..., 0] - y_B[..., 2]],
                     axis=-1).reshape(n, N, 4 * N_REGIONS)

        # OLS of each difference score on centered age; z holds the 20 slopes.
        sxx = (x ** 2).sum(axis=1)[:, None]
        z = np.einsum("ni,nik->nk", x, D) / sxx
        resid = D - D.mean(axis=1, keepdims=True) - z[:, None, :] * x[..., None]
        cov = np.einsum("nia,nib->nab", resid, resid) / (N - 2) / sxx[..., None]

        def stat(zz, cc, Q):
            # Estimate of (B's minus A's) disagreement, with the noise bias
            # tr(Q cov) removed, divided by its plug-in standard error.
            QC = Q @ cc
            Qz = zz @ Q
            est = np.einsum("na,na->n", Qz, zz) - np.trace(QC, axis1=1, axis2=2)
            var = (2 * np.einsum("nab,nba->n", QC, QC)
                   + 4 * np.einsum("na,nab,nb->n", Qz, cc, Qz))
            return est / np.sqrt(var)

        z_reg.append(np.stack([stat(z[:, 4 * r:4 * r + 4], cov[:, 4 * r:4 * r + 4,
                                                                4 * r:4 * r + 4], Q_REGION)
                               for r in range(N_REGIONS)], axis=1))
        z_all.append(stat(z, cov, Q_ALL))
    return np.concatenate(z_reg), np.concatenate(z_all)


def holm(z_reg, alpha=.05):
    """Which regions the Holm procedure declares, from one-sided z statistics."""
    p = stats.norm.sf(z_reg)
    order = np.argsort(p, axis=1)
    p_sorted = np.take_along_axis(p, order, axis=1)
    m = p.shape[1]
    ok = np.cumprod(p_sorted <= alpha / (m - np.arange(m)), axis=1).astype(bool)
    declared = np.zeros_like(ok)
    np.put_along_axis(declared, order, ok, axis=1)
    return declared


def power(N, alpha=.05, **kw):
    """Power of the pooled test, and each region's chance of being declared by
    the Holm-corrected secondary tests."""
    z_reg, z_all = z_stats(N, **kw)
    return (z_all > stats.norm.isf(alpha)).mean(), holm(z_reg, alpha).mean(axis=0)


def detectable_scale(N, target_power=.80, alpha=.05, gain_A=gain_A, gain_B=gain_B,
                     rng=None, tol=.01, **kw):
    """Smallest k such that B's gains gain_A + k (gain_B - gain_A) give the pooled
    test target_power, or inf if no k keeping every gain positive does. The same
    random draws are used for every k."""
    gain_A, gain_B = np.asarray(gain_A, float), np.asarray(gain_B, float)
    shrink = gain_B < gain_A                       # gains that fall as k grows
    k_max = (gain_A[shrink] / (gain_A - gain_B)[shrink]).min() if shrink.any() else 64.0

    def pw(k):
        z = z_stats(N, gain_A=gain_A, gain_B=gain_A + k * (gain_B - gain_A),
                    rng=rng, **kw)[1]
        return (z > stats.norm.isf(alpha)).mean()

    lo, hi = 0.0, min(1.0, 0.99 * k_max)
    while pw(hi) < target_power:
        if hi >= 0.99 * k_max:
            return np.inf
        lo, hi = hi, min(2 * hi, 0.99 * k_max)
    while hi - lo > tol * hi:
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if pw(mid) < target_power else (lo, mid)
    return hi


# ---- Power curve over N -----------------------------------------------------

if __name__ == "__main__":
    kw = dict(gain_A=gain_A, gain_B=gain_B, slope=age_slope, sd_between=sd_between,
              sd_scan=sd_scan, sd_noise=sd_noise, sd_site=sd_site, n_rep=n_rep,
              noise_corr=noise_corr, age_min=age_min, age_max=age_max, n_sims=n_sims)
    gap_B = largest_gap(gain_B, age_slope)
    rows = []
    for N in N_grid:
        # Seeded per N, so a row does not change when N_grid is edited.
        p_all, p_reg = power(N, alpha=alpha, rng=[seed, N], **kw)
        k = detectable_scale(N, target_power, alpha, rng=[seed, N, 1], **kw)
        rows.append((N, p_all, p_reg.min(), p_reg.mean(),
                     largest_gap(gain_A + k * (gain_B - gain_A), age_slope)
                     if np.isfinite(k) else np.nan))

    print(f"Disagreement (sum of squared between-site slope differences, ppb^2/year^2) "
          f"by region: A {np.round(disagreement(gain_A, age_slope), 3)}, "
          f"B {np.round(disagreement(gain_B, age_slope), 3)}")
    print(f"Largest between-site slope difference in B: {gap_B:.3g} ppb/year "
          f"({gap_B / age_slope:.0%} of the age slope)")
    print(f"Power of the pooled test, chance each region is declared (Holm), and the "
          f"smallest such largest slope difference in B (ppb/year) with "
          f"{target_power:.0%} power:")
    print("   N  Pooled power  Region (min)  Region (mean)  Smallest B difference")
    for N, p, rmin, rmean, d in rows:
        print(f"{N:4d}  {p:12.3f}  {rmin:12.3f}  {rmean:13.3f}  "
              + (f"{d:21.3g}" if np.isfinite(d) else f"{'not reached':>21}"))
    np.savetxt("harmonization_curve.csv", rows, delimiter=",",
               header="N,power_pooled,power_region_min,power_region_mean,"
                      "smallest_B_slope_difference",
               comments="", fmt=["%d", "%.4f", "%.4f", "%.4f", "%.6g"])

    scans = f"{n_rep} scan{'s' * (n_rep > 1)}"
    a_note = "; A perfectly harmonized" if not disagreement(gain_A).any() else ""
    plot_curves(
        [r[0] for r in rows], f"Subjects ({scans} at each of 3 sites)",
        [dict(series=[("Pooled over regions", [r[1] for r in rows]),
                      ("One region (Holm)", [r[3] for r in rows])],
              title=f"Power, B's slopes differ by up to {gap_B:g} ppb/year",
              ylim=(0, 1), ref=(target_power, f"{target_power:.0%} target")),
         dict(y=[r[4] for r in rows],
              title=f"Smallest such difference at {target_power:.0%} power (ppb/year)",
              ylim=(0, 1.08 * np.nanmax([r[4] for r in rows])))],
        f"A vs B: total between-site slope disagreement over 5 regions, one-sided "
        f"test at α = {alpha:g}{a_note}",
        "harmonization_curve.png")
    print("\nSaved harmonization_curve.csv and harmonization_curve.png")
