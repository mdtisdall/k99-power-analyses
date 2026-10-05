"""Smallest detectable group-by-region effect in a multi-site study.

See docs/group-effect-derivation.md for the derivation. Each subject is scanned
once, at one site. The analysis model is

    y ~ 0 + region + region:age_c + site + region:group + (1 | subject)

and the test is each region's group difference (each non-reference group vs the
reference group), Bonferroni-corrected across regions. This script reports the
smallest group difference detectable at the target power, in outcome units and
as Cohen's d.

All parameter values below are placeholders. Replace them with your design and
with estimates from prior data.
"""

import numpy as np
from scipy import optimize, stats

# ---- Parameters -------------------------------------------------------------

n_sites = 20                       # number of sites
subjects_per_site = 10             # subjects at each site, split evenly across groups
n_groups = 2                       # group 0 is the reference group
n_regions = 10                     # number of regions
sd_total = 1.0                     # SD of one region's value across subjects (same group, site, age)
icc = 0.5                          # correlation between two regions of the same subject
age_min, age_max = 25, 65          # age distribution: uniform(age_min, age_max)
alpha = 0.05                       # two-sided, family-wise across the Bonferroni family
bonferroni = True                  # False: test a single pre-specified region at alpha
target_power = 0.80
n_designs = 1000                   # random age draws to average power over
seed = 1

# ---- Design and standard error ----------------------------------------------


def make_design(n_sites, subjects_per_site, n_groups, age_min, age_max, rng):
    """Site, group, and centered age for every subject. Groups are allocated
    round-robin within each site, so they are as balanced as possible."""
    site = np.repeat(np.arange(n_sites), subjects_per_site)
    group = np.tile(np.arange(subjects_per_site) % n_groups, n_sites)
    age = rng.uniform(age_min, age_max, site.size)
    return site, group, age - age.mean()


def dummies(codes, n):
    """Indicator columns for levels 1..n-1 (level 0 is the reference)."""
    return (codes[:, None] == np.arange(1, n)).astype(float)


def group_se(site, group, age_c, n_sites, n_groups, n_regions, sd_total, icc):
    """SE and Satterthwaite df of each region's group-k-vs-reference difference,
    for k = 1..n_groups-1. Returns arrays of length n_groups - 1."""
    N = site.size
    tau2, sigma2 = icc * sd_total**2, (1 - icc) * sd_total**2
    G = dummies(group, n_groups)

    # Between-subject stratum: subject means on intercept, age, group, site.
    Zb = np.column_stack([np.ones(N), age_c, G, dummies(site, n_sites)])
    # Within-subject stratum: regional deviations on intercept, age, group.
    Zw = np.column_stack([np.ones(N), age_c, G])
    k = slice(2, 2 + n_groups - 1)
    cb = np.diag(np.linalg.inv(Zb.T @ Zb))[k]
    cw = np.diag(np.linalg.inv(Zw.T @ Zw))[k]

    vb = (tau2 + sigma2 / n_regions) * cb
    vw = sigma2 * (1 - 1 / n_regions) * cw
    dfb = N - Zb.shape[1]
    dfw = (n_regions - 1) * (N - Zw.shape[1])
    se2 = vb + vw
    df = se2**2 / (vb**2 / dfb + (vw**2 / dfw if n_regions > 1 else 0))
    return np.sqrt(se2), df


def power_curve_inputs(n_sites, subjects_per_site, n_groups, n_regions, sd_total,
                       icc, age_min, age_max, n_designs, rng):
    """SE and df for each contrast, for n_designs random age draws."""
    out = [group_se(*make_design(n_sites, subjects_per_site, n_groups,
                                 age_min, age_max, rng),
                    n_sites, n_groups, n_regions, sd_total, icc)
           for _ in range(n_designs)]
    se = np.array([o[0] for o in out])          # (n_designs, n_groups - 1)
    df = np.array([o[1] for o in out])
    return se, df


def power(effect, se, df, alpha_test):
    """Power of a two-sided t test of one region's group difference,
    averaged over designs."""
    tc = stats.t.ppf(1 - alpha_test / 2, df)
    ncp = effect / se
    return (stats.nct.sf(tc, df, ncp) + stats.nct.cdf(-tc, df, ncp)).mean()


def detectable_effect(se, df, alpha_test, target_power):
    """Smallest group difference with power >= target_power."""
    hi = 20 * se.max()
    return optimize.brentq(lambda e: power(e, se, df, alpha_test) - target_power,
                           0, hi, xtol=1e-10 * hi)


# ---- Detectable group-by-region effect --------------------------------------

if __name__ == "__main__":
    rng = np.random.default_rng(seed)
    n_tests = n_regions * (n_groups - 1) if bonferroni else 1
    alpha_test = alpha / n_tests
    se, df = power_curve_inputs(n_sites, subjects_per_site, n_groups, n_regions,
                                sd_total, icc, age_min, age_max, n_designs, rng)

    N = n_sites * subjects_per_site
    print(f"{n_sites} sites x {subjects_per_site} subjects = {N} subjects, "
          f"{n_groups} groups, {n_regions} regions")
    print(f"{n_tests} tests, each two-sided at alpha = {alpha_test:.4g}; "
          f"target power {target_power:.0%}")
    print()
    print("Contrast         SE (outcome units)  df      Detectable effect  Cohen's d")
    for k in range(n_groups - 1):
        effect = detectable_effect(se[:, k], df[:, k], alpha_test, target_power)
        print(f"group {k + 1} vs 0      {se[:, k].mean():<18.4g}  "
              f"{df[:, k].mean():<6.0f}  {effect:<17.4g}  {effect / sd_total:.3f}")
