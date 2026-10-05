"""Smallest detectable group-by-region effect in a multi-site study.

See docs/group-effect-derivation.md for the derivation. Each subject is scanned
once, at one site. The analysis model is

    y ~ 0 + region + region:age_c + site + region:group + (1 | subject)

and the test is each region's group difference (each non-reference group vs the
reference group), Bonferroni-corrected across regions. The power calculation
matches a Wald t test with Satterthwaite degrees of freedom, as reported by
lmerTest (R) for this model; z-based p-values are too liberal at small df.

Balanced-sites simplification: for the power calculation, every site recruits
the same number of subjects from each group. This is the best case for a given
total N; imbalance within sites raises the detectable effect.

For each number of subjects per site, this script reports the power to detect a
group difference of effect_d (Cohen's d), and the smallest group difference
detectable at the target power, in outcome units and as Cohen's d. It saves the
table to group_effect_curve.csv and plots both curves in group_effect_curve.png.

All parameter values below are placeholders. Replace them with your design and
with estimates from prior data.
"""

import numpy as np
from scipy import optimize, stats

from plotting import plot_curves

# ---- Parameters -------------------------------------------------------------

n_sites = 20                       # number of sites
subjects_per_site = [2, 4, 6, 8, 10, 12, 14, 16, 18, 20]  # values to evaluate; each split equally across groups
n_groups = 2                       # group 0 is the reference group
n_regions = 10                     # number of regions
sd_total = 1.0                     # SD of one region's value across subjects (same group, site, age)
region_corr = 0.5                  # correlation between two regions of the same subject
age_min, age_max = 25, 65          # age distribution: uniform(age_min, age_max)
alpha = 0.05                       # two-sided, family-wise across the Bonferroni family
bonferroni = True                  # False: one pre-specified test (one region, one group) at alpha
target_power = 0.80
effect_d = 0.5                     # group difference (Cohen's d) for the power curve
n_designs = 1000                   # random age draws to average power over
seed = 1

# ---- Design and standard error ----------------------------------------------


def make_design(n_sites, subjects_per_site, n_groups, age_min, age_max, rng):
    """Site, group, and centered age for every subject. Groups are allocated
    round-robin across all subjects, site by site: exactly balanced at every
    site (the balanced-sites simplification) when subjects_per_site is divisible
    by n_groups; otherwise each site is as even as possible and the extra
    subjects rotate between groups, so overall group sizes stay as equal as
    possible."""
    site = np.repeat(np.arange(n_sites), subjects_per_site)
    group = np.arange(site.size) % n_groups
    age = rng.uniform(age_min, age_max, site.size)
    return site, group, age - age.mean()


def check_design(n_sites, subjects_per_site, n_groups, n_regions, region_corr,
                 age_min, age_max):
    """Raise a clear error for designs the model cannot be fit to."""
    N = n_sites * subjects_per_site
    problems = []
    if n_sites < 1 or subjects_per_site < 1 or n_regions < 1:
        problems.append("n_sites, subjects_per_site, and n_regions must be at least 1")
    if n_groups < 2:
        problems.append("n_groups must be at least 2")
    if N - (n_groups + n_sites) < 1:
        problems.append(f"{n_sites} sites x {subjects_per_site} subjects leaves no "
                        f"degrees of freedom for {n_groups} groups and {n_sites} site "
                        f"effects (need N > n_groups + n_sites)")
    if not age_max > age_min:
        problems.append("age_max must be greater than age_min")
    if not 0 <= region_corr <= 1:
        problems.append("region_corr must be between 0 and 1")
    if problems:
        raise ValueError("; ".join(problems))


def dummies(codes, n):
    """Indicator columns for levels 1..n-1 (level 0 is the reference)."""
    return (codes[:, None] == np.arange(1, n)).astype(float)


def group_se(site, group, age_c, n_sites, n_groups, n_regions, sd_total,
             region_corr):
    """SE and Satterthwaite df of each region's group-k-vs-reference difference,
    for k = 1..n_groups-1. Returns arrays of length n_groups - 1."""
    N = site.size
    tau2, sigma2 = region_corr * sd_total**2, (1 - region_corr) * sd_total**2
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
                       region_corr, age_min, age_max, n_designs, rng):
    """SE and df for each contrast, for n_designs random age draws."""
    check_design(n_sites, subjects_per_site, n_groups, n_regions, region_corr,
                 age_min, age_max)
    out = [group_se(*make_design(n_sites, subjects_per_site, n_groups,
                                 age_min, age_max, rng),
                    n_sites, n_groups, n_regions, sd_total, region_corr)
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


# ---- Power curve over subjects per site -------------------------------------


def curve(n_sites, subjects_per_site, n_groups, n_regions, sd_total, region_corr,
          age_min, age_max, alpha_test, target_power, effect_d, n_designs, seed):
    """One row per subjects-per-site value. With more than 2 groups, each row
    reports the least powerful comparison with the reference group."""
    rows = []
    for n in subjects_per_site:
        # Seeded per value, so a row does not change when the list is edited.
        rng = np.random.default_rng([seed, n])
        se, df = power_curve_inputs(n_sites, n, n_groups, n_regions, sd_total,
                                    region_corr, age_min, age_max, n_designs, rng)
        k = np.argmax(se.mean(axis=0))
        rows.append(dict(
            subjects_per_site=n, N=n_sites * n, df=df[:, k].mean(),
            power=power(effect_d * sd_total, se[:, k], df[:, k], alpha_test),
            detectable=detectable_effect(se[:, k], df[:, k], alpha_test,
                                         target_power)))
    return rows


# ---- Run --------------------------------------------------------------------

if __name__ == "__main__":
    n_tests = n_regions * (n_groups - 1) if bonferroni else 1
    alpha_test = alpha / n_tests
    rows = curve(n_sites, subjects_per_site, n_groups, n_regions, sd_total, region_corr,
                 age_min, age_max, alpha_test, target_power, effect_d,
                 n_designs, seed)

    summary = (f"{n_sites} sites, {n_groups} groups, {n_regions} regions; "
               f"{n_tests} tests, each two-sided at alpha = {alpha_test:.4g}")
    print(summary)
    unbalanced = [n for n in subjects_per_site if n % n_groups]
    if unbalanced:
        print(f"Note: {unbalanced} subjects per site can't be split equally across "
              f"{n_groups} groups, so the balanced-sites simplification holds only "
              f"approximately for those rows (marked *): each site is as even as "
              f"possible, with the extra subjects rotating between groups.")
    else:
        print("Balanced-sites simplification: equal group sizes at every site")
    if n_groups > 2:
        print("Each row shows the least powerful comparison with the reference group.")
    print()
    print(f"Subjects/site     N    df  Power at d={effect_d:<5g}  "
          f"Detectable d  Detectable (outcome units)")
    for r in rows:
        mark = "*" if r["subjects_per_site"] in unbalanced else " "
        print(f"{r['subjects_per_site']:>12}{mark} {r['N']:>5} {r['df']:>5.0f}  "
              f"{r['power']:>15.3f}  {r['detectable'] / sd_total:>12.3f}  "
              f"{r['detectable']:>26.4g}")

    with open("group_effect_curve.csv", "w") as f:
        f.write("subjects_per_site,N,df,power_at_effect_d,detectable_d,"
                "detectable_outcome_units\n")
        for r in rows:
            f.write(f"{r['subjects_per_site']},{r['N']},{r['df']:.1f},"
                    f"{r['power']:.4f},{r['detectable'] / sd_total:.4f},"
                    f"{r['detectable']:.6g}\n")
    plot_curves(
        subjects_per_site, f"Subjects per site ({n_sites} sites)",
        [dict(y=[r["power"] for r in rows], title=f"Power to detect d = {effect_d:g}",
              ylim=(0, 1), ref=(target_power, f"{target_power:.0%} target")),
         dict(y=[r["detectable"] / sd_total for r in rows],
              title=f"Smallest detectable d at {target_power:.0%} power")],
        summary, "group_effect_curve.png")
    print("\nSaved group_effect_curve.csv and group_effect_curve.png")
