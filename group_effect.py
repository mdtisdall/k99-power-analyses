"""Smallest detectable group-by-region effect in a multi-site study.

See docs/group-effect-derivation.md for the derivation. Each subject is scanned
at one site, at every visit in visit_years (e.g. once a year for three years).
The analysis model is

    y ~ 0 + region + region:age_c + site + region:group
        + (1 | subject) + (1 | subject:visit) + (1 | subject:region)

and the test is each region's group difference (each non-reference group vs the
reference group, the same at every visit), Bonferroni-corrected across regions.
With one visit, the two extra random effects drop out. The power calculation
matches a Wald t test with Satterthwaite degrees of freedom, as reported by
lmerTest (R) for this model; z-based p-values are too liberal at small df.

Balanced-sites simplification: for the power calculation, every site recruits
the same number of subjects from each group. This is the best case for a given
total N; imbalance within sites raises the detectable effect.

For each number of subjects per site, this script reports the smallest group
difference detectable at the target power, as Cohen's d and, for each SD scenario
in sd_scenarios, in ppb; and the power to detect a group difference of effect_ppb
under each scenario. It saves the table to group_effect_curve.csv and plots both
curves, in ppb, in group_effect_curve.png.

The design values below (14 sites, 22 subjects per site, scanned once a year
for 3 years) are the planned study. The SD scenarios (10 and 30 ppb) are the low
and high values from Lancione et al. 2022, Naji et al. 2022, and Li et al. 2023
(see the README). The other noise values and effect_ppb are placeholders.
"""

import numpy as np
from scipy import optimize, stats

from plotting import plot_curves

# ---- Parameters -------------------------------------------------------------

n_sites = 14                       # number of sites
planned_per_site = 22              # marked on the plot (None for no mark)
subjects_per_site = [6, 10, 14, 18, 22, 26, 30]  # values to evaluate; each split as equally as possible
                                   # (22 per site = 308 subjects)
n_groups = 2                       # group 0 is the reference group
n_regions = 10                     # number of regions
visit_years = [0, 1, 2]            # when each subject is scanned, in years from baseline
                                   # (always at the same site); [0] for one scan
sd_scenarios = {"Low SD": 10.0,    # SD (ppb) of one region's value in one scan across
                "High SD": 30.0}   # subjects (same group, site, age); one curve each
region_corr = 0.5                  # correlation between two regions of the same scan
retest_corr = 0.8                  # correlation between two visits' values of one region in one subject
                                   # (after age); 1 - retest_corr is the visit-to-visit share of variance
visit_region_corr = 0.5            # correlation between two regions' visit-to-visit fluctuations
age_min, age_max = 25, 65          # baseline age distribution: uniform(age_min, age_max)
alpha = 0.05                       # two-sided, family-wise across the Bonferroni family
bonferroni = True                  # False: one pre-specified test (one region, one group) at alpha
target_power = 0.80
effect_ppb = 10.0                  # group difference (ppb) for the power curve
n_designs = 1000                   # random age draws to average power over
seed = 1

# ---- Design and standard error ----------------------------------------------


def make_design(n_sites, subjects_per_site, n_groups, age_min, age_max, rng):
    """Site, group, and centered baseline age for every subject. Groups are
    allocated round-robin across all subjects, site by site: exactly balanced at
    every site (the balanced-sites simplification) when subjects_per_site is
    divisible by n_groups; otherwise each site is as even as possible and the
    extra subjects rotate between groups, so overall group sizes stay as equal as
    possible. Every subject has the same visit schedule, so centered baseline age
    is also each subject's centered mean age over visits."""
    site = np.repeat(np.arange(n_sites), subjects_per_site)
    group = np.arange(site.size) % n_groups
    age = rng.uniform(age_min, age_max, site.size)
    return site, group, age - age.mean()


def variance_components(sd_total, region_corr, retest_corr, visit_region_corr):
    """Split one scan's variance sd_total**2 into four parts: the subject effect
    (shared by all regions and visits), the subject-by-region effect (stable
    across visits), the visit effect (shared by all regions in one scan), and
    the visit-by-region noise."""
    v = sd_total**2
    visit = (1 - retest_corr) * v
    visit_shared = visit_region_corr * visit
    visit_region = visit - visit_shared
    return dict(subject=region_corr * v - visit_shared,
                subject_region=(1 - region_corr) * v - visit_region,
                visit=visit_shared, visit_region=visit_region)


def check_design(n_sites, subjects_per_site, n_groups, n_regions, region_corr,
                 age_min, age_max, visit_years=(0,), retest_corr=1.0,
                 visit_region_corr=0.0):
    """Raise a clear error for designs the model cannot be fit to."""
    counts = dict(n_sites=n_sites, subjects_per_site=subjects_per_site,
                  n_groups=n_groups, n_regions=n_regions)
    problems = [f"{k} must be a whole number (got {v})"
                for k, v in counts.items() if int(v) != v]
    if not problems:
        N = n_sites * subjects_per_site
        if n_sites < 1 or subjects_per_site < 1 or n_regions < 1:
            problems.append("n_sites, subjects_per_site, and n_regions must be at least 1")
        elif n_groups < 2:
            problems.append("n_groups must be at least 2")
        elif N - (n_groups + n_sites) < 1:
            problems.append(f"{n_sites} sites x {subjects_per_site} subjects leaves no "
                            f"degrees of freedom for {n_groups} groups and {n_sites} "
                            f"site effects (need N > n_groups + n_sites)")
        else:
            # Each group must be comparable with the reference group within
            # sites; otherwise its effect is confounded with site.
            site, group, age_c = make_design(int(n_sites), int(subjects_per_site),
                                             int(n_groups), 0, 1,
                                             np.random.default_rng(0))
            Zb = np.column_stack([np.ones(N), age_c, dummies(group, n_groups),
                                  dummies(site, n_sites)])
            if np.linalg.matrix_rank(Zb) < Zb.shape[1]:
                problems.append(f"with {subjects_per_site} subjects per site and "
                                f"{n_groups} groups, some group effects can't be "
                                f"separated from site effects (a group never shares "
                                f"a site with the reference group); use more "
                                f"subjects per site")
    if not age_max > age_min:
        problems.append("age_max must be greater than age_min")
    if len(visit_years) < 1 or len(set(visit_years)) < len(visit_years):
        problems.append("visit_years must list at least one visit, with no repeats")
    corrs = dict(region_corr=region_corr, retest_corr=retest_corr,
                 visit_region_corr=visit_region_corr)
    bad = [k for k, v in corrs.items() if not 0 <= v <= 1]
    problems += [f"{k} must be between 0 and 1" for k in bad]
    if not bad and len(visit_years) > 1:
        vc = variance_components(1.0, region_corr, retest_corr, visit_region_corr)
        if vc["subject"] < -1e-12:
            problems.append("visit-to-visit variance shared by all regions, "
                            "visit_region_corr * (1 - retest_corr), can't exceed "
                            "region_corr")
        if vc["subject_region"] < -1e-12:
            problems.append("region-specific visit-to-visit variance, "
                            "(1 - visit_region_corr) * (1 - retest_corr), can't "
                            "exceed 1 - region_corr")
    if problems:
        raise ValueError("; ".join(problems))


def dummies(codes, n):
    """Indicator columns for levels 1..n-1 (level 0 is the reference)."""
    return (codes[:, None] == np.arange(1, n)).astype(float)


def _stratum_pair(Z, k, lam_mean, lam_dev, kappa):
    """Variance of the group coefficients k from a time-mean stratum (OLS of
    subject values on Z, error variance lam_mean) combined with the matching
    time-deviation stratum, which adds information kappa / lam_dev on the age
    slope (column 1 of Z) only. Also returns the derivatives of that variance
    with respect to lam_mean and lam_dev, for the Satterthwaite df."""
    A = np.linalg.inv(Z.T @ Z)
    c, a, b = np.diag(A)[k], A[1, k], A[1, 1]
    # Sherman-Morrison: [(Z'Z / lam_mean + (kappa / lam_dev) e e')^-1]_kk,
    # written so that kappa = 0 (one visit) and lam_dev = 0 are both safe.
    D = lam_dev + kappa * lam_mean * b
    r = np.divide(kappa * a**2, D**2, out=np.zeros_like(a), where=D > 0)
    v = lam_mean * c - lam_mean**2 * r * D
    d_mean = c - r * lam_mean * (2 * lam_dev + kappa * lam_mean * b)
    d_dev = r * lam_mean**2
    return v, d_mean, d_dev


def group_se(site, group, age_c, n_sites, n_groups, n_regions, sd_total,
             region_corr, visit_years=(0,), retest_corr=1.0, visit_region_corr=0.0):
    """SE and Satterthwaite df of each region's group-k-vs-reference difference,
    for k = 1..n_groups-1. Returns arrays of length n_groups - 1. age_c is each
    subject's centered mean age over visits."""
    N, R, T = site.size, n_regions, len(visit_years)
    vc = variance_components(sd_total, region_corr, retest_corr, visit_region_corr)
    # Error variances of the four strata (time mean or deviation x region mean
    # or deviation) of one subject's T x R values.
    lam = dict(mean_mean=vc["subject"] + vc["subject_region"] / R + vc["visit"] / T
               + vc["visit_region"] / (T * R),
               mean_dev=vc["subject_region"] + vc["visit_region"] / T,
               dev_mean=vc["visit"] + vc["visit_region"] / R,
               dev_dev=vc["visit_region"])
    kappa = N * np.var(visit_years) * T      # within-subject age sum of squares
    G = dummies(group, n_groups)
    k = np.arange(2, 2 + n_groups - 1)

    # Region means: subject means on intercept, age, group, site.
    Zb = np.column_stack([np.ones(N), age_c, G, dummies(site, n_sites)])
    vb, gb, gb_dev = _stratum_pair(Zb, k, lam["mean_mean"], lam["dev_mean"], kappa)
    # Regional deviations: on intercept, age, group, for each region.
    Zw = np.column_stack([np.ones(N), age_c, G])
    vw, gw, gw_dev = (x * (1 - 1 / R) for x in
                      _stratum_pair(Zw, k, lam["mean_dev"], lam["dev_dev"], kappa))

    # Satterthwaite: var(v_hat) = sum of (dv/dlam)^2 * 2 lam^2 / df over strata.
    terms = [(gb * lam["mean_mean"], N - Zb.shape[1]),
             (gb_dev * lam["dev_mean"], N * (T - 1)),
             (gw * lam["mean_dev"], (R - 1) * (N - Zw.shape[1])),
             (gw_dev * lam["dev_dev"], (R - 1) * N * (T - 1))]
    se2 = vb + vw
    df = se2**2 / sum(g**2 / dfs for g, dfs in terms if dfs > 0)
    return np.sqrt(se2), df


def power_curve_inputs(n_sites, subjects_per_site, n_groups, n_regions, sd_total,
                       region_corr, age_min, age_max, n_designs, rng,
                       visit_years=(0,), retest_corr=1.0, visit_region_corr=0.0):
    """SE and df for each contrast, for n_designs random age draws."""
    check_design(n_sites, subjects_per_site, n_groups, n_regions, region_corr,
                 age_min, age_max, visit_years, retest_corr, visit_region_corr)
    n_sites, subjects_per_site, n_groups, n_regions = map(
        int, (n_sites, subjects_per_site, n_groups, n_regions))
    out = [group_se(*make_design(n_sites, subjects_per_site, n_groups,
                                 age_min, age_max, rng),
                    n_sites, n_groups, n_regions, sd_total, region_corr,
                    visit_years, retest_corr, visit_region_corr)
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
    while power(hi, se, df, alpha_test) < target_power:   # very small df
        hi *= 2
    return optimize.brentq(lambda e: power(e, se, df, alpha_test) - target_power,
                           0, hi, xtol=1e-10 * hi)


# ---- Power curve over subjects per site -------------------------------------


def curve(n_sites, subjects_per_site, n_groups, n_regions, sd_total, region_corr,
          age_min, age_max, alpha_test, target_power, effect_d, n_designs, seed,
          visit_years=(0,), retest_corr=1.0, visit_region_corr=0.0):
    """One row per subjects-per-site value. With more than 2 groups, each row
    reports the least powerful comparison with the reference group."""
    rows = []
    for n in subjects_per_site:
        # Seeded per value, so a row does not change when the list is edited.
        rng = np.random.default_rng([seed, n])
        se, df = power_curve_inputs(n_sites, n, n_groups, n_regions, sd_total,
                                    region_corr, age_min, age_max, n_designs, rng,
                                    visit_years, retest_corr, visit_region_corr)
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
    # One curve per SD scenario. The designs are seeded identically, so the
    # detectable d is the same in every scenario; only its size in ppb differs.
    curves = {name: curve(n_sites, subjects_per_site, n_groups, n_regions, sd,
                          region_corr, age_min, age_max, alpha_test, target_power,
                          effect_ppb / sd, n_designs, seed, visit_years,
                          retest_corr, visit_region_corr)
              for name, sd in sd_scenarios.items()}
    rows = next(iter(curves.values()))
    labels = {name: f"{name} ({sd:g} ppb)" for name, sd in sd_scenarios.items()}

    planned = (planned_per_site, f"Planned (N = {n_sites * planned_per_site})"
               ) if planned_per_site else None
    n_visits = len(visit_years)
    summary = (f"{n_sites} sites, {n_groups} groups, {n_regions} regions, "
               f"{n_visits} visit{'s' if n_visits > 1 else ''} per subject; "
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
    sds = list(sd_scenarios.values())
    print(f"Subjects/site     N    df  Detectable d  " + "  ".join(
        f"{labels[k]:>32}" for k in curves))
    print(" " * 41 + "  ".join(f"{f'Power at {effect_ppb:g} ppb':>16}{'Detectable ppb':>16}"
                               for k in curves))
    for i, r in enumerate(rows):
        mark = "*" if r["subjects_per_site"] in unbalanced else " "
        print(f"{r['subjects_per_site']:>12}{mark} {r['N']:>5} {r['df']:>5.0f}  "
              f"{r['detectable'] / sds[0]:>12.3f}  " + "  ".join(
                  f"{c[i]['power']:>16.3f}{c[i]['detectable']:>16.2f}"
                  for c in curves.values()))

    with open("group_effect_curve.csv", "w") as f:
        f.write("subjects_per_site,N,df,detectable_d," + ",".join(
            f"sd_ppb_{k.lower().replace(' ', '_')},power_at_{effect_ppb:g}_ppb_"
            f"{k.lower().replace(' ', '_')},detectable_ppb_{k.lower().replace(' ', '_')}"
            for k in curves) + "\n")
        for i, r in enumerate(rows):
            f.write(f"{r['subjects_per_site']},{r['N']},{r['df']:.1f},"
                    f"{r['detectable'] / sds[0]:.4f}," + ",".join(
                        f"{sd:g},{c[i]['power']:.4f},{c[i]['detectable']:.4f}"
                        for sd, c in zip(sds, curves.values())) + "\n")
    plot_curves(
        subjects_per_site, f"Subjects per site ({n_sites} sites)",
        [dict(series=[(labels[k], [r["power"] for r in c]) for k, c in curves.items()],
              title=f"Power to detect a {effect_ppb:g} ppb group difference",
              ylim=(0, 1), ref=(target_power, f"{target_power:.0%} target"),
              legend_loc="center left", mark_x=planned),
         dict(series=[(labels[k], [r["detectable"] for r in c]) for k, c in curves.items()],
              title=f"Smallest detectable difference at {target_power:.0%} power",
              ylabel="ppb", legend_loc="center left", mark_x=planned)],
        summary, "group_effect_curve.png")
    print("\nSaved group_effect_curve.csv and group_effect_curve.png")
