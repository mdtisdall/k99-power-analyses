# k99-power-analyses

Two power analyses for multi-site MRI studies:

1. [**Cross-site equivalence of age effects**](#analysis-1-cross-site-equivalence-of-age-effects)
   ([`equivalence.py`](equivalence.py)): every subject is scanned at all 3 sites.
   How small a between-site difference in regional age effects can we rule out?
2. [**Detectable group-by-region effect**](#analysis-2-detectable-group-by-region-effect)
   ([`group_effect.py`](group_effect.py)): each subject is scanned at one of many
   sites. How large a group difference in a region can we detect?

## Quick start

1. Install the dependencies (numpy, scipy, and matplotlib for plots). Any Python 3, including the one
   built into macOS, works:

   ```bash
   python3 -m pip install -r requirements.txt
   ```

2. Open the script for the analysis you want and edit the **Parameters** block at
   the top (see below).

3. Run it:

   ```bash
   python3 equivalence.py
   ```

   ```bash
   python3 group_effect.py
   ```

   Each run takes a few seconds. On Windows, type `python` instead of `python3`.

(With Nix, `nix-shell --run "python3 equivalence.py"` does steps 1 and 3 together,
and likewise for `group_effect.py`.)

All parameter values in both scripts are **placeholders**. Replace them with your
design and with estimates from prior data.

---

## Analysis 1: cross-site equivalence of age effects

N subjects are each scanned at all 3 sites, and the outcome is measured in 5
regions. Each site's data are analyzed separately with

```
y ~ 0 + region + region:age_c + (1 | subject)
```

For a given number of subjects, [`equivalence.py`](equivalence.py) estimates the
**smallest difference in age effects between sites that we can rule out**, i.e. how
tightly we can show that the age effects are "the same" across sites.

### Parameters

Estimate the noise parameters from prior data, ideally traveling-subject or
test-retest scans.

| Parameter | What it is |
|---|---|
| `N_grid` | Number(s) of subjects to evaluate, e.g. `[5]` or `[5, 10, 20]`. |
| `target_power` | Required probability that the study shows equivalence. Usually 0.80. |
| `s_ss` | SD of a whole-scan offset: how much a subject's values shift together, across all regions, from one scan to another. |
| `sigma` | Measurement-noise SD for one region in one scan. The same for all sites and regions. |
| `age_min`, `age_max` | Age range of the sample. Ages are drawn uniformly from this range. |
| `dslope` | Optional true differences in age slope between sites (5 regions × 3 sites). Leave at 0 for the standard calculation. |
| `alpha` | Level of each one-sided test. Leave at 0.05 (gives 90% CIs). |
| `n_sims` | Simulations per sample size. 10000 gives Δ to about ±1%. |
| `seed` | Random seed, so results are reproducible. |

You don't need between-subject variance, regional means, the true age slopes, or
fixed site offsets, because they cancel out.

### Output: Δ

The script prints Δ for each number of subjects and saves the table to
`detectable_delta.csv`. Δ is the **equivalence margin**: a between-site difference
in age slope, in outcome units per year. With the given number of subjects, a study
has `target_power` chance of showing that every between-site slope difference lies
within ±Δ. Smaller is better. Compare Δ with the expected age slope: a Δ around 20%
of the slope would be a convincing "same", while a Δ as large as the slope means the
sites could differ by the whole age effect. When you analyze real data, choose the
margin before you see the data.

### How it works

To show that the age effects are the *same*, the script uses equivalence tests
(TOST) rather than a non-significant site × age interaction, which can't show
sameness. For each region (5) and each pair of sites (3), it regresses the
per-subject difference between the two sites on age. That slope is exactly the
difference between the two sites' age effects. The 90% CI for that slope must fall
inside (−Δ, +Δ), and the sites count as "the same" only if **all 15** comparisons
pass. The script simulates studies and reports the smallest Δ for which that
happens with probability `target_power`.

Full derivation and assumptions: [docs/equivalence-derivation.md](docs/equivalence-derivation.md).

---

## Analysis 2: detectable group-by-region effect

Subjects are recruited at many sites, and each subject is scanned **once, at one
site**. The outcome is measured in several regions, and subjects belong to one of
several groups (e.g. patients and controls). The analysis model is

```
y ~ 0 + region + region:age_c + site + region:group + (1 | subject)
```

The group-by-region effect is the difference between a group and the reference
group in one region. Each region's effect is tested separately, Bonferroni-corrected
across regions (and across groups, if there are more than 2).
For each number of subjects per site, [`group_effect.py`](group_effect.py) reports
the **power** to detect a given group difference and the **smallest group difference
detectable** at the target power.

> **Balanced-sites simplification.** For the power calculation, every site is
> assumed to recruit the same number of subjects from each group. This is a
> planning simplification, not a requirement of the analysis. It gives the best
> case for a given total N. Groups that are unbalanced within sites need a somewhat
> larger effect to reach the same power, typically a few percent larger. Use values
> of `subjects_per_site` divisible by `n_groups`; otherwise the script marks those
> rows and warns that the simplification holds only approximately.

### Parameters

| Parameter | Default | What it is |
|---|---|---|
| `n_sites` | 20 | Number of sites. |
| `subjects_per_site` | 2, 4, …, 20 | List of subjects-per-site values to evaluate (the x-axis of the power curve). Each is split equally across groups (balanced-sites simplification). |
| `n_groups` | 2 | Number of groups. Group 0 is the reference; each other group is compared with it. |
| `n_regions` | 10 | Number of regions (and of Bonferroni-corrected tests, with 2 groups). |
| `sd_total` | 1.0 | SD of one region's value across subjects with the same group, site, and age. Only scales the answer in outcome units; leave at 1 to get the answer in SD units. |
| `icc` | 0.5 | Correlation between two regions of the same subject. Has little effect under the balanced-sites simplification. |
| `age_min`, `age_max` | 25, 65 | Age range. Ages are drawn uniformly from this range, independently of group. |
| `alpha` | 0.05 | Family-wise two-sided level, split across the tests by Bonferroni. |
| `bonferroni` | `True` | Set to `False` to test a single pre-specified region at `alpha`. |
| `target_power` | 0.80 | Required power for each region's test. |
| `effect_d` | 0.5 | Group difference, as Cohen's d, at which the power curve is computed. |
| `n_designs` | 1000 | Random age draws averaged over. |
| `seed` | 1 | Random seed, so results are reproducible. |

To estimate `sd_total` and `icc` from prior data, fit the same model and take
`sd_total` = √(τ² + σ²) and `icc` = τ² / (τ² + σ²), where τ² is the subject variance
and σ² the residual variance.

### Output

For each value of `subjects_per_site`, the script prints the total N, the degrees
of freedom, the **power at `effect_d`**, and the **smallest detectable group
difference** (as Cohen's d and in outcome units). It saves the table to
`group_effect_curve.csv` and plots both curves in `group_effect_curve.png`:

- left: power to detect `effect_d` against subjects per site, with the target power
  marked;
- right: smallest detectable d at the target power against subjects per site.

With the defaults (20 sites, 2 groups, 10 regions, so each test is at α = 0.005):

| Subjects per site | N | Power at d = 0.5 | Smallest detectable d |
|---|---|---|---|
| 4 | 80 | 0.27 | 0.83 |
| 8 | 160 | 0.63 | 0.58 |
| 10 | 200 | 0.76 | 0.52 |
| 12 | 240 | 0.85 | 0.47 |
| 16 | 320 | 0.95 | 0.41 |
| 20 | 400 | 0.99 | 0.37 |

With more than 2 groups, each row reports the least powerful comparison with the
reference group.

### How it works

The script uses an exact closed-form standard error instead of simulating, so it
runs in seconds. With one scan per subject, the model's random subject intercept
splits each subject's data into a subject mean and deviations from that mean, and
the mixed model's estimate of each group-by-region effect combines an ordinary
least-squares fit to each.

A practical consequence of the balanced-sites simplification: the answer is
essentially that of a two-sample t test on one region (with two equal groups,
SE ≈ 2·SD/√N). The number of sites, the number of regions (except through
Bonferroni), and `icc` barely matter.

Full derivation, simulation check, and assumptions:
[docs/group-effect-derivation.md](docs/group-effect-derivation.md).
