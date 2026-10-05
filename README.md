# k99-power-analyses

For a given number of subjects, estimates the **smallest difference in age effects
between 3 MRI sites that we can rule out**, i.e. how tightly we can show that age
effects in 5 brain regions are "the same" across sites. Every subject is scanned at
all 3 sites.

## Quick start

1. Install the two dependencies (numpy and scipy). Any Python 3, including the one
   built into macOS, works:

   ```bash
   python3 -m pip install -r requirements.txt
   ```

2. Open [`power.py`](power.py) and edit the **Parameters** block at the top (see
   below).

3. Run it:

   ```bash
   python3 power.py
   ```

   It prints the smallest Δ that can be shown with the target power for each
   number of subjects, and saves the same table to `detectable_delta.csv`. A run
   takes a few seconds.

   On Windows, type `python` instead of `python3`.

(With Nix, `nix-shell --run "python3 power.py"` does steps 1 and 3 together.)

## Parameters to set

The values in `power.py` are **placeholders**. Replace them with estimates from
prior data, ideally traveling-subject or test-retest scans.

| Parameter | What it is |
|---|---|
| `N_grid` | Number(s) of subjects to evaluate, e.g. `[5]` or `[5, 10, 20]`. |
| `target_power` | Required probability that the study shows equivalence. Usually 0.80. |
| `s_ss` | SD of a whole-scan offset: how much a subject's values shift together, across all regions, from one scan to another. |
| `sigma` | Measurement-noise SD for one region in one scan, one value per site: `[site1, site2, site3]`. |
| `age_min`, `age_max` | Age range of the sample. Ages are drawn uniformly from this range. |
| `dslope` | Optional true differences in age slope between sites (5 regions × 3 sites). Leave at 0 for the standard calculation. |
| `alpha` | Test level. Leave at 0.05. |
| `n_sims` | Simulations per sample size. 10000 gives Δ to about ±1%. |
| `seed` | Random seed, so results are reproducible. |

You don't need between-subject variance, regional means, the true age slopes, or
fixed site offsets, because they cancel out (see below).

## Output: Δ

Δ is the **equivalence margin**: a between-site difference in age slope, in outcome
units per year. With the given number of subjects, a study has `target_power` chance
of showing that every between-site slope difference lies within ±Δ. Smaller is
better. To judge whether it's useful, compare Δ with the expected age slope. A Δ
around 20% of the slope would be a convincing "same". A Δ as large as the slope
itself means the sites could differ by the whole age effect.

When you analyze real data, choose the margin before you see the data.

## What it calculates

Each site's data are analyzed separately with the model

```
y ~ 0 + region + region:age_c + (1 | subject)
```

To show that the age effects are the *same*, the script uses equivalence tests
(TOST) rather than a non-significant site × age interaction, which can't show
sameness. For each region (5) and each pair of sites (3), it regresses the
per-subject difference between the two sites on age. That slope is exactly the
difference between the two sites' age effects. The 90% CI for that slope must fall
inside (−Δ, +Δ). The sites count as "the same" only if **all 15** comparisons pass,
and **power is the probability that all 15 pass**.

In each simulated study, the smallest Δ that study would pass is the largest
`|slope difference| + t × SE` across the 15 comparisons. The reported Δ is the
`target_power` quantile of that value across simulations.

Because each comparison uses differences within the same subject, anything stable
within a subject cancels. The detectable Δ depends only on scan-to-scan noise (`s_ss`,
`sigma`), N, and the spread of ages.

## Assumptions

- **Complete data**: every subject has all 5 regions at all 3 sites. If data will be
  missing, this shortcut isn't exact, and you would need a joint mixed model across
  sites instead.
- **Very small N**: with 5 subjects each test has only 3 degrees of freedom, so the
  result depends heavily on how spread out the subjects' ages are.
- **Scale of Δ**: scanner differences that multiply values also scale the slopes.
  Decide in advance whether Δ applies to raw, harmonized, or log-transformed values.
