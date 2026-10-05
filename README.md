# k99-power-analyses

Power analysis for showing that **region-specific age effects are equivalent across
3 MRI sites**, when the original single-site model is fit at each site separately.

## Design

- N subjects, each scanned at all 3 sites; no drop-out, no missing regions.
- Outcome measured in 5 brain regions per scan.
- Single-site analysis model (fit separately at each site), with `age_c` = centered age:

  ```r
  y ~ 0 + region + region:age_c + (1 | subject)
  ```

## Method

1. **"The same" is an equivalence claim.** Use TOST against a pre-specified margin Δ
   (outcome units per year). A non-significant site×age interaction does not show
   sameness.
2. **Difference-score regression.** With complete, balanced data (every subject has all
   5 regions and the same age at every site), the mixed model's region slopes equal
   per-region OLS slopes (GLS = OLS). So the difference between two sites' slopes for
   region *r* equals the slope from regressing the per-subject difference score
   `D_i = y[i, r, j] - y[i, r, k]` on `age_c`. That regression's SE is exact and
   accounts for the same subjects being scanned at both sites. Comparing SEs from the
   separate fits would ignore that pairing.
3. **Decision rule.** For each of 5 regions × 3 site pairs = 15 comparisons, the 90% CI
   (t, df = N − 2) of the D-on-age slope must lie inside (−Δ, +Δ). Sites are "the same"
   only if all 15 pass. This is an intersection-union test, so no α adjustment is
   needed, but **power = P(all 15 pass)**.
4. **What cancels.** Between-subject variance (τ²), stable subject-specific regional
   variance, region intercepts, true age slopes, and constant site offsets all drop out
   of the differences. Power depends only on scan-level variability, Δ, N, and
   Var(age). With equal noise at all sites:

   ```
   Var(slope diff) ≈ 2(σ²_subj:site + σ²) / (N · Var(age))
   ```

   Use this as a sanity check. Each test needs roughly |true diff| + t₀.₉₅ · SE < Δ,
   and all 15 must pass.

## Variables you need to supply

Set these in the **Parameters** block at the top of [`power.R`](power.R) or
[`power.py`](power.py); both use the same names. All the
current values are **placeholders**. Replace them with estimates from prior data,
ideally traveling-subject or test-retest scans.

| Variable | Meaning | Where it comes from |
|---|---|---|
| `Delta` | Equivalence margin Δ, in outcome units per year. **The most consequential choice.** | Justify before seeing the data, e.g. ±20% of the expected age slope. Decide whether it applies to the raw, harmonized, or log scale (see Caveats). |
| `s_ss` | σ_subj:site: SD of the scan-level global offset, shared across all regions within one scan. | Prior multi-site / test-retest data. |
| `sigma` | Length-3 vector: region-level measurement-noise SD at each site/scanner. | Prior data, per scanner. Alternatively, estimate the SD of between-site difference scores per region and pair directly. |
| `age_min`, `age_max` | Age distribution. The sim draws ages from uniform(`age_min`, `age_max`). Only Var(age) matters. | Expected recruitment range. If your age distribution isn't uniform, edit the `runif` line in `sim_once` (R) or the `rng.uniform` line in `power()` (Python). |
| `N_grid` | Candidate sample sizes to evaluate. | Feasible enrollment range. |
| `dslope` | Optional 5 × 3 matrix (regions × sites) of true site slope deviations from a common slope. Base case: all 0. | Set small non-zero values to get power when the sites truly differ slightly. |
| `alpha` | One-sided level of each TOST (0.05 gives a 90% CI). | Usually leave at 0.05. |
| `n_sims` | Monte Carlo replicates per N. | 2000 gives an MC SE of ≤ ~0.011 on power. |
| `seed` | RNG seed, for reproducibility. | |

You do **not** need τ², region intercepts, the true age slopes, or constant site
offsets, because they all cancel (see Method item 4).

## Running

There are two equivalent implementations. Pick whichever you prefer.

| | R: [`power.R`](power.R) | Python: [`power.py`](power.py) |
|---|---|---|
| Dependencies | base R only | numpy, scipy |
| How it fits | `lm()` per comparison, per simulation | closed-form OLS slope/SE, all simulations vectorized |
| Default run (10 N × 2000 sims) | ~15 s | ~1–2 s |
| Output | `power_curve.csv` | `power_curve_py.csv` |

Both print the power curve (`N`, `power` = P(all 15 TOSTs pass)) and write it to a
CSV, which is git-ignored. They agree to within Monte Carlo error.

With Nix, `shell.nix` provides both R and Python:

```bash
nix-shell --run "Rscript power.R"
```

```bash
nix-shell --run "python power.py"
```

Or without Nix, using a local R installation, or a Python environment with numpy
and scipy (`pip install numpy scipy`):

```bash
Rscript power.R
```

```bash
python power.py
```

To get power at a single N interactively:

```r
# R: source the functions, then
mean(replicate(2000, sim_once(N = 60, Delta = 0.003)))
```

```python
# Python
from power import power
power(N=60, Delta=0.003, n_sims=2000, rng=1)
```

## Caveats

- **Missing data.** If any data are missing, GLS = OLS no longer holds and the
  difference-score shortcut is no longer exact. Use a joint model instead, e.g.

  ```r
  y ~ 0 + region:site + region:site:age_c +
      (1 | subject) + (1 | subject:region) + (1 | subject:site)
  ```

  and get the slope-difference contrasts from `vcov()`. Alternatively, use a
  subject-level bootstrap that refits all three sites on the same resampled subjects.
  This needs `lme4` (and optionally `simr`). Under Nix, use `rWrapper` so the packages
  are visible to `Rscript`:

  ```bash
  nix-shell -p 'rWrapper.override { packages = with rPackages; [ lme4 simr ]; }' --run "Rscript your_script.R"
  ```

  On macOS (as of nixpkgs with R 4.4.2), `simr` fails to build because its
  `SparseM` dependency fails to compile with gfortran. `lme4` alone builds and loads fine:

  ```bash
  nix-shell -p 'rWrapper.override { packages = with rPackages; [ lme4 ]; }' --run "Rscript your_script.R"
  ```

- **Multiplicative scanner bias** scales slopes proportionally. Decide in advance
  whether Δ applies to the raw, harmonized, or log-transformed outcome.
