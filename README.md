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

| Parameter | Symbol | What it is |
|---|---|---|
| `N_grid` | $N$ | Number(s) of subjects to evaluate, e.g. `[5]` or `[5, 10, 20]`. |
| `target_power` | $p$ | Required probability that the study shows equivalence. Usually 0.80. |
| `s_ss` | $\sigma_w$ | SD of a whole-scan offset: how much a subject's values shift together, across all regions, from one scan to another. |
| `sigma` | $\sigma$ | Measurement-noise SD for one region in one scan. The same for all sites and regions. |
| `age_min`, `age_max` | $R$ = range | Age range of the sample. Ages are drawn uniformly from this range. |
| `dslope` | $\delta_{rs}$ | Optional true differences in age slope between sites (5 regions × 3 sites). Leave at 0 for the standard calculation. |
| `alpha` | $\alpha$ | Level of each one-sided test. Leave at 0.05 (gives 90% CIs). |
| `n_sims` | | Simulations per sample size. 10000 gives Δ to about ±1%. |
| `seed` | | Random seed, so results are reproducible. |

You don't need between-subject variance, regional means, the true age slopes, or
fixed site offsets, because they cancel out (see the [derivation](#derivation)).

## Output: Δ

Δ is the **equivalence margin**: a between-site difference in age slope, in outcome
units per year. With the given number of subjects, a study has `target_power` chance
of showing that every between-site slope difference lies within ±Δ. Smaller is
better. To judge whether it's useful, compare Δ with the expected age slope. A Δ
around 20% of the slope would be a convincing "same". A Δ as large as the slope
itself means the sites could differ by the whole age effect.

When you analyze real data, choose the margin before you see the data.

## What it calculates, in brief

Each site's data are analyzed separately with the model

```
y ~ 0 + region + region:age_c + (1 | subject)
```

To show that the age effects are the *same*, the script uses equivalence tests
(TOST) rather than a non-significant site × age interaction, which can't show
sameness. For each region (5) and each pair of sites (3), it regresses the
per-subject difference between the two sites on age. That slope is exactly the
difference between the two sites' age effects. The 90% CI for that slope must fall
inside (−Δ, +Δ). The sites count as "the same" only if **all 15** comparisons pass.
The script reports the smallest Δ for which that happens with probability
`target_power`.

## Derivation

### Notation and data-generating model

Subjects are $i = 1, \dots, N$, regions $r = 1, \dots, 5$, and sites
$s \in \lbrace 1, 2, 3 \rbrace$. Every subject is scanned once at every site, at
effectively the same age $a_i$. Let $x_i = a_i - \bar a$ be centered age, and
$S_{xx} = \sum_i x_i^2$.

The outcome for subject $i$, region $r$, site $s$ is

$$
y_{irs} = \mu_{rs} + \beta_{rs}\, x_i + u_i + v_{ir} + w_{is} + \varepsilon_{irs},
$$

where all random terms are independent with mean 0:

| Term | Meaning | Variance |
|---|---|---|
| $\mu_{rs}$ | region mean at site $s$, including any constant site offset | fixed |
| $\beta_{rs} = \beta_r + \delta_{rs}$ | age slope: common slope plus site deviation (`dslope`) | fixed |
| $u_i$ | subject effect | $\tau^2$ |
| $v_{ir}$ | stable subject-specific regional deviation | $\sigma_v^2$ |
| $w_{is}$ | whole-scan offset, shared by all regions in one scan (`s_ss`) | $\sigma_w^2$ |
| $\varepsilon_{irs}$ | region-level measurement noise (`sigma`) | $\sigma^2$ |

$w_{is}$ and $\varepsilon_{irs}$ are assumed normal. Write
$\lambda^2 = \sigma_w^2 + \sigma^2$ for the total scan-level variance of one
region's value.

### Step 1: Each site's mixed model gives per-region OLS slopes

At a single site, collect the five regional values of subject $i$ into
$\mathbf y_i \in \mathbb R^5$, and let $\mathbf z_i = (1, x_i)^\top$. The analysis
model `y ~ 0 + region + region:age_c + (1 | subject)` has design matrix
$X_i = I_5 \otimes \mathbf z_i^\top$ for the coefficients
$\boldsymbol\theta = (\alpha_1, \beta_1, \dots, \alpha_5, \beta_5)^\top$, and some
within-subject covariance $V$ ($5 \times 5$, the same for every subject). Given the
estimated variance components, a mixed model's fixed effects are the GLS estimate

$$
\hat{\boldsymbol\theta} = \Big(\sum_i X_i^\top V^{-1} X_i\Big)^{-1} \sum_i X_i^\top V^{-1} \mathbf y_i .
$$

Using the Kronecker mixed-product rule,

$$
X_i^\top V^{-1} X_i = V^{-1} \otimes \mathbf z_i \mathbf z_i^\top,
\qquad
X_i^\top V^{-1} \mathbf y_i = (V^{-1} \otimes I_2)(\mathbf y_i \otimes \mathbf z_i).
$$

With $M = \sum_i \mathbf z_i \mathbf z_i^\top$,

$$
\hat{\boldsymbol\theta}
= (V \otimes M^{-1})(V^{-1} \otimes I_2) \sum_i \mathbf y_i \otimes \mathbf z_i
= (I_5 \otimes M^{-1}) \sum_i \mathbf y_i \otimes \mathbf z_i .
$$

$V$ cancels. The block for region $r$ is $M^{-1} \sum_i \mathbf z_i\, y_{ir}$, which
is the OLS fit of that region alone on age. This is the classical result that GLS
equals OLS when every equation has the same regressors (Zellner, 1962;
Kruskal, 1968). It holds for **any** $V$, so it does not matter how well the
`(1 | subject)` term describes the true within-subject covariance, or what values
REML estimates. Because $x$ is centered, $M = \mathrm{diag}(N, S_{xx})$ and

$$
\hat\beta_{rs} = \frac{1}{S_{xx}} \sum_i x_i\, y_{irs}.
$$

This step needs complete data: every subject must have all 5 regions, so that every
region has the same regressors.

### Step 2: The slope difference is a difference-score regression

For sites $j$ and $k$, define the per-subject difference score
$D_{ir}^{jk} = y_{irj} - y_{irk}$. Because the same subjects, with the same $x_i$,
appear at both sites,

$$
\hat\gamma_r^{jk} \equiv \hat\beta_{rj} - \hat\beta_{rk}
= \frac{1}{S_{xx}} \sum_i x_i\, D_{ir}^{jk},
$$

which is exactly the OLS slope from regressing $D_{ir}^{jk}$ on $x_i$, with an
intercept. This is an algebraic identity, not an approximation. The point of using
it is that this regression's standard error correctly accounts for the same subjects
being scanned at both sites. Comparing the standard errors of the separate per-site
fits would ignore that pairing.

### Step 3: Distribution of each difference-score regression

Substituting the model into $D_{ir}^{jk}$:

$$
D_{ir}^{jk} = (\mu_{rj} - \mu_{rk}) + \gamma_r^{jk}\, x_i +
(w_{ij} - w_{ik}) + (\varepsilon_{irj} - \varepsilon_{irk}),
\qquad \gamma_r^{jk} = \delta_{rj} - \delta_{rk}.
$$

The subject effect $u_i$, the regional deviation $v_{ir}$, and the common slope
$\beta_r$ cancel exactly. Constant site offsets go into the intercept. The errors are
independent across subjects and distributed $N(0, 2\lambda^2)$. So, conditional on
the ages, this is a classical normal simple linear regression:

$$
\hat\gamma_r^{jk} \sim N\Big(\gamma_r^{jk},\ \frac{2\lambda^2}{S_{xx}}\Big),
\qquad
\widehat{\mathrm{SE}}^2 = \frac{\mathrm{RSS}}{(N-2)\, S_{xx}},
\qquad
T = \frac{\hat\gamma_r^{jk} - \gamma_r^{jk}}{\widehat{\mathrm{SE}}} \sim t_{N-2}.
$$

So the only quantities that affect power are $\lambda$ (through $\sigma_w$ and
$\sigma$), $N$, the ages (through $S_{xx}$), and any true site differences
$\gamma$. Neither $\tau^2$, $\sigma_v^2$, the regional means, nor the true age slopes
appear.

### Step 4: An equivalence test for one comparison

For a margin $\Delta > 0$, test
$H_0: \lvert\gamma\rvert \ge \Delta$ against $H_1: \lvert\gamma\rvert < \Delta$
with two one-sided tests (TOST; Schuirmann, 1987), each at level $\alpha$. With
$t^\ast = t_{1-\alpha,\,N-2}$, both reject when

$$
\hat\gamma - t^\ast\, \widehat{\mathrm{SE}} > -\Delta
\quad\text{and}\quad
\hat\gamma + t^\ast\, \widehat{\mathrm{SE}} < \Delta
\quad\Longleftrightarrow\quad
\lvert\hat\gamma\rvert + t^\ast\, \widehat{\mathrm{SE}} < \Delta .
$$

Equivalently, the $100(1-2\alpha)\%$ CI (90% for $\alpha = 0.05$) lies inside
$(-\Delta, \Delta)$. The test has size $\alpha$.

### Step 5: Combining the 15 comparisons

The claim "the age effects are the same at all three sites" means that
$\lvert\gamma_r^{jk}\rvert < \Delta$ for all 5 regions and all 3 site pairs. The
null hypothesis is the union of the 15 individual nulls, and we reject it only if
all 15 TOSTs reject. This is an intersection-union test, which has size at most
$\alpha$ with **no multiplicity adjustment** (Berger, 1982). The cost is in power,
which is the probability that all 15 pass:

$$
\text{power}(\Delta) = P\Big(\max_{r,\,(j,k)} \big[\lvert\hat\gamma_r^{jk}\rvert + t^\ast\, \widehat{\mathrm{SE}}_r^{jk}\big] < \Delta\Big).
$$

### Step 6: Why simulation is needed

The 15 comparisons are not independent. Conditional on the ages, the 15
estimates $\hat\gamma_r^{jk}$ are jointly normal, with these covariances (each
multiplied by $S_{xx}$):

| Pair of comparisons | $S_{xx}\,\mathrm{Cov}$ |
|---|---|
| same region, same site pair (variance) | $2\lambda^2$ |
| different regions, same site pair | $2\sigma_w^2$ |
| same region, pairs (1,2) and (1,3), or (1,3) and (2,3) | $+\lambda^2$ |
| same region, pairs (1,2) and (2,3) | $-\lambda^2$ |
| different regions, pairs sharing a site | $\pm\sigma_w^2$ (same signs as above) |

The three site pairs are also linearly dependent:
$\hat\gamma_r^{23} = \hat\gamma_r^{13} - \hat\gamma_r^{12}$. The 15 standard errors
are random and correlated, and all 15 share the same $S_{xx}$. The probability above
therefore has no convenient closed form, so `power.py` estimates it by Monte Carlo.
In each simulated study it:

1. draws $N$ ages uniformly on `[age_min, age_max]` and centers them;
2. draws $w_{is}$ and $\varepsilon_{irs}$, and adds $\delta_{rs}\, x_i$;
3. forms the 15 difference scores and computes each $\hat\gamma$ and
   $\widehat{\mathrm{SE}}$ in closed form (Step 3).

It leaves out $\mu$, $\beta_r$, $u_i$, and $v_{ir}$, because they cancel exactly
(Step 3), so including them would not change any result. Ages are redrawn in every
simulated study, so the result averages over possible samples of ages.

### Step 7: From power to the smallest detectable Δ

Define, for each simulated study, the smallest margin it would pass:

$$
\Delta_{\min} = \max_{r,\,(j,k)} \big[\lvert\hat\gamma_r^{jk}\rvert + t^\ast\, \widehat{\mathrm{SE}}_r^{jk}\big].
$$

The study shows equivalence at margin $\Delta$ exactly when $\Delta_{\min} < \Delta$,
so $\text{power}(\Delta) = F(\Delta)$, the distribution function of
$\Delta_{\min}$. This increases with $\Delta$. The smallest margin with power at
least $p$ is therefore its $p$-quantile,

$$
\Delta^\ast_p = F^{-1}(p),
$$

which `power.py` estimates as the empirical $p$-quantile of $\Delta_{\min}$ across
`n_sims` simulated studies. No search over $\Delta$ is needed.

### Step 8: How the answer scales

With no true site differences ($\delta = 0$) and ages uniform over a range of
width $R$ = `age_max − age_min`, every $\hat\gamma$ and $\widehat{\mathrm{SE}}$ is
proportional to $\lambda / R$. Hence

$$
\Delta^\ast_p = \frac{\sqrt{\sigma_w^2 + \sigma^2}}{R}\; g_p(N, \rho),
\qquad \rho = \frac{\sigma_w^2}{\sigma_w^2 + \sigma^2},
$$

where $g_p$ depends only on $N$, the split of the noise $\rho$, and $p$. By
simulation, for $N = 5$ and $p = 0.8$, $g \approx 17.0$, $17.0$, $16.6$, and $15.7$
for $\rho = 0$, $0.2$, $0.5$, and $0.8$. So the answer depends almost entirely on the
**total** scan-level noise relative to the age range, and hardly on how that noise
splits between whole-scan offsets and region-level noise. With the placeholder
values ($\lambda = 0.112$, $R = 40$, ages 25–65, $\rho = 0.2$), this gives
$\Delta^\ast_{0.8} \approx 0.112 / 40 \times 17.0 \approx 0.048$.

Why $g$ is so large at $N = 5$: for a single comparison, $S_{xx}$ is about
$(N-1)R^2/12$, so $\mathrm{SE} \approx \sqrt{24/(N-1)}\, \lambda/R \approx 2.4\,\lambda/R$.
With only 3 degrees of freedom, $t^\ast = 2.35$, so even a single comparison with
$\hat\gamma = 0$ needs $\Delta > 5.8\,\lambda/R$. Taking the worst of 15
comparisons, with noisy SEs and a randomly varying age spread, roughly triples
that.

## Assumptions

- **Complete data.** Every subject has all 5 regions at all 3 sites. Step 1 relies
  on this. If data will be missing, the per-site slopes are no longer per-region
  OLS, and you would need a joint mixed model across sites instead.
- **Same age at every site.** Each subject's three scans are close enough in time
  that age is the same. Step 2 relies on this.
- **Noise.** Whole-scan offsets and region-level noise are normal, independent
  across subjects, and have the same SD at every site and region. Normality makes
  the t distribution in Step 3 exact. Without it, the result is approximate.
- **Site differences are additive.** Constant site offsets cancel. A multiplicative
  scanner bias scales the slopes ($\beta_{rs} = c_s\, \beta_r$), which shows up as a
  true site difference $\gamma_r^{jk} = (c_j - c_k)\,\beta_r$. Model it with `dslope`,
  or decide in advance whether Δ applies to raw, harmonized, or log-transformed
  values.
- **Random ages.** Ages are redrawn in every simulated study. With 5 subjects, the
  result depends heavily on how spread out their ages are. If the actual subjects
  are already known, their real ages would give a more relevant answer.

## References

- Berger, R. L. (1982). Multiparameter hypothesis testing and acceptance sampling.
  *Technometrics*, 24(4), 295–300.
- Kruskal, W. (1968). When are Gauss–Markov and least squares estimators identical?
  A coordinate-free approach. *Annals of Mathematical Statistics*, 39(1), 70–75.
- Schuirmann, D. J. (1987). A comparison of the two one-sided tests procedure and
  the power approach for assessing the equivalence of average bioavailability.
  *Journal of Pharmacokinetics and Biopharmaceutics*, 15(6), 657–680.
- Zellner, A. (1962). An efficient method of estimating seemingly unrelated
  regressions and tests for aggregation bias. *Journal of the American Statistical
  Association*, 57(298), 348–368.
