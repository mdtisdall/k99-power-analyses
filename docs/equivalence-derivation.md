# Derivation: cross-site equivalence of age effects

This is the derivation behind [`equivalence.py`](../equivalence.py). See the
[README](../README.md#analysis-1-cross-site-equivalence-of-age-effects) for how to
run it and which parameters to set.

**The question.** The same subjects are scanned at 3 sites. Each site's data, on
their own, would be analyzed with

```
y ~ 0 + region + region:age_c + (1 | subject)
```

How small a between-site difference in regional age slopes (the margin Δ) can a
study rule out?

**The test.** The per-site mixed models are not used for testing. They only define
the slopes being compared (Step 1). The test is 15 ordinary regressions, one per
region and pair of sites, of each subject's between-site difference on age
(Step 2), each checked with an equivalence test (Steps 4–5).

## Notation and data-generating model

Subjects are $i = 1, \dots, N$, regions $r = 1, \dots, 5$, and sites
$s \in \lbrace 1, 2, 3 \rbrace$. Every subject is scanned once at every site, at
effectively the same age $a_i$. Let $x_i = a_i - \bar a$ be centered age, and
$S_{xx} = \sum_i x_i^2$.

The outcome for subject $i$, region $r$, site $s$ is

$$
y_{irs} = \mu_{rs} + \beta_{rs}\thinspace x_i + u_i + v_{ir} + w_{is} + \varepsilon_{irs},
$$

where all random terms are independent with mean 0:

| Term | Meaning | Variance | In the script |
|---|---|---|---|
| $\mu_{rs}$ | region mean at site $s$, including any constant site offset | fixed | not needed |
| $\beta_{rs} = \beta_r + \delta_{rs}$ | age slope: common slope plus site deviation | fixed | $\delta_{rs}$ = `true_dev` × Δ |
| $u_i$ | subject effect | $\tau^2$ | not needed |
| $v_{ir}$ | stable subject-specific regional deviation | $\sigma_v^2$ | not needed |
| $w_{is}$ | whole-scan offset, shared by all regions in one scan | $\sigma_w^2$ | `sd_scan` $= \sigma_w$ |
| $\varepsilon_{irs}$ | region-level measurement noise | $\sigma^2$ | `sd_noise` $= \sigma$ |

$w_{is}$ and $\varepsilon_{irs}$ are assumed normal, with the same variance at every
site and region. Write $\lambda^2 = \sigma_w^2 + \sigma^2$ for the total scan-level
variance of one region's value, and $\rho = \sigma_w^2 / \lambda^2$ for the share of
it that is shared across regions.

## Step 1: Each site's mixed model gives per-region least-squares slopes

At a single site, collect the five regional values of subject $i$ into
$\mathbf y_i \in \mathbb R^5$, and let $\mathbf z_i = (1, x_i)^\top$. The model
`y ~ 0 + region + region:age_c + (1 | subject)` has design matrix
$X_i = I_5 \otimes \mathbf z_i^\top$ for the coefficients
$\boldsymbol\theta = (\alpha_1, \beta_1, \dots, \alpha_5, \beta_5)^\top$, and some
within-subject covariance $V$ ($5 \times 5$, the same for every subject). Given the
estimated variance components, a mixed model's fixed effects are the generalized
least-squares (GLS) estimate

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

$V$ cancels, so each region's mixed-model slope equals an ordinary least-squares
(OLS) fit of that region alone on age, whatever the variance estimates. This is the
classical result for equations that share the same regressors (Zellner, 1962;
Kruskal, 1968). Because $x$ is centered, $M = \mathrm{diag}(N, S_{xx})$ and

$$
\hat\beta_{rs} = \frac{1}{S_{xx}} \sum_i x_i\thinspace y_{irs}.
$$

This step needs complete data: every subject must have all 5 regions.

## Step 2: The slope difference is a difference-score regression

For sites $j$ and $k$, define each subject's difference score
$D_{ir}^{jk} = y_{irj} - y_{irk}$. Because the same subjects, with the same $x_i$,
appear at both sites,

$$
\hat\gamma_r^{jk} \equiv \hat\beta_{rj} - \hat\beta_{rk}
= \frac{1}{S_{xx}} \sum_i x_i\thinspace D_{ir}^{jk},
$$

which is exactly the OLS slope from regressing $D_{ir}^{jk}$ on $x_i$, with an
intercept. This is an algebraic identity, not an approximation. Its standard error
correctly accounts for the same subjects being scanned at both sites; comparing the
standard errors of the separate per-site fits would ignore that pairing.

## Step 3: Distribution of each difference-score regression

Substituting the model into $D_{ir}^{jk}$:

$$
D_{ir}^{jk} = (\mu_{rj} - \mu_{rk}) + \gamma_r^{jk}\thinspace x_i +
(w_{ij} - w_{ik}) + (\varepsilon_{irj} - \varepsilon_{irk}),
\qquad \gamma_r^{jk} = \delta_{rj} - \delta_{rk}.
$$

The subject effect $u_i$, the regional deviation $v_{ir}$, and the common slope
$\beta_r$ cancel exactly. Constant site offsets go into the intercept. The errors are
independent across subjects and distributed $N(0, 2\lambda^2)$. So, conditional on
the ages, this is a standard normal simple linear regression:

$$
\hat\gamma_r^{jk} \sim N\Big(\gamma_r^{jk},\ \frac{2\lambda^2}{S_{xx}}\Big),
\qquad
\widehat{\mathrm{SE}}^2 = \frac{\mathrm{RSS}}{(N-2)\thinspace S_{xx}},
\qquad
T = \frac{\hat\gamma_r^{jk} - \gamma_r^{jk}}{\widehat{\mathrm{SE}}} \sim t_{N-2}.
$$

So power depends on $\lambda$, $N$, the ages (through $S_{xx}$), the true site
differences $\gamma$, and, slightly, on $\rho$ (through the correlation between
comparisons; Step 6). Neither $\tau^2$, $\sigma_v^2$, the regional means, nor the
true age slopes appear. Each comparison here uses its own standard error with
$N - 2$ degrees of freedom; see [the pooled-variance option](#pooled-variance-option)
for an alternative.

## Step 4: An equivalence test for one comparison

For a margin $\Delta > 0$, test
$H_0: \lvert\gamma\rvert \ge \Delta$ against $H_1: \lvert\gamma\rvert < \Delta$
with two one-sided tests (TOST; Schuirmann, 1987), each at level $\alpha$. With
$t^\ast = t_{1-\alpha,\thinspace N-2}$, both reject when

$$
\hat\gamma - t^\ast\thinspace \widehat{\mathrm{SE}} > -\Delta
\quad\text{and}\quad
\hat\gamma + t^\ast\thinspace \widehat{\mathrm{SE}} < \Delta .
$$

Equivalently, the $100(1-2\alpha)$% CI (90% for $\alpha = 0.05$) lies inside
$(-\Delta, \Delta)$. The test has level $\alpha$. It is conservative (rejects less
often than $\alpha$ when $\lvert\gamma\rvert = \Delta$) when the standard error is
large relative to $\Delta$, e.g. a rate of 0.026 instead of 0.05 at $N = 5$.

## Step 5: Combining the 15 comparisons

The claim "the age effects are the same at all three sites" means that
$\lvert\gamma_r^{jk}\rvert < \Delta$ for all 5 regions and all 3 site pairs. All 15
tests must pass, so no multiplicity correction is needed (Berger, 1982); the price
is lower power. Power is the probability that all 15 pass.

## Step 6: Why simulation is needed

The 15 comparisons are not independent. Conditional on the ages, the 15
estimates $\hat\gamma_r^{jk}$ are jointly normal, with these covariances (each
multiplied by $S_{xx}$):

| Pair of comparisons | $S_{xx}\thinspace \mathrm{Cov}$ |
|---|---|
| same region, same site pair (variance) | $2\lambda^2$ |
| different regions, same site pair | $2\sigma_w^2$ |
| same region, pairs (1,2) and (1,3), or (1,3) and (2,3) | $+\lambda^2$ |
| same region, pairs (1,2) and (2,3) | $-\lambda^2$ |
| different regions, pairs sharing a site | $\pm\sigma_w^2$ (same signs as above) |

The three site pairs are also linearly dependent:
$\hat\gamma_r^{23} = \hat\gamma_r^{13} - \hat\gamma_r^{12}$. The 15 standard errors
are random and correlated, and all 15 share the same $S_{xx}$. The probability that
all 15 pass therefore has no convenient closed form, so
[`equivalence.py`](../equivalence.py) estimates it by Monte Carlo. In each simulated
study it:

1. draws $N$ ages uniformly on `[age_min, age_max]` and centers them;
2. draws $w_{is}$ and $\varepsilon_{irs}$;
3. forms the 15 difference scores and computes each slope and
   $\widehat{\mathrm{SE}}$ in closed form (Step 3).

It leaves out $\mu$, $\beta_r$, $u_i$, and $v_{ir}$, because they cancel exactly
(Step 3). Ages are redrawn in every simulated study, so the result averages over
possible samples of ages.

## Step 7: Power and the smallest detectable Δ

**Planning scenario.** Power must be computed for some assumed true difference
between sites. Assuming the sites are exactly equal ($\gamma = 0$) is the best case,
and it is rarely true. The script instead assumes site deviations that are a
fixed fraction of the margin, $\delta_{rs} = f_{rs}\thinspace\Delta$ with
$f$ = `true_dev`. The default is $f_{r3} = 0.25$ for every region: site 3's age
slope differs from the other two by $0.25\Delta$. So $\gamma_r^{jk} = g_r^{jk}\Delta$
with $g = 0.25$ for pairs (1,3) and (2,3) and $g = 0$ for pair (1,2). This matters: at
$N = 20$, power at $\Delta = 0.015$ is 0.90 with no true difference, 0.70 with
$0.25\Delta$, and 0.15 with $0.5\Delta$.

**The smallest margin each study passes.** Write each slope estimate as
$\hat\gamma = \gamma + b$, where $b$ is its noise part. The true difference shifts
the estimate by exactly $\gamma$ and does not change $\widehat{\mathrm{SE}}$ (the
residuals are unaffected). So with $\gamma = g\Delta$ and $\lvert g \rvert < 1$, the
two conditions in Step 4 become

$$
(1 - g)\thinspace\Delta > b + t^\ast\thinspace\widehat{\mathrm{SE}}
\quad\text{and}\quad
(1 + g)\thinspace\Delta > t^\ast\thinspace\widehat{\mathrm{SE}} - b .
$$

Both are of the form "Δ above a threshold", so each simulated study has a smallest
margin it passes:

$$
\Delta_{\min} = \max_{r,\thinspace (j,k)} \max\Big(
\frac{b + t^\ast\thinspace\widehat{\mathrm{SE}}}{1 - g},\quad
\frac{t^\ast\thinspace\widehat{\mathrm{SE}} - b}{1 + g}\Big).
$$

With $g = 0$ this is $\lvert b\rvert + t^\ast\thinspace\widehat{\mathrm{SE}}$. The
study shows equivalence at margin $\Delta$ exactly when $\Delta_{\min} < \Delta$, so

$$
\text{power}(\Delta) = F(\Delta), \qquad \Delta^\ast_p = F^{-1}(p),
$$

where $F$ is the distribution function of $\Delta_{\min}$. The script estimates
both from the same `n_sims` simulated studies: power at the chosen `Delta` is the
fraction of studies with $\Delta_{\min} <$ `Delta`, and the smallest margin with
power $p$ is the empirical $p$-quantile of $\Delta_{\min}$. No search over $\Delta$
is needed.

## Step 8: How the answer scales

With true differences set as fractions of Δ and ages uniform over a range of width
$R$ = `age_max − age_min`, the distribution of every $\Delta_{\min}$ is
proportional to $\lambda / R$. Hence

$$
\Delta^\ast_p = \frac{\sqrt{\sigma_w^2 + \sigma^2}}{R}\thinspace g_p(N, \rho, f),
$$

where $g_p$ depends only on $N$, $\rho$, the planning scenario $f$, and $p$. By
simulation, for $N = 5$ and $p = 0.8$:

| $\rho$ | 0 | 0.2 | 0.5 | 0.8 |
|---|---|---|---|---|
| $g$, no true difference | 17.2 | 17.1 | 16.7 | 15.8 |
| $g$, default scenario ($0.25\Delta$) | 20.2 | 19.9 | 19.4 | 18.1 |

So the answer depends mostly on the **total** scan-level noise relative to the age
range, and little on how that noise splits between whole-scan offsets and
region-level noise. With the placeholder values ($\lambda = 0.112$, $R = 40$,
$\rho = 0.2$), this gives $\Delta^\ast_{0.8} \approx 0.112 / 40 \times 19.9 \approx 0.056$.

Why $g$ is so large at $N = 5$: for a single comparison, $S_{xx}$ is about
$(N-1)R^2/12$, so $\mathrm{SE} \approx \sqrt{24/(N-1)}\thinspace \lambda/R \approx 2.4\thinspace \lambda/R$.
With only 3 degrees of freedom, $t^\ast = 2.35$, so even a single comparison with
$\hat\gamma = 0$ needs $\Delta > 5.8\thinspace \lambda/R$. Taking the worst of 15
comparisons, with noisy standard errors and a randomly varying age spread, roughly
triples that.

## Pooled-variance option

Steps 3–4 give each of the 15 comparisons its own standard error, estimated from
that comparison's residuals with $N - 2$ degrees of freedom. If the noise really is
the same at every site and region (as the power calculation assumes), the analysis
can instead pool the residual variance across comparisons, for example by fitting
one mixed model to all three sites,

```
y ~ 0 + region:site + region:site:age_c + (1 | subject) + (1 | subject:region) + (1 | subject:site)
```

and testing the 15 slope-difference contrasts with Satterthwaite or Kenward–Roger
degrees of freedom.

**Why it helps.** With few subjects, separate standard errors cost power in two
ways. First, the critical value is large: $t^\ast = 2.35$ with 3 degrees of freedom
at $N = 5$, against about 1.70 with a pooled estimate. Second, each separate
standard error is itself very noisy, and because all 15 comparisons must pass, the
one that happens to be estimated largest decides the outcome. A pooled estimate is
close to the true standard error and removes both costs. In the default scenario
(smallest Δ at 80% power, with the pooled version approximated by averaging the
residual variances of the 10 independent comparisons):

| N | Separate SEs (script) | Pooled SE | Known variance (best possible) |
|---|---|---|---|
| 5 | 0.056 | 0.041 | 0.040 |
| 10 | 0.026 | 0.024 | 0.024 |
| 20 | 0.016 | 0.016 | 0.016 |

So pooling matters only for very small studies; from about 20 subjects it makes no
practical difference. The script reports the separate-SE version, which is
conservative.

**When not to pool.** Pooling is valid only if the noise is truly equal across
sites and regions, or the model allows for the differences. If one site or region
is noisier, a pooled standard error understates its uncertainty and the test
becomes too liberal for it. The comparisons are also correlated (Step 6), so a
pooled variance has fewer effective degrees of freedom than a simple count
suggests; the joint mixed model accounts for this, a simple average does not.

## Assumptions

- **Complete data.** Every subject has all 5 regions at all 3 sites. Step 1 relies
  on this. If data will be missing, the per-site slopes are no longer per-region
  OLS, and you would need a joint mixed model across sites instead.
- **Same age at every site.** Each subject's three scans are close enough in time
  that age is the same. Step 2 relies on this. There are also no scan-order or
  practice effects.
- **Noise.** Whole-scan offsets and region-level noise are normal, independent
  across subjects, and have the same SD at every site and region. Region-level
  noise is independent across regions, apart from the shared whole-scan offset.
  Normality makes the t distribution in Step 3 exact. Equal noise is optimistic:
  scanners differ, and the noisiest site or region will decide the result. Running
  the script with the largest plausible noise values gives a safer answer.
- **One margin for all regions.** The same absolute Δ applies to every region. That
  makes sense only if the regions are on a common scale with similar age slopes;
  otherwise Δ should be set per region.
- **Planning scenario.** Power is computed for the true site differences in
  `true_dev` (default: one site off by $0.25\Delta$ in every region). Larger true
  differences reduce power sharply (Step 7).
- **Site differences are additive.** Constant site offsets cancel. A multiplicative
  scanner bias scales the slopes ($\beta_{rs} = c_s\thinspace \beta_r$), which shows
  up as a true site difference $\gamma_r^{jk} = (c_j - c_k)\thinspace \beta_r$.
  Include it in `true_dev`, or decide in advance whether Δ applies to raw,
  harmonized, or log-transformed values.
- **Random ages.** Ages are redrawn in every simulated study. With few subjects, the
  result depends heavily on how spread out their ages are. If the actual subjects
  are already known, their real ages would give a more relevant answer.

## Numerical checks

[`checks/check_equivalence.py`](../checks/check_equivalence.py) verifies the
closed-form $\Delta_{\min}$ with true differences against brute-force simulation
(power 0.578 vs 0.580 at $N = 10$; 0.689 vs 0.679 at $N = 20$, Monte Carlo SE
0.009), that per-site mixed-model slopes equal per-region OLS slopes (to
$10^{-16}$), and the scaling in Step 8 (a factor of exactly 6 for 3 times the noise
and half the age range). The covariance table in Step 6 and the $t_{N-2}$
distribution in Step 3 were also checked by simulation during development.

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
