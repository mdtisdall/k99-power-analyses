# Derivation: comparing harmonizations by between-site disagreement in age slopes

This is the derivation behind [`harmonization.py`](../harmonization.py). See the
[README](../README.md#analysis-1-comparing-harmonizations-by-between-site-age-slope-disagreement)
for how to run it and which parameters to set.

**The question.** The same subjects are scanned at 3 sites, twice at each. Two
harmonization approaches, A and B, are applied to the same scans. After each
approach, each site's data, on their own, would be analyzed with

```
y ~ 0 + region + region:age_c + (1 | subject)
```

Does A leave less between-site disagreement in each region's age slope than B?

**The test.** The per-site mixed models only define the slopes being compared
(Step 1). Each between-site slope difference is a regression of each subject's
between-site difference on age (Step 2). An approach's disagreement in a region is
the sum of its squared slope differences (Step 3). The test estimates B's
disagreement minus A's, removes the bias that noise adds to squared estimates
(Step 4), and is one-sided (Step 5).

## Notation and data-generating model

Subjects are $i = 1, \dots, N$, regions $r = 1, \dots, 5$, sites
$s \in \lbrace 1, 2, 3 \rbrace$, scans $t = 1, \dots, n$ at each site, and approaches
$m \in \lbrace A, B \rbrace$. Every subject is scanned at every site at effectively
the same age $a_i$. Let $x_i = a_i - \bar a$ be centered age, and
$S_{xx} = \sum_i x_i^2$.

The subject's true value in region $r$ is
$\chi_{ir} = \mu_r + \beta_r\thinspace x_i + u_{ir}$, where $u_{ir}$ (SD
`sd_between`) is how the subject differs from others of the same age. After
approach $m$, scan $t$ at site $s$ gives

$$
y^{m}_{irst} = o^{m}_{rs} + c^{m}_{rs}\thinspace \chi_{ir} + e^{m}_{irst},
\qquad
e^{m}_{irst} = \eta^{m}_{irs} + w^{m}_{ist} + \varepsilon^{m}_{irst},
$$

| Term | Meaning | Variance | In the script |
|---|---|---|---|
| $o^{m}_{rs}$ | constant offset left by approach $m$ at site $s$ | fixed | not needed |
| $c^{m}_{rs}$ | gain (scale) left by approach $m$ at site $s$; 1 = harmonized | fixed | `gain_A`, `gain_B` |
| $\eta^{m}_{irs}$ | subject-by-site deviation, the same in every scan at that site | $\sigma_\eta^2$ | `sd_site` |
| $w^{m}_{ist}$ | whole-scan offset, shared by all regions in one scan | $\sigma_w^2$ | `sd_scan` |
| $\varepsilon^{m}_{irst}$ | region-level scan-to-scan noise | $\sigma^2$ | `sd_noise` |

The noise terms are normal, independent across subjects, and have the same
variance at every site and region. A and B are computed from the same scans, so
each of A's noise terms has correlation $\kappa$ (`noise_corr`) with the
corresponding term of B.

Each site's value is the average of its $n$ scans,
$\bar y^{m}_{irs} = o^{m}_{rs} + c^{m}_{rs}\thinspace\chi_{ir} + \bar e^{m}_{irs}$,
where $\bar e$ has variance $\sigma_\eta^2 + (\sigma_w^2 + \sigma^2)/n$. Repeats
reduce the scan-to-scan terms, not $\eta$.

## Step 1: Each site's mixed model gives per-region least-squares slopes

For one approach at one site, collect the five regional values of subject $i$ into
$\mathbf y_i \in \mathbb R^5$, and let $\mathbf z_i = (1, x_i)^\top$. The model has
design matrix $X_i = I_5 \otimes \mathbf z_i^\top$ and some within-subject
covariance $V$, the same for every subject. Its fixed effects are the generalized
least-squares estimate

$$
\hat{\boldsymbol\theta} = \Big(\sum_i X_i^\top V^{-1} X_i\Big)^{-1} \sum_i X_i^\top V^{-1} \mathbf y_i
= (I_5 \otimes M^{-1}) \sum_i \mathbf y_i \otimes \mathbf z_i ,
\qquad M = \sum_i \mathbf z_i \mathbf z_i^\top,
$$

using the Kronecker mixed-product rule. $V$ cancels, so each region's slope equals
an ordinary least-squares (OLS) fit of that region alone on age, whatever the
variance estimates (Zellner, 1962; Kruskal, 1968):
$\hat\beta^{m}_{rs} = S_{xx}^{-1} \sum_i x_i\thinspace \bar y^{m}_{irs}$. This needs
complete data: every subject has all 5 regions.

## Step 2: Each slope difference is a difference-score regression

For sites $j$ and $k$, let $D^{m,jk}_{ir} = \bar y^{m}_{irj} - \bar y^{m}_{irk}$.
Because the same subjects, with the same $x_i$, appear at both sites,

$$
\hat\gamma^{m,jk}_r \equiv \hat\beta^{m}_{rj} - \hat\beta^{m}_{rk}
= \frac{1}{S_{xx}} \sum_i x_i\thinspace D^{m,jk}_{ir},
$$

the OLS slope of $D^{m,jk}_{ir}$ on $x_i$. Substituting the model,

$$
D^{m,jk}_{ir} = \text{const} + \gamma^{m,jk}_r\thinspace x_i
+ (c^{m}_{rj} - c^{m}_{rk})\thinspace u_{ir} + \bar e^{m}_{irj} - \bar e^{m}_{irk},
\qquad \gamma^{m,jk}_r = (c^{m}_{rj} - c^{m}_{rk})\thinspace \beta_r .
$$

So the true between-site slope difference comes only from the gains: constant
offsets go into the intercept. A gain mismatch also lets between-subject
variation $u_{ir}$ leak into the difference scores, which adds noise of SD
$\lvert c_{rj} - c_{rk}\rvert\thinspace$`sd_between`. The script includes this.

## Step 3: Disagreement is a quadratic form

An approach's disagreement in region $r$ is

$$
\theta^{m}_r = \sum_{j<k} \big(\gamma^{m,jk}_r\big)^2
= 3 \sum_s \big(\beta^{m}_{rs} - \bar\beta^{m}_{r}\big)^2 ,
$$

three times the variance of that region's slope across sites. Only two of the
three differences are free ($\gamma^{23} = \gamma^{13} - \gamma^{12}$), so with
$\mathbf g^{m}_r = (\gamma^{m,12}_r, \gamma^{m,13}_r)^\top$,

$$
\theta^{m}_r = 2\big(\gamma^{m,12}_r\big)^2 + 2\big(\gamma^{m,13}_r\big)^2 -
2\thinspace\gamma^{m,12}_r\gamma^{m,13}_r
= \mathbf g^{m\top}_r K\thinspace \mathbf g^{m}_r ,
$$

where $K$ is the $2 \times 2$ matrix with 2 on the diagonal and −1 off it. Stack
$\mathbf z_r = (\mathbf g^{A}_r, \mathbf g^{B}_r)$, and let $Q$ be the
block-diagonal matrix with blocks $-K$ (for A) and $K$ (for B). The quantity
tested in region $r$ is

$$
\Theta_r = \theta^{B}_r - \theta^{A}_r = \mathbf z_r^\top Q\thinspace \mathbf z_r ,
$$

and the pooled quantity is $\Theta = \sum_r \Theta_r$, the same form for the
20-vector of all four slope differences in all five regions, with $Q$ repeated
along the diagonal.

## Step 4: Removing the noise bias, and the standard error

The four (or twenty) estimated slope differences are each an OLS slope of a
difference score on the same $x$, so, conditional on the ages, they are jointly
normal with mean $\mathbf z$ and covariance $\Sigma / S_{xx}$, where $\Sigma$ is
the covariance of the difference-score errors. It is estimated from the residuals
of the difference-score regressions,
$\hat\Sigma = \sum_i \hat{\mathbf r}_i \hat{\mathbf r}_i^\top / (N-2)$, which is
unbiased.

For a normal vector, $E[\hat{\mathbf z}^\top Q \hat{\mathbf z}] = \mathbf z^\top Q
\mathbf z + \mathrm{tr}(Q\Sigma)/S_{xx}$. The second term is the noise in the
squared estimates. It differs between the approaches whenever their noise
differs, so leaving it in would reward a less noisy approach even if it removed no
scanner bias. The estimate is

$$
\hat\Theta = \hat{\mathbf z}^\top Q\thinspace \hat{\mathbf z} - \mathrm{tr}(Q\hat\Sigma)/S_{xx},
$$

which is unbiased. Its variance is
$2\thinspace\mathrm{tr}\big((Q\Sigma_z)^2\big) + 4\thinspace\mathbf z^\top Q\Sigma_z Q\thinspace\mathbf z$
with $\Sigma_z = \Sigma/S_{xx}$. The script plugs in $\hat\Sigma$ and
$\hat{\mathbf z}$. Plugging in $\hat{\mathbf z}$ overstates the second term on
average (by $4\thinspace\mathrm{tr}((Q\Sigma_z)^2)$), so the test is conservative.
Removing that excess made the test too liberal when neither approach had
disagreement (rejection rates up to 0.19 at nominal 0.05), so the plug-in version
is used.

## Step 5: The tests

- **Primary (pooled).** Reject $H_0: \Theta \le 0$ in favor of $H_1: \Theta > 0$
  (A disagrees less than B, summed over regions) when
  $\hat\Theta / \widehat{\mathrm{SE}} > z_{1-\alpha}$.
- **Secondary (each region).** The same test for each $\Theta_r$, with Holm's
  correction across the 5 regions, says in which regions A helps. Requiring all 5
  regions to pass (an intersection-union test) is a much stronger claim and needs
  several times as many subjects.

Pooling assumes the scanner bias that B leaves is spread over regions. If B's
mismatch is confined to one region, the pooled test still detects it but is
diluted by the noise of the other four.

**Why not simpler tests.** A model-by-site-by-age interaction tests whether the
slope differences *change* between approaches, not whether they *shrink*: A could
reverse them or move them to another region. Comparing a site-by-age test from A
with one from B compares p-values, and ignores that A and B come from the same
scans. Comparing squared slope differences without the correction in Step 4
confounds noise with disagreement.

## Step 6: Simulation

The tests use estimated standard errors and, for the per-region tests, Holm's
procedure, so power is estimated by Monte Carlo. In each simulated study,
[`harmonization.py`](../harmonization.py):

1. draws $N$ ages uniformly on `[age_min, age_max]` and centers them, and draws
   $u_{ir}$;
2. draws B's averaged noise terms directly (sums of normals are normal), and A's
   with correlation $\kappa$ to B's;
3. forms the 20 difference scores, their slopes and residual covariance in closed
   form, and the pooled and per-region statistics.

Offsets and regional means are left out because they cancel. The **smallest
detectable mismatch** scales B's gain deviations from A,
$c^{B} = c^{A} + k\thinspace(c^{B}_{\text{scenario}} - c^{A})$, finds the smallest
$k$ with the target power by bisection (with the same random draws for every
$k$), and reports B's largest between-site slope difference at that $k$.

## Assumptions

- **Complete data**, the **same age at every site**, and no scan-order effects, as
  for Steps 1–2.
- **Scanner bias acts as a gain and an offset.** Offsets cancel. Bias that depends
  on age in some other way would add to $\gamma$ in the same way as a gain.
- **Noise.** Normal, the same SD at every site, region and approach, and
  independent across regions apart from the whole-scan offset. The between-subject
  terms $u_{ir}$ are independent across regions. A and B have the same noise
  level; if one approach adds noise, its slope differences become less precise,
  but the bias correction keeps the comparison fair.
- **Harmonization parameters estimated out of sample.** If A or B estimates its
  site parameters from the same traveling subjects, it fits their noise and its
  disagreement looks too small. Estimate the parameters for each subject from the
  other subjects (cross-fitting) or from separate calibration data. The remaining
  estimation error then acts as a small leftover gain mismatch for that approach.
- **Scale.** An approach that shrinks all values shrinks all slopes, and so its
  slope differences, without improving agreement. Check that the pooled slopes are
  similar under A and B; if not, compare relative disagreement
  ($\theta^{m}_r / \bar\beta^{m\thinspace 2}_r$), with a bootstrap over subjects
  for inference.
- **Large-sample test.** The statistic is compared with a normal distribution. For
  the real analysis, a bootstrap over subjects (resampling each subject's 3 sites,
  2 scans, 5 regions and both approaches together) is a good cross-check.

## Numerical checks

[`checks/check_harmonization.py`](../checks/check_harmonization.py) checks that:

1. the vectorized script matches a brute-force simulation that draws every scan
   (two per site), averages them, fits each site's slope separately, and confirms
   the Step 2 identity (pooled power 0.693 vs 0.677 at $N = 20$ in the default
   scenario, 0.873 vs 0.869 at $N = 40$ with A at ±5% and B at ±10% gains, Monte
   Carlo SE about 0.01), and that the estimate is unbiased (mean 1.180 vs a true
   1.200, Monte Carlo SE 0.012; 0.2226 vs 0.2250, SE 0.0018);
2. the tests keep their level at the null boundary. With equal disagreement in A
   and B, the rejection rate is at most 0.044 (pooled) and 0.039 (per region,
   before Holm) for $N$ from 10 to 100, and at most 0.005 when neither approach
   has any disagreement;
3. per-site mixed-model slopes equal per-region OLS slopes (to $10^{-16}$).

## References

- Holm, S. (1979). A simple sequentially rejective multiple test procedure.
  *Scandinavian Journal of Statistics*, 6(2), 65–70.
- Kruskal, W. (1968). When are Gauss–Markov and least squares estimators identical?
  A coordinate-free approach. *Annals of Mathematical Statistics*, 39(1), 70–75.
- Zellner, A. (1962). An efficient method of estimating seemingly unrelated
  regressions and tests for aggregation bias. *Journal of the American Statistical
  Association*, 57(298), 348–368.
