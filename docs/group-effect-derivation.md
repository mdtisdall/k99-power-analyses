# Derivation: detectable group-by-region effect

This is the derivation behind [`group_effect.py`](../group_effect.py). See the
[README](../README.md#analysis-2-detectable-group-by-region-effect) for how to run it
and which parameters to set.

## Design and model

There are $S$ sites with $n$ subjects each, so $N = Sn$ subjects in total. Each
subject is scanned **once, at one site**, and the outcome is measured in $R$
regions. Subjects belong to one of $K$ groups; group 0 is the reference group.
Within each site, subjects are split as evenly as possible across groups. For
subject $i$, let $s(i)$ be its site, $G_{ik}$ the indicator that it is in group $k$
($k = 1, \dots, K-1$), and $x_i$ its centered age.

The analysis model is

```
y ~ 0 + region + region:age_c + site + region:group + (1 | subject)
```

that is,

$$
y_{ir} = \alpha_r + \beta_r x_i + \sum_{k=1}^{K-1} \gamma_{rk} G_{ik} + \eta_{s(i)} + u_i + e_{ir},
\qquad
u_i \sim N(0, \tau^2), \quad e_{ir} \sim N(0, \sigma^2),
$$

with region intercepts $\alpha_r$, region-specific age slopes $\beta_r$, a site
effect $\eta_s$ shared by all regions ($\eta_1 = 0$), and the **group-by-region
effects** $\gamma_{rk}$: the difference between group $k$ and the reference group
in region $r$. The test of interest is $H_0: \gamma_{rk} = 0$ for each region and
each non-reference group, Bonferroni-corrected across the $m = R(K-1)$ tests.

The script's two noise parameters are `sd_total` $= \mathrm{SD}$ and `icc` $= \rho$,
where

$$
\mathrm{SD}^2 = \tau^2 + \sigma^2,
\qquad
\rho = \frac{\tau^2}{\tau^2 + \sigma^2}.
$$

$\mathrm{SD}$ is the SD of one region's value across subjects with the same group,
site, and age. $\rho$ is the correlation between two regions of the same subject.
Because each subject is scanned only once, whole-scan offsets are part of $u_i$, and
stable subject-specific regional deviations and measurement noise are part of
$e_{ir}$. So `(1 | subject)` describes the within-subject covariance correctly, as
long as all regions have the same SD and every pair of regions has the same
correlation. Both can be estimated by fitting this model to prior data:
$\mathrm{SD} = \sqrt{\hat\tau^2 + \hat\sigma^2}$ and
$\rho = \hat\tau^2 / (\hat\tau^2 + \hat\sigma^2)$.

## Step 1: The covariance splits into two strata

For subject $i$, collect the $R$ regional values into $\mathbf y_i$. Its covariance
is compound symmetric:

$$
V = \tau^2 J_R + \sigma^2 I_R ,
$$

where $J_R$ is the $R \times R$ matrix of ones. $V$ has two eigenspaces: the
direction $\mathbf 1_R$ (the subject's mean over regions), with eigenvalue
$R\tau^2 + \sigma^2$, and its orthogonal complement (deviations from that mean),
with eigenvalue $\sigma^2$. Accordingly, split each subject's data into

$$
\bar y_i = \frac{1}{R} \sum_r y_{ir}
\qquad\text{and}\qquad
d_{ir} = y_{ir} - \bar y_i .
$$

Substituting the model, with bars denoting averages over regions:

$$
\bar y_i = \bar\alpha + \bar\beta x_i + \sum_k \bar\gamma_k G_{ik} + \eta_{s(i)} + u_i + \bar e_i,
$$

$$
d_{ir} = (\alpha_r - \bar\alpha) + (\beta_r - \bar\beta) x_i + \sum_k (\gamma_{rk} - \bar\gamma_k) G_{ik} + (e_{ir} - \bar e_i).
$$

The site effects $\eta_s$ and the subject effects $u_i$ appear only in the subject
means. The errors of the two strata are uncorrelated (and, under normality,
independent), and the fixed effects split into two disjoint sets:
$(\bar\alpha, \bar\beta, \bar\gamma_k, \eta_s)$ for the means and the
region deviations $(\alpha_r - \bar\alpha, \beta_r - \bar\beta, \gamma_{rk} - \bar\gamma_k)$
for the deviations. GLS therefore separates into one fit per stratum. Within each
stratum the errors are homoscedastic, so GLS reduces to OLS:

- **Between subjects.** OLS of $\bar y_i$ on
  $Z_b = [\mathbf 1, x, G, \text{site indicators}]$, with error variance
  $\tau^2 + \sigma^2/R$ and $\mathrm{df}_b = N - (K + S)$ residual degrees of
  freedom.
- **Within subjects.** OLS of $d_{ir}$ on $Z_w = [\mathbf 1, x, G]$ for each
  region. Every region has the same regressors, so the joint GLS equals the
  per-region OLS fits (Zellner, 1962). Each fit has error variance
  $\sigma^2 (1 - 1/R)$, and the stratum has $\mathrm{df}_w = (R - 1)(N - K - 1)$
  residual degrees of freedom.

## Step 2: The standard error of a group-by-region effect

Because $\gamma_{rk} = \bar\gamma_k + (\gamma_{rk} - \bar\gamma_k)$, the estimate
combines one piece from each stratum, and the two pieces are independent:

$$
\mathrm{Var}(\hat\gamma_{rk}) = v_b + v_w,
\qquad
v_b = \Big(\tau^2 + \frac{\sigma^2}{R}\Big) c_b,
\qquad
v_w = \sigma^2 \Big(1 - \frac{1}{R}\Big) c_w,
$$

where $c_b$ and $c_w$ are the diagonal entries for group $k$ of
$(Z_b^\top Z_b)^{-1}$ and $(Z_w^\top Z_w)^{-1}$. These depend only on the design:
the group sizes, how groups are spread across sites, and the ages. The script
computes them exactly.

This closed form is exactly the GLS variance from the full mixed model with known
variance components. It was checked against the brute-force
$(\sum_i X_i^\top V^{-1} X_i)^{-1}$ to machine precision for several designs,
including 3 unequal groups and a single region.

**What this means in practice.** When groups are balanced within every site and
age is unrelated to group, the site indicators remove no group information, so
$c_b \approx c_w = c$ and

$$
\mathrm{Var}(\hat\gamma_{rk}) \approx (\tau^2 + \sigma^2)\thinspace c = \mathrm{SD}^2 c .
$$

For two equal groups, $c \approx 4/N$, so $\mathrm{SE} \approx 2\thinspace\mathrm{SD}/\sqrt N$.
That is the same as a two-sample t test on that one region: the correlation between
regions, the number of regions, and the site effects hardly matter. The random
subject effect does not reduce the error of a group difference because group varies
only between subjects, so each region's group comparison is subject to that region's
full between-subject variability. The correlation $\rho$ starts to matter when
groups are unbalanced within sites (the site adjustment then costs information,
$c_b > c_w$, and higher $\rho$ puts more weight on $c_b$). It would also matter for
a different question, whether the group difference *varies across regions*, which
uses only the within-subject stratum and its smaller variance $\sigma^2$.

## Step 3: The test and its degrees of freedom

Each effect is tested with the Wald statistic
$T = \hat\gamma_{rk} / \widehat{\mathrm{SE}}$, two-sided at level $\alpha / m$
(Bonferroni). $\widehat{\mathrm{SE}}^2$ combines variance estimates from both
strata, so its degrees of freedom follow from the Satterthwaite approximation:

$$
\mathrm{df} = \frac{(v_b + v_w)^2}{v_b^2 / \mathrm{df}_b + v_w^2 / \mathrm{df}_w}.
$$

When the true effect is $\gamma$, $T$ is approximately noncentral t with these df
and noncentrality $\gamma / \mathrm{SE}$, so

$$
\text{power}(\gamma) = P\big(T > t_c\big) + P\big(T < -t_c\big),
\qquad t_c = t_{1 - \alpha/(2m),\thinspace \mathrm{df}} .
$$

With the default design ($N = 200$), $\mathrm{df} \approx 550$, so this is
essentially a z test. The df matter only for small designs.

## Step 4: The smallest detectable effect

Power increases with $\gamma$, so the smallest detectable effect $\gamma^\ast$ is
the root of $\text{power}(\gamma) = p$, which the script finds numerically. Ages
are random, so $c_b$ and $c_w$ vary slightly between possible samples. The script
averages power over `n_designs` random age draws before solving. It reports
$\gamma^\ast$ in outcome units and as Cohen's $d = \gamma^\ast / \mathrm{SD}$. As
a check, the usual normal approximation gives

$$
\gamma^\ast \approx \big(z_{1-\alpha/(2m)} + z_p\big)\thinspace \mathrm{SE}.
$$

With the defaults ($m = 10$, $\alpha/m = 0.005$, $p = 0.8$, two equal groups of
100), this is $(2.81 + 0.84) \times 2/\sqrt{200} \approx 0.52$ SD, matching the
script's $d = 0.519$.

## Verification by simulation

The derivation was also checked against real mixed-model fits (statsmodels
`MixedLM`, REML). Data were simulated from the model with the defaults, with region
means, age slopes, and random site offsets added and an effect of $d = 0.519$ in one
region only. The fitted model above was then tested at $\alpha = 0.005$. In 1500 simulated
studies:

| Region | Rejection rate | Expected |
|---|---|---|
| with the effect (power) | 0.808 | 0.80 (Monte Carlo SE 0.010) |
| without an effect (type I error) | 0.0067 | 0.005 (Monte Carlo SE 0.0018) |

Both agree with the closed-form calculation within Monte Carlo error.

## Assumptions

- **One scan per subject, complete data.** Every subject has all $R$ regions.
  Step 1 relies on this. With missing regions the strata no longer separate, and
  power would need to be simulated with full mixed-model fits.
- **Compound-symmetric covariance.** All regions have the same SD, and every pair of
  regions has the same correlation. If some regions are noisier, the detectable
  effect for region $r$ scales with that region's SD. Run the script with each
  region's SD to see the range.
- **Additive site effects.** A site shifts all regions equally, as in the model.
  Region-specific scanner effects are not modeled. Balancing groups within sites
  protects the group comparison from them.
- **Balanced groups, comparable ages.** Groups are allocated evenly within each
  site, and age has the same distribution in every group. Groups that differ in age
  lose some information to the age adjustment.
- **Normal errors.** Needed for the t-based power calculation to be exact. Otherwise
  it is approximate, which is usually fine at these sample sizes.

## References

- Satterthwaite, F. E. (1946). An approximate distribution of estimates of variance
  components. *Biometrics Bulletin*, 2(6), 110–114.
- Zellner, A. (1962). An efficient method of estimating seemingly unrelated
  regressions and tests for aggregation bias. *Journal of the American Statistical
  Association*, 57(298), 348–368.
