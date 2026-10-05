# Derivation: detectable group-by-region effect

This is the derivation behind [`group_effect.py`](../group_effect.py). See the
[README](../README.md#analysis-2-detectable-group-by-region-effect) for how to run it
and which parameters to set.

## Design and model

There are $S$ sites with $n$ subjects each, so $N = Sn$ subjects in total. Each
subject is scanned **once, at one site**, and the outcome is measured in $R$
regions. Subjects belong to one of $K$ groups; group 0 is the reference group. For
subject $i$, let $s(i)$ be its site, $G_{ik}$ the indicator that it is in group $k$
($k = 1, \dots, K-1$), and $x_i$ its centered age.

### Balanced-sites simplification

**For the purposes of the power calculation, every site is assumed to recruit the
same number of subjects from each group** ($n/K$ per group per site). This is a
simplification for planning, not a requirement of the analysis: the model
adjusts for site whatever the allocation. Real recruitment is rarely exactly
balanced, so treat the result as the best case for a given total $N$ (see
Step 2 for what imbalance costs).

When $n$ is not divisible by $K$, the script splits each site as evenly as
possible, rotates the extra subjects between groups from site to site (so overall
group sizes stay as equal as possible), computes the SE exactly for that design,
and prints a note that the simplification holds only approximately.

### Analysis model

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

The power calculation matches a Wald t test with Satterthwaite degrees of freedom,
which is what lmerTest (R) reports for this model. Software that reports z-based
p-values (e.g. statsmodels `MixedLM`) is too liberal when the degrees of freedom
are small: at a nominal $\alpha = 0.005$, its true false-positive rate is 0.012
with 18 df and 0.007 with 56 df. With one region ($R = 1$) there is no random
effect to fit; use ordinary regression, for which the same formulas hold.

The script's two noise parameters are `sd_total` $= \mathrm{SD}$ and
`region_corr` $= \rho$, where

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

## Step 1: The data split into subject means and regional deviations

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
means. The errors of the two parts are uncorrelated (and, under normality,
independent), and the fixed effects split into two disjoint sets:
$(\bar\alpha, \bar\beta, \bar\gamma_k, \eta_s)$ for the means and the regional
deviations $(\alpha_r - \bar\alpha, \beta_r - \bar\beta, \gamma_{rk} - \bar\gamma_k)$
for the deviations. The generalized least-squares (GLS) fit of the mixed model
therefore separates into one fit for each part, and each reduces to ordinary least
squares (OLS):

- **Subject means.** OLS of $\bar y_i$ on
  $Z_b = [\mathbf 1, x, G, \text{site indicators}]$. There is one value per subject,
  with independent errors of variance $\tau^2 + \sigma^2/R$, so GLS is OLS. This fit
  has $\mathrm{df}_b = N - (K + S)$ residual degrees of freedom.
- **Regional deviations.** OLS of $d_{ir}$ on $Z_w = [\mathbf 1, x, G]$ for each
  region. The deviations of one subject are correlated with each other, but every
  region has the same regressors, so the joint GLS fit equals the separate OLS fits
  (Zellner, 1962). Each fit has error variance $\sigma^2 (1 - 1/R)$, and together
  they have $\mathrm{df}_w = (R - 1)(N - K - 1)$ residual degrees of freedom.

## Step 2: The standard error of a group-by-region effect

Because $\gamma_{rk} = \bar\gamma_k + (\gamma_{rk} - \bar\gamma_k)$, the estimate
combines one piece from each part, and the two pieces are independent:

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
variance components. It matches the brute-force $(\sum_i X_i^\top V^{-1} X_i)^{-1}$
to machine precision, including for unbalanced designs, 3 unequal groups, and a
single region ([checks](#numerical-checks)).

**What the balanced-sites simplification buys.** Under the simplification, group
is unrelated to site, so adjusting for site removes no group information. If age
is also unrelated to group, $c_b \approx c_w = c$ and

$$
\mathrm{Var}(\hat\gamma_{rk}) \approx (\tau^2 + \sigma^2)\thinspace c = \mathrm{SD}^2 c .
$$

With $K$ groups of $N/K$ subjects, $c \approx 2K/N$. For two groups that is $4/N$,
so $\mathrm{SE} \approx 2\thinspace\mathrm{SD}/\sqrt N$, the same as a two-sample t
test on that one region. For a fixed total $N$, the number of sites, the number of
regions (apart from the Bonferroni correction), and $\rho$ hardly matter. Group
differs only between subjects, so the random subject effect cannot remove
subject-to-subject noise from a group comparison.

**What imbalance costs.** When groups are unbalanced within sites, adjusting for
site costs information ($c_b > c_w$), and a higher $\rho$ puts more weight on
$c_b$. For example, with 20 sites of 3 subjects, split alternately 2:1 and 1:2 (30
per group overall), the SE is larger than the ideal $2\thinspace\mathrm{SD}/\sqrt N$
by about 2%, 4.5%, and 7% for $\rho = 0.1$, $0.5$, and $0.9$. The detectable effect
grows by the same proportion. A correlation between regions would also matter for a
different question, whether the group difference *varies across regions*, which
uses only the regional deviations and their smaller variance $\sigma^2$.

## Step 3: The test and its degrees of freedom

Each effect is tested with the Wald statistic
$T = \hat\gamma_{rk} / \widehat{\mathrm{SE}}$, two-sided at level $\alpha / m$
(Bonferroni). $\widehat{\mathrm{SE}}^2$ combines variance estimates from both
parts, so its degrees of freedom come from the Satterthwaite approximation:

$$
\mathrm{df} = \frac{(v_b + v_w)^2}{v_b^2 / \mathrm{df}_b + v_w^2 / \mathrm{df}_w}.
$$

This equals the Satterthwaite df that lmerTest computes at the true variance
components. When the true effect is $\gamma$, $T$ is approximately noncentral t
with these df and noncentrality $\gamma / \mathrm{SE}$, so

$$
\text{power}(\gamma) = P\big(T > t_c\big) + P\big(T < -t_c\big),
\qquad t_c = t_{1 - \alpha/(2m),\thinspace \mathrm{df}} .
$$

This is an approximation even for normal data, because the Satterthwaite df are
(it is exact only when $R = 1$ or $\rho = 1$). With the default design
($N = 200$), $\mathrm{df} \approx 550$, so it is essentially a z test; the df matter
only for small designs.

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

## Step 5: The power curve over subjects per site

The script repeats Steps 2–4 for each value of `subjects_per_site`, with the number
of sites fixed. For each value it reports $\text{power}(\gamma)$ at
$\gamma = d \cdot \mathrm{SD}$ (the parameter `effect_d`) and $\gamma^\ast$. Under
the balanced-sites simplification with two groups,
$\mathrm{SE} \approx 2\thinspace\mathrm{SD}/\sqrt{Sn}$, so the detectable effect
falls roughly as $1/\sqrt n$: doubling the subjects per site shrinks it by about
30%. At very small $n$ it falls faster, because the degrees of freedom are small
and the t critical value is larger.

## What the power refers to

`target_power` is the power of **one region's test**: the chance of detecting a
group difference of the given size in a particular region. If several regions have
true effects, the chance of detecting at least one of them is higher.

Bonferroni treats the $m$ tests as unrelated, but under the balanced-sites
simplification the region estimates are correlated, with correlation $\rho$. So
Bonferroni is conservative: with 10 regions its family-wise false-positive rate is
0.047, 0.039, and 0.023 for $\rho = 0.2$, $0.5$, and $0.8$, rather than 0.05. A
max-T correction (the critical value of the largest of the $m$ correlated
statistics, e.g. from `multcomp` in R) uses the correlation and would reduce the
detectable effect by about 0.5%, 2.5%, and 8% for those values of $\rho$. The script
uses Bonferroni, which is simpler and slightly conservative.

## Numerical checks

[`checks/check_group_effect.py`](../checks/check_group_effect.py) verifies the
closed-form SE against brute-force GLS for several designs and checks power against
real mixed-model fits (statsmodels `MixedLM`, REML). Data were simulated from the
model with the defaults, with region means, age slopes, and random site offsets
added and an effect of $d = 0.519$ in one region only, then tested at
$\alpha = 0.005$ using `MixedLM`'s z-based p-values (with about 550 df, z and t are
practically the same here). In 1500 simulated studies:

| Region | Rejection rate | Expected |
|---|---|---|
| with the effect (power) | 0.808 | 0.80 (Monte Carlo SE 0.010) |
| without an effect (false positives) | 0.0067 | 0.005 (Monte Carlo SE 0.0018) |

Both agree with the closed-form calculation within Monte Carlo error. An
independent check at $N = 40$, where the df matter more, gave power 0.796
(Monte Carlo SE 0.015) with t-based critical values.

## Assumptions

- **One scan per subject, complete data.** Every subject has all $R$ regions.
  Step 1 relies on this. With missing regions the two parts no longer separate, and
  power would need to be simulated with full mixed-model fits.
- **Equal SD and correlation across regions.** If some regions are noisier, the
  model with one residual variance is misspecified, and its tests are too liberal
  for the noisy regions and too conservative for the others. Changing `sd_total`
  does not show this: it only rescales the answer in outcome units, and $d$ is
  already relative to the SD. To handle unequal SDs, fit region-specific residual
  variances (e.g. `nlme::lme` with `weights = varIdent(form = ~ 1 | region)`) or
  analyze each region separately. The detectable $d$ then applies to each region
  relative to its own SD.
- **Site effects shift all regions equally.** If scanners affect regions
  differently (site-by-region effects), the model above leaves those effects in
  the residual. Under the balanced-sites simplification the group estimate stays
  unbiased, but the estimated residual variance is inflated, so the SE is larger
  and power is lower than calculated (e.g. an SE of 0.147 instead of 0.142 SD when
  site-by-region effects have an SD of 0.3). Adding `site:region` to the model
  removes this cost.
- **Balanced-sites simplification.** Every site recruits equally from each group
  (see [above](#balanced-sites-simplification)). This gives the smallest detectable
  effect for a given total $N$. Imbalance within sites makes it larger.
- **Comparable ages.** Age has the same distribution in every group. Groups that
  differ in age lose some information to the age adjustment.
- **Independent subjects.** Apart from the fixed site effects, subjects at the same
  site are independent.

## References

- Satterthwaite, F. E. (1946). An approximate distribution of estimates of variance
  components. *Biometrics Bulletin*, 2(6), 110–114.
- Zellner, A. (1962). An efficient method of estimating seemingly unrelated
  regressions and tests for aggregation bias. *Journal of the American Statistical
  Association*, 57(298), 348–368.
