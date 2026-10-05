# Derivation: detectable group-by-region effect

This is the derivation behind [`group_effect.py`](../group_effect.py). See the
[README](../README.md#analysis-2-detectable-group-by-region-effect) for how to run it
and which parameters to set.

## Design and model

There are $S$ sites with $n$ subjects each, so $N = Sn$ subjects in total. Each
subject is scanned at one site, and the outcome is measured in $R$ regions.
Subjects belong to one of $K$ groups; group 0 is the reference group. For subject
$i$, let $s(i)$ be its site, $G_{ik}$ the indicator that it is in group $k$
($k = 1, \dots, K-1$), and $x_i$ its centered age.

Steps 1–5 derive the calculation for **one scan per subject**. The planned study
scans each subject $T = 3$ times, a year apart, always at the same site;
[Repeated visits](#repeated-visits) extends the derivation to that design. With
one visit (`visit_years = [0]`), the script reduces exactly to Steps 1–5.

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
are small: at a nominal $\alpha = 0.005$, its true false-positive rate is 0.007
with 56 df (20 sites × 2 subjects, scanned once) and 0.012 with 18 df. With one region ($R = 1$) there is no random
effect to fit; use ordinary regression, for which the same formulas hold.

The script's two noise parameters are $\mathrm{SD}$ (each value in `sd_scenarios`)
and `region_corr` $= \rho$, where

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
$c_b$. In one example, 20 sites of 3 subjects split alternately 2:1 and 1:2 (30 per
group overall), the detectable effect is about 1%, 4%, and 6% larger than for a
balanced design of the same size (15 sites of 4) when $\rho = 0.1$, $0.5$, and
$0.9$. More severe imbalance, such as sites that recruit mostly one group, costs
more.

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

Even for normal data this is approximate, because the Satterthwaite df are. It is
exact only when one part carries all the variance ($R = 1$ or $\rho = 1$). With 20 sites × 10 subjects
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

For 20 sites × 10 subjects scanned once ($m = 10$, $\alpha/m = 0.005$, $p = 0.8$,
two equal groups of 100), this is $(2.81 + 0.84) \times 2/\sqrt{200} \approx 0.52$ SD, matching the
script's $d = 0.519$.

## Step 5: The power curve over subjects per site

The script repeats Steps 2–4 for each value of `subjects_per_site`, with the number
of sites fixed. For each value it reports $\text{power}(\gamma)$ at
the group difference `effect_ppb` and $\gamma^\ast$, for each SD in `sd_scenarios`. Under
the balanced-sites simplification with two groups,
$\mathrm{SE} \approx 2\thinspace\mathrm{SD}/\sqrt{Sn}$, so the detectable effect
falls roughly as $1/\sqrt n$: doubling the subjects per site shrinks it by about
30%. At very small $n$ it falls faster, because the degrees of freedom are small
and the t critical value is larger.

## Repeated visits

Now each subject is scanned at $T$ visits, at times $t_1, \dots, t_T$ (years after
the first scan, `visit_years`), always at the same site. Let $x_{iv}$ be the centered
age of subject $i$ at visit $v$. The analysis model is

```
y ~ 0 + region + region:age_c + site + region:group +
    (1 | subject) + (1 | subject:visit) + (1 | subject:region)
```

that is,

$$
y_{ivr} = \alpha_r + \beta_r x_{iv} + \sum_{k=1}^{K-1} \gamma_{rk} G_{ik} + \eta_{s(i)} +
u_i + p_{ir} + w_{iv} + e_{ivr},
$$

with four independent random terms:

| Term | Variance | What it is |
|---|---|---|
| $u_i$ | $\tau^2$ | Subject effect, shared by all regions and visits |
| $p_{ir}$ | $\psi^2$ | Subject-by-region effect: stable across visits, specific to one region |
| $w_{iv}$ | $\omega^2$ | Visit effect, shared by all regions in one scan (e.g. a whole-scan offset) |
| $e_{ivr}$ | $\sigma^2$ | Visit-by-region noise |

In a single scan, $u_i + w_{iv}$ acts as the subject effect and
$p_{ir} + e_{ivr}$ as the residual of Steps 1–5. The script's parameters are

$$
\mathrm{SD}^2 = \tau^2 + \omega^2 + \psi^2 + \sigma^2,
\qquad
\rho = \frac{\tau^2 + \omega^2}{\mathrm{SD}^2},
\qquad
r = \frac{\tau^2 + \psi^2}{\mathrm{SD}^2},
\qquad
\rho_v = \frac{\omega^2}{\omega^2 + \sigma^2},
$$

(each value in `sd_scenarios`, `region_corr`, `retest_corr`, `visit_region_corr`). $\mathrm{SD}$ and
$\rho$ mean the same as before, for one scan. $r$ is the correlation between one
region's values at two visits of the same subject, and $\rho_v$ is the correlation
between two regions' visit-to-visit changes. Cohen's $d$ is still
$\gamma / \mathrm{SD}$, relative to the single-scan SD.

### Four strata

For subject $i$, collect the $T \times R$ values. Their covariance is
$\tau^2 J_T \otimes J_R + \psi^2 J_T \otimes I_R + \omega^2 I_T \otimes J_R + \sigma^2 I_T \otimes I_R$,
whose eigenspaces are the four combinations of (mean over visits, deviations from
it) and (mean over regions, deviations from it). Each is a stratum with its own
error variance:

| | Mean over regions | Regional deviations |
|---|---|---|
| **Mean over visits** | $\lambda_1 = \tau^2 + \psi^2/R + \omega^2/T + \sigma^2/(TR)$ | $\lambda_2 = \psi^2 + \sigma^2/T$ |
| **Deviations from visit mean** | $\lambda_3 = \omega^2 + \sigma^2/R$ | $\lambda_4 = \sigma^2$ |

(as in Step 1, a single region's deviation has variance $\lambda (1 - 1/R)$).

- **Visit means** (top row). Averaging each subject's values over visits gives
  exactly the single-scan problem of Steps 1–2, with the subject's mean age
  $\bar x_i$ as the age and $\lambda_1$, $\lambda_2$ in place of
  $\tau^2 + \sigma^2/R$ and $\sigma^2$. Because every subject has the same visit
  schedule, $\bar x_i$ is the centered age at the first visit, shifted by a
  constant.
- **Deviations from the visit mean** (bottom row). Region intercepts, group, and
  site are constant across visits, so they drop out. What remains is age: every
  subject's ages deviate from their mean by the same $t_v - \bar t$. These strata
  carry information only about the age slopes: $\kappa / \lambda_3$ about
  $\bar\beta$ and $\kappa / \lambda_4$ about each $\beta_r - \bar\beta$, where
  $\kappa = N \sum_v (t_v - \bar t)^2$.

### Standard error

The group effects get their information from the visit means, but share the age
slopes with the deviation strata, which therefore help slightly. Let $c_b$, $a_b$,
and $b_b$ be the (group, group), (group, age), and (age, age) entries of
$(Z_b^\top Z_b)^{-1}$, with $Z_b$ as in Step 1 and $\bar x_i$ as the age, and
likewise $c_w$, $a_w$, $b_w$ for $Z_w$. Adding the age information to each fit and
inverting (Sherman–Morrison),

$$
\mathrm{Var}(\hat\gamma_{rk}) = v_b + v_w,
\qquad
v_b = \lambda_1 c_b - \frac{\kappa \lambda_1^2 a_b^2}{\lambda_3 + \kappa \lambda_1 b_b},
\qquad
v_w = \Big(1 - \frac{1}{R}\Big) \bigg[\lambda_2 c_w - \frac{\kappa \lambda_2^2 a_w^2}{\lambda_4 + \kappa \lambda_2 b_w}\bigg].
$$

With $T = 1$, $\kappa = 0$ and this is Step 2. It matches the brute-force GLS
variance with all four variance components to machine precision
([checks](#numerical-checks)).

The second terms are small: $a_b$ and $a_w$ measure how much group and age are
confounded, which is little when ages are unrelated to group. Dropping them, and
with the balanced-sites simplification ($c_b \approx c_w = c$),

$$
\mathrm{Var}(\hat\gamma_{rk}) \approx \Big(\tau^2 + \psi^2 + \frac{\omega^2 + \sigma^2}{T}\Big) c
= \mathrm{SD}^2 \Big(r + \frac{1 - r}{T}\Big) c .
$$

That is the variance of a two-sample comparison of each subject's average over
visits. Only the visit-to-visit share of the variance, $1 - r$, is averaged down.
With $r = 0.8$ and $T = 3$, the SE is $\sqrt{0.8 + 0.2/3} = 0.93$ times the
single-scan SE, and no number of visits can take it below $\sqrt{r} = 0.89$ times.
$\rho$ and $\rho_v$ barely matter.

### Degrees of freedom

The Satterthwaite approximation now involves four estimated variances. Each
$\hat\lambda_s$ has approximate variance $2 \lambda_s^2 / \mathrm{df}_s$, with

$$
\mathrm{df}_1 = N - K - S, \quad
\mathrm{df}_2 = (R - 1)(N - K - 1), \quad
\mathrm{df}_3 = N(T - 1), \quad
\mathrm{df}_4 = (R - 1) N (T - 1),
$$

so, writing $v = v_b + v_w$,

$$
\mathrm{df} = \frac{v^2}{\sum_{s=1}^{4} \big(\lambda_s \thinspace \partial v / \partial \lambda_s\big)^2 / \mathrm{df}_s} .
$$

With $T = 1$ this is Step 3. The age slopes are shared between strata, so the split
of degrees of freedom between them is approximate; this makes no visible
difference. Fitting the model with lmerTest to simulated data and evaluating the
formula at the REML estimates from the same fit:

| Design | lmerTest SE | Formula SE | lmerTest df | Formula df |
|---|---|---|---|---|
| 14 sites × 22 subjects × 3 visits, $R = 10$ (planned) | 0.104858 | 0.104858 | 934.2 | 934.7 |
| 14 sites × 6 subjects × 3 visits, $R = 5$ | 0.200384 | 0.200384 | 175.7 | 175.8 |
| 5 sites × 4 subjects × visits at 0, 1, 3 years, $R = 6$ | 0.363358 | 0.363358 | 27.7 | 28.4 |
| 8 sites × 2 subjects × 2 visits, $R = 3$ | 0.416810 | 0.416810 | 11.1 | 11.3 |

Steps 4 and 5 are unchanged.

## What the power refers to

`target_power` is the power of **one region's test**: the chance of detecting a
group difference of the given size in a particular region. If several regions have
true effects, the chance of detecting at least one of them is higher.

Bonferroni is valid whatever the correlation between tests, but it ignores that
correlation. Under the balanced-sites simplification the region estimates are
correlated, with correlation approximately $\rho$ for one scan per subject (with
repeated visits, approximately $(\tau^2 + \omega^2/T) / (\tau^2 + \psi^2 + (\omega^2 + \sigma^2)/T)$,
the correlation between two regions' averages over visits), so Bonferroni is
conservative.
With 2 groups, 10 regions, and large df, its family-wise false-positive rate is
0.047, 0.039, and 0.023 for $\rho = 0.2$, $0.5$, and $0.8$, rather than 0.05. A
max-T correction (based on the largest of the $m$ correlated test statistics) uses
the correlation and would reduce the detectable effect by about 0.5%, 2.5%, and 8%
for those values of $\rho$. In R, `multcomp` computes it, but for mixed models it
uses z statistics by default, which are too liberal at small df. The script uses
Bonferroni, which is simpler and slightly conservative.

## Numerical checks

[`checks/check_group_effect.py`](../checks/check_group_effect.py) verifies the
closed-form SE against brute-force GLS for several designs and checks power against
real mixed-model fits (statsmodels `MixedLM`, REML). Data are simulated from the
model with region means, age slopes, and random site offsets added and an effect,
in one region only, of the size the closed form says is detectable with 80% power.
Each fit's Wald statistic is compared with the t critical value at the closed-form
Satterthwaite df, as lmerTest would do, at $\alpha = 0.005$:

| Design | Studies | Power (expected 0.80) | False positives (expected 0.005) |
|---|---|---|---|
| 20 sites × 10 subjects, one scan, $d = 0.519$ | 1500 | 0.797 (MC SE 0.010) | 0.0027 (MC SE 0.0018) |
| 20 sites × 2 subjects, one scan (56 df), $d = 1.222$ | 700 | 0.787 (MC SE 0.015) | 0.0057 (MC SE 0.0027) |

To reproduce: `python3 checks/check_group_effect.py 1500 10` and
`python3 checks/check_group_effect.py 700 2` (about 10 and 4 minutes).

Both agree with the closed-form calculation within Monte Carlo error.

[`checks/check_group_effect_lmer.R`](../checks/check_group_effect_lmer.R) does the
same for repeated visits, with lmerTest's own Satterthwaite p-values. Data are
simulated from the repeated-visits model ($\rho = 0.5$, $r = 0.8$, $\rho_v = 0.5$,
$R = 10$, visits at 0, 1, 2 years), at $\alpha = 0.005$:

| Design | Studies | Power (expected 0.80) | False positives (expected 0.005) |
|---|---|---|---|
| 14 sites × 22 subjects × 3 visits (planned, 901 df), $d = 0.389$ | 1200 | 0.812 (MC SE 0.012) | 0.0042 (MC SE 0.0020) |
| 14 sites × 6 subjects × 3 visits (212 df), $d = 0.753$ | 1200 | 0.792 (MC SE 0.012) | 0.0067 (MC SE 0.0020) |
| 6 sites × 2 subjects × 3 visits, $R = 4$ (10 df), $d = 2.64$ | 1200 | 0.757 (MC SE 0.012) | 0.0083 (MC SE 0.0020) |

The first two agree within Monte Carlo error. The last, a deliberately tiny design
of 12 subjects, shows that with about 10 df the calculation is somewhat optimistic,
because the variance estimates are then very imprecise; this does not affect
designs of realistic size.

## Assumptions

- **Complete data.** Every subject has all $R$ regions at every visit. Step 1 and
  the four strata rely on this. With missing regions or missed visits the parts no
  longer separate, and power would need to be simulated with full mixed-model fits.
  As a rough guide, dropout removes visits, which costs little (see
  [Standard error](#standard-error)), but a subject lost entirely costs a whole
  subject.
- **Same visit schedule, same site.** Every subject is scanned at the same times
  (`visit_years`) and always at the same site. Visits a little early or late
  make almost no difference.
- **Visit-to-visit changes are independent.** The visit effects $w_{iv}$ and
  $e_{ivr}$ are independent across visits, and the group difference is the same at
  every visit. The model has no subject-specific rates of change (random age
  slopes). If subjects do change at different rates, that adds to the
  differences between their visit averages, which the variance of the visit-mean
  stratum (estimated from the data) absorbs, so the group test stays
  approximately valid. If a group difference in rate of change is of interest,
  that is a different test (`region:group:age_c`), not covered here.
- **Equal SD and correlation across regions.** If some regions are noisier, the
  model with one residual variance is misspecified, and its tests are too liberal
  for the noisy regions and too conservative for the others. Changing the SD
  scenarios does not show this: they only rescale the answer in ppb, and $d$ is
  already relative to the SD. To handle unequal SDs, fit region-specific residual
  variances (e.g. `nlme::lme` with `weights = varIdent(form = ~ 1 | region)`) or
  analyze each region separately. The detectable $d$ then applies to each region
  relative to its own SD. A separate analysis per region has only $N - K - S$ df
  (18 with 2 subjects per site at 20 sites, against the script's 56), so at small
  $n$ its power is a little lower than the script reports.
- **Site effects shift all regions equally.** If scanners affect regions
  differently (site-by-region effects), the model above leaves those effects in
  the residual. Under the balanced-sites simplification the group estimate stays
  unbiased, but the estimated residual variance is inflated, so the SE is larger
  and power is lower than calculated (e.g. an SE of 0.147 instead of 0.142 SD at
  20 sites × 10 subjects scanned once, with $\rho = 0.5$, when
  site-by-region effects have an SD of 0.3 SD). Adding `site:region` to the model
  removes this cost, though it also lowers the df to $(R - 1)(N - K - S)$ for the
  regional deviations, so at small $n$ the script slightly overstates the df.
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
