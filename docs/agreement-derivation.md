# Derivation: do scanners preserve regional differences as well as repeat scans?

This is the derivation behind [`agreement.py`](../agreement.py). See the
[README](../README.md#analysis-1-scanner-agreement-in-regional-differences) for how
to run it and which parameters to set.

**The question.** The same subjects are scanned twice at each of 3 sites. Within
each scan, take every pairwise difference between the 5 regions (10 differences).
Are a subject's regional differences as close between scanners as between two
repeat scans on one scanner? Between-scanner disagreement includes both fixed
scanner offsets (a scanner shifting one region relative to another in everyone)
and subject-specific scanner scatter.

## Notation and model

Subjects are $i = 1, \dots, N$, regions $r = 1, \dots, 5$, sites
$s = 1, 2, 3$, and repeats $t = 1, \dots, n$. The value of region $r$ in scan $t$ of
subject $i$ at site $s$ is

$$
y_{irst} = \mu_{ir} + o_{rs} + w_{ist} + \eta_{irs} + \varepsilon_{irst},
$$

| Term | Meaning | Variance |
|---|---|---|
| $\mu_{ir}$ | the subject's true value (including age effects and the subject's own regional pattern) | fixed per subject |
| $o_{rs}$ | fixed scanner offset in region $r$ at site $s$, the same for every subject | fixed (`offset_pattern`, scaled) |
| $w_{ist}$ | whole-scan offset, shared by all regions in one scan | any |
| $\eta_{irs}$ | subject-by-site deviation, the same in both repeats at that site | $\sigma_\eta^2$ (`sd_site`) |
| $\varepsilon_{irst}$ | region-level scan-to-scan noise | $\sigma^2$ (`sd_noise`) |

$\eta$ and $\varepsilon$ are normal, independent across subjects, regions, sites and
scans, with the same variance everywhere.

**Within-scan differences.** For a contrast vector $\mathbf m$ over regions whose
entries sum to 0 (for example $\mathbf m = e_a - e_b$ for the difference between
regions $a$ and $b$), the contrast of one scan is

$$
c_{ist} = \mathbf m^\top \mathbf y_{ist}
= \mathbf m^\top \boldsymbol\mu_{i} + \mathbf m^\top \mathbf o_{s}
+ \mathbf m^\top \boldsymbol\eta_{is} + \mathbf m^\top \boldsymbol\varepsilon_{ist}.
$$

The whole-scan offset $w_{ist}$ cancels, because $\mathbf m$ sums to 0. The noise
terms have variances $\lVert\mathbf m\rVert^2\sigma_\eta^2$ and
$\lVert\mathbf m\rVert^2\sigma^2$ ($\lVert\mathbf m\rVert^2 = 2$ for a pairwise
difference).

## The target ratio

For one contrast, the squared difference between two repeat scans at one site,
and between single scans at two sites $j$ and $k$, have expectations

$$
E\big[(c_{ist} - c_{ist'})^2\big] = 2\lVert\mathbf m\rVert^2\sigma^2,
\qquad
E\big[(c_{ijt} - c_{ikt})^2\big] = (\mathbf m^\top \mathbf o_j - \mathbf m^\top \mathbf o_k)^2
+ 2\lVert\mathbf m\rVert^2(\sigma_\eta^2 + \sigma^2).
$$

Averaging the second over the 3 pairs of sites, the ratio is

$$
R = 1 + \frac{\sigma_\eta^2 + \bar\Delta^2 / (2\lVert\mathbf m\rVert^2)}{\sigma^2},
\qquad
\bar\Delta^2 = \tfrac{1}{3}\sum_{j<k}(\mathbf m^\top \mathbf o_j - \mathbf m^\top \mathbf o_k)^2 .
$$

$R \ge 1$, because a between-scanner difference always includes the repeat noise.
$R = 1$ means scanners add nothing beyond repeat noise. The claim to show is
$R < R_0$ for a margin $R_0 > 1$; "the same size or smaller" cannot be shown
without a margin.

## Step 1: Mean squares

For one contrast, let $\bar c_{is}$ be the mean of subject $i$'s $n$ repeats at
site $s$, and $\bar c_{i}$ the mean of those over sites. Define

$$
\mathrm{MS}_{\text{site}} = \frac{n}{N(S-1)} \sum_{i}\sum_{s} (\bar c_{is} - \bar c_{i})^2,
\qquad
\mathrm{MS}_{\text{rep}} = \frac{1}{N S (n-1)} \sum_{i}\sum_{s}\sum_{t} (c_{ist} - \bar c_{is})^2 .
$$

$\mathrm{MS}_{\text{site}}$ is the variation between sites within subjects: in a
two-way ANOVA with subject, site and subject-by-site terms, it is the site plus the
interaction sum of squares, with $N(S-1)$ degrees of freedom. It does **not**
remove the site means, so fixed offsets count as disagreement.
$\mathrm{MS}_{\text{rep}}$ is the residual (repeat) mean square, with $NS(n-1)$
degrees of freedom. The subject's own pattern $\boldsymbol\mu_i$ cancels in both.

With $S = 3$, $\sum_s (o_s - \bar o)^2 = \bar\Delta^2$ for the contrast's offsets,
so

$$
E[\mathrm{MS}_{\text{rep}}] = \lVert\mathbf m\rVert^2\sigma^2,
\qquad
E[\mathrm{MS}_{\text{site}}] = \lVert\mathbf m\rVert^2\sigma^2 + n\big(\lVert\mathbf m\rVert^2\sigma_\eta^2 + \bar\Delta^2/2\big),
$$

and the ratio of expectations is $\tau = 1 + n(R - 1)$.

## Step 2: The test for one contrast

The two mean squares are independent: one uses the site means, the other the
deviations of repeats from them. $\mathrm{MS}_{\text{rep}}$ is a scaled
$\chi^2_{NS(n-1)}$. Without offsets, $\mathrm{MS}_{\text{site}}$ is a scaled
$\chi^2_{N(S-1)}$, so

$$
F = \mathrm{MS}_{\text{site}} / \mathrm{MS}_{\text{rep}} \sim \tau \cdot F_{N(S-1),\thinspace NS(n-1)} .
$$

Test $H_0: R \ge R_0$ against $H_1: R < R_0$: reject when
$F < \tau_0\thinspace q_\alpha$, with $\tau_0 = 1 + n(R_0 - 1)$ and $q_\alpha$ the
lower $\alpha$ quantile of $F_{N(S-1),\thinspace NS(n-1)}$. Equivalently, the
smallest margin a study passes is $1 + (F/q_\alpha - 1)/n$.

**With fixed offsets**, $\mathrm{MS}_{\text{site}}$ is a scaled noncentral
$\chi^2$. For a given expectation, the noncentral version is *less* spread than
the central one, so it falls below the critical value less often at the null
boundary: the test stays at or below its level whatever the mix of offsets and
scatter (checked by simulation; the rejection rate at $R = R_0$ was 0.034–0.051
with offsets making up 0–100% of $R - 1$).

## Step 3: Pooled and every-pair tests

**Pooled (primary).** The 10 pairwise differences span a 4-dimensional space of
contrasts. Take 4 orthonormal contrasts in that space. Because $\eta$ and
$\varepsilon$ are independent with equal variance across regions, the 4 contrasts
have independent noise, so their mean squares can be added: the pooled $F$ has
$4N(S-1)$ and $4NS(n-1)$ degrees of freedom and the same test applies, with $R$ the
average of the 4 contrasts' ratios. For any set of scans, the pooled sum of
squares over the orthonormal contrasts equals the sum over the 10 pairwise
differences divided by 5, so the pooled test is the test on all pairwise
differences together.

Power without offsets is exact:
$\Pr\big(F_{d_1, d_2} < (\tau_0/\tau)\thinspace q_\alpha\big)$. With offsets, the
script simulates.

**Every pair (secondary).** Apply the one-contrast test to each of the 10 pairwise
differences, with margin `R0_pair`; the claim holds only if all 10 pass. This is an
intersection-union test, so it needs no multiplicity correction (Berger, 1982).
It is decided by the worst pair, so fixed offsets concentrated in one region hurt
it much more than the pooled test.

## Step 4: Simulation

[`agreement.py`](../agreement.py) builds each planning scenario from a pooled $R$
and the share of $R - 1$ due to fixed offsets (`offset_share`): the offsets take the
shape of `offset_pattern`, scaled so that they contribute that share, and
$\sigma_\eta$ supplies the rest. Each simulated study draws $\eta$ and
$\varepsilon$ for $N$ subjects, adds the offsets, forms the contrasts, and records
the smallest margin the pooled and every-pair tests pass. Power at a margin is the
fraction of studies that pass it; the smallest margin with power $p$ is the
$p$-quantile. Subject means, age and whole-scan offsets are not simulated, because
they cancel.

## Assumptions

- **Complete, balanced data:** every subject has all 5 regions in $n$ scans at each
  of the 3 sites.
- **Equal noise across regions, sites and subjects, independent across regions.**
  The pooled $F$ is exact only then. If some regions are noisier, or region noise
  is correlated, the pooled test is approximate and the per-pair tests (each exact
  on its own) are the safer evidence. If one site is noisier, its between-site
  differences are larger and $R$ is larger, which the test correctly counts as
  worse agreement; but the repeat noise is then averaged over sites.
- **Normality.** Variance-ratio tests are sensitive to heavy tails; check the
  repeat differences for outliers (for example motion-corrupted scans).
- **Repeats are true repeats:** reposition the subject between the two scans at a
  site, so that $\varepsilon$ includes positioning effects. Otherwise the repeat
  noise is too small and $R$ looks worse than it is.

## Numerical checks

[`checks/check_agreement.py`](../checks/check_agreement.py) checks that:

1. the mean squares equal those of a two-way ANOVA fit with statsmodels (site plus
   subject-by-site sums of squares, and the residual);
2. the pooled sum of squares over orthonormal contrasts equals the sum over the 10
   pairwise differences divided by 5;
3. simulated pooled power without offsets matches the exact formula (0.797 vs 0.797,
   0.631 vs 0.630, 0.915 vs 0.914);
4. at $R = R_0$ the pooled test rejects at most 0.051 of the time (Monte Carlo SE
   0.0015) for $N$ = 6, 20 and 50, with offsets making up 0, 50 or 100% of
   $R - 1$, concentrated in one region or spread over several; the every-pair test
   at the worst pair's $R$ rejects at most 0.001.

## References

- Berger, R. L. (1982). Multiparameter hypothesis testing and acceptance sampling.
  *Technometrics*, 24(4), 295–300.
