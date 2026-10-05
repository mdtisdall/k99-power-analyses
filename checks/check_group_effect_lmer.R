# Power check of docs/group-effect-derivation.md (repeated visits) with lmerTest.
#
# Simulates studies from the repeated-visits model, fits
#   y ~ 0 + region + region:age_c + site + region:group +
#       (1 | subject) + (1 | subject:visit) + (1 | subject:region)
# with lmerTest (REML, Satterthwaite df), and counts how often region 0's group
# effect (true effect `effect`, in SD units) and region 1's (true effect 0) are
# significant at `alpha`. For the planned design, group_effect.py gives a
# detectable d of 0.3886 at 80% power, so the expected power is 0.80 and the
# expected false-positive rate 0.005.
#
# Run from the repository root (about 1.5 s per study for the planned design):
#   Rscript checks/check_group_effect_lmer.R S n R region_corr retest_corr \
#     visit_region_corr visit_years effect alpha n_sims seed
# e.g. Rscript checks/check_group_effect_lmer.R 14 22 10 0.5 0.8 0.5 0,1,2 0.3886 0.005 1200 1
# Needs R with lme4 and lmerTest, e.g. with Nix:
#   nix-shell -p 'rWrapper.override{packages=with rPackages;[lme4 lmerTest];}'
suppressMessages(library(lmerTest))
a <- commandArgs(TRUE)
S <- as.integer(a[1]); n <- as.integer(a[2]); R <- as.integer(a[3])
rho <- as.numeric(a[4]); rt <- as.numeric(a[5]); vrc <- as.numeric(a[6])
yrs <- as.numeric(strsplit(a[7], ",")[[1]]); eff <- as.numeric(a[8]); alpha <- as.numeric(a[9])
nsim <- as.integer(a[10]); set.seed(as.integer(a[11]))
T <- length(yrs); N <- S * n; tt <- yrs - mean(yrs)
visit <- (1 - rt); sw <- vrc * visit; se <- visit - sw; tu <- rho - sw; sb <- (1 - rho) - se
site <- rep(0:(S - 1), each = n); group <- (0:(N - 1)) %% 2
# long format: region fastest, then visit, then subject
I <- rep(0:(N - 1), each = T * R); Tt <- rep(rep(0:(T - 1), each = R), N); Rr <- rep(0:(R - 1), N * T)
hit <- fp <- 0
for (s in 1:nsim) {
  x <- runif(N, 25, 65); x <- x - mean(x)
  y <- sqrt(tu) * rnorm(N)[I + 1] + sqrt(sb) * rnorm(N * R)[I * R + Rr + 1] +
       sqrt(sw) * rnorm(N * T)[I * T + Tt + 1] + sqrt(se) * rnorm(N * T * R) +
       rnorm(R)[Rr + 1] + 0.4 * rnorm(S)[site[I + 1] + 1] +
       seq(-.02, .02, length.out = R)[Rr + 1] * (x[I + 1] + tt[Tt + 1]) +
       eff * (Rr == 0) * group[I + 1]
  d <- data.frame(y, region = factor(Rr), age_c = x[I + 1] + tt[Tt + 1], site = factor(site[I + 1]),
                  group = group[I + 1], subject = factor(I), visit = factor(Tt))
  f <- lmer(y ~ 0 + region + region:age_c + site + region:group + (1|subject) + (1|subject:visit) + (1|subject:region), d)
  co <- summary(f)$coefficients
  hit <- hit + (co["region0:group", "Pr(>|t|)"] < alpha)
  fp <- fp + (co["region1:group", "Pr(>|t|)"] < alpha)
}
cat(sprintf("S=%d n=%d R=%d years=%s effect=%.4f: power %.3f (MC SE %.3f), false positives %.4f (MC SE %.4f), %d fits\n",
            S, n, R, a[7], eff, hit / nsim, sqrt(.16 / nsim), fp / nsim, sqrt(alpha * (1 - alpha) / nsim), nsim))
