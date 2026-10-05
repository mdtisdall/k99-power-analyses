# Power for cross-site equivalence of region-specific age effects.
#
# Design: N subjects, each scanned at all 3 sites, 5 regions per scan, complete
# data. The single-site model fit at each site is
#   y ~ 0 + region + region:age_c + (1 | subject)
# With complete, balanced data its region slopes equal per-region OLS slopes, so
# the slope difference between two sites for a region is the slope of the
# per-subject difference score D on age_c. "Same across sites" requires all
# 5 regions x 3 site pairs = 15 TOSTs (90% CI inside (-Delta, +Delta)) to pass.
# Power = P(all 15 pass). See README.md.
#
# All parameter values below are placeholders. Replace them with estimates from
# prior (ideally traveling-subject / test-retest) data.

# ---- Parameters -------------------------------------------------------------

Delta   <- 0.003               # equivalence margin, outcome units per year
s_ss    <- 0.05                # SD of subject:site offset (shared across regions in a scan)
sigma   <- c(.10, .10, .10)    # region-level measurement-noise SD at each site
age_min <- 50                  # age distribution: uniform(age_min, age_max)
age_max <- 80
dslope  <- matrix(0, 5, 3)     # true slope deviation, regions (rows) x sites (cols)
alpha   <- 0.05                # one-sided level for each TOST (gives a 90% CI)
N_grid  <- seq(20, 200, by = 20)
n_sims  <- 2000
seed    <- 1

# ---- Simulation -------------------------------------------------------------

# Simulate only the terms that do not cancel in between-site differences:
# subject:site offsets, measurement noise, and site slope deviations.
sim_once <- function(N, Delta, s_ss = .05, sigma = c(.10, .10, .10),
                     dslope = matrix(0, 5, 3), alpha = .05,
                     age_min = 50, age_max = 80) {
  age_c <- runif(N, age_min, age_max); age_c <- age_c - mean(age_c)
  e <- array(0, c(N, 5, 3))
  for (k in 1:3)
    e[, , k] <- rnorm(N, 0, s_ss) +                     # subject:site offset
                matrix(rnorm(N * 5, 0, sigma[k]), N) +  # measurement noise
                outer(age_c, dslope[, k])               # true site slope deviation
  tq <- qt(1 - alpha, N - 2)
  for (r in 1:5) for (p in list(c(1, 2), c(1, 3), c(2, 3))) {
    D  <- e[, r, p[1]] - e[, r, p[2]]
    cf <- summary(lm(D ~ age_c))$coefficients["age_c", 1:2]
    if (!(cf[1] - tq * cf[2] > -Delta && cf[1] + tq * cf[2] < Delta)) return(FALSE)
  }
  TRUE
}

# ---- Power curve ------------------------------------------------------------

set.seed(seed)
power <- sapply(N_grid, function(N)
  mean(replicate(n_sims, sim_once(N, Delta, s_ss, sigma, dslope, alpha,
                                  age_min, age_max))))

result <- data.frame(N = N_grid, power = power)
print(result, row.names = FALSE)
write.csv(result, "power_curve.csv", row.names = FALSE)
