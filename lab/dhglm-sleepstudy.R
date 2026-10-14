################################################################################
#
# A double hierarchical model on the sleep study data (October 2026): random
# effects on the residual precision, which R-INLA cannot fit. Morales-Otero,
# Gomez-Rubio and Nunez-Anton (2022, Statistics and Computing 32(4),
# arXiv:2201.07493) fitted it with INLA inside adaptive multiple importance
# sampling over the 18 log-precisions: 21 rounds of 1,000 INLA fits, a first
# proposal whose results they discarded, and a final effective sample size of
# 2.0. Code and data: https://github.com/becarioprecario/DHGLM-INLA, folder
# sleepstudy_rcoef.
#
# Their model (equation 3), with reaction times in seconds:
#   y_ij ~ N(beta0 + beta_i day_ij, 1 / tau_i), log tau_i = gamma + u_i,
#   beta_i ~ N(0, 1 / tau_beta), u_i ~ N(0, 1 / tau_u),
#   beta0, gamma ~ N(0, precision 0.001), tau_beta ~ Gamma(1, 0.00005),
#   tau_u ~ Gamma(0.01, 0.01).
# The paper's text gives Gamma(1, 0.00005) for tau_u, but its JAGS and INLA
# code both use Gamma(0.01, 0.01), and so does this script. The random slopes
# have mean 0, so they also carry the average effect of the days.
#
# Both random effects integrate out inside logpost. Given tau_i, subject i's
# 10 responses are normal with covariance I / tau_i + d d' / tau_beta, a
# rank-one update. Then u_i (as v_i = log tau_i) integrates out by adaptive
# Gauss-Hermite quadrature, centred at the mode of its integrand and scaled
# by its curvature, since a fixed rule fails when tau_u is small. That leaves
# 4 parameters for inla_vanilla(): beta0, gamma, log tau_beta and log tau_u.
# The subject effects come from sampling(): at each draw, v_i from its
# conditional density on a grid, then beta_i from its normal given v_i.
# Reference: Stan (NUTS, non-centred) on all 40 parameters.
#
# Needs lme4 (for the data) and rstan. Run from the package root:
#   Rscript lab/dhglm-sleepstudy.R
#
# Results (October 2026; Stan 4 chains of 3,000 kept draws, R-hat 1.004, no
# divergences): inla_vanilla() matches Stan. The four outer parameters are
# within 0.05 posterior SDs in mean and 0.17 in the 2.5% and 97.5% quantiles,
# with SDs within 5%. The 36 subject effects (log tau_i and beta_i) are within
# 0.03 SDs in mean, with SDs within 4%. On the natural scale the mean of
# tau_u is 2.59 against Stan's 2.72 (the paper: 2.12 by AMIS, 2.66 by JAGS),
# though its SD, which its long right tail drives, is 1.9 against 3.1 (the
# paper: 0.93 and 2.3). The fit took 0.27 s (389 calls of logpost) and the
# subject draws 5 s, against 4.5 s of sampling for Stan after a minute of
# compiling: this model is small for Stan, and the contrast is with the
# paper's 21,000 INLA fits. Every check in summary() passes, and 41 quadrature
# nodes change logpost at the mode by 2e-8.
#
################################################################################

## ----- Configuration ---------------------------------------------------------

smoke <- FALSE
suppressMessages(devtools::load_all(".", quiet = TRUE))
data(sleepstudy, package = "lme4")
n_iter <- if (smoke) 1000 else 6000 # Stan iterations per chain, half warm-up
n_draws <- if (smoke) 2000 else 10000 # copula draws
n_nodes <- 15 # adaptive Gauss-Hermite nodes per subject
options(mc.cores = 4)

## =============================================================================
## Model and data
## =============================================================================

y <- sleepstudy$Reaction / 1000
day <- sleepstudy$Days
id <- as.integer(sleepstudy$Subject)
n_sub <- max(id)
n_obs <- tabulate(id)
dd <- sum(day[id == 1]^2) # d'd, the same for every subject (days 0 to 9)

## ----- Gauss-Hermite rule ----------------------------------------------------

# Nodes and weights for the standard normal, from the Golub-Welsch algorithm
gh_normal <- function(k) {
  j <- seq_len(k - 1)
  jacobi <- matrix(0, k, k)
  jacobi[cbind(j, j + 1)] <- sqrt(j)
  jacobi[cbind(j + 1, j)] <- sqrt(j)
  e <- eigen(jacobi, symmetric = TRUE)
  list(z = e$values, w = e$vectors[1, ]^2)
}
gh <- gh_normal(n_nodes)

## ----- Subject log-likelihood given log tau_i --------------------------------

# Sufficient statistics of the residuals y - beta0: S_i = sum r^2 and
# P_i = sum d r, per subject
subject_stats <- function(beta0) {
  r <- y - beta0
  list(s = rowsum(r^2, id)[, 1], p = rowsum(day * r, id)[, 1])
}

# log N(y_i; beta0, a I + b d d') at a = exp(-v), b = 1 / tau_beta, for a
# vector v (one value per subject, or one row of nodes per subject)
loglik_given_v <- function(v, st, tau_beta) {
  a <- exp(-v)
  b <- 1 / tau_beta
  -0.5 *
    n_obs *
    log(2 * pi * a) -
    0.5 * log1p(b * dd / a) -
    0.5 * (st$s - b * st$p^2 / (a + b * dd)) / a
}

# The integrand over v_i on the log scale: the subject's likelihood times the
# N(gamma, 1 / tau_u) density of v_i
log_integrand <- function(v, st, theta) {
  loglik_given_v(v, st, exp(theta[["log_tau_beta"]])) +
    dnorm(v, theta[["gamma"]], exp(-theta[["log_tau_u"]] / 2), log = TRUE)
}

# The mode and curvature scale of each subject's integrand, by Newton steps
# with central differences (a fixed number, so the result is smooth in theta)
integrand_mode <- function(st, theta, h = 1e-3, steps = 10) {
  v <- log(n_obs / st$s)
  for (k in seq_len(steps)) {
    f0 <- log_integrand(v, st, theta)
    fp <- log_integrand(v + h, st, theta)
    fm <- log_integrand(v - h, st, theta)
    d1 <- (fp - fm) / (2 * h)
    d2 <- (fp - 2 * f0 + fm) / h^2
    v <- v - d1 / pmin(d2, -1e-8)
  }
  f0 <- log_integrand(v, st, theta)
  d2 <- (log_integrand(v + h, st, theta) -
    2 * f0 +
    log_integrand(v - h, st, theta)) /
    h^2
  list(v = v, scale = 1 / sqrt(-pmin(d2, -1e-8)))
}

## ----- Log-posterior ---------------------------------------------------------

logpost <- function(theta) {
  st <- subject_stats(theta[["beta0"]])
  m <- integrand_mode(st, theta)
  # Adaptive Gauss-Hermite: nodes v_ik = mode_i + scale_i z_k
  v <- outer(m$v, rep(1, n_nodes)) + outer(m$scale, gh$z)
  st_k <- list(s = rep(st$s, n_nodes), p = rep(st$p, n_nodes))
  h <- matrix(log_integrand(as.vector(v), st_k, theta), n_sub)
  terms <- sweep(h, 2, log(gh$w) - dnorm(gh$z, log = TRUE), "+")
  top <- apply(terms, 1, max)
  loglik <- sum(log(m$scale) + top + log(rowSums(exp(terms - top))))
  tau_beta <- exp(theta[["log_tau_beta"]])
  tau_u <- exp(theta[["log_tau_u"]])
  loglik +
    dnorm(theta[["beta0"]], 0, sqrt(1000), log = TRUE) +
    dnorm(theta[["gamma"]], 0, sqrt(1000), log = TRUE) +
    dgamma(tau_beta, 1, 0.00005, log = TRUE) +
    theta[["log_tau_beta"]] +
    dgamma(tau_u, 0.01, 0.01, log = TRUE) +
    theta[["log_tau_u"]]
}

init <- c(
  beta0 = mean(y),
  gamma = log(1 / var(y)),
  log_tau_beta = log(1 / 0.01^2),
  log_tau_u = log(2)
)

## =============================================================================
## Fits
## =============================================================================

## ----- inla_vanilla() --------------------------------------------------------

t_fit <- system.time(fit <- inla_vanilla(logpost, init))[["elapsed"]]

# Quadrature check at the mode: 15 against 41 nodes
gh_15 <- gh
gh <- gh_normal(41)
n_nodes <- 41
lp_41 <- logpost(fit$mode)
gh <- gh_15
n_nodes <- 15
quad_error <- lp_41 - logpost(fit$mode)

# Subject effects: v_i from its conditional density on a grid of 201 points
# around its mode, then beta_i | v_i, which is normal
draw_subjects <- function(theta) {
  st <- subject_stats(theta[["beta0"]])
  m <- integrand_mode(st, theta)
  grid <- seq(-8, 8, length.out = 201)
  v <- outer(m$v, rep(1, 201)) + outer(m$scale, grid)
  st_g <- list(s = rep(st$s, 201), p = rep(st$p, 201))
  h <- matrix(log_integrand(as.vector(v), st_g, theta), n_sub)
  prob <- exp(h - apply(h, 1, max))
  pick <- apply(prob, 1, function(p) sample.int(201, 1, prob = p))
  step <- m$scale * (grid[2] - grid[1])
  v_i <- v[cbind(seq_len(n_sub), pick)] + step * (runif(n_sub) - 0.5)
  tau_i <- exp(v_i)
  tau_beta <- exp(theta[["log_tau_beta"]])
  prec <- tau_beta + tau_i * dd
  beta_i <- rnorm(n_sub, tau_i * st$p / prec, 1 / sqrt(prec))
  c(v_i, beta_i)
}
set.seed(1)
t_draws <- system.time(
  sub_iv <- sampling(fit, nsamp = n_draws, fun = draw_subjects)
)[["elapsed"]]
natural_iv <- sampling(fit, nsamp = n_draws, fun = function(theta) {
  c(
    tau_beta = exp(theta[["log_tau_beta"]]),
    tau_u = exp(theta[["log_tau_u"]])
  )
})

## ----- Stan ------------------------------------------------------------------

stan_code <- "
data {
  int<lower=1> N;
  int<lower=1> J;
  array[N] int<lower=1, upper=J> id;
  vector[N] day;
  vector[N] y;
}
parameters {
  real beta0;
  real gamma;
  real<lower=0> tau_beta;
  real<lower=0> tau_u;
  vector[J] z_beta;
  vector[J] z_u;
}
transformed parameters {
  vector[J] beta = z_beta / sqrt(tau_beta);
  vector[J] log_tau = gamma + z_u / sqrt(tau_u);
}
model {
  y ~ normal(beta0 + beta[id] .* day, exp(-0.5 * log_tau[id]));
  beta0 ~ normal(0, sqrt(1000));
  gamma ~ normal(0, sqrt(1000));
  tau_beta ~ gamma(1, 0.00005);
  tau_u ~ gamma(0.01, 0.01);
  z_beta ~ std_normal();
  z_u ~ std_normal();
}
generated quantities {
  real log_tau_beta = log(tau_beta);
  real log_tau_u = log(tau_u);
}
"
stan_mod <- rstan::stan_model(model_code = stan_code)
t_stan <- system.time(
  stan_fit <- rstan::sampling(
    stan_mod,
    data = list(N = length(y), J = n_sub, id = id, day = day, y = y),
    chains = 4,
    iter = n_iter,
    seed = 2026,
    refresh = 0,
    control = list(adapt_delta = 0.99)
  )
)[["elapsed"]]
stan_draws <- as.matrix(stan_fit)

## =============================================================================
## Comparison
## =============================================================================

## ----- Outer parameters ------------------------------------------------------

par_names <- names(init)
ref <- cbind(
  mean = colMeans(stan_draws[, par_names]),
  sd = apply(stan_draws[, par_names], 2, sd),
  t(apply(stan_draws[, par_names], 2, quantile, c(0.025, 0.975)))
)
s <- fit$summary[par_names, ]
cat("\n== Stan reference ==\n")
print(signif(ref, 3))
cat("\n== inla_vanilla() against Stan, in Stan's posterior SDs ==\n")
print(round(
  cbind(
    mean = (s[, "Mean"] - ref[, "mean"]) / ref[, "sd"],
    sd_ratio = s[, "SD"] / ref[, "sd"],
    q025 = (s[, "2.5%"] - ref[, "2.5%"]) / ref[, "sd"],
    q975 = (s[, "97.5%"] - ref[, "97.5%"]) / ref[, "sd"]
  ),
  3
))

# The precisions on their natural scale, next to the paper's Table 8 (AMIS
# with INLA, and their JAGS run)
natural <- function(draws) c(mean = mean(draws), sd = sd(draws))
cat("\n== tau_beta and tau_u: mean (SD) ==\n")
print(round(
  rbind(
    inla_vanilla = c(
      natural(natural_iv[, "tau_beta"]),
      natural(natural_iv[, "tau_u"])
    ),
    stan = c(
      natural(stan_draws[, "tau_beta"]),
      natural(stan_draws[, "tau_u"])
    ),
    paper_amis = c(8240.22, 2742.75, 2.1222, 0.9348),
    paper_jags = c(8002.25, 2898.70, 2.6565, 2.2836)
  ),
  3
))

## ----- Subject effects -------------------------------------------------------

sub_stan <- stan_draws[, c(
  paste0("log_tau[", 1:n_sub, "]"),
  paste0("beta[", 1:n_sub, "]")
)]
sd_ref <- apply(sub_stan, 2, sd)
gap <- rbind(
  mean = abs(colMeans(sub_iv) - colMeans(sub_stan)) / sd_ref,
  sd_ratio = apply(sub_iv, 2, sd) / sd_ref
)
cat("\n== Subject log tau_i and beta_i against Stan, range over the 18 ==\n")
print(round(
  rbind(
    log_tau = c(
      max_mean_gap = max(gap["mean", 1:n_sub]),
      range(gap["sd_ratio", 1:n_sub])
    ),
    beta = c(
      max_mean_gap = max(gap["mean", n_sub + 1:n_sub]),
      range(gap["sd_ratio", n_sub + 1:n_sub])
    )
  ),
  3
))

## ----- Cost and checks -------------------------------------------------------

cat("\n== Cost ==\n")
cat("inla_vanilla():", fit$counts[["logpost"]], "calls,", round(t_fit, 2))
cat(" s, then", round(t_draws, 2), "s for the subject draws\n")
cat("Stan:", round(t_stan, 1), "s of sampling\n")
cat(
  "Stan: max R-hat",
  round(max(rstan::summary(stan_fit)$summary[, "Rhat"]), 4)
)
cat(", divergences", rstan::get_num_divergent(stan_fit), "\n")
cat(
  "Quadrature: logpost at the mode, 41 nodes minus 15:",
  signif(quad_error, 2)
)
cat("\n\n")
summary(fit)
