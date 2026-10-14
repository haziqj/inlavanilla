################################################################################
#
# A stochastic production frontier with exponential inefficiency (October
# 2026): a model R-INLA cannot fit. On the R-INLA discussion group (7-8
# October 2020, "Stochastic production frontiers",
# https://groups.google.com/g/r-inla-discussion-group/c/7rd0KAncAlo) a user
# asked to move this model from Stan to INLA, and the answer was that INLA
# "cannot do an exponential random effect like this".
#
# The model, with the poster's priors and the sign of u corrected as Finn
# Lindgren read it in the thread:
#   log y_i = x_i' beta - u_i + v_i, v_i ~ N(0, sigma^2), u_i ~ Exp(lambda),
#   beta ~ N(0, 100^2), sigma ~ InvGamma(0.001, 0.001),
#   lambda ~ Exp(-log(0.875)),
# the last the prior of van den Broeck et al. (1994), with prior median
# efficiency 0.875. Each u_i integrates out in closed form (Aigner, Lovell and
# Schmidt, 1977): with e_i = log y_i - x_i' beta,
#   p(e_i) = lambda exp(lambda e_i + lambda^2 sigma^2 / 2)
#            Phi(-e_i / sigma - lambda sigma),
# so inla_vanilla() fits 6 parameters: beta (4) and log sigma and log lambda.
# Given e_i, u_i is N(-e_i - lambda sigma^2, sigma^2) truncated at 0, so the
# efficiencies exp(-u_i) come from sampling(), a draw of u_i at each draw of
# the parameters.
#
# Data: the Cobb-Douglas frontier of 344 rice harvests (43 Philippine farms,
# 8 years) in frontier::riceProdPhil (Coelli et al., 2005), pooled as in the
# thread, one u_i per observation. Reference: the poster's own Stan model,
# with the 344 u_i sampled.
#
# Needs frontier and rstan. Run from the package root:
#   Rscript lab/frontier-exponential.R
#
# Results (October 2026; Stan 4 chains of 3,000 kept draws, R-hat 1.003, no
# divergences): inla_vanilla() matches Stan. The six parameters are within
# 0.03 posterior SDs in mean and 0.07 in the 2.5% and 97.5% quantiles, with
# SDs within 2.5%. The efficiencies of the 344 harvests are within 0.006 of
# Stan's posterior means and 0.02 of its 95% limits, and the mean efficiency
# is 0.789 in both. The fit took 0.2 s (776 calls of logpost, no gradient)
# and the efficiency draws 0.4 s, against 21 s of sampling for Stan, after a
# minute of compiling. Every check in summary() passes.
#
################################################################################

## ----- Configuration ---------------------------------------------------------

smoke <- FALSE
suppressMessages(devtools::load_all(".", quiet = TRUE))
data(riceProdPhil, package = "frontier")
n_iter <- if (smoke) 1000 else 6000 # Stan iterations per chain, half warm-up
n_draws <- if (smoke) 2000 else 10000 # copula draws
options(mc.cores = 4)

## =============================================================================
## Model and data
## =============================================================================

y <- log(riceProdPhil$PROD)
X <- cbind(
  b0 = 1,
  b_area = log(riceProdPhil$AREA),
  b_labor = log(riceProdPhil$LABOR),
  b_npk = log(riceProdPhil$NPK)
)
n <- length(y)

## ----- Log-posterior ---------------------------------------------------------

# The marginal log-likelihood with each u_i integrated out, and the priors on
# the unconstrained scale with their Jacobians (log sigma and log lambda)
logpost <- function(theta, y, X) {
  beta <- theta[1:4]
  sigma <- exp(theta[["log_sigma"]])
  lambda <- exp(theta[["log_lambda"]])
  e <- y - drop(X %*% beta)
  loglik <- sum(
    log(lambda) +
      lambda * e +
      lambda^2 * sigma^2 / 2 +
      pnorm(-e / sigma - lambda * sigma, log.p = TRUE)
  )
  log_prior_sigma <- 0.001 *
    log(0.001) -
    lgamma(0.001) -
    1.001 * log(sigma) -
    0.001 / sigma +
    log(sigma)
  log_prior_lambda <- dexp(lambda, -log(0.875), log = TRUE) + log(lambda)
  loglik +
    sum(dnorm(beta, 0, 100, log = TRUE)) +
    log_prior_sigma +
    log_prior_lambda
}

ols <- lm.fit(X, y)
init <- c(
  ols$coefficients,
  log_sigma = log(sd(ols$residuals) / 2),
  log_lambda = log(3)
)

## =============================================================================
## Fits
## =============================================================================

## ----- inla_vanilla() --------------------------------------------------------

t_fit <- system.time(
  fit <- inla_vanilla(logpost, init, y = y, X = X)
)[["elapsed"]]

# Efficiencies: at each draw of the parameters, u_i from its truncated normal,
# drawn through the upper tail on the log scale for stability
efficiency <- function(theta) {
  sigma <- exp(theta[["log_sigma"]])
  lambda <- exp(theta[["log_lambda"]])
  e <- y - drop(X %*% theta[1:4])
  mu <- -e - lambda * sigma^2
  log_tail <- log(runif(n)) + pnorm(mu / sigma, log.p = TRUE)
  u <- mu + sigma * qnorm(log_tail, lower.tail = FALSE, log.p = TRUE)
  exp(-u)
}
set.seed(1)
t_eff <- system.time(
  eff_iv <- sampling(fit, nsamp = n_draws, fun = efficiency)
)[["elapsed"]]

## ----- Stan, the poster's model ----------------------------------------------

stan_code <- "
data {
  int<lower=1> N;
  int<lower=1> K;
  matrix[N, K] X;
  vector[N] y;
}
parameters {
  vector[K] beta;
  real<lower=0> sigma;
  real<lower=0> lambda;
  vector<lower=0>[N] u;
}
model {
  beta ~ normal(0, 100);
  sigma ~ inv_gamma(0.001, 0.001);
  lambda ~ exponential(-log(0.875));
  u ~ exponential(lambda);
  y ~ normal(X * beta - u, sigma);
}
generated quantities {
  real log_sigma = log(sigma);
  real log_lambda = log(lambda);
  vector[N] eff = exp(-u);
}
"
stan_mod <- rstan::stan_model(model_code = stan_code)
t_stan <- system.time(
  stan_fit <- rstan::sampling(
    stan_mod,
    data = list(N = n, K = ncol(X), X = X, y = y),
    chains = 4,
    iter = n_iter,
    seed = 2026,
    refresh = 0,
    control = list(adapt_delta = 0.95)
  )
)[["elapsed"]]
stan_draws <- as.matrix(stan_fit)
colnames(stan_draws)[1:4] <- colnames(X)

## =============================================================================
## Comparison
## =============================================================================

## ----- Parameters ------------------------------------------------------------

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

## ----- Efficiencies ----------------------------------------------------------

eff_stan <- stan_draws[, grep("^eff\\[", colnames(stan_draws))]
eff_summary <- function(draws) {
  cbind(
    mean = colMeans(draws),
    q025 = apply(draws, 2, quantile, 0.025),
    q975 = apply(draws, 2, quantile, 0.975)
  )
}
e_iv <- eff_summary(eff_iv)
e_stan <- eff_summary(eff_stan)
cat("\n== Efficiencies of the 344 harvests: largest gap to Stan ==\n")
print(round(apply(abs(e_iv - e_stan), 2, max), 4))
cat(
  "Mean efficiency: inla_vanilla()",
  round(mean(e_iv[, "mean"]), 4),
  "| Stan",
  round(mean(e_stan[, "mean"]), 4),
  "\n"
)

## ----- Cost and checks -------------------------------------------------------

cat("\n== Cost ==\n")
cat("inla_vanilla():", fit$counts[["logpost"]], "calls,", round(t_fit, 2))
cat(" s, then", round(t_eff, 2), "s for the efficiency draws\n")
cat("Stan:", round(t_stan, 1), "s of sampling\n")
cat(
  "Stan: max R-hat",
  round(max(rstan::summary(stan_fit)$summary[, "Rhat"]), 4)
)
cat(", divergences", rstan::get_num_divergent(stan_fit), "\n\n")
summary(fit)
