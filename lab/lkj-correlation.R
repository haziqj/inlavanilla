################################################################################
#
# Example: a correlation matrix with an LKJ prior. Four variables are drawn
# from a multivariate normal with means mu, SDs sigma and correlation matrix
# Omega, and Omega has the LKJ(eta) prior of Lewandowski, Kurowicka and Joe
# (2009). inla_vanilla() works on an unconstrained parameter vector, so Omega
# enters through its six canonical partial correlations (CPCs) z, on the atanh
# scale, as in Stan's cholesky_factor_corr type: the CPCs fill the Cholesky
# factor L of Omega row by row, and logpost() adds the log-Jacobian of the map
# from atanh(z) to L. sampling() then turns each draw back into Omega. If rstan
# is installed, the same model is fitted in Stan as a check. Step through it in
# an R session.
#
################################################################################

## ----- Configuration ---------------------------------------------------------

library(inlavanilla)
set.seed(1)
n <- 50 # observations
eta <- 2 # LKJ shape: eta > 1 shrinks Omega towards the identity
nsamp <- 4000 # draws from the fitted posterior
run_stan <- requireNamespace("rstan", quietly = TRUE)

mu_true <- c(0, 1, -1, 2)
sigma_true <- c(1, 2, 0.5, 1)
# fmt: skip
omega_true <- matrix(c(
   1.0, 0.6, 0.3, -0.2,
   0.6, 1.0, 0.5,  0.0,
   0.3, 0.5, 1.0,  0.4,
  -0.2, 0.0, 0.4,  1.0
), 4, 4)

## =============================================================================
## Model
## =============================================================================

## ----- From unconstrained values to a correlation matrix ---------------------

# Builds the Cholesky factor L of a d x d correlation matrix from a vector w of
# d(d - 1)/2 unconstrained values, as Stan does. Each w gives a CPC
# z = tanh(w) in (-1, 1). Row i of L has unit length, and its entries take,
# from left to right, the share z of the squared length still left, so what is
# left is a product of the (1 - z^2). Keeping it on the log scale keeps the
# diagonal and the Jacobian exact when |w| is large and tanh(w) rounds to 1.
cholesky_corr <- function(w, d) {
  z <- tanh(w)
  log_1mz2 <- 2 * (log(2) - abs(w) - log1p(exp(-2 * abs(w)))) # log(1 - z^2)
  L <- diag(d)
  log_diag <- numeric(d) # log L[i, i]
  log_jac <- sum(log_1mz2) # dz / dw
  k <- 0
  for (i in seq_len(d)[-1]) {
    log_left <- 0 # log of the squared length left in row i
    for (j in seq_len(i - 1)) {
      k <- k + 1
      L[i, j] <- z[k] * exp(log_left / 2)
      log_jac <- log_jac + log_left / 2 # dL[i, j] / dz
      log_left <- log_left + log_1mz2[k]
    }
    L[i, i] <- exp(log_left / 2)
    log_diag[i] <- log_left / 2
  }
  list(L = L, log_diag = log_diag, log_jac = log_jac)
}

# The row and column of L that each CPC fills, in the order of w
cpc_index <- function(d) {
  data.frame(
    row = rep(seq_len(d)[-1], seq_len(d - 1)),
    col = sequence(seq_len(d - 1))
  )
}

## ----- LKJ prior -------------------------------------------------------------

# The LKJ(eta) density of Omega is proportional to det(Omega)^(eta - 1). As a
# density of L it gains the Jacobian of L -> Omega and becomes
# prod_i L[i, i]^(d - i + 2 eta - 2), Stan's lkj_corr_cholesky. The normalising
# constant comes from the CPCs, which under LKJ(eta) are independent, with
# (z + 1) / 2 ~ Beta(b, b) and b = eta + (d - 1 - j) / 2 in column j of L
# (Lewandowski et al., 2009). With it, fit$logml is a proper log marginal
# likelihood.
lkj_cholesky_lpdf <- function(log_diag, eta) {
  d <- length(log_diag)
  b <- eta + (d - 1 - cpc_index(d)$col) / 2
  log_const <- -sum(lbeta(b, b) + (2 * b - 1) * log(2))
  sum((d - seq_len(d) + 2 * eta - 2) * log_diag) + log_const
}

## ----- Log-posterior ---------------------------------------------------------

# Splits theta into the means, the SDs and the Cholesky factor of Omega
unpack <- function(theta, d) {
  log_sigma <- theta[d + seq_len(d)]
  c(
    list(mu = theta[seq_len(d)], log_sigma = log_sigma, sigma = exp(log_sigma)),
    cholesky_corr(theta[-seq_len(2 * d)], d)
  )
}

# Each row is x_i = mu + diag(sigma) L e_i with e_i ~ N(0, I), so the
# log-likelihood is that of the e_i, less n log det(diag(sigma) L). Priors:
# mu ~ N(0, 10^2) and sigma ~ Exp(1), with the Jacobian of log sigma, and
# Omega ~ LKJ(eta), with the Jacobian of atanh(z) -> L.
logpost <- function(theta, x, eta) {
  d <- ncol(x)
  p <- unpack(theta, d)
  e <- forwardsolve(p$L, (t(x) - p$mu) / p$sigma)
  loglik <- sum(dnorm(e, log = TRUE)) -
    nrow(x) * (sum(p$log_sigma) + sum(p$log_diag))
  logprior <- sum(dnorm(p$mu, 0, 10, log = TRUE)) +
    sum(dexp(p$sigma, rate = 1, log = TRUE)) +
    sum(p$log_sigma) +
    lkj_cholesky_lpdf(p$log_diag, eta) +
    p$log_jac
  loglik + logprior
}

## ----- Check the Jacobian ----------------------------------------------------

# Under LKJ(eta) each w = atanh(z) has the density
# (1 - z^2)^b / (2^(2b - 1) B(b, b)), independently, so the LKJ density of L
# times the Jacobian of w -> L must equal their product at any w. The prior
# could be written that way directly. The LKJ form is the one that carries
# over to a prior on Omega that the CPCs do not factorise.
d <- ncol(omega_true)
w <- rnorm(choose(d, 2), sd = 2)
chol_w <- cholesky_corr(w, d)
b <- eta + (d - 1 - cpc_index(d)$col) / 2
stopifnot(all.equal(
  lkj_cholesky_lpdf(chol_w$log_diag, eta) + chol_w$log_jac,
  sum(b * log1p(-tanh(w)^2) - lbeta(b, b) - (2 * b - 1) * log(2))
))

## =============================================================================
## Fit
## =============================================================================

## ----- Simulate data ---------------------------------------------------------

sigma_mat <- diag(sigma_true) %*% omega_true %*% diag(sigma_true)
x <- matrix(rnorm(n * d), n, d) %*% chol(sigma_mat) + rep(mu_true, each = n)

## ----- Fit on the atanh scale ------------------------------------------------

# The means and SDs start at their sample values, and the CPCs at zero, the
# identity
cpc <- cpc_index(d)
init <- c(
  setNames(colMeans(x), paste0("mu", seq_len(d))),
  setNames(log(apply(x, 2, sd)), paste0("log_sigma", seq_len(d))),
  setNames(numeric(nrow(cpc)), paste0("atanh_z", cpc$row, cpc$col))
)
fit <- inla_vanilla(logpost, init, x = x, eta = eta)
summary(fit)

## ----- Back to the correlations ----------------------------------------------

# Each draw on the natural scale: the means, the SDs, and the correlations of
# Omega = L L' named after their row and column
natural <- function(theta, d) {
  p <- unpack(theta, d)
  omega <- tcrossprod(p$L)
  low <- lower.tri(omega)
  r <- omega[low]
  names(r) <- paste0("r", col(omega)[low], row(omega)[low])
  c(p$mu, setNames(p$sigma, paste0("sigma", seq_len(d))), r)
}
draws <- sampling(fit, nsamp = nsamp, fun = natural, d = d)

post <- data.frame(
  truth = c(mu_true, sigma_true, omega_true[lower.tri(omega_true)]),
  mean = colMeans(draws),
  sd = apply(draws, 2, sd),
  t(apply(draws, 2, quantile, c(0.025, 0.5, 0.975))),
  check.names = FALSE
)
round(post, 3)

# The posterior mean of Omega next to the sample correlations
r_cols <- grep("^r", colnames(draws))
omega_mean <- diag(d)
omega_mean[lower.tri(omega_mean)] <- colMeans(draws[, r_cols])
omega_mean[upper.tri(omega_mean)] <- t(omega_mean)[upper.tri(omega_mean)]
round(omega_mean, 3)
round(cor(x), 3)

## =============================================================================
## Check against Stan
## =============================================================================

if (run_stan) {
  lkj_stan <- "
  data {
    int<lower=1> N;
    int<lower=1> D;
    array[N] vector[D] x;
    real<lower=0> eta;
  }
  parameters {
    vector[D] mu;
    vector<lower=0>[D] sigma;
    cholesky_factor_corr[D] L;
  }
  model {
    mu ~ normal(0, 10);
    sigma ~ exponential(1);
    L ~ lkj_corr_cholesky(eta);
    x ~ multi_normal_cholesky(mu, diag_pre_multiply(sigma, L));
  }
  generated quantities {
    corr_matrix[D] Omega = multiply_lower_tri_self_transpose(L);
  }
  "
  stan_fit <- rstan::stan(
    model_code = lkj_stan,
    data = list(N = n, D = d, x = x, eta = eta),
    chains = 4,
    iter = 6000,
    warmup = 1000,
    seed = 1,
    refresh = 0
  )
  low <- lower.tri(diag(d))
  stan_pars <- c(
    sprintf("mu[%d]", seq_len(d)),
    sprintf("sigma[%d]", seq_len(d)),
    sprintf("Omega[%d,%d]", row(diag(d))[low], col(diag(d))[low])
  )
  stan_draws <- as.matrix(stan_fit, pars = c("mu", "sigma", "Omega"))
  stan_draws <- stan_draws[, stan_pars]
  colnames(stan_draws) <- colnames(draws)

  # Gaps in posterior SDs of the Stan draws: under 0.1 is a close match
  compare <- data.frame(
    mean = post$mean,
    stan_mean = colMeans(stan_draws),
    sd = post$sd,
    stan_sd = apply(stan_draws, 2, sd),
    row.names = colnames(draws)
  )
  compare$mean_gap <- (compare$mean - compare$stan_mean) / compare$stan_sd
  compare$sd_ratio <- compare$sd / compare$stan_sd
  print(round(compare, 3))
}

## =============================================================================
## Plot the correlations
## =============================================================================

# The posterior of each correlation from the draws, with Stan's dashed when it
# ran and the true value in grey
op <- par(mfrow = c(2, 3), mar = c(4, 3, 2, 1))
for (j in r_cols) {
  dens <- density(draws[, j], from = -1, to = 1)
  stan_dens <- if (run_stan) density(stan_draws[, j], from = -1, to = 1)
  plot(
    dens,
    main = colnames(draws)[j],
    xlab = "",
    ylim = c(0, max(dens$y, stan_dens$y))
  )
  if (run_stan) {
    lines(stan_dens, lty = 2)
  }
  abline(v = post$truth[j], col = "grey")
}
if (run_stan) {
  legend("topleft", c("inla_vanilla", "Stan"), lty = 1:2, bty = "n")
}
par(op)
