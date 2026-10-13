# HMC reference for the Political Democracy SEM of
# vignettes/articles/rcpp.qmd: the check behind the heavy-tail limitation in
# lab/DESIGN.md.
#
# Run the article's chunks first (interactively, in order), so that `Y`,
# `init`, `logpost_cpp()`, `grad_cpp()` and the fit with the gradient, `fit`,
# exist. Needs the posterior package for the summaries.
#
# Plain HMC in the whitened coordinates of the Laplace approximation
# (identity mass matrix there), with a jittered step size adapted in warm-up
# to an acceptance of 0.8 and 1 to 40 leapfrog steps. 4 chains of 2,000
# warm-up and 10,000 kept draws took about 30 s. Result (October 2026):
# R-hat 1.000, bulk ESS over 11,000, tail ESS over 8,000, 1 divergence.

L <- t(chol(fit$cov))
centre <- fit$mode
m <- length(init)
lp_u <- function(u) {
  v <- logpost_cpp(centre + drop(L %*% u), Y)
  if (is.finite(v)) v else -Inf
}
gr_u <- function(u) drop(crossprod(L, grad_cpp(centre + drop(L %*% u), Y)))

hmc_chain <- function(seed, n_warm = 2000, n_keep = 10000, l_max = 40) {
  set.seed(seed)
  u <- rnorm(m) * 0.5
  lp <- lp_u(u)
  g <- gr_u(u)
  log_eps <- log(0.2)
  out <- matrix(NA_real_, n_keep, m, dimnames = list(NULL, names(init)))
  for (it in seq_len(n_warm + n_keep)) {
    e <- exp(log_eps) * runif(1, 0.8, 1.2)
    p <- rnorm(m)
    h0 <- lp - 0.5 * sum(p^2)
    u1 <- u
    g1 <- g
    p1 <- p + 0.5 * e * g1
    ok <- TRUE
    n_steps <- sample.int(l_max, 1)
    for (s in seq_len(n_steps)) {
      u1 <- u1 + e * p1
      g1 <- gr_u(u1)
      if (!all(is.finite(g1))) {
        ok <- FALSE
        break
      }
      if (s < n_steps) p1 <- p1 + e * g1
    }
    a <- 0
    if (ok) {
      p1 <- p1 + 0.5 * e * g1
      lp1 <- lp_u(u1)
      h1 <- lp1 - 0.5 * sum(p1^2)
      if (is.finite(h1)) a <- min(1, exp(h1 - h0))
    }
    if (runif(1) < a) {
      u <- u1
      lp <- lp1
      g <- g1
    }
    if (it <= n_warm) {
      log_eps <- log_eps + (a - 0.8) / sqrt(it + 10)
    } else {
      out[it - n_warm, ] <- centre + drop(L %*% u)
    }
  }
  out
}

chains <- lapply(1:4, hmc_chain)
draws <- posterior::as_draws_array(
  aperm(simplify2array(chains), c(1, 3, 2))
)
sm <- posterior::summarise_draws(
  draws,
  "mean",
  "sd",
  ~ quantile(.x, c(0.025, 0.975)),
  "rhat",
  "ess_bulk",
  "ess_tail"
)
sn <- fit$summary[sm$variable, ]
cmp <- data.frame(
  par = sm$variable,
  sd_ratio = sn[, "SD"] / sm$sd,
  mean_err = (sn[, "Mean"] - sm$mean) / sm$sd,
  q025_err = (sn[, "2.5%"] - sm$`2.5%`) / sm$sd,
  q975_err = (sn[, "97.5%"] - sm$`97.5%`) / sm$sd,
  bend = fit$diagnostics$marginals[sm$variable, "bend"]
)
print(cmp[order(cmp$sd_ratio), ], digits = 3, row.names = FALSE)
