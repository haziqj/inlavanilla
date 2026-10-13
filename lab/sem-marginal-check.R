# Skew-normal marginals against full Laplace marginals on the Political
# Democracy SEM (39 parameters) of vignettes/articles/rcpp.qmd.
#
# Run the article's chunks first (interactively, in order), so that `Y`,
# `init`, `logpost_cpp()`, `grad_cpp()` and the fit with the gradient, `fit`,
# exist. The full Laplace marginal of parameter j maximises the log-posterior
# over the other parameters at each point of a grid in theta_j, and corrects
# by half the log-determinant of their Hessian there. It is about 10 times the
# cost of the skew-normal stage, and is the reference for the straight-scan
# limitation recorded in lab/DESIGN.md.
#
# Result (October 2026): the skew-normal SD is too small by about 25% for
# logv_dem65 (bend 0.44), 19% for logv_x2 (bend 0.11), 6% for logv_x1 (bend
# 0.14), 5.5% for b (bend 0.11) and 2.6% for the intercepts (bend about 0.05,
# not flagged). For logv_dem65 the 2.5% quantile is -3.60 against -4.65.

full_laplace <- function(j, ngrid = 31) {
  sdj <- sqrt(fit$cov[j, j])
  grid <- fit$mode[j] + seq(-4.5, 4.5, length.out = ngrid) * sdj
  out <- numeric(ngrid)
  for (k in seq_along(grid)) {
    f <- function(x) -logpost_cpp(append(x, grid[k], after = j - 1), Y)
    g <- function(x) -grad_cpp(append(x, grid[k], after = j - 1), Y)[-j]
    start <- fit$mode[-j] +
      fit$cov[-j, j] / fit$cov[j, j] * (grid[k] - fit$mode[j])
    o <- nlminb(start, f, g, control = list(iter.max = 500, rel.tol = 1e-12))
    H <- sapply(seq_along(o$par), function(i) {
      h <- 1e-5
      e <- replace(numeric(length(o$par)), i, h)
      (g(o$par + e) - g(o$par - e)) / (2 * h)
    })
    H <- (H + t(H)) / 2
    out[k] <- -o$objective - 0.5 * determinant(H, logarithm = TRUE)$modulus
  }
  fine <- seq(min(grid), max(grid), length.out = 2001)
  p <- exp(stats::spline(grid, out - max(out), xout = fine)$y)
  p <- p / sum(p)
  m <- sum(fine * p)
  s <- sqrt(sum((fine - m)^2 * p))
  q <- stats::approx(cumsum(p), fine, c(0.025, 0.5, 0.975), ties = mean)$y
  c(mean = m, sd = s, q = q)
}

parnames <- names(init)
pars <- c(
  "logv_dem65",
  "logv_x2",
  "logv_x1",
  "b",
  "z_y4_y8",
  "beta_dem65_ind60",
  "lam_x2",
  "nu_y1"
)
res <- t(sapply(pars, function(p) {
  fl <- full_laplace(match(p, parnames))
  sn <- fit$summary[p, c("Mean", "SD", "2.5%", "50%", "97.5%")]
  c(
    sd_ratio = sn[[2]] / fl[[2]],
    mean_diff_sd = (sn[[1]] - fl[[1]]) / fl[[2]],
    fl_q025 = fl[[3]],
    sn_q025 = sn[[3]],
    fl_q975 = fl[[5]],
    sn_q975 = sn[[5]]
  )
}))
print(round(res, 3))
print(fit$diagnostics$marginals[pars, c("bend", "nmad", "tails")])
