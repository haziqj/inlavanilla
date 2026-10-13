# Skew-normal marginals against full Laplace marginals on the Political
# Democracy SEM (39 parameters) of vignettes/articles/rcpp.qmd.
#
# Run the article's chunks first (interactively, in order), so that `Y`,
# `init`, `logpost_cpp()`, `grad_cpp()` and the fit with the gradient, `fit`,
# exist. The full Laplace marginal of parameter j maximises the log-posterior
# over the other parameters at each point of a grid in theta_j, and corrects
# by half the log-determinant of their Hessian there.
#
# The grid must be wide: for the two log residual variances below, the left
# tail reaches 12 Laplace SDs from the mode. An earlier run on +-4.5 SDs cut
# it off and understated the error (about 25% for logv_dem65).
#
# Result (October 2026), which agrees with HMC (lab/sem-hmc.R):
#               mean    sd    2.5%    50%   97.5%
# logv_x2     full Laplace -2.623 1.412 -6.707 -2.225 -1.257
#             HMC          -2.659 1.470 -6.984 -2.237 -1.254
#             skew-normal  -2.165 0.566 -3.461 -2.086 -1.271
# logv_dem65  full Laplace -2.243 1.839 -7.340 -1.703 -0.275
#             HMC          -2.251 1.864 -7.341 -1.707 -0.258
#             skew-normal  -1.598 0.861 -3.598 -1.460 -0.302
# The left tail of the log-density has slope 0.50 in both: the
# gamma(1, 0.5) prior on the SD with its Jacobian.

full_laplace <- function(j, lo, hi, ngrid = 81) {
  grid <- seq(lo, hi, length.out = ngrid)
  out <- rep(NA_real_, ngrid)
  start <- fit$mode[-j]
  for (k in seq_along(grid)) {
    f <- function(x) {
      v <- -logpost_cpp(append(x, grid[k], after = j - 1), Y)
      if (is.finite(v)) v else 1e10
    }
    g <- function(x) -grad_cpp(append(x, grid[k], after = j - 1), Y)[-j]
    o <- nlminb(
      start,
      f,
      g,
      control = list(iter.max = 1000, eval.max = 2000, rel.tol = 1e-12)
    )
    start <- o$par # warm start along the grid
    H <- sapply(seq_along(o$par), function(i) {
      e <- replace(numeric(length(o$par)), i, 1e-5)
      (g(o$par + e) - g(o$par - e)) / 2e-5
    })
    d <- determinant((H + t(H)) / 2, logarithm = TRUE)
    out[k] <- if (d$sign > 0) -o$objective - 0.5 * d$modulus else NA
  }
  ok <- is.finite(out)
  fine <- seq(lo, hi, length.out = 4001)
  p <- exp(stats::spline(grid[ok], out[ok] - max(out[ok]), xout = fine)$y)
  p <- p / sum(p)
  m <- sum(fine * p)
  s <- sqrt(sum((fine - m)^2 * p))
  q <- stats::approx(cumsum(p), fine, c(0.025, 0.5, 0.975), ties = mean)$y
  c(mean = m, sd = s, q025 = q[1], q50 = q[2], q975 = q[3])
}

for (p in c("logv_x2", "logv_dem65")) {
  j <- match(p, names(init))
  fl <- full_laplace(j, -16, fit$mode[[j]] + 6 * sqrt(fit$cov[j, j]))
  sn <- fit$summary[p, c("Mean", "SD", "2.5%", "50%", "97.5%")]
  print(round(rbind(full_laplace = fl, skew_normal = sn), 3))
}
