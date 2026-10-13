# Stage 2: skew-normal marginals, one per parameter, ported from INLAvaan.
#
# For parameter j the joint log-posterior is scanned along the direction of
# the conditional mode path (CMP) of the Laplace Gaussian,
#   theta(t) = mode + t v_j,  v_j = Sigma[, j] / sqrt(Sigma[j, j]),
# so theta_j(t) = mode_j + t sd_j and the other parameters sit at their
# Gaussian conditional means (Lemma 3.1). The log-heights are corrected for
# the change in the volume of the other parameters along the path by adding
# t gamma1_j (Corollary 3.3), and a skew-normal is fitted to the corrected
# log-heights on a grid of t.

# The volume correction gamma1_j: the rate of change, along v_j, of minus half
# the log-determinant of the conditional precision of the other parameters,
# taken by a forward difference of `d_outer` posterior SDs from the mode. In
# whitened coordinates (columns of L) the precision at the mode is the
# identity, so its trace is m and v_j' H v_j = 1 there.

# INLAvaan's "shortcut" correction, from the gradient: the whitened trace and
# the curvature along v_j at mode + d_outer v_j by central differences of the
# gradient, 2m + 2 gradient calls.
gamma1_grad <- function(j, gr, mode, Vscan, L, d_outer = 0.01, h = 1e-5) {
  m <- length(mode)
  vj <- Vscan[, j]
  th <- mode + vj * d_outer
  tr <- 0
  for (k in seq_len(m)) {
    Lk <- L[, k]
    tr <- tr + sum(Lk * (gr(th - Lk * h) - gr(th + Lk * h))) / (2 * h)
  }
  vHv <- sum(vj * (gr(th - vj * h) - gr(th + vj * h))) / (2 * h)
  -0.5 * (tr - m) / d_outer + 0.5 * (vHv - 1) / d_outer
}

# The same correction from the log-posterior alone. Each whitened curvature is
# a second difference of the log-posterior rather than a difference of two
# finite-difference gradients, so the cost is 2m + 5 calls per parameter
# instead of the 4m^2 + 4m that the gradient route would take with
# finite-difference gradients. Second differences carry a truncation error of
# their own, so the baselines at the mode (`tr0`, and the curvature along v_j
# computed here) are taken by the same rule rather than set to m and 1, and
# the error cancels in the forward difference.
gamma1_lp <- function(
  j,
  lp,
  mode,
  lp_max,
  Vscan,
  L,
  tr0,
  d_outer = 0.01,
  h = 1e-2
) {
  m <- length(mode)
  vj <- Vscan[, j]
  th <- mode + vj * d_outer
  f_th <- lp(th)
  tr <- 0
  for (k in seq_len(m)) {
    tr <- tr - fd_curvature(lp, th, L[, k], h, f0 = f_th)
  }
  vHv <- -fd_curvature(lp, th, vj, h, f0 = f_th)
  vHv0 <- -fd_curvature(lp, mode, vj, h, f0 = lp_max)
  -0.5 * (tr - tr0) / d_outer + 0.5 * (vHv - vHv0) / d_outer
}

# Whitened trace of the precision at the mode by the rule of gamma1_lp().
whitened_trace_lp <- function(lp, mode, lp_max, L, h = 1e-2) {
  tr <- 0
  for (k in seq_len(ncol(L))) {
    tr <- tr - fd_curvature(lp, mode, L[, k], h, f0 = lp_max)
  }
  tr
}

# Scan each parameter over +-4 Laplace SDs, as in INLAvaan, and fit its
# skew-normal. `tails` is FALSE when the corrected log-height at an end of the
# scan is still within 3 of its maximum (over 5% of the peak density): the
# Laplace SD has badly understated a heavy tail, and the skew-normal
# extrapolates it. Moderately skewed marginals, which the skew-normal fits
# well, end around -4 to -6 on their heavy side. Extending the scan was tried
# and rejected: it improved some heavy tails but pulled the fit away from the
# bulk of others (a Fisher-z correlation lost an order of magnitude in SD
# accuracy), so the problem is reported instead.
fit_marginals <- function(lp, gr, lap, ctrl, tick = function() NULL) {
  mode <- lap$mode
  Sigma <- lap$Sigma
  L <- lap$L
  m <- length(mode)
  sds <- sqrt(diag(Sigma))
  Vscan <- sweep(Sigma, 2, sds, "/")
  z <- seq(-4, 4, length.out = ctrl$n_grid)

  tr0 <- if (is.null(gr)) {
    whitened_trace_lp(lp, mode, lap$lp_max, L)
  } else {
    NULL
  }

  out <- matrix(
    NA_real_,
    m,
    7,
    dimnames = list(
      names(mode),
      c("xi", "omega", "alpha", "gamma1", "nmad", "tails", "bend")
    )
  )
  # Scan points for the bend check (see scan_bend())
  i_bend <- unique(c(which.min(abs(z + 2)), which.min(abs(z - 2))))

  for (j in seq_len(m)) {
    g1 <- if (is.null(gr)) {
      gamma1_lp(j, lp, mode, lap$lp_max, Vscan, L, tr0)
    } else {
      gamma1_grad(j, gr, mode, Vscan, L)
    }
    # Log-posterior along the scan; the mode's height is known
    raw <- vapply(
      z,
      function(t) {
        if (t == 0) lap$lp_max else lp(mode + Vscan[, j] * t)
      },
      numeric(1)
    )
    above <- function(h, y) {
      is.finite(h) && h - max(y[is.finite(y)]) > -3
    }

    # Corrected log-heights
    y <- raw + g1 * z
    tails <- !above(y[1], y) && !above(y[length(y)], y)
    bend <- scan_bend(j, lp, gr, lap, Vscan, z[i_bend], raw[i_bend])

    ok <- is.finite(y)
    fit <- fit_skew_normal(
      x = z[ok],
      y = y[ok] - max(y[ok]),
      threshold_log_drop = ctrl$log_threshold
    )

    # Back from the scan coordinate t to theta_j = mode_j + t sd_j
    out[j, ] <- c(
      mode[[j]] + fit$xi * sds[[j]],
      fit$omega * sds[[j]],
      fit$alpha,
      g1,
      fit$nmad,
      tails,
      bend
    )
    tick()
  }
  out
}

# The mass each fitted marginal puts outside the window that was scanned to
# fit it, mode +- 4 Laplace SDs: there the skew-normal is extrapolated, so
# when this is large the outer quantiles are the ones to distrust. This is
# INLAvaan's `scan_end_mass` (in its diagnostics()). A Gaussian marginal
# gives 2 pnorm(-4) = 6.3e-5; larger values mean more of the marginal is
# extrapolated. summary() flags it, and the fit warns, above 0.05.
scan_end_mass <- function(sn, mode, sds) {
  lo <- psnorm(mode - 4 * sds, sn[, "xi"], sn[, "omega"], sn[, "alpha"])
  hi <- psnorm(
    mode + 4 * sds,
    sn[, "xi"],
    sn[, "omega"],
    sn[, "alpha"],
    lower_tail = FALSE
  )
  unname(lo + hi)
}

# How far the straight scan path misses the conditional mode of the other
# parameters. At scan points t (near +-2 posterior SDs), one Newton step in
# the other parameters, with the Laplace curvature, would raise the
# log-posterior by 0.5 g' C g, where g is the gradient in the other parameters
# and C their conditional covariance given theta_j,
#   C = Sigma[-j, -j] - Sigma[-j, j] Sigma[j, -j] / Sigma[j, j].
# That is the height the scan misses there, in log-density units. When it is
# large, the other parameters' conditional mode bends away from the line (as
# a noise scale does for a regression coefficient at small n), and the
# marginal comes out too narrow. Returns the larger of the gains. The
# gradient is the user's, or one-sided differences in the other coordinates
# only (m - 1 calls per point, reusing the scan's height), with their bias
# removed by the Laplace curvature as in the VB stage.
scan_bend <- function(j, lp, gr, lap, Vscan, t, f_t) {
  m <- length(lap$mode)
  if (m < 2L) {
    return(0)
  }
  Sigma <- lap$Sigma
  sds <- sqrt(diag(Sigma))
  C <- Sigma[-j, -j, drop = FALSE] -
    tcrossprod(Sigma[-j, j]) / Sigma[j, j]
  others <- seq_len(m)[-j]
  h <- 1e-3 * sds
  gain <- numeric(length(t))
  for (k in seq_along(t)) {
    x <- lap$mode + Vscan[, j] * t[k]
    g <- if (is.null(gr)) {
      vapply(
        others,
        function(i) {
          xp <- x
          xp[i] <- x[i] + h[i]
          (lp(xp) - f_t[k]) / h[i] + 0.5 * h[i] * lap$H[i, i]
        },
        numeric(1)
      )
    } else {
      gr(x)[others]
    }
    gain[k] <- if (all(is.finite(g))) 0.5 * sum(g * (C %*% g)) else Inf
  }
  max(gain)
}
