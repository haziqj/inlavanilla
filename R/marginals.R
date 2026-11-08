# Stage 2: skew-normal marginals, one per parameter, ported from INLAvaan.
#
# For parameter j the joint log-posterior is scanned along the direction of
# the conditional mode path (CMP) of the Laplace Gaussian,
#   theta(t) = mode + t v_j,  v_j = Sigma[, j] / sqrt(Sigma[j, j]),
# so theta_j(t) = mode_j + t sd_j and the other parameters sit at their
# Gaussian conditional means (Jamil & Rue, 2026, Lemma 3.1). The log-heights
# are corrected for the change in the volume of the other parameters along
# the path by adding t gamma1_j (Jamil & Rue, 2026, Corollary 3.3), and a
# skew-normal is fitted to the corrected log-heights on a grid of t.
#
# The correction is the linear term of the width gamma_j(t) =
# -1/2 log|H_{-j}(theta(t))|, the log-volume of the other parameters. Two
# more accurate scans can replace it, each removing one error (the `width`
# argument of inla_vanilla()): "cmp" takes the width exactly along the CMP,
# removing the linearisation, and "tk" moves the other parameters to their
# conditional mode at each point and takes height and width there, the
# Laplace marginal of Tierney and Kadane (1986), removing the error of the
# straight path as well. The skew-normal fit is the same for all three.
#
# Jamil, H., & Rue, H. (2026). Approximate Bayesian inference for structural
# equation models using integrated nested Laplace approximations
# (arXiv:2603.25690 [stat.ME]). arXiv.
# https://doi.org/10.48550/arXiv.2603.25690
#
# Tierney, L., & Kadane, J. B. (1986). Accurate approximations for posterior
# moments and marginal densities. Journal of the American Statistical
# Association, 81(393), 82-86. https://doi.org/10.1080/01621459.1986.10478240

# The volume correction gamma1_j: the rate of change, along v_j, of minus half
# the log-determinant of the conditional precision of the other parameters. In
# whitened coordinates (columns of L) the precision at the mode is the
# identity, so its trace is m and v_j' H v_j = 1 there.

# The correction for every parameter at once, from the gradient. With T the
# third derivative of the log-posterior at the mode,
#   gamma1_j = 1/2 (sum_k T(A_k, A_k, v_j) - T(v_j, v_j, v_j)),
# and T(a, a, .) is the second difference of the gradient along a. So one
# vector, the sum over k, holds the trace part of every parameter, and each
# v_j adds one more second difference: 4m + 1 gradient calls for all m
# parameters. The step is `h` posterior SDs. Truncation grows as h^2, and noise
# in the gradient, such as that of an inner optimisation, as 1 / h^2. Any A
# with A A' = Sigma gives the same sum in the limit. fit_marginals() passes the
# principal axes, so that the truncation does not depend on the order of the
# parameters.
gamma1_grad <- function(gr, mode, Vscan, A, h = 0.1) {
  m <- length(mode)
  g0 <- gr(mode)
  second_diff <- function(a) {
    (gr(mode + h * a) + gr(mode - h * a) - 2 * g0) / h^2
  }
  u <- numeric(m)
  for (k in seq_len(m)) {
    u <- u + second_diff(A[, k])
  }
  vapply(
    seq_len(m),
    function(j) 0.5 * sum(Vscan[, j] * (u - second_diff(Vscan[, j]))),
    numeric(1)
  )
}

# The correction from the log-posterior alone, by a forward difference of
# `d_outer` posterior SDs from the mode. Each whitened curvature is a second
# difference of the log-posterior, 2m + 5 calls per parameter. Second
# differences carry a truncation error of their own, so the baselines at the
# mode (`tr0`, and the curvature along v_j computed here) are taken by the same
# rule rather than set to m and 1, and the error cancels in the forward
# difference.
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

# The width term without the linearisation: the volume of the other
# parameters, gamma_j(x) = -1/2 log|H_{-j}(x)|, from the negative Hessian
# H_{-j} of the log-posterior in the other parameters at the point `x`. It
# comes from one call of the user's Hessian, from central differences of the
# gradient in the other coordinates (2(m - 1) gradient calls), or from second
# differences of the log-posterior there ((m - 1)m calls, with `f0` the value
# at `x`). The steps are those of the Hessian at the mode in laplace_at(), so
# the volume at the mode, taken from lap$H, is on the same footing. NA where
# H_{-j} is not positive definite or cannot be evaluated.
cond_volume <- function(x, j, lp, gr, he, lap, f0 = NULL) {
  others <- seq_along(x)[-j]
  h <- lap$scales[others]
  Hb <- if (!is.null(he)) {
    -he(x)[others, others, drop = FALSE]
  } else if (!is.null(gr)) {
    gr_others <- function(y) {
      x[others] <- y
      gr(x)[others]
    }
    -fd_jacobian(gr_others, x[others], h = 1e-4 * h)
  } else {
    lp_others <- function(y) {
      x[others] <- y
      lp(x)
    }
    if (is.null(f0)) {
      f0 <- lp(x)
    }
    -fd_hessian(lp_others, x[others], h = 0.01 * h, f0 = f0)
  }
  half_log_det(Hb)
}

# Minus half the log-determinant of a negative Hessian, symmetrised, or NA
# when it is not positive definite.
half_log_det <- function(H) {
  R <- chol_or_null(0.5 * (H + t(H)))
  if (is.null(R)) NA_real_ else -sum(log(diag(R)))
}

# Level 2, "cmp": the exact width term along the CMP, gamma_j(t) - gamma_j(0)
# at each scan point, with the heights `raw` already known. The volume at
# the mode is that of lap$H. Points where the height is not finite are left
# NA, with no calls.
width_cmp <- function(j, lp, gr, he, lap, Vscan, z, raw) {
  g0 <- half_log_det(lap$H[-j, -j, drop = FALSE])
  vapply(
    seq_along(z),
    function(k) {
      if (z[k] == 0) {
        return(0)
      }
      if (!is.finite(raw[k])) {
        return(NA_real_)
      }
      x <- lap$mode + Vscan[, j] * z[k]
      cond_volume(x, j, lp, gr, he, lap, f0 = raw[k]) - g0
    },
    numeric(1)
  )
}

# Level 1, "tk": the Laplace marginal of Tierney and Kadane (1986). At each
# scan point, theta_j = mode_j + t sd_j, the other parameters are moved to
# their conditional mode, found by stats::nlminb() from the CMP point, in the
# whitened coordinates of their conditional Laplace covariance C_j (where the
# conditional posterior has unit curvature at the mode). The height is the
# log-posterior there and the width the volume there, relative to the mode.
# Returns both, with NA where the conditional mode or its volume could not be
# found. Needs m > 1.
scan_tk <- function(j, lp, gr, he, lap, Vscan, z) {
  others <- seq_along(lap$mode)[-j]
  height <- width <- rep(NA_real_, length(z))
  g0 <- half_log_det(lap$H[-j, -j, drop = FALSE])
  # u -> theta_{-j}: x_{-j} = start + Lc u, with Lc Lc' = H_{-j}^{-1}
  Rc <- chol_or_null(lap$H[others, others, drop = FALSE])
  Lc <- if (is.null(Rc)) {
    diag(lap$scales[others], length(others))
  } else {
    backsolve(Rc, diag(length(others)))
  }
  for (k in seq_along(z)) {
    if (z[k] == 0) {
      height[k] <- lap$lp_max
      width[k] <- 0
      next
    }
    x0 <- lap$mode + Vscan[, j] * z[k]
    to_x <- function(u) {
      x <- x0
      x[others] <- x0[others] + drop(Lc %*% u)
      x
    }
    ob <- function(u) {
      v <- lp(to_x(u))
      if (is.finite(v)) -v else Inf
    }
    gro <- if (is.null(gr)) {
      NULL
    } else {
      function(u) -drop(crossprod(Lc, gr(to_x(u))[others]))
    }
    heo <- if (is.null(he)) {
      NULL
    } else {
      function(u) {
        Hb <- -he(to_x(u))[others, others, drop = FALSE]
        crossprod(Lc, Hb %*% Lc)
      }
    }
    opt <- tryCatch(
      stats::nlminb(
        numeric(length(others)),
        ob,
        gradient = gro,
        hessian = heo,
        control = list(iter.max = 200L, eval.max = 400L)
      ),
      error = function(e) NULL
    )
    if (is.null(opt) || !is.finite(opt$objective)) {
      next
    }
    x <- to_x(opt$par)
    height[k] <- -opt$objective
    width[k] <- cond_volume(x, j, lp, gr, he, lap, f0 = height[k]) - g0
  }
  list(height = height, width = width)
}

# The slope of a width curve at the mode, per posterior SD: the central
# difference across t = 0 between the nearest scan points either side. For
# the linear trace it is gamma1 itself.
width_slope <- function(z, w) {
  l <- which(z < 0)
  r <- which(z > 0)
  if (length(l) == 0L || length(r) == 0L) {
    return(NA_real_)
  }
  l <- l[which.max(z[l])]
  r <- r[which.min(z[r])]
  (w[r] - w[l]) / (z[r] - z[l])
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
#
# The width term (the volume of the other parameters) comes from `width`:
# "trace" (the default) adds the linear term gamma1_j t, "cmp" the exact
# volume along the CMP (width_cmp()), and "tk" replaces the scan by the
# Laplace marginal of Tierney and Kadane (scan_tk()), heights included. With
# `tk_fallback`, a parameter whose ridge bend is flagged (above
# summary_levels$ridge_bend) is scanned by "tk" whatever `width` says. The
# linear volume correction takes gamma1_grad() with a gradient, whether or not
# there is a Hessian, because a Hessian per parameter costs far more than the
# gradients it would save, and gamma1_lp() otherwise. With one parameter there
# are no other parameters, and the two agree.
#
# Returns `sn`, a matrix with a row per parameter (the skew-normal on the
# theta scale, the diagnostics, and `logC`, the log scale factor that sets
# the fitted density's height against the scan normalised to a peak of 1),
# the scan grid `z`, `logpost`, the log-heights along each scan (a column per
# parameter: along the CMP, or at the conditional modes for "tk"), `width`,
# the width term added to them (zero at the mode), and `method`, the width
# method each parameter used. plot(type = "sn_fit") draws them as INLAvaan's
# visual_debug() does.
fit_marginals <- function(
  lp,
  gr,
  lap,
  ctrl,
  tick = function() NULL,
  he = NULL,
  width = "trace",
  tk_fallback = FALSE
) {
  mode <- lap$mode
  Sigma <- lap$Sigma
  L <- lap$L
  m <- length(mode)
  sds <- sqrt(diag(Sigma))
  Vscan <- sweep(Sigma, 2, sds, "/")
  z <- seq(-4, 4, length.out = ctrl$n_grid)

  # The trace at the mode for gamma1_lp(), only when some parameter uses it
  tr0 <- if (is.null(gr) && (width == "trace" || m == 1L)) {
    whitened_trace_lp(lp, mode, lap$lp_max, L)
  } else {
    NULL
  }
  # The linear width term from the gradient, for all the parameters at once
  g1_grad <- if (!is.null(gr) && (width == "trace" || m == 1L)) {
    gamma1_grad(gr, mode, Vscan, principal_axes(Sigma))
  } else {
    NULL
  }

  out <- matrix(
    NA_real_,
    m,
    8,
    dimnames = list(
      names(mode),
      c("xi", "omega", "alpha", "gamma1", "nmad", "tails", "ridge_bend", "logC")
    )
  )
  heights <- widths <- matrix(
    NA_real_,
    length(z),
    m,
    dimnames = list(NULL, names(mode))
  )
  method <- stats::setNames(rep(width, m), names(mode))
  # Scan points for the ridge bend (see scan_bend())
  i_bend <- unique(c(which.min(abs(z + 2)), which.min(abs(z - 2))))
  # Log-posterior along the CMP; the mode's height is known
  cmp_height <- function(j, t) {
    vapply(
      t,
      function(t) if (t == 0) lap$lp_max else lp(mode + Vscan[, j] * t),
      numeric(1)
    )
  }
  above <- function(h, y) {
    is.finite(h) && h - max(y[is.finite(y)]) > -3
  }

  for (j in seq_len(m)) {
    # The heights along the CMP, which "tk" does not need beyond the ridge
    # bend, and that only without a gradient
    raw <- if (width != "tk" || m == 1L) {
      cmp_height(j, z)
    } else {
      replace(
        rep(NA_real_, length(z)),
        i_bend,
        if (is.null(gr)) cmp_height(j, z[i_bend]) else NA_real_
      )
    }
    bend <- scan_bend(j, lp, gr, lap, Vscan, z[i_bend], raw[i_bend])
    if (tk_fallback && off_level(bend, summary_levels$ridge_bend)) {
      method[[j]] <- "tk"
    }

    # The width term, and the heights it is added to. With one parameter
    # there are no others, so the three methods are the same scan, taken by
    # the linear route as before.
    use <- if (m == 1L) "trace" else method[[j]]
    if (use == "trace") {
      g1 <- if (is.null(gr)) {
        gamma1_lp(j, lp, mode, lap$lp_max, Vscan, L, tr0)
      } else {
        g1_grad[[j]]
      }
      w <- g1 * z
    } else if (use == "cmp") {
      w <- width_cmp(j, lp, gr, he, lap, Vscan, z, raw)
      g1 <- width_slope(z, w)
    } else {
      tk <- scan_tk(j, lp, gr, he, lap, Vscan, z)
      raw <- tk$height
      w <- tk$width
      g1 <- width_slope(z, w)
    }

    # Corrected log-heights
    y <- raw + w
    tails <- !above(y[1], y) && !above(y[length(y)], y)

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
      bend,
      fit$logC
    )
    heights[, j] <- raw
    widths[, j] <- w
    tick()
  }
  list(sn = out, z = z, logpost = heights, width = widths, method = method)
}

# The mass each fitted marginal puts outside the window that was scanned to
# fit it, mode +- 4 Laplace SDs: there the skew-normal is extrapolated, so
# when this is large the outer quantiles are the ones to distrust. This is
# INLAvaan's `scan_end_mass` (in its diagnostics()), renamed. A Gaussian
# marginal gives 2 pnorm(-4) = 6.3e-5; larger values mean more of the
# marginal is extrapolated. summary() flags it, and the fit warns, above
# 0.05. INLAvaan's docs give 1e-3 to 1e-2 as the healthy range; its lower
# end fails for healthy fits here (see lab/DESIGN.md), so no range is given.
unscanned_mass <- function(sn, mode, sds) {
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
# parameters: the linear CMP is a first-order approximation to a possibly
# curved ridge (Jamil & Rue, 2026, Section 3.2.1), and this measures the
# curvature it ignores, a check of inlavanilla's own. At scan points t (near
# +-2 posterior SDs), one Newton step in the other parameters, with the
# Laplace curvature, would raise the log-posterior by 0.5 g' C g, where g is
# the gradient in the other parameters and C their conditional covariance
# given theta_j,
#   C = Sigma[-j, -j] - Sigma[-j, j] Sigma[j, -j] / Sigma[j, j].
# That is the height the scan misses there, in log-density units. When it is
# large, the other parameters' conditional mode bends away from the line (as
# a noise scale does for a regression coefficient at small n), and the
# marginal comes out too narrow. Returns the larger of the gains. The
# gradient is the user's, or one-sided differences in the other coordinates
# only (m - 1 calls per point, reusing the scan's height), with their bias,
# half a step times the curvature, removed with the Laplace curvature.
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
