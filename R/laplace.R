# Stage 1a: the joint Gaussian (Laplace) approximation N(mode, Sigma) with
# Sigma = H(mode)^-1, and the Laplace estimate of the log marginal likelihood.
#
# `lp` and `gr` are the wrapped log-posterior and gradient (`gr` is NULL when
# the user gave no gradient). Returns the mode, the log-posterior there, the
# upper Cholesky factor R of the precision (R'R = H), the covariance Sigma,
# a whitening factor L (L L' = Sigma), the optimiser output, the largest
# whitened gradient at the mode, and the Laplace log marginal likelihood.
laplace_fit <- function(lp, gr, init, ctrl) {
  ob <- function(x) {
    v <- lp(x)
    if (is.finite(v)) -v else Inf
  }
  gro <- if (is.null(gr)) NULL else function(x) -gr(x)
  opt <- stats::nlminb(
    start = unname(init),
    objective = ob,
    gradient = gro,
    control = ctrl$nlminb
  )
  cur <- laplace_at(lp, gr, opt$par, -opt$objective, names(init))

  # nlminb() works on the raw parameters, and its convergence tests assume
  # they are on a unit scale. Far from that (posterior SDs of 1e4, or 1e-6)
  # it can stop well short of the mode, reporting success or false
  # convergence. When the whitened gradient says the mode is off, restart
  # from there with each parameter scaled by its curvature.
  if (is.null(cur$R) || cur$grad_mode > 1e-3) {
    scale <- if (is.null(cur$R)) 1 / cur$scales else sqrt(diag(cur$H))
    opt2 <- stats::nlminb(
      start = opt$par,
      objective = ob,
      gradient = gro,
      scale = scale,
      control = ctrl$nlminb
    )
    if (opt2$objective <= opt$objective) {
      opt <- opt2
      cur <- laplace_at(lp, gr, opt$par, -opt$objective, names(init))
    }
  }

  if (is.null(cur$R)) {
    H <- cur$H
    problem <- if (all(is.finite(H))) {
      "is not positive definite"
    } else {
      "has values that are not finite"
    }
    abort(
      c(
        "The Hessian at the posterior mode {problem}.",
        "i" = "The optimiser may have stopped short of the mode, or the
               posterior may be improper or not identified.",
        "i" = "Optimiser message: {.val {opt$message}}."
      ),
    )
  }
  m <- length(init)
  c(
    cur,
    list(
      opt = opt,
      # log pi(y) ~ log pi(y, mode) - log N(mode | mode, Sigma), with
      # log det(Sigma) = -2 sum(log(diag(R)))
      logml = cur$lp_max + (m / 2) * log(2 * pi) - sum(log(diag(cur$R)))
    )
  )
}

# Hessian, its Cholesky factor, the covariance and the whitened gradient at
# `mode`. Finite-difference steps are set per coordinate from fd_scales(), so
# they suit each parameter's own scale: central differences of the gradient
# at 1e-4 of that scale, or, without a gradient, second differences of the
# log-posterior at a hundredth of it. R is NULL when the Hessian is not
# positive definite.
laplace_at <- function(lp, gr, mode, lp_max, parnames) {
  m <- length(mode)
  mode <- stats::setNames(mode, parnames)
  scales <- fd_scales(lp, mode, lp_max)
  H <- if (is.null(gr)) {
    -fd_hessian(lp, mode, h = 0.01 * scales, f0 = lp_max)
  } else {
    -fd_jacobian(gr, mode, h = 1e-4 * scales)
  }
  # A Hessian is symmetric; a Jacobian of the user's gradient that is not
  # (relative to the curvatures) means the gradient is wrong.
  d <- sqrt(abs(diag(H)))
  asymmetry <- if (all(is.finite(H)) && all(d > 0)) {
    max(abs(H - t(H)) / outer(d, d))
  } else {
    0
  }
  H <- 0.5 * (H + t(H))
  dimnames(H) <- list(parnames, parnames)
  R <- chol_or_null(H)
  out <- list(
    mode = mode,
    lp_max = lp_max,
    H = H,
    R = R,
    scales = scales,
    asymmetry = asymmetry,
    grad_mode = Inf
  )
  if (is.null(R)) {
    return(out)
  }
  out$L <- backsolve(R, diag(m)) # upper triangular, L L' = H^-1
  out$Sigma <- chol2inv(R)
  dimnames(out$Sigma) <- list(parnames, parnames)
  # The gradient at the mode in whitened units: how far, in posterior SDs, a
  # Newton step would still move the mode.
  g <- if (is.null(gr)) {
    fd_grad(lp, mode, h = 1e-3 * sqrt(diag(out$Sigma)))
  } else {
    gr(mode)
  }
  out$grad_mode <- max(abs(crossprod(out$L, g)))
  out
}

# A curvature scale for each coordinate: roughly the conditional posterior SD
# 1 / sqrt(H_ii), found by scaling a step by tens until the symmetric drop
# f(x) - (f(x + h e_i) + f(x - h e_i)) / 2 = H_ii h^2 / 2 lies between 1e-3
# and 1, a few calls per coordinate. A fixed step fails at the extremes: for
# a posterior SD of 1e4 the second difference at h = 1e-4 drowns in rounding
# error, and for an SD of 1e-6 it measures far beyond the quadratic region.
fd_scales <- function(fn, x, f0) {
  m <- length(x)
  s <- numeric(m)
  for (i in seq_len(m)) {
    h <- 1e-4 * max(1, abs(x[[i]]))
    drop <- NA_real_
    for (k in 1:30) {
      xp <- xm <- x
      xp[i] <- x[i] + h
      xm[i] <- x[i] - h
      drop <- f0 - 0.5 * (fn(xp) + fn(xm))
      if (!is.finite(drop) || drop > 1) {
        h <- h / 10 # too far: out of the quadratic region or the support
      } else if (drop < 1e-3) {
        h <- h * 10 # too close: lost in rounding error
      } else {
        break
      }
    }
    s[i] <- if (is.finite(drop) && drop > 0) h / sqrt(2 * drop) else h
  }
  s
}

chol_or_null <- function(H) {
  if (!all(is.finite(H))) {
    return(NULL)
  }
  tryCatch(chol(H), error = function(e) NULL)
}
