# Stage 1a: the joint Gaussian (Laplace) approximation N(mode, Sigma) with
# Sigma = H(mode)^-1, and the Laplace estimate of the log marginal likelihood.
#
# `lp`, `gr` and `he` are the wrapped log-posterior, gradient and Hessian
# (`gr` and `he` are NULL when the user gave none). With `he`, nlminb() takes
# Newton steps with it, and the Hessian at the mode is the user's. Returns
# the mode, the log-posterior there, the upper Cholesky factor R of the
# precision (R'R = H), the covariance Sigma, a whitening factor L
# (L L' = Sigma), the optimiser output, the gradient of the log-posterior at
# the mode with the checks built on it, the condition number of the Hessian
# in posterior-SD units, and the Laplace log marginal likelihood. For verbose
# output it also says how the Hessian was taken (`hessian`) and whether
# nlminb() was restarted (`restart`, NULL when it was not).
laplace_fit <- function(lp, gr, init, ctrl, he = NULL) {
  ob <- function(x) {
    v <- lp(x)
    if (is.finite(v)) -v else Inf
  }
  gro <- if (is.null(gr)) NULL else function(x) -gr(x)
  heo <- if (is.null(he)) NULL else function(x) -he(x)
  opt <- stats::nlminb(
    start = unname(init),
    objective = ob,
    gradient = gro,
    hessian = heo,
    control = ctrl$nlminb
  )
  cur <- laplace_at(lp, gr, opt$par, -opt$objective, names(init), he)

  # nlminb() works on the raw parameters, and its convergence tests assume
  # they are on a unit scale. Far from that (posterior SDs of 1e4, or 1e-6)
  # it can stop well short of the mode, reporting success or false
  # convergence. When the whitened gradient says the mode is off, restart
  # from there with each parameter scaled by its curvature.
  restart <- NULL
  if (is.null(cur$R) || cur$grad_mode > 1e-3) {
    scale <- if (is.null(cur$R)) 1 / cur$scales else sqrt(diag(cur$H))
    opt2 <- stats::nlminb(
      start = opt$par,
      objective = ob,
      gradient = gro,
      hessian = heo,
      scale = scale,
      control = ctrl$nlminb
    )
    restart <- list(
      grad_mode = cur$grad_mode,
      iterations = opt$iterations,
      kept = opt2$objective <= opt$objective
    )
    if (restart$kept) {
      opt <- opt2
      cur <- laplace_at(lp, gr, opt$par, -opt$objective, names(init), he)
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
      )
    )
  }
  m <- length(init)
  # The condition number of the Hessian with each parameter in its own
  # posterior-SD units, D^1/2 H D^1/2 with D = diag(Sigma): the eigenvalue
  # ratio of its inverse, the correlation matrix cov2cor(Sigma). INLAvaan
  # takes it on the raw scale, which a mere change of units can make large;
  # posterior SDs here run from 1e-6 to 1e4.
  ev <- eigen(stats::cov2cor(cur$Sigma), symmetric = TRUE, only.values = TRUE)
  ev <- ev$values
  hess_cond <- if (min(ev) > 0) max(ev) / min(ev) else Inf
  c(
    cur,
    list(
      opt = opt,
      restart = restart,
      hess_cond = hess_cond,
      # log pi(y) ~ log pi(y, mode) - log N(mode | mode, Sigma), with
      # log det(Sigma) = -2 sum(log(diag(R)))
      logml = cur$lp_max + (m / 2) * log(2 * pi) - sum(log(diag(cur$R)))
    )
  )
}

# Hessian, its Cholesky factor, the covariance and the whitened gradient at
# `mode`. With `he` the Hessian is the user's, checked against the gradient
# by hess_discrepancy(). Otherwise finite-difference steps are set per
# coordinate from fd_scales(), so they suit each parameter's own scale:
# central differences of the gradient at 1e-4 of that scale, or, without a
# gradient, second differences of the log-posterior at a hundredth of it. R
# is NULL when the Hessian is not positive definite. `hessian` records the
# method and its calls, those of fd_scales() counted apart.
laplace_at <- function(lp, gr, mode, lp_max, parnames, he = NULL) {
  m <- length(mode)
  mode <- stats::setNames(mode, parnames)
  if (is.null(he)) {
    scales <- fd_scales(lp, mode, lp_max)
    H <- if (is.null(gr)) {
      -fd_hessian(lp, mode, h = 0.01 * scales, f0 = lp_max)
    } else {
      -fd_jacobian(gr, mode, h = 1e-4 * scales)
    }
  } else {
    H <- -he(mode)
    # The curvature scales, which only rescale a restart of nlminb(), from
    # the Hessian's own diagonal where it allows
    d <- diag(H)
    scales <- if (all(is.finite(d) & d > 0)) {
      structure(1 / sqrt(d), ok = rep(TRUE, m), calls = 0L)
    } else {
      fd_scales(lp, mode, lp_max)
    }
  }
  # A Hessian is symmetric. A Jacobian of the user's gradient that is not
  # (relative to the curvatures) means the gradient is wrong. A user Hessian
  # is checked against the gradient instead, by hess_discrepancy().
  d <- sqrt(abs(diag(H)))
  asymmetry <- if (!is.null(he)) {
    NA_real_
  } else if (all(is.finite(H)) && all(d > 0)) {
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
    hess_error = NA_real_,
    grad_mode = Inf,
    hessian = list(
      by = if (!is.null(he)) {
        "hess"
      } else if (is.null(gr)) {
        "logpost"
      } else {
        "grad"
      },
      calls = if (!is.null(he)) {
        1L
      } else if (is.null(gr)) {
        m^2 + m
      } else {
        2 * m
      },
      step_calls = attr(scales, "calls"),
      check_calls = 0L
    )
  )
  if (is.null(R)) {
    return(out)
  }
  out$L <- backsolve(R, diag(m)) # upper triangular, L L' = H^-1
  out$Sigma <- chol2inv(R)
  dimnames(out$Sigma) <- list(parnames, parnames)
  if (!is.null(he)) {
    out$hess_error <- hess_discrepancy(gr, mode, out$L)
    out$hessian$check_calls <- attr(out$hess_error, "calls")
  }
  # The gradient of the log-posterior at the mode, and two measures of how
  # far a Newton step would still move the mode. `grad_mode`, the largest
  # whitened component, decides the restart in laplace_fit().
  # `mode_shift` is INLAvaan's per-parameter Newton step in posterior SDs,
  # |(Sigma g)_j| / sd_j, which decides the fit-time check; it does not
  # depend on the order of the parameters.
  sds <- sqrt(diag(out$Sigma))
  g <- if (is.null(gr)) {
    fd_grad(lp, mode, h = 1e-3 * sds)
  } else {
    gr(mode)
  }
  out$grad <- stats::setNames(as.numeric(g), parnames)
  out$grad_mode <- max(abs(crossprod(out$L, g)))
  out$mode_shift <- stats::setNames(
    abs(as.numeric(out$Sigma %*% g)) / sds,
    parnames
  )
  out
}

# How far a user Hessian is from central differences of the gradient at the
# mode, in whitened coordinates, where the Hessian is the identity: the
# largest error in L'(-J)L w, for J the Jacobian of the gradient along L w,
# over two fixed unit directions w (all ones, and alternating signs). It
# costs 4 gradient calls whatever m, against 2m for a check of each column,
# and a wrong Hessian passes only if its error happens to vanish along both
# directions. Not finite when the gradient is not finite there.
hess_discrepancy <- function(gr, mode, L, h = 1e-4) {
  m <- length(mode)
  W <- unique(cbind(rep(1, m), rep_len(c(1, -1), m)), MARGIN = 2) / sqrt(m)
  err <- 0
  for (k in seq_len(ncol(W))) {
    u <- drop(L %*% W[, k])
    jl <- (gr(mode + h * u) - gr(mode - h * u)) / (2 * h)
    err <- max(err, abs(drop(crossprod(L, -jl)) - W[, k]))
  }
  structure(err, calls = 2L * ncol(W))
}

# A curvature scale for each coordinate: roughly the conditional posterior SD
# 1 / sqrt(H_ii), found by scaling a step by tens until the symmetric drop
# f(x) - (f(x + h e_i) + f(x - h e_i)) / 2 = H_ii h^2 / 2 lies between 1e-3
# and 1, a few calls per coordinate. A fixed step fails at the extremes: for
# a posterior SD of 1e4 the second difference at h = 1e-4 drowns in rounding
# error, and for an SD of 1e-6 it measures far beyond the quadratic region.
# The attribute `ok` marks the coordinates where the search succeeded, and
# `calls` counts the calls of `fn`.
fd_scales <- function(fn, x, f0) {
  m <- length(x)
  s <- numeric(m)
  ok <- logical(m)
  calls <- 0L
  for (i in seq_len(m)) {
    h <- 1e-4 * max(1, abs(x[[i]]))
    drop <- NA_real_
    for (k in 1:30) {
      xp <- xm <- x
      xp[i] <- x[i] + h
      xm[i] <- x[i] - h
      drop <- f0 - 0.5 * (fn(xp) + fn(xm))
      calls <- calls + 2L
      if (!is.finite(drop) || drop > 1) {
        h <- h / 10 # too far: out of the quadratic region or the support
      } else if (drop < 1e-3) {
        h <- h * 10 # too close: lost in rounding error
      } else {
        break
      }
    }
    ok[i] <- is.finite(drop) && drop >= 1e-3 && drop <= 1
    s[i] <- if (is.finite(drop) && drop > 0) h / sqrt(2 * drop) else h
  }
  # `ok` is FALSE where no step gave a drop in the window, such as where the
  # log-posterior is not locally concave
  structure(s, ok = ok, calls = calls)
}

chol_or_null <- function(H) {
  if (!all(is.finite(H))) {
    return(NULL)
  }
  tryCatch(chol(H), error = function(e) NULL)
}
