# Stage 1b: the variational Bayes mean correction. The Laplace Gaussian is
# recentred at mode + shift, keeping its covariance, where the shift maximises
# the Gaussian expectation of the log-posterior. At the optimum the expected
# gradient vanishes, E[grad log pi(mode + shift + z)] = 0 for z ~ N(0, Sigma).
#
# Ported from INLAvaan: the same fixed-point iteration with Anderson(1)
# acceleration and an nlminb() fallback, using its deterministic Gauss-Hermite
# node set (see vb_nodes()).

# Node set for the Gaussian expectation: the three-point Gauss-Hermite rule
# along each principal axis of Sigma, sharing the centre, so 2m + 1 nodes.
# Nodes are the centre and +/- sqrt(3) times each principal axis, with weights
# 1 - m/3 and 1/6. In whitened coordinates every odd moment and every cross
# moment E[z_r z_s] vanishes, so the rule is exact for the Gaussian expectation
# of the gradient whenever the log-posterior is a quartic polynomial in
# whitened coordinates. It gives the same shift on every run. The centre
# weight is negative once m > 3; sums are taken as differences from the centre
# so the large terms do not cancel.
vb_nodes <- function(Sigma) {
  m <- ncol(Sigma)
  e <- eigen(Sigma, symmetric = TRUE)
  axes <- sweep(e$vectors, 2, sqrt(pmax(e$values, 0)), "*")
  list(
    nodes = rbind(0, sqrt(3) * t(axes), -sqrt(3) * t(axes)),
    weights = c(1 - m / 3, rep(1 / 6, 2 * m))
  )
}

vb_correct <- function(lp, gr, mode, Sigma, R, L, ctrl) {
  m <- length(mode)
  rule <- vb_nodes(Sigma)
  Z <- rule$nodes
  w <- rule$weights
  nn <- nrow(Z)
  sds <- sqrt(diag(Sigma))

  # Expected gradient at mode + shift: the weighted mean over the nodes,
  # written as differences from the centre gradient.
  score <- function(shift) {
    mu <- mode + shift
    g0 <- gr(mu)
    s <- g0
    for (b in 2:nn) {
      s <- s + w[b] * (gr(mu + Z[b, ]) - g0)
    }
    s
  }

  # Negative expected log-posterior, for the fallback optimiser. A failed node
  # fails the objective, since the negative centre weight would otherwise
  # reward a -Inf there.
  objective <- function(shift) {
    mu <- mode + shift
    lps <- numeric(nn)
    for (b in seq_len(nn)) {
      lps[b] <- lp(mu + Z[b, ])
    }
    if (!all(is.finite(lps))) {
      return(1e40)
    }
    -(lps[1] + sum(w[-1] * (lps[-1] - lps[1])))
  }

  # Fixed-point iteration shift <- shift + Sigma %*% score(shift): with the
  # curvature of the Laplace approximation, the Newton step for the
  # stationarity condition is one multiplication by Sigma. Convergence is
  # judged on the step in posterior SDs. A step over `maxstep` SDs, or a
  # non-finite score, means the quadratic premise has failed (a strongly
  # skewed posterior), and the solve falls back to nlminb().
  maxstep <- 1
  shift <- numeric(m)
  iter <- 0L
  move <- Inf
  fallback <- FALSE
  step_prev <- shift_prev <- NULL
  for (it in seq_len(ctrl$vb_maxit)) {
    s <- score(shift)
    if (!all(is.finite(s))) {
      fallback <- TRUE
      break
    }
    step <- as.numeric(Sigma %*% s)
    move <- max(abs(step) / sds)
    if (move > maxstep) {
      fallback <- TRUE
      break
    }
    if (move < ctrl$vb_tol) {
      shift <- shift + step
      iter <- it
      break
    }
    # Anderson(1) acceleration: a secant estimate from consecutive steps
    # removes the dominant error mode. The plain step is kept whenever the
    # secant is degenerate or would move too far.
    update <- step
    if (!is.null(step_prev)) {
      df <- step - step_prev
      dx <- shift - shift_prev
      denom <- sum(df^2)
      if (denom > 0) {
        gam <- sum(step * df) / denom
        cand <- step - gam * (dx + df)
        if (
          all(is.finite(cand)) &&
            max(abs(cand) / sds) <= min(maxstep, 2 * move)
        ) {
          update <- cand
        }
      }
    }
    step_prev <- step
    shift_prev <- shift
    shift <- shift + update
    iter <- it
  }
  if (move >= ctrl$vb_tol) {
    fallback <- TRUE
  }

  if (fallback) {
    # Optimise in whitened coordinates, where the problem is well conditioned
    nl <- stats::nlminb(
      start = numeric(m),
      objective = function(d) objective(as.numeric(L %*% d)),
      gradient = function(d) {
        as.numeric(crossprod(L, -score(as.numeric(L %*% d))))
      },
      control = list(rel.tol = 1e-8)
    )
    shift <- as.numeric(L %*% nl$par)
    iter <- nl$iterations
  }

  names(shift) <- names(mode)
  list(
    shift = shift,
    iterations = iter,
    fallback = fallback,
    # Joint KL divergence between the Laplace Gaussian and its shifted copy,
    # 0.5 shift' H shift. It is the location term the Laplace estimate of the
    # marginal likelihood misses, and is added to it.
    kld = 0.5 * sum((R %*% shift)^2)
  )
}
