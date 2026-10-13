# Stage 1a: the joint Gaussian (Laplace) approximation N(mode, Sigma) with
# Sigma = H(mode)^-1, and the Laplace estimate of the log marginal likelihood.
#
# `lp` and `gr` are the wrapped log-posterior and gradient (`gr` is NULL when
# the user gave no gradient). Returns the mode, the log-posterior there, the
# upper Cholesky factor R of the precision (R'R = H), the covariance Sigma,
# a whitening factor L (L L' = Sigma), the optimiser output and the Laplace
# log marginal likelihood.
laplace_fit <- function(lp, gr, init, ctrl) {
  m <- length(init)
  parnames <- names(init)

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
  mode <- stats::setNames(opt$par, parnames)
  lp_max <- -opt$objective

  # Hessian at the mode. With a gradient, central differences of the gradient
  # (2m calls). Without one, second differences of the log-posterior, in two
  # passes: a rough pass with a fixed step gives the posterior SDs, and the
  # second pass steps a hundredth of an SD along each axis. Second differences
  # are sensitive to the step, and an absolute step that suits one parameter
  # can be far too large or too small for another.
  if (!is.null(gr)) {
    H <- -fd_jacobian(gr, mode, h = 1e-5)
  } else {
    H <- -fd_hessian(lp, mode, h = 1e-4, f0 = lp_max)
    R1 <- chol_or_null(H)
    if (!is.null(R1)) {
      sd1 <- sqrt(diag(chol2inv(R1)))
      H <- -fd_hessian(lp, mode, h = 1e-2 * sd1, f0 = lp_max)
    }
  }
  H <- 0.5 * (H + t(H))
  R <- chol_or_null(H)
  if (is.null(R)) {
    problem <- if (all(is.finite(H))) {
      "is not positive definite"
    } else {
      "has values that are not finite"
    }
    cli::cli_abort(c(
      "The Hessian at the posterior mode {problem}.",
      "i" = "The optimiser may have stopped short of the mode, or the
             posterior may be improper or not identified.",
      "i" = "Optimiser message: {.val {opt$message}}."
    ))
  }
  L <- backsolve(R, diag(m)) # upper triangular, L L' = H^-1
  Sigma <- chol2inv(R)
  dimnames(Sigma) <- dimnames(H) <- list(parnames, parnames)

  list(
    mode = mode,
    lp_max = lp_max,
    H = H,
    R = R,
    Sigma = Sigma,
    L = L,
    opt = opt,
    # log pi(y) ~ log pi(y, mode) - log N(mode | mode, Sigma), with
    # log det(Sigma) = -2 sum(log(diag(R)))
    logml = lp_max + (m / 2) * log(2 * pi) - sum(log(diag(R)))
  )
}

chol_or_null <- function(H) {
  if (!all(is.finite(H))) {
    return(NULL)
  }
  tryCatch(chol(H), error = function(e) NULL)
}
