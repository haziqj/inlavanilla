# Tests of the derivatives and of the VB mean correction:
#
# (a) fits with an analytic gradient and with finite differences agree, up to
#     the error of the VB rule used without a gradient;
# (b) the two routes to the volume correction, gamma1_grad() and gamma1_lp(),
#     agree with each other and with the derivative of the conditional
#     log-determinant, and vanish for a Gaussian;
# (c) fd_grad(), fd_jacobian() and fd_hessian() against analytic derivatives,
#     including their exact error terms on polynomials;
# (d) vb_nodes() integrates Gaussian moments up to degree five exactly;
# (e) the VB shift against a reference VB solution computed here by tensor
#     Gauss-Hermite quadrature (nodes by the Golub-Welsch method);
# (f) vb = FALSE against vb = TRUE;
# (g) the nlminb() fallback of vb_correct().
#
# The reference VB solution solves E_q[grad log p(mu + L z)] = 0 for mu, with
# z ~ N(0, I) and L L' = fit$cov held fixed, as the package does. Location
# errors are in posterior SDs (sqrt(diag(fit$cov))). Tolerances on
# approximations were calibrated by printing the actual errors first; each
# comment states the observed value.

# Helpers ----------------------------------------------------------------------

# Gauss-Hermite rule for the standard normal by the Golub-Welsch method: the
# nodes are the eigenvalues of the Jacobi matrix of the probabilists' Hermite
# polynomials (zero diagonal, off-diagonal sqrt(k)), and the weights are the
# squared first components of the normalised eigenvectors.
gauss_hermite <- function(n) {
  off <- sqrt(seq_len(n - 1))
  jacobi <- matrix(0, n, n)
  jacobi[cbind(seq_len(n - 1), 2:n)] <- off
  jacobi[cbind(2:n, seq_len(n - 1))] <- off
  e <- eigen(jacobi, symmetric = TRUE)
  list(nodes = e$values, weights = e$vectors[1, ]^2)
}

# Tensor-product Gauss-Hermite rule for N(0, I_m), n nodes per axis.
tensor_gauss_hermite <- function(m, n) {
  rule <- gauss_hermite(n)
  idx <- as.matrix(expand.grid(rep(list(seq_len(n)), m)))
  list(
    nodes = matrix(rule$nodes[idx], ncol = m),
    weights = apply(matrix(rule$weights[idx], ncol = m), 1, prod)
  )
}

# Expected gradient E[grad(mu + x)] over a node set `x` (rows) with weights.
expected_grad <- function(grad, mu, x, weights) {
  m <- length(mu)
  g <- vapply(seq_len(nrow(x)), function(i) grad(mu + x[i, ]), numeric(m))
  colSums(weights * matrix(g, ncol = m, byrow = TRUE))
}

# Reference VB mean with the covariance fixed at `Sigma`: the root of the
# expected gradient under N(mu, Sigma), by tensor Gauss-Hermite quadrature and
# the iteration mu <- mu + Sigma E[grad], run to 1e-13 SDs.
vb_reference <- function(grad, start, Sigma, n_gh = 30, maxit = 1000) {
  m <- length(start)
  rule <- tensor_gauss_hermite(m, n_gh)
  x <- rule$nodes %*% chol(Sigma) # rows ~ N(0, Sigma)
  sds <- sqrt(diag(Sigma))
  mu <- start
  for (it in seq_len(maxit)) {
    step <- drop(Sigma %*% expected_grad(grad, mu, x, rule$weights))
    mu <- mu + step
    if (max(abs(step) / sds) < 1e-13) {
      return(mu)
    }
  }
  stop("The reference VB iteration did not converge.")
}

# Whitened expected gradient under the package's own rule, L' E[grad], for
# L L' = Sigma: zero at a stationary point of the expected log-posterior.
rule_stationarity <- function(grad, mu, Sigma) {
  rule <- vb_nodes(Sigma)
  g <- expected_grad(grad, mu, rule$nodes, rule$weights)
  drop(crossprod(t(chol(Sigma)), g))
}

# Largest absolute difference between two vectors (or row-wise matrices) in
# units of the posterior SDs `sds`.
max_sd_diff <- function(a, b, sds) {
  max(abs(a - b) / sds)
}

# All discrepancies between two fits of the same model, on comparable scales.
fit_discrepancy <- function(fit_a, fit_b) {
  sds <- sqrt(diag(fit_a$cov))
  c(
    mode = max_sd_diff(fit_a$mode, fit_b$mode, sds),
    cov = max(abs(fit_a$cov - fit_b$cov) / outer(sds, sds)),
    xi = max_sd_diff(fit_a$sn[, "xi"], fit_b$sn[, "xi"], sds),
    omega = max(abs(fit_a$sn[, "omega"] / fit_b$sn[, "omega"] - 1)),
    alpha = max(abs(fit_a$sn[, "alpha"] - fit_b$sn[, "alpha"])),
    summary = max_sd_diff(fit_a$summary, fit_b$summary, sds),
    summary_table = max_sd_diff(
      summary(fit_a)$table,
      summary(fit_b)$table,
      sds
    ),
    logml = abs(fit_a$logml - fit_b$logml),
    logml_laplace = abs(fit_a$logml_laplace - fit_b$logml_laplace),
    mean_vb = max_sd_diff(fit_a$mean_vb, fit_b$mean_vb, sds)
  )
}

# Models -----------------------------------------------------------------------

# y_i ~ N(mu, sigma^2), mu ~ N(0, 10^2), sigma ~ Exp(1), on (mu, log_sigma)
# with the Jacobian. Skewed in log_sigma.
normal_model <- function() {
  set.seed(1)
  y <- stats::rnorm(20, 5, 2)
  n <- length(y)
  list(
    init = c(mu = 0, log_sigma = 0),
    logpost = function(theta) {
      mu <- theta[["mu"]]
      log_sigma <- theta[["log_sigma"]]
      sum(stats::dnorm(y, mu, exp(log_sigma), log = TRUE)) +
        stats::dnorm(mu, 0, 10, log = TRUE) +
        stats::dexp(exp(log_sigma), 1, log = TRUE) +
        log_sigma
    },
    grad = function(theta) {
      mu <- theta[["mu"]]
      log_sigma <- theta[["log_sigma"]]
      s2 <- exp(2 * log_sigma)
      c(
        sum(y - mu) / s2 - mu / 100,
        -n + sum((y - mu)^2) / s2 - exp(log_sigma) + 1
      )
    },
    hessian = function(theta) {
      mu <- theta[["mu"]]
      log_sigma <- theta[["log_sigma"]]
      s2 <- exp(2 * log_sigma)
      h12 <- -2 * sum(y - mu) / s2
      matrix(
        c(
          -n / s2 - 1 / 100,
          h12,
          h12,
          -2 * sum((y - mu)^2) / s2 - exp(log_sigma)
        ),
        2
      )
    }
  )
}

# Logistic regression with an intercept and four covariates, beta ~ N(0, 2.5^2).
logistic_model <- function() {
  set.seed(12)
  n <- 80
  X <- cbind(1, matrix(stats::rnorm(n * 4), n))
  y <- stats::rbinom(n, 1, stats::plogis(X %*% c(-0.5, 1, -0.8, 0.5, 0)))
  list(
    init = stats::setNames(numeric(5), paste0("b", 0:4)),
    logpost = function(beta) {
      eta <- drop(X %*% beta)
      sum(y * eta - log1p(exp(eta))) +
        sum(stats::dnorm(beta, 0, 2.5, log = TRUE))
    },
    grad = function(beta) {
      eta <- drop(X %*% beta)
      drop(crossprod(X, y - stats::plogis(eta))) - beta / 2.5^2
    }
  )
}

# Poisson regression with few counts, y_i ~ Poisson(exp(x_i' beta)),
# beta ~ N(0, 2^2), and its analytic Hessian. `p` is 2 or 3.
poisson_regression_model <- function(p) {
  if (p == 2) {
    set.seed(4)
    x <- seq(-1, 1, length.out = 10)
    X <- cbind(1, x)
    y <- stats::rpois(10, exp(0.3 + 0.8 * x))
  } else {
    set.seed(3)
    X <- cbind(1, stats::rnorm(12), stats::rnorm(12))
    y <- stats::rpois(12, exp(X %*% c(0.2, 0.6, -0.4)))
  }
  list(
    init = stats::setNames(numeric(p), paste0("b", seq_len(p) - 1)),
    logpost = function(beta) {
      eta <- drop(X %*% beta)
      sum(y * eta - exp(eta)) - sum(beta^2) / 8
    },
    grad = function(beta) {
      eta <- drop(X %*% beta)
      drop(crossprod(X, y - exp(eta))) - beta / 4
    },
    hessian = function(beta) {
      eta <- drop(X %*% beta)
      -crossprod(X, X * exp(eta)) - diag(p) / 4
    }
  )
}

# An exactly Gaussian log-posterior N(mu0, P^-1) in three correlated
# parameters, with log normalising constant 1.5 log(2 pi) - 0.5 log det P.
gaussian_model <- function() {
  P <- matrix(c(2, 0.8, -0.3, 0.8, 1.5, 0.4, -0.3, 0.4, 0.7), 3)
  mu0 <- c(1, -2, 0.5)
  list(
    init = c(a = 0, b = 0, c = 0),
    precision = P,
    mean = mu0,
    logml = 1.5 * log(2 * pi) - 0.5 * determinant(P)$modulus[[1]],
    logpost = function(x) -0.5 * sum((x - mu0) * (P %*% (x - mu0))),
    grad = function(x) -drop(P %*% (x - mu0))
  )
}

# Poisson-gamma on theta = log(lambda): the log-posterior is
# A theta - B exp(theta) + const with A = a + sum(y) and B = b + n. The mode is
# log(A / B) and the Laplace variance 1 / A. Under N(mu, s^2) the expected
# gradient A - B exp(mu + s^2 / 2) vanishes at mu = log(A / B) - s^2 / 2 (exact
# VB). The three-point rule replaces exp(s^2 / 2) by its approximation
# 2 / 3 + cosh(sqrt(3) s) / 3, so its stationary point is also closed form.
poisson_gamma_model <- function(y, a, b) {
  A <- a + sum(y)
  B <- b + length(y)
  list(
    init = c(log_lambda = 0),
    A = A,
    B = B,
    logpost = function(theta) {
      lambda <- exp(theta[[1]])
      sum(stats::dpois(y, lambda, log = TRUE)) +
        stats::dgamma(lambda, a, b, log = TRUE) +
        theta[[1]]
    },
    grad = function(theta) A - B * exp(theta[[1]]),
    hessian = function(theta) matrix(-B * exp(theta[[1]])),
    vb_exact = function(s) log(A / B) - s^2 / 2,
    vb_rule = function(s) log(A / B) - log(2 / 3 + cosh(sqrt(3) * s) / 3)
  )
}

# (a) Analytic gradient against finite differences ----------------------------

test_that("fits with and without a gradient agree on a two-parameter normal model", {
  mod <- normal_model()
  fit_grad <- inla_vanilla(
    mod$logpost,
    mod$init,
    grad = mod$grad,
    silent = TRUE
  )
  fit_fd <- inla_vanilla(mod$logpost, mod$init, silent = TRUE)
  expect_true(fit_grad$has_grad)
  expect_false(fit_fd$has_grad)
  expect_identical(fit_fd$counts[["grad"]], 0L)

  d <- fit_discrepancy(fit_grad, fit_fd)
  # Laplace and marginal stages
  expect_lt(d[["mode"]], 5e-7) # observed 1.4e-7 SDs
  expect_lt(d[["cov"]], 7e-6) # observed 2.3e-6
  expect_lt(d[["logml_laplace"]], 1.2e-6) # observed 4.0e-7
  expect_lt(d[["summary"]], 7e-6) # observed 2.4e-6 SDs
  expect_lt(d[["summary_table"]], 7e-6) # observed 2.4e-6 SDs
  # The skew-normal parameters trade off against each other where the shape
  # is weakly identified (alpha = -0.3 for mu), so they move more than the
  # summaries they produce.
  expect_lt(d[["xi"]], 8e-5) # observed 2.7e-5 SDs
  expect_lt(d[["omega"]], 2e-5) # observed 6.5e-6
  expect_lt(d[["alpha"]], 1e-4) # observed 3.8e-5
  # Both volume-correction routes are exercised, and agree
  expect_lt(
    max(abs(
      fit_fd$diagnostics$marginals$gamma1 -
        fit_grad$diagnostics$marginals$gamma1
    )),
    1e-6 # observed 3.4e-7
  )

  # VB stage. Without a gradient the expected gradient comes from
  # log-posterior values alone, by Stein's identity on the simplex rule,
  # which is exact for log-densities of degree four rather than six. The two
  # VB means differ by the error of that rule, for a shift of 0.29 SDs, and
  # the log marginal likelihood through its VB term.
  expect_lt(d[["mean_vb"]], 0.015) # observed 5.1e-3 SDs
  expect_lt(d[["logml"]], 4.5e-3) # observed 1.5e-3
})

test_that("fits with and without a gradient agree on a five-parameter logistic regression", {
  mod <- logistic_model()
  fit_grad <- inla_vanilla(
    mod$logpost,
    mod$init,
    grad = mod$grad,
    silent = TRUE
  )
  fit_fd <- inla_vanilla(mod$logpost, mod$init, silent = TRUE)
  # A genuinely skewed posterior, so the comparison is not vacuous
  expect_gt(max(abs(fit_grad$sn[, "alpha"])), 1)

  d <- fit_discrepancy(fit_grad, fit_fd)
  # Without a gradient nlminb() differences the objective itself and stops
  # slightly further from the mode. These values do not move when the data
  # are perturbed by 1e-14 relative, so they are not rounding luck, but other
  # datasets (seeds 1 to 30) gave up to 7e-5 SDs in the mode and 2e-5 in the
  # summaries.
  expect_lt(d[["mode"]], 2e-5) # observed 5.7e-6 SDs
  expect_lt(d[["cov"]], 8.5e-6) # observed 2.8e-6
  expect_lt(d[["logml_laplace"]], 6e-6) # observed 1.9e-6
  expect_lt(d[["summary"]], 8e-6) # observed 2.6e-6 SDs
  expect_lt(d[["summary_table"]], 8e-6) # observed 2.6e-6 SDs
  expect_lt(d[["xi"]], 1.3e-6) # observed 4.2e-7 SDs
  expect_lt(d[["omega"]], 3e-6) # observed 1.0e-6
  expect_lt(d[["alpha"]], 9e-6) # observed 2.9e-6
  expect_lt(
    max(abs(
      fit_fd$diagnostics$marginals$gamma1 -
        fit_grad$diagnostics$marginals$gamma1
    )),
    4e-6 # observed 1.3e-6
  )
  # VB stage, from log-posterior values without a gradient, as above, for a
  # shift of 0.40 SDs
  expect_lt(d[["mean_vb"]], 7.5e-3) # observed 2.1e-3 SDs
  expect_lt(d[["logml"]], 2e-3) # observed 4.4e-4
})

# (b) The volume correction ----------------------------------------------------

test_that("gamma1_grad() and gamma1_lp() agree for every parameter of a skewed model", {
  mod <- poisson_regression_model(3)
  lap <- laplace_fit(mod$logpost, mod$grad, mod$init, vanilla_control(list()))
  m <- length(lap$mode)
  sds <- sqrt(diag(lap$Sigma))
  v_scan <- sweep(lap$Sigma, 2, sds, "/")
  tr0 <- whitened_trace_lp(mod$logpost, lap$mode, lap$lp_max, lap$L)

  # Independent reference: gamma1_j is the slope along the scan direction of
  # minus half the log-determinant of the conditional precision of the other
  # parameters, -0.5 log det(H[-j, -j]), with the analytic Hessian and a
  # central difference in t.
  log_vol <- function(j, t) {
    H <- -mod$hessian(lap$mode + t * v_scan[, j])
    -0.5 * determinant(H[-j, -j, drop = FALSE])$modulus[[1]]
  }

  g_grad <- g_lp <- g_ref <- g_grad_half <- g_lp_half <- numeric(m)
  for (j in seq_len(m)) {
    g_grad[j] <- gamma1_grad(j, mod$grad, lap$mode, v_scan, lap$L)
    g_lp[j] <- gamma1_lp(
      j,
      mod$logpost,
      lap$mode,
      lap$lp_max,
      v_scan,
      lap$L,
      tr0
    )
    # The same, with the forward difference over half the distance
    g_grad_half[j] <- gamma1_grad(
      j,
      mod$grad,
      lap$mode,
      v_scan,
      lap$L,
      d_outer = 0.005
    )
    g_lp_half[j] <- gamma1_lp(
      j,
      mod$logpost,
      lap$mode,
      lap$lp_max,
      v_scan,
      lap$L,
      tr0,
      d_outer = 0.005
    )
    g_ref[j] <- (log_vol(j, 1e-4) - log_vol(j, -1e-4)) / 2e-4
  }

  # The model is skewed enough for the correction to matter
  expect_gt(max(abs(g_ref)), 0.1)
  # The two routes agree with each other closely ...
  expect_lt(max(abs(g_grad - g_lp)), 1e-6) # observed 2.0e-7
  # ... and with the exact slope to within the O(d_outer) error of a forward
  # difference 0.01 SDs from the mode (observed 2.3e-4 to 2.4e-4 for every
  # parameter, both routes)
  expect_lt(max(abs(g_grad - g_ref)), 6e-4)
  expect_lt(max(abs(g_lp - g_ref)), 6e-4)
  # That error is the first-order truncation and nothing else: halving
  # d_outer halves it (observed ratios 2.000 to 2.001), and the Richardson
  # extrapolation 2 gamma1(d / 2) - gamma1(d) removes it, leaving 1.7e-7
  # (gradient route) and 7.4e-8 (log-posterior route) against values up to
  # 0.17. A constant bias in either route, such as a wrong baseline at the
  # mode, would survive the extrapolation.
  expect_lt(max(abs((g_grad - g_ref) / (g_grad_half - g_ref) - 2)), 0.01)
  expect_lt(max(abs((g_lp - g_ref) / (g_lp_half - g_ref) - 2)), 0.01)
  expect_lt(max(abs(2 * g_grad_half - g_grad - g_ref)), 5e-7)
  expect_lt(max(abs(2 * g_lp_half - g_lp - g_ref)), 3e-7)
})

test_that("the volume correction vanishes for an exactly Gaussian log-posterior", {
  mod <- gaussian_model()
  ctrl <- vanilla_control(list())
  for (gr in list(mod$grad, NULL)) {
    lap <- laplace_fit(mod$logpost, gr, mod$init, ctrl)
    sds <- sqrt(diag(lap$Sigma))
    v_scan <- sweep(lap$Sigma, 2, sds, "/")
    tr0 <- whitened_trace_lp(mod$logpost, lap$mode, lap$lp_max, lap$L)
    for (j in 1:3) {
      expect_lt(
        abs(gamma1_grad(j, mod$grad, lap$mode, v_scan, lap$L)),
        1e-8 # observed 1.7e-9, rounding in differences of the gradient
      )
      expect_lt(
        abs(gamma1_lp(
          j,
          mod$logpost,
          lap$mode,
          lap$lp_max,
          v_scan,
          lap$L,
          tr0
        )),
        1e-10 # observed 1e-12
      )
    }
  }

  fit_grad <- inla_vanilla(
    mod$logpost,
    mod$init,
    grad = mod$grad,
    silent = TRUE
  )
  fit_fd <- inla_vanilla(mod$logpost, mod$init, silent = TRUE)
  expect_lt(max(abs(fit_grad$diagnostics$marginals$gamma1)), 1e-8)
  expect_lt(max(abs(fit_fd$diagnostics$marginals$gamma1)), 1e-10)
})

# (c) Finite differences -------------------------------------------------------

# A smooth function of three coordinates, the third on a scale of hundreds,
# with quartic, cross and exponential terms, and its analytic derivatives.
smooth_fn <- function(x) {
  u <- x[3] / 100
  e <- exp(0.3 * x[1] - 0.2 * x[2] + 0.5 * u)
  0.25 * x[1]^4 + 0.5 * x[1]^2 * x[2] - 0.1 * x[2]^4 + e + u^4 + x[1] * u
}

smooth_grad <- function(x) {
  u <- x[3] / 100
  e <- exp(0.3 * x[1] - 0.2 * x[2] + 0.5 * u)
  c(
    x[1]^3 + x[1] * x[2] + 0.3 * e + u,
    0.5 * x[1]^2 - 0.4 * x[2]^3 - 0.2 * e,
    (0.5 * e + 4 * u^3 + x[1]) / 100
  )
}

smooth_hessian <- function(x) {
  u <- x[3] / 100
  e <- exp(0.3 * x[1] - 0.2 * x[2] + 0.5 * u)
  h13 <- (0.15 * e + 1) / 100
  h23 <- -0.1 * e / 100
  matrix(
    c(
      3 * x[1]^2 + x[2] + 0.09 * e,
      x[1] - 0.06 * e,
      h13,
      x[1] - 0.06 * e,
      -1.2 * x[2]^2 + 0.04 * e,
      h23,
      h13,
      h23,
      (0.25 * e + 12 * u^2) / 1e4
    ),
    3,
    3
  )
}

test_that("finite differences match the analytic derivatives of a smooth function", {
  x <- c(0.7, -0.4, 150)
  scale <- c(1, 1, 100)

  # Errors in units of the coordinate scales: derivatives times `scale`
  g_err <- (fd_grad(smooth_fn, x, h = 1e-5 * scale) - smooth_grad(x)) * scale
  expect_lt(max(abs(g_err)), 2e-9) # observed 7.0e-10

  H <- smooth_hessian(x)
  s2 <- outer(scale, scale)
  J <- fd_jacobian(smooth_grad, x, h = 1e-5 * scale)
  expect_lt(max(abs((J - H) * s2)), 1.5e-9) # observed 5.5e-10

  H_fd <- fd_hessian(smooth_fn, x, h = 1e-3 * scale)
  hess_err <- max(abs((H_fd - H) * s2))
  expect_lt(hess_err, 6e-6) # observed 2.0e-6
  expect_identical(H_fd, t(H_fd))
  # Second order: steps ten times larger give an error 100 times larger
  # (observed 100.0; rounding takes over below steps of about 2e-4)
  H_fd_coarse <- fd_hessian(smooth_fn, x, h = 1e-2 * scale)
  ratio <- max(abs((H_fd_coarse - H) * s2)) / hess_err
  expect_gt(ratio, 80)
  expect_lt(ratio, 120)
})

test_that("finite-difference errors take their exact form on polynomials, step by step", {
  # Distinct steps per coordinate, so a step applied to the wrong coordinate
  # changes the error. On polynomials the truncation errors are exact:
  #   central difference of c x^3:     3 c x^2 + c h^2,
  #   second difference of c x^4:     12 c x^2 + 2 c h^2,
  #   mixed difference of d x1^2 x2^2: 4 d x1 x2 + d h1 h2,
  # and a bilinear term e x2 x3 is differentiated exactly.
  x <- c(0.6, -1.2, 0.9)
  h <- c(0.1, 0.02, 0.3)
  cc <- c(1, -2, 0.5)
  d <- 0.7
  e <- -1.3

  cubic <- function(x) sum(cc * x^3)
  expect_equal(
    fd_grad(cubic, x, h),
    3 * cc * x^2 + cc * h^2,
    tolerance = 1e-12
  )

  # Jacobian of the gradient of sum(cc * x^4)
  quartic_grad <- function(x) 4 * cc * x^3
  expect_equal(
    fd_jacobian(quartic_grad, x, h),
    diag(12 * cc * x^2 + 4 * cc * h^2),
    tolerance = 1e-12
  )

  poly4 <- function(x) sum(cc * x^4) + d * x[1]^2 * x[2]^2 + e * x[2] * x[3]
  exact <- diag(12 * cc * x^2 + c(2 * d * x[2]^2, 2 * d * x[1]^2, 0))
  exact[1, 2] <- exact[2, 1] <- 4 * d * x[1] * x[2]
  exact[2, 3] <- exact[3, 2] <- e
  truncation <- diag(2 * cc * h^2)
  truncation[1, 2] <- truncation[2, 1] <- d * h[1] * h[2]
  expect_equal(
    fd_hessian(poly4, x, h),
    exact + truncation,
    tolerance = 1e-10
  )
})

test_that("fd_jacobian() puts the derivatives with respect to x_j in column j", {
  # A map whose Jacobian is far from symmetric
  fn <- function(x) c(x[1] * x[2], x[2]^3, exp(x[1]) + 5 * x[3])
  x <- c(0.4, 1.5, -2)
  jac <- rbind(
    c(x[2], x[1], 0),
    c(0, 3 * x[2]^2, 0),
    c(exp(x[1]), 0, 5)
  )
  expect_equal(fd_jacobian(fn, x, h = 1e-6), jac, tolerance = 1e-8)
})

test_that("finite differences recycle a scalar step and make the documented number of calls", {
  x <- c(0.7, -0.4, 150)
  expect_identical(
    fd_grad(smooth_fn, x, h = 1e-4),
    fd_grad(smooth_fn, x, h = rep(1e-4, 3))
  )
  expect_identical(
    fd_hessian(smooth_fn, x, h = 1e-3),
    fd_hessian(smooth_fn, x, h = rep(1e-3, 3))
  )

  for (m in 1:5) {
    calls <- 0L
    counted <- function(x) {
      calls <<- calls + 1L
      sum(exp(x / 2)) + prod(x)
    }
    counted_grad <- function(x) {
      calls <<- calls + 1L
      exp(x / 2) / 2
    }
    x <- seq_len(m) / 10

    calls <- 0L
    H <- fd_hessian(counted, x, h = 1e-3)
    expect_identical(calls, as.integer(m^2 + m + 1))
    expect_identical(H, t(H))

    f0 <- counted(x)
    calls <- 0L
    expect_identical(fd_hessian(counted, x, h = 1e-3, f0 = f0), H)
    expect_identical(calls, as.integer(m^2 + m))

    calls <- 0L
    fd_grad(counted, x, h = 1e-4)
    expect_identical(calls, 2L * m)

    calls <- 0L
    fd_jacobian(counted_grad, x, h = 1e-4)
    expect_identical(calls, 2L * m)
  }
})

# (d) The degree-5 cubature rule -----------------------------------------------

# E[prod z_i^k_i] under N(0, I): the product of (k_i - 1)!! over the k_i,
# or zero when any k_i is odd.
gaussian_moment <- function(k) {
  if (any(k %% 2 == 1)) {
    return(0)
  }
  prod(vapply(k, function(ki) prod(seq(1, max(ki - 1, 1), by = 2)), 1))
}

rule_moment <- function(z, weights, k) {
  mono <- rep(1, nrow(z))
  for (i in seq_along(k)) {
    mono <- mono * z[, i]^k[i]
  }
  sum(weights * mono)
}

# A random covariance matrix with scales spread over two orders of magnitude
random_covariance <- function(m) {
  A <- matrix(stats::rnorm(m * m), m)
  S <- crossprod(A) + diag(0.1, m)
  D <- diag(exp(stats::runif(m, -2.5, 2.5)), m)
  D %*% S %*% D
}

# Every exponent vector of m variables with total degree <= max_degree, built
# up one variable at a time
exponents_up_to <- function(m, max_degree) {
  if (m == 1) {
    return(matrix(0:max_degree, ncol = 1))
  }
  do.call(
    rbind,
    lapply(0:max_degree, function(k) {
      cbind(k, exponents_up_to(m - 1, max_degree - k))
    })
  )
}

test_that("vb_nodes() has (m + 1)(m + 2) + 1 nodes in antipodal pairs, within 3 SDs", {
  # Fewer for m <= 3, where some directions coincide, and at m = 7, where the
  # vertex nodes have zero weight
  n_nodes <- c(3, 7, 15, 31, 43, 57, 57, 91, 111)
  expect_identical(n_nodes[4:9][-4], ((4:9 + 1) * (4:9 + 2) + 1)[-4])
  set.seed(7)
  for (m in 1:9) {
    Sigma <- random_covariance(m)
    rule <- vb_nodes(Sigma)
    expect_identical(dim(rule$nodes), as.integer(c(n_nodes[m], m)))
    expect_equal(sum(rule$weights), 1, tolerance = 1e-14)
    # The centre, then one node of each antipodal pair, then the negatives
    n_pairs <- (n_nodes[m] - 1) / 2
    plus <- 1 + seq_len(n_pairs)
    expect_identical(rule$nodes[1, ], rep(0, m))
    expect_identical(rule$nodes[plus + n_pairs, ], -rule$nodes[plus, ])
    expect_identical(rule$weights[plus + n_pairs], rule$weights[plus])
    # In whitened coordinates no node is more than 3 SDs out
    z <- rule$nodes %*% solve(chol(Sigma))
    expect_lte(max(sqrt(rowSums(z^2))), 3 + 1e-12)
  }
})

test_that("vb_nodes() integrates every Gaussian moment up to degree five exactly", {
  # Every weight is positive up to m = 7, and the vertex weights are negative
  # from m = 8, so both regimes are covered
  set.seed(8)
  for (m in c(1:9, 12)) {
    Sigma <- random_covariance(m)
    rule <- vb_nodes(Sigma)
    expect_identical(min(rule$weights) > 0, m <= 7)
    z <- rule$nodes %*% solve(chol(Sigma))
    exps <- exponents_up_to(m, 5)
    err <- vapply(
      seq_len(nrow(exps)),
      function(r) {
        k <- exps[r, ]
        rule_moment(z, rule$weights, k) - gaussian_moment(k)
      },
      1
    )
    # Rounding, which grows with the condition number of Sigma: observed up
    # to 3.4e-11, at m = 4 with condition number 5e4
    expect_lt(max(abs(err)), 1e-10)
    expect_equal(
      crossprod(rule$nodes, rule$weights * rule$nodes),
      Sigma,
      tolerance = 1e-12
    )
    # Degree six is beyond the rule
    expect_gt(abs(rule_moment(z, rule$weights, c(6, rep(0, m - 1))) - 15), 1)
  }
})

test_that("vb_nodes() gives the classical rules in one, two and three dimensions", {
  # m = 1: the three-point Gauss-Hermite rule, as vb_nodes() is
  rule <- vb_nodes(matrix(2.5^2))
  ord <- order(rule$nodes[, 1])
  expect_equal(
    rule$nodes[ord, 1],
    2.5 * c(-sqrt(3), 0, sqrt(3)),
    tolerance = 1e-14
  )
  expect_equal(rule$weights[ord], c(1, 4, 1) / 6, tolerance = 1e-14)

  # m = 2: Radon's seven-point rule, the centre with weight 1/2 and a regular
  # hexagon of radius 2 with weights 1/12
  rule <- vb_nodes(diag(2))
  expect_equal(sqrt(rowSums(rule$nodes^2)), c(0, rep(2, 6)), tolerance = 1e-14)
  expect_equal(rule$weights, c(1 / 2, rep(1 / 12, 6)), tolerance = 1e-14)
  angles <- sort(atan2(rule$nodes[-1, 2], rule$nodes[-1, 1]) %% (2 * pi))
  expect_equal(diff(angles), rep(pi / 3, 5), tolerance = 1e-12)

  # m = 3: the corners of a cube and the vertices of an octahedron, all at
  # radius sqrt(5), with weights 0.045 and 0.04, and 0.4 at the centre
  rule <- vb_nodes(diag(3))
  expect_equal(
    sqrt(rowSums(rule$nodes^2)),
    c(0, rep(sqrt(5), 14)),
    tolerance = 1e-14
  )
  expect_equal(rule$weights[1], 0.4, tolerance = 1e-14)
  expect_equal(sort(rule$weights[-1]), rep(c(0.04, 0.045), c(6, 8)))
  # Octahedron vertices are orthogonal or opposite to each other, cube
  # corners are at angles with cosine +-1/3 or opposite
  oct <- rule$nodes[-1, ][abs(rule$weights[-1] - 0.04) < 1e-12, ] / sqrt(5)
  cube <- rule$nodes[-1, ][abs(rule$weights[-1] - 0.045) < 1e-12, ] / sqrt(5)
  cos_oct <- tcrossprod(oct)[upper.tri(diag(6))]
  cos_cube <- tcrossprod(cube)[upper.tri(diag(8))]
  expect_true(all(abs(cos_oct) < 1e-12 | abs(abs(cos_oct) - 1) < 1e-12))
  expect_true(all(
    abs(abs(cos_cube) - 1 / 3) < 1e-12 | abs(abs(cos_cube) - 1) < 1e-12
  ))
})

# (e) The VB shift ------------------------------------------------------------

test_that("the Gauss-Hermite reference is exact where it can be checked", {
  # Golub-Welsch nodes reproduce N(0, 1) moments to degree 2n - 1
  gh <- gauss_hermite(5)
  expect_equal(sum(gh$weights), 1, tolerance = 1e-14)
  expect_equal(sum(gh$weights * gh$nodes^4), 3, tolerance = 1e-13)
  expect_equal(sum(gh$weights * gh$nodes^8), 105, tolerance = 1e-13)

  # The reference VB mean of the Poisson-gamma model is closed form
  mod <- poisson_gamma_model(y = c(2, 0, 1), a = 0.5, b = 1)
  s2 <- 1 / mod$A
  ref <- vb_reference(mod$grad, log(mod$A / mod$B), matrix(s2))
  expect_equal(ref[[1]], mod$vb_exact(sqrt(s2)), tolerance = 1e-12)
})

test_that("the VB shift is zero for an exactly Gaussian posterior", {
  mod <- gaussian_model()
  sds <- sqrt(diag(solve(mod$precision)))

  # At the exact mode and covariance the shift vanishes to rounding
  R <- chol(mod$precision)
  vb <- vb_correct(
    mod$logpost,
    mod$grad,
    mod$mean,
    chol2inv(R),
    R,
    backsolve(R, diag(3)),
    vanilla_control(list())
  )
  expect_lt(max(abs(vb$shift) / sds), 1e-14) # observed 3e-17
  expect_false(vb$fallback)
  expect_lt(vb$kld_global, 1e-28)

  # In a full fit the shift only absorbs where nlminb() stopped, so the VB
  # mean is the exact mean, and the log marginal likelihood is exact.
  for (gr in list(mod$grad, NULL)) {
    fit <- inla_vanilla(mod$logpost, mod$init, grad = gr, silent = TRUE)
    expect_lt(max_sd_diff(fit$mean_vb, fit$mode, sds), 1e-8) # observed 3e-10
    expect_lt(max_sd_diff(fit$mean_vb, mod$mean, sds), 1e-12) # observed 3e-14
    expect_equal(fit$logml, mod$logml, tolerance = 1e-10)
    expect_equal(fit$logml, fit$logml_laplace, tolerance = 1e-14)
    expect_false(fit$diagnostics$vb$fallback)
  }
})

test_that("the VB shift matches a Gauss-Hermite reference for a skewed 1D posterior", {
  # Beta-binomial on the logit scale, y = 2 of n = 12, Beta(1, 1) prior
  A <- 3
  B <- 11
  logpost <- function(theta) {
    A *
      stats::plogis(theta[[1]], log.p = TRUE) +
      B * stats::plogis(-theta[[1]], log.p = TRUE)
  }
  grad <- function(theta) A - (A + B) * stats::plogis(theta[[1]])
  fit <- inla_vanilla(logpost, c(t = 0), grad = grad, silent = TRUE)
  sds <- sqrt(fit$cov[1, 1])
  ref <- vb_reference(grad, fit$mode, fit$cov)
  # The reference is converged in the number of nodes
  expect_equal(
    vb_reference(grad, fit$mode, fit$cov, n_gh = 60),
    ref,
    tolerance = 1e-12
  )

  vb_err <- max_sd_diff(fit$mean_vb, ref, sds)
  laplace_err <- max_sd_diff(fit$mode, ref, sds)
  expect_lt(vb_err, 1.5e-3) # observed 5.3e-4 SDs, for a shift of 0.17 SDs
  expect_lt(vb_err, laplace_err / 100)

  # Without a gradient: Stein's identity on the three-point rule, exact for
  # log-densities of degree four rather than six
  fit_fd <- inla_vanilla(logpost, c(t = 0), silent = TRUE)
  vb_err_fd <- max_sd_diff(fit_fd$mean_vb, ref, sds)
  expect_lt(vb_err_fd, 0.02) # observed 6.6e-3 SDs
  expect_lt(vb_err_fd, laplace_err / 10)
})

test_that("the VB shift solves the three-point rule exactly for a Poisson-gamma posterior", {
  # In 1D the degree-5 rule is the three-point Gauss-Hermite rule, so the shift
  # has a closed form under the rule as well as under exact VB.
  mod <- poisson_gamma_model(y = c(2, 0, 1), a = 0.5, b = 1)
  fit <- inla_vanilla(
    mod$logpost,
    mod$init,
    grad = mod$grad,
    silent = TRUE,
    control = list(vb_tol = 1e-12)
  )
  s <- sqrt(fit$cov[1, 1])
  expect_equal(s^2, 1 / mod$A, tolerance = 1e-8)
  expect_lt(abs(fit$mean_vb[[1]] - mod$vb_rule(s)) / s, 1e-10)
  # and the rule's answer is close to exact VB: observed 3.4e-4 SDs, against
  # a shift of 0.27 SDs from the mode
  expect_lt(abs(fit$mean_vb[[1]] - mod$vb_exact(s)) / s, 1e-3)
})

test_that("the VB shift matches a Gauss-Hermite reference for a skewed 2D posterior", {
  mod <- poisson_regression_model(2)
  fit <- inla_vanilla(mod$logpost, mod$init, grad = mod$grad, silent = TRUE)
  sds <- sqrt(diag(fit$cov))
  expect_lt(stats::cov2cor(fit$cov)[1, 2], -0.5) # correlated parameters
  # The reference holds q's covariance at fit$cov, so check that against the
  # inverse of the analytic Hessian at the mode first (observed 1.6e-10
  # relative to the SDs): the reference must not inherit an error from it.
  sigma_exact <- solve(-mod$hessian(fit$mode))
  expect_lt(max(abs(fit$cov - sigma_exact) / outer(sds, sds)), 5e-10)
  ref <- vb_reference(mod$grad, fit$mode, fit$cov)

  vb_err <- max_sd_diff(fit$mean_vb, ref, sds)
  laplace_err <- max_sd_diff(fit$mode, ref, sds)
  expect_lt(vb_err, 7.5e-4) # observed 2.6e-4 SDs, for a shift of 0.25 SDs
  expect_lt(vb_err, laplace_err / 100)

  # Without a gradient, by Stein's identity on the simplex rule
  fit_fd <- inla_vanilla(mod$logpost, mod$init, silent = TRUE)
  vb_err_fd <- max_sd_diff(fit_fd$mean_vb, ref, sds)
  expect_lt(vb_err_fd, 0.015) # observed 4.6e-3 SDs
  expect_lt(vb_err_fd, laplace_err / 10)
})

# inla_vanilla() with only the gradient check at init silenced. That check
# steps 1e-5 whatever the scale of the parameter, so it warns about a correct
# gradient when the posterior SDs are small (see the test after the next
# one). Every other warning still gets through: the fit runs at the default
# output level, since silent = TRUE would muffle them, with its stage lines
# suppressed.
fit_past_init_check <- function(logpost, init, grad) {
  withCallingHandlers(
    suppressMessages(inla_vanilla(logpost, init, grad = grad)),
    inlavanilla_warning = function(w) {
      if (grepl("does not match finite differences", conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

test_that("the VB shift in posterior SDs does not depend on the scale of the parameters", {
  # theta = scale * beta for the 2D Poisson regression, so the posterior SDs
  # run from 3e-5 to 4e3. The fixed-point step, its convergence test and the
  # fallback trigger are all measured in posterior SDs, and the
  # finite-difference steps scale with them, so the shift in SDs should not
  # change and the iteration should not fall back.
  mod <- poisson_regression_model(2)
  fit <- inla_vanilla(mod$logpost, mod$init, grad = mod$grad, silent = TRUE)
  fit_fd <- inla_vanilla(mod$logpost, mod$init, silent = TRUE)
  shift_sd <- function(f) (f$mean_vb - f$mode) / sqrt(diag(f$cov))
  ref_sd <- (vb_reference(mod$grad, fit$mode, fit$cov) - fit$mode) /
    sqrt(diag(fit$cov))

  for (scale in list(c(1e-4, 1e-4), c(1e4, 1e4), c(1e-4, 1e4))) {
    lp_scaled <- function(theta) mod$logpost(theta / scale)
    gr_scaled <- function(theta) mod$grad(theta / scale) / scale
    for (gr in list(gr_scaled, NULL)) {
      fit_scaled <- fit_past_init_check(lp_scaled, mod$init, grad = gr)
      expect_false(fit_scaled$diagnostics$vb$fallback)
      # Against the exact VB shift, as for the unscaled model: observed 2.0e-4
      # to 2.6e-4 SDs with the gradient, and 4.2e-3 to 4.6e-3 without, from
      # Stein's identity
      tol <- if (is.null(gr)) 0.015 else 7.5e-4
      expect_lt(max(abs(shift_sd(fit_scaled) - ref_sd)), tol)
      if (scale[[1]] == scale[[2]]) {
        # A common scale changes nothing but rounding: observed 7.3e-7 SDs
        # with the gradient, 1.6e-7 without
        base <- if (is.null(gr)) fit_fd else fit
        expect_lt(max(abs(shift_sd(fit_scaled) - shift_sd(base))), 2e-6)
        expect_identical(
          fit_scaled$diagnostics$vb$iterations,
          base$diagnostics$vb$iterations
        )
      }
      # With scales 1e8 apart the shift moves by 9.6e-5 SDs and takes one
      # more iteration: the Anderson secant uses inner products of the raw
      # steps, so the path, and where it stops within vb_tol = 1e-3 SDs,
      # depends on the relative scales. Hence no tighter check there.
    }
  }
})

test_that("a correct gradient passes the check at init when the posterior SD is small", {
  # The Poisson-gamma posterior on theta = 1e-4 log(lambda), whose posterior
  # SD is 5.3e-5, with its exact gradient. The check at init takes central
  # differences with a step of 1e-5 (0.19 posterior SDs here) whatever the
  # scale, and their truncation error, 1.3%, exceeds its threshold of 0.1%,
  # so it warns that a correct gradient is wrong ("grad gives -5000, finite
  # differences give -5067"). The Laplace stage handles this scale (see the
  # test above); the check should too.
  mod <- poisson_gamma_model(y = c(2, 0, 1), a = 0.5, b = 1)
  lp_scaled <- function(theta) mod$logpost(theta / 1e-4)
  gr_scaled <- function(theta) mod$grad(theta / 1e-4) / 1e-4
  expect_equal(
    gr_scaled(c(t = 2e-5)),
    (lp_scaled(c(t = 2e-5 + 1e-9)) - lp_scaled(c(t = 2e-5 - 1e-9))) / 2e-9,
    tolerance = 1e-6
  )
  # silent = TRUE muffles warnings, but records them in fit$warnings
  fit <- inla_vanilla(lp_scaled, c(t = 0), grad = gr_scaled, silent = TRUE)
  expect_identical(fit$warnings, character())
})

test_that("the degree-5 rule gives the exact VB shift for a degree-6 polynomial log-posterior", {
  # A Gaussian plus cubic, quartic and sextic terms, bounded above. Its gradient
  # is a polynomial of degree five, which the rule integrates exactly against
  # any Gaussian, so the shift must match the reference to rounding error.
  P <- matrix(c(1.5, 0.4, 0.2, 0.4, 1, -0.3, 0.2, -0.3, 0.8), 3)
  logpost <- function(t) {
    -0.5 *
      sum(t * (P %*% t)) +
      0.15 * t[[1]]^3 -
      0.05 * t[[1]]^4 -
      0.03 * t[[1]]^2 * t[[2]]^2 +
      0.02 * t[[1]] * t[[2]]^2 +
      0.04 * t[[2]] * t[[3]]^2 -
      0.002 * t[[1]]^6 -
      0.01 * t[[3]]^4
  }
  grad <- function(t) {
    t1 <- t[[1]]
    t2 <- t[[2]]
    t3 <- t[[3]]
    -drop(P %*% t) +
      c(
        0.45 *
          t1^2 -
          0.2 * t1^3 -
          0.06 * t1 * t2^2 +
          0.02 * t2^2 -
          0.012 * t1^5,
        -0.06 * t1^2 * t2 + 0.04 * t1 * t2 + 0.04 * t3^2,
        0.08 * t2 * t3 - 0.04 * t3^3
      )
  }
  init <- c(a = 0, b = 0, c = 0)
  expect_equal(
    grad(c(0.3, -0.7, 1.1)),
    fd_grad(logpost, c(0.3, -0.7, 1.1), 1e-6)
  )

  fit <- inla_vanilla(
    logpost,
    init,
    grad = grad,
    silent = TRUE,
    control = list(vb_tol = 1e-12, vb_maxit = 200)
  )
  expect_false(fit$diagnostics$vb$fallback)
  sds <- sqrt(diag(fit$cov))
  # Eight nodes per axis are exact up to degree 15
  ref <- vb_reference(grad, fit$mode, fit$cov, n_gh = 8)
  expect_gt(max_sd_diff(ref, fit$mode, sds), 0.2) # a real shift
  # Solved to vb_tol = 1e-12 SDs with an exact rule: observed 4e-14 SDs
  expect_lt(max_sd_diff(fit$mean_vb, ref, sds), 1e-12)

  # With the default tolerance the iteration stops a little short
  fit_default <- inla_vanilla(logpost, init, grad = grad, silent = TRUE)
  expect_lt(max_sd_diff(fit_default$mean_vb, ref, sds), 1e-4) # observed 2e-5
})

test_that("without a gradient the VB shift is exact for a degree-4 polynomial log-posterior", {
  # The example above without its sextic term. Without a gradient the
  # expected gradient comes from Stein's identity on the simplex rule, which
  # integrates z f(mu + z) exactly when f is a polynomial of degree four.
  P <- matrix(c(1.5, 0.4, 0.2, 0.4, 1, -0.3, 0.2, -0.3, 0.8), 3)
  logpost <- function(t) {
    -0.5 *
      sum(t * (P %*% t)) +
      0.15 * t[[1]]^3 -
      0.05 * t[[1]]^4 -
      0.03 * t[[1]]^2 * t[[2]]^2 +
      0.02 * t[[1]] * t[[2]]^2 +
      0.04 * t[[2]] * t[[3]]^2 -
      0.01 * t[[3]]^4
  }
  grad <- function(t) {
    t1 <- t[[1]]
    t2 <- t[[2]]
    t3 <- t[[3]]
    -drop(P %*% t) +
      c(
        0.45 * t1^2 - 0.2 * t1^3 - 0.06 * t1 * t2^2 + 0.02 * t2^2,
        -0.06 * t1^2 * t2 + 0.04 * t1 * t2 + 0.04 * t3^2,
        0.08 * t2 * t3 - 0.04 * t3^3
      )
  }
  init <- c(a = 0, b = 0, c = 0)
  expect_equal(
    grad(c(0.3, -0.7, 1.1)),
    fd_grad(logpost, c(0.3, -0.7, 1.1), 1e-6)
  )

  # The gradient serves only the reference
  ctrl <- list(vb_tol = 1e-12, vb_maxit = 200)
  fit <- inla_vanilla(logpost, init, silent = TRUE, control = ctrl)
  expect_false(fit$has_grad)
  expect_false(fit$diagnostics$vb$fallback)
  sds <- sqrt(diag(fit$cov))
  ref <- vb_reference(grad, fit$mode, fit$cov, n_gh = 8)
  expect_gt(max_sd_diff(ref, fit$mode, sds), 0.2) # a real shift
  expect_lt(max_sd_diff(fit$mean_vb, ref, sds), 1e-11) # observed 2.1e-13 SDs

  # A sextic term is beyond it, though the gradient rule takes it exactly
  # (see the test above)
  logpost6 <- function(t) logpost(t) - 0.002 * t[[1]]^6
  grad6 <- function(t) grad(t) - c(0.012 * t[[1]]^5, 0, 0)
  fit6 <- inla_vanilla(logpost6, init, silent = TRUE, control = ctrl)
  sds6 <- sqrt(diag(fit6$cov))
  ref6 <- vb_reference(grad6, fit6$mode, fit6$cov, n_gh = 8)
  expect_gt(max_sd_diff(fit6$mean_vb, ref6, sds6), 1e-3) # observed 5.4e-3
})

test_that("without a gradient each VB sweep costs one log-posterior call per node", {
  mod <- logistic_model()
  fit <- inla_vanilla(mod$logpost, mod$init, silent = TRUE)
  m <- length(fit$mode)
  R <- chol(solve(fit$cov))
  calls <- 0L
  counted <- function(x) {
    calls <<- calls + 1L
    mod$logpost(x)
  }
  vb <- vb_correct(
    counted,
    NULL,
    fit$mode,
    fit$cov,
    R,
    backsolve(R, diag(m)),
    vanilla_control(list())
  )
  expect_false(vb$fallback)
  # (m + 1)(m + 2) = 42 calls per sweep, no centre: observed 4 sweeps
  n_nodes <- nrow(vb_nodes(fit$cov)$nodes) - 1L
  expect_identical(n_nodes, as.integer((m + 1) * (m + 2)))
  expect_identical(calls, vb$iterations * n_nodes)
  expect_equal(vb$shift, fit$mean_vb - fit$mode, tolerance = 1e-12)
})

# (f) vb = FALSE against vb = TRUE ---------------------------------------------

test_that("vb = FALSE keeps the mean at the mode and the Laplace marginal likelihood", {
  mod <- poisson_regression_model(2)
  fit <- inla_vanilla(
    mod$logpost,
    mod$init,
    grad = mod$grad,
    vb = FALSE,
    silent = TRUE
  )
  expect_false(fit$vb)
  expect_identical(fit$mean_vb, fit$mode)
  expect_identical(fit$logml, fit$logml_laplace)
  expect_null(fit$diagnostics$vb)
  expect_null(summary(fit)$vb)
  # The VB diagnostics did not run, so they are NA, as in INLAvaan's
  # diagnostics(); the mean gap needs a VB mean too
  d <- fit$diagnostics$marginals
  for (col in c(
    "kld",
    "vb_shift",
    "vb_shift_sigma",
    "vb_step_sigma",
    "mean_gap"
  )) {
    expect_identical(d[[col]], c(NA_real_, NA_real_))
  }
  g <- diagnostics(fit)
  expect_identical(g[["vb_applied"]], 0)
  for (nm in c(
    "vb_shift_max",
    "vb_kld_global",
    "kld_max",
    "vb_fallback",
    "vb_step_max"
  )) {
    expect_identical(g[[nm]], NA_real_)
  }
})

test_that("vb = TRUE adds 0.5 shift' H shift to the log marginal likelihood", {
  mod <- poisson_regression_model(2)
  fit_vb <- inla_vanilla(
    mod$logpost,
    mod$init,
    grad = mod$grad,
    silent = TRUE
  )
  fit_lap <- inla_vanilla(
    mod$logpost,
    mod$init,
    grad = mod$grad,
    vb = FALSE,
    silent = TRUE
  )
  shift <- fit_vb$mean_vb - fit_vb$mode
  H <- solve(fit_vb$cov)
  kld <- 0.5 * drop(t(shift) %*% H %*% shift)
  expect_gt(kld, 0.01) # observed 0.033
  expect_equal(fit_vb$logml - fit_vb$logml_laplace, kld, tolerance = 1e-10)
  expect_equal(fit_vb$diagnostics$vb$kld_global, kld, tolerance = 1e-10)
  expect_identical(fit_vb$logml_laplace, fit_lap$logml)

  # The shift on the scale of logpost, and in posterior SDs, and the KL
  # divergence of each shifted marginal, as INLAvaan's diagnostics() has them
  sds <- sqrt(diag(fit_vb$cov))
  d <- fit_vb$diagnostics$marginals
  expect_equal(d$vb_shift, unname(shift), tolerance = 1e-12)
  expect_equal(d$vb_shift_sigma, unname(shift / sds), tolerance = 1e-12)
  expect_equal(d$kld, 0.5 * d$vb_shift_sigma^2, tolerance = 1e-12)
  g <- diagnostics(fit_vb)
  expect_equal(g[["vb_kld_global"]], kld, tolerance = 1e-10)
  expect_equal(g[["kld_max"]], max(d$kld))
  expect_equal(g[["vb_shift_max"]], max(abs(d$vb_shift_sigma)))
  expect_equal(
    fit_vb$diagnostics$marginals$mean_gap,
    unname((fit_vb$summary[, "Mean"] - fit_vb$mean_vb) / sds),
    tolerance = 1e-12
  )

  # The marginals are scanned around the mode, so VB does not move them
  expect_identical(fit_vb$mode, fit_lap$mode)
  expect_identical(fit_vb$cov, fit_lap$cov)
  expect_identical(fit_vb$sn, fit_lap$sn)
  expect_identical(fit_vb$summary, fit_lap$summary)
})

test_that("the VB term of the log marginal likelihood is positive and uses the exact Hessian", {
  # The term is 0.5 shift' H shift, a KL divergence, so it is never negative.
  # Here H is the analytic Hessian at the mode rather than the package's
  # finite differences, on models that take the fixed point and the fallback
  # (the Poisson-gamma), with and without a gradient.
  models <- list(
    normal_model(),
    poisson_regression_model(2),
    poisson_regression_model(3),
    poisson_gamma_model(y = 0, a = 0.5, b = 1)
  )
  for (mod in models) {
    for (gr in list(mod$grad, NULL)) {
      fit <- inla_vanilla(mod$logpost, mod$init, grad = gr, silent = TRUE)
      shift <- fit$mean_vb - fit$mode
      kld <- 0.5 * drop(shift %*% -mod$hessian(fit$mode) %*% shift)
      expect_gt(kld, 0.03) # observed 0.033 to 0.23
      expect_gt(fit$logml - fit$logml_laplace, 0)
      # Observed relative errors up to 3.3e-9 with the gradient and 1.7e-5
      # without, where the Hessian comes from second differences
      expect_equal(
        fit$logml - fit$logml_laplace,
        kld,
        tolerance = if (is.null(gr)) 5e-5 else 1e-8
      )
    }
  }
})

# (g) The nlminb() fallback ---------------------------------------------------

test_that("a strongly skewed Poisson posterior falls back to nlminb() at a stationary point", {
  # One zero count and a Gamma(0.5, 1) prior on the rate: on the log scale
  # the posterior is log-gamma with shape 0.5 and the Laplace SD is sqrt(2).
  # The first Newton step is over one SD, which triggers the fallback.
  mod <- poisson_gamma_model(y = 0, a = 0.5, b = 1)
  fit <- inla_vanilla(mod$logpost, mod$init, grad = mod$grad, silent = TRUE)
  expect_true(fit$diagnostics$vb$fallback)
  s <- sqrt(fit$cov[1, 1])
  expect_equal(s, sqrt(2), tolerance = 1e-6)

  # A stationary point of the rule's expected log-posterior ...
  expect_lt(
    max(abs(rule_stationarity(mod$grad, fit$mean_vb, fit$cov))),
    6e-7 # observed 2.0e-7
  )
  expect_gt(max(abs(rule_stationarity(mod$grad, fit$mode, fit$cov))), 0.5)
  # ... which in 1D is closed form
  expect_lt(abs(fit$mean_vb[[1]] - mod$vb_rule(s)) / s, 6e-7) # observed 2.0e-7
  # and close to exact VB, unlike the mode: observed 0.028 SDs, against 0.71
  vb_err <- abs(fit$mean_vb[[1]] - mod$vb_exact(s)) / s
  expect_lt(vb_err, 0.06)
  expect_lt(vb_err, abs(fit$mode[[1]] - mod$vb_exact(s)) / s / 10)

  # The marginal likelihood still takes the VB term
  shift <- fit$mean_vb - fit$mode
  expect_equal(
    fit$logml - fit$logml_laplace,
    0.5 * sum(shift^2) / fit$cov[1, 1],
    tolerance = 1e-10
  )
})

test_that("a strongly skewed 2D posterior falls back to nlminb() at a stationary point", {
  # theta_1 as above, theta_2 | theta_1 ~ N(0.8 theta_1, 0.5^2). Exact VB has
  # mu_1 = log(A / B) - Sigma_11 / 2 and mu_2 = 0.8 mu_1.
  pg <- poisson_gamma_model(y = 0, a = 0.5, b = 1)
  logpost <- function(t) {
    pg$logpost(t[1]) + stats::dnorm(t[[2]], 0.8 * t[[1]], 0.5, log = TRUE)
  }
  grad <- function(t) {
    r <- (t[[2]] - 0.8 * t[[1]]) / 0.25
    c(pg$grad(t[1]) + 0.8 * r, -r)
  }
  init <- c(a = 0, b = 0)
  fit <- inla_vanilla(logpost, init, grad = grad, silent = TRUE)
  expect_true(fit$diagnostics$vb$fallback)
  sds <- sqrt(diag(fit$cov))
  expect_lt(max(abs(rule_stationarity(grad, fit$mean_vb, fit$cov))), 6e-7) # obs. 1.9e-7
  # theta_2 enters linearly, so mu_2 = 0.8 mu_1 holds under any rule
  expect_equal(fit$mean_vb[[2]], 0.8 * fit$mean_vb[[1]], tolerance = 1e-10)
  mu1_exact <- log(pg$A / pg$B) - fit$cov[1, 1] / 2
  mode_err <- abs(fit$mode[[1]] - mu1_exact) / sds[[1]]
  expect_lt(abs(fit$mean_vb[[1]] - mu1_exact) / sds[[1]], 0.06) # obs. 0.026

  # Without a gradient the expected gradient comes from Stein's identity on
  # the simplex rule, whose first step here is under one SD, so there is no
  # fallback. The rule is exact for log-densities of degree four only, and
  # this posterior is far from that: observed 0.096 SDs from exact VB, against
  # 0.026 with the gradient and 0.71 for the mode.
  fit_fd <- inla_vanilla(logpost, init, silent = TRUE)
  sds <- sqrt(diag(fit_fd$cov))
  vb_err_fd <- abs(fit_fd$mean_vb[[1]] - mu1_exact) / sds[[1]]
  expect_lt(vb_err_fd, 0.2)
  expect_lt(vb_err_fd, mode_err / 5)
  # mu_2 = 0.8 mu_1 to the tolerance of the fixed-point iteration, which
  # stops when a step is under 1e-3 SDs: observed 4.7e-3 SDs
  expect_lt(
    abs(fit_fd$mean_vb[[2]] - 0.8 * fit_fd$mean_vb[[1]]) / sds[[2]],
    0.015
  )
})

test_that("a forced fallback agrees with the fixed-point solution", {
  mod <- poisson_regression_model(2)
  fit_fixed <- inla_vanilla(
    mod$logpost,
    mod$init,
    grad = mod$grad,
    silent = TRUE,
    control = list(vb_tol = 1e-10)
  )
  # One iteration cannot converge, so the solve falls back to nlminb()
  fit_forced <- inla_vanilla(
    mod$logpost,
    mod$init,
    grad = mod$grad,
    silent = TRUE,
    control = list(vb_maxit = 1)
  )
  expect_false(fit_fixed$diagnostics$vb$fallback)
  expect_true(fit_forced$diagnostics$vb$fallback)
  sds <- sqrt(diag(fit_fixed$cov))
  expect_lt(
    max_sd_diff(fit_forced$mean_vb, fit_fixed$mean_vb, sds),
    2e-7 # observed 6.1e-8
  )
  # Both settled: the fixed point's last step was below vb_tol, and the step
  # where the fallback landed is as small as its distance from the fixed
  # point
  expect_lt(max(fit_fixed$diagnostics$marginals$vb_step_sigma), 1e-10)
  expect_lt(max(fit_forced$diagnostics$marginals$vb_step_sigma), 2e-7)
  dg <- diagnostics_internal(
    fit_forced$diagnostics$optim,
    fit_forced$diagnostics$vb,
    fit_forced$diagnostics$marginals
  )
  expect_false(check_flags(dg, gaussian = FALSE)[["vb_step"]])
  expect_lt(
    max(abs(rule_stationarity(mod$grad, fit_forced$mean_vb, fit_forced$cov))),
    2e-7 # observed 6.4e-8
  )
  # observed 4.5e-11 relative
  expect_equal(fit_forced$logml, fit_fixed$logml, tolerance = 1e-10)
})

test_that("the VB step is the step still to go at the VB mean, in posterior SDs", {
  # The Poisson regression with each coefficient scaled by 100, so that the
  # posterior SDs are far from 1 and a step left in the parameters' own units
  # would be off by that factor. The step is checked against the rule's own
  # fixed-point step, Sigma E[grad], taken here at the VB mean.
  mod <- poisson_regression_model(2)
  s <- 100
  lp_s <- function(theta) mod$logpost(theta / s)
  gr_s <- function(theta) mod$grad(theta / s) / s
  step_at_mean <- function(fit) {
    rule <- vb_nodes(fit$cov)
    g <- expected_grad(gr_s, fit$mean_vb, rule$nodes, rule$weights)
    unname(abs(drop(fit$cov %*% g)) / sqrt(diag(fit$cov)))
  }
  fit_s <- function(...) {
    inla_vanilla(lp_s, s * mod$init, grad = gr_s, silent = TRUE, ...)
  }

  # A vb_tol above the settling level is met by the first step, of more
  # than 0.1 SDs, which says little about the step still to go: that is
  # taken afresh, and the correction settles
  loose <- fit_s(control = list(vb_tol = 0.5))
  d <- loose$diagnostics$marginals
  expect_false(loose$diagnostics$vb$fallback)
  expect_identical(loose$diagnostics$vb$iterations, 1L)
  expect_gt(max(abs(d$vb_shift_sigma)), 0.1)
  expect_equal(d$vb_step_sigma, step_at_mean(loose), tolerance = 1e-6)
  expect_lt(max(d$vb_step_sigma), 0.1)
  expect_length(loose$warnings, 0)

  # After the fallback, the step where it landed
  forced <- fit_s(control = list(vb_maxit = 1))
  expect_true(forced$diagnostics$vb$fallback)
  expect_equal(
    forced$diagnostics$marginals$vb_step_sigma,
    step_at_mean(forced),
    tolerance = 1e-4
  )
})

test_that("without a gradient a forced fallback lands close to the fixed point", {
  # Without a gradient nlminb() minimises the rule's negative expected
  # log-posterior, with the expected gradient by Stein's identity as its
  # gradient. The two agree only to the error of the rule, so the fallback
  # stops near the fixed point of the iteration rather than on it.
  mod <- poisson_regression_model(2)
  fit_fixed <- inla_vanilla(
    mod$logpost,
    mod$init,
    silent = TRUE,
    control = list(vb_tol = 1e-10)
  )
  fit_forced <- inla_vanilla(
    mod$logpost,
    mod$init,
    silent = TRUE,
    control = list(vb_maxit = 1)
  )
  expect_false(fit_fixed$diagnostics$vb$fallback)
  expect_true(fit_forced$diagnostics$vb$fallback)
  sds <- sqrt(diag(fit_fixed$cov))
  expect_lt(
    max_sd_diff(fit_forced$mean_vb, fit_fixed$mean_vb, sds),
    1.5e-3 # observed 4.5e-4
  )
  # The step where the fallback landed measures that distance: one more
  # fixed-point step would close most of it, far below the 0.1 SD level
  step <- fit_forced$diagnostics$marginals$vb_step_sigma
  expect_lt(max(step), 1.5e-3)
  expect_gt(
    max(step),
    0.2 * max_sd_diff(fit_forced$mean_vb, fit_fixed$mean_vb, sds)
  )
  # and as close to exact VB as the iteration is (both observed 4e-3 SDs)
  ref <- vb_reference(mod$grad, fit_fixed$mode, fit_fixed$cov)
  expect_lt(max_sd_diff(fit_forced$mean_vb, ref, sds), 0.015)
  expect_lt(abs(fit_forced$logml - fit_fixed$logml), 5e-4) # observed 1.2e-4
})
