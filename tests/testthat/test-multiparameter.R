# End-to-end tests of inla_vanilla() on models with two or three parameters,
# against exact posteriors or numerical quadrature written here:
#
# (a) normal data with unknown mean and variance under the conjugate
#     normal-inverse-gamma prior, on (mu, log_sigma): Student-t and log
#     inverse-gamma marginals, closed-form marginal likelihood;
# (b) a bivariate normal with an unknown correlation, on (mu1, mu2, z) with
#     z = atanh(rho): the means are integrated out analytically and z's
#     marginal is found by 1D quadrature;
# (c) a Weibull model on (log_shape, log_scale) with normal priors, which has
#     no closed form: marginals by brute-force 2D grid quadrature;
# (d) invariance of the marginals to the order of the parameters.
#
# The Laplace step (mode, covariance, logml_laplace) is also checked against
# an exact or independently differenced Hessian in (a) to (c), and in (a) the
# VB mean against the exact VB fixed point, which has a closed form there.
#
# Every logpost carries all of its normalising constants and the Jacobian of
# the reparameterisation, so fit$logml is comparable with the exact log
# marginal likelihood. Location errors are in reference posterior SDs, SD
# errors are relative. Tolerances were calibrated by printing the actual
# errors first; each comment states the observed value.

probs <- c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975)
q_cols <- paste0(100 * probs, "%")

# Helpers ----------------------------------------------------------------------

# Trapezoid weights for an equally spaced grid of k points with spacing h
trapezoid_weights <- function(k, h) {
  h * c(0.5, rep(1, k - 2), 0.5)
}

# Central differences written here, so that the reference checks do not lean
# on the package's own finite differences (fd_grad(), fd_jacobian()).
num_grad <- function(fn, x, h = 1e-5) {
  vapply(
    seq_along(x),
    function(j) {
      e <- replace(numeric(length(x)), j, h)
      (fn(x + e) - fn(x - e)) / (2 * h)
    },
    numeric(1)
  )
}

num_jacobian <- function(gr, x, h = 1e-5) {
  vapply(
    seq_along(x),
    function(j) {
      e <- replace(numeric(length(x)), j, h)
      (gr(x + e) - gr(x - e)) / (2 * h)
    },
    numeric(length(x))
  )
}

# Checks of the Laplace step against an independent mode and Hessian: the
# exact gradient `grad_ref` at fit$mode in whitened units (how many posterior
# SDs the mode is off, to first order), the covariance against the inverse of
# `hess_ref` at fit$mode (relative to the SDs), and logml_laplace against
# lp(mode) + m/2 log(2 pi) + 1/2 log det(Sigma) from that covariance.
laplace_errors <- function(fit, logpost, grad_ref, hess_ref) {
  S <- solve(-hess_ref(fit$mode))
  sds <- sqrt(diag(S))
  m <- length(fit$mode)
  list(
    # chol(S) = U with U'U = S, so U g is the gradient in whitened units
    grad = max(abs(chol(S) %*% grad_ref(fit$mode))),
    cov = max(abs(fit$cov - S) / outer(sds, sds)),
    logml_laplace = fit$logml_laplace -
      (logpost(fit$mode) + m / 2 * log(2 * pi) + 0.5 * log(det(S)))
  )
}

# A reference marginal from log-density values `log_f` (up to a constant) on
# an equally spaced grid `x` that covers all of the mass: the normalised
# density, the mean and SD by the trapezoid rule, quantiles by linear
# interpolation of the trapezoid CDF, and the log of the normalising constant.
grid_marginal <- function(x, log_f) {
  h <- x[2] - x[1]
  w <- trapezoid_weights(length(x), h)
  f <- exp(log_f - max(log_f))
  norm_const <- sum(w * f)
  f <- f / norm_const
  cdf <- cumsum(c(0, (f[-1] + f[-length(f)]) / 2 * h))
  mu <- sum(w * x * f)
  list(
    x = x,
    density = f,
    log_norm = max(log_f) + log(norm_const),
    mean = mu,
    sd = sqrt(sum(w * (x - mu)^2 * f)),
    quantile = function(p) approx(cdf, x, p, ties = base::mean)$y,
    # the density at the ends of the grid, relative to its peak
    edge = max(f[c(1, length(f))]) / max(f)
  )
}

# Errors of the fitted marginal of parameter `par` against a reference
# marginal `ref` (a list with mean, sd and quantile()), for the skew-normal
# and for the Laplace Gaussian N(mode, cov[par, par]). Quantile errors are
# signed, at `probs`, in reference SDs.
marginal_errors <- function(fit, par, ref) {
  tab <- summary(fit, probs = probs)$table
  lap_sd <- sqrt(fit$cov[par, par])
  q_exact <- ref$quantile(probs)
  list(
    quantile_sn = (tab[par, q_cols] - q_exact) / ref$sd,
    quantile_laplace = (qnorm(probs, fit$mode[[par]], lap_sd) - q_exact) /
      ref$sd,
    mean_sn = (tab[par, "Mean"] - ref$mean) / ref$sd,
    mean_laplace = (fit$mode[[par]] - ref$mean) / ref$sd,
    mean_vb = (fit$mean_vb[[par]] - ref$mean) / ref$sd,
    sd_sn = tab[par, "SD"] / ref$sd - 1,
    sd_laplace = lap_sd / ref$sd - 1
  )
}

# Log density of N(0, V) at x
mvn_logdens <- function(x, V) {
  U <- chol(V)
  r <- backsolve(U, x, transpose = TRUE)
  -length(x) / 2 * log(2 * pi) - sum(log(diag(U))) - 0.5 * sum(r^2)
}

# log(1 - tanh(z)^2) = 2 log(sech z), stable for large |z|
log_sech2 <- function(z) {
  log(4) - 2 * abs(z) - 2 * log1p(exp(-2 * abs(z)))
}

# Models -----------------------------------------------------------------------

# (a) y_i ~ N(mu, sigma^2), mu | sigma^2 ~ N(m0, sigma^2 / k0),
# sigma^2 ~ IG(a0, b0), on theta = (mu, log_sigma). The Jacobian of
# sigma^2 = exp(2 log_sigma) is 2 sigma^2. With
#   kn = k0 + n, mn = (k0 m0 + n ybar) / kn, an = a0 + n / 2,
#   bn = b0 + SS / 2 + k0 n (ybar - m0)^2 / (2 kn),
# the posterior is mu | sigma^2 ~ N(mn, sigma^2 / kn), sigma^2 ~ IG(an, bn).
# Marginally mu = mn + sqrt(bn / (an kn)) T with T ~ t(2 an), and
# 1 / sigma^2 ~ Gamma(an, rate = bn), so log_sigma = -log(tau) / 2 has
# E = -(digamma(an) - log(bn)) / 2, Var = trigamma(an) / 4, mode
# log(bn / an) / 2, and P(log_sigma <= q) = P(tau >= exp(-2 q)).
#
# On theta the log-posterior is const - (2 an + 1) s - B(mu) exp(-2 s), with
# s = log_sigma and B(mu) = bn + kn (mu - mn)^2 / 2, so the joint mode is
# (mn, log(2 bn / (2 an + 1)) / 2) and the Hessian there is diagonal with
# entries kn (2 an + 1) / (2 bn) and 2 (2 an + 1). The marginal likelihood is
#   p(y) = (2 pi)^(-n/2) sqrt(k0 / kn) b0^a0 Gamma(an) / (Gamma(a0) bn^an).
#
# The VB mean correction keeps the Laplace covariance diag(v1, v2) and moves
# the mean to the m where the expected gradient under N(m, diag(v1, v2))
# vanishes. With E[exp(-2 s)] = exp(-2 m2 + 2 v2), that gives m1 = mn and
#   m2 = v2 + log((2 bn + kn v1) / (2 an + 1)) / 2.
nig_model <- function(y, m0, k0, a0, b0) {
  n <- length(y)
  ybar <- mean(y)
  kn <- k0 + n
  mn <- (k0 * m0 + n * ybar) / kn
  an <- a0 + n / 2
  bn <- b0 + 0.5 * sum((y - ybar)^2) + k0 * n * (ybar - m0)^2 / (2 * kn)
  nu <- 2 * an
  t_scale <- sqrt(bn / (an * kn))

  list(
    logpost = function(theta) {
      mu <- theta[["mu"]]
      s <- theta[["log_sigma"]]
      sigma2 <- exp(2 * s)
      sum(dnorm(y, mu, sqrt(sigma2), log = TRUE)) +
        dnorm(mu, m0, sqrt(sigma2 / k0), log = TRUE) +
        a0 * log(b0) -
        lgamma(a0) -
        (a0 + 1) * log(sigma2) -
        b0 / sigma2 +
        log(2) +
        2 * s
    },
    grad = function(theta) {
      mu <- theta[["mu"]]
      e2s <- exp(-2 * theta[["log_sigma"]])
      b_mu <- bn + kn / 2 * (mu - mn)^2
      c(-kn * (mu - mn) * e2s, -(2 * an + 1) + 2 * b_mu * e2s)
    },
    mu = list(
      density = function(x) dt((x - mn) / t_scale, nu) / t_scale,
      mean = mn,
      sd = t_scale * sqrt(nu / (nu - 2)),
      mode = mn,
      quantile = function(p) mn + t_scale * qt(p, nu)
    ),
    log_sigma = list(
      density = function(s) {
        dgamma(exp(-2 * s), an, rate = bn) * 2 * exp(-2 * s)
      },
      mean = -0.5 * (digamma(an) - log(bn)),
      sd = 0.5 * sqrt(trigamma(an)),
      mode = 0.5 * log(bn / an),
      quantile = function(p) {
        -0.5 * log(qgamma(p, an, rate = bn, lower.tail = FALSE))
      }
    ),
    joint_mode = c(mu = mn, log_sigma = 0.5 * log(2 * bn / (2 * an + 1))),
    joint_cov = diag(c(2 * bn / (kn * (2 * an + 1)), 1 / (2 * (2 * an + 1)))),
    vb_mean = local({
      v1 <- 2 * bn / (kn * (2 * an + 1))
      v2 <- 1 / (2 * (2 * an + 1))
      c(mu = mn, log_sigma = v2 + 0.5 * log((2 * bn + kn * v1) / (2 * an + 1)))
    }),
    logml = lgamma(an) -
      lgamma(a0) +
      a0 * log(b0) -
      an * log(bn) +
      0.5 * log(k0 / kn) -
      n / 2 * log(2 * pi),
    an = an
  )
}

# (b) y_i ~ N2(mu, S(rho)) with unit variances and correlation rho,
# mu_k ~ N(0, s0^2) independently and rho ~ U(-1, 1), on
# theta = (mu1, mu2, z) with rho = tanh(z). The Jacobian of rho = tanh(z) is
# 1 - rho^2, so the prior density of z is (1 - rho^2) / 2.
#
# With ybar the sample mean, the likelihood factors as
#   (2 pi)^-n |S|^(-n/2) exp(-Q(rho) / 2)
#     exp(-n (ybar - mu)' S^-1 (ybar - mu) / 2),
# where Q(rho) is the centred quadratic form, so integrating mu out against
# its prior gives
#   p(y | z) = (2 pi)^-n |S|^(-n/2) exp(-Q / 2) (2 pi) |S / n|^(1/2)
#              N2(ybar | 0, S / n + s0^2 I),
# and mu | z, y ~ N(m(z), P(z)^-1) with P = n S^-1 + I / s0^2 and
# m = P^-1 n S^-1 ybar.
bvn_model <- function(Y, s0) {
  n <- nrow(Y)
  ybar <- colMeans(Y)
  C <- crossprod(sweep(Y, 2, ybar)) # centred sums of squares and products

  # log p(y, z) with the means integrated out, vectorised over z
  log_joint_z <- function(z) {
    vapply(
      z,
      function(zz) {
        rho <- tanh(zz)
        S <- matrix(c(1, rho, rho, 1), 2)
        q <- (C[1, 1] - 2 * rho * C[1, 2] + C[2, 2]) * cosh(zz)^2
        loglik_centred <- -n * log(2 * pi) - n / 2 * log_sech2(zz) - q / 2
        # log of the integral over mu, with log det(S / n) = log(1 - rho^2)
        # - 2 log(n)
        log_det_s_n <- log_sech2(zz) - 2 * log(n)
        V <- S / n + s0^2 * diag(2)
        log_int_mu <- log(2 * pi) + 0.5 * log_det_s_n + mvn_logdens(ybar, V)
        loglik_centred + log_int_mu + log(0.5) + log_sech2(zz)
      },
      numeric(1)
    )
  }

  # Mean and variances of mu given z, as a 4 x length(z) matrix
  mu_given_z <- function(z) {
    vapply(
      z,
      function(zz) {
        rho <- tanh(zz)
        S_inv <- solve(matrix(c(1, rho, rho, 1), 2))
        P_inv <- solve(n * S_inv + diag(2) / s0^2)
        c(P_inv %*% (n * S_inv %*% ybar), diag(P_inv))
      },
      numeric(4)
    )
  }

  list(
    logpost = function(theta) {
      mu <- c(theta[["mu1"]], theta[["mu2"]])
      z <- theta[["z"]]
      rho <- tanh(z)
      e1 <- Y[, 1] - mu[1]
      e2 <- Y[, 2] - mu[2]
      q <- sum(e1^2 - 2 * rho * e1 * e2 + e2^2) * cosh(z)^2 # 1 / (1 - rho^2)
      loglik <- -n * log(2 * pi) - n / 2 * log_sech2(z) - q / 2
      loglik + sum(dnorm(mu, 0, s0, log = TRUE)) + log(0.5) + log_sech2(z)
    },
    log_joint_z = log_joint_z,
    mu_given_z = mu_given_z,
    ybar = ybar
  )
}

# (c) y_i ~ Weibull(shape k, scale lambda), log k ~ N(0, 1),
# log lambda ~ N(0, 2^2), on theta = (log_shape, log_scale). The priors are
# on theta itself, so there is no Jacobian. `log_joint` is the same
# log-posterior written out by hand and vectorised over (a, b) for the grid.
weibull_model <- function(y, prior_sd = c(1, 2)) {
  n <- length(y)
  ly <- log(y)
  log_joint <- function(a, b) {
    k <- exp(a)
    s <- rowSums(exp(outer(k, ly) - k * b)) # sum_i (y_i / lambda)^k
    loglik <- n * a - n * k * b + (k - 1) * sum(ly) - s
    loglik +
      dnorm(a, 0, prior_sd[1], log = TRUE) +
      dnorm(b, 0, prior_sd[2], log = TRUE)
  }
  list(
    logpost = function(theta) {
      a <- theta[["log_shape"]]
      b <- theta[["log_scale"]]
      sum(dweibull(y, shape = exp(a), scale = exp(b), log = TRUE)) +
        dnorm(a, 0, prior_sd[1], log = TRUE) +
        dnorm(b, 0, prior_sd[2], log = TRUE)
    },
    grad = function(theta) {
      a <- theta[["log_shape"]]
      b <- theta[["log_scale"]]
      k <- exp(a)
      e <- exp(k * (ly - b))
      c(
        n - n * k * b + k * sum(ly) - k * sum(e * (ly - b)) - a / prior_sd[1]^2,
        -n * k + k * sum(e) - b / prior_sd[2]^2
      )
    },
    # With u_i = log y_i - b and e_i = exp(k u_i), the gradient is
    # (n + k sum(u) - k sum(e u) - a / sd_a^2, -n k + k sum(e) - b / sd_b^2)
    hess = function(theta) {
      a <- theta[["log_shape"]]
      b <- theta[["log_scale"]]
      k <- exp(a)
      u <- ly - b
      e <- exp(k * u)
      h_ab <- -n * k + k * sum(e) + k^2 * sum(e * u)
      matrix(
        c(
          k * sum(u) - k * sum(e * u) - k^2 * sum(e * u^2) - 1 / prior_sd[1]^2,
          h_ab,
          h_ab,
          -k^2 * sum(e) - 1 / prior_sd[2]^2
        ),
        2
      )
    },
    log_joint = log_joint
  )
}

# Reference marginals of the Weibull posterior on a 2D grid of
# `n_grid` x `n_grid` points over +/- `half_width` SDs around the mode. The
# mode and SDs come from optim(), not from the code under test. The marginal
# densities are trapezoid sums over the other parameter; the normalising
# constant gives the log marginal likelihood.
weibull_reference <- function(model, half_width = 12, n_grid = 601) {
  fn <- function(x) -model$log_joint(x[1], x[2])
  opt <- optim(c(0, 0), fn, method = "BFGS", control = list(reltol = 1e-12))
  sds <- sqrt(diag(solve(optimHess(opt$par, fn))))
  u <- seq(-half_width, half_width, length.out = n_grid)
  a <- opt$par[1] + u * sds[1]
  b <- opt$par[2] + u * sds[2]
  wa <- trapezoid_weights(n_grid, a[2] - a[1])
  wb <- trapezoid_weights(n_grid, b[2] - b[1])
  ab <- expand.grid(a = a, b = b)
  L <- matrix(model$log_joint(ab$a, ab$b), n_grid, n_grid)
  L_max <- max(L)
  P <- exp(L - L_max)
  list(
    log_shape = grid_marginal(a, log(as.numeric(P %*% wb))),
    log_scale = grid_marginal(b, log(as.numeric(wa %*% P))),
    logml = L_max + log(sum(P * outer(wa, wb))),
    mode = opt$par,
    sds = sds
  )
}

# Data -------------------------------------------------------------------------

set.seed(3)
nig <- nig_model(rnorm(10, 1, 2), m0 = 0, k0 = 1, a0 = 2, b0 = 2)

set.seed(11)
bvn_y <- local({
  n <- 15
  rho <- 0.6
  Z <- matrix(rnorm(2 * n), n)
  cbind(Z[, 1] + 1, rho * Z[, 1] + sqrt(1 - rho^2) * Z[, 2] - 1)
})
bvn <- bvn_model(bvn_y, s0 = 5)

set.seed(5)
weib <- weibull_model(rweibull(15, shape = 1.5, scale = 2))

# Reference marginals of z and of the means for (b), on a z grid wide enough
# to hold all of the mass (the density at its ends is below 1e-30 of the
# peak).
bvn_z_grid <- seq(-1.5, 3.5, length.out = 4001)
bvn_log_joint <- bvn$log_joint_z(bvn_z_grid)
bvn_ref_z <- grid_marginal(bvn_z_grid, bvn_log_joint)

# (a) Normal-inverse-gamma ----------------------------------------------------

test_that("the normal-inverse-gamma references agree with quadrature", {
  lp_grid <- function(mu, s) {
    vapply(
      mu,
      function(m) nig$logpost(c(mu = m, log_sigma = s)),
      numeric(1)
    )
  }
  s_lims <- nig$log_sigma$mean + c(-15, 15) * nig$log_sigma$sd

  # The log marginal likelihood as the log of the 2D integral of exp(logpost)
  lp_max <- nig$logpost(nig$joint_mode)
  inner <- function(s) {
    vapply(
      s,
      function(si) {
        integrate(
          function(mu) exp(lp_grid(mu, si) - lp_max),
          -Inf,
          Inf,
          rel.tol = 1e-10
        )$value
      },
      numeric(1)
    )
  }
  int <- integrate(inner, s_lims[1], s_lims[2], rel.tol = 1e-10)$value
  expect_equal(lp_max + log(int), nig$logml, tolerance = 1e-8)

  # The joint mode and the Hessian there
  expect_equal(
    nig$grad(nig$joint_mode),
    c(0, 0),
    tolerance = 1e-12
  )
  H <- -num_jacobian(nig$grad, nig$joint_mode, h = 1e-6)
  expect_equal(solve(H), nig$joint_cov, tolerance = 1e-7)
  # The analytic gradient matches finite differences of logpost away from
  # the mode
  th <- c(mu = 0.3, log_sigma = 0.9)
  expect_equal(
    nig$grad(th),
    num_grad(nig$logpost, th, h = 1e-6),
    tolerance = 1e-7
  )

  # Integrating the joint density over the other parameter gives the closed
  # form marginal densities
  for (mu in nig$mu$mean + c(-3, -1, 0, 2) * nig$mu$sd) {
    f <- integrate(
      function(s) {
        vapply(
          s,
          function(si) exp(lp_grid(mu, si) - nig$logml),
          numeric(1)
        )
      },
      s_lims[1],
      s_lims[2],
      rel.tol = 1e-10
    )$value
    expect_equal(f, nig$mu$density(mu), tolerance = 1e-7)
  }
  for (s in nig$log_sigma$mean + c(-2, 0, 1, 3) * nig$log_sigma$sd) {
    f <- integrate(
      function(mu) exp(lp_grid(mu, s) - nig$logml),
      -Inf,
      Inf,
      rel.tol = 1e-10
    )$value
    expect_equal(f, nig$log_sigma$density(s), tolerance = 1e-7)
  }

  # The expected gradient at the VB mean, by nested quadrature over the
  # Laplace Gaussian, is zero
  sds <- sqrt(diag(nig$joint_cov))
  expected_grad <- function(k) {
    inner <- function(z2) {
      vapply(
        z2,
        function(zz2) {
          integrate(
            function(z1) {
              vapply(
                z1,
                function(zz1) {
                  nig$grad(nig$vb_mean + sds * c(zz1, zz2))[k] * dnorm(zz1)
                },
                numeric(1)
              )
            },
            -10,
            10,
            rel.tol = 1e-10
          )$value *
            dnorm(zz2)
        },
        numeric(1)
      )
    }
    integrate(inner, -10, 10, rel.tol = 1e-10)$value
  }
  expect_lt(abs(expected_grad(1)) * sds[1], 1e-8)
  expect_lt(abs(expected_grad(2)) * sds[2], 1e-8)

  # Moments, mode and quantiles of the log_sigma density
  dens <- nig$log_sigma$density
  int_s <- function(fn) {
    integrate(fn, s_lims[1], s_lims[2], rel.tol = 1e-10)$value
  }
  expect_equal(int_s(dens), 1, tolerance = 1e-8)
  m1 <- int_s(function(s) s * dens(s))
  expect_equal(m1, nig$log_sigma$mean, tolerance = 1e-8)
  expect_equal(
    sqrt(int_s(function(s) (s - m1)^2 * dens(s))),
    nig$log_sigma$sd,
    tolerance = 1e-8
  )
  mode_q <- optimize(
    function(s) log(dens(s)),
    nig$log_sigma$mode + c(-1, 1) * nig$log_sigma$sd,
    maximum = TRUE,
    tol = 1e-10
  )$maximum
  expect_equal(mode_q, nig$log_sigma$mode, tolerance = 1e-6)
  for (p in probs) {
    cdf <- integrate(
      dens,
      s_lims[1],
      nig$log_sigma$quantile(p),
      rel.tol = 1e-10
    )$value
    expect_equal(cdf, p, tolerance = 1e-8)
  }
})

test_that("normal data with unknown mean and variance give the exact Laplace step", {
  grads <- list(finite_differences = NULL, analytic = nig$grad)
  for (gr in grads) {
    fit <- inla_vanilla(
      nig$logpost,
      init = c(mu = 0, log_sigma = 0),
      grad = gr,
      silent = TRUE
    )
    sds <- sqrt(diag(nig$joint_cov))
    # Mode: observed error 6e-8 SDs (the optimiser's tolerance)
    expect_lt(max(abs(fit$mode - nig$joint_mode) / sds), 1e-6)
    # Covariance, relative to the SDs: observed error 3.3e-6 with finite
    # differences, 1.4e-8 with the gradient
    cov_err <- max(abs(fit$cov - nig$joint_cov) / outer(sds, sds))
    expect_lt(cov_err, if (is.null(gr)) 1e-5 else 5e-8)
    # The Laplace log marginal likelihood from the exact mode and covariance
    logml_laplace <- nig$logpost(nig$joint_mode) +
      log(2 * pi) +
      0.5 * log(det(nig$joint_cov))
    # Observed error 6e-7 with finite differences, 1e-9 with the gradient
    expect_lt(abs(fit$logml_laplace - logml_laplace), 2e-6)
  }
})

test_that("normal data with unknown mean and variance give accurate marginals", {
  grads <- list(finite_differences = NULL, analytic = nig$grad)
  for (gr in grads) {
    fit <- inla_vanilla(
      nig$logpost,
      init = c(mu = 0, log_sigma = 0),
      grad = gr,
      silent = TRUE
    )

    # log_sigma: a log inverse-gamma, skewed to the right (skewness 0.39)
    err <- marginal_errors(fit, "log_sigma", nig$log_sigma)
    # Quantiles: observed worst SN error 0.0067 SDs, Laplace 0.68 SDs
    expect_lt(max(abs(err$quantile_sn)), 0.02)
    expect_lt(max(abs(err$quantile_sn)), max(abs(err$quantile_laplace)))
    expect_lt(mean(abs(err$quantile_sn)), mean(abs(err$quantile_laplace)))
    # Mean: observed SN error 0.0026 SDs, Laplace (the mode) 0.36 SDs
    expect_lt(abs(err$mean_sn), 0.008)
    expect_lt(abs(err$mean_sn), abs(err$mean_laplace))
    # SD: observed relative error 0.0022 for SN, 0.068 for Laplace
    expect_lt(abs(err$sd_sn), 0.007)
    expect_lt(abs(err$sd_sn), abs(err$sd_laplace))
    # Mode of the skew-normal: observed error 0.017 SDs
    mode_err <- (fit$summary["log_sigma", "Mode"] - nig$log_sigma$mode) /
      nig$log_sigma$sd
    expect_lt(abs(mode_err), 0.05)
    # Right skew, and the VB mean moves right: observed shift 0.36 SDs, VB
    # mean error 0.028 SDs against the posterior mean
    expect_gt(fit$sn["log_sigma", "alpha"], 1)
    expect_gt(fit$mean_vb[["log_sigma"]], fit$mode[["log_sigma"]])
    expect_lt(abs(err$mean_vb), 0.08)
    # Against the exact VB solution, the VB mean is off only by the error of
    # the cubature rule: observed 5.2e-5 SDs for log_sigma, 3e-5 (finite
    # differences) or 5e-9 (gradient) for mu. The VB term of the log marginal
    # likelihood, 0.5 shift' H shift, is off by 1.9e-5 against the exact
    # shift and the exact H.
    sds <- sqrt(diag(nig$joint_cov))
    expect_lt(max(abs(fit$mean_vb - nig$vb_mean) / sds), 2e-4)
    shift <- nig$vb_mean - nig$joint_mode
    expect_lt(
      abs(fit$logml - fit$logml_laplace - 0.5 * sum(shift^2 / sds^2)),
      6e-5
    )

    # mu: a Student-t with 2 an = 14 degrees of freedom, symmetric about mn.
    # The covariance at the joint mode is diagonal, so the scan for mu holds
    # log_sigma at its mode, the scanned log-heights are exactly quadratic in
    # mu, and the exact volume correction at the mode is zero by symmetry.
    # The skew-normal is then the Laplace Gaussian
    # N(mn, 2 bn / (kn (2 an + 1))), whose SD is
    # sqrt(2 (an - 1) / (2 an + 1)) = 0.894 times the Student-t SD. The
    # linear scan cannot see the extra spread from the uncertainty in sigma,
    # which is an O(1 / n) limitation of the method, not of this code.
    err <- marginal_errors(fit, "mu", nig$mu)
    ratio <- sqrt(2 * (nig$an - 1) / (2 * nig$an + 1))
    expect_equal(err$sd_sn + 1, ratio, tolerance = 1e-4)
    # The whole skew-normal is that exact Gaussian. Its shape alpha is not
    # checked: near alpha = 0 a skew-normal differs from a shifted Gaussian
    # only at third order, so alpha is weakly identified (it came out between
    # -0.026 and 9e-6 on equivalent runs, a skewness below 1e-5). Observed:
    # SD equal to 1e-8 relative; quantiles shifted by -3.3e-4 Laplace SDs
    # with or without the gradient, which is the volume correction gamma1:
    # its forward difference of 0.01 SDs picks up the curvature of a
    # log-determinant that is flat at the mode, so it is O(0.01), not 0.
    lap_sd <- sqrt(nig$joint_cov[1, 1])
    tab <- summary(fit, probs = probs)$table
    expect_lt(abs(tab["mu", "SD"] / lap_sd - 1), 1e-6)
    q_lap <- qnorm(probs, nig$joint_mode[["mu"]], lap_sd)
    expect_lt(max(abs(tab["mu", q_cols] - q_lap)) / lap_sd, 1e-3)
    # Mean and median against the t: the same -3e-4 SDs
    expect_lt(abs(err$mean_sn), 1e-3)
    expect_lt(abs(err$quantile_sn[[4]]), 1e-3)
    # Quantiles: observed worst SN error 0.233 SDs at the 2.5% and 97.5%
    # points, the same as Laplace, and 0.038 SDs at the quartiles. These
    # errors are fixed by the SD ratio checked above (they are the
    # narrower Gaussian against the heavier-tailed t), not noise.
    expect_lt(max(abs(err$quantile_sn)), 0.35)
    expect_lt(max(abs(err$quantile_sn[3:5])), 0.08)

    # Log marginal likelihood: observed error 0.0011 with the VB term,
    # -0.063 for plain Laplace
    expect_lt(abs(fit$logml - nig$logml), 0.003)
    expect_lt(
      abs(fit$logml - nig$logml),
      abs(fit$logml_laplace - nig$logml)
    )
  }
})

# (b) Bivariate normal with an unknown correlation ----------------------------

test_that("the bivariate normal references agree with quadrature", {
  # The logpost is the bivariate normal likelihood plus the priors and the
  # Jacobian, written another way
  th <- c(mu1 = 0.4, mu2 = -1.2, z = 0.9)
  S <- matrix(c(1, tanh(0.9), tanh(0.9), 1), 2)
  expect_equal(
    bvn$logpost(th),
    sum(apply(bvn_y, 1, function(y) mvn_logdens(y - th[1:2], S))) +
      sum(dnorm(th[1:2], 0, 5, log = TRUE)) +
      log((1 - tanh(0.9)^2) / 2),
    tolerance = 1e-12
  )

  # The edges of the z grid carry no mass
  expect_lt(bvn_ref_z$edge, 1e-20)

  # Integrating the means out analytically agrees with 2D quadrature
  for (z in c(0.3, 1.4)) {
    lp_max <- bvn$logpost(c(mu1 = bvn$ybar[[1]], mu2 = bvn$ybar[[2]], z = z))
    inner <- function(mu1) {
      vapply(
        mu1,
        function(m1) {
          integrate(
            function(mu2) {
              vapply(
                mu2,
                function(m2) {
                  exp(bvn$logpost(c(mu1 = m1, mu2 = m2, z = z)) - lp_max)
                },
                numeric(1)
              )
            },
            bvn$ybar[[2]] - 4,
            bvn$ybar[[2]] + 4,
            rel.tol = 1e-10
          )$value
        },
        numeric(1)
      )
    }
    int <- integrate(
      inner,
      bvn$ybar[[1]] - 4,
      bvn$ybar[[1]] + 4,
      rel.tol = 1e-10
    )$value
    expect_equal(lp_max + log(int), bvn$log_joint_z(z), tolerance = 1e-8)
  }

  # The grid marginal of z against adaptive quadrature
  log_max <- max(bvn_log_joint)
  dens <- function(z) exp(bvn$log_joint_z(z) - log_max)
  norm_const <- integrate(dens, -1.5, 3.5, rel.tol = 1e-10)$value
  for (p in c(0.025, 0.5, 0.975)) {
    cdf <- integrate(
      dens,
      -1.5,
      bvn_ref_z$quantile(p),
      rel.tol = 1e-10
    )$value
    # Observed CDF errors up to 9e-7
    expect_lt(abs(cdf / norm_const - p), 5e-6)
  }
  m1 <- integrate(function(z) z * dens(z), -1.5, 3.5, rel.tol = 1e-10)$value
  expect_equal(m1 / norm_const, bvn_ref_z$mean, tolerance = 1e-8)
  expect_equal(log_max + log(norm_const), bvn_ref_z$log_norm, tolerance = 1e-10)
})

test_that("the marginal of a correlation is recovered on the Fisher z scale", {
  fit <- inla_vanilla(
    bvn$logpost,
    init = c(mu1 = 0, mu2 = 0, z = 0),
    silent = TRUE
  )

  # The Laplace step against a Hessian by nested central differences written
  # here, which agree to 1e-7 across steps from 2e-5 to 5e-4. Observed:
  # whitened gradient at the mode 3e-7, covariance error 1.4e-5 relative to
  # the SDs (the package's second differences of logpost), logml_laplace
  # error 9e-7.
  lp <- function(x) bvn$logpost(stats::setNames(x, c("mu1", "mu2", "z")))
  lap <- laplace_errors(
    fit,
    lp,
    grad_ref = function(x) num_grad(lp, x),
    hess_ref = function(x) {
      num_jacobian(function(y) num_grad(lp, y, 1e-4), x, 1e-4)
    }
  )
  expect_lt(lap$grad, 1e-6)
  expect_lt(lap$cov, 5e-5)
  expect_lt(abs(lap$logml_laplace), 3e-6)

  # z is strongly skewed to the left (rho near 0.8 with n = 15)
  err <- marginal_errors(fit, "z", bvn_ref_z)
  # Quantiles: observed worst SN error 0.047 SDs (at 10%), Laplace 1.02 SDs
  # (at 2.5%); averages 0.017 and 0.52
  expect_lt(max(abs(err$quantile_sn)), 0.1)
  expect_lt(max(abs(err$quantile_sn)), max(abs(err$quantile_laplace)))
  expect_lt(mean(abs(err$quantile_sn)), mean(abs(err$quantile_laplace)))
  # Mean: observed SN error 0.012 SDs, Laplace 0.44 SDs
  expect_lt(abs(err$mean_sn), 0.03)
  expect_lt(abs(err$mean_sn), abs(err$mean_laplace))
  # SD: observed relative error 2.4e-4 for SN, 0.14 for Laplace
  expect_lt(abs(err$sd_sn), 1e-3)
  expect_lt(abs(err$sd_sn), abs(err$sd_laplace))
  # Left skew, and the VB mean moves left: observed shift -0.44 SDs, VB mean
  # error 0.066 SDs
  expect_lt(fit$sn["z", "alpha"], -1)
  expect_lt(fit$mean_vb[["z"]], fit$mode[["z"]])
  expect_lt(abs(err$mean_vb), 0.15)

  # The means: mixtures over z of the Gaussians mu | z, y, close to Gaussian
  mz <- bvn$mu_given_z(bvn_z_grid)
  w <- trapezoid_weights(length(bvn_z_grid), diff(bvn_z_grid[1:2])) *
    bvn_ref_z$density
  for (k in 1:2) {
    par <- paste0("mu", k)
    mu_mean <- sum(w * mz[k, ])
    mu_sd <- sqrt(sum(w * (mz[k + 2, ] + mz[k, ]^2)) - mu_mean^2)
    mix_cdf <- function(x) sum(w * pnorm(x, mz[k, ], sqrt(mz[k + 2, ])))
    q_exact <- vapply(
      probs,
      function(p) {
        uniroot(
          function(x) mix_cdf(x) - p,
          mu_mean + c(-6, 6) * mu_sd,
          tol = 1e-10
        )$root
      },
      numeric(1)
    )
    ref <- list(
      mean = mu_mean,
      sd = mu_sd,
      quantile = function(p) q_exact[match(p, probs)]
    )
    err <- marginal_errors(fit, par, ref)
    # Observed errors: mean 1.6e-4 SDs, SD 1e-4 relative, quantiles 3.7e-4
    # SDs
    expect_lt(abs(err$mean_sn), 6e-4)
    expect_lt(abs(err$sd_sn), 4e-4)
    expect_lt(max(abs(err$quantile_sn)), 1e-3)
  }

  # Log marginal likelihood, the log integral of p(y, z): observed error
  # -0.0013 with the VB term, -0.098 for plain Laplace
  logml <- bvn_ref_z$log_norm
  expect_lt(abs(fit$logml - logml), 0.004)
  expect_lt(abs(fit$logml - logml), abs(fit$logml_laplace - logml))
})

# (c) Weibull, no closed form -------------------------------------------------

test_that("the Weibull grid reference is converged and matches the logpost", {
  # The vectorised log-posterior of the grid is the logpost given to the fit
  a <- c(-0.5, 0.2, 0.6, 1.4)
  b <- c(0.1, 1.2, 0.7, -0.3)
  expect_equal(
    weib$log_joint(a, b),
    vapply(
      seq_along(a),
      function(i) weib$logpost(c(log_shape = a[i], log_scale = b[i])),
      numeric(1)
    ),
    tolerance = 1e-12
  )
  th <- c(log_shape = 0.3, log_scale = 0.4)
  expect_equal(
    weib$grad(th),
    num_grad(weib$logpost, th, h = 1e-6),
    tolerance = 1e-7
  )
  # The analytic Hessian matches differences of the analytic gradient
  expect_equal(
    weib$hess(th),
    num_jacobian(weib$grad, th, h = 1e-6),
    tolerance = 1e-7
  )

  ref <- weibull_reference(weib)
  # Nothing is lost off the grid: the marginal densities at its edges are
  # 2e-13 (log_shape) and 5e-8 (log_scale, which has heavy tails) of their
  # peaks
  expect_lt(ref$log_shape$edge, 1e-6)
  expect_lt(ref$log_scale$edge, 1e-6)
  # A coarser, narrower grid gives the same answers: observed differences
  # 7e-4 SDs in the quantiles and 8e-7 in the log marginal likelihood
  coarse <- weibull_reference(weib, half_width = 9, n_grid = 301)
  for (par in c("log_shape", "log_scale")) {
    expect_lt(
      max(abs(coarse[[par]]$quantile(probs) - ref[[par]]$quantile(probs))) /
        ref[[par]]$sd,
      3e-3
    )
  }
  expect_lt(abs(coarse$logml - ref$logml), 5e-6)
  # The grid's log marginal likelihood against nested adaptive quadrature
  lims_a <- ref$mode[1] + c(-12, 12) * ref$sds[1]
  lims_b <- ref$mode[2] + c(-12, 12) * ref$sds[2]
  L_max <- weib$log_joint(ref$mode[1], ref$mode[2])
  inner <- function(a) {
    vapply(
      a,
      function(ai) {
        integrate(
          function(b) exp(weib$log_joint(rep(ai, length(b)), b) - L_max),
          lims_b[1],
          lims_b[2],
          rel.tol = 1e-10
        )$value
      },
      numeric(1)
    )
  }
  int <- integrate(inner, lims_a[1], lims_a[2], rel.tol = 1e-10)$value
  # Observed difference 4e-12
  expect_lt(abs(L_max + log(int) - ref$logml), 1e-8)
})

test_that("a non-conjugate Weibull posterior is recovered against 2D quadrature", {
  ref <- weibull_reference(weib)
  grads <- list(finite_differences = NULL, analytic = weib$grad)
  for (gr in grads) {
    fit <- inla_vanilla(
      weib$logpost,
      init = c(log_shape = 0, log_scale = 0),
      grad = gr,
      silent = TRUE
    )
    err_shape <- marginal_errors(fit, "log_shape", ref$log_shape)
    err_scale <- marginal_errors(fit, "log_scale", ref$log_scale)

    # The Laplace step against the analytic gradient and Hessian. Observed:
    # whitened gradient at the mode 1e-7; covariance error 1.8e-6 relative to
    # the SDs with finite differences, 2.5e-10 with the gradient;
    # logml_laplace error 1.2e-6 and 1e-10.
    lap <- laplace_errors(fit, weib$logpost, weib$grad, weib$hess)
    expect_lt(lap$grad, 1e-6)
    expect_lt(lap$cov, if (is.null(gr)) 1e-5 else 1e-8)
    expect_lt(abs(lap$logml_laplace), if (is.null(gr)) 5e-6 else 1e-8)

    # log_shape is skewed to the left (skewness -0.27). Quantiles: observed
    # worst SN error 0.0044 SDs, Laplace 0.56 SDs
    expect_lt(fit$sn["log_shape", "alpha"], -1)
    expect_lt(max(abs(err_shape$quantile_sn)), 0.015)
    expect_lt(
      max(abs(err_shape$quantile_sn)),
      max(abs(err_shape$quantile_laplace))
    )
    # Mean: observed SN error 1.5e-4 SDs, Laplace 0.34 SDs
    expect_lt(abs(err_shape$mean_sn), 5e-4)
    expect_lt(abs(err_shape$mean_sn), abs(err_shape$mean_laplace))
    # SD: observed relative error 1.2e-4 for SN, 0.047 for Laplace
    expect_lt(abs(err_shape$sd_sn), 4e-4)
    expect_lt(abs(err_shape$sd_sn), abs(err_shape$sd_laplace))

    # log_scale is nearly symmetric (skewness -0.06) but heavy-tailed (excess
    # kurtosis 0.72): given the shape k its conditional SD is about
    # 1 / (k sqrt(n)), so the uncertainty in k spreads it out. As for mu in
    # (a), the linear scan with its linear volume correction cannot see this,
    # and both the skew-normal and the Laplace SDs are about 11% too small.
    # While calibrating, a full Laplace marginal (exact conditional modes and
    # curvatures) had an SD error of 0.1%, so the gap is the linearisation in
    # the method, not an error in this code. Quantiles: observed worst SN
    # error 0.28 SDs at 97.5% (Laplace 0.26), and 0.044 SDs at the quartiles
    # (Laplace 0.078)
    expect_lt(max(abs(err_scale$quantile_sn)), 0.5)
    expect_lt(max(abs(err_scale$quantile_sn[3:5])), 0.1)
    # Mean: observed SN error 0.010 SDs, Laplace 0.032 SDs
    expect_lt(abs(err_scale$mean_sn), 0.025)
    expect_lt(abs(err_scale$mean_sn), abs(err_scale$mean_laplace))
    # SD: observed relative error -0.110 for SN, -0.118 for Laplace
    expect_lt(abs(err_scale$sd_sn), 0.13)
    expect_lt(abs(err_scale$sd_sn), abs(err_scale$sd_laplace))

    # The VB means, against the reference means: observed errors 0.020 SDs
    # (log_shape, moved by -0.34 SDs from the mode) and 0.004 SDs (log_scale)
    expect_lt(abs(err_shape$mean_vb), 0.05)
    expect_lt(abs(err_scale$mean_vb), 0.05)

    # Averaged over the quantiles, the skew-normal beats Laplace for each
    # parameter (observed 0.003 against 0.37 SDs for log_shape, 0.110 against
    # 0.124 for log_scale; the same ordering on seeds 1 to 8) and by a
    # factor of 4 overall
    avg <- function(e) mean(abs(e))
    expect_lt(avg(err_shape$quantile_sn), avg(err_shape$quantile_laplace))
    expect_lt(avg(err_scale$quantile_sn), avg(err_scale$quantile_laplace))
    expect_lt(
      avg(c(err_shape$quantile_sn, err_scale$quantile_sn)),
      0.5 * avg(c(err_shape$quantile_laplace, err_scale$quantile_laplace))
    )

    # Log marginal likelihood: observed error 0.0030 with the VB term, -0.063
    # for plain Laplace
    expect_lt(abs(fit$logml - ref$logml), 0.01)
    expect_lt(abs(fit$logml - ref$logml), abs(fit$logml_laplace - ref$logml))
  }
})

# (d) Invariance to the order of the parameters --------------------------------

# Compare two fits of the same posterior with the parameters in different
# orders, matching parameters by name. Location differences are in posterior
# SDs. The Laplace and VB stages (mode, covariance, VB mean, logml) are held
# to `tol_laplace`, and everything that passes through a skew-normal fit to
# `tol_sn`, since nlminb() stops the fit at a tolerance of its own. The
# skew-normal parameters themselves are not compared: for a near-Gaussian
# marginal (the means in (b)) xi and alpha are weakly identified, and they
# differ by 1e-4 between the orders while the fitted distributions agree to
# 1e-8 SDs.
expect_same_marginals <- function(fit1, fit2, tol_laplace, tol_sn) {
  pars <- names(fit1$mode)
  expect_setequal(names(fit2$mode), pars)
  sds <- sqrt(diag(fit1$cov))
  expect_lt(max(abs(fit2$mode[pars] - fit1$mode) / sds), tol_laplace)
  expect_lt(max(abs(fit2$mean_vb[pars] - fit1$mean_vb) / sds), tol_laplace)
  expect_lt(
    max(abs(fit2$cov[pars, pars] - fit1$cov) / outer(sds, sds)),
    tol_laplace
  )
  expect_lt(abs(fit2$logml - fit1$logml), tol_laplace)
  tab1 <- summary(fit1, probs = probs)$table
  tab2 <- summary(fit2, probs = probs)$table[pars, ]
  expect_lt(max(abs(tab2 - tab1) / sds), tol_sn)
  d1 <- summary(fit1)$diagnostics
  d2 <- summary(fit2)$diagnostics[pars, ]
  expect_lt(max(abs(d2$vb_shift_sigma - d1$vb_shift_sigma)), tol_laplace)
  expect_lt(max(abs(d2$sn_vs_vb - d1$sn_vs_vb)), tol_sn)
}

test_that("permuting the parameters permutes the marginals (finite differences)", {
  fit1 <- inla_vanilla(
    bvn$logpost,
    init = c(mu1 = 0, mu2 = 0, z = 0),
    silent = TRUE
  )
  # The same posterior with the parameters in the order (mu2, z, mu1), a
  # cycle rather than a reversal, so that the permutation is not its own
  # inverse
  logpost_perm <- function(theta) {
    bvn$logpost(c(mu1 = theta[[3]], mu2 = theta[[1]], z = theta[[2]]))
  }
  fit2 <- inla_vanilla(
    logpost_perm,
    init = c(mu2 = 0, z = 0, mu1 = 0),
    silent = TRUE
  )
  expect_equal(rownames(fit2$summary), c("mu2", "z", "mu1"))
  # Observed differences: 3.4e-8 SDs in the summaries, up to 5e-10 in the
  # mode, VB mean and covariance, and 3e-10 in the log marginal likelihood
  expect_same_marginals(fit1, fit2, tol_laplace = 1e-8, tol_sn = 1e-7)
})

test_that("permuting the parameters permutes the marginals (analytic gradient)", {
  fit1 <- inla_vanilla(
    weib$logpost,
    init = c(log_shape = 0, log_scale = 0),
    grad = weib$grad,
    silent = TRUE
  )
  logpost_perm <- function(theta) {
    weib$logpost(c(log_shape = theta[[2]], log_scale = theta[[1]]))
  }
  grad_perm <- function(theta) {
    rev(weib$grad(c(log_shape = theta[[2]], log_scale = theta[[1]])))
  }
  fit2 <- inla_vanilla(
    logpost_perm,
    init = c(log_scale = 0, log_shape = 0),
    grad = grad_perm,
    silent = TRUE
  )
  expect_equal(rownames(fit2$summary), c("log_scale", "log_shape"))
  # Observed differences: 4e-9 SDs in the summaries, 1e-15 or less in the
  # Laplace and VB stages
  expect_same_marginals(fit1, fit2, tol_laplace = 1e-8, tol_sn = 1e-7)
})
