# End-to-end tests of inla_vanilla() on one-parameter conjugate models, whose
# posteriors and marginal likelihoods are known in closed form. Every logpost
# carries all of its normalising constants (and the Jacobian of the
# reparameterisation), so fit$logml is directly comparable with the exact log
# marginal likelihood. Location errors are measured in exact posterior SDs and
# SD errors as relative errors.
#
# Each stage is checked against its own exact target, so that a bug in one
# stage cannot hide inside the approximation error of another:
# * Laplace: the exact mode, the exact curvature there and the exact Laplace
#   log marginal likelihood (all in closed form), tightly.
# * VB: the exact VB mean, the root of E[d logpost / d theta] = 0 under
#   N(mu, 1 / H) (in closed form or by quadrature), and the exact VB term
#   0.5 H (mu - mode)^2. Only the error of the cubature rule is allowed.
# * Skew-normal: the summary against the fitted skew-normal (by quadrature),
#   and the fitted skew-normal against the exact posterior, with tolerances
#   that reflect the approximation error and an assertion that it beats the
#   Laplace Gaussian.
#
# The tolerances on the non-Gaussian cases were calibrated by printing the
# actual errors first; each comment states the observed value. The skew-normal
# results with finite differences and with an analytic gradient agree to about
# 1e-7 SDs, so the same approximation tolerances serve both.

probs <- c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975)

# Models ----------------------------------------------------------------------
#
# Each builder returns the logpost (a function of the named vector theta, as
# inla_vanilla() passes it), its gradient, the exact posterior of theta
# (density, mean, SD, mode and quantile function), the exact log marginal
# likelihood, the exact curvature H = -d2 logpost / d theta2 at the mode and
# the exact VB mean.

# The exact Laplace estimate of the log marginal likelihood, and the exact VB
# estimate (Laplace plus the location term 0.5 H (mu_vb - mode)^2).
add_laplace_and_vb <- function(model) {
  model$laplace_logml <- model$logpost(c(theta = model$mode)) +
    0.5 * log(2 * pi) -
    0.5 * log(model$laplace_H)
  model$vb_logml <- model$laplace_logml +
    0.5 * model$laplace_H * (model$vb_mean - model$mode)^2
  model
}

# y_i ~ N(mu, sigma^2) with sigma known, mu ~ N(m0, s0^2). The posterior of mu
# is Gaussian, and marginally y ~ N_n(m0 1, sigma^2 I + s0^2 1 1'). For a
# Gaussian the Laplace and VB answers are exact.
normal_model <- function(y, sigma, m0, s0) {
  n <- length(y)
  prec <- 1 / s0^2 + n / sigma^2
  post_mean <- (m0 / s0^2 + sum(y) / sigma^2) / prec
  post_sd <- sqrt(1 / prec)

  V <- sigma^2 * diag(n) + s0^2 * matrix(1, n, n)
  U <- chol(V)
  r <- backsolve(U, y - m0, transpose = TRUE)
  logml <- -n / 2 * log(2 * pi) - sum(log(diag(U))) - 0.5 * sum(r^2)

  add_laplace_and_vb(list(
    logpost = function(theta) {
      mu <- theta[[1]]
      sum(dnorm(y, mu, sigma, log = TRUE)) + dnorm(mu, m0, s0, log = TRUE)
    },
    grad = function(theta) {
      mu <- theta[[1]]
      sum(y - mu) / sigma^2 - (mu - m0) / s0^2
    },
    density = function(t) dnorm(t, post_mean, post_sd),
    mean = post_mean,
    sd = post_sd,
    mode = post_mean,
    quantile = function(p) qnorm(p, post_mean, post_sd),
    logml = logml,
    laplace_H = prec,
    vb_mean = post_mean
  ))
}

# y ~ Binomial(n, p), p ~ Beta(a, b), on theta = logit(p). The Jacobian is
# dp / dtheta = p (1 - p). Exactly p | y ~ Beta(A, B) with A = a + y and
# B = b + n - y, so theta | y is a logit-beta: E = digamma(A) - digamma(B),
# Var = trigamma(A) + trigamma(B), mode log(A / B), quantiles
# qlogis(qbeta(...)). The marginal likelihood is the beta-binomial pmf.
# d logpost / d theta = A - (A + B) p, so H = A B / (A + B) at the mode, and
# the VB mean solves E[plogis(mu + Z / sqrt(H))] = A / (A + B), Z ~ N(0, 1).
beta_binomial_model <- function(y, n, a, b) {
  A <- a + y
  B <- b + n - y
  H <- A * B / (A + B)
  expected_p <- function(mu) {
    integrate(
      function(z) plogis(mu + z / sqrt(H)) * dnorm(z),
      -12,
      12,
      rel.tol = 1e-13
    )$value
  }
  vb_mean <- uniroot(
    function(mu) expected_p(mu) - A / (A + B),
    log(A / B) + c(-1, 1),
    tol = 1e-13
  )$root
  add_laplace_and_vb(list(
    logpost = function(theta) {
      p <- plogis(theta[[1]])
      dbinom(y, n, p, log = TRUE) +
        dbeta(p, a, b, log = TRUE) +
        plogis(theta[[1]], log.p = TRUE) +
        plogis(-theta[[1]], log.p = TRUE)
    },
    grad = function(theta) A - (A + B) * plogis(theta[[1]]),
    density = function(t) dbeta(plogis(t), A, B) * plogis(t) * plogis(-t),
    mean = digamma(A) - digamma(B),
    sd = sqrt(trigamma(A) + trigamma(B)),
    mode = log(A / B),
    quantile = function(p) qlogis(qbeta(p, A, B)),
    logml = lchoose(n, y) + lbeta(A, B) - lbeta(a, b),
    laplace_H = H,
    vb_mean = vb_mean,
    natural_quantile = function(p) qbeta(p, A, B),
    natural_mean = A / (A + B),
    natural_sd = sqrt(A * B / ((A + B)^2 * (A + B + 1)))
  ))
}

# y_i ~ Poisson(lambda), lambda ~ Gamma(shape, rate), on theta = log(lambda).
# The Jacobian is dlambda / dtheta = lambda. Exactly lambda | y ~ Gamma(S, R)
# with S = shape + sum(y) and R = rate + n, so theta | y is a log-gamma:
# E = digamma(S) - log(R), Var = trigamma(S), mode log(S / R), quantiles
# log(qgamma(...)). The marginal likelihood:
#   p(y) = int prod_i [exp(-lambda) lambda^y_i / y_i!]
#            rate^shape lambda^(shape - 1) exp(-rate lambda) / Gamma(shape)
#            dlambda
#        = rate^shape / (Gamma(shape) prod_i y_i!)
#            int lambda^(S - 1) exp(-R lambda) dlambda
#        = rate^shape Gamma(S) / (Gamma(shape) prod_i y_i! R^S),
# a negative-binomial-type closed form, checked by quadrature below.
# d logpost / d theta = S - R exp(theta), so H = S at the mode. Under
# N(mu, 1 / S), E[exp(theta)] = exp(mu + 1 / (2 S)), so the VB mean is
# log(S / R) - 1 / (2 S).
gamma_poisson_model <- function(y, shape, rate) {
  n <- length(y)
  S <- shape + sum(y)
  R <- rate + n
  log_prior_const <- shape * log(rate) - lgamma(shape)
  logml <- log_prior_const - sum(lfactorial(y)) + lgamma(S) - S * log(R)
  add_laplace_and_vb(list(
    logpost = function(theta) {
      lambda <- exp(theta[[1]])
      sum(dpois(y, lambda, log = TRUE)) +
        dgamma(lambda, shape = shape, rate = rate, log = TRUE) +
        theta[[1]]
    },
    grad = function(theta) S - R * exp(theta[[1]]),
    density = function(t) dgamma(exp(t), shape = S, rate = R) * exp(t),
    mean = digamma(S) - log(R),
    sd = sqrt(trigamma(S)),
    mode = log(S / R),
    quantile = function(p) log(qgamma(p, shape = S, rate = R)),
    logml = logml,
    laplace_H = S,
    vb_mean = log(S / R) - 1 / (2 * S),
    natural_quantile = function(p) qgamma(p, shape = S, rate = R),
    natural_mean = S / R,
    natural_sd = sqrt(S) / R
  ))
}

set.seed(1)
normal <- normal_model(y = rnorm(12, 3, 1.5), sigma = 1.5, m0 = 1, s0 = 2)
# Skewness of theta about -0.46 (A = 3, B = 14)
bb_skewed <- beta_binomial_model(y = 2, n = 15, a = 1, b = 1)
# Skewness of theta about -0.03 (A = 20, B = 24)
bb_symmetric <- beta_binomial_model(y = 18, n = 40, a = 2, b = 2)
# Skewness of theta about -0.47 (S = 5, R = 6)
gp_skewed <- gamma_poisson_model(y = c(0, 1, 0, 2, 1), shape = 1, rate = 1)
# Skewness of theta about -0.09 (S = 120, R = 20.5)
gp_symmetric <- gamma_poisson_model(
  y = c(5, 7, 4, 6, 9, 3, 6, 8, 5, 7, 6, 4, 7, 5, 8, 6, 5, 7, 6, 4),
  shape = 2,
  rate = 0.5
)

# Helpers ---------------------------------------------------------------------

# logpost on a vector of scalar values, for quadrature
logpost_vec <- function(model, t) {
  vapply(t, function(ti) model$logpost(c(theta = ti)), numeric(1))
}

# Integrate fn over the exact mode +/- 15 exact SDs, where all the mass is
integrate_posterior <- function(model, fn) {
  lims <- model$mode + c(-15, 15) * model$sd
  integrate(fn, lims[1], lims[2], rel.tol = 1e-10)$value
}

# The skew-normal density, written out here so that the checks below do not
# rely on the package's own dsnorm()
dsn_ref <- function(x, xi, omega, alpha) {
  z <- (x - xi) / omega
  2 / omega * dnorm(z) * pnorm(alpha * z)
}

# Largest discrepancy, in exact posterior SDs, between the default summary
# table and the skew-normal it claims to summarise, all by quadrature: mean,
# SD, the CDF at each quantile (converted to a quantile error through the
# density) and the mode.
summary_vs_sn <- function(fit, model) {
  sn <- fit$sn[1, ]
  tab <- fit$summary[1, ]
  f <- function(x) dsn_ref(x, sn[["xi"]], sn[["omega"]], sn[["alpha"]])
  lims <- sn[["xi"]] + c(-20, 20) * sn[["omega"]]
  int <- function(fn, upper = lims[2]) {
    integrate(fn, lims[1], upper, rel.tol = 1e-12, subdivisions = 500L)$value
  }
  m1 <- int(function(x) x * f(x))
  s1 <- sqrt(int(function(x) (x - m1)^2 * f(x)))
  q_names <- c("2.5%", "50%", "97.5%")
  q_err <- vapply(
    seq_along(q_names),
    function(i) {
      q <- tab[[q_names[i]]]
      (int(f, q) - c(0.025, 0.5, 0.975)[i]) / f(q)
    },
    numeric(1)
  )
  mode <- optimize(
    function(x) log(f(x)),
    tab[["Mode"]] + c(-0.5, 0.5) * s1,
    maximum = TRUE,
    tol = 1e-12
  )$maximum
  max(abs(c(
    tab[["Mean"]] - m1,
    tab[["SD"]] - s1,
    q_err,
    tab[["Mode"]] - mode
  ))) /
    model$sd
}

# Errors of a fit against the exact posterior, for the skew-normal marginal
# and for the Laplace Gaussian N(mode, cov): mean, mode and the largest
# quantile error over `probs` in exact SDs; SD as a relative error; the L1
# distance between densities; and the log marginal likelihood (VB-corrected
# and plain Laplace). Then each stage against its own exact target: the
# Laplace curvature (relative) and log marginal likelihood, the VB mean (in
# Laplace SDs) and the VB log marginal likelihood.
approx_errors <- function(fit, model) {
  tab <- fit$summary
  qq <- summary(fit, probs = probs)$table
  q_sn <- qq[1, paste0(100 * probs, "%")]
  lap_sd <- sqrt(fit$cov[1, 1])
  q_laplace <- qnorm(probs, fit$mode[[1]], lap_sd)
  q_exact <- model$quantile(probs)

  tt <- seq(
    model$mode - 12 * model$sd,
    model$mode + 12 * model$sd,
    length.out = 4001
  )
  h <- tt[2] - tt[1]
  f_exact <- model$density(tt)
  f_sn <- dsn_ref(tt, fit$sn[1, "xi"], fit$sn[1, "omega"], fit$sn[1, "alpha"])
  f_laplace <- dnorm(tt, fit$mode[[1]], lap_sd)

  list(
    mean_sn = (tab[1, "Mean"] - model$mean) / model$sd,
    mean_laplace = (fit$mode[[1]] - model$mean) / model$sd,
    mean_vb = (fit$mean_vb[[1]] - model$mean) / model$sd,
    sd_sn = tab[1, "SD"] / model$sd - 1,
    sd_laplace = lap_sd / model$sd - 1,
    quantile_sn = max(abs(q_sn - q_exact)) / model$sd,
    quantile_laplace = max(abs(q_laplace - q_exact)) / model$sd,
    mode_sn = (tab[1, "Mode"] - model$mode) / model$sd,
    mode_laplace = (fit$mode[[1]] - model$mode) / model$sd,
    l1_sn = sum(abs(f_sn - f_exact)) * h,
    l1_laplace = sum(abs(f_laplace - f_exact)) * h,
    logml = fit$logml - model$logml,
    logml_laplace = fit$logml_laplace - model$logml,
    curvature = 1 / (fit$cov[1, 1] * model$laplace_H) - 1,
    laplace_vs_exact = fit$logml_laplace - model$laplace_logml,
    vb_vs_exact = (fit$mean_vb[[1]] - model$vb_mean) * sqrt(model$laplace_H),
    vb_logml_vs_exact = fit$logml - model$vb_logml
  )
}

# The ground truth itself ------------------------------------------------------

test_that("the exact references agree with numerical quadrature", {
  models <- list(normal, bb_skewed, bb_symmetric, gp_skewed, gp_symmetric)
  for (model in models) {
    # logpost - logml is the exact log density of theta, pointwise, so the
    # logpost includes every constant and the Jacobian
    tt <- model$mode + seq(-4, 4, by = 0.5) * model$sd
    expect_equal(
      exp(logpost_vec(model, tt) - model$logml),
      model$density(tt),
      tolerance = 1e-10
    )

    # The log marginal likelihood as the log integral of exp(logpost)
    lp_max <- model$logpost(c(theta = model$mode))
    int_lp <- integrate_posterior(
      model,
      function(t) exp(logpost_vec(model, t) - lp_max)
    )
    expect_equal(lp_max + log(int_lp), model$logml, tolerance = 1e-8)

    # Moments, mode and quantiles of the exact density
    expect_equal(integrate_posterior(model, model$density), 1, tolerance = 1e-8)
    mean_q <- integrate_posterior(model, function(t) t * model$density(t))
    expect_equal(mean_q, model$mean, tolerance = 1e-8)
    var_q <- integrate_posterior(
      model,
      function(t) (t - model$mean)^2 * model$density(t)
    )
    expect_equal(sqrt(var_q), model$sd, tolerance = 1e-8)
    mode_q <- optimize(
      function(t) log(model$density(t)),
      model$mode + c(-1, 1) * model$sd,
      maximum = TRUE,
      tol = 1e-10
    )$maximum
    expect_equal(mode_q, model$mode, tolerance = 1e-6)
    lower <- model$mode - 15 * model$sd
    for (p in probs) {
      cdf <- integrate(
        model$density,
        lower,
        model$quantile(p),
        rel.tol = 1e-10
      )$value
      expect_equal(cdf, p, tolerance = 1e-8)
    }

    # The analytic gradient, and the curvature at the mode, by central
    # differences of logpost
    for (t in tt[c(3, 9, 15)]) {
      h <- 1e-5 * model$sd
      fd <- (model$logpost(c(theta = t + h)) -
        model$logpost(c(theta = t - h))) /
        (2 * h)
      expect_equal(model$grad(c(theta = t)), fd, tolerance = 1e-7)
    }
    h <- 1e-3 * model$sd
    curv <- -(model$logpost(c(theta = model$mode + h)) -
      2 * lp_max +
      model$logpost(c(theta = model$mode - h))) /
      h^2
    expect_equal(model$laplace_H, curv, tolerance = 1e-5)

    # The VB mean: the expected gradient under N(vb_mean, 1 / H) vanishes
    s <- 1 / sqrt(model$laplace_H)
    expected_grad <- integrate(
      function(z) {
        vapply(
          z,
          function(zi) model$grad(c(theta = model$vb_mean + s * zi)),
          1
        ) *
          dnorm(z)
      },
      -12,
      12,
      rel.tol = 1e-12
    )$value
    expect_lt(abs(expected_grad) * s, 1e-9)
  }
})

test_that("the gamma-Poisson marginal likelihood matches quadrature on the natural scale", {
  for (pars in list(
    list(y = c(0, 1, 0, 2, 1), shape = 1, rate = 1),
    list(
      y = c(5, 7, 4, 6, 9, 3, 6, 8, 5, 7, 6, 4, 7, 5, 8, 6, 5, 7, 6, 4),
      shape = 2,
      rate = 0.5
    )
  )) {
    model <- gamma_poisson_model(pars$y, pars$shape, pars$rate)
    joint <- function(lambda) {
      loglik <- vapply(
        lambda,
        function(l) sum(dpois(pars$y, l, log = TRUE)),
        numeric(1)
      )
      # Scaled by the closed form, so the integral is one when it is right
      exp(
        loglik +
          dgamma(lambda, shape = pars$shape, rate = pars$rate, log = TRUE) -
          model$logml
      )
    }
    upper <- model$natural_quantile(1 - 1e-12)
    int <- integrate(joint, 0, upper, rel.tol = 1e-10)$value
    expect_equal(int, 1, tolerance = 1e-8)
  }
})

# (a) Normal mean, known variance --------------------------------------------

# Everything is exact for a Gaussian posterior, up to numerical precision:
# `mode_tol` (in SDs) for the optimiser and `laplace_tol` for the curvature,
# the log marginal likelihoods and the VB mean. The skew-normal fitted to an
# exactly quadratic scan is scale-free, so its tolerances are fixed.
expect_exact_gaussian <- function(fit, model, mode_tol, laplace_tol = 1e-10) {
  err <- approx_errors(fit, model)
  diag <- summary(fit)$diagnostics

  # Laplace: observed mode error 6e-7 SDs with finite differences and 1e-15
  # with the gradient; curvature relative error 4e-12 with the gradient and,
  # with finite differences, 4e-11 on x86-64 and 1e-10 on macOS ARM64, whose
  # rounding differs; log marginal likelihood 2e-11.
  expect_lt(abs(err$mode_laplace), mode_tol)
  expect_lt(abs(err$curvature), laplace_tol)
  expect_lt(abs(err$logml_laplace), laplace_tol)

  # The skew-normal reproduces the Gaussian. Observed errors over the data
  # sets and scales below: mean 5e-10 and SD up to 2.7e-9 (relative),
  # quantiles 2e-8 and mode 1.6e-8 (SDs), L1 distance 8e-9, NMAD 7.5e-9.
  expect_lt(abs(err$mean_sn), 1e-8)
  expect_lt(abs(err$sd_sn), 1e-8)
  expect_lt(err$quantile_sn, 5e-8)
  expect_lt(abs(err$mode_sn), 5e-8)
  expect_lt(err$l1_sn, 2e-8)
  expect_lt(diag[1, "nmad"], 2e-8)
  # alpha is only weakly identified near 0, since SN(xi, omega, alpha) with
  # small alpha differs from a shifted Gaussian at order alpha^3; observed
  # |alpha| up to 0.0054 across several data sets, scales and priors.
  expect_lt(abs(fit$sn[1, "alpha"]), 0.01)
  # With one parameter there is no other parameter whose volume changes along
  # the scan, so the volume correction vanishes (observed exactly 0).
  expect_lt(abs(diag[1, "gamma1"]), 1e-8)

  # The VB mean is the exact mean, so the shift only undoes the optimiser's
  # error in the mode. Observed VB mean error up to 6e-12 SDs.
  expect_lt(abs(diag[1, "vb_shift_sigma"]), mode_tol)
  expect_lt(abs(err$mean_vb), laplace_tol)
  expect_lt(abs(diag[1, "sn_vs_vb"]), 1e-8)
  # Laplace is exact for a Gaussian, and the VB term is 0.5 shift^2
  expect_lt(abs(err$logml), laplace_tol)
}

test_that("a Gaussian posterior is recovered exactly with finite differences", {
  fit <- inla_vanilla(normal$logpost, init = c(mu = 0), silent = TRUE)
  expect_exact_gaussian(fit, normal, mode_tol = 2e-6, laplace_tol = 1e-9)
})

test_that("a Gaussian posterior is recovered exactly with an analytic gradient", {
  fit <- inla_vanilla(
    normal$logpost,
    init = c(mu = 0),
    grad = normal$grad,
    silent = TRUE
  )
  expect_exact_gaussian(fit, normal, mode_tol = 1e-12)
})

test_that("a Gaussian posterior is recovered exactly at small and large scales", {
  # The same data in units a thousand times smaller and larger: posterior SDs
  # of 4e-4 and 4e2. Observed mode errors up to 1.9e-6 SDs, and curvature and
  # log marginal likelihood errors up to 2e-8.
  set.seed(1)
  z <- rnorm(12, 3, 1.5)
  for (k in c(1e-3, 1e3)) {
    model <- normal_model(z * k, sigma = 1.5 * k, m0 = k, s0 = 2 * k)
    for (gr in list(NULL, model$grad)) {
      fit <- inla_vanilla(
        model$logpost,
        init = c(mu = k),
        grad = gr,
        silent = TRUE
      )
      expect_exact_gaussian(fit, model, mode_tol = 5e-6, laplace_tol = 5e-8)
    }
  }
})

# House prices in dollars: sigma = 1e5 known and a N(2e5, 2e5^2) prior, so the
# posterior SD is 2.9e4, and the search starts from the prior mean, 4.3
# posterior SDs below the mode. Nothing here is extreme for data in natural
# units. A fit with a well-scaled mode search (control$nlminb$scale = 1 / 3e4)
# has mode error 2e-15 SDs, quantile error 2e-8 SDs, curvature error 2e-6 and
# log marginal likelihood error 9e-7, within the tolerances used below.
house_y <- 100 *
  c(2374, 3184, 2164, 4595, 3330, 2180, 3487, 3738, 3576, 2695, 4512, 3390)
house_prices <- normal_model(house_y, sigma = 1e5, m0 = 2e5, s0 = 2e5)

# The mode and the summary first, so that a failure reads clearly, then
# everything else once those are right
expect_house_prices_exact <- function(fit) {
  err <- approx_errors(fit, house_prices)
  expect_lt(abs(err$mode_laplace), 1e-5)
  expect_lt(err$quantile_sn, 1e-6)
  if (abs(err$mode_laplace) < 1e-5 && err$quantile_sn < 1e-6) {
    expect_exact_gaussian(
      fit,
      house_prices,
      mode_tol = 1e-5,
      laplace_tol = 1e-5
    )
  }
}

test_that("a Gaussian posterior in dollars is recovered exactly with finite differences", {
  # Fails (package bug): the first pass of the finite-difference Hessian steps
  # a fixed 1e-4, where the second difference of logpost is pure rounding
  # error (H h^2 is 1.2e-17 against |logpost| eps of 3.7e-14), so the rough
  # Hessian is often not positive definite and inla_vanilla() aborts with
  # "The Hessian at the posterior mode is not positive definite".
  fit <- NULL
  expect_no_error(
    fit <- inla_vanilla(
      house_prices$logpost,
      init = c(mu = 2e5),
      silent = TRUE
    )
  )
  if (!is.null(fit)) {
    expect_house_prices_exact(fit)
  }
})

test_that("a Gaussian posterior in dollars is recovered exactly with an analytic gradient", {
  # Fails (package bug): nlminb() starts from an identity Hessian, so its first
  # step is the raw gradient, 1.5e-4 dollars against a posterior SD of 2.9e4,
  # and it stops at once with "X-convergence (3)". inla_vanilla() reports
  # success with no warning, although the Newton step at the "mode"
  # (diagnostics$marginals$mode_shift_sigma) is 4.3 SDs. The scan is then
  # centred 4.3 SDs from the mode, and the summary quantiles are off by 2
  # SDs.
  fit <- inla_vanilla(
    house_prices$logpost,
    init = c(mu = 2e5),
    grad = house_prices$grad,
    silent = TRUE
  )
  expect_house_prices_exact(fit)
})

test_that("the default summary of a Gaussian posterior matches its exact quantiles", {
  fit <- inla_vanilla(normal$logpost, init = c(mu = 0), silent = TRUE)
  exact <- c(
    normal$mean,
    normal$sd,
    normal$quantile(c(0.025, 0.5, 0.975)),
    normal$mode
  )
  expect_equal(
    colnames(fit$summary),
    c("Mean", "SD", "2.5%", "50%", "97.5%", "Mode")
  )
  expect_equal(rownames(fit$summary), "mu")
  expect_equal(unname(fit$summary[1, ]), exact, tolerance = 1e-7)
})

test_that("without the VB correction the log marginal likelihood is the Laplace one", {
  fit <- inla_vanilla(
    bb_skewed$logpost,
    init = c(theta = 0),
    vb = FALSE,
    silent = TRUE
  )
  expect_identical(fit$logml, fit$logml_laplace)
  expect_equal(fit$mean_vb, fit$mode)
  # No VB correction, so no VB diagnostics, as in INLAvaan's diagnostics()
  vb_cols <- c("kld", "vb_shift", "vb_shift_sigma", "sn_vs_vb")
  expect_true(all(is.na(summary(fit)$diagnostics[vb_cols])))
  # The skew-normal scan is centred at the mode, so the VB correction does not
  # move the marginals
  fit_vb <- inla_vanilla(
    bb_skewed$logpost,
    init = c(theta = 0),
    silent = TRUE
  )
  expect_identical(fit$sn, fit_vb$sn)
  expect_identical(fit$summary, fit_vb$summary)
  expect_identical(fit$logml_laplace, fit_vb$logml_laplace)
})

# (b) and (c): skewed posteriors ---------------------------------------------

# The Laplace and VB stages against their exact targets. Observed, over the
# four models below:
# * mode error up to 3e-6 SDs with finite differences and 8e-9 with the
#   gradient;
# * curvature relative error up to 1.7e-6 and 1.1e-9;
# * Laplace log marginal likelihood error up to 8.3e-7 and 5.4e-10.
# `vb_tol` bounds the VB mean's error in Laplace SDs, which comes from the
# cubature rule (for the log-gamma the degree-5 rule's relative error in
# E[exp(theta)] is s^6 / 120 for Laplace SD s, so 6e-5 at S = 5, and 1.4e-4
# SDs in the VB mean). Without a gradient the VB mean comes from Stein's
# identity on the three-point rule, exact for log-densities of degree four
# rather than six, and its error is larger. `vb_logml_tol` bounds the error
# of fit$logml against the exact Laplace plus VB term.
expect_stages_exact <- function(fit, model, gradient, vb_tol, vb_logml_tol) {
  err <- approx_errors(fit, model)
  if (gradient) {
    expect_lt(abs(err$mode_laplace), 3e-8)
    expect_lt(abs(err$curvature), 3e-9)
    expect_lt(abs(err$laplace_vs_exact), 1.5e-9)
  } else {
    expect_lt(abs(err$mode_laplace), 1e-5)
    expect_lt(abs(err$curvature), 5e-6)
    expect_lt(abs(err$laplace_vs_exact), 2.5e-6)
  }
  expect_lt(abs(err$vb_vs_exact), vb_tol)
  expect_lt(abs(err$vb_logml_vs_exact), vb_logml_tol)
  # The VB term is 0.5 shift' H shift with the fit's own H
  expect_equal(
    fit$logml - fit$logml_laplace,
    0.5 * (fit$mean_vb[[1]] - fit$mode[[1]])^2 / fit$cov[1, 1],
    tolerance = 1e-8
  )
  # No volume correction with one parameter (observed at most 1e-14)
  expect_lt(abs(summary(fit)$diagnostics[1, "gamma1"]), 1e-8)
  # The summary is that of the fitted skew-normal (observed at most 2.7e-8 SDs)
  expect_lt(summary_vs_sn(fit, model), 1e-7)
  err
}

fit_both_ways <- function(model, init) {
  list(
    finite_differences = inla_vanilla(
      model$logpost,
      init = init,
      silent = TRUE
    ),
    analytic = inla_vanilla(
      model$logpost,
      init = init,
      grad = model$grad,
      silent = TRUE
    )
  )
}

# (b) Beta-binomial on the logit scale ---------------------------------------

test_that("a skewed beta-binomial posterior is recovered on the logit scale", {
  fits <- fit_both_ways(bb_skewed, c(theta = 0))
  for (way in names(fits)) {
    fit <- fits[[way]]
    # VB mean: observed error 3.5e-5 Laplace SDs with the gradient and 5.3e-3
    # without; VB log marginal likelihood: observed 6.9e-6 and 1.1e-3
    analytic <- way == "analytic"
    err <- expect_stages_exact(
      fit,
      bb_skewed,
      gradient = analytic,
      vb_tol = if (analytic) 2e-4 else 0.015,
      vb_logml_tol = if (analytic) 4e-5 else 3e-3
    )

    # Mean: observed SN error 0.0022 SDs, Laplace (the mode) 0.20 SDs
    expect_lt(abs(err$mean_sn), 0.007)
    expect_lt(abs(err$mean_sn), abs(err$mean_laplace))
    # SD: observed relative error 0.0070 for SN, 0.071 for Laplace
    expect_lt(abs(err$sd_sn), 0.02)
    expect_lt(abs(err$sd_sn), abs(err$sd_laplace))
    # Quantiles 2.5% to 97.5%: observed worst SN error 0.031 SDs (at 2.5%),
    # Laplace 0.55 SDs
    expect_lt(err$quantile_sn, 0.08)
    expect_lt(err$quantile_sn, err$quantile_laplace)
    # Mode of the skew-normal: observed error 0.021 SDs
    expect_lt(abs(err$mode_sn), 0.06)
    # Density: observed L1 distance 0.011 for SN, 0.12 for Laplace
    expect_lt(err$l1_sn, 0.03)
    expect_lt(err$l1_sn, err$l1_laplace)

    # The skew-normal is skewed the right way (left skew; observed alpha
    # -1.81), and the VB mean moves the right way: observed shift -0.19 SDs,
    # VB mean error against the posterior mean 0.023 SDs
    expect_lt(fit$sn[1, "alpha"], -1)
    expect_lt(fit$mean_vb[[1]], fit$mode[[1]])
    expect_lt(abs(err$mean_vb), 0.06)
    expect_lt(abs(err$mean_vb), abs(err$mean_laplace))

    # Log marginal likelihood: observed error -0.0097 with the VB term,
    # -0.029 for plain Laplace
    expect_lt(abs(err$logml), 0.025)
    expect_lt(abs(err$logml), abs(err$logml_laplace))
  }
})

test_that("a near-symmetric beta-binomial posterior is recovered on the logit scale", {
  fits <- fit_both_ways(bb_symmetric, c(theta = 0))
  for (way in names(fits)) {
    fit <- fits[[way]]
    # VB mean: observed error 2.1e-5 Laplace SDs with the gradient and 2.3e-4
    # without; VB log marginal likelihood: observed 2.8e-7 and 3.2e-6
    analytic <- way == "analytic"
    err <- expect_stages_exact(
      fit,
      bb_symmetric,
      gradient = analytic,
      vb_tol = if (analytic) 1e-4 else 7e-4,
      vb_logml_tol = if (analytic) 1e-6 else 1e-5
    )

    # Mean: observed SN error 0.00012 SDs, Laplace 0.014 SDs
    expect_lt(abs(err$mean_sn), 4e-4)
    expect_lt(abs(err$mean_sn), abs(err$mean_laplace))
    # SD: observed relative error 0.0006 for SN, 0.012 for Laplace
    expect_lt(abs(err$sd_sn), 0.002)
    expect_lt(abs(err$sd_sn), abs(err$sd_laplace))
    # Quantiles: observed worst SN error 0.0053 SDs, Laplace 0.053 SDs
    expect_lt(err$quantile_sn, 0.015)
    expect_lt(err$quantile_sn, err$quantile_laplace)
    # Mode of the skew-normal: observed error 0.0016 SDs
    expect_lt(abs(err$mode_sn), 0.005)
    # Density: observed L1 distance 0.0042 for SN, 0.011 for Laplace
    expect_lt(err$l1_sn, 0.012)
    expect_lt(err$l1_sn, err$l1_laplace)
    # VB mean: observed shift -0.013 SDs, error 0.0005 SDs
    expect_lt(abs(err$mean_vb), 0.0015)
    expect_lt(abs(err$mean_vb), abs(err$mean_laplace))
    # Log marginal likelihood: observed error -0.0057 with the VB term and
    # -0.0057 for plain Laplace. The location term is negligible here, and the
    # remaining error comes from the tails, which neither corrects.
    expect_lt(abs(err$logml), 0.015)
  }
})

test_that("the fit does not depend on the starting value", {
  # Starts 3.6 SDs below and 8.1 SDs above the mode. Observed: the summaries
  # agree with the start at 0 to 4.2e-8 SDs, the log marginal likelihoods to
  # 1.5e-8
  ref <- inla_vanilla(bb_skewed$logpost, init = c(theta = 0), silent = TRUE)
  for (init in c(-4, 4)) {
    fit <- inla_vanilla(
      bb_skewed$logpost,
      init = c(theta = init),
      silent = TRUE
    )
    expect_lt(
      max(abs(fit$summary - ref$summary)) / bb_skewed$sd,
      1e-7
    )
    expect_lt(abs(fit$logml - ref$logml), 5e-8)
  }
})

test_that("copula draws reproduce the exact beta posterior of p", {
  fit <- inla_vanilla(bb_skewed$logpost, init = c(theta = 0), silent = TRUE)
  set.seed(42)
  p <- sampling(fit, nsamp = 20000, fun = function(theta) plogis(theta))
  expect_equal(dim(p), c(20000L, 1L))
  expect_true(all(p > 0 & p < 1))

  # The tolerances here are set by the Monte Carlo noise rather than by the
  # error at this seed, so that a change of seed or of RNG does not break them.
  # Mean: Monte Carlo SE about 0.0006; observed error 0.0004 at seed 42, and
  # 0.0014 (99th percentile) and 0.0017 (max) over seeds 1 to 300.
  expect_lt(abs(mean(p) - bb_skewed$natural_mean), 0.003)
  # Quantiles in exact SDs of p: observed worst error 0.0011 at seed 42, and
  # 0.068 (99th percentile) and 0.084 (max) over seeds 1 to 300, nearly all
  # Monte Carlo noise, whose SE is about 0.02 at 97.5%. The plain Laplace
  # quantiles mapped to p are off by up to 0.49.
  pp <- c(0.025, 0.5, 0.975)
  q_draws <- quantile(p, pp, names = FALSE)
  q_exact <- bb_skewed$natural_quantile(pp)
  q_laplace <- plogis(qnorm(pp, fit$mode[[1]], sqrt(fit$cov[1, 1])))
  err_draws <- max(abs(q_draws - q_exact)) / bb_skewed$natural_sd
  err_laplace <- max(abs(q_laplace - q_exact)) / bb_skewed$natural_sd
  expect_lt(err_draws, 0.08)
  expect_lt(err_draws, err_laplace)
})

# (c) Gamma-Poisson on the log scale -----------------------------------------

test_that("a skewed gamma-Poisson posterior is recovered on the log scale", {
  fits <- fit_both_ways(gp_skewed, c(log_lambda = 0))
  for (way in names(fits)) {
    fit <- fits[[way]]
    # VB mean: observed error 1.4e-4 Laplace SDs with the gradient, the error
    # of the cubature rule, and 4.3e-3 without; VB log marginal likelihood:
    # observed 3e-5 and 9.6e-4
    analytic <- way == "analytic"
    err <- expect_stages_exact(
      fit,
      gp_skewed,
      gradient = analytic,
      vb_tol = if (analytic) 4e-4 else 0.013,
      vb_logml_tol = if (analytic) 1e-4 else 3e-3
    )

    # Mean: observed SN error 0.0024 SDs, Laplace 0.22 SDs
    expect_lt(abs(err$mean_sn), 0.007)
    expect_lt(abs(err$mean_sn), abs(err$mean_laplace))
    # SD: observed relative error 0.0033 for SN, 0.049 for Laplace
    expect_lt(abs(err$sd_sn), 0.01)
    expect_lt(abs(err$sd_sn), abs(err$sd_laplace))
    # Quantiles: observed worst SN error 0.012 SDs, Laplace 0.53 SDs
    expect_lt(err$quantile_sn, 0.035)
    expect_lt(err$quantile_sn, err$quantile_laplace)
    # Mode of the skew-normal: observed error 0.024 SDs
    expect_lt(abs(err$mode_sn), 0.06)
    # Density: observed L1 distance 0.011 for SN, 0.12 for Laplace
    expect_lt(err$l1_sn, 0.03)
    expect_lt(err$l1_sn, err$l1_laplace)

    # Left skew (observed alpha -1.90), and a VB mean that moves left:
    # observed shift -0.22 SDs, VB mean error against the posterior mean
    # 0.0072 SDs
    expect_lt(fit$sn[1, "alpha"], -1)
    expect_lt(fit$mean_vb[[1]], fit$mode[[1]])
    expect_lt(abs(err$mean_vb), 0.02)
    expect_lt(abs(err$mean_vb), abs(err$mean_laplace))
    # The skew-normal and VB means agree (observed gap 0.005 SDs)
    expect_lt(abs(summary(fit)$diagnostics[1, "sn_vs_vb"]), 0.015)

    # Log marginal likelihood: observed error 0.0083 with the VB term (which
    # overshoots), -0.017 for plain Laplace
    expect_lt(abs(err$logml), 0.025)
    expect_lt(abs(err$logml), abs(err$logml_laplace))
  }
})

test_that("a near-symmetric gamma-Poisson posterior is recovered on the log scale", {
  fits <- fit_both_ways(gp_symmetric, c(log_lambda = 0))
  for (way in names(fits)) {
    fit <- fits[[way]]
    # VB mean: observed error 5.2e-8 Laplace SDs with the gradient (the
    # cubature rule is nearly exact at S = 120) and 3.8e-5 without; VB log
    # marginal likelihood: observed 3.6e-8 and 1.8e-6
    analytic <- way == "analytic"
    err <- expect_stages_exact(
      fit,
      gp_symmetric,
      gradient = analytic,
      vb_tol = if (analytic) 1.5e-7 else 1.2e-4,
      vb_logml_tol = if (analytic) 1e-7 else 6e-6
    )

    # Mean: observed SN error 8e-5 SDs, Laplace 0.046 SDs
    expect_lt(abs(err$mean_sn), 3e-4)
    expect_lt(abs(err$mean_sn), abs(err$mean_laplace))
    # SD: observed relative error 0.0002 for SN, 0.0021 for Laplace
    expect_lt(abs(err$sd_sn), 7e-4)
    expect_lt(abs(err$sd_sn), abs(err$sd_laplace))
    # Quantiles: observed worst SN error 0.0026 SDs, Laplace 0.093 SDs
    expect_lt(err$quantile_sn, 0.008)
    expect_lt(err$quantile_sn, err$quantile_laplace)
    # Mode of the skew-normal: observed error 0.0004 SDs
    expect_lt(abs(err$mode_sn), 0.0012)
    # Density: observed L1 distance 0.0025 for SN, 0.024 for Laplace
    expect_lt(err$l1_sn, 0.0075)
    expect_lt(err$l1_sn, err$l1_laplace)
    # VB mean: observed shift -0.045 SDs, error 6e-5 SDs
    expect_lt(abs(err$mean_vb), 2e-4)
    expect_lt(abs(err$mean_vb), abs(err$mean_laplace))
    # Log marginal likelihood: observed error 0.00035 with the VB term,
    # -0.00069 for plain Laplace
    expect_lt(abs(err$logml), 0.001)
    expect_lt(abs(err$logml), abs(err$logml_laplace))
  }
})

test_that("copula draws reproduce the exact gamma posterior of lambda", {
  fit <- inla_vanilla(
    gp_skewed$logpost,
    init = c(log_lambda = 0),
    silent = TRUE
  )
  set.seed(42)
  lambda <- sampling(fit, nsamp = 20000, fun = function(theta) exp(theta))

  # In exact SDs of lambda, with tolerances set by the Monte Carlo noise. Mean:
  # Monte Carlo SE about 0.007; observed error 0.0013 at seed 42, and 0.017
  # (99th percentile) and 0.021 (max) over seeds 1 to 300.
  err_mean <- (mean(lambda) - gp_skewed$natural_mean) / gp_skewed$natural_sd
  expect_lt(abs(err_mean), 0.03)
  # Quantiles: observed worst error 0.017 at seed 42, and 0.086 (99th
  # percentile) and 0.104 (max) over seeds 1 to 300, nearly all Monte Carlo
  # noise, whose SE is about 0.03 at 97.5% in the long right tail. The plain
  # Laplace quantiles mapped to lambda are off by up to 0.79.
  pp <- c(0.025, 0.5, 0.975)
  q_draws <- quantile(lambda, pp, names = FALSE)
  q_exact <- gp_skewed$natural_quantile(pp)
  q_laplace <- exp(qnorm(pp, fit$mode[[1]], sqrt(fit$cov[1, 1])))
  err_draws <- max(abs(q_draws - q_exact)) / gp_skewed$natural_sd
  err_laplace <- max(abs(q_laplace - q_exact)) / gp_skewed$natural_sd
  expect_lt(err_draws, 0.12)
  expect_lt(err_draws, err_laplace)
})
