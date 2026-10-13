# End-to-end tests of inla_vanilla() on one-parameter conjugate models, whose
# posteriors and marginal likelihoods are known in closed form. Every logpost
# carries all of its normalising constants (and the Jacobian of the
# reparameterisation), so fit$logml is directly comparable with the exact log
# marginal likelihood. Location errors are measured in exact posterior SDs and
# SD errors as relative errors.
#
# The tolerances on the non-Gaussian cases were calibrated by printing the
# actual errors first; each comment states the observed value. The results
# with finite differences and with an analytic gradient agree to about 1e-7,
# so the same tolerances serve both.

probs <- c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975)

# Models ----------------------------------------------------------------------
#
# Each builder returns the logpost (a function of the named vector theta, as
# inla_vanilla() passes it), its gradient, the exact posterior of theta
# (density, mean, SD, mode and quantile function) and the exact log marginal
# likelihood.

# y_i ~ N(mu, sigma^2) with sigma known, mu ~ N(m0, s0^2). The posterior of mu
# is Gaussian, and marginally y ~ N_n(m0 1, sigma^2 I + s0^2 1 1').
normal_model <- function(y, sigma, m0, s0) {
  n <- length(y)
  prec <- 1 / s0^2 + n / sigma^2
  post_mean <- (m0 / s0^2 + sum(y) / sigma^2) / prec
  post_sd <- sqrt(1 / prec)

  V <- sigma^2 * diag(n) + s0^2 * matrix(1, n, n)
  U <- chol(V)
  r <- backsolve(U, y - m0, transpose = TRUE)
  logml <- -n / 2 * log(2 * pi) - sum(log(diag(U))) - 0.5 * sum(r^2)

  list(
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
    logml = logml
  )
}

# y ~ Binomial(n, p), p ~ Beta(a, b), on theta = logit(p). The Jacobian is
# dp / dtheta = p (1 - p). Exactly p | y ~ Beta(A, B) with A = a + y and
# B = b + n - y, so theta | y is a logit-beta: E = digamma(A) - digamma(B),
# Var = trigamma(A) + trigamma(B), mode log(A / B), quantiles
# qlogis(qbeta(...)). The marginal likelihood is the beta-binomial pmf.
beta_binomial_model <- function(y, n, a, b) {
  A <- a + y
  B <- b + n - y
  list(
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
    natural_quantile = function(p) qbeta(p, A, B),
    natural_mean = A / (A + B),
    natural_sd = sqrt(A * B / ((A + B)^2 * (A + B + 1)))
  )
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
gamma_poisson_model <- function(y, shape, rate) {
  n <- length(y)
  S <- shape + sum(y)
  R <- rate + n
  log_prior_const <- shape * log(rate) - lgamma(shape)
  logml <- log_prior_const - sum(lfactorial(y)) + lgamma(S) - S * log(R)
  list(
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
    natural_quantile = function(p) qgamma(p, shape = S, rate = R),
    natural_mean = S / R,
    natural_sd = sqrt(S) / R
  )
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

# Errors of a fit against the exact posterior, for the skew-normal marginal
# and for the Laplace Gaussian N(mode, cov): mean, mode and the largest
# quantile error over `probs` in exact SDs; SD as a relative error; the L1
# distance between densities; and the log marginal likelihood (VB-corrected
# and plain Laplace).
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
  f_sn <- dsnorm(tt, fit$sn[1, "xi"], fit$sn[1, "omega"], fit$sn[1, "alpha"])
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
    logml_laplace = fit$logml_laplace - model$logml
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

expect_exact_gaussian <- function(fit, model) {
  err <- approx_errors(fit, model)
  # The Laplace step is exact for a Gaussian posterior, up to the precision of
  # the optimiser (observed mode error 6e-7 SDs) and of the Hessian (observed
  # relative error 1e-11).
  expect_lt(abs(err$mode_laplace), 1e-4)
  expect_equal(fit$cov[1, 1], model$sd^2, tolerance = 1e-6)
  # The skew-normal reproduces the Gaussian (observed errors below 2e-8 SDs)
  expect_lt(abs(err$mean_sn), 1e-6)
  expect_lt(abs(err$sd_sn), 1e-6)
  expect_lt(err$quantile_sn, 1e-6)
  expect_lt(abs(err$mode_sn), 1e-6)
  expect_lt(err$l1_sn, 1e-6)
  # alpha is only weakly identified near 0, since SN(xi, omega, alpha) with
  # small alpha differs from a shifted Gaussian at order alpha^3; observed
  # |alpha| about 0.005 across several data sets and priors.
  expect_lt(abs(fit$sn[1, "alpha"]), 0.015)
  # The VB mean is the exact mean, so the shift only undoes the optimiser's
  # error in the mode (observed 6e-7 SDs).
  shift <- summary(fit)$diagnostics[1, "vb_shift"]
  expect_lt(abs(shift), 1e-4)
  expect_lt(abs(err$mean_vb), 1e-6)
  expect_lt(abs(summary(fit)$diagnostics[1, "sn_vs_vb"]), 1e-6)
  # Laplace is exact for a Gaussian, and the VB term is 0.5 shift^2, below
  # 1e-8 (observed errors 1e-11).
  expect_lt(abs(err$logml_laplace), 1e-6)
  expect_lt(abs(err$logml), 1e-6)
}

test_that("a Gaussian posterior is recovered exactly with finite differences", {
  fit <- inla_vanilla(normal$logpost, init = c(mu = 0), verbose = FALSE)
  expect_exact_gaussian(fit, normal)
})

test_that("a Gaussian posterior is recovered exactly with an analytic gradient", {
  fit <- inla_vanilla(
    normal$logpost,
    init = c(mu = 0),
    grad = normal$grad,
    verbose = FALSE
  )
  expect_exact_gaussian(fit, normal)
})

test_that("the default summary of a Gaussian posterior matches its exact quantiles", {
  fit <- inla_vanilla(normal$logpost, init = c(mu = 0), verbose = FALSE)
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
    verbose = FALSE
  )
  expect_identical(fit$logml, fit$logml_laplace)
  expect_equal(fit$mean_vb, fit$mode)
  expect_true(all(summary(fit)$diagnostics$vb_shift == 0))
})

# (b) Beta-binomial on the logit scale ---------------------------------------

test_that("a skewed beta-binomial posterior is recovered on the logit scale", {
  grads <- list(finite_differences = NULL, analytic = bb_skewed$grad)
  for (gr in grads) {
    fit <- inla_vanilla(
      bb_skewed$logpost,
      init = c(theta = 0),
      grad = gr,
      verbose = FALSE
    )
    err <- approx_errors(fit, bb_skewed)

    # The Laplace mode is the exact mode (observed error 8e-7 SDs)
    expect_lt(abs(err$mode_laplace), 1e-4)

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

    # The skew-normal is skewed the right way (left skew), and the VB mean
    # moves the right way: observed shift -0.19 SDs, VB mean error 0.023 SDs
    expect_lt(fit$sn[1, "alpha"], -1)
    expect_lt(fit$mean_vb[[1]], fit$mode[[1]])
    expect_lt(abs(err$mean_vb), 0.06)
    expect_lt(abs(err$mean_vb), abs(err$mean_laplace))

    # Log marginal likelihood: observed error -0.0097 with the VB term,
    # -0.029 for plain Laplace
    expect_lt(abs(err$logml), 0.025)
    expect_lt(abs(err$logml), abs(err$logml_laplace))
    # The VB term is 0.5 shift' H shift
    expect_equal(
      fit$logml - fit$logml_laplace,
      0.5 * (fit$mean_vb[[1]] - fit$mode[[1]])^2 / fit$cov[1, 1],
      tolerance = 1e-8
    )
  }
})

test_that("a near-symmetric beta-binomial posterior is recovered on the logit scale", {
  fit <- inla_vanilla(
    bb_symmetric$logpost,
    init = c(theta = 0),
    verbose = FALSE
  )
  err <- approx_errors(fit, bb_symmetric)

  # Observed mode error 5e-9 SDs
  expect_lt(abs(err$mode_laplace), 1e-4)
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
})

test_that("copula draws reproduce the exact beta posterior of p", {
  fit <- inla_vanilla(bb_skewed$logpost, init = c(theta = 0), verbose = FALSE)
  set.seed(42)
  p <- sampling(fit, nsamp = 20000, fun = function(theta) plogis(theta))
  expect_equal(dim(p), c(20000L, 1L))
  expect_true(all(p > 0 & p < 1))

  # Monte Carlo SE of the mean is about 0.0006 (observed error 0.0003; up to
  # 0.0009 over other seeds)
  expect_lt(abs(mean(p) - bb_skewed$natural_mean), 0.003)
  # Quantiles in exact SDs of p: observed worst error 0.0011 (up to 0.042
  # over other seeds), nearly all Monte Carlo noise, whose SE is about 0.02 at
  # 97.5%. The plain Laplace quantiles mapped to p are off by up to 0.49.
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
  fit <- inla_vanilla(
    gp_skewed$logpost,
    init = c(log_lambda = 0),
    verbose = FALSE
  )
  err <- approx_errors(fit, gp_skewed)

  # Observed mode error 5e-9 SDs
  expect_lt(abs(err$mode_laplace), 1e-4)
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

  # Left skew, and a VB mean that moves left: observed shift -0.22 SDs, VB
  # mean error 0.0072 SDs
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
})

test_that("a near-symmetric gamma-Poisson posterior is recovered on the log scale", {
  fit <- inla_vanilla(
    gp_symmetric$logpost,
    init = c(log_lambda = 0),
    grad = gp_symmetric$grad,
    verbose = FALSE
  )
  err <- approx_errors(fit, gp_symmetric)

  # Observed mode error 8e-9 SDs with the gradient
  expect_lt(abs(err$mode_laplace), 1e-4)
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
})

test_that("copula draws reproduce the exact gamma posterior of lambda", {
  fit <- inla_vanilla(
    gp_skewed$logpost,
    init = c(log_lambda = 0),
    verbose = FALSE
  )
  set.seed(42)
  lambda <- sampling(fit, nsamp = 20000, fun = function(theta) exp(theta))

  # In exact SDs of lambda. Monte Carlo SE of the mean is about 0.007
  # (observed error 0.0013; up to 0.0086 over other seeds)
  err_mean <- (mean(lambda) - gp_skewed$natural_mean) / gp_skewed$natural_sd
  expect_lt(abs(err_mean), 0.03)
  # Quantiles: observed worst error 0.017 (up to 0.068 over other seeds),
  # nearly all Monte Carlo noise, whose SE is about 0.03 at 97.5% in the long
  # right tail. The plain Laplace quantiles mapped to lambda are off by up to
  # 0.79.
  pp <- c(0.025, 0.5, 0.975)
  q_draws <- quantile(lambda, pp, names = FALSE)
  q_exact <- gp_skewed$natural_quantile(pp)
  q_laplace <- exp(qnorm(pp, fit$mode[[1]], sqrt(fit$cov[1, 1])))
  err_draws <- max(abs(q_draws - q_exact)) / gp_skewed$natural_sd
  err_laplace <- max(abs(q_laplace - q_exact)) / gp_skewed$natural_sd
  expect_lt(err_draws, 0.12)
  expect_lt(err_draws, err_laplace)
})
