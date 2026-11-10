# Natural-scale reporting: a parameter named log_x, logit_x or atanh_x is
# reported as x by print(), summary(), coef(), plot() and sampling().
#
# The engine is untouched, so there are three things to check. (a) The names
# are read as documented. (b) Given a skew-normal for theta, the summary of
# psi = inverse(theta) is exact: the references are quadratures of the
# skew-normal density written in this file. (c) End to end, the natural-scale
# summaries match posteriors with a closed form on that scale: gamma-Poisson
# (lambda is gamma), beta-binomial (p is beta), an exactly Gaussian log
# parameter (log-normal) and a correlation by one-dimensional quadrature.
# Their tolerances are those of the skew-normal approximation itself, which
# the other test files measure on the scale of logpost.

# References on the natural scale ---------------------------------------------

sn_density_ref <- function(x, xi, omega, alpha) {
  z <- (x - xi) / omega
  2 * dnorm(z) * pnorm(alpha * z) / omega
}

# E h(theta) for theta ~ SN(xi, omega, alpha), over 12 SDs either side, in
# pieces so that integrate() cannot miss a narrow integrand
sn_expect <- function(h, xi, omega, alpha) {
  cuts <- xi + omega * seq(-12, 12, by = 1)
  sum(vapply(
    seq_len(length(cuts) - 1L),
    function(k) {
      integrate(
        function(t) h(t) * sn_density_ref(t, xi, omega, alpha),
        cuts[k],
        cuts[k + 1L],
        rel.tol = 1e-12
      )$value
    },
    numeric(1)
  ))
}

maps <- list(
  log = list(inverse = exp, forward = log, lower = 0, upper = Inf),
  logit = list(inverse = plogis, forward = qlogis, lower = 0, upper = 1),
  atanh = list(inverse = tanh, forward = atanh, lower = -1, upper = 1)
)

# Skew-normals from nearly Gaussian to strongly skewed, narrow and wide
sn_cases <- list(
  c(xi = 0.3, omega = 0.8, alpha = 3),
  c(xi = -1.2, omega = 0.4, alpha = -6),
  c(xi = 1, omega = 0.05, alpha = 20),
  c(xi = 0.5, omega = 1.5, alpha = 0),
  c(xi = -0.4, omega = 0.2, alpha = 0.5)
)

# (a) Reading the names --------------------------------------------------------

test_that("a prefix followed by an underscore and a name is a transform", {
  tr <- detect_transforms(c(
    "mu",
    "log_sigma",
    "logit_p",
    "atanh_rho",
    "log_tau[2]",
    "log_a_b"
  ))
  expect_identical(
    tr$transform,
    c("identity", "log", "logit", "atanh", "log", "log")
  )
  expect_identical(tr$name, c("mu", "sigma", "p", "rho", "tau[2]", "a_b"))
  expect_identical(tr$parameter, rownames(tr))
})

test_that("names that only resemble a prefix are left alone", {
  nms <- c(
    "logistic",
    "log",
    "log_",
    "logit_",
    "Log_sigma",
    "LOG_sigma",
    "log.sigma",
    "logsigma",
    "my_log_sigma",
    "tanh_rho",
    "logv_x"
  )
  tr <- expect_silent(detect_transforms(nms))
  expect_true(all(tr$transform == "identity"))
  expect_identical(tr$name, nms)
})

test_that("a name that would clash is left on the unconstrained scale", {
  expect_warning(
    tr <- detect_transforms(c("log_sigma", "sigma", "logit_p")),
    "unconstrained scale",
    class = "inlavanilla_warning"
  )
  expect_identical(tr$transform, c("identity", "identity", "logit"))
  expect_identical(tr$name, c("log_sigma", "sigma", "p"))

  # two transforms of one name clash with each other
  expect_warning(tr <- detect_transforms(c("log_x", "logit_x")), "log_x")
  expect_true(all(tr$transform == "identity"))

  # undoing one clash can make another
  expect_warning(tr <- detect_transforms(c("log_log_x", "log_x", "x")))
  expect_true(all(tr$transform == "identity"))
  expect_false(anyDuplicated(tr$name) > 0L)
})

test_that("detection can be switched off", {
  tr <- detect_transforms(c("mu", "log_sigma"), detect = FALSE)
  expect_true(all(tr$transform == "identity"))
  expect_identical(tr$name, c("mu", "log_sigma"))
})

# (b) The dictionary and the summary of a transformed skew-normal -------------

test_that("each transform inverts, increases and has the right Jacobian", {
  theta <- c(-4, -1.3, 0, 0.4, 2.5, 5)
  for (key in names(transform_dict)) {
    tf <- transform_dict[[key]]
    expect_equal(tf$forward(tf$inverse(theta)), theta, tolerance = 1e-9)
    expect_true(all(diff(tf$inverse(theta)) > 0))
    expect_equal(tf$inverse(theta), maps[[key]]$inverse(theta))
    # the log-derivative of the inverse, by central differences
    h <- 1e-5
    slope <- (tf$inverse(theta + h) - tf$inverse(theta - h)) / (2 * h)
    expect_equal(tf$log_jac(theta), log(slope), tolerance = 1e-7)
  }
  # far out, where 1 - tanh(theta)^2 underflows to zero in double precision
  expect_equal(transform_dict$atanh$log_jac(40), log(4) - 80)
  expect_equal(transform_dict$logit$log_jac(-800), -800)
})

test_that("the mean and SD of a transformed skew-normal are exact", {
  # Observed errors: 2e-15 at most, relative, for the closed form of the log
  # and for the trapezoid rule of the other two
  for (key in names(maps)) {
    for (sn in sn_cases) {
      inverse <- maps[[key]]$inverse
      m1 <- sn_expect(inverse, sn[["xi"]], sn[["omega"]], sn[["alpha"]])
      v <- sn_expect(
        function(t) (inverse(t) - m1)^2,
        sn[["xi"]],
        sn[["omega"]],
        sn[["alpha"]]
      )
      got <- natural_summary(
        sn[["xi"]],
        sn[["omega"]],
        sn[["alpha"]],
        key,
        probs = 0.5
      )
      expect_equal(got[1], m1, tolerance = 1e-10)
      expect_equal(got[2], sqrt(v), tolerance = 1e-8)
    }
  }
})

test_that("the SD of exp(theta) keeps its accuracy for a tiny scale", {
  # To first order the SD of exp(theta) is exp(xi) times that of theta, and
  # a difference of moments would lose it all at omega = 1e-9
  for (omega in c(1e-3, 1e-6, 1e-9)) {
    for (alpha in c(0, 4)) {
      delta <- alpha / sqrt(1 + alpha^2)
      sd_theta <- omega * sqrt(1 - 2 * delta^2 / pi)
      got <- lognormal_moments(2, omega, alpha)
      expect_equal(got$sd, exp(2) * sd_theta, tolerance = 10 * omega + 1e-7)
    }
  }
})

test_that("the quantiles of a transformed skew-normal map straight through", {
  probs <- c(0.01, 0.25, 0.5, 0.9, 0.999)
  for (key in names(maps)) {
    for (sn in sn_cases) {
      got <- natural_summary(
        sn[["xi"]],
        sn[["omega"]],
        sn[["alpha"]],
        key,
        probs
      )
      q <- got[2 + seq_along(probs)]
      # the probability below each, from the density of theta
      below <- vapply(
        q,
        function(qi) {
          sn_expect(
            function(t) as.numeric(t <= maps[[key]]$forward(qi)),
            sn[["xi"]],
            sn[["omega"]],
            sn[["alpha"]]
          )
        },
        numeric(1)
      )
      expect_equal(below, probs, tolerance = 1e-5)
    }
  }
})

test_that("the mode is that of the natural-scale density", {
  for (key in names(maps)) {
    for (sn in sn_cases) {
      xi <- sn[["xi"]]
      omega <- sn[["omega"]]
      alpha <- sn[["alpha"]]
      got <- natural_summary(xi, omega, alpha, key, probs = 0.5)[4]
      # The density of psi by the change of variables, with the derivative of
      # the forward map by central differences, maximised on a fine grid
      forward <- maps[[key]]$forward
      dens <- function(psi) {
        h <- 1e-6 * pmin(psi - maps[[key]]$lower, maps[[key]]$upper - psi, 1)
        sn_density_ref(forward(psi), xi, omega, alpha) *
          (forward(psi + h) - forward(psi - h)) /
          (2 * h)
      }
      lims <- maps[[key]]$inverse(xi + c(-9, 9) * omega)
      grid <- seq(lims[1], lims[2], length.out = 4e5 + 2)
      grid <- grid[-c(1, length(grid))]
      best <- grid[which.max(dens(grid))]
      expect_equal(got, best, tolerance = 2 * diff(grid[1:2]) / abs(best))
      # and it is not the image of the mode of theta, unless the Jacobian is
      # flat there
      expect_gt(dens(got), dens(maps[[key]]$inverse(xi)) * (1 - 1e-9))
    }
  }
})

test_that("a wide logit posterior has its mode at an end of the interval", {
  # theta ~ N(0, 3^2): the density of plogis(theta) is U-shaped, with its two
  # peaks where theta = +/- 9, and the search must not stop at the dip at 1/2
  got <- natural_summary(0, 3, 0, "logit", probs = 0.5)[4]
  expect_lt(min(got, 1 - got), plogis(-8.9))
})

test_that("the natural-scale density integrates to one", {
  for (key in names(maps)) {
    sn <- sn_cases[[1]]
    tf <- transform_dict[[key]]
    dens <- function(psi) {
      natural_density(
        psi,
        function(t) sn_density_ref(t, sn[["xi"]], sn[["omega"]], sn[["alpha"]]),
        tf
      )
    }
    lims <- tf$inverse(sn[["xi"]] + c(-9, 9) * sn[["omega"]])
    expect_equal(
      integrate(dens, lims[1], lims[2], rel.tol = 1e-10)$value,
      1,
      tolerance = 1e-8
    )
  }
})

# (c) End to end against closed forms -----------------------------------------

# Gamma-Poisson: y_i ~ Poisson(lambda), lambda ~ Gamma(shape, rate), so
# lambda | y ~ Gamma(S, R) with S = shape + sum(y) and R = rate + n
gp_y <- c(5, 7, 4, 6, 9, 3, 6, 8, 5, 7)
gp_S <- 2 + sum(gp_y)
gp_R <- 0.5 + length(gp_y)
gp_logpost <- function(theta) {
  lambda <- exp(theta[["log_lambda"]])
  sum(dpois(gp_y, lambda, log = TRUE)) +
    dgamma(lambda, shape = 2, rate = 0.5, log = TRUE) +
    theta[["log_lambda"]]
}
gp_fit <- inla_vanilla(gp_logpost, c(log_lambda = 0), silent = TRUE)

# Beta-binomial: y ~ Binomial(n, p), p ~ Beta(a, b), so p | y ~ Beta(A, B)
bb_A <- 1 + 6
bb_B <- 1 + 19
bb_logpost <- function(theta) {
  p <- plogis(theta[["logit_p"]])
  dbinom(6, 25, p, log = TRUE) +
    dbeta(p, 1, 1, log = TRUE) +
    log(p) +
    log1p(-p)
}
bb_fit <- inla_vanilla(bb_logpost, c(logit_p = 0), silent = TRUE)

test_that("lambda of the gamma-Poisson model matches its gamma posterior", {
  expect_identical(rownames(gp_fit$summary), "lambda")
  sd_exact <- sqrt(gp_S) / gp_R
  tab <- gp_fit$summary
  # In exact SDs of lambda. Observed: mean 1.7e-4, SD 0.05% too large,
  # quantiles 4.1e-3 at most, mode 3.6e-3
  expect_lt(abs(tab[, "Mean"] - gp_S / gp_R) / sd_exact, 1e-3)
  expect_equal(tab[, "SD"], sd_exact, tolerance = 2e-3, ignore_attr = TRUE)
  q_exact <- qgamma(c(0.025, 0.5, 0.975), gp_S, gp_R)
  expect_lt(max(abs(tab[, 3:5] - q_exact)) / sd_exact, 0.01)
  expect_lt(abs(tab[, "Mode"] - (gp_S - 1) / gp_R) / sd_exact, 0.01)
})

test_that("p of the beta-binomial model matches its beta posterior", {
  expect_identical(rownames(bb_fit$summary), "p")
  mean_exact <- bb_A / (bb_A + bb_B)
  sd_exact <- sqrt(mean_exact * (1 - mean_exact) / (bb_A + bb_B + 1))
  tab <- bb_fit$summary
  # In exact SDs of p. Observed: mean 1.1e-4, SD 0.02% too small, quantiles
  # 3.0e-3 at most, mode 6.2e-3
  expect_lt(abs(tab[, "Mean"] - mean_exact) / sd_exact, 1e-3)
  expect_equal(tab[, "SD"], sd_exact, tolerance = 2e-3, ignore_attr = TRUE)
  q_exact <- qbeta(c(0.025, 0.5, 0.975), bb_A, bb_B)
  expect_lt(max(abs(tab[, 3:5] - q_exact)) / sd_exact, 0.01)
  mode_exact <- (bb_A - 1) / (bb_A + bb_B - 2)
  expect_lt(abs(tab[, "Mode"] - mode_exact) / sd_exact, 0.02)
})

test_that("an exactly Gaussian log parameter gives the log-normal", {
  m0 <- 0.7
  s0 <- 0.45
  lp <- function(theta) dnorm(theta[["log_x"]], m0, s0, log = TRUE)
  # Observed relative errors: 8e-9 at most, the skew-normal fit's
  for (marg in c("skewnormal", "gaussian")) {
    fit <- inla_vanilla(lp, c(log_x = 0), marginals = marg, silent = TRUE)
    tab <- fit$summary["x", ]
    expect_equal(tab[["Mean"]], exp(m0 + s0^2 / 2), tolerance = 1e-6)
    expect_equal(
      tab[["SD"]],
      exp(m0 + s0^2 / 2) * sqrt(expm1(s0^2)),
      tolerance = 1e-6
    )
    expect_equal(
      unname(tab[3:5]),
      qlnorm(c(0.025, 0.5, 0.975), m0, s0),
      tolerance = 1e-6
    )
    expect_equal(tab[["Mode"]], exp(m0 - s0^2), tolerance = 1e-6)
  }
})

test_that("a correlation matches one-dimensional quadrature", {
  # Standardised pairs with known means and variances and a uniform prior on
  # the correlation, written on Fisher's z with its Jacobian
  set.seed(3)
  n <- 40
  u <- rnorm(n)
  v <- 0.6 * u + sqrt(1 - 0.6^2) * rnorm(n)
  loglik <- function(rho) {
    -n *
      log(2 * pi) -
      n / 2 * log1p(-rho^2) -
      (sum(u^2) - 2 * rho * sum(u * v) + sum(v^2)) / (2 * (1 - rho^2))
  }
  lp <- function(theta) {
    rho <- tanh(theta[["atanh_rho"]])
    loglik(rho) + log(0.5) + log1p(-rho^2)
  }
  fit <- inla_vanilla(lp, c(atanh_rho = 0), silent = TRUE)
  expect_identical(rownames(fit$summary), "rho")

  peak <- optimize(loglik, c(-0.999, 0.999), maximum = TRUE)$objective
  dens <- function(rho) exp(vapply(rho, loglik, numeric(1)) - peak)
  norm_const <- integrate(dens, -1, 1, rel.tol = 1e-10)$value
  m1 <- integrate(function(r) r * dens(r), -1, 1, rel.tol = 1e-10)$value /
    norm_const
  m2 <- integrate(function(r) r^2 * dens(r), -1, 1, rel.tol = 1e-10)$value /
    norm_const
  sd_exact <- sqrt(m2 - m1^2)
  cdf <- function(q) integrate(dens, -1, q, rel.tol = 1e-10)$value / norm_const
  tab <- fit$summary
  # Observed: mean 1.3e-3 exact SDs, SD 0.4% too small, and the probability
  # below each quantile within 1.4e-3 of its own
  expect_lt(abs(tab[, "Mean"] - m1) / sd_exact, 5e-3)
  expect_equal(tab[, "SD"], sd_exact, tolerance = 0.01, ignore_attr = TRUE)
  expect_lt(
    max(abs(vapply(tab[, 3:5], cdf, numeric(1)) - c(0.025, 0.5, 0.975))),
    5e-3
  )
  expect_true(all(abs(tab[, 3:5]) < 1))
})

# (d) The interface ------------------------------------------------------------

set.seed(1)
nm_y <- rnorm(30, 5, 2)
nm_logpost <- function(theta) {
  sigma <- exp(theta[["log_sigma"]])
  sum(dnorm(nm_y, theta[["mu"]], sigma, log = TRUE)) +
    dnorm(theta[["mu"]], 0, 10, log = TRUE) +
    dexp(sigma, log = TRUE) +
    theta[["log_sigma"]]
}
nm_init <- c(mu = 0, log_sigma = 0)
nm_fit <- inla_vanilla(nm_logpost, nm_init, silent = TRUE)
nm_fit_off <- inla_vanilla(
  nm_logpost,
  nm_init,
  transform = FALSE,
  silent = TRUE
)

test_that("the transform only changes what is reported", {
  for (part in c("sn", "mode", "cov", "mean_vb", "logml", "counts", "scan")) {
    expect_identical(nm_fit[[part]], nm_fit_off[[part]])
  }
  expect_identical(nm_fit$diagnostics, nm_fit_off$diagnostics)
  expect_identical(rownames(nm_fit$sn), c("mu", "log_sigma"))
  expect_identical(diagnostics(nm_fit, "param")$names, c("mu", "log_sigma"))

  expect_identical(nm_fit$transform$transform, c("identity", "log"))
  expect_identical(nm_fit$transform$name, c("mu", "sigma"))
  expect_true(all(nm_fit_off$transform$transform == "identity"))
})

test_that("the summaries are on the natural scale, and can be taken off it", {
  expect_identical(rownames(nm_fit$summary), c("mu", "sigma"))
  expect_identical(rownames(nm_fit_off$summary), c("mu", "log_sigma"))
  # an untransformed parameter is unchanged
  expect_identical(nm_fit$summary["mu", ], nm_fit_off$summary["mu", ])
  # the quantiles of sigma are exp() of those of log_sigma
  expect_equal(
    nm_fit$summary["sigma", 3:5],
    exp(nm_fit_off$summary["log_sigma", 3:5])
  )
  # and the mean is above exp() of the mean, by Jensen's inequality
  expect_gt(
    nm_fit$summary["sigma", "Mean"],
    exp(nm_fit_off$summary["log_sigma", "Mean"])
  )

  s <- summary(nm_fit, probs = c(0.1, 0.9))
  expect_identical(rownames(s$table), c("mu", "sigma"))
  expect_identical(colnames(s$table), c("Mean", "SD", "10%", "90%", "Mode"))
  expect_identical(rownames(s$diagnostics), c("mu", "log_sigma"))
  expect_identical(
    summary(nm_fit, transform = FALSE)$table,
    summary(nm_fit_off)$table
  )
  expect_null(summary(nm_fit, transform = FALSE)$transform)

  expect_identical(names(coef(nm_fit)), c("mu", "sigma"))
  expect_identical(coef(nm_fit)[["sigma"]], nm_fit$summary[["sigma", "Mean"]])
  expect_identical(names(coef(nm_fit_off)), c("mu", "log_sigma"))
})

test_that("print() and summary() say which rows are on the natural scale", {
  note <- "Natural scale: sigma = exp\\(log_sigma\\)\\."
  expect_output(print(nm_fit), note)
  expect_output(print(summary(nm_fit)), note)
  expect_no_match(capture.output(print(nm_fit_off)), "Natural scale")
  expect_no_match(
    capture.output(print(summary(nm_fit, transform = FALSE))),
    "Natural scale"
  )
  # no wider than the other output
  expect_true(all(nchar(capture.output(print(nm_fit))) <= 80))
})

test_that("many transformed parameters are counted, not listed", {
  nms <- c(paste0("log_s", 1:4), "atanh_r", "b")
  lp <- function(theta) sum(dnorm(theta, 0, 0.3, log = TRUE))
  fit <- inla_vanilla(
    lp,
    stats::setNames(numeric(6), nms),
    marginals = "gaussian",
    vb = "off",
    silent = TRUE
  )
  out <- capture.output(print(fit))
  expect_match(
    out[length(out)],
    "Natural scale: exp\\(\\) of 4 log_ parameters, tanh\\(\\) of 1 atanh_ parameter\\.$"
  )
  expect_true(all(nchar(out) <= 80))
  # four or fewer are listed in full
  fit <- inla_vanilla(
    lp,
    stats::setNames(numeric(5), nms[-1]),
    marginals = "gaussian",
    vb = "off",
    silent = TRUE
  )
  note <- gsub(" +", " ", paste(capture.output(print(fit)), collapse = " "))
  expect_match(note, "s4 = exp\\(log_s4\\), r = tanh\\(atanh_r\\)\\.")
})

test_that("sampling() returns natural-scale draws, and fun still gets theta", {
  set.seed(8)
  natural <- sampling(nm_fit, nsamp = 200)
  set.seed(8)
  working <- sampling(nm_fit, nsamp = 200, transform = FALSE)
  expect_identical(colnames(natural), c("mu", "sigma"))
  expect_identical(colnames(working), c("mu", "log_sigma"))
  expect_identical(natural[, "mu"], working[, "mu"])
  expect_identical(natural[, "sigma"], exp(working[, "log_sigma"]))
  set.seed(8)
  expect_identical(sampling(nm_fit_off, nsamp = 200), working)

  # fun receives what logpost receives, whatever `transform`
  set.seed(8)
  seen <- sampling(nm_fit, nsamp = 200, fun = function(theta) theta)
  expect_identical(seen, working)
  set.seed(8)
  cv <- sampling(
    nm_fit,
    nsamp = 200,
    fun = function(theta) exp(theta[["log_sigma"]]) / theta[["mu"]]
  )
  expect_equal(cv[, 1], natural[, "sigma"] / natural[, "mu"])
})

test_that("the draws agree with the natural-scale summary", {
  set.seed(21)
  draws <- sampling(nm_fit, nsamp = 4e4)
  tab <- nm_fit$summary
  # Monte Carlo: the SE of a mean is SD / 200, and four of them are allowed
  # (observed 1.3). The SDs were within 0.7% and the quantiles of sigma
  # within 0.1%.
  expect_lt(
    max(abs(colMeans(draws) - tab[, "Mean"]) / tab[, "SD"]),
    4 / 200
  )
  expect_equal(apply(draws, 2, sd), tab[, "SD"], tolerance = 0.03)
  expect_equal(
    unname(quantile(draws[, "sigma"], c(0.025, 0.5, 0.975))),
    unname(tab["sigma", 3:5]),
    tolerance = 0.01
  )
})

test_that("plot() draws the natural scale and takes either name", {
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit({
    grDevices::dev.off()
    unlink(path)
  })
  expect_no_error(plot(nm_fit))
  expect_no_error(plot(nm_fit, pars = "sigma"))
  expect_no_error(plot(nm_fit, pars = "log_sigma"))
  expect_no_error(plot(nm_fit, pars = 2, transform = FALSE))
  expect_no_error(plot(nm_fit, type = "sn_fit"))
  expect_no_error(plot(bb_fit))
  expect_error(plot(nm_fit, pars = "tau"), "Unknown parameter")
  expect_error(plot(nm_fit, transform = NA), "transform")
})

test_that("a clash of names warns once, and the fit records it", {
  lp <- function(theta) sum(dnorm(theta, 0, 1, log = TRUE))
  init <- c(log_sigma = 0, sigma = 0)
  expect_warning(
    fit <- inla_vanilla(lp, init, marginals = "gaussian", vb = "off"),
    "unconstrained scale"
  ) |>
    suppressMessages()
  expect_identical(rownames(fit$summary), c("log_sigma", "sigma"))
  expect_length(fit$warnings, 1L)
  expect_match(fit$warnings, "log_sigma")

  fit <- expect_silent(inla_vanilla(
    lp,
    init,
    marginals = "gaussian",
    vb = "off",
    silent = TRUE
  ))
  expect_length(fit$warnings, 1L)
  # switched off, there is nothing to clash
  fit <- expect_silent(inla_vanilla(
    lp,
    init,
    marginals = "gaussian",
    vb = "off",
    transform = FALSE,
    silent = TRUE
  ))
  expect_length(fit$warnings, 0L)
})

test_that("`transform` must be TRUE or FALSE", {
  expect_error(
    inla_vanilla(nm_logpost, nm_init, transform = "log"),
    "transform",
    class = "inlavanilla_error"
  )
  expect_error(summary(nm_fit, transform = NA), class = "inlavanilla_error")
  expect_error(
    sampling(nm_fit, nsamp = 2, transform = 1),
    class = "inlavanilla_error"
  )
})

test_that("a fit saved before transforms existed still works", {
  old <- nm_fit_off
  old$transform <- NULL
  expect_identical(summary(old)$table, summary(nm_fit_off)$table)
  expect_output(print(old), "log_sigma")
  set.seed(4)
  draws <- sampling(old, nsamp = 5)
  expect_identical(colnames(draws), c("mu", "log_sigma"))
})
