# Tests of the user interface of inla_vanilla() and of its S3 methods: input
# validation, the check of an analytic gradient, parameter names, arguments
# passed through `...`, log-posteriors that are not finite away from the mode,
# improper posteriors, progress messages, and the print(), summary() and
# plot() methods.
#
# Ground truth, where a test needs one, is closed form: Gaussian posteriors,
# where the skew-normal marginals should be exact, a gamma posterior written
# on its natural scale, the gamma-Poisson model on the log scale, and a
# truncated normal. Tolerances on the approximate (non-Gaussian) cases were
# calibrated by printing the actual errors first; each comment states the
# observed value.

# Models ----------------------------------------------------------------------

# A Gaussian log-posterior with every mean equal to 1 and unit variances. The
# mode is kept away from zero: nlminb() without a gradient can report "false
# convergence" when it lands exactly on a mode whose differenced gradient is
# exactly zero (as at the origin of a centred quadratic).
lp_gauss <- function(theta) -0.5 * sum((theta - 1)^2)
gr_gauss <- function(theta) -(theta - 1)

# y_i ~ Poisson(lambda), lambda ~ Gamma(a0, b0), on theta = log(lambda), with
# the Jacobian. Exactly lambda | y ~ Gamma(a0 + sum(y), b0 + n), so the mode of
# theta is log(A / B), its mean digamma(A) - log(B) and its SD sqrt(trigamma(A)).
# The data reach logpost and grad through `...`.
pois_y <- c(2, 0, 3, 1, 4, 2, 1)
pois_a0 <- 2
pois_b0 <- 1
pois_A <- pois_a0 + sum(pois_y)
pois_B <- pois_b0 + length(pois_y)

lp_pois <- function(theta, y, a0, b0) {
  lambda <- exp(theta[[1]])
  sum(dpois(y, lambda, log = TRUE)) +
    dgamma(lambda, a0, b0, log = TRUE) +
    theta[[1]]
}
gr_pois <- function(theta, y, a0, b0) {
  lambda <- exp(theta[[1]])
  sum(y) - length(y) * lambda + a0 - b0 * lambda
}

# Two independent, skewed parameters: the Poisson rate above (log scale) and a
# binomial probability with a uniform prior (logit scale). Data are fixed, so
# the fit is deterministic. Used for the printed output and the plots.
lp_two <- function(theta) {
  lambda <- exp(theta[["log_rate"]])
  p <- plogis(theta[["logit_p"]])
  sum(dpois(pois_y, lambda, log = TRUE)) +
    dgamma(lambda, 2, 1, log = TRUE) +
    theta[["log_rate"]] +
    dbinom(3, 12, p, log = TRUE) +
    log(p) +
    log1p(-p)
}
init_two <- c(log_rate = 0.1, logit_p = 0.2)

fit_two <- inla_vanilla(lp_two, init_two, verbose = FALSE)

# Helpers ---------------------------------------------------------------------

# Make printed output stable between runs. Timings (12ms, 1.2s) vary, and the
# rule above print(fit) is padded with dashes to a fixed width, so the number
# of dashes varies with the width of the time. The numbers of calls are
# deterministic on one machine, but they follow the optimiser's path, which
# floating-point differences between platforms could change; the gradient at
# the mode is at the level of finite-difference noise.
scrub_output <- function(x) {
  x <- gsub("[0-9.]+(ms|s)\\b", "<time>", x)
  x <- gsub("-+ <time> --", "-- <time> --", x)
  x <- gsub("[0-9]+ (logpost|grad)\\b", "<n> \\1", x)
  x <- gsub("gradient [-+0-9.e]+ SDs", "gradient <g> SDs", x)
  x
}

# Expect an error raised by the package itself with cli::cli_abort(), rather
# than one from base R deep inside the code. Unlike expect_error(class = ), a
# base R error is reported as a failure, so the rest of the test still runs.
expect_package_error <- function(object) {
  err <- tryCatch(
    {
      object
      NULL
    },
    error = identity
  )
  got <- if (is.null(err)) {
    "no error"
  } else {
    paste0("<", class(err)[1], "> ", conditionMessage(err))
  }
  expect(
    inherits(err, "rlang_error"),
    paste0("Expected an error from the package (<rlang_error>), got ", got)
  )
  invisible(err)
}

# Graphical parameters that every high-level plot sets, whatever the method
# does: the user coordinates and the axis tick positions.
par_set_by_any_plot <- c("usr", "xaxp", "yaxp")

# (a) Input validation ---------------------------------------------------------

test_that("logpost and grad must be functions", {
  expect_error(
    inla_vanilla(1, c(a = 0), verbose = FALSE),
    "`logpost` must be a function",
    class = "rlang_error"
  )
  expect_error(
    inla_vanilla("lp_gauss", c(a = 0), verbose = FALSE),
    "`logpost` must be a function",
    class = "rlang_error"
  )
  expect_error(
    inla_vanilla(lp_gauss, c(a = 0), grad = "gr", verbose = FALSE),
    "`grad` must be a function",
    class = "rlang_error"
  )
  expect_error(
    inla_vanilla(lp_gauss, c(a = 0), grad = 1, verbose = FALSE),
    "`grad` must be a function",
    class = "rlang_error"
  )
})

test_that("init must be a non-empty, finite numeric vector", {
  for (bad in list("a", numeric(0), NULL, list(a = 1), NA, TRUE)) {
    expect_error(
      inla_vanilla(lp_gauss, bad, verbose = FALSE),
      "`init` must be a non-empty numeric vector",
      class = "rlang_error"
    )
  }
  for (bad in list(c(a = NA_real_), c(a = 0, b = NaN), c(a = Inf), -Inf)) {
    expect_error(
      inla_vanilla(lp_gauss, bad, verbose = FALSE),
      "`init` must contain only finite values",
      class = "rlang_error"
    )
  }
})

test_that("init must be fully named, or not named at all, without duplicates", {
  expect_error(
    inla_vanilla(lp_gauss, c(a = 0, 1), verbose = FALSE),
    "fully named, or not named at all",
    class = "rlang_error"
  )
  expect_error(
    inla_vanilla(
      lp_gauss,
      stats::setNames(c(0, 1), c("a", NA)),
      verbose = FALSE
    ),
    "fully named, or not named at all",
    class = "rlang_error"
  )
  expect_error(
    inla_vanilla(lp_gauss, c(a = 0, b = 1, a = 2), verbose = FALSE),
    "duplicated names",
    class = "rlang_error"
  )
})

test_that("logpost must return a single finite number at init", {
  not_single <- list(
    vector = function(theta) theta,
    string = function(theta) "-1",
    logical_na = function(theta) NA,
    empty = function(theta) numeric(0)
  )
  for (lp in not_single) {
    expect_error(
      inla_vanilla(lp, c(a = 0.5, b = 0.5), verbose = FALSE),
      "`logpost` must return a single number",
      class = "rlang_error"
    )
  }
  not_finite <- list(
    na = function(theta) NA_real_,
    nan = function(theta) NaN,
    minus_inf = function(theta) -Inf,
    plus_inf = function(theta) Inf
  )
  for (lp in not_finite) {
    expect_error(
      inla_vanilla(lp, c(a = 0.5), verbose = FALSE),
      "`logpost` is not finite at `init`",
      class = "rlang_error"
    )
  }
  # Finite near the mode but not at init: the check is about init itself
  lp_hole <- function(theta) if (theta[[1]] > 5) -Inf else lp_gauss(theta)
  expect_error(
    inla_vanilla(lp_hole, c(a = 6), verbose = FALSE),
    "not finite at `init`",
    class = "rlang_error"
  )
})

test_that("grad must return a finite vector of the right length at init", {
  init <- c(a = 0.5, b = 0.5)
  for (gr in list(
    function(theta) 1,
    function(theta) c(1, 2, 3),
    function(theta) "a",
    function(theta) NULL
  )) {
    expect_error(
      inla_vanilla(lp_gauss, init, grad = gr, verbose = FALSE),
      "`grad` must return a numeric vector of length 2",
      class = "rlang_error"
    )
  }
  for (gr in list(
    function(theta) c(NA, 1),
    function(theta) c(NaN, 1),
    function(theta) c(1, Inf)
  )) {
    expect_error(
      inla_vanilla(lp_gauss, init, grad = gr, verbose = FALSE),
      "`grad` is not finite at `init`",
      class = "rlang_error"
    )
  }
})

test_that("vb and verbose must be a single TRUE or FALSE", {
  for (bad in list(NA, "yes", 1, c(TRUE, FALSE), logical(0), NULL)) {
    expect_error(
      inla_vanilla(lp_gauss, c(a = 0.5), vb = bad, verbose = FALSE),
      "`vb` must be `TRUE` or `FALSE`",
      class = "rlang_error"
    )
    expect_error(
      inla_vanilla(lp_gauss, c(a = 0.5), verbose = bad),
      "`verbose` must be `TRUE` or `FALSE`",
      class = "rlang_error"
    )
  }
})

test_that("control must be a list of known entries", {
  for (bad in list(1, "n_grid", c(n_grid = 21))) {
    expect_error(
      inla_vanilla(lp_gauss, c(a = 0.5), control = bad, verbose = FALSE),
      "`control` must be a list",
      class = "rlang_error"
    )
  }
  expect_error(
    inla_vanilla(lp_gauss, c(a = 0.5), control = list(foo = 1), verbose = FALSE),
    "Unknown `control` entry: \"foo\"",
    class = "rlang_error"
  )
  expect_error(
    inla_vanilla(
      lp_gauss,
      c(a = 0.5),
      control = list(n_grid = 21, foo = 1, ngrid = 2),
      verbose = FALSE
    ),
    "Unknown `control` entries: \"foo\" and \"ngrid\"",
    class = "rlang_error"
  )
  # Unnamed entries cannot be matched to anything
  expect_error(
    inla_vanilla(lp_gauss, c(a = 0.5), control = list(21), verbose = FALSE),
    "Unknown `control`",
    class = "rlang_error"
  )
  expect_error(
    inla_vanilla(
      lp_gauss,
      c(a = 0.5),
      control = list(n_grid = 21, 1),
      verbose = FALSE
    ),
    "Unknown `control`",
    class = "rlang_error"
  )
})

test_that("control$n_grid must be a single number of at least 7", {
  for (bad in list(6, 0, -21, "21", c(21, 31), integer(0))) {
    expect_error(
      inla_vanilla(
        lp_gauss,
        c(a = 0.5),
        control = list(n_grid = bad),
        verbose = FALSE
      ),
      "`control$n_grid` must be a single number of at least 7",
      fixed = TRUE,
      class = "rlang_error"
    )
  }
})

test_that("a missing control$n_grid gets the package's own error", {
  # Observed: NA fails inside the check with base R's "missing value where
  # TRUE/FALSE needed".
  expect_package_error(
    inla_vanilla(
      lp_gauss,
      c(a = 0.5),
      control = list(n_grid = NA_real_),
      verbose = FALSE
    )
  )
})

test_that("an infinite control$n_grid gets the package's own error", {
  # Observed: Inf passes the check, then fails in seq() with base R's
  # "'length.out' must be a non-negative number".
  expect_package_error(
    inla_vanilla(
      lp_gauss,
      c(a = 0.5),
      control = list(n_grid = Inf),
      verbose = FALSE
    )
  )
})

test_that("control$log_threshold must be a single negative number", {
  for (bad in list(0, 1, "-6", c(-6, -8))) {
    expect_error(
      inla_vanilla(
        lp_gauss,
        c(a = 0.5),
        control = list(log_threshold = bad),
        verbose = FALSE
      ),
      "`control$log_threshold` must be a single negative number",
      fixed = TRUE,
      class = "rlang_error"
    )
  }
})

test_that("a missing control$log_threshold gets the package's own error", {
  # Observed: NA fails inside the check with base R's "missing value where
  # TRUE/FALSE needed".
  expect_package_error(
    inla_vanilla(
      lp_gauss,
      c(a = 0.5),
      control = list(log_threshold = NA_real_),
      verbose = FALSE
    )
  )
})

test_that("valid control entries are accepted and keep a Gaussian exact", {
  ctrl <- list(
    n_grid = 31,
    log_threshold = -8,
    vb_maxit = 50,
    vb_tol = 1e-4,
    nlminb = list(rel.tol = 1e-12)
  )
  fit <- inla_vanilla(
    lp_gauss,
    c(a = 0.3, b = 1.4),
    control = ctrl,
    verbose = FALSE
  )
  expect_equal(unname(fit$summary[, "Mean"]), c(1, 1), tolerance = 1e-6)
  expect_equal(unname(fit$summary[, "SD"]), c(1, 1), tolerance = 1e-6)
})

test_that("any accepted n_grid gives exact marginals for a Gaussian posterior", {
  # N(3, 2^2): a skew-normal with alpha = 0 is exact, so the summary should
  # be exact for every grid. Odd grids are (to ~1e-8). With an even n_grid
  # the grid has no point at the mode, yet fit_marginals() gives the point
  # nearest to it the height of the mode: observed for n_grid = 8, a mean
  # off by 0.074 (0.037 SDs), an SD 2.1% too small and a mode off by 0.19
  # (0.095 SDs); for n_grid = 20, an SD 0.18% too small.
  lp <- function(theta) dnorm(theta[[1]], 3, 2, log = TRUE)
  for (n_grid in c(7, 8, 20, 21)) {
    fit <- inla_vanilla(
      lp,
      c(a = 0.5),
      control = list(n_grid = n_grid),
      verbose = FALSE
    )
    expect_equal(
      unname(fit$summary[1, c("Mean", "SD", "Mode")]),
      c(3, 2, 3),
      tolerance = 1e-6,
      label = paste0("Mean, SD and Mode with n_grid = ", n_grid)
    )
  }
})

test_that("validation errors are informative", {
  expect_snapshot(error = TRUE, {
    inla_vanilla(function(theta) 0, c(a = 0, b = 1, a = 2))
  })
  expect_snapshot(error = TRUE, {
    inla_vanilla("lp", c(a = 0))
  })
  expect_snapshot(error = TRUE, {
    inla_vanilla(function(theta) theta, c(a = 0.5, b = 0.5))
  })
  expect_snapshot(error = TRUE, {
    inla_vanilla(function(theta) -Inf, c(a = 0.5))
  })
  expect_snapshot(error = TRUE, {
    inla_vanilla(lp_gauss, c(a = 0.5, b = 0.5), grad = function(theta) 1)
  })
  expect_snapshot(error = TRUE, {
    inla_vanilla(lp_gauss, c(a = 0.5), control = list(foo = 1, ngrid = 21))
  })
})

test_that("errors are reported as coming from inla_vanilla()", {
  # The user called inla_vanilla(), so the error should name it, as r-lib
  # packages do by passing `call = caller_env()` through their checkers.
  # Observed: the snapshots above show the internal helper instead, e.g.
  # "Error in `check_init()`" and "Error in `vanilla_control()`".
  errs <- list(
    function_check = function() inla_vanilla(1, c(a = 0)),
    init_check = function() inla_vanilla(lp_gauss, "a"),
    logpost_check = function() inla_vanilla(function(theta) NaN, c(a = 0.5)),
    grad_check = function() {
      inla_vanilla(lp_gauss, c(a = 0.5), grad = function(theta) NA)
    },
    control_check = function() {
      inla_vanilla(lp_gauss, c(a = 0.5), control = list(foo = 1))
    }
  )
  called <- vapply(
    errs,
    function(f) {
      err <- tryCatch(f(), error = identity)
      if (inherits(err, "rlang_error")) deparse(err$call[[1]]) else NA
    },
    ""
  )
  expect_identical(
    called,
    c(
      function_check = "inla_vanilla",
      init_check = "inla_vanilla",
      logpost_check = "inla_vanilla",
      grad_check = "inla_vanilla",
      control_check = "inla_vanilla"
    )
  )
})

# (b) The gradient check ------------------------------------------------------

test_that("a wrong analytic gradient triggers the mismatch warning", {
  # Right for `a`, twice too large for `b`: the warning names `b`. (A sign
  # error would also be caught, but then the optimiser climbs away from the
  # mode and the fit fails later with the Hessian error.)
  gr_wrong <- function(theta) c(-(theta[[1]] - 1), -2 * (theta[[2]] - 1))
  expect_warning(
    inla_vanilla(
      lp_gauss,
      c(a = 0.3, b = 0.2),
      grad = gr_wrong,
      verbose = FALSE
    ),
    "does not match finite differences",
    class = "rlang_warning"
  )
  expect_snapshot(
    invisible(inla_vanilla(
      lp_gauss,
      c(a = 0.3, b = 0.2),
      grad = gr_wrong,
      verbose = FALSE
    ))
  )
  # A gradient off by a constant factor is also caught
  expect_warning(
    inla_vanilla(
      lp_pois,
      0.3,
      grad = function(theta, ...) 2 * gr_pois(theta, ...),
      y = pois_y,
      a0 = pois_a0,
      b0 = pois_b0,
      verbose = FALSE
    ),
    "does not match finite differences",
    class = "rlang_warning"
  )
})

test_that("a correct analytic gradient passes the check silently", {
  expect_silent(
    inla_vanilla(
      lp_gauss,
      c(a = 0.3, b = 0.2),
      grad = gr_gauss,
      verbose = FALSE
    )
  )
  expect_silent(
    inla_vanilla(
      lp_pois,
      0.3,
      grad = gr_pois,
      y = pois_y,
      a0 = pois_a0,
      b0 = pois_b0,
      verbose = FALSE
    )
  )
  # A 1 x m matrix, as TMB's obj$gr() returns, is accepted
  fit <- inla_vanilla(
    lp_gauss,
    c(a = 0.3, b = 0.2),
    grad = function(theta) matrix(gr_gauss(theta), nrow = 1),
    verbose = FALSE
  )
  expect_equal(unname(fit$mode), c(1, 1), tolerance = 1e-8)
})

# (c) Names, `...`, non-finite regions, one parameter -------------------------

test_that("unnamed init gets the names theta[1], theta[2], ...", {
  seen_lp <- NULL
  seen_gr <- NULL
  lp <- function(theta) {
    seen_lp <<- union(seen_lp, list(names(theta)))
    lp_gauss(theta)
  }
  gr <- function(theta) {
    seen_gr <<- union(seen_gr, list(names(theta)))
    gr_gauss(theta)
  }
  fit <- inla_vanilla(lp, c(0.1, 0.2, 0.3), grad = gr, verbose = FALSE)

  nms <- c("theta[1]", "theta[2]", "theta[3]")
  # Every call of logpost and grad saw the same names
  expect_identical(seen_lp, list(nms))
  expect_identical(seen_gr, list(nms))
  expect_identical(names(fit$mode), nms)
  expect_identical(names(fit$mean_vb), nms)
  expect_identical(rownames(fit$sn), nms)
  expect_identical(colnames(fit$sn), c("xi", "omega", "alpha"))
  expect_identical(dimnames(fit$cov), list(nms, nms))
  expect_identical(rownames(fit$summary), nms)
  expect_identical(
    colnames(fit$summary),
    c("Mean", "SD", "2.5%", "50%", "97.5%", "Mode")
  )
  expect_identical(rownames(summary(fit)$diagnostics), nms)

  set.seed(1)
  draws <- sampling(fit, nsamp = 5)
  expect_identical(colnames(draws), nms)
  seen_fun <- NULL
  set.seed(1)
  out <- sampling(fit, nsamp = 5, fun = function(theta) {
    seen_fun <<- union(seen_fun, list(names(theta)))
    theta[["theta[2]"]]
  })
  expect_identical(seen_fun, list(nms))
  expect_equal(out[, 1], draws[, "theta[2]"])
})

test_that("named init keeps its names, and integer init is accepted", {
  fit <- inla_vanilla(lp_gauss, c(alpha = 0L, beta = 2L), verbose = FALSE)
  expect_identical(rownames(fit$summary), c("alpha", "beta"))
  expect_equal(unname(fit$mode), c(1, 1), tolerance = 1e-6)
})

test_that("arguments in ... reach logpost and grad", {
  # The mode of log(lambda) is log(A / B) exactly; it depends on all of y, a0
  # and b0, so the fit can only get it right if they arrive.
  fit_fd <- inla_vanilla(
    lp_pois,
    0.3,
    y = pois_y,
    a0 = pois_a0,
    b0 = pois_b0,
    verbose = FALSE
  )
  fit_gr <- inla_vanilla(
    lp_pois,
    0.3,
    grad = gr_pois,
    y = pois_y,
    a0 = pois_a0,
    b0 = pois_b0,
    verbose = FALSE
  )
  exact_mode <- log(pois_A / pois_B)
  expect_equal(fit_fd$mode[[1]], exact_mode, tolerance = 1e-6)
  expect_equal(fit_gr$mode[[1]], exact_mode, tolerance = 1e-7)
  expect_gt(fit_gr$counts[["grad"]], 0)
  expect_identical(fit_fd$counts[["grad"]], 0L)

  # Different data through ... give a different, still exact, mode
  y2 <- c(10, 12, 9)
  fit2 <- inla_vanilla(
    lp_pois,
    0.3,
    grad = gr_pois,
    y = y2,
    a0 = pois_a0,
    b0 = pois_b0,
    verbose = FALSE
  )
  expect_equal(
    fit2$mode[[1]],
    log((pois_a0 + sum(y2)) / (pois_b0 + length(y2))),
    tolerance = 1e-7
  )

  # Without them logpost fails
  expect_error(inla_vanilla(lp_pois, 0.3, verbose = FALSE), "missing")
})

test_that("NaN outside the support, away from the mode, does not break the fit", {
  # x ~ Gamma(10, 2) on its natural scale, written so that logpost is NaN for
  # x <= 0 (and -Inf in a second version). The Laplace SD is 1.5, so the
  # boundary is 3 SDs below the mode 4.5 and the bottom of the scan, down to
  # 4 SDs, falls outside the support.
  shape <- 10
  rate <- 2
  n_outside <- 0
  lp_nan <- function(theta) {
    x <- theta[[1]]
    if (x <= 0) {
      n_outside <<- n_outside + 1
      return(NaN)
    }
    (shape - 1) * log(x) - rate * x
  }
  lp_inf <- function(theta) {
    x <- theta[[1]]
    if (x <= 0) -Inf else (shape - 1) * log(x) - rate * x
  }

  expect_silent(fit <- inla_vanilla(lp_nan, c(x = 4), verbose = FALSE))
  expect_gt(n_outside, 0)
  expect_true(all(is.finite(fit$summary)))
  expect_true(all(is.finite(fit$sn)))
  expect_true(is.finite(fit$logml))

  exact_mean <- shape / rate
  exact_sd <- sqrt(shape) / rate
  exact_q <- qgamma(c(0.025, 0.5, 0.975), shape, rate)
  # Observed: mean 0.0055 SDs off, SD 0.22% small, quantiles at most 0.016
  # SDs off. The Laplace mean (the mode) is 0.32 SDs off.
  err_mean <- abs(fit$summary[1, "Mean"] - exact_mean) / exact_sd
  expect_lt(err_mean, 0.02)
  expect_lt(abs(fit$summary[1, "SD"] / exact_sd - 1), 0.01)
  expect_lt(
    max(abs(fit$summary[1, c("2.5%", "50%", "97.5%")] - exact_q)) / exact_sd,
    0.04
  )
  expect_lt(err_mean, abs(fit$mode[[1]] - exact_mean) / exact_sd)

  # Starting near the boundary, where the optimiser's first steps can leave
  # the support, finds the same mode
  fit_near <- inla_vanilla(lp_nan, c(x = 0.1), verbose = FALSE)
  expect_equal(fit_near$mode[[1]], (shape - 1) / rate, tolerance = 1e-5)

  # -Inf and NaN are treated alike
  fit_inf <- inla_vanilla(lp_inf, c(x = 4), verbose = FALSE)
  expect_equal(fit_inf$summary, fit$summary)
})

test_that("-Inf inside the scan range of one parameter does not break the fit", {
  # a ~ N(1, 1) truncated to a > -1.5 (2.5 SDs below the mode), and
  # b | a ~ N(a / 2, 0.75). The scan of `a` reaches 4 SDs, so its lowest grid
  # points are outside the region, where the Gaussian log-height would still
  # be above the threshold of the fit.
  n_outside <- 0
  lp <- function(theta) {
    a <- theta[["a"]]
    if (a < -1.5) {
      n_outside <<- n_outside + 1
      return(-Inf)
    }
    dnorm(a, 1, 1, log = TRUE) +
      dnorm(theta[["b"]], 0.5 * a, sqrt(0.75), log = TRUE)
  }
  expect_silent(fit <- inla_vanilla(lp, c(a = 0.3, b = -0.2), verbose = FALSE))
  expect_gt(n_outside, 0)
  expect_true(all(is.finite(fit$summary)))

  # Exact truncated-normal moments of `a`
  lo <- -2.5 # truncation point in standard units
  z_mass <- 1 - pnorm(lo)
  lambda <- dnorm(lo) / z_mass
  exact_mean <- 1 + lambda
  exact_sd <- sqrt(1 + lo * lambda - lambda^2)
  # Observed: mean 0.018 SDs off and SD 2.3% too large; the skew-normal sees
  # the Gaussian inside the region and cannot represent the hard edge.
  expect_lt(abs(fit$summary["a", "Mean"] - exact_mean) / exact_sd, 0.05)
  expect_lt(abs(fit$summary["a", "SD"] / exact_sd - 1), 0.05)
})

test_that("a one-parameter model works end to end", {
  fit <- inla_vanilla(
    lp_pois,
    0.3,
    grad = gr_pois,
    y = pois_y,
    a0 = pois_a0,
    b0 = pois_b0,
    verbose = FALSE
  )
  expect_s3_class(fit, "inla_vanilla")
  expect_identical(dim(fit$summary), c(1L, 6L))
  expect_identical(dim(fit$sn), c(1L, 3L))
  expect_identical(dim(fit$cov), c(1L, 1L))
  expect_identical(rownames(fit$summary), "theta[1]")

  # Against the exact log-gamma posterior. Observed: mean 0.00066 SDs off,
  # SD 0.011% off.
  exact_mean <- digamma(pois_A) - log(pois_B)
  exact_sd <- sqrt(trigamma(pois_A))
  expect_lt(abs(fit$summary[1, "Mean"] - exact_mean) / exact_sd, 0.005)
  expect_lt(abs(fit$summary[1, "SD"] / exact_sd - 1), 0.002)

  out <- capture_output(print(fit))
  expect_match(out, "1 parameter[^s]")
  out <- capture_output(print(summary(fit)))
  expect_match(out, "theta[1]", fixed = TRUE)

  grDevices::pdf(NULL)
  dev <- grDevices::dev.cur()
  on.exit(grDevices::dev.off(dev), add = TRUE)
  expect_no_error(plot(fit))
  expect_no_error(plot(fit, laplace = FALSE))

  set.seed(1)
  draws <- sampling(fit, nsamp = 2000)
  expect_identical(dim(draws), c(2000L, 1L))
  expect_identical(colnames(draws), "theta[1]")
  expect_true(all(is.finite(draws)))
  # Within 4 Monte Carlo standard errors of the skew-normal mean
  expect_lt(
    abs(mean(draws) - fit$summary[1, "Mean"]),
    4 * fit$summary[1, "SD"] / sqrt(2000)
  )
  set.seed(1)
  lambda <- sampling(fit, nsamp = 2000, fun = function(theta) {
    c(lambda = exp(theta[[1]]))
  })
  expect_identical(dim(lambda), c(2000L, 1L))
  expect_identical(colnames(lambda), "lambda")
  expect_equal(lambda[, 1], exp(draws[, 1]))
})

# (d) Improper or non-identified posteriors ----------------------------------

test_that("a posterior flat in one direction fails with the Hessian error", {
  lp_flat <- function(theta) -0.5 * (theta[["a"]] - 1)^2
  gr_flat <- function(theta) c(-(theta[["a"]] - 1), 0)
  expect_error(
    inla_vanilla(lp_flat, c(a = 0.3, b = 0.2), verbose = FALSE),
    "Hessian at the posterior mode is not positive definite",
    class = "rlang_error"
  )
  expect_error(
    inla_vanilla(lp_flat, c(a = 0.3, b = 0.2), grad = gr_flat, verbose = FALSE),
    "Hessian at the posterior mode is not positive definite",
    class = "rlang_error"
  )
})

test_that("a non-identified posterior fails with the Hessian error", {
  # Only a + b is identified
  lp_sum <- function(theta) -0.5 * (theta[["a"]] + theta[["b"]] - 1)^2
  gr_sum <- function(theta) rep(-(theta[["a"]] + theta[["b"]] - 1), 2)
  expect_error(
    inla_vanilla(lp_sum, c(a = 0.3, b = 0.2), verbose = FALSE),
    "Hessian at the posterior mode is not positive definite",
    class = "rlang_error"
  )
  expect_error(
    inla_vanilla(lp_sum, c(a = 0.3, b = 0.2), grad = gr_sum, verbose = FALSE),
    "Hessian at the posterior mode is not positive definite",
    class = "rlang_error"
  )
})

test_that("an improper posterior with no mode fails with the Hessian error", {
  # Increasing without bound in b, and flat (logpost = 0) in one dimension
  lp_up <- function(theta) -0.5 * (theta[["a"]] - 1)^2 + theta[["b"]]
  expect_error(
    inla_vanilla(lp_up, c(a = 0.3, b = 0.2), verbose = FALSE),
    "Hessian at the posterior mode",
    class = "rlang_error"
  )
  expect_error(
    inla_vanilla(function(theta) 0, c(a = 0.3), verbose = FALSE),
    "Hessian at the posterior mode",
    class = "rlang_error"
  )
  expect_snapshot(
    error = TRUE,
    transform = function(x) gsub("message: .*$", "message: <msg>", x),
    inla_vanilla(function(theta) 0, c(a = 0.3), verbose = FALSE)
  )
})

# (e) Progress messages and printed output -----------------------------------

test_that("verbose = TRUE reports each stage, verbose = FALSE nothing", {
  msgs <- capture_messages(
    fit <- inla_vanilla(lp_two, init_two, verbose = TRUE)
  )
  i_mode <- grep("Found the posterior mode in [0-9]+ iterations", msgs)
  i_vb <- grep("VB mean correction, max shift", msgs)
  i_sn <- grep("Fitted 2 skew-normal marginals", msgs)
  expect_length(i_mode, 1)
  expect_length(i_vb, 1)
  expect_length(i_sn, 1)
  expect_true(i_mode < i_vb && i_vb < i_sn)
  # Each line ends with the time it took
  expect_match(msgs[c(i_mode, i_vb, i_sn)], "\\[[0-9.]+m?s\\]\\s*$")
  # The result does not depend on verbose
  expect_equal(fit$summary, fit_two$summary)

  # No VB line without the VB correction
  msgs <- capture_messages(
    inla_vanilla(lp_two, init_two, vb = FALSE, verbose = TRUE)
  )
  expect_false(any(grepl("VB", msgs)))
  expect_length(grep("posterior mode|skew-normal", msgs), 2)

  # One parameter: singular
  msgs <- capture_messages(inla_vanilla(lp_gauss, c(a = 0.3), verbose = TRUE))
  expect_match(msgs, "Fitted 1 skew-normal marginal\\b", all = FALSE)

  expect_silent(inla_vanilla(lp_two, init_two, verbose = FALSE))
  expect_silent(
    inla_vanilla(lp_two, init_two, vb = FALSE, verbose = FALSE)
  )
})

test_that("print() shows the summary table and returns the fit invisibly", {
  testthat::local_reproducible_output(width = 80, crayon = FALSE)
  expect_snapshot(print(fit_two), transform = scrub_output)

  res <- NULL
  out <- capture_output(res <- withVisible(print(fit_two)))
  expect_false(res$visible)
  expect_identical(res$value, fit_two)
  expect_match(out, "2 parameters")
  expect_match(out, "VB-corrected")

  fit_novb <- inla_vanilla(lp_two, init_two, vb = FALSE, verbose = FALSE)
  out <- capture_output(print(fit_novb))
  expect_no_match(out, "VB-corrected")
  expect_match(out, "skew-normal marginals")
})

test_that("print(summary()) shows marginals, diagnostics and the evidence", {
  testthat::local_reproducible_output(width = 80, crayon = FALSE)
  expect_snapshot(print(summary(fit_two)), transform = scrub_output)

  fit_novb <- inla_vanilla(lp_two, init_two, vb = FALSE, verbose = FALSE)
  expect_snapshot(print(summary(fit_novb)), transform = scrub_output)

  fit_grad <- inla_vanilla(
    lp_pois,
    0.3,
    grad = gr_pois,
    y = pois_y,
    a0 = pois_a0,
    b0 = pois_b0,
    verbose = FALSE
  )
  out <- capture_output(print(summary(fit_grad)))
  expect_match(out, "[0-9]+ logpost and [0-9]+ grad calls")

  s <- summary(fit_two)
  res <- NULL
  out <- capture_output(res <- withVisible(print(s)))
  expect_false(res$visible)
  expect_identical(res$value, s)
})

test_that("summary() returns the table and the diagnostics", {
  s <- summary(fit_two)
  expect_s3_class(s, "summary.inla_vanilla")
  expect_identical(
    colnames(s$table),
    c("Mean", "SD", "2.5%", "25%", "50%", "75%", "97.5%", "Mode")
  )
  expect_identical(rownames(s$table), names(init_two))
  # The default columns agree with fit$summary
  shared <- colnames(fit_two$summary)
  expect_equal(s$table[, shared], fit_two$summary)
  expect_s3_class(s$diagnostics, "data.frame")
  expect_identical(
    names(s$diagnostics),
    c("vb_shift", "sn_vs_vb", "alpha", "gamma1", "nmad")
  )
  expect_identical(rownames(s$diagnostics), names(init_two))
  expect_identical(s$logml, fit_two$logml)
  expect_identical(s$logml_laplace, fit_two$logml_laplace)
})

test_that("summary(probs = ) gives the requested quantile columns", {
  s <- summary(fit_two, probs = c(0.1, 0.9))
  expect_identical(colnames(s$table), c("Mean", "SD", "10%", "90%", "Mode"))
  expect_equal(s$table[, c("Mean", "SD", "Mode")], fit_two$summary[, c(
    "Mean",
    "SD",
    "Mode"
  )])

  # The columns are quantiles of the fitted skew-normals: check the CDF at
  # them by integrating the skew-normal density written out here.
  for (j in seq_len(nrow(fit_two$sn))) {
    xi <- fit_two$sn[j, "xi"]
    omega <- fit_two$sn[j, "omega"]
    alpha <- fit_two$sn[j, "alpha"]
    dens <- function(x) {
      z <- (x - xi) / omega
      2 / omega * dnorm(z) * pnorm(alpha * z)
    }
    for (p in c(0.1, 0.9)) {
      q <- s$table[j, paste0(100 * p, "%")]
      cdf <- integrate(dens, -Inf, q, rel.tol = 1e-10)$value
      expect_equal(cdf, p, tolerance = 1e-6)
    }
  }

  # Odd probabilities are labelled as quantile() labels them
  s <- summary(fit_two, probs = c(0.005, 1 / 3))
  expect_identical(
    colnames(s$table),
    c("Mean", "SD", names(quantile(1, c(0.005, 1 / 3))), "Mode")
  )
  expect_output(print(s), "33.33333%", fixed = TRUE)
})

test_that("summary() rejects probabilities outside (0, 1)", {
  for (bad in list(0, 1, c(0.5, 1), -0.1, 1.5, NA, c(0.5, NA), "0.5", NULL)) {
    expect_error(
      summary(fit_two, probs = bad),
      "`probs` must be probabilities strictly between 0 and 1",
      class = "rlang_error"
    )
  }
  expect_snapshot(error = TRUE, summary(fit_two, probs = c(0.5, 1)))
})

test_that("summary() with no probabilities does not fail internally", {
  # Observed: probs = numeric(0) passes the check and then fails in
  # sn_summary() with base R's "length of 'dimnames' [2] not equal to array
  # extent". Either a table without quantile columns or an informative error
  # would do; this asks for the informative error, as for other bad probs.
  expect_package_error(summary(fit_two, probs = numeric(0)))
})

# (f) plot() -----------------------------------------------------------------

test_that("plot() draws all, named and indexed parameters and returns the fit", {
  grDevices::pdf(NULL)
  dev <- grDevices::dev.cur()
  on.exit(grDevices::dev.off(dev), add = TRUE)

  expect_no_error(res <- withVisible(plot(fit_two)))
  expect_false(res$visible)
  expect_identical(res$value, fit_two)

  expect_no_error(plot(fit_two, pars = "logit_p"))
  expect_no_error(plot(fit_two, pars = c("logit_p", "log_rate")))
  expect_no_error(plot(fit_two, pars = 2))
  expect_no_error(plot(fit_two, pars = c(2, 1)))
  expect_no_error(plot(fit_two, laplace = FALSE))
  expect_no_error(plot(fit_two, pars = "log_rate", laplace = FALSE))
  # Graphical parameters the method does not set itself pass through
  expect_no_error(plot(fit_two, col.axis = "grey40", cex.axis = 0.8))

  res <- withVisible(plot(fit_two, pars = 1, laplace = FALSE))
  expect_false(res$visible)
  expect_identical(res$value, fit_two)

  # More parameters than fit on one row
  fit5 <- inla_vanilla(lp_gauss, c(0.1, 0.2, 0.3, 0.4, 0.5), verbose = FALSE)
  expect_no_error(plot(fit5))
  expect_no_error(plot(fit5, pars = c("theta[2]", "theta[5]")))
})

test_that("plot() rejects parameters that are not in the fit", {
  grDevices::pdf(NULL)
  dev <- grDevices::dev.cur()
  on.exit(grDevices::dev.off(dev), add = TRUE)

  expect_error(
    plot(fit_two, pars = "sigma"),
    "Unknown parameter: \"sigma\"",
    class = "rlang_error"
  )
  expect_error(
    plot(fit_two, pars = c("sigma", "logit_p", "tau")),
    "Unknown parameters: \"sigma\" and \"tau\"",
    class = "rlang_error"
  )
  for (bad in list(0, 3, -1, c(1, 3), TRUE)) {
    expect_error(
      plot(fit_two, pars = bad),
      "`pars` must name or index parameters of the fit",
      class = "rlang_error"
    )
  }
  expect_snapshot(error = TRUE, plot(fit_two, pars = c("sigma", "tau")))
})

test_that("plot() restores the graphical parameters", {
  grDevices::pdf(NULL)
  dev <- grDevices::dev.cur()
  on.exit(grDevices::dev.off(dev), add = TRUE)

  before <- graphics::par(no.readonly = TRUE)
  keep <- setdiff(names(before), par_set_by_any_plot)
  plot(fit_two)
  expect_identical(graphics::par(no.readonly = TRUE)[keep], before[keep])
  plot(fit_two, pars = 2, laplace = FALSE)
  expect_identical(graphics::par(no.readonly = TRUE)[keep], before[keep])
  expect_error(plot(fit_two, pars = "sigma"))
  expect_identical(graphics::par(no.readonly = TRUE)[keep], before[keep])
})

test_that("plot() restores graphical parameters the caller has changed", {
  grDevices::pdf(NULL)
  dev <- grDevices::dev.cur()
  on.exit(grDevices::dev.off(dev), add = TRUE)

  # Observed: `cex` comes back as 1, not 0.8. Setting `mfrow` resets `cex`,
  # and restoring the old `mfrow` resets it again, but `cex` itself is not
  # among the parameters that plot() saves.
  graphics::par(mar = c(1, 2, 3, 4), las = 2, cex = 0.8, mgp = c(3, 1, 0))
  before <- graphics::par(no.readonly = TRUE)
  keep <- setdiff(names(before), par_set_by_any_plot)
  plot(fit_two)
  after <- graphics::par(no.readonly = TRUE)
  expect_identical(after[keep], before[keep])
})

test_that("plot() accepts axis labels and titles in ...", {
  grDevices::pdf(NULL)
  dev <- grDevices::dev.cur()
  on.exit(grDevices::dev.off(dev), add = TRUE)

  # The documentation says `...` holds further graphical parameters passed to
  # graphics::plot(). Observed: `xlab` and `main` fail with "formal argument
  # "xlab" matched by multiple actual arguments", since plot() sets both.
  expect_no_error(plot(fit_two, xlab = "unconstrained scale"))
  expect_no_error(plot(fit_two, pars = 1, main = "Rate"))
})
