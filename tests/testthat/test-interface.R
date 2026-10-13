# Tests of the user interface of inla_vanilla() and of its S3 methods: input
# validation, the checks of an analytic gradient, parameter names, arguments
# passed through `...`, log-posteriors that fail or are not finite away from
# the mode, improper posteriors and an unreached mode, progress messages, and
# the print(), summary() and plot() methods.
#
# The package raises its errors and warnings with abort() and warn()
# (R/utils.R), which give them the classes "inlavanilla_error" and
# "inlavanilla_warning"; the tests ask for those classes, so that an error
# from base R deep inside the code does not pass for a validation message.
#
# Ground truth, where a test needs one, is closed form: Gaussian posteriors,
# where the skew-normal marginals should be exact, a gamma posterior written
# on its natural scale, the gamma-Poisson model on the log scale, and a
# truncated normal; the Poisson regression behind the printed output is
# checked against a 2D grid. Tolerances on the approximate (non-Gaussian)
# cases were calibrated by printing the actual errors first; each comment
# states the observed value. The plots are checked through the display list
# that grDevices::recordPlot() keeps of the low-level graphics calls.

# Models ----------------------------------------------------------------------

# A Gaussian log-posterior with every mean equal to 1 and unit variances.
lp_gauss <- function(theta) -0.5 * sum((theta - 1)^2)
gr_gauss <- function(theta) -(theta - 1)

# y_i ~ Poisson(lambda), lambda ~ Gamma(a0, b0), on theta = log(lambda), with
# the Jacobian. Exactly lambda | y ~ Gamma(A, B) with A = a0 + sum(y) and
# B = b0 + n, so the mode of theta is log(A / B), its mean digamma(A) - log(B),
# its SD sqrt(trigamma(A)) and its quantiles log(qgamma(p, A, B)). The data
# reach logpost and grad through `...`.
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

# A Poisson regression with fixed data: two correlated (r = -0.56), skewed
# parameters, so that every diagnostic in summary() is well away from
# rounding noise (the slope's scan-path bend is flagged, at 0.14). The fit is
# deterministic. Used for the printed output and the plots.
reg_x <- c(-1.75, -1.25, -0.75, -0.25, 0.25, 0.75, 1.25, 1.75)
reg_y <- c(1, 0, 2, 1, 3, 4, 2, 6)
lp_reg <- function(theta) {
  eta <- theta[["intercept"]] + theta[["slope"]] * reg_x
  sum(dpois(reg_y, exp(eta), log = TRUE)) +
    dnorm(theta[["intercept"]], 0, 2, log = TRUE) +
    dnorm(theta[["slope"]], 0, 2, log = TRUE)
}
init_reg <- c(intercept = 0.2, slope = 0.1)
gr_reg <- function(theta) {
  eta <- theta[["intercept"]] + theta[["slope"]] * reg_x
  res <- reg_y - exp(eta)
  c(sum(res), sum(reg_x * res)) - theta / 4
}

fit_reg <- inla_vanilla(lp_reg, init_reg, silent = TRUE)
fit_reg_novb <- inla_vanilla(lp_reg, init_reg, vb = FALSE, silent = TRUE)

# Helpers ---------------------------------------------------------------------

# Make printed output stable between runs and platforms. Timings (12ms, 1.2s)
# vary. The numbers of calls and iterations follow the optimiser's path, which
# floating-point differences between platforms could change, and the gradient
# at the mode is at the level of finite-difference noise.
scrub_output <- function(x) {
  x <- gsub("[0-9.]+(ms|s)\\b", "<time>", x)
  x <- gsub("[0-9]+ (logpost|grad)\\b", "<n> \\1", x)
  x <- gsub("[0-9]+ iterations?\\b", "<n> iteration(s)", x)
  x <- gsub(
    "Newton step [-+0-9.e]+ posterior SDs",
    "Newton step <g> posterior SDs",
    x
  )
  x
}

# For the detail lines of verbose = TRUE, which show many more numbers, some
# at the level of rounding or finite-difference noise: every number with a
# decimal point or an exponent is scrubbed too (but not the thresholds, after
# "above"), along with what follows the optimiser's path (nlminb's message
# and evaluations, the calls that set the finite-difference steps, the list
# of VB steps and whether Anderson acceleration was used), and the runs of
# spaces that align the numbers in the table of marginals.
scrub_numbers <- function(x) {
  x <- scrub_output(x)
  x <- gsub("nlminb: .*, after", "nlminb: <message>, after", x)
  x <- gsub(
    "[0-9]+ function and [0-9]+ gradient",
    "<n> function and <n> gradient",
    x
  )
  x <- gsub("and [0-9]+ more", "and <n> more", x)
  x <- gsub("(Fixed-point steps[^:]*): .*$", "\\1: <steps>", x)
  x <- gsub(
    ", (Anderson acceleration at [0-9]+ of them|without Anderson acceleration)",
    ", <anderson>",
    x
  )
  x <- gsub(
    "(?<![A-Za-z0-9_.])(?<!above )-?([0-9]+\\.[0-9]+(e[-+]?[0-9]+)?|[0-9]+e[-+]?[0-9]+)",
    "#",
    x,
    perl = TRUE
  )
  x <- gsub("(\\S)  +", "\\1 ", x)
  x
}

# The diagnostics of a fit, as its warning and the flags of summary() read
# them
fit_dg <- function(fit) {
  diagnostics_internal(
    fit$diagnostics$optim,
    fit$diagnostics$vb,
    fit$diagnostics$marginals,
    pars = fit$mean_vb
  )
}

# Expect an error raised by the package itself, rather than one from base R
# deep inside the code. Unlike expect_error(class = ), a base R error is
# reported as a failure, so the rest of the test still runs.
expect_package_error <- function(object, regexp = NULL) {
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
    inherits(err, "inlavanilla_error") &&
      (is.null(regexp) || grepl(regexp, conditionMessage(err), fixed = TRUE)),
    paste0(
      "Expected an error from the package (<inlavanilla_error>)",
      if (!is.null(regexp)) paste0(" matching '", regexp, "'"),
      ", got ",
      got
    )
  )
  invisible(err)
}

# Run `expr` on a null pdf device with its display list recorded, and return
# the low-level graphics calls it made, each as list(name, args).
record_plot <- function(expr) {
  grDevices::pdf(NULL)
  dev <- grDevices::dev.cur()
  on.exit(grDevices::dev.off(dev), add = TRUE)
  grDevices::dev.control("enable")
  force(expr)
  dl <- grDevices::recordPlot()[[1]]
  lapply(dl, function(el) {
    args <- as.list(el[[2]])
    list(name = args[[1]]$name, args = args[-1])
  })
}

# The short strings drawn in a recorded plot: titles, axis labels, legend.
plot_strings <- function(calls) {
  strs <- character(0)
  walk <- function(x) {
    if (is.character(x)) {
      strs <<- c(strs, x)
    } else if (is.list(x)) {
      for (el in x) {
        walk(el)
      }
    }
  }
  walk(lapply(calls, `[[`, "args"))
  unique(strs)
}

# The skew-normal density, written out here independently of dsnorm().
sn_density <- function(x, xi, omega, alpha) {
  z <- (x - xi) / omega
  2 / omega * dnorm(z) * pnorm(alpha * z)
}

# Graphical parameters that every high-level plot sets, whatever the method
# does: the user coordinates and the axis tick positions.
par_set_by_any_plot <- c("usr", "xaxp", "yaxp")

# (a) Input validation ---------------------------------------------------------

test_that("logpost and grad must be functions", {
  expect_error(
    inla_vanilla(1, c(a = 0), silent = TRUE),
    "`logpost` must be a function, not a number",
    class = "inlavanilla_error"
  )
  expect_error(
    inla_vanilla(lp_gauss, c(a = 0), grad = "gr", silent = TRUE),
    "`grad` must be a function, not a string",
    class = "inlavanilla_error"
  )
})

test_that("init must be a non-empty, finite numeric vector", {
  for (bad in list("a", numeric(0), list(a = 1))) {
    expect_error(
      inla_vanilla(lp_gauss, bad, silent = TRUE),
      "`init` must be a non-empty numeric vector",
      class = "inlavanilla_error"
    )
  }
  for (bad in list(c(a = 0, b = NA), c(a = -Inf))) {
    expect_error(
      inla_vanilla(lp_gauss, bad, silent = TRUE),
      "`init` must contain only finite values",
      class = "inlavanilla_error"
    )
  }
})

test_that("init must be fully named, or not named at all", {
  expect_error(
    inla_vanilla(lp_gauss, c(a = 0, 1), silent = TRUE),
    "fully named, or not named at all",
    class = "inlavanilla_error"
  )
  expect_error(
    inla_vanilla(
      lp_gauss,
      stats::setNames(c(0, 1), c("a", NA)),
      silent = TRUE
    ),
    "fully named, or not named at all",
    class = "inlavanilla_error"
  )
})

test_that("repeated names in init are numbered, as TMB's obj$par needs", {
  seen <- NULL
  lp <- function(theta) {
    seen <<- union(seen, list(names(theta)))
    lp_gauss(theta)
  }
  init <- c(a = 0.2, b = 0.3, a = 0.4)
  nms <- c("a[1]", "b", "a[2]")
  expect_silent(fit <- inla_vanilla(lp, init, grad = gr_gauss, silent = TRUE))
  expect_identical(fit$warnings, character())
  expect_identical(seen, list(nms))
  expect_identical(rownames(fit$summary), nms)
  expect_identical(names(fit$mode), nms)
  expect_equal(unname(fit$mode), c(1, 1, 1), tolerance = 1e-6)

  # Unless silent = TRUE, a note says which names were changed
  msgs <- capture_messages(inla_vanilla(lp_gauss, init, grad = gr_gauss))
  expect_match(
    msgs,
    "Numbered the repeated names in `init`: \"a[1]\" and \"a[2]\"",
    fixed = TRUE,
    all = FALSE
  )

  # Names that clash after numbering are an error
  expect_error(
    inla_vanilla(lp_gauss, c(a = 0, a = 1, `a[1]` = 2), silent = TRUE),
    "`init` has names that cannot be made unique",
    class = "inlavanilla_error"
  )
})

test_that("logpost must return a single finite number at init", {
  not_single <- list(
    vector = function(theta) theta,
    string = function(theta) "-1",
    logical_na = function(theta) NA
  )
  for (lp in not_single) {
    expect_error(
      inla_vanilla(lp, c(a = 0.5, b = 0.5), silent = TRUE),
      "`logpost` must return a single number",
      class = "inlavanilla_error"
    )
  }
  not_finite <- list(
    nan = function(theta) NaN,
    minus_inf = function(theta) -Inf
  )
  for (lp in not_finite) {
    expect_error(
      inla_vanilla(lp, c(a = 0.5), silent = TRUE),
      "`logpost` is not finite at `init`",
      class = "inlavanilla_error"
    )
  }
  # Finite near the mode but not at init: the check is about init itself
  lp_hole <- function(theta) if (theta[[1]] > 5) -Inf else lp_gauss(theta)
  expect_error(
    inla_vanilla(lp_hole, c(a = 6), silent = TRUE),
    "not finite at `init`",
    class = "inlavanilla_error"
  )
  # A 1 x 1 matrix, such as a quadratic form, is a single number
  lp_form <- function(theta) -0.5 * crossprod(theta - 1)
  fit <- inla_vanilla(lp_form, c(a = 0.3, b = 1.6), silent = TRUE)
  expect_equal(unname(fit$summary[, "SD"]), c(1, 1), tolerance = 1e-6)
})

test_that("grad must return a finite vector of the right length at init", {
  init <- c(a = 0.5, b = 0.5)
  for (gr in list(function(theta) 1, function(theta) "a")) {
    expect_error(
      inla_vanilla(lp_gauss, init, grad = gr, silent = TRUE),
      "`grad` must return a numeric vector of length 2",
      class = "inlavanilla_error"
    )
  }
  for (gr in list(function(theta) c(NA, 1))) {
    expect_error(
      inla_vanilla(lp_gauss, init, grad = gr, silent = TRUE),
      "`grad` is not finite at `init`",
      class = "inlavanilla_error"
    )
  }
})

test_that("vb, verbose and silent must be a single TRUE or FALSE", {
  for (bad in list(NA, c(TRUE, FALSE))) {
    expect_error(
      inla_vanilla(lp_gauss, c(a = 0.5), vb = bad, silent = TRUE),
      "`vb` must be `TRUE` or `FALSE`",
      class = "inlavanilla_error"
    )
    expect_error(
      inla_vanilla(lp_gauss, c(a = 0.5), verbose = bad),
      "`verbose` must be `TRUE` or `FALSE`",
      class = "inlavanilla_error"
    )
    expect_error(
      inla_vanilla(lp_gauss, c(a = 0.5), silent = bad),
      "`silent` must be `TRUE` or `FALSE`",
      class = "inlavanilla_error"
    )
  }
})

test_that("verbose and silent cannot both be TRUE", {
  expect_error(
    inla_vanilla(lp_gauss, c(a = 0.5), verbose = TRUE, silent = TRUE),
    "`verbose` and `silent` cannot both be `TRUE`",
    class = "inlavanilla_error"
  )
  expect_snapshot(
    error = TRUE,
    inla_vanilla(lp_gauss, c(a = 0.5), verbose = TRUE, silent = TRUE)
  )
})

test_that("a bad verbose or silent is reported even when init has repeated names", {
  # check_init() numbers repeated names and says so unless silent is TRUE,
  # so it must not use the flags before they have been checked.
  # Observed: base R errors from `if (verbose)` inside check_init():
  # "missing value where TRUE/FALSE needed" for NA, "argument is not
  # interpretable as logical" for "yes", and "the condition has length > 1"
  # for c(TRUE, TRUE).
  for (bad in list(NA, "yes", c(TRUE, TRUE))) {
    expect_package_error(
      inla_vanilla(lp_gauss, c(a = 0, a = 1), verbose = bad),
      "`verbose` must be `TRUE` or `FALSE`"
    )
    expect_package_error(
      inla_vanilla(lp_gauss, c(a = 0, a = 1), silent = bad),
      "`silent` must be `TRUE` or `FALSE`"
    )
  }
})

test_that("control must be a list of named, known entries", {
  for (bad in list("n_grid", c(n_grid = 21))) {
    expect_error(
      inla_vanilla(lp_gauss, c(a = 0.5), control = bad, silent = TRUE),
      "`control` must be a list",
      class = "inlavanilla_error"
    )
  }
  expect_error(
    inla_vanilla(
      lp_gauss,
      c(a = 0.5),
      control = list(foo = 1),
      silent = TRUE
    ),
    "Unknown `control` entry: \"foo\"",
    class = "inlavanilla_error"
  )
  expect_error(
    inla_vanilla(
      lp_gauss,
      c(a = 0.5),
      control = list(n_grid = 21, foo = 1, ngrid = 2),
      silent = TRUE
    ),
    "Unknown `control` entries: \"foo\" and \"ngrid\"",
    class = "inlavanilla_error"
  )
  # Unnamed entries cannot be matched to anything
  for (bad in list(list(21), list(n_grid = 21, 1))) {
    expect_error(
      inla_vanilla(lp_gauss, c(a = 0.5), control = bad, silent = TRUE),
      "Every `control` entry must be named",
      class = "inlavanilla_error"
    )
  }
})

test_that("each numeric control entry is checked", {
  bad_values <- list(
    n_grid = list(6, "21", NA_real_, Inf),
    log_threshold = list(0, NA_real_),
    vb_maxit = list(0, "a"),
    vb_tol = list(-1, Inf)
  )
  messages <- c(
    n_grid = "`control$n_grid` must be a single number of at least 7",
    log_threshold = "`control$log_threshold` must be a single negative number",
    vb_maxit = "`control$vb_maxit` must be a single number of at least 1",
    vb_tol = "`control$vb_tol` must be a single positive number"
  )
  for (entry in names(bad_values)) {
    for (bad in bad_values[[entry]]) {
      ctrl <- stats::setNames(list(bad), entry)
      expect_package_error(
        inla_vanilla(lp_gauss, c(a = 0.5), control = ctrl, silent = TRUE),
        messages[[entry]]
      )
    }
  }
  expect_error(
    inla_vanilla(
      lp_gauss,
      c(a = 0.5),
      control = list(nlminb = 5),
      silent = TRUE
    ),
    "`control$nlminb` must be a list",
    fixed = TRUE,
    class = "inlavanilla_error"
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
    silent = TRUE
  )
  expect_equal(unname(fit$summary[, "Mean"]), c(1, 1), tolerance = 1e-6)
  expect_equal(unname(fit$summary[, "SD"]), c(1, 1), tolerance = 1e-6)
})

test_that("control entries reach the stages they tune", {
  # n_grid: each parameter's scan takes n_grid - 1 calls (the mode's height
  # is known), so twenty more points cost exactly 20 m more logpost calls
  # and no more gradient calls.
  counts <- function(n_grid) {
    inla_vanilla(
      lp_gauss,
      c(a = 0.3, b = 1.4),
      grad = gr_gauss,
      control = list(n_grid = n_grid),
      silent = TRUE
    )$counts
  }
  expect_identical(counts(41) - counts(21), c(logpost = 40L, grad = 0L))

  # vb_tol and vb_maxit: the VB correction of the skewed gamma-Poisson
  # posterior takes 3 fixed-point iterations at the default tolerance.
  fit_vb <- function(...) {
    inla_vanilla(
      lp_pois,
      0.3,
      y = pois_y,
      a0 = pois_a0,
      b0 = pois_b0,
      control = list(...),
      silent = TRUE
    )$diagnostics$vb
  }
  default <- fit_vb()
  expect_false(default$fallback)
  expect_lt(fit_vb(vb_tol = 0.5)$iterations, default$iterations)
  # One iteration is not enough, so the fallback optimiser runs
  expect_true(fit_vb(vb_maxit = 1)$fallback)

  # log_threshold: a threshold of -1 keeps only the top of the scan, which
  # changes the skew-normal fitted to a skewed posterior
  fit_thr <- function(thr) {
    inla_vanilla(
      lp_pois,
      0.3,
      y = pois_y,
      a0 = pois_a0,
      b0 = pois_b0,
      control = list(log_threshold = thr),
      silent = TRUE
    )$sn
  }
  expect_gt(max(abs(fit_thr(-1) - fit_thr(-6))), 1e-3)

  # nlminb: see the test of an unreached mode, where iter.max = 1 stops the
  # optimiser after one iteration
})

test_that("any accepted n_grid gives exact marginals for a Gaussian posterior", {
  # N(3, 2^2): a skew-normal with alpha = 0 is exact, so the summary should
  # be exact for every grid, odd or even (an even grid has no point at the
  # mode, and once gave the nearest point the mode's height), and for a
  # non-integer n_grid, which is rounded (7.5 to 8). Observed: exact to
  # ~1e-8.
  lp <- function(theta) dnorm(theta[[1]], 3, 2, log = TRUE)
  for (n_grid in c(7, 7.5, 20)) {
    fit <- inla_vanilla(
      lp,
      c(a = 0.5),
      control = list(n_grid = n_grid),
      silent = TRUE
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
    inla_vanilla("lp", c(a = 0))
  })
  expect_snapshot(error = TRUE, {
    inla_vanilla(function(theta) theta, c(a = 0.5, b = 0.5))
  })
  expect_snapshot(error = TRUE, {
    inla_vanilla(function(theta) -Inf, c(a = 0.5))
  })
  expect_snapshot(error = TRUE, {
    inla_vanilla(lp_gauss, c(a = 0.5), control = list(foo = 1, ngrid = 21))
  })
})

test_that("errors do not name the package's internal helpers", {
  # Users call inla_vanilla(); an error should not point them at
  # check_init() or vanilla_control(). The call is left out, or is
  # inla_vanilla() itself.
  errs <- list(
    function() inla_vanilla(lp_gauss, "a"),
    function() inla_vanilla(lp_gauss, c(a = 0), vb = NA),
    function() inla_vanilla(lp_gauss, c(a = 0.5), grad = function(theta) NA),
    function() inla_vanilla(function(theta) 0, c(a = 0.3), silent = TRUE),
    function() summary(fit_reg, probs = 2),
    function() {
      grDevices::pdf(NULL)
      on.exit(grDevices::dev.off())
      plot(fit_reg, pars = "sigma")
    }
  )
  for (f in errs) {
    err <- tryCatch(f(), error = identity)
    expect_s3_class(err, "inlavanilla_error")
    call <- conditionCall(err)
    expect_true(
      is.null(call) || identical(call, quote(inla_vanilla())),
      label = paste("the call", deparse(call), "is NULL or inla_vanilla()")
    )
  }
})

test_that("an error inside logpost or grad says where the fit stopped", {
  # logpost fails beyond 3, which only the scan of the marginal reaches (the
  # mode is 1 and the scan goes to 1 + 4 SDs)
  lp_fails <- function(theta) {
    if (theta[[1]] > 3) {
      stop("cannot evaluate here")
    }
    lp_gauss(theta)
  }
  err <- expect_error(
    inla_vanilla(lp_fails, c(a = 0.5), silent = TRUE),
    "while fitting the skew-normal marginals",
    class = "inlavanilla_error"
  )
  msg <- conditionMessage(err)
  # The value reported is one where logpost fails
  at <- regmatches(
    msg,
    regexpr("tried was a = -?[0-9]+(\\.[0-9]+)?(e[-+]?[0-9]+)?", msg)
  )
  expect_length(at, 1)
  expect_gt(as.numeric(sub("tried was a = ", "", at)), 3)
  expect_match(msg, "return `-Inf` there instead of failing", fixed = TRUE)
  expect_match(msg, "cannot evaluate here", fixed = TRUE)
  # The original error is kept
  expect_identical(conditionMessage(err$parent), "cannot evaluate here")

  # At init, for logpost and for grad
  expect_error(
    inla_vanilla(function(theta) stop("bad data"), c(a = 0.5), silent = TRUE),
    "while evaluating logpost at init",
    class = "inlavanilla_error"
  )
  expect_error(
    inla_vanilla(
      lp_gauss,
      c(a = 0.5),
      grad = function(theta) stop("bad data"),
      silent = TRUE
    ),
    "while evaluating grad at init",
    class = "inlavanilla_error"
  )

  # The package's own errors pass through unwrapped
  err <- expect_error(
    inla_vanilla(function(theta) 0, c(a = 0.3), silent = TRUE),
    class = "inlavanilla_error"
  )
  expect_no_match(conditionMessage(err), "The fit stopped")
})

# (b) The gradient checks ------------------------------------------------------

test_that("a wrong analytic gradient triggers the mismatch warning", {
  # Right for `a`, twice too large for `b`: the warning names `b`. (A sign
  # error would also be caught, but then the optimiser climbs away from the
  # mode and the fit fails later with the Hessian error.)
  # (Warnings are tested at the default output level, without silent = TRUE,
  # which would muffle them; its stage lines are suppressed.)
  gr_wrong <- function(theta) c(-(theta[[1]] - 1), -2 * (theta[[2]] - 1))
  expect_snapshot(
    invisible(suppressMessages(
      inla_vanilla(lp_gauss, c(a = 0.3, b = 0.2), grad = gr_wrong)
    ))
  )
  # A gradient off by a constant factor is also caught
  expect_warning(
    suppressMessages(inla_vanilla(
      lp_pois,
      0.3,
      grad = function(theta, ...) 2 * gr_pois(theta, ...),
      y = pois_y,
      a0 = pois_a0,
      b0 = pois_b0
    )),
    "does not match finite differences",
    class = "inlavanilla_warning"
  )
})

test_that("the gradient check works whatever the scale of the posterior", {
  # theta ~ N(s, s^2). A gradient twice too large has the right zero, so the
  # mode is right and nothing later notices, but the Hessian taken from it is
  # twice too large: the Laplace SD is off by sqrt(2), and with it fit$cov,
  # the copula and the log marginal likelihood. The check once compared
  # gradients with an absolute tolerance of 1e-3 (1 + |g|), which hid the
  # mismatch at s = 1000, where the gradient is below 1e-3: the Laplace SD
  # came out as 707 instead of 1000, with no warning.
  # The right gradient passes silently at every scale, and the Laplace SD
  # is then exact.
  for (s in c(1e-3, 1, 1000)) {
    lp <- function(theta) -0.5 * (theta[[1]] / s - 1)^2
    gr_right <- function(theta) -(theta[[1]] / s - 1) / s
    gr_twice <- function(theta) 2 * gr_right(theta)
    expect_warning(
      suppressMessages(inla_vanilla(lp, c(a = 0.5 * s), grad = gr_twice)),
      "does not match finite differences",
      class = "inlavanilla_warning",
      label = paste("a gradient twice too large at posterior SD", s)
    )
    # silent = TRUE muffles warnings, but records them in fit$warnings
    expect_silent(
      fit <- inla_vanilla(lp, c(a = 0.5 * s), grad = gr_right, silent = TRUE)
    )
    expect_identical(fit$warnings, character())
    expect_equal(sqrt(fit$cov[[1]]), s, tolerance = 1e-6)
  }
})

test_that("a gradient that is not a gradient field is caught at the mode", {
  # A rotation added to the right gradient: zero at the mode (1, 1), where
  # the fit starts, so the check at init passes, but it is wrong one
  # curvature scale away, and its Jacobian is not symmetric.
  gr_rot <- function(theta) {
    d <- theta - 1
    c(-d[[1]] + 0.5 * d[[2]], -d[[2]] - 0.5 * d[[1]])
  }
  fit <- inla_vanilla(lp_gauss, c(a = 1, b = 1), grad = gr_rot, silent = TRUE)
  expect_length(fit$warnings, 2)
  expect_match(
    fit$warnings[[1]],
    "does not match finite differences of `logpost` one curvature scale"
  )
  expect_match(
    fit$warnings[[2]],
    "does not look like the gradient of `logpost`"
  )
})

test_that("a correct analytic gradient passes the checks silently", {
  # silent = TRUE muffles warnings, but records them in fit$warnings
  fit <- inla_vanilla(
    lp_gauss,
    c(a = 0.3, b = 0.2),
    grad = gr_gauss,
    silent = TRUE
  )
  expect_identical(fit$warnings, character())
  fit <- inla_vanilla(
    lp_pois,
    0.3,
    grad = gr_pois,
    y = pois_y,
    a0 = pois_a0,
    b0 = pois_b0,
    silent = TRUE
  )
  expect_identical(fit$warnings, character())
  # Also for a Student-t posterior (5 df) with scale 1e-6, far from
  # quadratic beyond a few scale units. The check once took differences at a
  # fixed step of 1e-5, ten scale units here, which gave a gradient 20 times
  # too small at init and flagged this correct one.
  s <- 1e-6
  lp_t <- function(theta) -3 * log1p(((theta[[1]] - 1) / s)^2 / 5)
  gr_t <- function(theta) {
    z <- (theta[[1]] - 1) / s
    -6 * z / (5 * s * (1 + z^2 / 5))
  }
  fit <- inla_vanilla(lp_t, c(a = 1 + 0.5 * s), grad = gr_t, silent = TRUE)
  # No gradient warning. (A skew-normal cannot follow the t tails, so the
  # fit-time diagnostics warning reports its NMAD; that is tested elsewhere.)
  expect_false(any(grepl(
    "does not match finite differences|does not look like the gradient",
    fit$warnings
  )))
  # A 1 x m matrix, as TMB's obj$gr() returns, is accepted
  fit <- inla_vanilla(
    lp_gauss,
    c(a = 0.3, b = 0.2),
    grad = function(theta) matrix(gr_gauss(theta), nrow = 1),
    silent = TRUE
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
  fit <- inla_vanilla(lp, c(0.1, 0.2, 0.3), grad = gr, silent = TRUE)

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
  fit <- inla_vanilla(lp_gauss, c(alpha = 0L, beta = 2L), silent = TRUE)
  expect_identical(rownames(fit$summary), c("alpha", "beta"))
  expect_equal(unname(fit$mode), c(1, 1), tolerance = 1e-6)
})

test_that("arguments in ... reach logpost and grad", {
  # The mode of log(lambda) is log(A / B) exactly; it depends on all of y, a0
  # and b0, so the fit can only get it right if they arrive. Observed: the
  # mode is exact to 3e-9 without a gradient and 2e-10 with one.
  fit_fd <- inla_vanilla(
    lp_pois,
    0.3,
    y = pois_y,
    a0 = pois_a0,
    b0 = pois_b0,
    silent = TRUE
  )
  fit_gr <- inla_vanilla(
    lp_pois,
    0.3,
    grad = gr_pois,
    y = pois_y,
    a0 = pois_a0,
    b0 = pois_b0,
    silent = TRUE
  )
  exact_mode <- log(pois_A / pois_B)
  expect_equal(fit_fd$mode[[1]], exact_mode, tolerance = 1e-7)
  expect_equal(fit_gr$mode[[1]], exact_mode, tolerance = 1e-8)
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
    silent = TRUE
  )
  expect_equal(
    fit2$mode[[1]],
    log((pois_a0 + sum(y2)) / (pois_b0 + length(y2))),
    tolerance = 1e-8
  )

  # Without them logpost fails, and the error says so
  expect_error(
    inla_vanilla(lp_pois, 0.3, silent = TRUE),
    "while evaluating logpost at init",
    class = "inlavanilla_error"
  )
})

test_that("data arguments named like vb, verbose or control reach logpost", {
  # The arguments after `...` must be named in full, so v, verb, sil, con
  # and contr are data, not abbreviations of vb, verbose, silent and control
  # (verbose and silent are left at their defaults on purpose). The mode is
  # v + verb + sil + con + contr = 15.
  lp <- function(theta, v, verb, sil, con, contr) {
    -0.5 * (theta[[1]] - v - verb - sil - con - contr)^2
  }
  fit <- suppressMessages(
    inla_vanilla(lp, c(a = 0.5), v = 1, verb = 2, sil = 5, con = 3, contr = 4)
  )
  expect_equal(fit$mode[["a"]], 15, tolerance = 1e-6)
  expect_true(fit$vb)
  expect_false(fit$verbose)
  expect_false(fit$silent)
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

  expect_silent(fit <- inla_vanilla(lp_nan, c(x = 4), silent = TRUE))
  expect_identical(fit$warnings, character())
  expect_gt(n_outside, 0)
  expect_true(all(is.finite(fit$summary)))
  expect_true(all(is.finite(fit$sn)))
  expect_true(is.finite(fit$logml))

  exact_mean <- shape / rate
  exact_sd <- sqrt(shape) / rate
  exact_q <- qgamma(c(0.025, 0.5, 0.975), shape, rate)
  lap_sd <- sqrt(fit$cov[1, 1])
  lap_q <- qnorm(c(0.025, 0.5, 0.975), fit$mode[[1]], lap_sd)
  # Observed: mean 0.0055 SDs off, SD 0.23% small, quantiles at most 0.016
  # SDs off. Laplace: mean (the mode) 0.32 SDs off, SD 5.1% small, quantiles
  # up to 0.70 SDs off.
  err_mean <- abs(fit$summary[1, "Mean"] - exact_mean) / exact_sd
  err_sd <- abs(fit$summary[1, "SD"] / exact_sd - 1)
  err_q <- max(abs(fit$summary[1, c("2.5%", "50%", "97.5%")] - exact_q)) /
    exact_sd
  expect_lt(err_mean, 0.015)
  expect_lt(err_sd, 0.006)
  expect_lt(err_q, 0.04)
  expect_lt(err_mean, abs(fit$mode[[1]] - exact_mean) / exact_sd)
  expect_lt(err_sd, abs(lap_sd / exact_sd - 1))
  expect_lt(err_q, max(abs(lap_q - exact_q)) / exact_sd)

  # Starting near the boundary, where the optimiser's first steps can leave
  # the support, finds the same mode. Observed: exact to 4e-9.
  fit_near <- inla_vanilla(lp_nan, c(x = 0.1), silent = TRUE)
  expect_equal(fit_near$mode[[1]], (shape - 1) / rate, tolerance = 1e-7)

  # -Inf and NaN are treated alike
  fit_inf <- inla_vanilla(lp_inf, c(x = 4), silent = TRUE)
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
  expect_silent(fit <- inla_vanilla(lp, c(a = 0.3, b = -0.2), silent = TRUE))
  expect_identical(fit$warnings, character())
  expect_gt(n_outside, 0)
  expect_true(all(is.finite(fit$summary)))
  expect_true(all(fit$diagnostics$marginals$tails))

  # Exact truncated-normal moments of `a`, and of b = a / 2 + e with
  # e ~ N(0, 0.75) independent of a
  lo <- -2.5 # truncation point in standard units
  lambda <- dnorm(lo) / (1 - pnorm(lo))
  mean_a <- 1 + lambda
  sd_a <- sqrt(1 + lo * lambda - lambda^2)
  mean_b <- mean_a / 2
  sd_b <- sqrt(sd_a^2 / 4 + 0.75)
  # Observed: a's mean 0.018 SDs off and SD 2.3% too large, b's mean 0.0089
  # SDs off and SD 0.56% too large. The skew-normal sees only the Gaussian
  # inside the region (the points outside are dropped), so it matches the
  # Laplace approximation and cannot represent the hard edge; the errors are
  # the truncation's own effect.
  expect_lt(abs(fit$summary["a", "Mean"] - mean_a) / sd_a, 0.05)
  expect_lt(abs(fit$summary["a", "SD"] / sd_a - 1), 0.05)
  expect_lt(abs(fit$summary["b", "Mean"] - mean_b) / sd_b, 0.025)
  expect_lt(abs(fit$summary["b", "SD"] / sd_b - 1), 0.015)
})

test_that("a one-parameter model works end to end", {
  fit <- inla_vanilla(
    lp_pois,
    0.3,
    grad = gr_pois,
    y = pois_y,
    a0 = pois_a0,
    b0 = pois_b0,
    silent = TRUE
  )
  expect_s3_class(fit, "inla_vanilla")
  expect_identical(dim(fit$summary), c(1L, 6L))
  expect_identical(dim(fit$sn), c(1L, 3L))
  expect_identical(dim(fit$cov), c(1L, 1L))
  expect_identical(rownames(fit$summary), "theta[1]")

  # Against the exact log-gamma posterior. Observed: mean 0.00066 SDs off,
  # SD 0.011% off, quantiles at most 0.0051 SDs off; the Laplace mean is
  # 0.13 SDs off. The log marginal likelihood is 0.0028 off the exact
  # log p(y) = -sum(log y!) + a0 log b0 - log Gamma(a0) + log Gamma(A)
  # - A log B.
  exact_mean <- digamma(pois_A) - log(pois_B)
  exact_sd <- sqrt(trigamma(pois_A))
  exact_q <- log(qgamma(c(0.025, 0.5, 0.975), pois_A, pois_B))
  expect_lt(abs(fit$summary[1, "Mean"] - exact_mean) / exact_sd, 0.002)
  expect_lt(abs(fit$summary[1, "SD"] / exact_sd - 1), 3e-4)
  expect_lt(
    max(abs(fit$summary[1, c("2.5%", "50%", "97.5%")] - exact_q)) / exact_sd,
    0.015
  )
  exact_logml <- -sum(lgamma(pois_y + 1)) +
    pois_a0 * log(pois_b0) -
    lgamma(pois_a0) +
    lgamma(pois_A) -
    pois_A * log(pois_B)
  expect_lt(abs(fit$logml - exact_logml), 0.008)

  out <- capture_output(print(fit))
  expect_match(out, "theta[1]", fixed = TRUE)
  out <- capture_output(print(summary(fit)))
  expect_match(out, "theta[1]", fixed = TRUE)
  # No scan-path bend with one parameter
  expect_identical(summary(fit)$diagnostics$bend, 0)

  calls <- record_plot(plot(fit))
  expect_true("theta[1]" %in% plot_strings(calls))
  expect_no_error(record_plot(plot(fit, laplace = FALSE)))

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

# (d) Improper or non-identified posteriors, and an unreached mode ------------

test_that("a posterior flat in one direction fails with the Hessian error", {
  lp_flat <- function(theta) -0.5 * (theta[["a"]] - 1)^2
  gr_flat <- function(theta) c(-(theta[["a"]] - 1), 0)
  expect_error(
    inla_vanilla(lp_flat, c(a = 0.3, b = 0.2), silent = TRUE),
    "Hessian at the posterior mode is not positive definite",
    class = "inlavanilla_error"
  )
  expect_error(
    inla_vanilla(lp_flat, c(a = 0.3, b = 0.2), grad = gr_flat, silent = TRUE),
    "Hessian at the posterior mode is not positive definite",
    class = "inlavanilla_error"
  )
})

test_that("a non-identified posterior fails with the Hessian error", {
  # Only a + b is identified
  lp_sum <- function(theta) -0.5 * (theta[["a"]] + theta[["b"]] - 1)^2
  gr_sum <- function(theta) rep(-(theta[["a"]] + theta[["b"]] - 1), 2)
  expect_error(
    inla_vanilla(lp_sum, c(a = 0.3, b = 0.2), silent = TRUE),
    "Hessian at the posterior mode is not positive definite",
    class = "inlavanilla_error"
  )
  expect_error(
    inla_vanilla(lp_sum, c(a = 0.3, b = 0.2), grad = gr_sum, silent = TRUE),
    "Hessian at the posterior mode is not positive definite",
    class = "inlavanilla_error"
  )
})

test_that("an improper posterior with no mode fails with the Hessian error", {
  # Increasing without bound in b, and flat (logpost = 0) in one dimension
  lp_up <- function(theta) -0.5 * (theta[["a"]] - 1)^2 + theta[["b"]]
  expect_error(
    inla_vanilla(lp_up, c(a = 0.3, b = 0.2), silent = TRUE),
    "Hessian at the posterior mode",
    class = "inlavanilla_error"
  )
  expect_snapshot(
    error = TRUE,
    transform = function(x) gsub("message: .*$", "message: <msg>", x),
    inla_vanilla(function(theta) 0, c(a = 0.3), silent = TRUE)
  )
})

test_that("an unreached mode gives a warning, a flag and a yellow stage line", {
  # One nlminb() iteration from init = 3, 9 posterior SDs above the mode
  # log(A / B) = 0.81, cannot reach it. The mode comes first in the one
  # diagnostics warning, with nlminb()'s message and a hint; the bad mode
  # also spoils the marginal, which the same warning reports.
  msgs <- NULL
  w <- expect_warning(
    msgs <- capture_messages(
      fit <- inla_vanilla(
        lp_pois,
        3,
        y = pois_y,
        a0 = pois_a0,
        b0 = pois_b0,
        control = list(nlminb = list(iter.max = 1))
      )
    ),
    class = "inlavanilla_diagnostics_warning"
  )
  text <- cli::ansi_strip(conditionMessage(w))
  expect_match(
    text,
    "x The gradient at the posterior mode is not zero \\(max \\|grad\\| = [0-9.]+\\): a Newton step"
  )
  expect_match(text, "Raise `iter.max` or `eval.max`")
  # Observed 3.4 posterior SDs
  expect_gt(max(fit$diagnostics$marginals$mode_shift_sigma), 1)
  expect_output(
    print(summary(fit)),
    "Mode search stopped short (\"iteration limit",
    fixed = TRUE
  )
  expect_match(msgs, "! Stopped near the posterior mode", all = FALSE)
})

test_that("nlminb's false convergence at an exact mode gives no warning", {
  # Starting at the mode of a centred quadratic, nlminb() without a gradient
  # reports "false convergence (8)" (observed here; it may differ between
  # platforms). The Newton step there is zero, so the package must not warn,
  # whatever nlminb() says; diagnostics() still reports nlminb()'s code as
  # INLAvaan defines it.
  expect_no_warning(
    fit <- suppressMessages(
      inla_vanilla(function(theta) -0.5 * sum(theta^2), c(a = 0, b = 0))
    )
  )
  expect_identical(fit$warnings, character())
  expect_equal(unname(fit$summary[, "SD"]), c(1, 1), tolerance = 1e-6)
  expect_equal(unname(fit$summary[, "Mean"]), c(0, 0), tolerance = 1e-6)
})

# (e) Progress messages and printed output -----------------------------------

test_that("by default each stage prints one line, and nothing else", {
  msgs <- capture_messages(fit <- inla_vanilla(lp_reg, init_reg))
  i_mode <- grep("Found the posterior mode in [0-9]+ iterations", msgs)
  i_vb <- grep("VB mean correction, max shift", msgs)
  i_sn <- grep("Fitted 2 skew-normal marginals", msgs)
  expect_length(msgs, 3)
  expect_length(i_mode, 1)
  expect_length(i_vb, 1)
  expect_length(i_sn, 1)
  expect_true(i_mode < i_vb && i_vb < i_sn)
  # Each line ends with the time it took
  expect_match(msgs, "\\[[0-9.]+m?s\\]\\s*$")
  # The lines are messages: nothing goes to standard output
  expect_identical(
    capture_output(suppressMessages(inla_vanilla(lp_reg, init_reg))),
    ""
  )
  expect_false(fit$verbose)
  expect_false(fit$silent)

  # No VB line without the VB correction
  msgs <- capture_messages(inla_vanilla(lp_reg, init_reg, vb = FALSE))
  expect_false(any(grepl("VB", msgs)))
  expect_length(msgs, 2)
  expect_length(grep("Found the posterior mode|Fitted 2 skew-normal", msgs), 2)

  # One parameter: singular
  msgs <- capture_messages(inla_vanilla(lp_gauss, c(a = 0.3)))
  expect_match(msgs, "Fitted 1 skew-normal marginal\\b", all = FALSE)
})

test_that("progress bars follow the user's cli.progress_show_after", {
  # A knitr or Quarto setup chunk sets the option to keep bars, which cannot
  # redraw in place there, out of the document, so the package must not
  # override it: it sets half a second only when the option is unset, and
  # restores the option afterwards. cli reads the delay when a bar is
  # created, so that is when it is recorded here.
  seen <- NULL
  create_bar <- cli::cli_progress_bar
  local_mocked_bindings(
    cli_progress_bar = function(...) {
      seen <<- c(seen, getOption("cli.progress_show_after", NA))
      create_bar(...)
    },
    .package = "cli"
  )
  # add = TRUE: a plain on.exit() would replace the exit handler that undoes
  # the mock, which then leaks into later tests
  old <- options(cli.progress_show_after = NULL)
  on.exit(options(old), add = TRUE)

  # Unset: half a second for each of the three stages, then unset again
  suppressMessages(inla_vanilla(lp_reg, init_reg))
  expect_identical(seen, rep(0.5, 3))
  expect_null(getOption("cli.progress_show_after"))

  # Set by the user: kept
  seen <- NULL
  options(cli.progress_show_after = 2)
  suppressMessages(inla_vanilla(lp_reg, init_reg))
  expect_identical(seen, rep(2, 3))
  expect_identical(getOption("cli.progress_show_after"), 2)

  # No bars at all under silent = TRUE
  seen <- NULL
  inla_vanilla(lp_reg, init_reg, silent = TRUE)
  expect_null(seen)
})

test_that("verbose = TRUE adds detail lines under the stage lines", {
  stage_text <- function(x) sub(" \\[[0-9.]+m?s\\]\\s*$", "", x)
  default <- stage_text(capture_messages(inla_vanilla(lp_reg, init_reg)))
  msgs <- capture_messages(
    fit <- inla_vanilla(lp_reg, init_reg, verbose = TRUE)
  )
  expect_true(fit$verbose)
  expect_false(fit$silent)
  # The stage lines are the default ones, in the same order. Every other
  # line is an indented detail line, or a line of its own marked with an i.
  is_stage <- stage_text(msgs) %in% default
  expect_identical(stage_text(msgs)[is_stage], default)
  expect_match(msgs[!is_stage], "^(  |i )")
  expect_gt(sum(!is_stage), 15)

  text <- paste(msgs, collapse = "")
  for (pattern in c(
    "i 2 parameters: intercept and slope",
    "Derivatives by finite differences of logpost",
    "logpost at init -?[0-9.]+",
    "nlminb: .+, after [0-9]+ iterations? \\([0-9]+ function and [0-9]+ gradient evaluations\\)",
    "logpost at the mode -?[0-9.]+, whitened gradient",
    "Hessian by second differences of logpost: 6 calls",
    "Posterior correlation matrix: largest correlation -0.5[0-9]* \\(intercept, slope\\), smallest eigenvalue 0.4",
    "Laplace log marginal likelihood -?[0-9.]+",
    "Fixed-point steps, largest over the parameters \\(SD\\): ",
    "KL term [0-9.e-]+, giving a log marginal likelihood of -?[0-9.]+",
    "Largest shift -0.2[0-9]* SD, for intercept",
    "Shape +NMAD +Bend +Tails +End mass +gamma1",
    "marks what summary\\(\\) flags",
    "i Done: [0-9]+ logpost calls in [0-9.]+m?s"
  )) {
    expect_match(text, pattern)
  }
  # The log marginal likelihoods shown are the fit's own
  expect_match(
    text,
    paste(
      "Laplace log marginal likelihood",
      format(round(fit$logml_laplace, 2), nsmall = 2)
    )
  )
  expect_match(
    text,
    paste(
      "log marginal likelihood of",
      format(round(fit$logml, 2), nsmall = 2)
    )
  )
  # The slope's bend is flagged, as in summary(), and nothing else
  rows <- grep("^  (intercept|slope) ", msgs, value = TRUE)
  expect_length(rows, 2)
  expect_match(rows[[2]], "^  slope .*[0-9]!")
  expect_no_match(rows[[1]], "!")
  # The calls of the stages add up to the total
  per_stage <- as.integer(sub(
    "^  ([0-9]+) logpost calls\n$",
    "\\1",
    grep("^  [0-9]+ logpost calls\n$", msgs, value = TRUE)
  ))
  expect_length(per_stage, 3)
  expect_identical(sum(per_stage), fit$counts[["logpost"]])
  expect_match(text, paste("Done:", fit$counts[["logpost"]], "logpost calls"))

  # With a gradient: its check at init and the asymmetry of its Jacobian
  text <- paste(
    capture_messages(
      inla_vanilla(lp_reg, init_reg, grad = gr_reg, verbose = TRUE)
    ),
    collapse = ""
  )
  expect_match(
    text,
    "Gradient from grad, checked at init: relative discrepancy [0-9.e-]+ \\(warning above 0.001\\), [0-9]+ logpost calls"
  )
  expect_match(text, "Hessian by central differences of grad: 4 calls")
  expect_match(text, "Jacobian asymmetry [0-9.e-]+ \\(warning above 0.001\\)")
  expect_match(text, "Done: [0-9]+ logpost and [0-9]+ grad calls")

  # Without VB or scans, it says so
  msgs <- capture_messages(inla_vanilla(
    lp_reg,
    init_reg,
    vb = FALSE,
    marginals = "gaussian",
    verbose = TRUE
  ))
  expect_match(msgs, "^i No VB correction \\(vb = FALSE\\)", all = FALSE)
  expect_match(msgs, "^i Gaussian marginals .*: no scans", all = FALSE)
  expect_no_match(msgs, "Fixed-point|Shape")

  # One parameter: no correlations
  msgs <- capture_messages(inla_vanilla(lp_gauss, c(a = 0.3), verbose = TRUE))
  expect_match(msgs, "^i 1 parameter: a\n$", all = FALSE)
  expect_no_match(msgs, "correlation")
})

test_that("verbose output keeps its layout", {
  testthat::local_reproducible_output(width = 80, crayon = FALSE)
  expect_snapshot(
    invisible(inla_vanilla(lp_reg, init_reg, verbose = TRUE)),
    transform = scrub_numbers
  )
  expect_snapshot(
    invisible(inla_vanilla(lp_reg, init_reg, grad = gr_reg, verbose = TRUE)),
    transform = scrub_numbers
  )
})

test_that("verbose = TRUE reports a restart, a VB fallback and heavy tails", {
  # The unreached mode of the test above: nlminb() is restarted with the
  # parameters rescaled, the VB correction falls back to nlminb(), and the
  # scan ends well above the threshold.
  msgs <- capture_messages(suppressWarnings(inla_vanilla(
    lp_pois,
    3,
    y = pois_y,
    a0 = pois_a0,
    b0 = pois_b0,
    control = list(nlminb = list(iter.max = 1)),
    verbose = TRUE
  )))
  text <- paste(msgs, collapse = "")
  expect_match(text, "Restarted nlminb with each parameter scaled")
  expect_match(
    text,
    "whitened gradient [0-9.e-]+ SD, Newton step [0-9.e-]+ SD! \\(warning above"
  )
  expect_match(text, "! Fell back to nlminb, as a step was over 1 posterior SD")
  expect_match(text, "heavy!")
  # the unreached mode and what it leads to, in one diagnostics warning
  expect_match(text, "Done: .*, with 1 warning\n")
})

test_that("silent = TRUE prints nothing at all", {
  expect_silent(fit <- inla_vanilla(lp_reg, init_reg, silent = TRUE))
  expect_true(fit$silent)
  expect_false(fit$verbose)
  expect_identical(fit$warnings, character())
  # Nor the note on repeated names, nor the messages and warnings of logpost
  lp_chatty <- function(theta) {
    message("evaluating logpost")
    warning("logpost is chatty")
    lp_gauss(theta)
  }
  expect_silent(fit <- inla_vanilla(lp_chatty, c(a = 0, a = 1), silent = TRUE))
  expect_identical(names(fit$mode), c("a[1]", "a[2]"))
  expect_identical(fit$warnings, "logpost is chatty")
  # Errors are still raised
  expect_error(
    inla_vanilla(function(theta) stop("bad data"), c(a = 0.5), silent = TRUE),
    "while evaluating logpost at init",
    class = "inlavanilla_error"
  )
})

test_that("silent = TRUE records each warning once, as plain text", {
  # An unreached mode (as above), and a logpost that warns at every call
  fit_short <- function() {
    inla_vanilla(
      lp_pois,
      3,
      y = pois_y,
      a0 = pois_a0,
      b0 = pois_b0,
      control = list(nlminb = list(iter.max = 1)),
      silent = TRUE
    )
  }
  expect_silent(fit <- fit_short())
  # The unreached mode and what it leads to, in one diagnostics warning
  expect_length(fit$warnings, 1)
  expect_match(
    fit$warnings[1],
    "^Fit diagnostics flagged [0-9] potential issues:\n"
  )
  expect_match(
    fit$warnings[1],
    "Raise `iter.max` or `eval.max` in `control$nlminb`.",
    fixed = TRUE
  )
  # The diagnostics show the same unreached mode as without silent
  expect_gt(max(fit$diagnostics$marginals$mode_shift_sigma), 1)

  # Plain text, whatever the console: no colours, no cli markup and no line
  # wrapping, so the headline, each bullet and the closing pointer are a line
  # each
  local({
    testthat::local_reproducible_output(width = 30, crayon = TRUE)
    # testthat leaves conditions unwrapped; wrap them at 30 characters
    old <- options(cli.condition_width = 30)
    on.exit(options(old), add = TRUE)
    fit <- fit_short()
    expect_false(any(grepl("\033", fit$warnings, fixed = TRUE)))
    expect_false(any(grepl("{.", fit$warnings, fixed = TRUE)))
    lines <- strsplit(fit$warnings[1], "\n")[[1]]
    expect_match(lines[-1], "^[xi] ")
    expect_length(lines, length(diagnostic_problems(fit_dg(fit), FALSE)) + 2)
    # whereas the warning as signalled is styled and wrapped
    w <- tryCatch(
      suppressMessages(inla_vanilla(
        lp_pois,
        3,
        y = pois_y,
        a0 = pois_a0,
        b0 = pois_b0,
        control = list(nlminb = list(iter.max = 1))
      )),
      inlavanilla_warning = identity
    )
    expect_true(grepl("\033", conditionMessage(w), fixed = TRUE))
    expect_gt(length(strsplit(conditionMessage(w), "\n")[[1]]), length(lines))
  })

  n_calls <- 0
  lp_warns <- function(theta) {
    n_calls <<- n_calls + 1
    warning("logpost is chatty")
    lp_gauss(theta)
  }
  fit <- inla_vanilla(lp_warns, c(a = 0.3, b = 0.2), silent = TRUE)
  expect_gt(n_calls, 50)
  expect_identical(fit$warnings, "logpost is chatty")
})

test_that("without silent = TRUE, warnings are signalled and also recorded", {
  w <- expect_warning(
    fit <- suppressMessages(inla_vanilla(
      lp_pois,
      3,
      y = pois_y,
      a0 = pois_a0,
      b0 = pois_b0,
      control = list(nlminb = list(iter.max = 1))
    )),
    "The gradient at the posterior mode is not zero",
    class = "inlavanilla_diagnostics_warning"
  )
  # The record is the warning's message, as plain text
  squish <- function(x) gsub("\\s+", " ", cli::ansi_strip(x))
  expect_length(fit$warnings, 1)
  expect_identical(squish(fit$warnings[1]), squish(conditionMessage(w)))

  # Every warning of logpost is signalled, and recorded once
  n_calls <- 0
  lp_warns <- function(theta) {
    n_calls <<- n_calls + 1
    warning("logpost is chatty")
    lp_gauss(theta)
  }
  n_signalled <- 0
  fit <- withCallingHandlers(
    suppressMessages(inla_vanilla(lp_warns, c(a = 0.3, b = 0.2))),
    warning = function(w) {
      n_signalled <<- n_signalled + 1
      invokeRestart("muffleWarning")
    }
  )
  expect_identical(n_signalled, n_calls)
  expect_identical(fit$warnings, "logpost is chatty")
})

test_that("the output level does not change the fit", {
  # Only the record of the settings, the timing and the call differ
  strip <- function(fit) {
    fit[c("verbose", "silent", "timing", "call")] <- NULL
    fit
  }
  settings <- list(
    list(),
    list(grad = gr_reg),
    list(vb = FALSE, marginals = "gaussian")
  )
  for (s in settings) {
    fit_at <- function(...) {
      suppressMessages(do.call(
        inla_vanilla,
        c(list(lp_reg, init_reg), s, list(...))
      ))
    }
    base <- strip(fit_at(silent = TRUE))
    expect_identical(strip(fit_at()), base)
    expect_identical(strip(fit_at(verbose = TRUE)), base)
  }
})

test_that("print() shows the summary table and returns the fit invisibly", {
  testthat::local_reproducible_output(width = 80, crayon = FALSE)
  expect_snapshot(print(fit_reg), transform = scrub_output)

  res <- NULL
  out <- capture_output(res <- withVisible(print(fit_reg)))
  expect_false(res$visible)
  expect_identical(res$value, fit_reg)
  # The rule names the approximation, not the time it took
  expect_match(out, "^-- INLAvanilla -+ Skew-normal \\+ VB --\n")
  expect_no_match(out, "[0-9](ms|s)\\b|parameters")
  # The log marginal likelihood shown is fit$logml, to 2 decimals
  expect_match(out, format(round(fit_reg$logml, 2), nsmall = 2), fixed = TRUE)

  out <- capture_output(print(fit_reg_novb))
  expect_match(out, "^-- INLAvanilla -+ Skew-normal --\n")
  expect_match(
    capture_output(print(summary(fit_reg_novb))),
    "-- Skew-normal --"
  )
})

test_that("print() points to the fit's warnings and to flagged checks", {
  testthat::local_reproducible_output(width = 80, crayon = FALSE)
  # A Gaussian posterior: no warnings and no flags, so the output ends with
  # the log marginal likelihood
  fit_ok <- inla_vanilla(lp_gauss, c(a = 0.3, b = 0.2), silent = TRUE)
  out <- capture_output(print(fit_ok))
  expect_match(out, "Log marginal likelihood: [-0-9.]+$")

  fit <- fit_ok
  fit$warnings <- c("first", "second")
  out <- capture_output(print(fit))
  expect_match(
    out,
    "! 2 warnings during the fit, kept in `$warnings`.",
    fixed = TRUE
  )
  expect_no_match(out, "flagged")

  # One line for any number of flagged checks, counted as summary() flags
  # them: the count is the number of "!" lines of summary(), for each check
  # alone and for all of them at once
  n_counted <- function(fit) {
    out <- capture_output(print(fit))
    n <- regmatches(
      out,
      regexpr("[0-9]+(?= checks? flagged by)", out, perl = TRUE)
    )
    if (length(n) == 0L) 0L else as.integer(n)
  }
  n_flag_lines <- function(fit) {
    lines <- strsplit(capture_output(print(summary(fit))), "\n")[[1]]
    length(grep("^! ", lines))
  }
  trip <- list(
    mode = function(f) {
      f$diagnostics$marginals$mode_shift_sigma[1] <- 0.5
      f
    },
    hess_cond = function(f) {
      f$diagnostics$optim$hess_cond <- 1e9
      f
    },
    vb = function(f) {
      f$diagnostics$marginals$vb_shift_sigma[2] <- -1.5
      f
    },
    nmad = function(f) {
      f$diagnostics$marginals$nmad[2] <- 0.2
      f
    },
    scan_end_mass = function(f) {
      f$diagnostics$marginals$scan_end_mass[1] <- 0.3
      f
    },
    tails = function(f) {
      f$diagnostics$marginals$tails[2] <- FALSE
      f
    },
    bend = function(f) {
      f$diagnostics$marginals$bend[1] <- 0.5
      f
    },
    sn_vs_vb = function(f) {
      f$diagnostics$marginals$sn_vs_vb[2] <- -0.3
      f
    }
  )
  # Every check that summary() flags
  expect_identical(names(trip), names(check_flags(fit_dg(fit_ok), FALSE)))
  expect_identical(n_counted(fit_ok), 0L)
  expect_identical(n_flag_lines(fit_ok), 0L)
  for (check in names(trip)) {
    fit <- trip[[check]](fit_ok)
    expect_identical(names(which(check_flags(fit_dg(fit), FALSE))), check)
    expect_identical(n_counted(fit), 1L, label = check)
    expect_identical(n_flag_lines(fit), 1L, label = check)
  }
  out <- capture_output(print(trip$bend(fit_ok)))
  expect_match(out, "! 1 check flagged by `summary()`.", fixed = TRUE)
  expect_no_match(out, "warning")
  fit <- fit_ok
  for (f in trip) {
    fit <- f(fit)
  }
  out <- capture_output(print(fit))
  expect_match(out, "! 8 checks flagged by `summary()`.", fixed = TRUE)
  expect_identical(n_flag_lines(fit), 8L)

  # The small Poisson regression: summary() flags the bend of the slope
  out <- capture_output(print(fit_reg))
  expect_match(out, "! 1 check flagged by `summary()`.", fixed = TRUE)
  expect_identical(n_flag_lines(fit_reg), 1L)

  # A heavy-tailed posterior: its warning and its flags
  lp <- function(theta) dt(theta[[1]], 1, log = TRUE)
  fit <- inla_vanilla(lp, c(theta = 0.1), silent = TRUE)
  out <- capture_output(print(fit))
  expect_match(out, "! 1 warning during the fit", fixed = TRUE)
  expect_match(out, "checks flagged by `summary()`", fixed = TRUE)
  expect_identical(n_counted(fit), n_flag_lines(fit))
})

test_that("print() rounds each row to three significant figures of its SD", {
  # Two independent Gaussians on very different scales: N(100, 0.002^2) and
  # N(1000, 300^2). Each row gets the decimals its own SD needs, so the
  # first shows five decimals and the second none.
  testthat::local_reproducible_output(width = 80, crayon = FALSE)
  lp <- function(theta) {
    dnorm(theta[[1]], 100, 0.002, log = TRUE) +
      dnorm(theta[[2]], 1000, 300, log = TRUE)
  }
  fit <- inla_vanilla(lp, c(a = 99.999, b = 500), silent = TRUE)
  expect_equal(unname(fit$summary[, "SD"]), c(0.002, 300), tolerance = 1e-6)
  out <- strsplit(capture_output(print(fit)), "\n")[[1]]
  row_a <- strsplit(trimws(grep("^a ", out, value = TRUE)), " +")[[1]]
  row_b <- strsplit(trimws(grep("^b ", out, value = TRUE)), " +")[[1]]
  expect_identical(
    row_a,
    c(
      "a",
      "100.00000",
      "0.00200",
      "99.99608",
      "100.00000",
      "100.00392",
      "100.00000"
    )
  )
  expect_identical(row_b, c("b", "1000", "300", "412", "1000", "1588", "1000"))
})

test_that("print(summary()) shows marginals, diagnostics and the evidence", {
  testthat::local_reproducible_output(width = 80, crayon = FALSE)
  expect_snapshot(print(summary(fit_reg)), transform = scrub_output)

  expect_snapshot(print(summary(fit_reg_novb)), transform = scrub_output)

  fit_grad <- inla_vanilla(
    lp_pois,
    0.3,
    grad = gr_pois,
    y = pois_y,
    a0 = pois_a0,
    b0 = pois_b0,
    silent = TRUE
  )
  out <- capture_output(print(summary(fit_grad)))
  expect_match(
    out,
    paste0(
      fit_grad$counts[["logpost"]],
      " logpost and ",
      fit_grad$counts[["grad"]],
      " grad calls"
    )
  )

  s <- summary(fit_reg)
  res <- NULL
  out <- capture_output(res <- withVisible(print(s)))
  expect_false(res$visible)
  expect_identical(res$value, s)
})

test_that("printed output is plain ASCII unless the console shows UTF-8", {
  fit_text <- function() {
    c(
      capture_output(print(fit_reg)),
      capture_output(print(summary(fit_reg))),
      capture_messages(inla_vanilla(lp_reg, init_reg)),
      capture_messages(inla_vanilla(lp_reg, init_reg, verbose = TRUE)),
      capture_messages(inla_vanilla(
        lp_reg,
        init_reg,
        grad = gr_reg,
        vb = FALSE,
        marginals = "gaussian",
        verbose = TRUE
      ))
    )
  }
  local({
    testthat::local_reproducible_output(unicode = FALSE)
    out <- fit_text()
    expect_false(any(grepl("[^\001-\177]", out, useBytes = TRUE)))
    expect_match(out, "-- INLAvanilla -+ Skew-normal \\+ VB --", all = FALSE)
    expect_match(out, "max shift [0-9.e-]+ SD", all = FALSE)
    expect_match(out, "Largest shift [0-9.e-]+ SD", all = FALSE)
  })
  local({
    testthat::local_reproducible_output(unicode = TRUE)
    out <- fit_text()
    expect_match(out, "\u2500 INLAvanilla \u2500+ Skew-normal", all = FALSE)
    expect_match(out, "max shift [0-9.e-]+\u03c3", all = FALSE)
    expect_match(out, "Largest shift [0-9.e-]+\u03c3", all = FALSE)
  })
})

test_that("summary() returns the table and the diagnostics", {
  s <- summary(fit_reg)
  expect_s3_class(s, "summary.inla_vanilla")
  expect_identical(
    colnames(s$table),
    c("Mean", "SD", "2.5%", "25%", "50%", "75%", "97.5%", "Mode")
  )
  expect_identical(rownames(s$table), names(init_reg))
  # The default columns agree with fit$summary
  shared <- colnames(fit_reg$summary)
  expect_equal(s$table[, shared], fit_reg$summary)
  expect_s3_class(s$diagnostics, "data.frame")
  expect_identical(
    names(s$diagnostics),
    c(
      "mode_shift_sigma",
      "kld",
      "vb_shift",
      "vb_shift_sigma",
      "nmad",
      "alpha",
      "scan_end_mass",
      "sn_vs_vb",
      "gamma1",
      "tails",
      "bend"
    )
  )
  expect_identical(rownames(s$diagnostics), names(init_reg))
  expect_identical(s$logml, fit_reg$logml)
  expect_identical(s$logml_laplace, fit_reg$logml_laplace)

  # The diagnostics are what they say: shifts and gaps in posterior SDs
  sds <- sqrt(diag(fit_reg$cov))
  expect_equal(s$diagnostics$vb_shift, unname(fit_reg$mean_vb - fit_reg$mode))
  expect_equal(
    s$diagnostics$vb_shift_sigma,
    unname((fit_reg$mean_vb - fit_reg$mode) / sds)
  )
  expect_equal(
    s$diagnostics$sn_vs_vb,
    unname((fit_reg$summary[, "Mean"] - fit_reg$mean_vb) / sds)
  )
  expect_equal(s$diagnostics$alpha, unname(fit_reg$sn[, "alpha"]))
})

test_that("the printed model's summary agrees with 2D quadrature", {
  # The snapshots above freeze the summary of fit_reg; check that what they
  # freeze is right, against a 201 x 201 grid over +-7 Laplace SDs. Observed:
  # skew-normal means 0.0022 and 0.0059 SDs off (Laplace 0.23 and 0.084),
  # the intercept's SD 0.37% too large, and the slope's SD 3.3% too small:
  # the small-n narrowing of the straight scan path, which the bend
  # diagnostic flags for the slope alone (0.14, against 0.037 for the
  # intercept). The log marginal likelihood is 0.017 too large.
  sds <- sqrt(diag(fit_reg$cov))
  ga <- fit_reg$mode[["intercept"]] + seq(-7, 7, length.out = 201) * sds[[1]]
  gb <- fit_reg$mode[["slope"]] + seq(-7, 7, length.out = 201) * sds[[2]]
  grid <- expand.grid(a = ga, b = gb)
  eta <- outer(grid$a, rep(1, length(reg_x))) + outer(grid$b, reg_x)
  lp <- rowSums(sweep(eta, 2, reg_y, "*") - exp(eta)) -
    sum(lgamma(reg_y + 1)) +
    dnorm(grid$a, 0, 2, log = TRUE) +
    dnorm(grid$b, 0, 2, log = TRUE)
  w <- exp(lp - max(lp))
  w <- w / sum(w)
  exact_mean <- c(sum(w * grid$a), sum(w * grid$b))
  exact_sd <- sqrt(c(
    sum(w * (grid$a - exact_mean[1])^2),
    sum(w * (grid$b - exact_mean[2])^2)
  ))
  exact_logml <- max(lp) +
    log(sum(exp(lp - max(lp))) * diff(ga[1:2]) * diff(gb[1:2]))

  err_mean <- abs(fit_reg$summary[, "Mean"] - exact_mean) / exact_sd
  expect_true(all(err_mean < 0.015))
  expect_true(all(err_mean < abs(fit_reg$mode - exact_mean) / exact_sd))
  sd_ratio <- fit_reg$summary[, "SD"] / exact_sd
  expect_lt(abs(sd_ratio[["intercept"]] - 1), 0.01)
  expect_gt(sd_ratio[["slope"]], 0.94)
  expect_lt(sd_ratio[["slope"]], 1)
  bend <- summary(fit_reg)$diagnostics$bend
  expect_identical(bend > 0.1, c(FALSE, TRUE))
  expect_lt(abs(fit_reg$logml - exact_logml), 0.035)
})

test_that("summary(probs = ) gives the requested quantile columns", {
  s <- summary(fit_reg, probs = c(0.1, 0.9))
  expect_identical(colnames(s$table), c("Mean", "SD", "10%", "90%", "Mode"))
  expect_equal(
    s$table[, c("Mean", "SD", "Mode")],
    fit_reg$summary[, c("Mean", "SD", "Mode")]
  )

  # The columns are quantiles of the fitted skew-normals: check the CDF at
  # them by integrating the skew-normal density written out here. Observed:
  # agreement to 1e-9.
  for (j in seq_len(nrow(fit_reg$sn))) {
    sn <- fit_reg$sn[j, ]
    dens <- function(x) sn_density(x, sn[["xi"]], sn[["omega"]], sn[["alpha"]])
    for (p in c(0.1, 0.9)) {
      q <- s$table[j, paste0(100 * p, "%")]
      cdf <- integrate(dens, -Inf, q, rel.tol = 1e-10)$value
      expect_equal(cdf, p, tolerance = 1e-6)
    }
  }

  # Probabilities in any order give columns in that order
  s_rev <- summary(fit_reg, probs = c(0.9, 0.1))
  expect_identical(colnames(s_rev$table), c("Mean", "SD", "90%", "10%", "Mode"))
  expect_identical(s_rev$table[, "90%"], s$table[, "90%"])
  expect_identical(s_rev$table[, "10%"], s$table[, "10%"])

  # Odd probabilities are labelled as quantile() labels them
  s <- summary(fit_reg, probs = c(0.005, 1 / 3))
  expect_identical(
    colnames(s$table),
    c("Mean", "SD", names(quantile(1, c(0.005, 1 / 3))), "Mode")
  )
  expect_output(print(s), "33.33333%", fixed = TRUE)
})

test_that("summary() rejects probabilities outside (0, 1)", {
  for (bad in list(0, c(0.5, 1), c(0.5, NA), "0.5")) {
    expect_error(
      summary(fit_reg, probs = bad),
      "`probs` must be probabilities strictly between 0 and 1",
      class = "inlavanilla_error"
    )
  }
  expect_snapshot(error = TRUE, summary(fit_reg, probs = c(0.5, 1)))
})

test_that("summary() with no probabilities gives a table without quantiles", {
  # probs = numeric(0) used to pass the check and then fail in sn_summary()
  # with base R's "length of 'dimnames' [2] not equal to array extent",
  # because format_probs(numeric(0)) returned "%".
  res <- tryCatch(summary(fit_reg, probs = numeric(0)), error = identity)
  expect_s3_class(res, "summary.inla_vanilla")
  expect_identical(colnames(res$table), c("Mean", "SD", "Mode"))
  expect_equal(res$table, fit_reg$summary[, c("Mean", "SD", "Mode")])
  expect_output(print(res), "Posterior marginals")
})

test_that("summary() flags each diagnostic at its documented level", {
  # The levels documented in ?summary.inla_vanilla, those of the fit-time
  # warning: the Newton step at the mode 0.1 SDs, the Hessian condition
  # number 1e8, the VB shift 1 SD, NMAD 0.1, the scan-end mass 0.05 and
  # tails; and, flagged by summary() alone, the bend 0.1 nats and the SN - VB
  # gap 0.1 SDs. Every diagnostic is first set 2% below its level, where
  # nothing may be flagged, and then 2% above it, where each must be, with
  # the warning's sentence where it has one.
  testthat::local_reproducible_output(width = 200, crayon = FALSE)
  below <- summary(fit_reg)
  below$optim$hess_cond <- 0.98e8
  below$diagnostics$mode_shift_sigma <- c(0.098, 0.098)
  below$diagnostics$vb_shift_sigma <- c(0.98, -0.98)
  below$diagnostics$bend <- c(0.098, 0.098)
  below$diagnostics$nmad <- c(0.098, 0.098)
  below$diagnostics$scan_end_mass <- c(0.049, 0.049)
  below$diagnostics$sn_vs_vb <- c(0.098, -0.098)
  below$diagnostics$tails <- c(TRUE, TRUE)
  below$vb$fallback <- FALSE
  # Fixed, so the lines do not depend on the optimiser's path
  below$optim$message <- "relative convergence (4)"
  below$optim$iterations <- 9L
  below$vb$iterations <- 3L
  flagged <- function(s) {
    lines <- strsplit(capture_output(print(s)), "\n")[[1]]
    grep("^! ", lines, value = TRUE)
  }

  expect_length(flagged(below), 0)
  out <- capture_output(print(below))
  for (line in c(
    "v Mode found in 9 iterations (Newton step 0.098 posterior SDs).",
    "v Hessian condition number 9.8e+07, in posterior-SD units.",
    "v VB correction in 3 iterations: max shift 0.98 posterior SDs.",
    "v Skew-normal fits: max NMAD 0.098.",
    "v Skew-normal and VB means agree to within 0.098 posterior SDs."
  )) {
    expect_match(out, line, fixed = TRUE)
  }

  above <- below
  above$optim$hess_cond <- 1.02e8
  above$diagnostics$mode_shift_sigma <- c(0.098, 0.102)
  above$diagnostics$vb_shift_sigma <- c(0.98, -1.02)
  above$diagnostics$bend <- c(0.098, 0.102)
  above$diagnostics$nmad <- c(0.098, 0.102)
  above$diagnostics$scan_end_mass <- c(0.049, 0.051)
  above$diagnostics$sn_vs_vb <- c(0.098, -0.102)
  above$diagnostics$tails <- c(TRUE, FALSE)
  expect_identical(
    flagged(above),
    c(
      "! Mode search stopped short (\"relative convergence (4)\"): the Newton step is more than 0.1 posterior SDs for `slope` (0.10).",
      "! The Hessian at the mode is near-singular (condition number = 1.02e+08); some parameters may be weakly identified.",
      "! VB correction in 3 iterations: the shift is more than 1 posterior SD for `slope` (1.02).",
      "! The fitted marginal deviates from the scanned posterior (NMAD > 0.1) for `slope` (0.10).",
      "! The fitted marginal puts more than 0.05 of its mass beyond the scanned window (4 posterior SDs either side of the mode) for `slope` (0.05); its credible limits rely on extrapolation.",
      "! Heavy tails for `slope`: the density is still over 5% of its peak at the ends of the scanned window (4 posterior SDs either side of the mode), and the skew-normal extrapolates beyond it.",
      "! The straight scan misses more than 0.1 nats at +-2 posterior SDs for `slope` (0.10): the marginal may be too narrow.",
      "! Skew-normal and VB means differ by more than 0.1 posterior SDs for `slope` (0.10)."
    )
  )

  # Every flagged parameter is listed, worst first, unlike the warning,
  # which shows at most three
  both <- below
  both$diagnostics$nmad <- c(0.12, 0.25)
  both$diagnostics$vb_shift_sigma <- c(-1.5, 1.2)
  expect_identical(
    flagged(both),
    c(
      "! VB correction in 3 iterations: the shift is more than 1 posterior SD for `intercept` (1.50), `slope` (1.20).",
      "! The fitted marginal deviates from the scanned posterior (NMAD > 0.1) for `slope` (0.25), `intercept` (0.12)."
    )
  )

  # The fallback optimiser is reported but not flagged, as in INLAvaan
  fallback <- below
  fallback$vb$fallback <- TRUE
  expect_length(flagged(fallback), 0)
  expect_match(
    capture_output(print(fallback)),
    "v VB correction in 3 iterations (via the fallback optimiser): max shift 0.98 posterior SDs.",
    fixed = TRUE
  )

  # A value that is not finite is flagged, not an error
  for (col in c("mode_shift_sigma", "nmad", "bend")) {
    odd <- below
    odd$diagnostics[[col]][1] <- NaN
    expect_length(flagged(odd), 1)
  }
  odd <- below
  odd$optim$hess_cond <- Inf
  expect_identical(
    flagged(odd),
    "! The Hessian at the mode is near-singular (condition number = Inf); some parameters may be weakly identified."
  )
})

test_that("a heavy-tailed posterior is flagged in summary()", {
  # A Cauchy posterior: 4 Laplace SDs out (2.8 scale units) the density is
  # still 11% of its peak, and no skew-normal fits it. Observed: tails FALSE
  # and NMAD 0.18.
  fit <- inla_vanilla(
    function(theta) -log1p((theta[[1]] - 1)^2),
    c(x = 0.5),
    silent = TRUE
  )
  d <- summary(fit)$diagnostics
  expect_false(d$tails)
  expect_gt(d$nmad, 0.1)
  out <- capture_output(print(summary(fit)))
  expect_match(out, "Heavy tails for `x`")
  expect_match(
    out,
    "deviates from the scanned posterior (NMAD > 0.1) for `x`",
    fixed = TRUE
  )
})

# (f) plot() -----------------------------------------------------------------

test_that("plot() draws the requested parameters and returns the fit", {
  grDevices::pdf(NULL)
  dev <- grDevices::dev.cur()
  on.exit(grDevices::dev.off(dev), add = TRUE)

  expect_no_error(res <- withVisible(plot(fit_reg)))
  expect_false(res$visible)
  expect_identical(res$value, fit_reg)
  res <- withVisible(plot(fit_reg, pars = 1, laplace = FALSE))
  expect_false(res$visible)
  expect_identical(res$value, fit_reg)

  # Each panel is titled by its parameter, and the legend names both curves
  strs <- plot_strings(record_plot(plot(fit_reg)))
  expect_true(all(c("intercept", "slope", "Skew-normal", "Laplace") %in% strs))
  for (pars in list("slope", 2)) {
    strs <- plot_strings(record_plot(plot(fit_reg, pars = pars)))
    expect_true("slope" %in% strs)
    expect_false("intercept" %in% strs)
  }
  strs <- plot_strings(record_plot(plot(fit_reg, pars = c(2, 1))))
  expect_true(all(c("intercept", "slope") %in% strs))
  strs <- plot_strings(record_plot(plot(fit_reg, laplace = FALSE)))
  expect_false(any(c("Skew-normal", "Laplace") %in% strs))

  # More parameters than fit on one row
  fit5 <- inla_vanilla(lp_gauss, c(0.1, 0.2, 0.3, 0.4, 0.5), silent = TRUE)
  strs <- plot_strings(record_plot(plot(fit5)))
  expect_true(all(paste0("theta[", 1:5, "]") %in% strs))
  strs <- plot_strings(record_plot(plot(
    fit5,
    pars = c("theta[2]", "theta[5]")
  )))
  expect_identical(
    intersect(paste0("theta[", 1:5, "]"), strs),
    c("theta[2]", "theta[5]")
  )
})

test_that("plot() draws the fitted skew-normal and Laplace densities", {
  # The filled curve of the `slope` panel is the skew-normal density of
  # fit$sn["slope", ] (written out here), the dashed curve the Gaussian
  # N(mode, cov), and the vertical line marks the skew-normal mean.
  calls <- record_plot(plot(fit_reg, pars = "slope"))
  names <- vapply(calls, `[[`, "", "name")
  sn <- fit_reg$sn["slope", ]

  poly <- calls[[which(names == "C_polygon")]]$args
  n <- length(poly[[1]]) / 2
  x <- poly[[1]][seq_len(n)]
  expect_equal(
    poly[[2]][seq_len(n)],
    sn_density(x, sn[["xi"]], sn[["omega"]], sn[["alpha"]]),
    tolerance = 1e-10
  )
  # The range covers the Laplace and skew-normal 0.05% and 99.95% points
  cdf <- function(q) {
    integrate(
      function(t) sn_density(t, sn[["xi"]], sn[["omega"]], sn[["alpha"]]),
      -Inf,
      q,
      rel.tol = 1e-10
    )$value
  }
  expect_lte(cdf(min(x)), 0.0005 + 1e-8)
  expect_gte(cdf(max(x)), 0.9995 - 1e-8)
  lap_sd <- sqrt(fit_reg$cov["slope", "slope"])
  expect_lte(min(x), qnorm(0.0005, fit_reg$mode[["slope"]], lap_sd) + 1e-8)
  expect_gte(max(x), qnorm(0.9995, fit_reg$mode[["slope"]], lap_sd) - 1e-8)

  # The lines over the same grid: the solid skew-normal curve and the
  # dashed Laplace curve. In the display list, a line's arguments are its
  # coordinates, type, pch and line type.
  lines <- Filter(
    function(a) identical(a[[2]], "l"),
    lapply(calls[names == "C_plotXY"], `[[`, "args")
  )
  solid <- Filter(function(a) identical(a[[4]], "solid"), lines)
  dashed <- Filter(function(a) !identical(a[[4]], "solid"), lines)
  expect_length(solid, 1)
  expect_length(dashed, 1)
  expect_identical(solid[[1]][[1]]$x, x)
  expect_identical(dashed[[1]][[1]]$x, x)
  expect_equal(
    solid[[1]][[1]]$y,
    sn_density(x, sn[["xi"]], sn[["omega"]], sn[["alpha"]]),
    tolerance = 1e-10
  )
  expect_equal(
    dashed[[1]][[1]]$y,
    dnorm(x, fit_reg$mode[["slope"]], lap_sd),
    tolerance = 1e-10
  )

  # The mean line, from 0 up to the curve
  seg <- unname(unlist(calls[[which(names == "C_segments")[1]]]$args[1:4]))
  mean_sn <- fit_reg$summary["slope", "Mean"]
  expect_equal(
    seg,
    c(
      mean_sn,
      0,
      mean_sn,
      sn_density(mean_sn, sn[["xi"]], sn[["omega"]], sn[["alpha"]])
    ),
    tolerance = 1e-10
  )
})

test_that("plot() rejects parameters that are not in the fit", {
  grDevices::pdf(NULL)
  dev <- grDevices::dev.cur()
  on.exit(grDevices::dev.off(dev), add = TRUE)

  expect_error(
    plot(fit_reg, pars = "sigma"),
    "Unknown parameter: \"sigma\"",
    class = "inlavanilla_error"
  )
  expect_error(
    plot(fit_reg, pars = c("sigma", "slope", "tau")),
    "Unknown parameters: \"sigma\" and \"tau\"",
    class = "inlavanilla_error"
  )
  for (bad in list(0, c(1, 3), -1, TRUE)) {
    expect_error(
      plot(fit_reg, pars = bad),
      "`pars` must name or index parameters of the fit",
      class = "inlavanilla_error"
    )
  }
  expect_snapshot(error = TRUE, plot(fit_reg, pars = c("sigma", "tau")))
})

test_that("plot() rejects a missing or empty selection of parameters", {
  grDevices::pdf(NULL)
  dev <- grDevices::dev.cur()
  on.exit(grDevices::dev.off(dev), add = TRUE)

  # These once failed with base R errors: c(1, NA) inside the check, with
  # "missing value where TRUE/FALSE needed", and an empty selection in par(),
  # with "invalid value specified for graphical parameter "mfrow"".
  for (bad in list(c(1, NA), NA_real_, integer(0), character(0))) {
    expect_package_error(plot(fit_reg, pars = bad), "`pars`")
  }
})

test_that("plot() checks laplace as inla_vanilla() checks its flags", {
  grDevices::pdf(NULL)
  dev <- grDevices::dev.cur()
  on.exit(grDevices::dev.off(dev), add = TRUE)

  # Observed: base R errors from `if (laplace)`, "missing value where
  # TRUE/FALSE needed" for NA and "argument is not interpretable as logical"
  # for "no", after the device has been set up.
  for (bad in list(NA, "no")) {
    expect_package_error(plot(fit_reg, laplace = bad), "`laplace`")
  }
})

test_that("plot() restores the graphical parameters", {
  grDevices::pdf(NULL)
  dev <- grDevices::dev.cur()
  on.exit(grDevices::dev.off(dev), add = TRUE)

  before <- graphics::par(no.readonly = TRUE)
  keep <- setdiff(names(before), par_set_by_any_plot)
  plot(fit_reg)
  expect_identical(graphics::par(no.readonly = TRUE)[keep], before[keep])
  plot(fit_reg, pars = 2, laplace = FALSE)
  expect_identical(graphics::par(no.readonly = TRUE)[keep], before[keep])
  expect_error(plot(fit_reg, pars = "sigma"))
  expect_identical(graphics::par(no.readonly = TRUE)[keep], before[keep])
  for (type in c("sn_fit", "sn_fit_log")) {
    plot(fit_reg, type = type)
    expect_identical(graphics::par(no.readonly = TRUE)[keep], before[keep])
  }
})

test_that("plot() restores graphical parameters the caller has changed", {
  grDevices::pdf(NULL)
  dev <- grDevices::dev.cur()
  on.exit(grDevices::dev.off(dev), add = TRUE)

  # Setting `mfrow` resets `cex`, and so does restoring the old `mfrow`; plot()
  # once left `cex` at 1 here instead of the caller's 0.8.
  graphics::par(mar = c(1, 2, 3, 4), las = 2, cex = 0.8, mgp = c(3, 1, 0))
  before <- graphics::par(no.readonly = TRUE)
  keep <- setdiff(names(before), par_set_by_any_plot)
  plot(fit_reg)
  after <- graphics::par(no.readonly = TRUE)
  expect_identical(after[keep], before[keep])

  # The caller's own layout of panels, which plot() replaces while it draws
  graphics::par(mfrow = c(2, 2), cex = 0.7)
  before <- graphics::par(no.readonly = TRUE)
  plot(fit_reg, pars = 1)
  expect_identical(graphics::par(no.readonly = TRUE)[keep], before[keep])
})

test_that("plot() passes axis labels, titles and other arguments to each panel", {
  strs <- plot_strings(
    record_plot(plot(fit_reg, xlab = "log scale", main = "Rate", col.axis = 2))
  )
  expect_true(all(c("log scale", "Rate") %in% strs))
  # main replaces the default title, the parameter name
  expect_false(any(c("intercept", "slope") %in% strs))
  expect_no_error(record_plot(plot(fit_reg, pars = 1, ylab = "density")))
})

test_that("plot(type = \"sn_fit\") draws each skew-normal fit against its scan", {
  # INLAvaan's sn_fit plot: on the scan coordinate, the raw and the
  # corrected scan, each scaled to a peak of 1, and the fitted skew-normal
  # at the height it was fitted at, from the scans kept in the fit
  calls <- record_plot(
    res <- withVisible(plot(fit_reg, pars = "slope", type = "sn_fit"))
  )
  expect_false(res$visible)
  expect_identical(res$value, fit_reg)
  strs <- plot_strings(calls)
  expect_true(all(c("slope", "Raw", "Corrected", "Skew-normal fit") %in% strs))
  expect_false("intercept" %in% strs)
  curves <- function(calls) {
    names <- vapply(calls, `[[`, "", "name")
    lines <- lapply(calls[names == "C_plotXY"], `[[`, "args")
    lapply(Filter(function(a) identical(a[[2]], "l"), lines), `[[`, 1)
  }
  lines <- curves(calls)
  expect_length(lines, 3)
  z <- fit_reg$scan$z
  raw <- fit_reg$scan$logpost[, "slope"]
  y <- raw + fit_reg$diagnostics$marginals$gamma1[2] * z
  expect_identical(lines[[1]]$x, z)
  expect_equal(lines[[1]]$y, unname(exp(raw - max(raw))))
  expect_equal(lines[[2]]$y, unname(exp(y - max(y))))
  # The fit, drawn on 201 points, of which every tenth is a scan point: the
  # largest gap there is the NMAD of the fit
  expect_length(lines[[3]]$x, 201)
  at_z <- seq(1, 201, by = 10)
  expect_equal(lines[[3]]$x[at_z], z)
  expect_equal(
    max(abs(exp(y - max(y)) - lines[[3]]$y[at_z])),
    fit_reg$diagnostics$marginals$nmad[2],
    tolerance = 1e-8
  )

  # The same in base-10 logs
  lines <- curves(record_plot(plot(fit_reg, pars = 2, type = "sn_fit_log")))
  expect_equal(lines[[2]]$y, unname(log10(exp(y - max(y)))))

  # Both parameters, and the type checked as inla_vanilla() checks its
  # choices
  strs <- plot_strings(record_plot(plot(fit_reg, type = "sn_fit")))
  expect_true(all(c("intercept", "slope") %in% strs))
  expect_error(
    plot(fit_reg, type = "density"),
    "`type` must be one of \"marg_pdf\", \"sn_fit\", or \"sn_fit_log\"",
    class = "inlavanilla_error"
  )
})

test_that("scan_end_mass is the fitted mass outside the scanned window", {
  # The fitted skew-normal integrated outside mode +- 4 Laplace SDs, with the
  # density written out independently of the package
  d <- fit_reg$diagnostics$marginals
  sds <- sqrt(diag(fit_reg$cov))
  for (j in 1:2) {
    f <- function(x) {
      sn_density(
        x,
        fit_reg$sn[j, "xi"],
        fit_reg$sn[j, "omega"],
        fit_reg$sn[j, "alpha"]
      )
    }
    lo <- fit_reg$mode[[j]] - 4 * sds[[j]]
    hi <- fit_reg$mode[[j]] + 4 * sds[[j]]
    ref <- integrate(f, -Inf, lo, rel.tol = 1e-10)$value +
      integrate(f, hi, Inf, rel.tol = 1e-10)$value
    expect_equal(d$scan_end_mass[j], ref, tolerance = 1e-6)
  }
  expect_true(all(d$scan_end_mass > 0 & d$scan_end_mass < 0.05))

  # A Gaussian posterior: the fit is Gaussian, so the mass is 2 pnorm(-4)
  lp <- function(theta) sum(dnorm(theta, c(1, -2), c(0.5, 3), log = TRUE))
  fit <- inla_vanilla(lp, c(a = 0, b = 0), silent = TRUE)
  expect_equal(
    fit$diagnostics$marginals$scan_end_mass,
    rep(2 * pnorm(-4), 2),
    tolerance = 1e-3
  )

  # Shown in the table, and not flagged
  out <- capture_output(print(summary(fit_reg)))
  expect_match(out, "End mass")
  expect_no_match(out, "beyond the scanned window")
})

test_that("a gradient wrong by a factor is caught when init is at the mode", {
  # At the mode both gradients are near zero, so the check at init alone
  # passes a gradient that forgot a factor of 1 / sigma^2; the second check,
  # one curvature scale away, catches it.
  set.seed(1)
  n <- 5000
  x <- rnorm(n)
  y <- 1 + 2 * x + rnorm(n, 0, 2)
  X <- cbind(1, x)
  lp <- function(b) {
    sum(dnorm(y, drop(X %*% b), 2, log = TRUE)) +
      sum(dnorm(b, 0, 10, log = TRUE))
  }
  gr_ok <- function(b) drop(crossprod(X, y - X %*% b)) / 4 - b / 100
  gr_bad <- function(b) 4 * gr_ok(b)
  init <- stats::setNames(coef(lm(y ~ x)), c("a", "b"))
  fit <- inla_vanilla(lp, init, grad = gr_bad, silent = TRUE)
  expect_match(
    fit$warnings,
    "does not match finite differences of `logpost` one curvature scale",
    all = FALSE
  )
  fit <- inla_vanilla(lp, init, grad = gr_ok, silent = TRUE)
  expect_identical(fit$warnings, character())
})

test_that("silent = TRUE copes with conditions that have no muffle restart", {
  # signalCondition() sets up no muffle restart. The handlers once called
  # invokeRestart(), which made such a condition an error that stopped the
  # fit. No handler can muffle it, so it carries on to testthat's, which
  # leaves warnings alone under a negative `warn`; at top level, where
  # nothing handles it, nothing prints.
  lp <- function(theta) {
    signalCondition(simpleMessage("note\n"))
    signalCondition(simpleWarning("bare warning"))
    -0.5 * sum(theta^2)
  }
  old <- options(warn = -1)
  on.exit(options(old), add = TRUE)
  fit <- inla_vanilla(lp, c(a = 0.3), silent = TRUE)
  options(old)
  expect_s3_class(fit, "inla_vanilla")
  expect_identical(fit$warnings, "bare warning")
})
