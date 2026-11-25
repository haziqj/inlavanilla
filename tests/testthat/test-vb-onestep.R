# The one-step VB correction, the default of inla_vanilla(): the `vb`
# argument, the axis rule, the cost of the step and its accuracy against the
# true posterior mean. With a gradient, the check of the step from a second
# sweep, and the solve and the warning when it fails. Without one, the gain
# ratio and the guard that hands a large first step to the iteration.

# The log-density of theta = log(lambda) for lambda ~ Gamma(a, 1). Mode log(a),
# Laplace variance 1 / a, true mean digamma(a).
log_gamma_model <- function(a) {
  force(a)
  list(
    logpost = function(theta) a * theta[[1]] - exp(theta[[1]]),
    grad = function(theta) a - exp(theta[[1]]),
    init = c(eta = log(a) + 0.1)
  )
}

# The same with no density below `cut`, where the gradient is not finite
cut_log_gamma_model <- function(a, cut) {
  force(a)
  force(cut)
  list(
    logpost = function(theta) {
      if (theta[[1]] < cut) -Inf else a * theta[[1]] - exp(theta[[1]])
    },
    grad = function(theta) {
      if (theta[[1]] < cut) NaN else a - exp(theta[[1]])
    },
    init = c(eta = log(a) + 0.1)
  )
}

# The three-point rule's first step on that model, in Laplace SDs
log_gamma_first_step <- function(a) {
  s <- 1 / sqrt(a)
  sqrt(a) * (1 - (2 / 3 + (exp(s * sqrt(3)) + exp(-s * sqrt(3))) / 6))
}

# A logistic regression with a normal prior on its coefficients
logistic_model <- function(n, m, seed) {
  set.seed(seed)
  X <- cbind(1, matrix(stats::rnorm(n * (m - 1)), n))
  beta <- c(-0.5, stats::rnorm(m - 1, 0, 0.5))
  y <- stats::rbinom(n, 1, stats::plogis(X %*% beta))
  list(
    logpost = function(theta) {
      eta <- as.numeric(X %*% theta)
      sum(y * eta - log1p(exp(eta))) +
        sum(stats::dnorm(theta, 0, 2.5, log = TRUE))
    },
    grad = function(theta) {
      as.numeric(crossprod(X, y - stats::plogis(as.numeric(X %*% theta)))) -
        theta / 6.25
    },
    init = stats::setNames(rep(0, m), paste0("b", seq_len(m)))
  )
}

skew_normal_model <- function(alpha) {
  force(alpha)
  list(
    logpost = function(theta) {
      stats::dnorm(theta[[1]], log = TRUE) +
        stats::pnorm(alpha * theta[[1]], log.p = TRUE)
    },
    grad = function(theta) {
      x <- theta[[1]]
      -x +
        alpha *
          exp(
            stats::dnorm(alpha * x, log = TRUE) -
              stats::pnorm(alpha * x, log.p = TRUE)
          )
    },
    init = c(x = 0.3)
  )
}

fit_model <- function(model, ...) {
  inla_vanilla(model$logpost, model$init, grad = model$grad, ..., silent = TRUE)
}

fit_no_grad <- function(model, ...) {
  inla_vanilla(model$logpost, model$init, ..., silent = TRUE)
}

# Independent log-gamma terms with shapes `a` on mixed coordinates,
# eta = M theta, so that the parameters are correlated. A small shape makes
# the posterior far from Gaussian along one direction.
mixed_gamma_model <- function(a, names = paste0("t", seq_along(a))) {
  m <- length(a)
  M <- diag(m) + 0.3
  list(
    logpost = function(theta) {
      eta <- as.numeric(M %*% theta)
      sum(a * eta - exp(eta))
    },
    grad = function(theta) {
      eta <- as.numeric(M %*% theta)
      as.numeric(crossprod(M, a - exp(eta)))
    },
    init = stats::setNames(solve(M, log(a)) + 0.1, names)
  )
}

# x ~ Gamma(shape, 1) for x = theta + 1, with no density at or below `edge`,
# where the gradient is not finite. Mode shape - 1, Laplace variance shape - 1.
edge_gamma_model <- function(shape, edge = 0) {
  force(shape)
  force(edge)
  list(
    logpost = function(theta) {
      x <- theta[[1]] + 1
      if (x <= edge) -Inf else (shape - 1) * log(x) - x
    },
    grad = function(theta) {
      x <- theta[[1]] + 1
      if (x <= edge) NaN else (shape - 1) / x - 1
    },
    init = c(a = shape - 1.5)
  )
}

# The Laplace stage of a model and its one-step VB correction, called
# directly, with the gradient calls of the correction counted
vb_direct <- function(mod, budget = Inf, fixed = integer()) {
  ctrl <- vanilla_control()
  lap <- laplace_fit(mod$logpost, mod$grad, mod$init, ctrl)
  calls <- 0L
  gr <- function(theta) {
    calls <<- calls + 1L
    mod$grad(theta)
  }
  vb <- vb_correct(
    mod$logpost,
    gr,
    lap$mode,
    lap$Sigma,
    lap$R,
    lap$A,
    ctrl,
    method = "onestep",
    budget = budget,
    fixed = fixed
  )
  list(vb = vb, lap = lap, calls = calls)
}

# The one step as it was before the check existed: Sigma times the expected
# gradient on axis_nodes(), with the step of the fixed parameters at zero
plain_onestep <- function(mod, lap, fixed = integer()) {
  rule <- axis_nodes(lap$Sigma, lap$A)
  g0 <- mod$grad(lap$mode)
  s <- g0
  for (b in 2:nrow(rule$nodes)) {
    s <- s + rule$weights[b] * (mod$grad(lap$mode + rule$nodes[b, ]) - g0)
  }
  step <- as.numeric(lap$Sigma %*% s)
  step[fixed] <- 0
  stats::setNames(step, names(lap$mode))
}

# The whitened gradient of log_gamma_model(a), and the residual of the VB
# equation on the three-point rule at a whitened shift d
log_gamma_slope <- function(a, z) sqrt(a) * (1 - exp(z / sqrt(a)))
log_gamma_resid <- function(a, d) {
  2 /
    3 *
    log_gamma_slope(a, d) +
    (log_gamma_slope(a, d + sqrt(3)) + log_gamma_slope(a, d - sqrt(3))) / 6
}

# The warning of a flagged one step, in plain text
check_warning <- function(name) {
  paste0(
    "VB one-step check flagged that the posterior is far from Gaussian along ",
    "the mean correction, most of all in `",
    name,
    "`.\n",
    "i Treat posterior means, mloglik, looic and waic with care.\n",
    "i This reflects the model, the data, or the prior, not a numerical ",
    "failure. If these numbers matter, compare with MCMC."
  )
}

mod2 <- log_gamma_model(2)

# The argument --------------------------------------------------------------

test_that("vb_correction takes a method, an abbreviation, TRUE or FALSE", {
  vb_of <- function(...) fit_model(mod2, ...)$vb_correction
  expect_identical(vb_of(), "onestep")
  expect_identical(vb_of(vb_correction = "onestep"), "onestep")
  expect_identical(vb_of(vb_correction = "full"), "full")
  expect_identical(vb_of(vb_correction = "none"), "none")
  expect_identical(vb_of(vb_correction = TRUE), "onestep")
  expect_identical(vb_of(vb_correction = FALSE), "none")
  expect_identical(vb_of(vb_correction = "o"), "onestep")
  expect_identical(vb_of(vb_correction = "FULL"), "full")
  expect_identical(vb_of(vb_correction = "n"), "none")
  for (bad in list("off", "fast", NA, c(TRUE, FALSE), 1)) {
    expect_error(
      fit_model(mod2, vb_correction = bad),
      "`vb_correction` must be one of",
      class = "inlavanilla_error"
    )
  }
})

test_that("vb, the old name of vb_correction, still works, with a warning", {
  rm(list = ls(renamed), envir = renamed)
  # Not silent, which would muffle the warning; the stage lines are messages
  fit_loud <- function(...) {
    suppressMessages(
      inla_vanilla(mod2$logpost, mod2$init, grad = mod2$grad, ...)
    )
  }
  new <- fit_model(mod2, vb_correction = "none")
  expect_warning(
    old <- fit_loud(vb = "off"),
    "`vb` is deprecated",
    class = "inlavanilla_deprecated"
  )
  expect_identical(old$vb_correction, "none")
  expect_identical(old$summary, new$summary)
  # Once per session
  expect_no_warning(fit_loud(vb = FALSE))
  for (pair in list(
    list(TRUE, "onestep"),
    list(FALSE, "none"),
    list("full", "full"),
    list("one", "onestep")
  )) {
    expect_identical(fit_model(mod2, vb = pair[[1]])$vb_correction, pair[[2]])
  }
  err <- expect_error(
    fit_model(mod2, vb = "none"),
    "`vb` must be one of",
    class = "inlavanilla_error"
  )
  expect_match(
    cli::ansi_strip(conditionMessage(err)),
    "`vb` is now `vb_correction`, which takes \"onestep\""
  )
  expect_error(
    fit_model(mod2, vb = "off", vb_correction = "none"),
    "Give `vb_correction` alone",
    class = "inlavanilla_error"
  )
  # A silent fit records the warning, and leaves the next fit to show it
  rm(list = ls(renamed), envir = renamed)
  fit <- fit_model(mod2, vb = "off")
  expect_match(fit$warnings, "`vb` is deprecated", all = FALSE)
  expect_warning(fit_loud(vb = "off"), class = "inlavanilla_deprecated")
  rm(list = ls(renamed), envir = renamed)
})

test_that("vb_correction = \"none\" is the Laplace centre, and FALSE is the same fit", {
  off <- fit_model(mod2, vb_correction = "none")
  no <- fit_model(mod2, vb_correction = FALSE)
  expect_identical(off$mean_vb, off$mode)
  expect_identical(off$logml, off$logml_laplace)
  expect_null(off$diagnostics$vb)
  expect_identical(no$mean_vb, off$mean_vb)
  expect_identical(no$logml, off$logml)
})

# The rule ------------------------------------------------------------------

test_that("axis_nodes() has 2m + 1 antipodal nodes at sqrt(3) SDs", {
  Sigma <- matrix(c(4, 1, 0.5, 1, 2, 0.2, 0.5, 0.2, 1), 3)
  rule <- axis_nodes(Sigma)
  expect_identical(nrow(rule$nodes), 7L)
  expect_equal(unname(rule$nodes[1, ]), c(0, 0, 0))
  expect_equal(rule$nodes[2:4, ], -rule$nodes[5:7, ])
  expect_equal(sum(rule$weights), 1)
  # Whitened, every node is sqrt(3) from the centre
  white <- rule$nodes %*% t(chol(solve(Sigma)))
  expect_equal(unname(sqrt(rowSums(white[-1, ]^2))), rep(sqrt(3), 6))
})

test_that("axis_nodes() integrates Gaussian moments up to degree three", {
  Sigma <- matrix(c(4, 1, 1, 2), 2)
  rule <- axis_nodes(Sigma)
  x <- rule$nodes
  w <- rule$weights
  expect_equal(sum(w * x[, 1]), 0)
  expect_equal(sum(w * x[, 1] * x[, 2]), Sigma[1, 2])
  expect_equal(sum(w * x[, 1]^2), Sigma[1, 1])
  expect_equal(sum(w * x[, 2]^2), Sigma[2, 2])
  expect_equal(sum(w * x[, 1]^2 * x[, 2]), 0)
  expect_equal(sum(w * x[, 1]^3), 0)
})

test_that("in one dimension the axis rule is the three-point rule of vb_nodes()", {
  Sigma <- matrix(2.5^2)
  expect_equal(sort(axis_nodes(Sigma)$nodes), sort(vb_nodes(Sigma)$nodes))
  expect_equal(sort(axis_nodes(Sigma)$weights), sort(vb_nodes(Sigma)$weights))
})

# One step ------------------------------------------------------------------

test_that("one step is the first step of the rule, with no iteration", {
  for (a in c(1, 2, 5, 20)) {
    fit <- fit_model(log_gamma_model(a))
    shift <- (fit$mean_vb[[1]] - fit$mode[[1]]) / sqrt(fit$cov[1, 1])
    expect_equal(shift, log_gamma_first_step(a), tolerance = 1e-5)
    expect_identical(fit$diagnostics$vb$iterations, 1L)
    expect_true(fit$diagnostics$vb$one_pass)
    expect_false(fit$diagnostics$vb$solve)
    expect_false(fit$diagnostics$vb$fallback)
  }
})

test_that("with a gradient one step costs 4m + 2 gradient calls", {
  mod <- logistic_model(80, 6, 1)
  on <- fit_model(mod, marginals = "gaussian")
  off <- fit_model(mod, vb_correction = "none", marginals = "gaussian")
  # The sweep of the step and the sweep of its check, and no calls of logpost
  expect_identical(on$counts[["grad"]] - off$counts[["grad"]], 26L)
  expect_identical(on$counts[["logpost"]], off$counts[["logpost"]])
})

test_that("without a gradient one step is one pass of logpost calls", {
  mod <- logistic_model(80, 4, 2)
  on <- fit_no_grad(mod, marginals = "gaussian")
  off <- fit_no_grad(mod, vb_correction = "none", marginals = "gaussian")
  n_nodes <- nrow(vb_nodes(on$cov)$nodes) - 1L
  # The pass, and 4m + 2 calls for the gain ratio
  expect_identical(
    on$counts[["logpost"]] - off$counts[["logpost"]],
    n_nodes + 18L
  )
  expect_identical(on$diagnostics$vb$iterations, 1L)
  # The check of the step needs a gradient
  expect_true(is.na(on$diagnostics$vb$solve))
  expect_true(all(is.na(on$diagnostics$vb$check)))
})

test_that("one step has no settling check, which is not a flag", {
  fit <- fit_model(mod2)
  expect_true(all(is.na(fit$diagnostics$marginals$vb_step_sigma)))
  expect_identical(fit$warnings, character())
  flags <- check_flags(
    diagnostics_internal(
      fit$diagnostics$optim,
      fit$diagnostics$vb,
      fit$diagnostics$marginals
    ),
    gaussian = FALSE
  )
  expect_false(flags[["vb_step"]])
  expect_false(flags[["vb_check"]])
  expect_match(capture_output(print(fit)), "Converged in")
})

test_that("the log marginal likelihood adds the step's KL term", {
  fit <- fit_model(mod2)
  d <- (fit$mean_vb - fit$mode) / sqrt(diag(fit$cov))
  expect_equal(fit$logml - fit$logml_laplace, 0.5 * sum(d^2), tolerance = 1e-6)
})

test_that("the labels name the VB method", {
  testthat::local_reproducible_output(width = 80, crayon = FALSE)
  expect_match(
    capture_output(print(fit_model(mod2))),
    "Skew-normal + one-step VB --",
    fixed = TRUE
  )
  expect_match(
    capture_output(print(fit_model(mod2, vb_correction = "full"))),
    "Skew-normal + full VB --",
    fixed = TRUE
  )
  expect_match(
    capture_output(print(summary(fit_model(mod2)))),
    "VB correction in one step: max shift",
    fixed = TRUE
  )
  expect_match(
    capture_output(print(summary(fit_model(mod2, vb_correction = "full")))),
    "VB correction settled in",
    fixed = TRUE
  )
})

# Accuracy ------------------------------------------------------------------

test_that("one step is about as close to the true mean as the full solve", {
  for (a in c(1, 2, 5, 20)) {
    mod <- log_gamma_model(a)
    on <- fit_model(mod)
    full <- fit_model(mod, vb_correction = "full")
    s <- sqrt(on$cov[1, 1])
    truth <- digamma(a)
    # Observed 0.0005 to 0.061 SDs for one step, against 0.0009 to 0.084
    expect_lt(abs(on$mean_vb[[1]] - truth) / s, 0.07)
    expect_lt(
      abs(on$mean_vb[[1]] - truth) / s,
      abs(full$mean_vb[[1]] - truth) / s + 0.005
    )
  }
})

test_that("one step and the full solve agree to the order the theory gives", {
  # The gap shrinks like a^-3/2, in SDs, between step 1 and exact VB
  gap <- function(a) {
    mod <- log_gamma_model(a)
    on <- fit_model(mod)
    full <- fit_model(mod, vb_correction = "full")
    abs(on$mean_vb[[1]] - full$mean_vb[[1]]) / sqrt(on$cov[1, 1])
  }
  expect_lt(gap(50) / gap(5), 3 * (5 / 50)^1.5)
})

test_that("in several dimensions one step moves most of the way", {
  mod <- logistic_model(150, 8, 3)
  on <- fit_model(mod)
  full <- fit_model(mod, vb_correction = "full")
  sds <- sqrt(diag(on$cov))
  shift_on <- (on$mean_vb - on$mode) / sds
  shift_full <- (full$mean_vb - full$mode) / sds
  # Observed: the largest shift is 0.33 SDs, and the two differ by 0.03
  expect_lt(max(abs(shift_on - shift_full)), 0.05)
  expect_gt(stats::cor(shift_on, shift_full), 0.99)
})

# The check of the step, with a gradient -------------------------------------

test_that("a step that passes the check is the plain one step, to the bit", {
  for (mod in list(
    log_gamma_model(2),
    logistic_model(80, 6, 1),
    mixed_gamma_model(c(2, 5, 20))
  )) {
    out <- vb_direct(mod)
    plain <- plain_onestep(mod, out$lap)
    expect_false(out$vb$solve)
    expect_true(out$vb$one_pass)
    expect_identical(out$vb$shift, plain)
    expect_identical(
      out$vb$kld_global,
      0.5 * sum((out$lap$R %*% plain)^2)
    )
    expect_true(all(is.na(out$vb$step_sigma)))
    expect_identical(out$calls, 2L * (2L * length(plain) + 1L))
  }
})

test_that("the check is the contraction ratio and the slope at both ends", {
  # On the three-point rule the whitened residual and its slope are known
  for (a in c(0.4, 2)) {
    fit <- fit_model(log_gamma_model(a))
    d1 <- log_gamma_resid(a, 0)
    slope <- function(d) {
      1 +
        (log_gamma_slope(a, d + sqrt(3)) - log_gamma_slope(a, d - sqrt(3))) /
          (2 * sqrt(3))
    }
    expect_equal(
      fit$diagnostics$vb$check,
      c(
        lhat = abs(log_gamma_resid(a, d1)) / abs(d1),
        j0 = slope(0),
        j1 = slope(d1)
      ),
      tolerance = 1e-5
    )
    # The slope at the mode in closed form
    x <- sqrt(3 / a)
    expect_equal(
      fit$diagnostics$vb$check[["j0"]],
      1 - sinh(x) / x,
      tolerance = 1e-5
    )
  }
})

test_that("the check passes where the map contracts, and fails where not", {
  # The slope at the mode of the log-gamma crosses -1 at a = 0.633
  solved <- function(a) fit_model(log_gamma_model(a))$diagnostics$vb$solve
  expect_false(solved(2))
  expect_false(solved(0.65))
  expect_true(solved(0.62))
  expect_true(solved(0.4))
})

test_that("a failed check solves the VB equation of the rule", {
  a <- 0.4
  mod <- log_gamma_model(a)
  fit <- fit_model(mod)
  vb <- fit$diagnostics$vb
  expect_true(vb$solve)
  expect_false(vb$one_pass)
  expect_false(vb$fallback)
  expect_gt(vb$iterations, 1L)
  # The root of the three-point rule's equation, which is not the one step
  root <- stats::uniroot(
    function(d) log_gamma_resid(a, d),
    c(-3, 0),
    tol = 1e-12
  )$root
  s <- sqrt(fit$cov[1, 1])
  shift <- (fit$mean_vb[[1]] - fit$mode[[1]]) / s
  expect_equal(shift, root, tolerance = 2e-3)
  expect_gt(abs(log_gamma_first_step(a) - root), 0.5)
  # In one dimension the full solve uses the same rule
  full <- fit_model(mod, vb_correction = "full")
  expect_lt(abs(fit$mean_vb[[1]] - full$mean_vb[[1]]) / s, 2e-3)
  # The step still to go is recorded, below vb_tol where the solve converged
  expect_lt(fit$diagnostics$marginals$vb_step_sigma, 1e-3)
  # The evidence takes the KL term at the landing
  expect_equal(fit$logml - fit$logml_laplace, 0.5 * shift^2, tolerance = 1e-10)
})

test_that("a flagged step warns once, in the agreed words", {
  testthat::local_reproducible_output(width = 80, crayon = FALSE)
  mod <- log_gamma_model(0.4)
  fit <- fit_model(mod)
  expect_identical(
    fit$warnings[startsWith(fit$warnings, "VB one-step check")],
    check_warning("eta")
  )
  # Raised as a warning of its own, apart from the diagnostics warning
  seen <- list()
  withCallingHandlers(
    suppressMessages(inla_vanilla(mod$logpost, mod$init, grad = mod$grad)),
    warning = function(w) {
      seen[[length(seen) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )
  is_check <- vapply(
    seen,
    function(w) grepl("VB one-step check flagged", conditionMessage(w)),
    NA
  )
  expect_identical(sum(is_check), 1L)
  expect_s3_class(seen[[which(is_check)]], "inlavanilla_warning")
  expect_false(inherits(
    seen[[which(is_check)]],
    "inlavanilla_diagnostics_warning"
  ))
  # A step that passes raises nothing
  expect_identical(fit_model(log_gamma_model(2))$warnings, character())
})

test_that("the warning names the parameter the one step moves furthest", {
  mod <- mixed_gamma_model(c(20, 0.3, 5), names = c("tight", "wide", "mid"))
  out <- vb_direct(mod)
  expect_true(out$vb$solve)
  plain <- abs(plain_onestep(mod, out$lap)) / sqrt(diag(out$lap$Sigma))
  expect_identical(out$vb$check_worst, names(plain)[which.max(plain)])
  fit <- fit_model(mod)
  expect_identical(
    fit$warnings[startsWith(fit$warnings, "VB one-step check")],
    check_warning(out$vb$check_worst)
  )
})

test_that("the solve never takes the VB calls past the budget", {
  for (mod in list(
    log_gamma_model(0.4),
    mixed_gamma_model(c(20, 0.3, 5))
  )) {
    sweep <- 2L * length(mod$init) + 1L
    togo <- function(out) max(out$vb$step_sigma)
    last <- Inf
    for (budget in sweep * c(0, 2, 2.9, 3, 4, 6, 10)) {
      out <- vb_direct(mod, budget = budget)
      expect_true(out$vb$solve)
      # The two sweeps of the check are always taken
      expect_lte(out$calls, max(budget, 2L * sweep))
      expect_identical(out$calls, out$vb$history$calls)
      expect_true(all(is.finite(out$vb$shift)))
      expect_lte(togo(out), last + 1e-12)
      last <- togo(out)
    }
    # With no room for a sweep of the solve the mode stands, and the step
    # still to go is the one step
    none <- vb_direct(mod, budget = 2.9 * sweep)
    plain <- plain_onestep(mod, none$lap)
    expect_identical(none$calls, 2L * sweep)
    expect_identical(unname(none$vb$shift), numeric(length(plain)))
    expect_identical(none$vb$iterations, 0L)
    expect_identical(none$vb$kld_global, 0)
    expect_equal(
      none$vb$step_sigma,
      abs(plain) / sqrt(diag(none$lap$Sigma))
    )
  }
})

test_that("a whole fit keeps the VB stage within the mode and Hessian stage", {
  mod <- mixed_gamma_model(c(20, 0.3, 5))
  on <- fit_model(mod, marginals = "gaussian")
  off <- fit_model(mod, vb_correction = "none", marginals = "gaussian")
  expect_true(on$diagnostics$vb$solve)
  expect_identical(on$counts[["logpost"]], off$counts[["logpost"]])
  expect_lte(
    on$counts[["grad"]] - off$counts[["grad"]],
    off$counts[["logpost"]] + off$counts[["grad"]]
  )
})

test_that("the checked step does not depend on the order of the parameters", {
  ctrl <- vanilla_control()
  perm <- c(3L, 1L, 2L)
  back <- order(perm)
  # One model that passes the check and one that is solved
  for (a in list(c(2, 5, 20), c(20, 0.3, 5))) {
    mod <- mixed_gamma_model(a)
    out <- vb_direct(mod, budget = 70)
    lap <- out$lap
    Sigma <- lap$Sigma[perm, perm]
    moved <- vb_correct(
      function(theta) mod$logpost(theta[back]),
      function(theta) mod$grad(theta[back])[perm],
      lap$mode[perm],
      Sigma,
      chol(crossprod(lap$R)[perm, perm]),
      cov_root(Sigma),
      ctrl,
      method = "onestep",
      budget = 70
    )
    expect_identical(moved$solve, out$vb$solve)
    expect_identical(moved$iterations, out$vb$iterations)
    expect_identical(moved$history$sweeps, out$vb$history$sweeps)
    expect_identical(moved$check_worst, out$vb$check_worst)
    expect_equal(moved$check, out$vb$check, tolerance = 1e-8)
    expect_equal(moved$shift, out$vb$shift[perm], tolerance = 1e-8)
    expect_equal(moved$step_sigma, out$vb$step_sigma[perm], tolerance = 1e-6)
  }
})

test_that("the shift of a fixed parameter stays at zero", {
  # A step that passes: the plain one step with that entry at zero
  mod <- mixed_gamma_model(c(2, 5, 20))
  out <- vb_direct(mod, fixed = 2L)
  expect_false(out$vb$solve)
  expect_identical(out$vb$shift, plain_onestep(mod, out$lap, fixed = 2L))
  expect_identical(out$vb$shift[[2]], 0)
  # A step that is solved
  mod <- mixed_gamma_model(c(20, 0.15, 5))
  out <- vb_direct(mod, fixed = 3L)
  expect_true(out$vb$solve)
  expect_identical(out$vb$shift[[3]], 0)
  expect_identical(out$vb$step_sigma[[3]], 0)
  expect_true(all(out$vb$shift[-3] != 0))
})

test_that("a gradient not finite at a node of the second sweep is handled", {
  # Shape 0.7 passes the check when the support is the whole line. Here the
  # lowest node of the sweep at the mode is inside the support, and the one
  # step moves the mean down, so the lowest node of the second sweep is not.
  cut <- -3.3
  mod <- cut_log_gamma_model(0.7, cut)
  expect_false(fit_model(log_gamma_model(0.7))$diagnostics$vb$solve)
  out <- vb_direct(mod, budget = 60)
  vb <- out$vb
  low <- function(shift) {
    out$lap$mode[[1]] + shift - sqrt(3 * out$lap$Sigma[1, 1])
  }
  plain <- plain_onestep(mod, out$lap)
  expect_gt(low(0), cut)
  expect_true(is.finite(plain))
  expect_lt(low(plain[[1]]), cut)
  # The residual at the step is not finite, so the check fails: a value that
  # is not finite, and not the NA of a check that did not run
  expect_true(is.nan(vb$check[["lhat"]]))
  expect_true(is.nan(vb$check[["j1"]]))
  expect_true(is.finite(vb$check[["j0"]]))
  expect_true(vb$solve)
  # The solve lands where every node is finite, at the root of the rule
  expect_true(is.finite(vb$shift))
  expect_gt(low(vb$shift[[1]]), cut)
  expect_lt(vb$step_sigma[[1]], 1e-3)
  root <- stats::uniroot(
    function(d) log_gamma_resid(0.7, d),
    c(-3, 0),
    tol = 1e-12
  )$root
  expect_equal(
    vb$shift[[1]] / sqrt(out$lap$Sigma[1, 1]),
    root,
    tolerance = 2e-3
  )
  expect_lte(out$calls, 60L)
  # And a whole fit warns as for any flagged step
  testthat::local_reproducible_output(width = 80, crayon = FALSE)
  fit <- fit_model(mod)
  expect_true(fit$diagnostics$vb$solve)
  expect_true(check_warning("eta") %in% fit$warnings)
})

test_that("a node that touches the edge of the support is handled", {
  # The lowest node at the mode is at x = 0.025, where the gradient is finite
  # and large: the one step leaves the scale of the posterior
  mod <- edge_gamma_model(4.05)
  out <- vb_direct(mod, budget = 60)
  plain <- plain_onestep(mod, out$lap)
  sd <- sqrt(out$lap$Sigma[1, 1])
  expect_gt(abs(plain[[1]]) / sd, 5)
  expect_lt(out$vb$check[["j0"]], -1)
  expect_true(out$vb$solve)
  # Gamma(4.05, 1) has its mean 1 above its mode, 0.57 SDs
  expect_gt(out$vb$shift[[1]] / sd, 0)
  expect_lt(out$vb$shift[[1]] / sd, 1)
})

test_that("a sweep at the mode that is not finite goes to the fallback", {
  # a + 1 ~ Gamma(4, 1), with no density at or below 0.05: the lowest node at
  # the mode is at x = 0, so there is no step to check
  mod <- edge_gamma_model(4, edge = 0.05)
  on <- fit_model(mod)
  full <- fit_model(mod, vb_correction = "full")
  vb <- on$diagnostics$vb
  expect_true(vb$fallback)
  expect_false(vb$one_pass)
  expect_true(is.na(vb$solve))
  expect_true(all(is.na(vb$check)))
  sds <- sqrt(diag(on$cov))
  expect_lt(abs(on$mean_vb[[1]] - full$mean_vb[[1]]) / sds[[1]], 0.02)
  expect_false(any(startsWith(on$warnings, "VB one-step check")))
})

test_that("a step too small to check is taken with one sweep", {
  # A Gaussian posterior: the first step is zero
  lp <- function(theta) -0.5 * sum(theta^2)
  gr <- function(theta) -theta
  mod <- list(logpost = lp, grad = gr, init = c(a = 0.5, b = -0.2))
  out <- vb_direct(mod)
  expect_identical(out$calls, 5L)
  expect_true(out$vb$one_pass)
  expect_true(is.na(out$vb$solve))
  expect_true(all(is.na(out$vb$check)))
  fit <- fit_model(mod)
  expect_identical(fit$warnings, character())
  expect_true(all(is.na(
    diagnostics(fit)[c("vb_lhat", "vb_j0", "vb_j1", "vb_solve")]
  )))
})

test_that("diagnostics(), print() and summary() report the check", {
  testthat::local_reproducible_output(width = 80, crayon = FALSE)
  # A step that passes
  fit <- fit_model(mod2)
  dg <- diagnostics(fit)
  expect_identical(
    unname(dg[c("vb_lhat", "vb_j0", "vb_j1")]),
    unname(fit$diagnostics$vb$check)
  )
  expect_identical(dg[["vb_solve"]], 0)
  expect_true(is.na(dg[["vb_gain_ratio"]]))
  expect_true(is.na(dg[["vb_step_max"]]))
  # A step that is solved: the step still to go is a number, and not a flag
  mod <- log_gamma_model(0.4)
  short <- inla_vanilla(
    mod$logpost,
    mod$init,
    grad = mod$grad,
    control = list(vb_maxit = 1),
    silent = TRUE
  )
  dg <- diagnostics(short)
  expect_identical(dg[["vb_solve"]], 1)
  expect_lt(dg[["vb_j0"]], -1)
  expect_gt(dg[["vb_step_max"]], 0.1)
  expect_identical(
    diagnostics(short, type = "param")$vb_step_sigma,
    dg[["vb_step_max"]]
  )
  flags <- check_flags(
    diagnostics_internal(
      short$diagnostics$optim,
      short$diagnostics$vb,
      short$diagnostics$marginals
    ),
    gaussian = FALSE
  )
  expect_true(flags[["vb_check"]])
  expect_false(flags[["vb_step"]])
  # No bullet for the step still to go, and the fit counts as converged
  expect_false(any(grepl("did not settle", short$warnings)))
  expect_match(capture_output(print(short)), "Converged in")
  text <- gsub("\\s+", " ", capture_output(print(summary(short))))
  expect_match(
    text,
    "! VB one-step check flagged, most of all in `eta`",
    fixed = TRUE
  )
  expect_match(text, "solved in 1 iteration:", fixed = TRUE)
  expect_no_match(text, "VB correction in one step", fixed = TRUE)
})

# The gain ratio, without a gradient -----------------------------------------

test_that("the gain ratio is the rise delivered over the rise promised", {
  # Whitened: slope sqrt(a) (1 - exp(sigma z)), so the expected log-posterior
  # is known at the nodes, and the ratio can be written out
  a <- 2
  fit <- fit_no_grad(log_gamma_model(a))
  s <- 1 / sqrt(a)
  d <- (fit$mean_vb[[1]] - fit$mode[[1]]) / sqrt(fit$cov[1, 1])
  l <- function(z) sqrt(a) * z - a * (exp(s * z) - 1)
  expected <- function(d) {
    2 / 3 * l(d) + (l(d + sqrt(3)) + l(d - sqrt(3))) / 6
  }
  ratio <- (expected(d) - expected(0)) / (0.5 * d^2)
  expect_equal(fit$diagnostics$vb$gain_ratio, ratio, tolerance = 1e-4)
  expect_equal(
    unname(diagnostics(fit)[["vb_gain_ratio"]]),
    ratio,
    tolerance = 1e-4
  )
})

test_that("the gain ratio is near 1 with data and falls with skewness", {
  gain <- function(a) fit_no_grad(log_gamma_model(a))$diagnostics$vb$gain_ratio
  # Observed 0.997, 0.972, 0.836 and 0.425
  ratios <- vapply(c(50, 5, 1, 0.3), gain, 0)
  expect_true(all(diff(ratios) < 0))
  expect_gt(ratios[1], 0.95)
  expect_lt(ratios[4], 0.5)
})

test_that("a low gain ratio warns, and points to the full solve", {
  testthat::local_reproducible_output(width = 80, crayon = FALSE)
  mod <- log_gamma_model(0.3)
  w <- expect_warning(
    fit <- inla_vanilla(mod$logpost, mod$init, verbose = FALSE),
    class = "inlavanilla_diagnostics_warning"
  ) |>
    suppressMessages()
  text <- gsub("\\s+", " ", cli::ansi_strip(conditionMessage(w)))
  expect_match(text, "One-step VB delivered", fixed = TRUE)
  expect_match(text, "Compare with `vb_correction = \"full\"`", fixed = TRUE)
  expect_match(
    gsub("\\s+", " ", capture_output(print(summary(fit)))),
    "VB correction in one step delivered",
    fixed = TRUE
  )
})

test_that("a healthy gain ratio raises nothing and is printed", {
  testthat::local_reproducible_output(width = 80, crayon = FALSE)
  fit <- fit_no_grad(log_gamma_model(5))
  expect_identical(fit$warnings, character())
  expect_match(
    capture_output(print(summary(fit))),
    "gain ratio 0.972",
    fixed = TRUE
  )
})

test_that("the gain ratio runs only for one step without a gradient", {
  # Not when the correction iterated
  on <- fit_no_grad(skew_normal_model(10))
  expect_false(on$diagnostics$vb$one_pass)
  expect_true(is.na(on$diagnostics$vb$gain_ratio))
  full <- fit_no_grad(log_gamma_model(2), vb_correction = "full")
  expect_true(is.na(full$diagnostics$vb$gain_ratio))
  # With a gradient the check of the step takes its place
  expect_true(is.na(fit_model(log_gamma_model(2))$diagnostics$vb$gain_ratio))
})

test_that("a step too small to check has no gain ratio", {
  # A Gaussian posterior: the first step is zero
  lp <- function(theta) -0.5 * sum(theta^2)
  fit <- inla_vanilla(lp, c(a = 0.5, b = -0.2), silent = TRUE)
  expect_true(is.na(fit$diagnostics$vb$gain_ratio))
  expect_identical(fit$warnings, character())
})

test_that("below_level() flags low and non-finite values, not NA", {
  expect_identical(
    below_level(c(0.2, 0.5, 0.9, NA, NaN, Inf, -Inf), 0.5),
    c(TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE)
  )
})

# The guard, without a gradient ----------------------------------------------

test_that("a first step over the guard goes on to iterate", {
  # The first step of this strongly skewed posterior is 3.3 SDs
  mod <- skew_normal_model(10)
  on <- fit_no_grad(mod)
  full <- fit_no_grad(mod, vb_correction = "full")
  expect_false(on$diagnostics$vb$one_pass)
  expect_gt(on$diagnostics$vb$iterations, 1L)
  sds <- sqrt(diag(on$cov))
  expect_lt(abs(on$mean_vb[[1]] - full$mean_vb[[1]]) / sds[[1]], 0.02)
  expect_lt(max(on$diagnostics$marginals$vb_step_sigma), 0.1)
})

test_that("a first step under the guard is taken whole", {
  # Log-gamma with shape 0.2: strongly skewed, with a first step of 2.3 SDs
  on <- fit_no_grad(log_gamma_model(0.2))
  expect_true(on$diagnostics$vb$one_pass)
  expect_identical(on$diagnostics$vb$iterations, 1L)
  expect_true(all(is.na(on$diagnostics$marginals$vb_step_sigma)))
})

test_that("the iteration after the guard stops within the budget", {
  mod <- skew_normal_model(10)
  ctrl <- vanilla_control()
  lap <- laplace_fit(mod$logpost, NULL, mod$init, ctrl)
  calls <- 0L
  lp <- function(theta) {
    calls <<- calls + 1L
    mod$logpost(theta)
  }
  # Two calls is a pass of the three-point rule, so six is three passes: two
  # steps, and one pass kept for the step still to go
  vb <- vb_correct(
    lp,
    NULL,
    lap$mode,
    lap$Sigma,
    lap$R,
    lap$A,
    ctrl,
    method = "onestep",
    budget = 6
  )
  expect_identical(vb$iterations, 2L)
  expect_true(vb$history$ran_out)
  expect_false(vb$fallback)
  expect_identical(calls, 6L)
  expect_gt(max(vb$step_sigma), 0.1)
})
