# The one-step VB correction, the default of inla_vanilla(): the `vb`
# argument, the axis rule, the cost of the step and its accuracy against the
# true posterior mean. Then the check of the step from a second sweep, and the
# solve and the warning when it fails, with a gradient (sweeps of gradients on
# the axis rule) and without one (Stein passes on the simplex rule).

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
# directly, with or without the gradient, and with the calls of the correction
# counted (of the gradient with one, of the log-posterior without)
vb_direct <- function(mod, budget = Inf, fixed = integer(), grad = TRUE) {
  ctrl <- vanilla_control()
  lap <- laplace_fit(mod$logpost, if (grad) mod$grad, mod$init, ctrl)
  calls <- 0L
  lp <- function(theta) {
    if (!grad) {
      calls <<- calls + 1L
    }
    mod$logpost(theta)
  }
  gr <- function(theta) {
    calls <<- calls + 1L
    mod$grad(theta)
  }
  vb <- vb_correct(
    lp,
    if (grad) gr,
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

# The one step as it was before the check existed. With a gradient: Sigma
# times the expected gradient on axis_nodes(). Without one: a pass of Stein's
# identity on vb_nodes(). The step of the fixed parameters is zero.
plain_onestep <- function(mod, lap, fixed = integer(), grad = TRUE) {
  if (grad) {
    rule <- axis_nodes(lap$Sigma, lap$A)
    g0 <- mod$grad(lap$mode)
    s <- g0
    for (b in 2:nrow(rule$nodes)) {
      s <- s + rule$weights[b] * (mod$grad(lap$mode + rule$nodes[b, ]) - g0)
    }
    step <- as.numeric(lap$Sigma %*% s)
  } else {
    rule <- vb_nodes(lap$Sigma, lap$A)
    n_pairs <- (nrow(rule$nodes) - 1L) %/% 2L
    plus <- seq_len(n_pairs)
    v <- vapply(
      2:nrow(rule$nodes),
      function(b) mod$logpost(lap$mode + rule$nodes[b, ]),
      0
    )
    step <- as.numeric(crossprod(
      rule$nodes[1L + plus, , drop = FALSE],
      rule$weights[1L + plus] * (v[plus] - v[n_pairs + plus])
    ))
  }
  step[fixed] <- 0
  stats::setNames(step, names(lap$mode))
}

# The calls of one sweep: 2m + 1 gradients, or the values of a Stein pass
sweep_calls <- function(lap, grad) {
  if (grad) {
    2L * length(lap$mode) + 1L
  } else {
    nrow(vb_nodes(lap$Sigma, lap$A)$nodes) - 1L
  }
}

# The whitened gradient of log_gamma_model(a), and the residual of the VB
# equation at a whitened shift d: on the three-point rule of the gradient, or,
# without a gradient, by Stein's identity on the values of the whitened
# log-posterior at the same nodes
log_gamma_slope <- function(a, z) sqrt(a) * (1 - exp(z / sqrt(a)))
log_gamma_resid <- function(a, d, grad = TRUE) {
  if (grad) {
    centre <- log_gamma_slope(a, d)
    ends <- log_gamma_slope(a, d + sqrt(3)) + log_gamma_slope(a, d - sqrt(3))
    (4 * centre + ends) / 6
  } else {
    l <- function(z) sqrt(a) * z - a * (exp(z / sqrt(a)) - 1)
    sqrt(3) / 6 * (l(d + sqrt(3)) - l(d - sqrt(3)))
  }
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

test_that("without a gradient one step costs two Stein passes and 8m calls", {
  mod <- logistic_model(80, 4, 2)
  on <- fit_no_grad(mod, marginals = "gaussian")
  off <- fit_no_grad(mod, vb_correction = "none", marginals = "gaussian")
  n_nodes <- nrow(vb_nodes(on$cov)$nodes) - 1L
  # The pass of the step, the pass of its check, and four values per axis
  # for the slope along the step at each end
  expect_identical(
    on$counts[["logpost"]] - off$counts[["logpost"]],
    2L * n_nodes + 32L
  )
  expect_identical(on$diagnostics$vb$iterations, 1L)
  expect_false(on$diagnostics$vb$solve)
  expect_true(all(is.finite(on$diagnostics$vb$check)))
})

test_that("one step has no settling check, which is not a flag", {
  for (fit in list(fit_model(mod2), fit_no_grad(mod2))) {
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
  }
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
  for (fit in list(fit_model(mod2), fit_no_grad(mod2))) {
    expect_match(
      capture_output(print(summary(fit))),
      "VB correction in one step: max shift [0-9.]+ SD\\."
    )
  }
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

# The check of the step -------------------------------------------------------

test_that("a step that passes the check is the plain one step, to the bit", {
  for (grad in c(TRUE, FALSE)) {
    for (mod in list(
      log_gamma_model(2),
      logistic_model(80, 6, 1),
      mixed_gamma_model(c(2, 5, 20))
    )) {
      out <- vb_direct(mod, grad = grad)
      plain <- plain_onestep(mod, out$lap, grad = grad)
      expect_false(out$vb$solve)
      expect_true(out$vb$one_pass)
      expect_identical(out$vb$shift, plain)
      expect_identical(
        out$vb$kld_global,
        0.5 * sum((out$lap$R %*% plain)^2)
      )
      expect_true(all(is.na(out$vb$step_sigma)))
      expect_identical(out$vb$history$sweeps, 2L)
      # Two sweeps, and without a gradient four values per axis for each of
      # the two slopes
      expect_identical(
        out$calls,
        2L * sweep_calls(out$lap, grad) + if (grad) 0L else 8L * length(plain)
      )
    }
  }
})

test_that("the check is the contraction ratio and the slope at both ends", {
  # On the three-point rule the whitened residual and its slope are known.
  # The slope is that of the gradient's rule on both routes: without a
  # gradient it comes from differences of the log-posterior, not from Stein.
  for (grad in c(TRUE, FALSE)) {
    # Without a gradient the mode and the Laplace SD are less exact, and the
    # numbers of the check with them. Observed 6e-5.
    tol <- if (grad) 1e-5 else 5e-4
    for (a in c(0.4, 2)) {
      mod <- log_gamma_model(a)
      fit <- if (grad) fit_model(mod) else fit_no_grad(mod)
      d1 <- log_gamma_resid(a, 0, grad)
      slope <- function(d) {
        1 +
          (log_gamma_slope(a, d + sqrt(3)) - log_gamma_slope(a, d - sqrt(3))) /
            (2 * sqrt(3))
      }
      expect_equal(
        fit$diagnostics$vb$check,
        c(
          lhat = abs(log_gamma_resid(a, d1, grad)) / abs(d1),
          j0 = slope(0),
          j1 = slope(d1)
        ),
        tolerance = tol
      )
      # The slope at the mode in closed form
      x <- sqrt(3 / a)
      expect_equal(
        fit$diagnostics$vb$check[["j0"]],
        1 - sinh(x) / x,
        tolerance = tol
      )
    }
  }
})

test_that("the check passes where the map contracts, and fails where not", {
  # The slope at the mode of the log-gamma crosses -1 at a = 0.633
  for (fit_fun in list(fit_model, fit_no_grad)) {
    solved <- function(a) fit_fun(log_gamma_model(a))$diagnostics$vb$solve
    expect_false(solved(2))
    expect_false(solved(0.65))
    expect_true(solved(0.62))
    expect_true(solved(0.4))
  }
})

test_that("a failed check solves the VB equation of the rule", {
  a <- 0.4
  mod <- log_gamma_model(a)
  for (grad in c(TRUE, FALSE)) {
    fit_fun <- if (grad) fit_model else fit_no_grad
    fit <- fit_fun(mod)
    vb <- fit$diagnostics$vb
    expect_true(vb$solve)
    expect_false(vb$one_pass)
    expect_false(vb$fallback)
    expect_gt(vb$iterations, 1L)
    # The root of the rule's equation, which is not the one step
    root <- stats::uniroot(
      function(d) log_gamma_resid(a, d, grad),
      c(-3, 0),
      tol = 1e-12
    )$root
    s <- sqrt(fit$cov[1, 1])
    shift <- (fit$mean_vb[[1]] - fit$mode[[1]]) / s
    expect_equal(shift, root, tolerance = 2e-3)
    expect_gt(abs(log_gamma_resid(a, 0, grad) - root), 0.4)
    # In one dimension the full solve uses the same rule
    full <- fit_fun(mod, vb_correction = "full")
    expect_lt(abs(fit$mean_vb[[1]] - full$mean_vb[[1]]) / s, 2e-3)
    # The step still to go is recorded, below vb_tol where the solve converged
    expect_lt(fit$diagnostics$marginals$vb_step_sigma, 1e-3)
    # The evidence takes the KL term at the landing
    expect_equal(
      fit$logml - fit$logml_laplace,
      0.5 * shift^2,
      tolerance = 1e-10
    )
  }
})

test_that("a flagged step warns once, in the agreed words", {
  testthat::local_reproducible_output(width = 80, crayon = FALSE)
  mod <- log_gamma_model(0.4)
  for (grad in c(TRUE, FALSE)) {
    fit <- if (grad) fit_model(mod) else fit_no_grad(mod)
    expect_identical(
      fit$warnings[startsWith(fit$warnings, "VB one-step check")],
      check_warning("eta")
    )
    # Raised as a warning of its own, apart from the diagnostics warning
    seen <- list()
    withCallingHandlers(
      suppressMessages(
        inla_vanilla(mod$logpost, mod$init, grad = if (grad) mod$grad)
      ),
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
    expect_s3_class(seen[[which(is_check)]], "inlavanilla_vb_check_warning")
    expect_s3_class(seen[[which(is_check)]], "inlavanilla_warning")
    expect_false(inherits(
      seen[[which(is_check)]],
      "inlavanilla_diagnostics_warning"
    ))
  }
  # A step that passes raises nothing
  expect_identical(fit_model(log_gamma_model(2))$warnings, character())
  expect_identical(fit_no_grad(log_gamma_model(2))$warnings, character())
})

test_that("the warning names the parameter the one step moves furthest", {
  mod <- mixed_gamma_model(c(20, 0.3, 5), names = c("tight", "wide", "mid"))
  for (grad in c(TRUE, FALSE)) {
    out <- vb_direct(mod, grad = grad)
    expect_true(out$vb$solve)
    plain <- abs(plain_onestep(mod, out$lap, grad = grad)) /
      sqrt(diag(out$lap$Sigma))
    expect_identical(out$vb$check_worst, names(plain)[which.max(plain)])
    fit <- if (grad) fit_model(mod) else fit_no_grad(mod)
    expect_identical(
      fit$warnings[startsWith(fit$warnings, "VB one-step check")],
      check_warning(out$vb$check_worst)
    )
  }
})

test_that("the sweeps are capped by twice the calls of the mode and Hessian", {
  expect_identical(vb_max_sweeps(104, 3), 29)
  expect_identical(vb_max_sweeps(20, 3), 5)
  expect_identical(vb_max_sweeps(3, 3), 0)
  expect_identical(vb_max_sweeps(Inf, 3), Inf)
  # `vb_budget` is the factor: 1 is the cap of 26.11.25-1, and Inf is no cap
  expect_identical(vb_max_sweeps(104, 3, 1), 14)
  expect_identical(vb_max_sweeps(104, 3, 0), 0)
  expect_identical(vb_max_sweeps(104, 3, Inf), Inf)
  expect_identical(vanilla_control()$vb_budget, 2)
  expect_identical(vanilla_control(vb_budget = Inf)$vb_budget, Inf)
  for (grad in c(TRUE, FALSE)) {
    for (mod in list(
      log_gamma_model(0.4),
      mixed_gamma_model(c(20, 0.3, 5))
    )) {
      m <- length(mod$init)
      togo <- function(out) max(out$vb$step_sigma)
      last <- Inf
      # The budget is in calls, and the cap is twice that budget over 2m + 1,
      # on both routes
      for (cap in c(0, 2, 2.9, 3, 4, 6, 10)) {
        out <- vb_direct(mod, budget = cap * (2 * m + 1) / 2, grad = grad)
        expect_true(out$vb$solve)
        # The two sweeps of the check are always taken
        expect_lte(out$vb$history$sweeps, max(floor(cap), 2))
        expect_identical(out$vb$history$max_sweeps, floor(cap))
        expect_identical(
          out$calls,
          out$vb$history$sweeps *
            sweep_calls(out$lap, grad) +
            if (grad) 0L else 8L * m
        )
        expect_true(all(is.finite(out$vb$shift)))
        expect_lte(togo(out), last + 1e-12)
        last <- togo(out)
      }
      # With no sweep left for the solve the mode stands, and the step still
      # to go is the one step
      none <- vb_direct(mod, budget = 2.9 * (2 * m + 1) / 2, grad = grad)
      plain <- plain_onestep(mod, none$lap, grad = grad)
      expect_identical(none$vb$history$sweeps, 2L)
      expect_identical(unname(none$vb$shift), numeric(m))
      expect_identical(none$vb$iterations, 0L)
      expect_identical(none$vb$kld_global, 0)
      expect_equal(
        none$vb$step_sigma,
        abs(plain) / sqrt(diag(none$lap$Sigma))
      )
    }
  }
})

test_that("vb_budget sets the cap of a whole fit", {
  mod <- mixed_gamma_model(c(20, 0.3, 5))
  m <- length(mod$init)
  fit_at <- function(b) {
    fit_model(mod, marginals = "gaussian", control = list(vb_budget = b))
  }
  off <- fit_model(mod, vb_correction = "none", marginals = "gaussian")
  budget <- off$counts[["logpost"]] + off$counts[["grad"]]
  sweeps <- function(fit) {
    (fit$counts[["grad"]] - off$counts[["grad"]]) / (2 * m + 1)
  }
  # No cap: the solve runs until it converges
  free <- fit_at(Inf)
  expect_true(free$diagnostics$vb$solve)
  expect_gt(sweeps(free), 2)
  expect_lt(diagnostics(free)[["vb_step_max"]], 1e-3)
  # The default is 2
  expect_identical(
    fit_model(mod, marginals = "gaussian")$mean_vb,
    fit_at(2)$mean_vb
  )
  for (b in c(0, 0.5, 1, 2)) {
    fit <- fit_at(b)
    cap <- floor(b * budget / (2 * m + 1))
    expect_true(fit$diagnostics$vb$solve)
    expect_lte(sweeps(fit), max(cap, 2))
    # A cap that the free solve fits in changes nothing
    if (cap >= sweeps(free)) {
      expect_identical(fit$mean_vb, free$mean_vb)
    }
  }
  # With no sweep for the solve a flagged fit stays at its mode
  expect_identical(unname(fit_at(0)$mean_vb), unname(off$mode))
})

test_that("a whole fit keeps the VB stage within twice the mode and Hessian", {
  mod <- mixed_gamma_model(c(20, 0.3, 5))
  on <- fit_model(mod, marginals = "gaussian")
  off <- fit_model(mod, vb_correction = "none", marginals = "gaussian")
  expect_true(on$diagnostics$vb$solve)
  expect_identical(on$counts[["logpost"]], off$counts[["logpost"]])
  expect_lte(
    on$counts[["grad"]] - off$counts[["grad"]],
    2 * (off$counts[["logpost"]] + off$counts[["grad"]])
  )
})

test_that("a fit with no gradient at all is checked, solved and capped", {
  # The mode and the Hessian are found without a gradient too
  mod <- mixed_gamma_model(c(20, 0.3, 5))
  m <- length(mod$init)
  on <- fit_no_grad(mod, marginals = "gaussian")
  off <- fit_no_grad(mod, vb_correction = "none", marginals = "gaussian")
  with <- fit_model(mod, marginals = "gaussian")
  expect_identical(on$counts[["grad"]], 0L)
  expect_true(on$diagnostics$vb$solve)
  # The same verdict as with a gradient, from slopes that agree. Observed
  # -1.302 against -1.290: the two fits differ a little in mode and Hessian.
  expect_equal(
    on$diagnostics$vb$check[["j0"]],
    with$diagnostics$vb$check[["j0"]],
    tolerance = 0.03
  )
  # The sweeps are Stein passes, capped by twice the calls of this fit's own
  # mode and Hessian stage over 2m + 1. The slopes are not sweeps.
  n_nodes <- nrow(vb_nodes(on$cov)$nodes) - 1L
  sweeps <- (on$counts[["logpost"]] - off$counts[["logpost"]] - 8L * m) /
    n_nodes
  expect_identical(sweeps, round(sweeps))
  expect_gt(sweeps, 2)
  expect_lte(sweeps, floor(2 * off$counts[["logpost"]] / (2 * m + 1)))
  # Its mean is close to the mean of the solve with a gradient
  sds <- sqrt(diag(with$cov))
  expect_lt(max(abs(on$mean_vb - with$mean_vb) / sds), 0.15)
  # A model that passes takes its Stein step, as before the check
  mod <- mixed_gamma_model(c(2, 5, 20))
  on <- fit_no_grad(mod)
  expect_false(on$diagnostics$vb$solve)
  expect_identical(on$warnings, character())
})

test_that("the checked step does not depend on the order of the parameters", {
  ctrl <- vanilla_control()
  perm <- c(3L, 1L, 2L)
  back <- order(perm)
  for (grad in c(TRUE, FALSE)) {
    # One model that passes the check and one that is solved
    for (a in list(c(2, 5, 20), c(20, 0.3, 5))) {
      mod <- mixed_gamma_model(a)
      out <- vb_direct(mod, budget = 70, grad = grad)
      lap <- out$lap
      Sigma <- lap$Sigma[perm, perm]
      moved <- vb_correct(
        function(theta) mod$logpost(theta[back]),
        if (grad) function(theta) mod$grad(theta[back])[perm],
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
      expect_equal(moved$check, out$vb$check, tolerance = 1e-7)
      expect_equal(moved$shift, out$vb$shift[perm], tolerance = 1e-7)
      expect_equal(
        moved$step_sigma,
        out$vb$step_sigma[perm],
        tolerance = 1e-5
      )
    }
  }
})

test_that("the shift of a fixed parameter stays at zero", {
  for (grad in c(TRUE, FALSE)) {
    # A step that passes: the plain one step with that entry at zero
    mod <- mixed_gamma_model(c(2, 5, 20))
    out <- vb_direct(mod, fixed = 2L, grad = grad)
    expect_false(out$vb$solve)
    expect_identical(
      out$vb$shift,
      plain_onestep(mod, out$lap, fixed = 2L, grad = grad)
    )
    expect_identical(out$vb$shift[[2]], 0)
    # A step that is solved
    mod <- mixed_gamma_model(c(20, 0.15, 5))
    out <- vb_direct(mod, fixed = 3L, grad = grad)
    expect_true(out$vb$solve)
    expect_identical(out$vb$shift[[3]], 0)
    expect_identical(out$vb$step_sigma[[3]], 0)
    expect_true(all(out$vb$shift[-3] != 0))
  }
})

test_that("a node of the second sweep that is not finite is handled", {
  # Shape 0.7 passes the check when the support is the whole line. Here the
  # lowest node of the sweep at the mode is inside the support, and the one
  # step moves the mean down, so the lowest node of the second sweep is not:
  # the gradient is not finite there, and the log-posterior is -Inf.
  testthat::local_reproducible_output(width = 80, crayon = FALSE)
  cut <- -3.3
  mod <- cut_log_gamma_model(0.7, cut)
  for (grad in c(TRUE, FALSE)) {
    fit_fun <- if (grad) fit_model else fit_no_grad
    expect_false(fit_fun(log_gamma_model(0.7))$diagnostics$vb$solve)
    out <- vb_direct(mod, budget = 60, grad = grad)
    vb <- out$vb
    low <- function(shift) {
      out$lap$mode[[1]] + shift - sqrt(3 * out$lap$Sigma[1, 1])
    }
    plain <- plain_onestep(mod, out$lap, grad = grad)
    expect_gt(low(0), cut)
    expect_true(is.finite(plain))
    expect_lt(low(plain[[1]]), cut)
    # The residual at the step is not finite, so the check fails: a value
    # that is not finite, and not the NA of a check that did not run
    expect_false(is.finite(vb$check[["lhat"]]))
    expect_false(is.na(vb$check[["lhat"]]) && !is.nan(vb$check[["lhat"]]))
    expect_true(is.nan(vb$check[["j1"]]))
    expect_true(is.finite(vb$check[["j0"]]))
    expect_true(vb$solve)
    # The solve lands where every node is finite, at the root of the rule
    expect_true(is.finite(vb$shift))
    expect_gt(low(vb$shift[[1]]), cut)
    expect_lt(vb$step_sigma[[1]], 1e-3)
    root <- stats::uniroot(
      function(d) log_gamma_resid(0.7, d, grad),
      c(-3, 0),
      tol = 1e-12
    )$root
    expect_equal(
      vb$shift[[1]] / sqrt(out$lap$Sigma[1, 1]),
      root,
      tolerance = 2e-3
    )
    expect_lte(vb$history$sweeps, 20L)
    # And a whole fit warns as for any flagged step
    fit <- fit_fun(mod)
    expect_true(fit$diagnostics$vb$solve)
    expect_true(check_warning("eta") %in% fit$warnings)
  }
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
  for (fit_fun in list(fit_model, fit_no_grad)) {
    on <- fit_fun(mod)
    full <- fit_fun(mod, vb_correction = "full")
    vb <- on$diagnostics$vb
    expect_true(vb$fallback)
    expect_false(vb$one_pass)
    expect_true(is.na(vb$solve))
    expect_true(all(is.na(vb$check)))
    sds <- sqrt(diag(on$cov))
    expect_lt(abs(on$mean_vb[[1]] - full$mean_vb[[1]]) / sds[[1]], 0.02)
    expect_false(any(startsWith(on$warnings, "VB one-step check")))
  }
})

test_that("a step too small to check is taken with one sweep", {
  # A Gaussian posterior: the first step is zero
  lp <- function(theta) -0.5 * sum(theta^2)
  gr <- function(theta) -theta
  mod <- list(logpost = lp, grad = gr, init = c(a = 0.5, b = -0.2))
  for (grad in c(TRUE, FALSE)) {
    out <- vb_direct(mod, grad = grad)
    expect_identical(out$calls, sweep_calls(out$lap, grad))
    expect_identical(out$vb$history$sweeps, 1L)
    expect_true(out$vb$one_pass)
    expect_true(is.na(out$vb$solve))
    expect_true(all(is.na(out$vb$check)))
    fit <- if (grad) fit_model(mod) else fit_no_grad(mod)
    expect_identical(fit$warnings, character())
    expect_true(all(is.na(
      diagnostics(fit)[c("vb_lhat", "vb_j0", "vb_j1", "vb_solve")]
    )))
  }
})

test_that("diagnostics(), print() and summary() report the check", {
  testthat::local_reproducible_output(width = 80, crayon = FALSE)
  # A step that passes, with the same entries with and without a gradient
  fit <- fit_model(mod2)
  dg <- diagnostics(fit)
  expect_identical(
    unname(dg[c("vb_lhat", "vb_j0", "vb_j1")]),
    unname(fit$diagnostics$vb$check)
  )
  expect_identical(dg[["vb_solve"]], 0)
  expect_true(is.na(dg[["vb_step_max"]]))
  dg_no_grad <- diagnostics(fit_no_grad(mod2))
  expect_identical(names(dg_no_grad), names(dg))
  expect_identical(dg_no_grad[["vb_solve"]], 0)
  expect_false("vb_gain_ratio" %in% names(dg))
  # A step that is solved: the step still to go is a number, and not a flag
  mod <- log_gamma_model(0.4)
  for (grad in c(TRUE, FALSE)) {
    short <- inla_vanilla(
      mod$logpost,
      mod$init,
      grad = if (grad) mod$grad,
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
  }
})
