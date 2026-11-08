# The one-step VB correction, the default of inla_vanilla(): the `vb`
# argument, the axis rule, the cost of the step, its accuracy against the true
# posterior mean, and the guard that hands a large first step to the iteration.

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

mod2 <- log_gamma_model(2)

# The argument --------------------------------------------------------------

test_that("vb takes a method, an abbreviation, TRUE or FALSE", {
  expect_identical(fit_model(mod2)$vb, "onestep")
  expect_identical(fit_model(mod2, vb = "onestep")$vb, "onestep")
  expect_identical(fit_model(mod2, vb = "full")$vb, "full")
  expect_identical(fit_model(mod2, vb = "off")$vb, "off")
  expect_identical(fit_model(mod2, vb = TRUE)$vb, "onestep")
  expect_identical(fit_model(mod2, vb = FALSE)$vb, "off")
  expect_identical(fit_model(mod2, vb = "one")$vb, "onestep")
  expect_identical(fit_model(mod2, vb = "FULL")$vb, "full")
  expect_identical(fit_model(mod2, vb = "of")$vb, "off")
  # "o" could be onestep or off
  for (bad in list("o", "fast", NA, c(TRUE, FALSE), 1)) {
    expect_error(
      fit_model(mod2, vb = bad),
      "`vb` must be one of",
      class = "inlavanilla_error"
    )
  }
})

test_that("vb = \"off\" is the Laplace centre, and FALSE is the same fit", {
  off <- fit_model(mod2, vb = "off")
  no <- fit_model(mod2, vb = FALSE)
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
    expect_false(fit$diagnostics$vb$fallback)
  }
})

test_that("with a gradient one step costs 2m + 1 gradient calls", {
  mod <- logistic_model(80, 6, 1)
  on <- fit_model(mod, marginals = "gaussian")
  off <- fit_model(mod, vb = "off", marginals = "gaussian")
  expect_identical(on$counts[["grad"]] - off$counts[["grad"]], 13L)
  expect_identical(on$counts[["logpost"]], off$counts[["logpost"]])
})

test_that("without a gradient one step is one pass of logpost calls", {
  mod <- logistic_model(80, 4, 2)
  on <- inla_vanilla(
    mod$logpost,
    mod$init,
    marginals = "gaussian",
    silent = TRUE
  )
  off <- inla_vanilla(
    mod$logpost,
    mod$init,
    vb = "off",
    marginals = "gaussian",
    silent = TRUE
  )
  n_nodes <- nrow(vb_nodes(on$cov)$nodes) - 1L
  expect_identical(on$counts[["logpost"]] - off$counts[["logpost"]], n_nodes)
  expect_identical(on$diagnostics$vb$iterations, 1L)
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
    capture_output(print(fit_model(mod2, vb = "full"))),
    "Skew-normal + full VB --",
    fixed = TRUE
  )
  expect_match(
    capture_output(print(summary(fit_model(mod2)))),
    "VB correction in one step: max shift",
    fixed = TRUE
  )
  expect_match(
    capture_output(print(summary(fit_model(mod2, vb = "full")))),
    "VB correction settled in",
    fixed = TRUE
  )
})

# Accuracy ------------------------------------------------------------------

test_that("one step is about as close to the true mean as the full solve", {
  for (a in c(1, 2, 5, 20)) {
    mod <- log_gamma_model(a)
    on <- fit_model(mod)
    full <- fit_model(mod, vb = "full")
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
    full <- fit_model(mod, vb = "full")
    abs(on$mean_vb[[1]] - full$mean_vb[[1]]) / sqrt(on$cov[1, 1])
  }
  expect_lt(gap(50) / gap(5), 3 * (5 / 50)^1.5)
})

test_that("in several dimensions one step moves most of the way", {
  mod <- logistic_model(150, 8, 3)
  on <- fit_model(mod)
  full <- fit_model(mod, vb = "full")
  sds <- sqrt(diag(on$cov))
  shift_on <- (on$mean_vb - on$mode) / sds
  shift_full <- (full$mean_vb - full$mode) / sds
  # Observed: the largest shift is 0.33 SDs, and the two differ by 0.03
  expect_lt(max(abs(shift_on - shift_full)), 0.05)
  expect_gt(stats::cor(shift_on, shift_full), 0.99)
})

# The guard ------------------------------------------------------------------

test_that("a first step over the guard goes on to iterate", {
  # The first step of this strongly skewed posterior is 2.9 SDs
  mod <- skew_normal_model(10)
  on <- fit_model(mod)
  full <- fit_model(mod, vb = "full")
  expect_false(on$diagnostics$vb$one_pass)
  expect_gt(on$diagnostics$vb$iterations, 1L)
  sds <- sqrt(diag(on$cov))
  expect_lt(abs(on$mean_vb[[1]] - full$mean_vb[[1]]) / sds[[1]], 0.02)
  expect_lt(max(on$diagnostics$marginals$vb_step_sigma), 0.1)
})

test_that("a first step under the guard is taken whole", {
  # Log-gamma with shape 0.5: strongly skewed, with a first step of 1.14 SDs
  mod <- log_gamma_model(0.5)
  on <- fit_model(mod)
  expect_true(on$diagnostics$vb$one_pass)
  shift <- (on$mean_vb[[1]] - on$mode[[1]]) / sqrt(on$cov[1, 1])
  expect_equal(shift, log_gamma_first_step(0.5), tolerance = 1e-5)
})

test_that("a node on a boundary is handled by the guard", {
  # a + 1 ~ Gamma(4, 1): the lowest node at the mode sits on the barrier
  lp <- function(theta) {
    x <- theta[[1]] + 1
    if (x <= 0) -Inf else 3 * log(x) - x
  }
  gr <- function(theta) {
    x <- theta[[1]] + 1
    if (x <= 0) NaN else 3 / x - 1
  }
  on <- inla_vanilla(lp, c(a = 1.5), grad = gr, silent = TRUE)
  full <- inla_vanilla(lp, c(a = 1.5), grad = gr, vb = "full", silent = TRUE)
  expect_false(on$diagnostics$vb$one_pass)
  sds <- sqrt(diag(on$cov))
  expect_lt(abs(on$mean_vb[[1]] - full$mean_vb[[1]]) / sds[[1]], 0.02)
})

test_that("the iteration after the guard stops within the budget", {
  mod <- skew_normal_model(10)
  ctrl <- vanilla_control(list())
  lap <- laplace_fit(mod$logpost, mod$grad, mod$init, ctrl)
  calls <- 0L
  gr <- function(theta) {
    calls <<- calls + 1L
    mod$grad(theta)
  }
  # Six calls is two passes of the three-point rule
  vb <- vb_correct(
    mod$logpost,
    gr,
    lap$mode,
    lap$Sigma,
    lap$R,
    lap$L,
    ctrl,
    method = "onestep",
    budget = 6
  )
  expect_identical(vb$iterations, 2L)
  expect_true(vb$history$ran_out)
  expect_false(vb$fallback)
  # The two passes, and one more for the step still to go
  expect_identical(calls, 9L)
  expect_gt(max(vb$step_sigma), 0.1)
})
