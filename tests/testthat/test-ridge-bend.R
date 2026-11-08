# The ridge bend without a gradient. The calls that gamma1_lp() already makes
# give the direction in which the ridge bends, so a parameter whose ridge
# barely bends costs two more calls, and only a bend that can reach the flag
# takes the full gradient. For a linear regression with an unknown SD and flat
# priors the ridge bend of every coefficient is exactly 4 / n, all of it along
# log sigma, which has no bend of its own (lab/bend-threshold.md).

## ----- A linear regression with an exact ridge bend --------------------------

reg_model <- function(n, p = 5L) {
  set.seed(1)
  x <- cbind(1, matrix(rnorm(n * (p - 1L)), n, p - 1L))
  y <- drop(x %*% rnorm(p)) + rnorm(n, 0, 1.5)
  list(
    logpost = function(theta) {
      eta <- drop(x %*% theta[seq_len(p)])
      sum(dnorm(y, eta, exp(theta[[p + 1L]]), log = TRUE))
    },
    init = stats::setNames(
      numeric(p + 1L),
      c(paste0("b", seq_len(p)), "ls")
    )
  )
}

# The calls of logpost that fit_marginals() makes without a gradient, and what
# it returns
marginals_no_grad <- function(mod, ...) {
  ctrl <- vanilla_control(list())
  lap <- laplace_fit(mod$logpost, NULL, mod$init, ctrl)
  calls <- 0L
  lp <- function(x) {
    calls <<- calls + 1L
    mod$logpost(x)
  }
  out <- fit_marginals(lp, NULL, lap, ctrl, ...)
  list(
    calls = calls,
    bend = out$sn[, "ridge_bend"],
    lap = lap,
    out = out,
    lp = lp
  )
}

# With n_grid = 21: 2m calls at the mode, then for each parameter 20 scan
# points, 2m + 5 calls of gamma1_lp(), and 2 of the screen
calls_screened <- function(m) 2L * m + m * (20L + 2L * m + 5L + 2L)

## ----- The screen ------------------------------------------------------------

test_that("a ridge bend under the screen costs two calls and is exact here", {
  res <- marginals_no_grad(reg_model(100))
  m <- 6L
  expect_identical(res$calls, calls_screened(m))
  expect_equal(unname(res$bend[1:5]), rep(4 / 100, 5), tolerance = 1e-3)
  expect_lt(res$bend[["ls"]], 1e-6)
})

test_that("a ridge bend above the screen takes the full gradient", {
  res <- marginals_no_grad(reg_model(20))
  m <- 6L
  # The five coefficients pay for 2(m - 1) more calls, log sigma does not
  expect_identical(res$calls, calls_screened(m) + 5L * 2L * (m - 1L))
  expect_equal(unname(res$bend[1:5]), rep(4 / 20, 5), tolerance = 1e-3)
  expect_lt(res$bend[["ls"]], 1e-6)

  # The values reported are those of the route without the screen
  lap <- res$lap
  v_scan <- sweep(lap$Sigma, 2, sqrt(diag(lap$Sigma)), "/")
  z <- res$out$z
  i <- c(which.min(abs(z + 2)), which.min(abs(z - 2)))
  full <- scan_bend(1L, res$lp, NULL, lap, v_scan, z[i], res$out$logpost[i, 1L])
  expect_identical(res$bend[["b1"]], full)
})

test_that("the gain along the direction is a lower bound on the full gain", {
  set.seed(3)
  n <- 60
  x <- cbind(1, rnorm(n), rnorm(n))
  y <- rbinom(n, 1, plogis(drop(x %*% c(-0.5, 1.2, -0.8))))
  mod <- list(
    logpost = function(theta) {
      eta <- drop(x %*% theta)
      sum(y * eta - log1p(exp(eta))) + sum(dnorm(theta, 0, 2.5, log = TRUE))
    },
    init = c(b0 = 0, b1 = 0, b2 = 0)
  )
  res <- marginals_no_grad(mod)
  lap <- res$lap
  v_scan <- sweep(lap$Sigma, 2, sqrt(diag(lap$Sigma)), "/")
  z <- res$out$z
  i <- c(which.min(abs(z + 2)), which.min(abs(z - 2)))
  at_mode <- whitened_lp(res$lp, lap$mode, lap$lp_max, lap$L)
  for (j in 1:3) {
    d <- gamma1_lp(j, res$lp, lap$mode, lap$lp_max, v_scan, lap$L, at_mode)
    f_t <- res$out$logpost[i, j]
    along <- bend_along(j, res$lp, lap, v_scan, z[i], f_t, d$bend_dir)
    full <- scan_bend(j, res$lp, NULL, lap, v_scan, z[i], f_t)
    # Equal to the rounding of the one-sided differences, and most of the gain
    expect_lte(along, full * (1 + 1e-3))
    expect_gt(along, 0.8 * full)
  }
})

test_that("gamma1_lp() gives a unit direction across the scan", {
  mod <- reg_model(100)
  lap <- laplace_fit(mod$logpost, NULL, mod$init, vanilla_control(list()))
  v_scan <- sweep(lap$Sigma, 2, sqrt(diag(lap$Sigma)), "/")
  at_mode <- whitened_lp(mod$logpost, lap$mode, lap$lp_max, lap$L)
  res <- gamma1_lp(
    1L,
    mod$logpost,
    lap$mode,
    lap$lp_max,
    v_scan,
    lap$L,
    at_mode
  )
  a <- lap$L[1L, ] / sqrt(sum(lap$L[1L, ]^2))
  expect_equal(sum(res$bend_dir^2), 1)
  expect_lt(abs(sum(a * res$bend_dir)), 1e-10)
  # On the parameter scale it points along log sigma
  d <- drop(lap$L %*% res$bend_dir)
  expect_gt(abs(d[[6L]]) / sqrt(sum(d^2)), 0.999)
})

test_that("with no finite direction the ridge bend takes the full route", {
  mod <- reg_model(100)
  lap <- laplace_fit(mod$logpost, NULL, mod$init, vanilla_control(list()))
  sds <- sqrt(diag(lap$Sigma))
  v_scan <- sweep(lap$Sigma, 2, sds, "/")
  at_mode <- whitened_lp(mod$logpost, lap$mode, lap$lp_max, lap$L)
  # A wall just past the point where gamma1_lp() starts: the point itself is
  # finite, and one side of the first whitened axis is not
  edge <- lap$mode[[1L]] + 0.0105 * sds[[1L]]
  walled <- function(x) if (x[[1L]] > edge) -Inf else mod$logpost(x)
  res <- gamma1_lp(1L, walled, lap$mode, lap$lp_max, v_scan, lap$L, at_mode)
  expect_null(res$bend_dir)
})

## ----- Other routes through the check ----------------------------------------

test_that("tk_fallback rescans a flagged parameter without a gradient", {
  set.seed(4)
  y <- rnorm(12, 10, 2)
  logpost <- function(theta) {
    sum(dnorm(y, theta[["mu"]], exp(theta[["ls"]]), log = TRUE))
  }
  fit <- inla_vanilla(
    logpost,
    c(mu = mean(y), ls = log(sd(y))),
    tk_fallback = TRUE,
    silent = TRUE
  )
  # 4 / n, above the level of 0.25
  expect_equal(
    fit$diagnostics$marginals["mu", "ridge_bend"],
    4 / 12,
    tolerance = 1e-3
  )
  expect_identical(fit$fallback, "mu")
})

test_that("the other widths keep the full route without a gradient", {
  mod <- reg_model(100, p = 2L)
  trace_fit <- marginals_no_grad(mod)
  cmp_fit <- marginals_no_grad(mod, width = "cmp")
  expect_equal(cmp_fit$bend, trace_fit$bend, tolerance = 1e-3)
  expect_equal(cmp_fit$bend[["b1"]], 4 / 100, tolerance = 1e-3)
})
