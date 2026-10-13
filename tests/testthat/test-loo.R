# loo() and waic(): the deterministic second-order Taylor approximation.
#
# (a) Exact cases. For a linear regression with known noise SD and a Gaussian
#     prior, each unit's log-likelihood is quadratic and the posterior is
#     Gaussian, so the second-order expansion is exact: the log CPO, the lpd
#     and the WAIC penalty all have closed forms.
# (b) A skewed case, a Poisson regression, against the exact LOO by 2D grid
#     quadrature.
# (c) The fallback to first order, the loo package's layout, the interface
#     and the printed output.

set.seed(20)
lin_n <- 25
lin_X <- cbind(1, rnorm(lin_n), runif(lin_n))
lin_sigma <- 0.7
lin_tau <- 2
lin_y <- drop(lin_X %*% c(1, -0.5, 2) + rnorm(lin_n, 0, lin_sigma))

ll_lin <- function(theta, X, y) {
  dnorm(y, drop(X %*% theta), lin_sigma, log = TRUE)
}
lpost_lin <- function(theta, X, y) {
  sum(ll_lin(theta, X, y)) + sum(dnorm(theta, 0, lin_tau, log = TRUE))
}
fit_lin <- inla_vanilla(
  lpost_lin,
  c(a = 0, b = 0, c = 0),
  X = lin_X,
  y = lin_y,
  silent = TRUE
)

# Conjugate posterior of the coefficients given rows `rows`
lin_post <- function(rows) {
  X <- lin_X[rows, , drop = FALSE]
  Q <- crossprod(X) / lin_sigma^2 + diag(3) / lin_tau^2
  V <- solve(Q)
  list(m = drop(V %*% crossprod(X, lin_y[rows])) / lin_sigma^2, V = V)
}
# log p(y_i | posterior given rows): the predictive is normal
lin_pred <- function(i, post) {
  x <- lin_X[i, ]
  dnorm(
    lin_y[i],
    sum(x * post$m),
    sqrt(lin_sigma^2 + drop(t(x) %*% post$V %*% x)),
    log = TRUE
  )
}

test_that("LOO is exact when the expansion is exact", {
  res <- loo(fit_lin, ll_lin, X = lin_X, y = lin_y)
  full <- lin_post(seq_len(lin_n))
  cpo <- vapply(
    seq_len(lin_n),
    function(i) lin_pred(i, lin_post(-i)),
    numeric(1)
  )
  lpd <- vapply(seq_len(lin_n), function(i) lin_pred(i, full), numeric(1))
  # observed 1.5e-10
  expect_equal(unname(res$pointwise[, "elpd_loo"]), cpo, tolerance = 1e-8)
  expect_equal(
    res$estimates["elpd_loo", "Estimate"],
    sum(cpo),
    tolerance = 1e-8
  )
  expect_equal(
    res$estimates["p_loo", "Estimate"],
    sum(lpd - cpo),
    tolerance = 1e-7
  )
  expect_equal(
    res$estimates["looic", "Estimate"],
    -2 * sum(cpo),
    tolerance = 1e-8
  )
  expect_equal(
    res$estimates["elpd_loo", "SE"],
    sqrt(lin_n * var(cpo)),
    tolerance = 1e-7
  )
  expect_identical(res$order, 2L)
  # m^2 + m + 1 calls for m = 3
  expect_identical(res$calls, 13L)
  expect_identical(res$theta, fit_lin$mean_vb)
})

test_that("WAIC is exact when the expansion is exact", {
  res <- waic(fit_lin, ll_lin, X = lin_X, y = lin_y)
  full <- lin_post(seq_len(lin_n))
  lpd <- vapply(seq_len(lin_n), function(i) lin_pred(i, full), numeric(1))
  # Var of l_i = c - r^2 / (2 sigma^2), with r ~ N(mu, v) under the posterior
  p_waic <- vapply(
    seq_len(lin_n),
    function(i) {
      x <- lin_X[i, ]
      mu <- lin_y[i] - sum(x * full$m)
      v <- drop(t(x) %*% full$V %*% x)
      (4 * mu^2 * v + 2 * v^2) / (4 * lin_sigma^4)
    },
    numeric(1)
  )
  expect_equal(unname(res$pointwise[, "p_waic"]), p_waic, tolerance = 1e-7)
  expect_equal(
    res$estimates["elpd_waic", "Estimate"],
    sum(lpd - p_waic),
    tolerance = 1e-8
  )
  expect_equal(
    res$estimates["waic", "Estimate"],
    -2 * sum(lpd - p_waic),
    tolerance = 1e-8
  )
  expect_identical(colnames(res$pointwise), c("elpd_waic", "p_waic", "waic"))
  expect_s3_class(res, c("inla_vanilla_waic", "waic", "loo"), exact = TRUE)
})

# (b) Poisson regression against quadrature ----------------------------------

pois_x <- c(-1.75, -1.25, -0.75, -0.25, 0.25, 0.75, 1.25, 1.75)
pois_y <- c(1, 0, 2, 1, 3, 4, 2, 6)
ll_pois <- function(theta, x, y) {
  dpois(y, exp(theta[["a"]] + theta[["b"]] * x), log = TRUE)
}
lpost_pois <- function(theta, x, y) {
  sum(ll_pois(theta, x, y)) + sum(dnorm(theta, 0, 2, log = TRUE))
}

test_that("LOO is close to the exact LOO for a skewed posterior", {
  # A data argument named `x` reaches loglik: the generic's first argument
  # is `object`
  fit <- inla_vanilla(
    lpost_pois,
    c(a = 0, b = 0),
    x = pois_x,
    y = pois_y,
    silent = TRUE
  )
  res <- loo(fit, ll_pois, x = pois_x, y = pois_y)
  g <- expand.grid(
    a = seq(-2.5, 2.5, length.out = 301),
    b = seq(-1.5, 2.5, length.out = 301)
  )
  L <- vapply(
    seq_along(pois_y),
    function(i) dpois(pois_y[i], exp(g$a + g$b * pois_x[i]), log = TRUE),
    numeric(nrow(g))
  )
  total <- rowSums(L) +
    dnorm(g$a, 0, 2, log = TRUE) +
    dnorm(g$b, 0, 2, log = TRUE)
  cpo <- vapply(
    seq_along(pois_y),
    function(i) {
      w <- total - L[, i]
      mx <- max(w)
      log(sum(exp(w - mx + L[, i]))) - log(sum(exp(w - mx)))
    },
    numeric(1)
  )
  # Eight observations: observed errors up to 0.016 per unit and 0.054 in
  # total, from the Gaussian posterior and the truncated expansion
  expect_lt(max(abs(res$pointwise[, "elpd_loo"] - cpo)), 0.03)
  expect_lt(abs(res$estimates["elpd_loo", "Estimate"] - sum(cpo)), 0.1)
  expect_true(all(res$pointwise[, "k"] > 0 & res$pointwise[, "k"] < 1))
})

test_that("the expansion is at the VB mean, or at the mode without VB", {
  fit <- inla_vanilla(
    lpost_pois,
    c(a = 0, b = 0),
    x = pois_x,
    y = pois_y,
    vb = FALSE,
    silent = TRUE
  )
  res <- loo(fit, ll_pois, x = pois_x, y = pois_y)
  expect_identical(res$theta, fit$mode)
})

# (c) Fallback, layout, interface, output -----------------------------------

test_that("a unit with no second-order term sends every estimate to first order", {
  tp <- list(
    n = 3L,
    units = c("u1", "u2", "u3"),
    l0 = c(-1, -2, -3),
    quad1 = c(0.1, 0.2, 0.3),
    ssq = c(0.01, 0.02, 0.03),
    log_cpo_1 = c(-1.05, -2.1, -3.15),
    lpd_1 = c(-0.95, -1.9, -2.85),
    log_cpo_2 = c(-1.06, NA, -3.16),
    lpd_2 = c(-0.94, -1.89, NA),
    k_max = c(0.2, 1.3, 0.1),
    k_min = c(0, 0, -1.2),
    theta = c(a = 0),
    calls = 3L
  )
  expect_warning(
    res <- loo_from_taylor(tp),
    "first-order LOO",
    class = "inlavanilla_warning"
  )
  expect_identical(res$order, 1L)
  expect_identical(res$first_order_units, "u2")
  expect_equal(unname(res$pointwise[, "elpd_loo"]), tp$log_cpo_1)
  expect_equal(unname(res$pointwise[, "p_loo"]), tp$quad1)
  expect_match(res$warnings, "first-order LOO")

  expect_warning(
    res <- waic_from_taylor(tp),
    "first-order WAIC",
    class = "inlavanilla_warning"
  )
  expect_identical(res$order, 1L)
  expect_identical(res$first_order_units, "u3")
  # At first order the WAIC equals the first-order LOO
  expect_equal(unname(res$pointwise[, "elpd_waic"]), tp$log_cpo_1)

  # A missing lpd alone keeps the LOO at second order; p_loo takes the
  # unit's first-order term
  tp$log_cpo_2[2] <- -2.11
  expect_no_warning(res <- loo_from_taylor(tp))
  expect_identical(res$order, 2L)
  expect_equal(
    unname(res$pointwise[, "p_loo"]),
    c(-0.94 + 1.06, -1.89 + 2.11, 0.3)
  )

  # Silent: no warning, but recorded
  tp$log_cpo_2[2] <- NA
  expect_no_warning(res <- loo_from_taylor(tp, silent = TRUE))
  expect_match(res$warnings, "first-order LOO")
})

test_that("the results have the loo package's layout", {
  res <- loo(fit_lin, ll_lin, X = lin_X, y = lin_y)
  expect_s3_class(res, c("inla_vanilla_loo", "loo"), exact = TRUE)
  expect_identical(rownames(res$estimates), c("elpd_loo", "p_loo", "looic"))
  expect_identical(colnames(res$estimates), c("Estimate", "SE"))
  expect_identical(
    colnames(res$pointwise),
    c("elpd_loo", "p_loo", "looic", "k")
  )
  expect_equal(nrow(res$pointwise), lin_n)
  # One column whose name starts with "elpd", as loo::loo_compare() needs
  expect_length(grep("^elpd", colnames(res$pointwise)), 1)
})

test_that("loo::loo_compare() compares the results", {
  skip_if_not_installed("loo")
  fit_a <- inla_vanilla(
    lpost_pois,
    c(a = 0, b = 0),
    x = pois_x,
    y = pois_y,
    silent = TRUE
  )
  ll_const <- function(theta, x, y) dpois(y, exp(theta[["a"]]), log = TRUE)
  lpost_const <- function(theta, x, y) {
    sum(ll_const(theta, x, y)) + dnorm(theta[["a"]], 0, 2, log = TRUE)
  }
  fit_b <- inla_vanilla(
    lpost_const,
    c(a = 0),
    x = pois_x,
    y = pois_y,
    silent = TRUE
  )
  a <- loo(fit_a, ll_pois, x = pois_x, y = pois_y)
  b <- loo(fit_b, ll_const, x = pois_x, y = pois_y)
  comp <- loo::loo_compare(list(trend = a, constant = b))
  d <- b$pointwise[, "elpd_loo"] - a$pointwise[, "elpd_loo"]
  # loo 2.10 returns a data frame with a `model` column; earlier versions a
  # matrix with the models as row names
  if (is.data.frame(comp)) {
    row_b <- comp[comp$model == "constant", c("elpd_diff", "se_diff")]
  } else {
    row_b <- as.list(comp["constant", c("elpd_diff", "se_diff")])
  }
  expect_equal(row_b$elpd_diff, sum(d), tolerance = 1e-12)
  expect_equal(row_b$se_diff, sqrt(length(d)) * sd(d), tolerance = 1e-12)
})

test_that("other objects go to the loo package's loo() and waic()", {
  skip_if_not_installed("loo")
  set.seed(3)
  ll <- matrix(rnorm(4000, -1, 0.1), 400, 10)
  res <- suppressWarnings(loo(ll))
  expect_s3_class(res, "psis_loo")
  expect_equal(res, suppressWarnings(loo::loo(ll)))
  expect_s3_class(suppressWarnings(waic(ll)), "waic")
})

test_that("warnings from loglik are recorded, and muffled for a silent fit", {
  chatty <- function(theta, X, y) {
    warning("loglik is chatty")
    ll_lin(theta, X, y)
  }
  # fit_lin is silent: nothing is signalled, the warning is recorded once
  expect_silent(res <- loo(fit_lin, chatty, X = lin_X, y = lin_y))
  expect_identical(res$warnings, "loglik is chatty")
  expect_silent(res <- waic(fit_lin, chatty, X = lin_X, y = lin_y))
  expect_identical(res$warnings, "loglik is chatty")
  # A fit that is not silent: signalled at every call, and recorded once
  fit_loud <- fit_lin
  fit_loud$silent <- FALSE
  n_signalled <- 0
  res <- withCallingHandlers(
    suppressMessages(loo(fit_loud, chatty, X = lin_X, y = lin_y)),
    warning = function(w) {
      n_signalled <<- n_signalled + 1
      invokeRestart("muffleWarning")
    }
  )
  expect_identical(n_signalled, 13)
  expect_identical(res$warnings, "loglik is chatty")
})

test_that("the methods are reachable through the loo package's generics", {
  skip_if_not_installed("loo")
  res <- loo::loo(fit_lin, ll_lin, X = lin_X, y = lin_y)
  expect_s3_class(res, "inla_vanilla_loo")
  res <- loo::waic(fit_lin, ll_lin, X = lin_X, y = lin_y)
  expect_s3_class(res, "inla_vanilla_waic")
})

test_that("loglik is checked", {
  expect_error(loo(fit_lin), "`loglik` is missing", class = "inlavanilla_error")
  expect_error(
    loo(fit_lin, 1),
    "`loglik` must be a function",
    class = "inlavanilla_error"
  )
  # The summed log-likelihood instead of its terms
  expect_error(
    loo(
      fit_lin,
      function(theta, X, y) sum(ll_lin(theta, X, y)),
      X = lin_X,
      y = lin_y
    ),
    "single value",
    class = "inlavanilla_error"
  )
  expect_error(
    loo(fit_lin, function(theta) "a"),
    "numeric vector",
    class = "inlavanilla_error"
  )
  # A length that changes with theta
  calls <- 0
  vary <- function(theta) {
    calls <<- calls + 1
    rep(-1, if (calls == 1) 5 else 4)
  }
  expect_error(
    loo(fit_lin, vary),
    "one term per unit",
    class = "inlavanilla_error"
  )
  expect_error(
    loo(fit_lin, function(theta) c(-1, -Inf, -2)),
    "not finite for unit 2",
    class = "inlavanilla_error"
  )
  # An error inside loglik names the parameter value and keeps the cause
  err <- expect_error(
    loo(fit_lin, function(theta) stop("boom")),
    "`loglik` failed",
    class = "inlavanilla_error"
  )
  expect_match(conditionMessage(err), "boom")
  expect_match(conditionMessage(err), "a = ")
})

test_that("unit names label the pointwise results", {
  named <- function(theta, X, y) {
    stats::setNames(ll_lin(theta, X, y), paste0("obs", seq_along(y)))
  }
  res <- loo(fit_lin, named, X = lin_X, y = lin_y)
  expect_identical(rownames(res$pointwise), paste0("obs", seq_len(lin_n)))
  res <- loo(fit_lin, ll_lin, X = lin_X, y = lin_y)
  expect_null(rownames(res$pointwise))
})

test_that("print() shows the estimates and the leverage", {
  testthat::local_reproducible_output(width = 80, crayon = FALSE)
  res <- loo(fit_lin, ll_lin, X = lin_X, y = lin_y)
  expect_snapshot(print(res))
  expect_snapshot(print(waic(fit_lin, ll_lin, X = lin_X, y = lin_y)))
  out <- NULL
  expect_output(out <- withVisible(print(res)))
  expect_false(out$visible)
})
