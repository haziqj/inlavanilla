# Tests of rtmb_latent(), which draws the random effects of an RTMB or TMB
# object from their Gaussian approximation given the parameters.
#
# The reference is a Gaussian in closed form. Most tests need no RTMB: a
# stand-in object has the parts rtmb_latent() reads (fn, env$last.par,
# env$random, env$spHess), with a mode and a precision that are known
# functions of the parameters. The tests with RTMB use random intercepts,
# whose distribution given the parameters is Gaussian with a known mean and
# variance, so the Laplace approximation is exact.
#
# Monte Carlo tolerances come from standard errors: a mean of n draws has SE
# sqrt(v / n), and a variance of n normal draws has SE v sqrt(2 / (n - 1)).
# Each comparison allows z_mc standard errors, and the seeds are fixed.

z_mc <- 4

# A stand-in for an RTMB object with two parameters and four random effects.
# The precision is tridiagonal, so its Cholesky factor is not diagonal.
stand_in <- function(sparse) {
  k <- 4
  tri <- diag(2, k)
  tri[cbind(1:(k - 1), 2:k)] <- -1
  tri[cbind(2:k, 1:(k - 1))] <- -1
  precision <- function(theta) exp(theta[1]) * tri + diag(exp(theta[2]), k)
  centre <- function(theta) solve(precision(theta), seq_len(k))
  env <- new.env()
  env$random <- 3:6
  env$last.par <- c(a = 0, b = 0, u = 0, u = 0, u = 0, u = 0)
  env$spHess <- function(par, random = TRUE) {
    H <- precision(par[1:2])
    if (sparse) Matrix::forceSymmetric(Matrix::Matrix(H, sparse = TRUE)) else H
  }
  list(
    par = c(a = 0, b = 0),
    fn = function(theta) {
      env$last.par[] <- c(theta, centre(theta))
      0
    },
    env = env,
    precision = precision,
    centre = centre
  )
}

expect_gaussian_draws <- function(draws, mean, variance) {
  n <- nrow(draws)
  z_mean <- (colMeans(draws) - mean) / sqrt(variance / n)
  z_var <- (apply(draws, 2, var) - variance) / (variance * sqrt(2 / (n - 1)))
  expect_lt(max(abs(z_mean)), z_mc)
  expect_lt(max(abs(z_var)), z_mc)
}

# A correlation of n normal draws has SE of about 1 / sqrt(n)
expect_correlations <- function(draws, covariance) {
  gap <- cov2cor(cov(draws)) - cov2cor(covariance)
  expect_lt(max(abs(gap)) * sqrt(nrow(draws)), z_mc)
}

test_that("precision_root() maps standard normals to covariance H^-1", {
  set.seed(1)
  n <- 30
  # A banded precision in a shuffled order, so that the permutation matters
  band <- diag(3, n)
  band[cbind(1:(n - 1), 2:n)] <- -1
  band[cbind(2:n, 1:(n - 1))] <- -1
  band[cbind(1:(n - 2), 3:n)] <- 0.5
  band[cbind(3:n, 1:(n - 2))] <- 0.5
  shuffle <- sample(n)
  H <- band[shuffle, shuffle]

  root <- function(H) {
    vapply(
      seq_len(n),
      function(i) {
        precision_root(H, as.numeric(seq_len(n) == i))
      },
      numeric(n)
    )
  }
  dense <- root(H)
  expect_equal(tcrossprod(dense), solve(H), tolerance = 1e-10)

  skip_if_not_installed("Matrix")
  sparse <- root(Matrix::forceSymmetric(Matrix::Matrix(H, sparse = TRUE)))
  expect_equal(tcrossprod(sparse), solve(H), tolerance = 1e-10)
})

test_that("the draws follow N(mode, H^-1) with a base matrix", {
  obj <- stand_in(sparse = FALSE)
  fun <- rtmb_latent(obj)
  theta <- c(a = 0.3, b = -0.5)
  set.seed(2)
  draws <- t(replicate(2e4, fun(theta)))
  expect_identical(colnames(draws), sprintf("u[%d]", 1:4))
  expect_gaussian_draws(
    draws,
    obj$centre(theta),
    diag(solve(obj$precision(theta)))
  )
  expect_correlations(draws, solve(obj$precision(theta)))
})

test_that("the draws follow N(mode, H^-1) with a sparse matrix", {
  skip_if_not_installed("Matrix")
  obj <- stand_in(sparse = TRUE)
  fun <- rtmb_latent(obj)
  theta <- c(a = -0.4, b = 0.2)
  set.seed(3)
  draws <- t(replicate(2e4, fun(theta)))
  expect_gaussian_draws(
    draws,
    obj$centre(theta),
    diag(solve(obj$precision(theta)))
  )
  expect_correlations(draws, solve(obj$precision(theta)))
})

test_that("latent_names() numbers the elements of a vector parameter", {
  expect_identical(
    latent_names(c("alpha", "u", "u", "v", "u"), 5),
    c("alpha", "u[1]", "u[2]", "v", "u[3]")
  )
  expect_identical(latent_names(NULL, 2), c("latent[1]", "latent[2]"))
})

test_that("rtmb_latent() rejects what is not an object with random effects", {
  expect_error(rtmb_latent(function(theta) 0), class = "inlavanilla_error")
  expect_error(rtmb_latent(list(fn = identity)), "MakeADFun")
  obj <- stand_in(sparse = FALSE)
  obj$env$random <- NULL
  expect_error(rtmb_latent(obj), "no random effects")
})

test_that("a draw that cannot be made stops with the package's error", {
  obj <- stand_in(sparse = FALSE)
  fun <- rtmb_latent(obj)
  expect_error(fun(c(1, 2, 3)), "same order", class = "inlavanilla_error")

  failed <- stand_in(sparse = FALSE)
  failed$fn <- function(theta) NaN
  expect_error(
    rtmb_latent(failed)(c(0, 0)),
    "no Gaussian approximation",
    class = "inlavanilla_error"
  )

  indefinite <- stand_in(sparse = FALSE)
  indefinite$env$spHess <- function(par, random = TRUE) -diag(4)
  expect_error(
    rtmb_latent(indefinite)(c(0, 0)),
    "positive definite",
    class = "inlavanilla_error"
  )
})

# Random intercepts in RTMB ----------------------------------------------------
#
# y_ij = mu + u_i + e_ij with u_i ~ N(0, tau^2) and e_ij ~ N(0, sigma^2), five
# observations per group. Given the parameters, u_i is Gaussian with variance
# v = 1 / (1 / tau^2 + 5 / sigma^2) and mean v * sum_j (y_ij - mu) / sigma^2.

intercepts <- function() {
  set.seed(4)
  group <- rep(1:6, each = 5)
  y <- 1 + rnorm(6)[group] + rnorm(30, sd = 0.7)
  nll <- function(p) {
    -sum(RTMB::dnorm(y, p$mu + p$u[group], exp(p$log_sigma), log = TRUE)) -
      sum(RTMB::dnorm(p$u, 0, exp(p$log_tau), log = TRUE)) -
      RTMB::dnorm(p$mu, 0, 10, log = TRUE) -
      RTMB::dnorm(p$log_tau, 0, 2, log = TRUE) -
      RTMB::dnorm(p$log_sigma, 0, 2, log = TRUE)
  }
  obj <- RTMB::MakeADFun(
    nll,
    list(mu = 0, log_tau = 0, log_sigma = 0, u = numeric(6)),
    random = "u",
    silent = TRUE
  )
  closed_form <- function(theta) {
    theta <- unname(theta)
    tau <- exp(theta[2])
    sigma <- exp(theta[3])
    v <- 1 / (1 / tau^2 + 5 / sigma^2)
    list(
      mean = unname(v * rowsum(y - theta[1], group)[, 1] / sigma^2),
      variance = rep(v, 6)
    )
  }
  list(obj = obj, closed_form = closed_form)
}

test_that("an RTMB object gives the closed-form Gaussian at each draw", {
  skip_if_not_installed("RTMB")
  skip_if_not_installed("Matrix")
  model <- intercepts()
  fun <- rtmb_latent(model$obj)

  # A first call at other parameters, whose precision is 9 times as large.
  # The draws that follow must use their own factor, not that of this call.
  fun(c(0.5, log(0.5), log(0.3)))
  theta <- c(1.2, log(1.5), log(0.9))
  set.seed(5)
  draws <- t(replicate(4000, fun(theta)))
  expect_identical(colnames(draws), sprintf("u[%d]", 1:6))
  ref <- model$closed_form(theta)
  expect_gaussian_draws(draws, ref$mean, ref$variance)
})

test_that("sampling() mixes the draws over the parameters of a fit", {
  skip_if_not_installed("RTMB")
  skip_if_not_installed("Matrix")
  model <- intercepts()
  obj <- model$obj
  fit <- inla_vanilla(
    function(theta) -obj$fn(theta),
    obj$par,
    grad = function(theta) -as.vector(obj$gr(theta)),
    silent = TRUE
  )
  set.seed(6)
  u <- sampling(fit, nsamp = 4000, fun = rtmb_latent(obj))
  expect_identical(dim(u), c(4000L, 6L))

  # The same parameter draws through the closed form: the mixture has the
  # mean of the conditional means, and the law of total variance gives its
  # variance.
  set.seed(6)
  parts <- sampling(fit, nsamp = 4000, fun = function(theta) {
    ref <- model$closed_form(theta)
    c(ref$mean, ref$variance)
  })
  total_mean <- colMeans(parts[, 1:6])
  total_var <- colMeans(parts[, 7:12]) + apply(parts[, 1:6], 2, var)
  expect_lt(max(abs(colMeans(u) - total_mean) / sqrt(total_var / 4000)), z_mc)
  expect_equal(unname(apply(u, 2, var)), total_var, tolerance = 0.1)
})
