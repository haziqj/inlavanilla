# Tests of sampling(), the Gaussian-copula draws from a fit.
#
# sampling() is defined by the fit it draws from: column j follows the
# skew-normal marginal SN(xi_j, omega_j, alpha_j) in fit$sn, and the columns
# are coupled by a Gaussian copula with correlation R = cov2cor(fit$cov). The
# references are therefore the fit's own skew-normals, evaluated here by
# quadrature of the skew-normal density (written in this file, independent of
# the package's psnorm() and qsnorm()), and R. The end-to-end checks against
# exact posteriors use an exactly Gaussian posterior (Bayesian linear
# regression with known noise) and the gamma-Poisson model on log(lambda).
#
# Monte Carlo tolerances come from standard errors, not guesses:
# * an empirical p-quantile from n draws has SE sqrt(p (1 - p) / n) / f(q_p),
#   the binomial SE of the empirical CDF divided by the density at q_p;
# * a Pearson correlation of bivariate normal draws has
#   Var(atanh(r)) = 1 / (n - 3);
# * a Spearman correlation of bivariate normal draws has
#   Var(atanh(r_s)) ~ (1 + r_s^2 / 2) / (n - 3) (Bonett and Wright, 2000,
#   Psychometrika 65, 23-28);
# and each comparison allows `z_mc` standard errors. With about a hundred
# comparisons in the file, z_mc = 4 keeps the chance of a false failure for
# any seed below about 1%. The seeds are fixed, so the results are
# deterministic; the observed largest standardised deviations are noted.

z_mc <- 4
probs <- c(0.01, 0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975, 0.99)
n_big <- 2e4

# Skew-normal references -------------------------------------------------------
#
# The standardised density is 2 phi(z) Phi(alpha z). Tail probabilities are
# integrated on the side of the smaller tail, so they keep full relative
# accuracy far out.

sn_std_density <- function(z, alpha) 2 * dnorm(z) * pnorm(alpha * z)

sn_density_ref <- function(x, xi, omega, alpha) {
  sn_std_density((x - xi) / omega, alpha) / omega
}

# Lower and upper tail probabilities of SN(0, 1, alpha) at a scalar z
sn_std_tails <- function(z, alpha) {
  f <- function(t) sn_std_density(t, alpha)
  if (z <= 0) {
    lower <- integrate(f, -Inf, z, rel.tol = 1e-11, abs.tol = 0)$value
    c(lower = lower, upper = 1 - lower)
  } else {
    upper <- integrate(f, z, Inf, rel.tol = 1e-11, abs.tol = 0)$value
    c(lower = 1 - upper, upper = upper)
  }
}

sn_cdf_ref <- function(q, xi, omega, alpha) {
  vapply(
    (q - xi) / omega,
    function(z) sn_std_tails(z, alpha)[["lower"]],
    numeric(1)
  )
}

sn_quantile_ref <- function(p, xi, omega, alpha) {
  vapply(
    p,
    function(pp) {
      root <- uniroot(
        function(z) sn_std_tails(z, alpha)[["lower"]] - pp,
        c(-12, 12),
        tol = 1e-13
      )
      xi + omega * root$root
    },
    numeric(1)
  )
}

# Normal scores qnorm(F(x)) of skew-normal values x: exact at knots spanning
# the range of x (from the smaller tail, through qnorm's matching tail), and
# a cubic spline between them. The map is smooth and nearly linear, so 301
# knots reproduce it to about 1e-10 (checked against psnorm() while writing
# these tests), far below both Monte Carlo error and the map's own tolerance.
sn_normal_score_ref <- function(x, xi, omega, alpha, n_knots = 301) {
  z <- (x - xi) / omega
  zk <- seq(min(z) - 1e-3, max(z) + 1e-3, length.out = n_knots)
  sk <- vapply(
    zk,
    function(zz) {
      tl <- sn_std_tails(zz, alpha)
      if (zz <= 0) {
        qnorm(tl[["lower"]])
      } else {
        qnorm(tl[["upper"]], lower.tail = FALSE)
      }
    },
    numeric(1)
  )
  splinefun(zk, sk)(z)
}

# SD of SN(xi, omega, alpha): omega sqrt(1 - 2 delta^2 / pi), with
# delta = alpha / sqrt(1 + alpha^2)
sn_sd_ref <- function(omega, alpha) {
  omega * sqrt(1 - 2 * alpha^2 / ((1 + alpha^2) * pi))
}

sn_par <- function(fit, j) {
  list(
    xi = fit$sn[j, "xi"],
    omega = fit$sn[j, "omega"],
    alpha = fit$sn[j, "alpha"]
  )
}

# Monte Carlo SE of the empirical p-quantile of n draws with density f there
quantile_se <- function(p, n, f) sqrt(p * (1 - p) / n) / f

# Models -----------------------------------------------------------------------

# A skewed, correlated posterior: Poisson regression with small counts, an
# uncentred covariate (so b0 and b1 are strongly correlated) and a binary one,
# with N(0, 2^2) priors on the coefficients. The skew-normal alphas are about
# -1.3, 1.2 and 0.6, and the correlations about -0.82, -0.08 and -0.39.
pois_x <- c(0.2, 0.5, 0.9, 1.3, 1.6, 2.0, 2.4, 2.9, 3.3, 3.8)
pois_d <- c(1, 0, 1, 1, 0, 0, 1, 0, 1, 1)
pois_y <- c(0, 1, 0, 2, 1, 3, 2, 4, 6, 5)
pois_X <- cbind(1, pois_x, pois_d, deparse.level = 0)
pois_logpost <- function(theta) {
  eta <- drop(pois_X %*% theta)
  sum(dpois(pois_y, exp(eta), log = TRUE)) +
    sum(dnorm(theta, 0, 2, log = TRUE))
}
pois_grad <- function(theta) {
  eta <- drop(pois_X %*% theta)
  drop(crossprod(pois_X, pois_y - exp(eta))) - theta / 4
}
pois_fit <- inla_vanilla(
  pois_logpost,
  c(b0 = 0, b1 = 0, b2 = 0),
  grad = pois_grad,
  verbose = FALSE
)

# An exactly Gaussian posterior: y = X b + e, e ~ N(0, sigma^2) with sigma
# known, b ~ N(0, tau^2 I). Then b | y ~ N(mu, V) with
# V = (X'X / sigma^2 + I / tau^2)^-1 and mu = V X'y / sigma^2.
lm_x <- c(0.3, 0.8, 1.1, 1.7, 2.2, 2.6, 3.1, 3.5, 4.0, 4.6, 5.2, 5.9)
lm_d <- c(1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0)
lm_X <- cbind(1, lm_x, lm_d, deparse.level = 0)
lm_sigma <- 0.7
lm_tau <- 3
set.seed(3)
lm_y <- drop(lm_X %*% c(1, 0.5, -0.8)) + rnorm(12, 0, lm_sigma)
lm_V <- solve(crossprod(lm_X) / lm_sigma^2 + diag(3) / lm_tau^2)
lm_mu <- drop(lm_V %*% crossprod(lm_X, lm_y)) / lm_sigma^2
lm_logpost <- function(theta) {
  sum(dnorm(lm_y, drop(lm_X %*% theta), lm_sigma, log = TRUE)) +
    sum(dnorm(theta, 0, lm_tau, log = TRUE))
}
lm_grad <- function(theta) {
  drop(crossprod(lm_X, lm_y - drop(lm_X %*% theta))) /
    lm_sigma^2 -
    theta / lm_tau^2
}
lm_fit <- inla_vanilla(
  lm_logpost,
  c(a = 0, b = 0, c = 0),
  grad = lm_grad,
  verbose = FALSE
)

# Gamma-Poisson on theta = log(lambda): y_i ~ Poisson(lambda),
# lambda ~ Gamma(1, 1), Jacobian included. Exactly lambda | y ~ Gamma(S, R)
# with S = 1 + sum(y) = 5 and R = 1 + n = 6; theta has skewness about -0.47.
gp_y <- c(0, 1, 0, 2, 1)
gp_S <- 1 + sum(gp_y)
gp_R <- 1 + length(gp_y)
gp_logpost <- function(theta) {
  lambda <- exp(theta[[1]])
  sum(dpois(gp_y, lambda, log = TRUE)) +
    dgamma(lambda, shape = 1, rate = 1, log = TRUE) +
    theta[[1]]
}
gp_grad <- function(theta) gp_S - gp_R * exp(theta[[1]])
gp_fit <- inla_vanilla(
  gp_logpost,
  c(log_lambda = 0),
  grad = gp_grad,
  verbose = FALSE
)

# Large draws shared by the distributional tests
set.seed(101)
pois_draws <- sampling(pois_fit, nsamp = n_big)
set.seed(102)
lm_draws <- sampling(lm_fit, nsamp = n_big)

# (a) Shape, names and reproducibility ------------------------------------------

test_that("sampling() returns an nsamp x m matrix named by the parameters", {
  set.seed(1)
  draws <- sampling(pois_fit, nsamp = 250)
  expect_true(is.matrix(draws))
  expect_type(draws, "double")
  expect_identical(dim(draws), c(250L, 3L))
  expect_identical(colnames(draws), c("b0", "b1", "b2"))
  expect_null(rownames(draws))
  expect_true(all(is.finite(draws)))

  # one parameter still gives a matrix
  draws_1 <- sampling(gp_fit, nsamp = 40)
  expect_identical(dim(draws_1), c(40L, 1L))
  expect_identical(colnames(draws_1), "log_lambda")

  # a single draw is a 1 x m matrix, not a vector
  expect_identical(dim(sampling(pois_fit, nsamp = 1)), c(1L, 3L))
})

test_that("sampling() draws 1000 by default and accepts a double nsamp", {
  expect_identical(nrow(sampling(pois_fit)), 1000L)
  set.seed(5)
  from_double <- sampling(pois_fit, nsamp = 300)
  set.seed(5)
  from_integer <- sampling(pois_fit, nsamp = 300L)
  expect_identical(from_double, from_integer)
  expect_identical(nrow(sampling(pois_fit, nsamp = 1e3)), 1000L)
})

test_that("sampling() is reproducible under set.seed()", {
  set.seed(11)
  first <- sampling(pois_fit, nsamp = 100)
  set.seed(11)
  second <- sampling(pois_fit, nsamp = 100)
  expect_identical(first, second)

  set.seed(12)
  other <- sampling(pois_fit, nsamp = 100)
  expect_false(isTRUE(all.equal(first, other)))
})

test_that("sampling() has no method for objects that are not fits", {
  expect_error(
    sampling(list(sn = pois_fit$sn, cov = pois_fit$cov)),
    "no applicable method"
  )
})

test_that("sampling() never calls the log-posterior", {
  # Documented: draws never call logpost, so a large nsamp is cheap
  n_calls <- 0
  counting_logpost <- function(theta) {
    n_calls <<- n_calls + 1
    gp_logpost(theta)
  }
  fit <- inla_vanilla(counting_logpost, c(log_lambda = 0), verbose = FALSE)
  n_calls <- 0
  set.seed(7)
  sampling(fit, nsamp = 5000)
  sampling(fit, nsamp = 50, fun = function(theta) exp(theta))
  expect_identical(n_calls, 0)
})

# (b) Marginals -----------------------------------------------------------------

test_that("the quadrature references agree with the fit's quantiles and moments", {
  # Anchors the references used below: fit$summary and qsnorm() on fit$sn
  # against the skew-normal quantiles solved here by quadrature, and the
  # summary's Mean and SD against the moments by quadrature. Observed: about
  # 5e-14 (quantiles) and 1e-16 (moments); the tolerance is set by the
  # quadrature's own rel.tol of 1e-11, not by the observed error.
  tab_probs <- c(0.025, 0.5, 0.975)
  for (j in rownames(pois_fit$sn)) {
    sn <- sn_par(pois_fit, j)
    q_ref <- sn_quantile_ref(probs, sn$xi, sn$omega, sn$alpha)
    expect_equal(
      qsnorm(probs, sn$xi, sn$omega, sn$alpha),
      q_ref,
      tolerance = 1e-10
    )
    expect_equal(
      unname(pois_fit$summary[j, c("2.5%", "50%", "97.5%")]),
      sn_quantile_ref(tab_probs, sn$xi, sn$omega, sn$alpha),
      tolerance = 1e-10
    )

    dens <- function(x) sn_density_ref(x, sn$xi, sn$omega, sn$alpha)
    mean_ref <- integrate(
      function(x) x * dens(x),
      -Inf,
      Inf,
      rel.tol = 1e-12
    )$value
    var_ref <- integrate(
      function(x) (x - mean_ref)^2 * dens(x),
      -Inf,
      Inf,
      rel.tol = 1e-12
    )$value
    expect_equal(pois_fit$summary[j, "Mean"], mean_ref, tolerance = 1e-10)
    expect_equal(pois_fit$summary[j, "SD"], sqrt(var_ref), tolerance = 1e-10)
  }
})

test_that("each column's empirical quantiles match its skew-normal marginal", {
  # Observed: largest deviation 2.3 Monte Carlo SEs over the 27 quantiles.
  for (j in rownames(pois_fit$sn)) {
    sn <- sn_par(pois_fit, j)
    q_ref <- sn_quantile_ref(probs, sn$xi, sn$omega, sn$alpha)
    se <- quantile_se(
      probs,
      n_big,
      sn_density_ref(q_ref, sn$xi, sn$omega, sn$alpha)
    )
    q_emp <- unname(quantile(pois_draws[, j], probs))
    expect_lt(max(abs(q_emp - q_ref) / se), z_mc)
  }
})

test_that("the draws reproduce fit$summary's quantiles, mean and SD", {
  # The summary's mean and SD are the skew-normal moments; the SE of a
  # sample SD is SD sqrt((kurtosis - 1) / (4 n)), with the skew-normal
  # kurtosis 3 + 2 (pi - 3) (b delta)^4 / (1 - (b delta)^2)^2, b = sqrt(2/pi).
  # Observed: largest deviations 2.1 SEs (quantiles) and 1.6 SEs (moments).
  tab_probs <- c(0.025, 0.5, 0.975)
  for (j in rownames(pois_fit$sn)) {
    sn <- sn_par(pois_fit, j)
    tab <- pois_fit$summary[j, ]
    q_tab <- unname(tab[c("2.5%", "50%", "97.5%")])
    se_q <- quantile_se(
      tab_probs,
      n_big,
      sn_density_ref(q_tab, sn$xi, sn$omega, sn$alpha)
    )
    q_emp <- unname(quantile(pois_draws[, j], tab_probs))
    expect_lt(max(abs(q_emp - q_tab) / se_q), z_mc)

    bd <- sqrt(2 / pi) * sn$alpha / sqrt(1 + sn$alpha^2)
    kurt <- 3 + 2 * (pi - 3) * bd^4 / (1 - bd^2)^2
    se_mean <- tab[["SD"]] / sqrt(n_big)
    se_sd <- tab[["SD"]] * sqrt((kurt - 1) / (4 * n_big))
    expect_lt(abs(mean(pois_draws[, j]) - tab[["Mean"]]) / se_mean, z_mc)
    expect_lt(abs(sd(pois_draws[, j]) - tab[["SD"]]) / se_sd, z_mc)
  }
})

test_that("each column passes a Kolmogorov-Smirnov test against its marginal", {
  # Normal scores through the reference CDF should be standard normal.
  # Observed p-values: 0.18, 0.87 and 0.26.
  for (j in rownames(pois_fit$sn)) {
    sn <- sn_par(pois_fit, j)
    s <- sn_normal_score_ref(pois_draws[, j], sn$xi, sn$omega, sn$alpha)
    expect_gt(ks.test(s, "pnorm")$p.value, 1e-3)
  }
})

test_that("a one-parameter draw is the skew-normal quantile of a normal", {
  # With m = 1 the copula is the identity, so draw i is Q(Phi(z_i)) for the
  # i-th standard normal z_i that rnorm() gives under the same seed. This
  # checks the quantile map, row by row, to far below Monte Carlo error.
  # Observed: largest error 4.85e-8 posterior SDs, which is the spline's
  # largest error anywhere on |z| <= 7 for this alpha (-1.9), between the
  # knots near z = -0.55. sn_score_map() documents 1e-6 SDs for any alpha
  # (5e-7 at |alpha| = 20); halving its knots would give about 8e-7 here.
  sn <- sn_par(gp_fit, 1)
  set.seed(21)
  z <- rnorm(2000)
  set.seed(21)
  draws <- sampling(gp_fit, nsamp = 2000)[, 1]
  s <- sn_normal_score_ref(draws, sn$xi, sn$omega, sn$alpha)
  # the error in z, mapped to the parameter scale through dQ/dz
  dq_dz <- dnorm(z) / sn_density_ref(draws, sn$xi, sn$omega, sn$alpha)
  sd_sn <- sn_sd_ref(sn$omega, sn$alpha)
  expect_lt(max(abs(s - z) * dq_dz) / sd_sn, 1.5e-7)
})

test_that("each column is its skew-normal quantile of a linear map of iid normals", {
  # A deterministic version of the copula checks, free of Monte Carlo error.
  # The normal scores S = qnorm(F_j(draws)) of the draws should be Z A for
  # the iid normals Z of the same seed and some square root A of R. Any
  # square root passes, so the test does not depend on the use of chol().
  # It does assume that sampling() takes its normals as one rnorm(nsamp * m)
  # call, filled column by column, as R/sampling.R does now; if that ever
  # changes, rebuild Z here to match. Observed over five seeds: largest
  # residual 1.5e-8 posterior SDs (the spline map's error at alpha = -1.35),
  # and A'A within 5e-9 of R (the same smooth map error, projected onto the
  # regression). A fitted A is chol(R) to six decimals.
  n <- 2000
  R <- cov2cor(pois_fit$cov)
  set.seed(61)
  z_iid <- matrix(rnorm(n * 3), n, 3)
  set.seed(61)
  draws <- sampling(pois_fit, nsamp = n)
  s <- draws
  dq_dz_sd <- draws
  for (j in 1:3) {
    sn <- sn_par(pois_fit, j)
    s[, j] <- sn_normal_score_ref(draws[, j], sn$xi, sn$omega, sn$alpha)
    dq_dz_sd[, j] <- dnorm(s[, j]) /
      sn_density_ref(draws[, j], sn$xi, sn$omega, sn$alpha) /
      sn_sd_ref(sn$omega, sn$alpha)
  }
  A <- qr.solve(z_iid, s)
  resid_sd <- abs(s - z_iid %*% A) * dq_dz_sd
  expect_lt(max(resid_sd), 5e-8)
  expect_lt(max(abs(crossprod(A) - R)), 1.5e-8)
})

# (c) Dependence ----------------------------------------------------------------

test_that("Spearman correlations match those of the Gaussian copula", {
  # For a Gaussian copula with correlation r, Spearman's rho is
  # (6 / pi) asin(r / 2), whatever the marginals. Observed: largest deviation
  # 0.8 SEs on the Fisher scale.
  R <- cov2cor(pois_fit$cov)
  rho_s <- (6 / pi) * asin(R / 2)
  emp <- cor(pois_draws, method = "spearman")
  lower <- lower.tri(R)
  se <- sqrt((1 + rho_s[lower]^2 / 2) / (n_big - 3))
  dev <- abs(atanh(emp[lower]) - atanh(rho_s[lower])) / se
  expect_lt(max(dev), z_mc)
  # the copula carries real dependence here, so the test has teeth
  expect_gt(max(abs(R[lower])), 0.8)
})

test_that("normal scores of the draws are N(0, R)", {
  # Mapping each column back through its marginal CDF recovers the latent
  # normals, whose Pearson correlations should be R and whose squared
  # Mahalanobis distances under R should be chi-squared with m degrees of
  # freedom. Observed: largest correlation deviation 0.4 SEs on the Fisher
  # scale; the KS p-value of the distances is 0.089.
  R <- cov2cor(pois_fit$cov)
  s <- vapply(
    rownames(pois_fit$sn),
    function(j) {
      sn <- sn_par(pois_fit, j)
      sn_normal_score_ref(pois_draws[, j], sn$xi, sn$omega, sn$alpha)
    },
    numeric(n_big)
  )
  lower <- lower.tri(R)
  dev <- abs(atanh(cor(s)[lower]) - atanh(R[lower])) * sqrt(n_big - 3)
  expect_lt(max(dev), z_mc)

  d2 <- rowSums((s %*% solve(R)) * s)
  expect_gt(ks.test(d2, "pchisq", df = 3)$p.value, 1e-3)
})

test_that("for a Gaussian posterior the draws are the exact normal", {
  # The Laplace step is exact here, so fit$mode is mu, fit$cov is V and the
  # skew-normals are (near) normal. Observed: mode error 4.2e-7 posterior SDs
  # (nlminb's stopping point; up to 4.4e-6 over 100 other data sets), cov
  # error 3.9e-12 relative (up to 1.2e-11), and |alpha| at most 0.0053 on
  # every data set, where the skew-normal fit stops on a shape that is flat
  # to sixth order.
  sd_exact <- sqrt(diag(lm_V))
  expect_lt(max(abs(lm_fit$mode - lm_mu) / sd_exact), 1e-5)
  expect_lt(max(abs(lm_fit$cov - lm_V)) / max(abs(lm_V)), 5e-11)
  expect_lt(max(abs(lm_fit$sn[, "alpha"])), 0.015)

  # Pearson correlations against R and against the exact correlations.
  # Observed: largest deviation 0.3 SEs on the Fisher scale.
  R <- cov2cor(lm_fit$cov)
  R_exact <- cov2cor(lm_V)
  lower <- lower.tri(R)
  emp <- cor(lm_draws)
  expect_lt(
    max(abs(atanh(emp[lower]) - atanh(R[lower])) * sqrt(n_big - 3)),
    z_mc
  )
  expect_lt(
    max(abs(atanh(emp[lower]) - atanh(R_exact[lower])) * sqrt(n_big - 3)),
    z_mc
  )
  expect_gt(max(abs(R[lower])), 0.75)

  # Means and covariances against the exact posterior. The SE of a sample
  # covariance of normal draws is sqrt((V_jj V_kk + V_jk^2) / n). Observed:
  # largest deviations 1.1 SEs (means) and 0.9 SEs (covariances).
  expect_lt(
    max(abs(colMeans(lm_draws) - lm_mu) / (sd_exact / sqrt(n_big))),
    z_mc
  )
  se_cov <- sqrt((outer(diag(lm_V), diag(lm_V)) + lm_V^2) / n_big)
  expect_lt(max(abs(cov(lm_draws) - lm_V) / se_cov), z_mc)

  # Quantiles against the exact normal quantiles. Observed: largest deviation
  # 2.4 SEs.
  for (j in 1:3) {
    q_exact <- qnorm(probs, lm_mu[j], sd_exact[j])
    se <- quantile_se(probs, n_big, dnorm(q_exact, lm_mu[j], sd_exact[j]))
    q_emp <- unname(quantile(lm_draws[, j], probs))
    expect_lt(max(abs(q_emp - q_exact) / se), z_mc)
  }

  # Squared Mahalanobis distances under the exact posterior are chi-squared
  # with 3 degrees of freedom. Observed p-value 0.39.
  centred <- sweep(lm_draws, 2, lm_mu)
  d2 <- rowSums((centred %*% solve(lm_V)) * centred)
  expect_gt(ks.test(d2, "pchisq", df = 3)$p.value, 1e-3)
})

# (d) fun -----------------------------------------------------------------------

test_that("fun receives each draw once, in order, as a named numeric vector", {
  seen <- list()
  record <- function(theta) {
    seen[[length(seen) + 1L]] <<- theta
    theta[["b1"]]
  }
  set.seed(31)
  draws <- sampling(pois_fit, nsamp = 5)
  set.seed(31)
  out <- sampling(pois_fit, nsamp = 5, fun = record)
  # one call per draw, so a fun with side effects (or a costly one) is not
  # run more often than asked
  expect_length(seen, 5L)
  for (theta in seen) {
    expect_type(theta, "double")
    expect_null(dim(theta))
    expect_named(theta, c("b0", "b1", "b2"))
  }
  # call i receives row i of the draws
  expect_identical(do.call(rbind, seen), draws)

  # with one parameter, too
  seen <- list()
  sampling(gp_fit, nsamp = 3, fun = function(theta) {
    seen[[length(seen) + 1L]] <<- theta
    exp(theta[["log_lambda"]])
  })
  for (theta in seen) {
    expect_type(theta, "double")
    expect_named(theta, "log_lambda")
  }
})

test_that("fun's results equal fun applied to the draws under the same seed", {
  f <- function(theta) {
    c(
      rate_ratio = exp(theta[["b2"]]),
      slope_sq = theta[["b1"]]^2,
      total = sum(theta)
    )
  }
  set.seed(41)
  draws <- sampling(pois_fit, nsamp = 60)
  set.seed(41)
  out <- sampling(pois_fit, nsamp = 60, fun = f)
  expected <- do.call(
    rbind,
    lapply(seq_len(nrow(draws)), function(i) f(draws[i, ]))
  )
  expect_equal(out, expected)
})

test_that("a scalar fun gives an nsamp x 1 matrix", {
  set.seed(42)
  draws <- sampling(pois_fit, nsamp = 30)
  set.seed(42)
  out <- sampling(pois_fit, nsamp = 30, fun = function(theta) theta[["b1"]])
  expect_true(is.matrix(out))
  expect_identical(dim(out), c(30L, 1L))
  expect_null(colnames(out))
  expect_equal(out[, 1], unname(draws[, "b1"]))

  # a named scalar names its column
  set.seed(42)
  named <- sampling(pois_fit, nsamp = 30, fun = function(theta) {
    c(slope = theta[["b1"]])
  })
  expect_identical(colnames(named), "slope")
  expect_equal(named[, "slope"], unname(draws[, "b1"]))

  # also for a single draw
  expect_identical(
    dim(sampling(pois_fit, nsamp = 1, fun = function(theta) theta[[1]])),
    c(1L, 1L)
  )
})

test_that("a vector-valued fun gives named columns", {
  set.seed(43)
  draws <- sampling(pois_fit, nsamp = 25)
  set.seed(43)
  out <- sampling(pois_fit, nsamp = 25, fun = function(theta) {
    c(mean_at_0 = exp(theta[["b0"]]), mean_at_1 = exp(sum(theta[1:2])))
  })
  expect_identical(dim(out), c(25L, 2L))
  expect_identical(colnames(out), c("mean_at_0", "mean_at_1"))
  expect_equal(out[, "mean_at_0"], unname(exp(draws[, "b0"])))
  expect_equal(out[, "mean_at_1"], unname(exp(draws[, "b0"] + draws[, "b1"])))

  # a matrix-valued fun (say a covariance matrix) is flattened column-major
  set.seed(43)
  flat <- sampling(pois_fit, nsamp = 25, fun = function(theta) {
    outer(theta[1:2], theta[1:2])
  })
  expect_identical(dim(flat), c(25L, 4L))
  expect_equal(flat[, 2], unname(draws[, "b1"] * draws[, "b0"]))
  expect_equal(flat[, 4], unname(draws[, "b1"]^2))
})

test_that("extra arguments in ... are passed on to fun", {
  shift_scale <- function(theta, scale, shift = 0) {
    theta[["b1"]] * scale + shift
  }
  set.seed(44)
  draws <- sampling(pois_fit, nsamp = 20)
  set.seed(44)
  out <- sampling(pois_fit, nsamp = 20, fun = shift_scale, scale = 2, shift = 1)
  expect_equal(out[, 1], unname(2 * draws[, "b1"] + 1))

  # positionally, as sampling(fit, nsamp, fun, ...)
  set.seed(44)
  out_pos <- sampling(pois_fit, 20, shift_scale, scale = 2)
  expect_equal(out_pos[, 1], unname(2 * draws[, "b1"]))

  # a missing extra argument is fun's own error
  expect_error(sampling(pois_fit, 5, fun = shift_scale), "scale")
})

test_that("a logical fun gives 0/1 indicators that estimate probabilities", {
  # P(b1 > 0.5) from the draws matches 1 - F(0.5) of b1's skew-normal within
  # binomial error. Observed: 0.4 SEs.
  sn <- sn_par(pois_fit, "b1")
  set.seed(45)
  ind <- sampling(pois_fit, nsamp = n_big, fun = function(theta) {
    theta[["b1"]] > 0.5
  })
  expect_true(all(ind %in% c(0, 1)))
  p_ref <- 1 - sn_cdf_ref(0.5, sn$xi, sn$omega, sn$alpha)
  se <- sqrt(p_ref * (1 - p_ref) / n_big)
  expect_lt(abs(mean(ind) - p_ref) / se, z_mc)
})

test_that("a fun returning a non-numeric or empty value errors clearly", {
  expect_error(
    sampling(pois_fit, 10, fun = function(theta) "a"),
    "must return a numeric vector"
  )
  expect_error(
    sampling(pois_fit, 10, fun = function(theta) list(theta[[1]])),
    "must return a numeric vector"
  )
  expect_error(
    sampling(pois_fit, 10, fun = function(theta) NULL),
    "must return a numeric vector"
  )
  expect_error(
    sampling(pois_fit, 10, fun = function(theta) numeric(0)),
    "empty"
  )
  # a length that changes between draws is the package's own clear error,
  # not silent recycling or vapply()'s error (the message wraps, hence \\s+)
  set.seed(46)
  expect_error(
    sampling(pois_fit, 50, fun = function(theta) theta[theta > 0]),
    "same\\s+length",
    class = "inlavanilla_error"
  )
})

# (e) Input validation ----------------------------------------------------------

test_that("an invalid nsamp is an error", {
  bad <- list(0, -5, NA, NA_real_, NaN, Inf, c(10, 20), "100", TRUE, NULL)
  for (nsamp in bad) {
    expect_error(
      sampling(pois_fit, nsamp = nsamp),
      "must be a single positive number"
    )
  }
})

test_that("a fun that is not a function is an error", {
  expect_error(sampling(pois_fit, 10, fun = "exp"), "must be a function")
  expect_error(sampling(pois_fit, 10, fun = 1), "must be a function")
})

# (f) A derived natural-scale quantity ------------------------------------------

# lambda | y ~ Gamma(5, 6) exactly. The fit's skew-normal for theta, pushed
# through exp(), is the distribution that sampling(fun = exp) draws from, so
# the checks split the two error sources: that distribution against the exact
# gamma (approximation error, no noise), and the draws against it (Monte Carlo
# error only); then the draws against the exact gamma, within both, and
# against the Laplace lognormal.
gp_sd_exact <- sqrt(gp_S) / gp_R
gp_q_exact <- qgamma(probs, shape = gp_S, rate = gp_R)
gp_sn <- sn_par(gp_fit, 1)
gp_q_sn <- exp(sn_quantile_ref(probs, gp_sn$xi, gp_sn$omega, gp_sn$alpha))
gp_q_laplace <- exp(qnorm(probs, gp_fit$mode, sqrt(gp_fit$cov[1, 1])))
# Quantile errors allowed for the skew-normal approximation, in exact SDs
gp_tol_sn <- ifelse(probs > 0.98, 0.1, 0.03)

# E exp(t X) for X ~ SN(xi, omega, alpha) by quadrature on the theta scale
gp_sn_exp_moment <- function(t) {
  integrate(
    function(x) {
      exp(t * x) * sn_density_ref(x, gp_sn$xi, gp_sn$omega, gp_sn$alpha)
    },
    gp_sn$xi - 40 * gp_sn$omega,
    gp_sn$xi + 20 * gp_sn$omega,
    rel.tol = 1e-12
  )$value
}
gp_moments_sn <- vapply(1:4, gp_sn_exp_moment, numeric(1))
gp_mean_sn <- gp_moments_sn[1]
gp_sd_sn <- sqrt(gp_moments_sn[2] - gp_mean_sn^2)
# Kurtosis of exp(theta) under the skew-normal, from the raw moments: 4.52,
# against 3 + 6 / S = 4.2 for the exact gamma. The Monte Carlo SEs of the
# draws use it, since the draws come from this distribution.
gp_kurt_sn <- (gp_moments_sn[4] -
  4 * gp_mean_sn * gp_moments_sn[3] +
  6 * gp_mean_sn^2 * gp_moments_sn[2] -
  3 * gp_mean_sn^4) /
  gp_sd_sn^4

test_that("the skew-normal for log(lambda) gives near-exact gamma quantiles and moments", {
  # No Monte Carlo here. Observed: quantile errors at most 0.013 exact SDs
  # from 1% to 97.5% and 0.046 at 99% (tolerances 0.03 and 0.1); the Laplace
  # lognormal is off by 0.15 to 1.14 SDs. Mean and SD errors 0.0022 and
  # 0.0046 SDs (tolerances 0.005 and 0.012); Laplace is off by 0.24 and 0.16.
  err_sn <- abs(gp_q_sn - gp_q_exact) / gp_sd_exact
  err_laplace <- abs(gp_q_laplace - gp_q_exact) / gp_sd_exact
  expect_true(all(err_sn < gp_tol_sn))
  expect_true(all(err_sn < err_laplace))

  mean_exact <- gp_S / gp_R
  expect_lt(abs(gp_mean_sn - mean_exact) / gp_sd_exact, 0.005)
  expect_lt(abs(gp_sd_sn - gp_sd_exact) / gp_sd_exact, 0.012)
  lap_var <- gp_fit$cov[1, 1]
  lap_mean <- exp(gp_fit$mode[[1]] + lap_var / 2)
  lap_sd <- lap_mean * sqrt(exp(lap_var) - 1)
  expect_lt(abs(gp_mean_sn - mean_exact), abs(lap_mean - mean_exact))
  expect_lt(abs(gp_sd_sn - gp_sd_exact), abs(lap_sd - gp_sd_exact))
})

test_that("exp() of gamma-Poisson draws matches the exact gamma posterior", {
  set.seed(51)
  lambda <- sampling(gp_fit, nsamp = n_big, fun = function(theta) exp(theta))
  expect_identical(dim(lambda), c(as.integer(n_big), 1L))
  q_emp <- unname(quantile(lambda[, 1], probs))

  # Against the skew-normal target, within Monte Carlo error only. The
  # density of lambda = exp(theta) at q is f_sn(log q) / q. Observed: largest
  # deviation 1.3 SEs.
  f_sn <- sn_density_ref(log(gp_q_sn), gp_sn$xi, gp_sn$omega, gp_sn$alpha) /
    gp_q_sn
  se_sn <- quantile_se(probs, n_big, f_sn)
  expect_lt(max(abs(q_emp - gp_q_sn) / se_sn), z_mc)

  # Against the exact gamma, within the approximation tolerance plus Monte
  # Carlo error, and closer than the Laplace lognormal at every quantile.
  # Observed: largest error 0.027 exact SDs.
  err_draws <- abs(q_emp - gp_q_exact)
  expect_true(all(err_draws < gp_tol_sn * gp_sd_exact + z_mc * se_sn))
  expect_true(all(err_draws < abs(gp_q_laplace - gp_q_exact)))

  # Mean and SD of lambda. The Monte Carlo SE of the mean is SD / sqrt(n),
  # and of the SD, SD sqrt((kurtosis - 1) / (4 n)), both from the
  # skew-normal's moments (over 1000 seeds these z-scores have SD 1.00 and
  # 1.01). Observed: 1.2 and 1.2 SEs from the skew-normal's mean and SD, 0.9
  # and 0.5 SEs from the exact values, and closer to them than the Laplace
  # lognormal's (0.24 and 0.16 SDs off).
  mean_exact <- gp_S / gp_R
  se_mean <- gp_sd_sn / sqrt(n_big)
  se_sd <- gp_sd_sn * sqrt((gp_kurt_sn - 1) / (4 * n_big))
  expect_lt(abs(mean(lambda) - gp_mean_sn), z_mc * se_mean)
  expect_lt(abs(sd(lambda) - gp_sd_sn), z_mc * se_sd)
  expect_lt(
    abs(mean(lambda) - mean_exact),
    0.005 * gp_sd_exact + z_mc * se_mean
  )
  expect_lt(abs(sd(lambda) - gp_sd_exact), 0.012 * gp_sd_exact + z_mc * se_sd)

  lap_var <- gp_fit$cov[1, 1]
  lap_mean <- exp(gp_fit$mode[[1]] + lap_var / 2)
  lap_sd <- lap_mean * sqrt(exp(lap_var) - 1)
  expect_lt(abs(mean(lambda) - mean_exact), abs(lap_mean - mean_exact))
  expect_lt(abs(sd(lambda) - gp_sd_exact), abs(lap_sd - gp_sd_exact))
})
