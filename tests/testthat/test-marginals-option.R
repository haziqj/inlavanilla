# The `marginals` argument and its independence from `vb`.
#
# marginals = "gaussian" skips the skew-normal scans and uses the Gaussian
# marginals N(mean_vb_j, Sigma_jj), centred at the VB mean when vb = TRUE and
# at the mode when vb = FALSE. Every other stage must be unchanged: the mode,
# the covariance, the VB correction and the log marginal likelihood are the
# same numbers whichever marginals are asked for.

# A Poisson regression with fixed data: two correlated, skewed parameters, so
# the Gaussian and skew-normal marginals differ visibly.
reg_x <- c(-1.75, -1.25, -0.75, -0.25, 0.25, 0.75, 1.25, 1.75)
reg_y <- c(1, 0, 2, 1, 3, 4, 2, 6)
lp_reg <- function(theta) {
  eta <- theta[["intercept"]] + theta[["slope"]] * reg_x
  sum(dpois(reg_y, exp(eta), log = TRUE)) +
    dnorm(theta[["intercept"]], 0, 2, log = TRUE) +
    dnorm(theta[["slope"]], 0, 2, log = TRUE)
}
init_reg <- c(intercept = 0.2, slope = 0.1)

fit_sn_vb <- inla_vanilla(lp_reg, init_reg, silent = TRUE)
fit_sn <- inla_vanilla(lp_reg, init_reg, vb = FALSE, silent = TRUE)
fit_g_vb <- inla_vanilla(
  lp_reg,
  init_reg,
  marginals = "gaussian",
  silent = TRUE
)
fit_g <- inla_vanilla(
  lp_reg,
  init_reg,
  vb = FALSE,
  marginals = "gaussian",
  silent = TRUE
)

scrub_output <- function(x) {
  x <- gsub("[0-9.]+(ms|s)\\b", "<time>", x)
  x <- gsub("-+ <time> --", "-- <time> --", x)
  x <- gsub("[0-9]+ (logpost|grad)\\b", "<n> \\1", x)
  x <- gsub("[0-9]+ iterations?\\b", "<n> iteration(s)", x)
  x <- gsub("gradient [-+0-9.e]+ SDs", "gradient <g> SDs", x)
  x
}

test_that("marginals takes one of two values, or an abbreviation", {
  expect_identical(fit_sn_vb$marginals, "skewnormal")
  expect_identical(fit_g_vb$marginals, "gaussian")
  # Unambiguous abbreviations, in any case, as with match.arg()
  choices <- c("skewnormal", "gaussian")
  for (ok in c("skew", "s", "SKEW", "Skewnormal")) {
    expect_identical(check_choice(ok, choices, "marginals"), "skewnormal")
  }
  for (ok in c("gaus", "gauss", "g", "Gaussian")) {
    expect_identical(check_choice(ok, choices, "marginals"), "gaussian")
  }
  fit <- inla_vanilla(lp_reg, init_reg, marginals = "gauss", silent = TRUE)
  expect_identical(fit$marginals, "gaussian")
  expect_identical(fit$summary, fit_g_vb$summary)
  for (bad in list(
    "laplace",
    "normal",
    "",
    "skewed",
    NA_character_,
    1,
    TRUE,
    NULL
  )) {
    expect_error(
      inla_vanilla(lp_reg, init_reg, marginals = bad, silent = TRUE),
      "`marginals` must be one of \"skewnormal\" or \"gaussian\"",
      class = "inlavanilla_error"
    )
  }
  # The full vector in another order is not the default
  expect_error(
    inla_vanilla(
      lp_reg,
      init_reg,
      marginals = c("gaussian", "skewnormal"),
      silent = TRUE
    ),
    class = "inlavanilla_error"
  )
  expect_snapshot(
    error = TRUE,
    inla_vanilla(lp_reg, init_reg, marginals = "laplace", silent = TRUE)
  )
})

test_that("a data argument named like marginals reaches logpost", {
  # `marginals` comes after `...`, so `marg` cannot partially match it
  lp <- function(theta, marg) dnorm(theta[[1]], marg, 1, log = TRUE)
  fit <- inla_vanilla(lp, c(a = 0), marg = 3, silent = TRUE)
  expect_equal(fit$mode[["a"]], 3, tolerance = 1e-6)
  expect_identical(fit$marginals, "skewnormal")
})

test_that("the choice of marginals leaves the other stages unchanged", {
  # Same mode, covariance, VB mean and evidence: only the marginals differ
  for (pair in list(list(fit_g_vb, fit_sn_vb), list(fit_g, fit_sn))) {
    a <- pair[[1]]
    b <- pair[[2]]
    expect_identical(a$mode, b$mode)
    expect_identical(a$cov, b$cov)
    expect_identical(a$mean_vb, b$mean_vb)
    expect_identical(a$logml, b$logml)
    expect_identical(a$logml_laplace, b$logml_laplace)
    expect_identical(
      a$diagnostics$marginals$vb_shift,
      b$diagnostics$marginals$vb_shift
    )
  }
  # vb = FALSE: the mean is the mode and the evidence is Laplace's
  expect_identical(fit_g$mean_vb, fit_g$mode)
  expect_identical(fit_g$logml, fit_g$logml_laplace)
  expect_false(isTRUE(all.equal(fit_g_vb$logml, fit_g_vb$logml_laplace)))
})

test_that("Gaussian marginals are N(mean_vb, Sigma_jj)", {
  probs <- c(0.01, 0.025, 0.25, 0.5, 0.75, 0.975, 0.99)
  for (fit in list(fit_g_vb, fit_g)) {
    centre <- fit$mean_vb
    sds <- sqrt(diag(fit$cov))
    expect_identical(unname(fit$sn[, "alpha"]), c(0, 0))
    tab <- summary(fit, probs = probs)$table
    expect_equal(tab[, "Mean"], centre, tolerance = 1e-12)
    expect_equal(tab[, "SD"], sds, tolerance = 1e-12)
    expect_equal(tab[, "Mode"], centre, tolerance = 1e-8)
    for (j in 1:2) {
      expect_equal(
        unname(tab[j, 3:9]),
        qnorm(probs, centre[[j]], sds[[j]]),
        tolerance = 1e-10
      )
    }
  }
  # The VB-centred Gaussian sits where the VB correction put it, which for
  # this skewed posterior is away from the mode
  expect_gt(
    max(abs(fit_g_vb$mean_vb - fit_g_vb$mode) / sqrt(diag(fit_g_vb$cov))),
    0.05
  )
})

test_that("an exactly Gaussian posterior is recovered exactly", {
  # theta ~ N(mu, S) with correlation 0.6: the Laplace approximation is exact,
  # the VB shift is zero, and the log-density is normalised, so log p(y) = 0
  mu <- c(1, -2)
  S <- matrix(c(4, 0.6 * 2 * 0.5, 0.6 * 2 * 0.5, 0.25), 2)
  Q <- solve(S)
  lp <- function(theta) {
    d <- theta - mu
    -log(2 * pi) - 0.5 * log(det(S)) - 0.5 * sum(d * (Q %*% d))
  }
  for (vb in c(TRUE, FALSE)) {
    fit <- inla_vanilla(
      lp,
      c(a = 0, b = 0),
      vb = vb,
      marginals = "gaussian",
      silent = TRUE
    )
    expect_equal(unname(fit$summary[, "Mean"]), mu, tolerance = 1e-6)
    expect_equal(unname(fit$cov), S, tolerance = 1e-6)
    expect_equal(fit$logml, 0, tolerance = 1e-6)
  }
})

test_that("Gaussian marginals skip the scans", {
  # With skew-normal marginals each extra grid point costs one call per
  # parameter; with Gaussian marginals the grid is never used.
  count <- function(...) {
    n <- 0L
    lp <- function(theta) {
      n <<- n + 1L
      lp_reg(theta)
    }
    inla_vanilla(lp, init_reg, silent = TRUE, ...)
    n
  }
  for (vb in c(TRUE, FALSE)) {
    n_g21 <- count(vb = vb, marginals = "gaussian")
    n_g41 <- count(vb = vb, marginals = "gaussian", control = list(n_grid = 41))
    n_sn21 <- count(vb = vb)
    n_sn41 <- count(vb = vb, control = list(n_grid = 41))
    expect_identical(n_g21, n_g41)
    expect_identical(n_sn41 - n_sn21, 40L)
    expect_lt(n_g21, n_sn21)
  }
  expect_identical(unname(fit_g$counts[["grad"]]), 0L)
})

test_that("sampling() draws from the joint Gaussian", {
  n <- 50
  for (fit in list(fit_g_vb, fit_g)) {
    set.seed(11)
    draws <- sampling(fit, nsamp = n)
    set.seed(11)
    z <- matrix(rnorm(n * 2), n, 2) %*% chol(cov2cor(fit$cov))
    expected <- sweep(
      sweep(z, 2, sqrt(diag(fit$cov)), "*"),
      2,
      fit$mean_vb,
      "+"
    )
    expect_equal(unname(draws), unname(expected), tolerance = 1e-12)
    expect_identical(colnames(draws), names(init_reg))
  }
  set.seed(2)
  draws <- sampling(fit_g_vb, nsamp = 20000)
  sds <- sqrt(diag(fit_g_vb$cov))
  # Monte Carlo error of the mean: 1/sqrt(20000) = 0.007 SDs
  expect_lt(max(abs(colMeans(draws) - fit_g_vb$mean_vb) / sds), 0.03)
  expect_equal(cor(draws)[1, 2], cov2cor(fit_g_vb$cov)[1, 2], tolerance = 0.02)
})

test_that("the stage lines have no skew-normal line", {
  msgs <- capture_messages(
    inla_vanilla(lp_reg, init_reg, marginals = "gaussian")
  )
  expect_length(msgs, 2)
  expect_match(msgs[1], "Found the posterior mode")
  expect_match(msgs[2], "VB mean correction")
  msgs <- capture_messages(
    inla_vanilla(lp_reg, init_reg, vb = FALSE, marginals = "gaussian")
  )
  expect_length(msgs, 1)
  expect_match(msgs, "Found the posterior mode")

  # verbose = TRUE says that no scans ran
  msgs <- capture_messages(
    inla_vanilla(lp_reg, init_reg, marginals = "gaussian", verbose = TRUE)
  )
  expect_match(msgs, "Gaussian marginals .*: no scans", all = FALSE)
  expect_no_match(msgs, "skew-normal|Shape")
})

test_that("print() and summary() describe Gaussian marginals", {
  testthat::local_reproducible_output(width = 80, crayon = FALSE)
  out <- capture_output(print(fit_g_vb))
  expect_match(out, "Gaussian marginals")
  expect_match(out, "VB-corrected")
  expect_no_match(out, "skew-normal")
  out <- capture_output(print(fit_g))
  expect_match(out, "Gaussian marginals")
  expect_no_match(out, "VB-corrected")

  expect_snapshot(print(fit_g_vb), transform = scrub_output)
  expect_snapshot(print(summary(fit_g_vb)), transform = scrub_output)
  expect_snapshot(print(summary(fit_g)), transform = scrub_output)

  # No skew-normal diagnostics or warnings: the bend that the skew-normal fit
  # of this model flags is not computed
  out <- capture_output(print(summary(fit_g_vb)))
  expect_no_match(out, "Shape|NMAD|Bend|SN - VB|straight scan")
  expect_match(out, "VB shift")
  out <- capture_output(print(summary(fit_g)))
  expect_no_match(out, "VB shift")
  expect_true(all(is.na(summary(fit_g)$diagnostics$nmad)))
})

test_that("plot() labels the curves by the marginals", {
  record_strings <- function(expr) {
    grDevices::pdf(NULL)
    dev <- grDevices::dev.cur()
    on.exit(grDevices::dev.off(dev), add = TRUE)
    grDevices::dev.control("enable")
    force(expr)
    dl <- grDevices::recordPlot()[[1]]
    unique(unlist(lapply(dl, function(el) {
      Filter(is.character, unlist(as.list(el[[2]])[-1], recursive = TRUE))
    })))
  }
  strs <- record_strings(plot(fit_sn_vb))
  expect_true(all(c("Skew-normal", "Laplace") %in% strs))
  strs <- record_strings(plot(fit_g_vb))
  expect_true(all(c("VB Gaussian", "Laplace") %in% strs))
  expect_false("Skew-normal" %in% strs)
  # Laplace marginals have nothing to compare against: no legend
  strs <- record_strings(plot(fit_g))
  expect_false(any(c("Skew-normal", "VB Gaussian", "Laplace") %in% strs))
  expect_true(all(c("intercept", "slope") %in% strs))
})
