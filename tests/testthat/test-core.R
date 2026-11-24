# vanilla_core(), the stages of a fit that inla_vanilla() and INLAvaan both
# run: the plain list it returns, what it prints and raises, and the hooks that
# only INLAvaan passes (hess_mode and gaussian_marginals).

## ----- A logistic regression with every derivative --------------------------

set.seed(11)
core_n <- 60
core_X <- cbind(1, rnorm(core_n), rnorm(core_n))
core_y <- rbinom(core_n, 1, plogis(core_X %*% c(-0.5, 1, -0.7)))
lp_core <- function(theta) {
  eta <- drop(core_X %*% theta)
  sum(core_y * eta - log1p(exp(eta))) + sum(dnorm(theta, 0, 2, log = TRUE))
}
gr_core <- function(theta) {
  p <- plogis(drop(core_X %*% theta))
  drop(crossprod(core_X, core_y - p)) - theta / 4
}
he_core <- function(theta) {
  p <- plogis(drop(core_X %*% theta))
  -crossprod(core_X, core_X * (p * (1 - p))) - diag(3) / 4
}
init_core <- c(b0 = 0, b1 = 0, b2 = 0)

core_fit <- function(...) {
  vanilla_core(lp_core, gr_core, init_core, ..., silent = TRUE)
}

## ----- The returned list ----------------------------------------------------

test_that("vanilla_core() returns a plain list, the fit of inla_vanilla()", {
  core <- core_fit()
  expect_null(attr(core, "class"))
  expect_true(all(
    c(
      "mode",
      "mean_vb",
      "cov",
      "sn",
      "logml",
      "logml_laplace",
      "lp_max",
      "opt",
      "restart",
      "diagnostics",
      "scan",
      "timing",
      "counts"
    ) %in%
      names(core)
  ))
  m <- length(init_core)
  expect_identical(names(core$mode), names(init_core))
  expect_identical(names(core$mean_vb), names(init_core))
  expect_identical(dim(core$cov), c(m, m))
  expect_identical(
    colnames(core$sn),
    c(
      "xi",
      "omega",
      "alpha",
      "gamma1",
      "nmad",
      "tails",
      "ridge_bend",
      "logC",
      "unscanned_mass"
    )
  )
  expect_identical(rownames(core$sn), names(init_core))
  expect_identical(names(core$diagnostics), c("optim", "vb", "marginals"))
  expect_identical(
    names(core$timing),
    c("init", "laplace", "vb", "marginals", "diagnostics", "total")
  )
  expect_identical(names(core$counts), c("logpost", "grad", "hess"))
  expect_equal(core$lp_max, lp_core(core$mode))
  expect_identical(core$opt$par, unname(core$mode))

  fit <- inla_vanilla(lp_core, init_core, gr_core, silent = TRUE)
  expect_identical(fit$sn, core$sn[, c("xi", "omega", "alpha")])
  expect_identical(fit$mode, core$mode)
  expect_identical(fit$cov, core$cov)
  expect_identical(fit$mean_vb, core$mean_vb)
  expect_identical(fit$logml, core$logml)
  expect_identical(fit$diagnostics, core$diagnostics)
  expect_identical(fit$counts, core$counts)
})

test_that("vanilla_core() prints nothing when silent, and lets warnings through", {
  expect_silent(core_fit())
  # A wrong gradient: the core warns about it, does not record or muffle the
  # warning, and leaves the diagnostics warning to its caller
  gr_wrong <- function(theta) c(1, 2, 1) * gr_core(theta)
  w <- character()
  core <- withCallingHandlers(
    vanilla_core(lp_core, gr_wrong, init_core, silent = TRUE),
    warning = function(cond) {
      w <<- c(w, class(cond)[1])
      invokeRestart("muffleWarning")
    }
  )
  expect_true("inlavanilla_warning" %in% w)
  expect_false("inlavanilla_diagnostics_warning" %in% w)
  expect_gt(max(core$diagnostics$marginals$grad_gap), 0.1)
  # The problems it would raise are left to diagnostic_problems()
  dg <- diagnostics_internal(
    core$diagnostics$optim,
    core$diagnostics$vb,
    core$diagnostics$marginals
  )
  expect_true(check_flags(dg, gaussian = FALSE)[["grad_gap"]])
})

test_that("vanilla_core() checks its own arguments", {
  expect_error(
    vanilla_core(lp_core, gr_core, c(0, 0, 0), silent = TRUE),
    "`init` must be a named numeric vector",
    class = "inlavanilla_error"
  )
  expect_error(
    vanilla_core(lp_core, NULL, init_core, hess_mode = he_core),
    "`hess_mode` needs `gr`",
    class = "inlavanilla_error"
  )
  expect_error(
    core_fit(vb_correction = "off"),
    "`vb_correction` must be one of",
    class = "inlavanilla_error"
  )
})

## ----- hess_mode -------------------------------------------------------------

test_that("hess_mode gives the Hessian at the mode, and nothing else", {
  ref <- core_fit()
  core <- core_fit(hess_mode = he_core)
  # nlminb() never sees it, so the mode search is the same
  expect_identical(core$mode, ref$mode)
  expect_identical(core$opt$iterations, ref$opt$iterations)
  # The covariance is its inverse at the mode, which the differences of the
  # gradient approach
  expect_equal(
    unname(core$cov),
    solve(-he_core(core$mode)),
    tolerance = 1e-10
  )
  expect_equal(core$cov, ref$cov, tolerance = 1e-6)
  expect_lt(max(abs(core$sn[, 1:3] - ref$sn[, 1:3])), 1e-5)
  # One call, at the mode: not in the VB correction or the scans, whatever
  # the scan
  expect_identical(core$counts[["hess"]], 1L)
  exact <- core_fit(hess_mode = he_core, marginal_scan = "exact")
  expect_identical(exact$counts[["hess"]], 1L)
  # It saves the 2m gradient calls of the differences at the mode, less the 4
  # of its check against the gradient, and the gradient and its diagonal are
  # still compared with logpost there (a correct one reads about 0.003)
  expect_identical(
    ref$counts[["grad"]] - core$counts[["grad"]],
    2L * length(init_core) - 4L
  )
  expect_true(all(is.finite(core$diagnostics$marginals$grad_gap)))
  expect_lt(max(core$diagnostics$marginals$grad_gap), 0.02)
})

test_that("a wrong hess_mode is caught at the mode", {
  he_wrong <- function(theta) {
    H <- he_core(theta)
    H[1, 2] <- H[2, 1] <- 0
    H
  }
  expect_warning(
    core_fit(hess_mode = he_wrong),
    "The Hessian at the mode does not look like the Hessian of the",
    class = "inlavanilla_warning"
  )
  he_asym <- function(theta) {
    H <- he_core(theta)
    H[lower.tri(H)] <- 0
    H
  }
  w <- character()
  withCallingHandlers(
    core_fit(hess_mode = he_asym),
    warning = function(cond) {
      w <<- c(w, cli::ansi_strip(conditionMessage(cond)))
      invokeRestart("muffleWarning")
    }
  )
  expect_match(w, "not symmetric at the mode", all = FALSE)
})

test_that("verbose output says where the Hessian at the mode came from", {
  msgs <- cli::ansi_strip(capture_messages(
    vanilla_core(
      lp_core,
      gr_core,
      init_core,
      hess_mode = he_core,
      verbose = TRUE
    )
  ))
  text <- paste(msgs, collapse = "")
  expect_match(text, "Hessian at the mode from hess_mode, checked there")
  expect_match(text, "Hessian from hess_mode: 1 call")
})
