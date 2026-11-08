# The `width` and `tk_fallback` arguments: three levels of the width term of
# each skew-normal scan, gamma_j(t) = -1/2 log|H_{-j}|.
#
# "trace" (the default) is the linear term gamma_j' t, "cmp" the exact width
# term along the same straight scan, and "tk" the Laplace marginal of
# Tierney and Kadane, with the other parameters at their conditional mode. On
# a Gaussian posterior the three are the same. On a skewed one "cmp" stays
# close to "trace" near the mode, and "tk" is the closest to the exact
# marginal. tk_fallback rescans with "tk" only the parameters whose ridge
# bend is flagged.

## ----- A Gaussian posterior --------------------------------------------------

gauss_mean <- c(a = 1, b = -2, c = 0.5)
gauss_cov <- matrix(
  c(1, 0.6, -0.3, 0.6, 2, 0.4, -0.3, 0.4, 0.5),
  3,
  3,
  dimnames = list(names(gauss_mean), names(gauss_mean))
)
gauss_prec <- solve(gauss_cov)
lp_gauss <- function(theta) {
  d <- theta - gauss_mean
  -0.5 * sum(d * (gauss_prec %*% d))
}
gr_gauss <- function(theta) -as.numeric(gauss_prec %*% (theta - gauss_mean))
he_gauss <- function(theta) -gauss_prec
init_gauss <- c(a = 0, b = 0, c = 0)

## ----- A skewed posterior with an exact reference ----------------------------

# A normal mean with an unknown SD at n = 12, on the scale (mu, log sigma).
# The marginal of mu is close to a Student t, which the straight scan
# narrows by about 11%; its ridge bend, 0.30, is flagged.
nm_y <- c(
  3.75,
  5.37,
  3.33,
  8.19,
  5.66,
  3.36,
  5.97,
  6.48,
  6.15,
  4.39,
  8.02,
  5.78
)
nm_n <- length(nm_y)
lp_nm <- function(theta) {
  sigma <- exp(theta[["ls"]])
  sum(dnorm(nm_y, theta[["mu"]], sigma, log = TRUE)) +
    dnorm(theta[["mu"]], 0, 10, log = TRUE) +
    dexp(sigma, 1, log = TRUE) +
    theta[["ls"]]
}
gr_nm <- function(theta) {
  mu <- theta[["mu"]]
  ls <- theta[["ls"]]
  r <- nm_y - mu
  s2 <- exp(2 * ls)
  c(sum(r) / s2 - mu / 100, -nm_n + sum(r^2) / s2 - exp(ls) + 1)
}
init_nm <- c(mu = 0, ls = 0)

# Exact marginal means and SDs, from the log-posterior on a fine grid
nm_exact <- local({
  mu <- seq(0, 11, length.out = 1001)
  ls <- seq(-1.5, 2.5, length.out = 1001)
  ss <- sum((nm_y - mean(nm_y))^2)
  lp <- outer(mu, ls, function(m, l) {
    s <- exp(l)
    -nm_n *
      l -
      (ss + nm_n * (mean(nm_y) - m)^2) / (2 * s^2) +
      dnorm(m, 0, 10, log = TRUE) -
      s +
      l
  })
  p <- exp(lp - max(lp))
  p <- p / sum(p)
  moments <- function(x, w) {
    m <- sum(x * w)
    c(mean = m, sd = sqrt(sum((x - m)^2 * w)))
  }
  rbind(mu = moments(mu, rowSums(p)), ls = moments(ls, colSums(p)))
})

fits_nm <- lapply(
  c(trace = "trace", cmp = "cmp", tk = "tk"),
  function(w) inla_vanilla(lp_nm, init_nm, gr_nm, width = w, silent = TRUE)
)

## ----- Tests -----------------------------------------------------------------

test_that("width takes one of three values, or an abbreviation", {
  expect_identical(fits_nm$trace$width, "trace")
  expect_identical(fits_nm$cmp$width, "cmp")
  expect_identical(fits_nm$tk$width, "tk")
  choices <- c("trace", "cmp", "tk")
  expect_identical(check_choice("TK", choices, "width"), "tk")
  expect_identical(check_choice("c", choices, "width"), "cmp")
  expect_identical(check_choice("tr", choices, "width"), "trace")
  expect_error(
    inla_vanilla(lp_nm, init_nm, width = "laplace", silent = TRUE),
    class = "inlavanilla_error"
  )
  expect_error(
    inla_vanilla(lp_nm, init_nm, tk_fallback = NA, silent = TRUE),
    class = "inlavanilla_error"
  )
})

test_that("the default is the linear trace, and unchanged", {
  fit <- inla_vanilla(lp_nm, init_nm, gr_nm, silent = TRUE)
  expect_identical(fit$summary, fits_nm$trace$summary)
  expect_identical(fit$width, "trace")
  expect_false(fit$tk_fallback)
  expect_identical(fit$fallback, character())
  expect_identical(fit$diagnostics$marginals$width, c("trace", "trace"))
  # The width term kept for the plot is the linear one
  expect_equal(
    unname(fit$scan$width),
    unname(outer(fit$scan$z, fit$diagnostics$marginals$gamma1))
  )
})

test_that("on a Gaussian posterior the three widths agree", {
  routes <- list(
    logpost = list(grad = NULL, hess = NULL),
    grad = list(grad = gr_gauss, hess = NULL),
    hess = list(grad = gr_gauss, hess = he_gauss)
  )
  for (route in names(routes)) {
    fits <- lapply(c("trace", "cmp", "tk"), function(w) {
      inla_vanilla(
        lp_gauss,
        init_gauss,
        routes[[route]]$grad,
        routes[[route]]$hess,
        width = w,
        silent = TRUE
      )
    })
    sds <- sqrt(diag(gauss_cov))
    for (k in 2:3) {
      # The width term is flat, the conditional mode lies on the CMP, and so
      # the marginals are the same
      expect_lt(max(abs(fits[[k]]$scan$width)), 1e-4)
      expect_lt(
        max(abs(fits[[k]]$summary - fits[[1]]$summary) / sds),
        1e-4
      )
    }
    # And they are the exact Gaussian marginals
    expect_lt(max(abs(fits[[3]]$summary[, "Mean"] - gauss_mean) / sds), 1e-3)
    expect_lt(max(abs(fits[[3]]$summary[, "SD"] / sds - 1)), 1e-3)
  }
})

test_that("on a skewed posterior tk is the closest to the exact marginal", {
  err <- vapply(
    fits_nm,
    function(f) f$summary[, "SD"] / nm_exact[, "sd"] - 1,
    numeric(2)
  )
  # The straight scan narrows mu by about 11% (trace) to 14% (cmp); the
  # Laplace marginal is within 1%
  expect_lt(abs(err["mu", "tk"]), 0.01)
  expect_lt(abs(err["mu", "tk"]), abs(err["mu", "trace"]))
  expect_lt(abs(err["mu", "tk"]), abs(err["mu", "cmp"]))
  expect_gt(abs(err["mu", "trace"]), 0.08)
  # log sigma, whose ridge does not bend, is right whichever the width
  expect_lt(max(abs(err["ls", ])), 0.005)
  mean_err <- vapply(
    fits_nm,
    function(f) (f$summary[, "Mean"] - nm_exact[, "mean"]) / nm_exact[, "sd"],
    numeric(2)
  )
  expect_lt(max(abs(mean_err)), 0.01)
})

test_that("near the mode the exact width along the CMP is the linear one", {
  z <- fits_nm$cmp$scan$z
  near <- abs(z) <= 1
  w_trace <- fits_nm$trace$scan$width
  w_cmp <- fits_nm$cmp$scan$width
  expect_lt(max(abs(w_cmp[near, ] - w_trace[near, ])), 0.05)
  # The two slopes at the mode agree to the error of the forward difference
  # behind the trace
  expect_equal(
    fits_nm$cmp$diagnostics$marginals$gamma1,
    fits_nm$trace$diagnostics$marginals$gamma1,
    tolerance = 0.02
  )
  # Away from the mode they part: the volume of log sigma curves along the
  # scan of mu, which the linear term misses
  expect_gt(max(abs(w_cmp[, "mu"] - w_trace[, "mu"])), 0.2)
})

test_that("tk takes the heights at the conditional modes", {
  tk <- fits_nm$tk
  z <- tk$scan$z
  # A profile: at least as high as the straight scan everywhere, and the
  # same at the mode
  cmp <- fits_nm$cmp$scan$logpost
  expect_true(all(tk$scan$logpost >= cmp - 1e-8))
  expect_equal(tk$scan$logpost[z == 0, ], cmp[z == 0, ])
  # Check one point by brute force: the conditional mode of log sigma given
  # mu, and its volume
  k <- which.min(abs(z - 3))
  mu <- tk$mode[["mu"]] + z[k] * sqrt(tk$cov[1, 1])
  opt <- optimize(
    function(ls) lp_nm(c(mu = mu, ls = ls)),
    c(-3, 3),
    maximum = TRUE,
    tol = 1e-10
  )
  expect_equal(
    unname(tk$scan$logpost[k, "mu"]),
    opt$objective,
    tolerance = 1e-8
  )
  # The volume relative to that at the joint mode
  curv <- function(mu, ls) {
    h <- 1e-4
    g <- function(ls) gr_nm(c(mu = mu, ls = ls))[2]
    -(g(ls + h) - g(ls - h)) / (2 * h)
  }
  expect_equal(
    unname(tk$scan$width[k, "mu"]),
    -0.5 *
      log(curv(mu, opt$maximum)) +
      0.5 * log(curv(tk$mode[["mu"]], tk$mode[["ls"]])),
    tolerance = 1e-5
  )
})

test_that("the routes to the width term agree", {
  for (w in c("cmp", "tk")) {
    fit_lp <- inla_vanilla(lp_nm, init_nm, width = w, silent = TRUE)
    sds <- fits_nm[[w]]$summary[, "SD"]
    expect_lt(
      max(abs(fit_lp$summary - fits_nm[[w]]$summary) / sds),
      1e-3
    )
  }
  for (w in c("cmp", "tk")) {
    fit_he <- inla_vanilla(
      lp_gauss,
      init_gauss,
      gr_gauss,
      he_gauss,
      width = w,
      silent = TRUE
    )
    expect_gt(fit_he$counts[["hess"]], 1)
  }
})

test_that("the costlier widths make more calls", {
  n_gr <- vapply(fits_nm, function(f) f$counts[["grad"]], numeric(1))
  expect_lt(n_gr[["trace"]], n_gr[["cmp"]])
  expect_lt(n_gr[["cmp"]], n_gr[["tk"]])
  # "cmp" replaces the 4m + 1 calls of the trace, shared by all the parameters,
  # by 2(m - 1) calls at each scan point but the mode, per parameter
  m <- 2
  n_grid <- length(fits_nm$cmp$scan$z)
  expect_identical(
    n_gr[["cmp"]] - n_gr[["trace"]],
    m * (n_grid - 1) * 2 * (m - 1) - (4 * m + 1)
  )
})

test_that("with one parameter the three widths are the same", {
  lp1 <- function(theta) dgamma(exp(theta), 3, 2, log = TRUE) + theta
  fits <- lapply(c("trace", "cmp", "tk"), function(w) {
    inla_vanilla(lp1, c(x = 0), width = w, silent = TRUE)
  })
  expect_identical(fits[[2]]$summary, fits[[1]]$summary)
  expect_identical(fits[[3]]$summary, fits[[1]]$summary)
})

test_that("tk_fallback rescans the flagged parameters with tk", {
  fit <- inla_vanilla(lp_nm, init_nm, gr_nm, tk_fallback = TRUE, silent = TRUE)
  d <- fit$diagnostics$marginals
  # mu has a flagged ridge bend, log sigma does not
  expect_gt(d["mu", "ridge_bend"], summary_levels$ridge_bend)
  expect_lt(d["ls", "ridge_bend"], summary_levels$ridge_bend)
  expect_identical(fit$fallback, "mu")
  expect_identical(d$width, c("tk", "trace"))
  expect_identical(fit$width, "trace")
  expect_true(fit$tk_fallback)
  # Each marginal is that of its own width
  expect_identical(fit$summary["mu", ], fits_nm$tk$summary["mu", ])
  expect_identical(fit$summary["ls", ], fits_nm$trace$summary["ls", ])
  # The ridge bend of a rescanned parameter is reported but not flagged, and
  # summary() says which parameters were rescanned
  dg <- diagnostics_internal(fit$diagnostics$optim, fit$diagnostics$vb, d)
  expect_false(check_flags(dg, gaussian = FALSE)[["ridge_bend"]])
  out <- cli::ansi_strip(capture.output(print(summary(fit))))
  expect_true(any(grepl("TK fallback", out)))
  expect_false(any(grepl("Ridge bend >", out)))
  # Without the fallback the flag stands
  out <- cli::ansi_strip(capture.output(print(summary(fits_nm$trace))))
  expect_true(any(grepl("Ridge bend >", out)))
  expect_false(any(grepl("TK fallback|Width term", out)))
})

test_that("tk_fallback leaves a fit with no flagged ridge bend unchanged", {
  fit <- inla_vanilla(lp_gauss, init_gauss, tk_fallback = TRUE, silent = TRUE)
  ref <- inla_vanilla(lp_gauss, init_gauss, silent = TRUE)
  expect_identical(fit$fallback, character())
  expect_identical(fit$summary, ref$summary)
  expect_identical(fit$counts, ref$counts)
})

test_that("tk_fallback works without a gradient and with cmp", {
  fit <- inla_vanilla(
    lp_nm,
    init_nm,
    width = "cmp",
    tk_fallback = TRUE,
    silent = TRUE
  )
  expect_identical(fit$fallback, "mu")
  expect_identical(fit$diagnostics$marginals$width, c("tk", "cmp"))
  expect_lt(abs(fit$summary["mu", "SD"] / nm_exact["mu", "sd"] - 1), 0.01)
})

test_that("the width is ignored for Gaussian marginals", {
  fit <- inla_vanilla(
    lp_nm,
    init_nm,
    marginals = "gaussian",
    width = "tk",
    tk_fallback = TRUE,
    silent = TRUE
  )
  ref <- inla_vanilla(lp_nm, init_nm, marginals = "gaussian", silent = TRUE)
  expect_identical(fit$summary, ref$summary)
  expect_identical(fit$width, NA_character_)
  expect_identical(fit$fallback, character())
})

test_that("summary(), verbose output and plots show the width", {
  out <- cli::ansi_strip(capture.output(print(summary(fits_nm$tk))))
  expect_true(any(grepl("Width term: Laplace", out)))
  # The ridge bend of a tk scan is not flagged
  expect_false(any(grepl("Ridge bend >", out)))
  out <- cli::ansi_strip(capture.output(print(summary(fits_nm$cmp))))
  expect_true(any(grepl("Width term: exact along", out)))
  msgs <- cli::ansi_strip(capture_messages(
    inla_vanilla(lp_nm, init_nm, gr_nm, width = "cmp", verbose = TRUE)
  ))
  expect_true(any(grepl("Width term by width", msgs)))
  msgs <- cli::ansi_strip(capture_messages(
    inla_vanilla(lp_nm, init_nm, gr_nm, tk_fallback = TRUE, verbose = TRUE)
  ))
  expect_true(any(grepl("TK fallback", msgs)))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_invisible(plot(fits_nm$tk, type = "sn_fit"))
  expect_invisible(plot(fits_nm$cmp, type = "sn_fit_log"))
})
