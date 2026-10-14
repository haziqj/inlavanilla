# The diagnostics of a fit (R/diagnostics.R): diagnostics(), the fit-time
# warning, and the stored scans.
#
# One warning, of class inlavanilla_diagnostics_warning, in INLAvaan's words
# and order: the mode shift above 0.1 posterior SDs, NMAD above 0.1, the
# scan-end mass above 0.05, heavy tails (inlavanilla's own), the VB shift
# above 1 SD and the Hessian condition number above 1e8. A healthy fit stays
# silent. diagnostics() returns INLAvaan's global and per-parameter
# diagnostics, computed by the same function as the warning.

# A Poisson regression with two mildly skewed, correlated parameters: healthy
reg_x <- c(-1.75, -1.25, -0.75, -0.25, 0.25, 0.75, 1.25, 1.75)
reg_y <- c(1, 0, 2, 1, 3, 4, 2, 6)
lp_reg <- function(theta) {
  eta <- theta[["intercept"]] + theta[["slope"]] * reg_x
  sum(dpois(reg_y, exp(eta), log = TRUE)) +
    dnorm(theta[["intercept"]], 0, 2, log = TRUE) +
    dnorm(theta[["slope"]], 0, 2, log = TRUE)
}
init_reg <- c(intercept = 0.2, slope = 0.1)

# A Cauchy posterior: far heavier tails than any skew-normal
lp_cauchy <- function(theta) dt(theta[[1]], 1, log = TRUE)

# The parts of the diagnostics of a fit with two healthy parameters, to
# perturb one check at a time
healthy_optim <- list(
  convergence = 0L,
  message = "relative convergence (4)",
  iterations = 5L,
  grad = c(a = 1e-7, b = -2e-7),
  grad_mode = 2e-7,
  hess_cond = 3
)
healthy_vb <- list(iterations = 3L, fallback = FALSE, kld_global = 0.03)
healthy <- data.frame(
  mode_shift_sigma = c(1e-7, 2e-7),
  kld = 0.5 * c(0.1, -0.2)^2,
  vb_shift = c(0.01, -0.02),
  vb_shift_sigma = c(0.1, -0.2),
  nmad = c(0.01, 0.02),
  alpha = c(0.5, -0.3),
  scan_end_mass = c(1e-3, 2e-3),
  sn_vs_vb = c(0.01, 0.02),
  gamma1 = c(0, 0),
  tails = c(TRUE, TRUE),
  ridge_bend = c(0.01, 0.02),
  row.names = c("a", "b")
)

dg_of <- function(d = healthy, optim = healthy_optim, vb = healthy_vb) {
  diagnostics_internal(optim, vb, d)
}
problems <- function(..., gaussian = FALSE) {
  diagnostic_problems(dg_of(...), gaussian)
}
# The warning's message as plain text (no colours, no wrapping), or "" when
# there is none
warning_text <- function(..., gaussian = FALSE) {
  w <- tryCatch(
    warn_diagnostics(dg_of(...), gaussian),
    inlavanilla_diagnostics_warning = identity
  )
  if (inherits(w, "condition")) w[["text"]] else ""
}

test_that("a healthy fit raises no diagnostics warning", {
  expect_no_warning(
    fit <- suppressMessages(inla_vanilla(lp_reg, init_reg)),
    class = "inlavanilla_diagnostics_warning"
  )
  expect_identical(fit$warnings, character())
  expect_length(problems(), 0)
  expect_false(warn_diagnostics(dg_of(), gaussian = FALSE))
})

test_that("each check adds its own line to the warning, in INLAvaan's words", {
  # The mode, with nlminb()'s message and a hint
  o <- healthy_optim
  o$grad <- c(a = 0.5, b = 0.01)
  d <- healthy
  d$mode_shift_sigma[2] <- 0.15
  text <- warning_text(d, o)
  expect_match(
    text,
    paste0(
      "x The gradient at the posterior mode is not zero (max |grad| = 0.5): ",
      "a Newton step would move `b` by 0.15 posterior SDs.\n",
      "i `nlminb()` says: \"relative convergence (4)\".\n",
      "i Try other starting values, or supply `grad`."
    ),
    fixed = TRUE
  )
  o$message <- "iteration limit reached without convergence (10)"
  expect_match(
    warning_text(d, o),
    "Raise `iter.max` or `eval.max` in `control$nlminb`.",
    fixed = TRUE
  )
  # strictly above the level
  d$mode_shift_sigma[2] <- 0.1
  expect_length(problems(d, o), 0)

  d <- healthy
  d$nmad[1] <- 0.12
  expect_match(
    warning_text(d),
    paste0(
      "x The fitted marginal deviates from the scanned posterior ",
      "(NMAD > 0.1) for `a` (0.12)."
    ),
    fixed = TRUE
  )
  d$nmad[1] <- 0.1
  expect_length(problems(d), 0)

  d <- healthy
  d$scan_end_mass[2] <- 0.06
  expect_match(
    warning_text(d),
    paste0(
      "x The fitted marginal puts more than 0.05 of its mass beyond the ",
      "scanned window (4 posterior SDs either side of the mode) for `b` ",
      "(0.06); its credible limits rely on extrapolation."
    ),
    fixed = TRUE
  )

  d <- healthy
  d$tails[2] <- FALSE
  expect_match(
    warning_text(d),
    paste0(
      "x Heavy tails for `b`: the density is still over 5% of its peak at ",
      "the ends of the scanned window (4 posterior SDs either side of the ",
      "mode), and the skew-normal extrapolates beyond it."
    ),
    fixed = TRUE
  )

  # The VB shift names the worst parameter only, as INLAvaan's does
  d <- healthy
  d$vb_shift_sigma <- c(1.2, -3)
  text <- warning_text(d)
  expect_match(
    text,
    paste0(
      "x The VB correction shifted `b` by 3 posterior SDs; the Gaussian ",
      "approximation at the mode may be inaccurate."
    ),
    fixed = TRUE
  )
  expect_no_match(text, "`a`")

  o <- healthy_optim
  o$hess_cond <- 2e10
  expect_match(
    warning_text(optim = o),
    paste0(
      "^Fit diagnostics flagged 1 potential issue:\n",
      "x The Hessian at the mode is near-singular \\(condition number = ",
      "2e\\+10\\); some parameters may be weakly identified\\.\n"
    )
  )
  o$hess_cond <- 1e8
  expect_length(problems(optim = o), 0)
})

test_that("the VB fallback is recorded, not warned about, as in INLAvaan", {
  vb_fell <- modifyList(healthy_vb, list(fallback = TRUE))
  expect_length(problems(vb = vb_fell), 0)
  expect_identical(dg_of(vb = vb_fell)$global[["vb_fallback"]], 1)
  # A large shift is reported whichever route found it
  d <- healthy
  d$vb_shift_sigma[2] <- -1.5
  bullets <- problems(d, vb = vb_fell)
  expect_length(bullets, 1)
  expect_match(bullets, "The VB correction shifted")
  # No VB correction, no VB check
  d[c("kld", "vb_shift", "vb_shift_sigma", "sn_vs_vb")] <- NA_real_
  expect_length(problems(d, vb = NULL), 0)
})

test_that("lists name at most three parameters, worst first, with values", {
  # Five parameters, one of them named with braces, which must not be read
  # as cli markup
  nm <- c("p1", "b{1}", "p3", "p4", "p5")
  d <- healthy[rep(1, 5), ]
  rownames(d) <- nm
  o <- modifyList(healthy_optim, list(grad = setNames(rep(1e-7, 5), nm)))
  d$nmad <- c(0.2, 0.5, 0.01, 0.3, 0.15)
  d$scan_end_mass <- c(0.06, 0.01, 0.01, 0.01, 0.07)
  d$tails <- c(TRUE, FALSE, TRUE, FALSE, TRUE)
  text <- warning_text(d, o)
  expect_match(
    text,
    "for `b{1}` (0.50), `p4` (0.30), `p1` (0.20) and 1 other.",
    fixed = TRUE
  )
  expect_match(text, "for `p5` (0.07), `p1` (0.06); its", fixed = TRUE)
  expect_match(text, "Heavy tails for `b{1}`, `p4`: the", fixed = TRUE)
  d$nmad <- rep(0.2, 5)
  expect_match(
    warning_text(d, o),
    "for `p1` (0.20), `b{1}` (0.20), `p3` (0.20) and 2 others.",
    fixed = TRUE
  )
})

test_that("every problem at once gives one warning, in INLAvaan's order", {
  o <- healthy_optim
  o$grad <- c(a = 2, b = 0)
  o$hess_cond <- 1e9
  d <- healthy
  d$mode_shift_sigma[1] <- 0.5
  d$tails[1] <- FALSE
  d$nmad[1] <- 0.2
  d$scan_end_mass[1] <- 0.3
  d$vb_shift_sigma[1] <- 2
  bullets <- problems(d, o)
  expect_identical(names(bullets), c("x", "i", "i", "x", "x", "x", "x", "x"))
  text <- warning_text(d, o)
  lines <- strsplit(text, "\n")[[1]]
  expect_identical(lines[1], "Fit diagnostics flagged 6 potential issues:")
  # each line shows its own value
  for (pattern in c(
    "^x The gradient at the posterior mode .* by 0.5 posterior SDs",
    "^i `nlminb\\(\\)` says",
    "^i Try other starting values",
    "^x The fitted marginal deviates .* `a` \\(0.20\\)",
    "^x The fitted marginal puts .* `a` \\(0.30\\)",
    "^x Heavy tails for `a`",
    "^x The VB correction shifted `a` by 2 posterior SDs",
    "^x The Hessian at the mode is near-singular \\(condition number = 1e\\+09\\)",
    "^i Inspect with `diagnostics\\(fit\\)` and `diagnostics\\(fit, type = \"param\"\\)`; `vignette\\(\"diagnostics\", package = \"inlavanilla\"\\)` explains each check\\.$"
  )) {
    expect_length(grep(pattern, lines), 1)
  }
  expect_identical(
    grep("^[xi] ", lines),
    2:10
  )
})

test_that("a value that is not finite is flagged; a check that did not run is not", {
  d <- healthy
  d$nmad[2] <- NaN
  expect_match(warning_text(d), "(NMAD > 0.1) for `b` (NaN)", fixed = TRUE)
  d <- healthy
  d$mode_shift_sigma[1] <- NA
  expect_match(warning_text(d), "would move `a` by NA posterior SDs")

  # Gaussian marginals: no scans, so only the mode, Hessian and VB checks
  d <- healthy
  d[c("nmad", "alpha", "scan_end_mass", "gamma1", "tails", "ridge_bend")] <- NA
  d$sn_vs_vb <- NA_real_
  expect_length(problems(d, gaussian = TRUE), 0)
  d$vb_shift_sigma[1] <- 3
  expect_length(problems(d, gaussian = TRUE), 1)
})

test_that("a heavy-tailed posterior warns once at the end of the fit", {
  w <- expect_warning(
    fit <- suppressMessages(inla_vanilla(lp_cauchy, c(theta = 0.1))),
    class = "inlavanilla_diagnostics_warning"
  )
  expect_s3_class(w, "inlavanilla_warning")
  text <- cli::ansi_strip(conditionMessage(w))
  expect_match(text, "Heavy tails for `theta`")
  expect_match(text, "deviates from the scanned posterior \\(NMAD > 0.1\\)")
  expect_length(fit$warnings, 1)
  expect_match(fit$warnings, "^Fit diagnostics flagged 2 potential issues:")

  # silent = TRUE: recorded, not signalled
  expect_silent(fit <- inla_vanilla(lp_cauchy, c(theta = 0.1), silent = TRUE))
  expect_identical(
    gsub("\\s+", " ", fit$warnings),
    gsub("\\s+", " ", text)
  )
})

test_that("an ill-conditioned Hessian warns, whatever the units", {
  # y ~ N(a + b, 1) with N(0, tau^2) priors: only a + b is identified, and
  # the condition number of the Hessian, with each parameter in posterior-SD
  # units, is 2 n tau^2 + 1. Observed: 2.000003e10 against 2e10 + 1.
  y <- qnorm(ppoints(100))
  lp <- function(theta, tau) {
    sum(dnorm(y, theta[[1]] + theta[[2]], 1, log = TRUE)) +
      sum(dnorm(theta, 0, tau, log = TRUE))
  }
  gr <- function(theta, tau) sum(y - theta[[1]] - theta[[2]]) - theta / tau^2
  w <- expect_warning(
    fit <- suppressMessages(
      inla_vanilla(lp, c(a = 0, b = 0), grad = gr, tau = 1e4)
    ),
    class = "inlavanilla_diagnostics_warning"
  )
  expect_match(conditionMessage(w), "The Hessian at the mode is near-singular")
  expect_equal(
    fit$diagnostics$optim$hess_cond,
    2 * 100 * 1e4^2 + 1,
    tolerance = 1e-4
  )
  expect_identical(
    diagnostics(fit)[["hess_cond"]],
    fit$diagnostics$optim$hess_cond
  )
  fit <- inla_vanilla(lp, c(a = 0, b = 0), grad = gr, tau = 1, silent = TRUE)
  expect_equal(fit$diagnostics$optim$hess_cond, 201, tolerance = 1e-4)

  # Independent parameters with posterior SDs 1e-5 and 1e3: INLAvaan's raw
  # condition number would be 1e16, but in posterior-SD units it is 1
  lp2 <- function(theta) {
    dnorm(theta[[1]], 0, 1e-5, log = TRUE) +
      dnorm(theta[[2]], 0, 1e3, log = TRUE)
  }
  expect_no_warning(
    fit <- inla_vanilla(lp2, c(a = 1e-6, b = 10), silent = TRUE)
  )
  expect_identical(fit$warnings, character())
  expect_equal(fit$diagnostics$optim$hess_cond, 1, tolerance = 1e-6)
})

test_that("the mode shift is INLAvaan's Newton step in posterior SDs", {
  # A correlated Gaussian, rho = 0.999, with the "mode" off by 1e-3 along
  # the narrow direction: INLAvaan's |Sigma g| / sd is small, while the
  # whitened gradient, which decides the restart, is not
  rho <- 0.999
  S <- matrix(c(1, rho, rho, 1), 2)
  P <- solve(S)
  lp <- function(x) -0.5 * sum(x * (P %*% x))
  gr <- function(x) -as.numeric(P %*% x)
  off <- c(1, -1) * 1e-3 / sqrt(2)
  lap <- laplace_at(lp, gr, off, lp(off), c("a", "b"))
  g <- gr(off)
  expect_equal(
    unname(lap$mode_shift),
    abs(as.numeric(lap$Sigma %*% g)) / unname(sqrt(diag(lap$Sigma))),
    tolerance = 1e-8
  )
  expect_equal(unname(lap$mode_shift), abs(off), tolerance = 1e-4)
  expect_gt(lap$grad_mode, 10 * max(lap$mode_shift))

  # With one parameter the two agree
  lp1 <- function(x) dnorm(x[[1]], 1, 2, log = TRUE)
  lap <- laplace_at(lp1, NULL, 1.3, lp1(1.3), "a")
  expect_equal(unname(lap$mode_shift), lap$grad_mode, tolerance = 1e-8)

  # At an accurate mode it is close to 0 (observed 7e-7, from finite
  # differences)
  fit <- inla_vanilla(lp_reg, init_reg, silent = TRUE)
  expect_lt(max(fit$diagnostics$marginals$mode_shift_sigma), 1e-4)
})

test_that("diagnostics() returns INLAvaan's diagnostics, in its order", {
  fit <- inla_vanilla(lp_reg, init_reg, silent = TRUE)
  g <- diagnostics(fit)
  expect_s3_class(g, "inla_vanilla_diagnostics")
  expect_identical(
    names(g),
    c(
      "npar",
      "converged",
      "iterations",
      "grad_inf",
      "grad_inf_rel",
      "grad_l2",
      "mode_shift_max",
      "hess_cond",
      "vb_applied",
      "vb_shift_max",
      "vb_kld_global",
      "kld_max",
      "kld_mean",
      "nmad_max",
      "nmad_mean",
      "scan_end_mass_max",
      "vb_fallback",
      "sn_vs_vb_max",
      "ridge_bend_max",
      "n_heavy_tails"
    )
  )
  expect_identical(diagnostics(fit, type = "global"), g)
  p <- diagnostics(fit, type = "param")
  expect_s3_class(p, "inla_vanilla_diagnostics_param")
  expect_identical(
    names(p),
    c(
      "names",
      "grad",
      "grad_abs",
      "grad_rel",
      "mode_shift_sigma",
      "kld",
      "vb_shift",
      "vb_shift_sigma",
      "nmad",
      "alpha",
      "scan_end_mass",
      "sn_vs_vb",
      "gamma1",
      "tails",
      "ridge_bend"
    )
  )
  expect_identical(p$names, names(init_reg))
  expect_identical(diagnostics(fit, type = "par"), p)

  # The definitions: INLAvaan's, from what the fit stores
  d <- fit$diagnostics$marginals
  grad <- fit$diagnostics$optim$grad
  expect_identical(g[["npar"]], 2)
  expect_identical(
    g[["converged"]],
    as.numeric(fit$diagnostics$optim$convergence == 0)
  )
  expect_identical(g[["grad_inf"]], max(abs(grad)))
  expect_identical(
    g[["grad_inf_rel"]],
    max(abs(grad) / (abs(fit$mean_vb) + 1e-6))
  )
  expect_identical(g[["grad_l2"]], sqrt(sum(grad^2)))
  expect_identical(g[["mode_shift_max"]], max(d$mode_shift_sigma))
  expect_identical(g[["vb_kld_global"]], fit$diagnostics$vb$kld_global)
  expect_equal(g[["vb_kld_global"]], fit$logml - fit$logml_laplace)
  expect_identical(g[["kld_mean"]], mean(d$kld))
  expect_identical(g[["nmad_max"]], max(d$nmad))
  expect_identical(g[["scan_end_mass_max"]], max(d$scan_end_mass))
  expect_identical(g[["ridge_bend_max"]], max(d$ridge_bend))
  expect_identical(g[["n_heavy_tails"]], 0)
  expect_equal(p$grad, unname(grad))
  for (col in setdiff(names(p), c("names", "grad", "grad_abs", "grad_rel"))) {
    expect_identical(p[[col]], d[[col]])
  }
  # The gradient is that of minus logpost at the mode: finite differences
  # here, the user's gradient when there is one
  gr_reg <- function(theta) {
    eta <- theta[["intercept"]] + theta[["slope"]] * reg_x
    res <- reg_y - exp(eta)
    c(sum(res), sum(reg_x * res)) - theta / 4
  }
  fit_g <- inla_vanilla(lp_reg, init_reg, grad = gr_reg, silent = TRUE)
  expect_equal(
    unname(fit_g$diagnostics$optim$grad),
    -unname(gr_reg(fit_g$mode)),
    tolerance = 1e-12
  )

  expect_error(
    diagnostics(fit, type = "everything"),
    "`type` must be one of \"global\" or \"param\"",
    class = "inlavanilla_error"
  )
})

test_that("diagnostics() prints in INLAvaan's formats", {
  testthat::local_reproducible_output(width = 80)
  fit <- inla_vanilla(lp_reg, init_reg, silent = TRUE)
  g <- diagnostics(fit)
  out <- capture_output(res <- withVisible(print(g)))
  expect_false(res$visible)
  expect_identical(res$value, g)
  # Counts as integers, the gradient norms and wide-ranging quantities in
  # scientific notation, the rest to four decimals
  fmt <- function(nm) {
    lines <- strsplit(out, "\n")[[1]]
    i <- grep(paste0("(^| )", nm, "( |$)"), lines)
    heads <- strsplit(trimws(lines[i]), " +")[[1]]
    vals <- strsplit(trimws(lines[i + 1]), " +")[[1]]
    vals[match(nm, heads)]
  }
  expect_identical(fmt("npar"), "2")
  expect_identical(
    fmt("iterations"),
    as.character(fit$diagnostics$optim$iterations)
  )
  expect_match(fmt("grad_inf"), "^[0-9]\\.[0-9]{2}e[-+][0-9]{2}$")
  expect_match(fmt("hess_cond"), "^[0-9]\\.[0-9]{2}e[-+][0-9]{2}$")
  expect_match(fmt("scan_end_mass_max"), "^[0-9]\\.[0-9]{2}e[-+][0-9]{2}$")
  expect_match(fmt("nmad_max"), "^0\\.[0-9]{4}$")
  expect_identical(fmt("n_heavy_tails"), "0")

  # The table rounds to four decimals, but keeps three significant figures
  # of the scan-end mass, which is below the fourth decimal when healthy
  p <- diagnostics(fit, type = "param")
  out <- capture_output(res <- withVisible(print(p)))
  expect_identical(res$value, p)
  expect_match(out, "scan_end_mass")
  expect_match(out, format(signif(p$scan_end_mass[1], 3)), fixed = TRUE)

  # Missing values print as NA
  fit <- inla_vanilla(lp_reg, init_reg, vb = FALSE, silent = TRUE)
  out <- capture_output(print(diagnostics(fit)))
  expect_match(out, "NA")
})

test_that("diagnostics() of other objects is an error, or INLAvaan's", {
  expect_error(
    diagnostics(1:3),
    "There is no `diagnostics\\(\\)` method for an integer vector",
    class = "inlavanilla_error"
  )
  # An S4 object is passed on to INLAvaan's diagnostics() only when INLAvaan
  # is loaded
  skip_if(isNamespaceLoaded("INLAvaan"))
  expect_error(
    diagnostics(asS4(structure(list(), class = "not_inlavaan"))),
    "There is no `diagnostics\\(\\)` method",
    class = "inlavanilla_error"
  )
})

test_that("the stored scans reproduce the NMAD of each skew-normal fit", {
  # NMAD is the largest gap between the corrected scan, scaled to a peak of
  # 1, and the fitted skew-normal on the scan coordinate (INLAvaan's
  # sn_fit plot draws both); recomputed here from fit$scan
  fit <- inla_vanilla(lp_reg, init_reg, silent = TRUE)
  expect_identical(fit$scan$z, seq(-4, 4, length.out = 21))
  expect_identical(dim(fit$scan$logpost), c(21L, 2L))
  expect_identical(colnames(fit$scan$logpost), names(init_reg))
  expect_identical(names(fit$scan$logC), names(init_reg))
  # The scan passes through the mode
  expect_identical(
    unname(fit$scan$logpost[11, 1]),
    lp_reg(fit$mode)
  )
  nmad_from_scan <- function(fit) {
    d <- fit$diagnostics$marginals
    sds <- sqrt(diag(fit$cov))
    z <- fit$scan$z
    vapply(
      seq_len(nrow(d)),
      function(j) {
        y <- fit$scan$logpost[, j] + d$gamma1[j] * z
        y <- y - max(y)
        fitted <- dsnorm(
          z,
          (fit$sn[j, "xi"] - fit$mode[[j]]) / sds[[j]],
          fit$sn[j, "omega"] / sds[[j]],
          fit$sn[j, "alpha"],
          logC = fit$scan$logC[[j]]
        )
        max(abs(exp(y) - fitted))
      },
      numeric(1)
    )
  }
  expect_equal(
    nmad_from_scan(fit),
    fit$diagnostics$marginals$nmad,
    tolerance = 1e-10
  )
  # The same with another grid, and with a gradient
  gr_reg <- function(theta) {
    eta <- theta[["intercept"]] + theta[["slope"]] * reg_x
    res <- reg_y - exp(eta)
    c(sum(res), sum(reg_x * res)) - theta / 4
  }
  fit <- inla_vanilla(
    lp_reg,
    init_reg,
    grad = gr_reg,
    silent = TRUE,
    control = list(n_grid = 9)
  )
  expect_identical(fit$scan$z, seq(-4, 4, by = 1))
  expect_equal(
    nmad_from_scan(fit),
    fit$diagnostics$marginals$nmad,
    tolerance = 1e-10
  )
})
