# The fit-time diagnostics warning (R/diagnostics.R) and the scan-end mass
# flag in summary().
#
# One warning, of class inlavanilla_diagnostics_warning, when a marginal is
# clearly off: heavy tails, NMAD above 0.1, scan-end mass above 0.05, a VB
# shift above 1 SD or the VB fallback. A healthy fit stays silent.

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

# Diagnostics of two healthy parameters, to perturb one check at a time
healthy <- data.frame(
  vb_shift = c(0.1, -0.2),
  sn_vs_vb = c(0.01, 0.02),
  alpha = c(0.5, -0.3),
  gamma1 = c(0, 0),
  nmad = c(0.01, 0.02),
  tails = c(TRUE, TRUE),
  bend = c(0.01, 0.02),
  scan_end_mass = c(1e-3, 2e-3),
  row.names = c("a", "b")
)
vb_ok <- list(iterations = 3L, fallback = FALSE, kld = 0.01)

warning_text <- function(diag, vbd) {
  w <- tryCatch(
    warn_diagnostics(diag, vbd),
    inlavanilla_diagnostics_warning = identity
  )
  if (inherits(w, "condition")) cli::ansi_strip(conditionMessage(w)) else ""
}

test_that("a healthy fit raises no diagnostics warning", {
  expect_no_warning(
    fit <- suppressMessages(inla_vanilla(lp_reg, init_reg)),
    class = "inlavanilla_diagnostics_warning"
  )
  expect_identical(fit$warnings, character())
  expect_length(diagnostic_problems(healthy, vb_ok), 0)
  expect_false(warn_diagnostics(healthy, vb_ok))
})

test_that("each check adds its own line to the warning", {
  d <- healthy
  d$tails[2] <- FALSE
  expect_length(diagnostic_problems(d, vb_ok), 1)
  expect_match(warning_text(d, vb_ok), "Heavy tails for b:")

  d <- healthy
  d$nmad[1] <- 0.12
  expect_match(
    warning_text(d, vb_ok),
    "Poor skew-normal fit for a: largest NMAD 0.12, above 0.1"
  )
  # strictly above the level
  d$nmad[1] <- 0.1
  expect_length(diagnostic_problems(d, vb_ok), 0)

  d <- healthy
  d$scan_end_mass[2] <- 0.06
  expect_match(
    warning_text(d, vb_ok),
    "The skew-normal for b puts up to 0.06 of its mass beyond the scanned"
  )

  d <- healthy
  d$vb_shift[2] <- -1.5
  expect_match(
    warning_text(d, vb_ok),
    "moved the mean of b by up to 1.5 posterior SDs"
  )
  # The fallback is reported instead of the shift
  vb_fell <- modifyList(vb_ok, list(fallback = TRUE))
  bullets <- diagnostic_problems(d, vb_fell)
  expect_length(bullets, 1)
  expect_match(bullets, "fallback optimiser")
  # No VB correction, no VB check
  expect_length(diagnostic_problems(d, NULL), 0)

  # Everything at once: one warning, four lines, then the pointer to help
  d <- healthy
  d$tails[1] <- FALSE
  d$nmad[1] <- 0.2
  d$scan_end_mass[1] <- 0.3
  d$vb_shift[1] <- 2
  expect_length(diagnostic_problems(d, vb_ok), 4)
  text <- warning_text(d, vb_ok)
  # each line shows its own value
  expect_match(text, "largest NMAD 0.2,")
  expect_match(text, "puts up to 0.3 of its mass")
  expect_match(text, "by up to 2 posterior SDs")
  expect_match(text, "^The approximation may be poor for some parameters")
  expect_match(text, 'vignette\\("diagnostics", package = "inlavanilla"\\)')
})

test_that("Gaussian marginals are checked for the VB correction only", {
  d <- healthy
  d[c(
    "sn_vs_vb",
    "alpha",
    "gamma1",
    "nmad",
    "tails",
    "bend",
    "scan_end_mass"
  )] <- NA
  expect_length(diagnostic_problems(d, vb_ok), 0)
  d$vb_shift[1] <- 3
  expect_length(diagnostic_problems(d, vb_ok), 1)
})

test_that("a heavy-tailed posterior warns once at the end of the fit", {
  w <- expect_warning(
    fit <- suppressMessages(inla_vanilla(lp_cauchy, c(theta = 0.1))),
    class = "inlavanilla_diagnostics_warning"
  )
  expect_s3_class(w, "inlavanilla_warning")
  text <- cli::ansi_strip(conditionMessage(w))
  expect_match(text, "Heavy tails for theta")
  expect_match(text, "Poor skew-normal fit for theta")
  expect_length(fit$warnings, 1)
  expect_match(fit$warnings, "^The approximation may be poor")

  # silent = TRUE: recorded, not signalled
  expect_silent(fit <- inla_vanilla(lp_cauchy, c(theta = 0.1), silent = TRUE))
  expect_identical(
    gsub("\\s+", " ", fit$warnings),
    gsub("\\s+", " ", text)
  )
})

test_that("summary() flags a scan-end mass above 0.05", {
  testthat::local_reproducible_output(width = 200, crayon = FALSE)
  fit <- inla_vanilla(lp_reg, init_reg, silent = TRUE)
  out <- capture_output(print(summary(fit)))
  expect_no_match(out, "beyond the scanned window")
  fit$diagnostics$marginals$scan_end_mass[2] <- 0.2
  out <- capture_output(print(summary(fit)))
  expect_match(
    out,
    "The skew-normal for slope puts up to 0.2 of its mass beyond the scanned window"
  )
})
