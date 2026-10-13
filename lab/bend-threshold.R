# Run from the package root: Rscript lab/bend-threshold.R
# Reproduce the analytic calibration and counterexample in bend-threshold.md.
# Research only: does not change package thresholds or write output files.
devtools::load_all(".", quiet = TRUE)

# Normal data with mean zero and sum of squares n; prior flat in (mu, log_sd).
# The posterior is proper for n > 1 and its mu variance exists for n > 3.
# Fractional n is allowed as an interpolation of this analytic family.
normal_scale <- function(n) {
  lp <- function(th) {
    -n * th[2] - 0.5 * n * (1 + th[1]^2) * exp(-2 * th[2])
  }
  gr <- function(th) {
    c(-n * th[1] * exp(-2 * th[2]), -n + n * (1 + th[1]^2) * exp(-2 * th[2]))
  }
  fit <- inla_vanilla(
    lp,
    c(mu = 0, log_sd = 0),
    grad = gr,
    vb = FALSE,
    silent = TRUE
  )
  # Exact marginal: sqrt(n - 1) * mu ~ t_(n - 1).
  marginal <- function(x) sqrt(n - 1) * dt(x * sqrt(n - 1), df = n - 1)
  sd_quad <- sqrt(
    integrate(function(x) x^2 * marginal(x), -Inf, Inf, rel.tol = 1e-10)$value
  )
  sd_exact <- 1 / sqrt(n - 3)
  b <- unname(fit$diagnostics$marginals["mu", "bend"])
  sd_fit <- unname(fit$summary["mu", "SD"])
  stopifnot(
    abs(b - 4 / n) < 1e-5,
    abs(sd_quad / sd_exact - 1) < 1e-8,
    abs(sd_fit * sqrt(n) - 1) < 1e-4
  )
  data.frame(
    n = n,
    bend = b,
    sd_shortfall = 1 - sd_fit / sd_exact,
    analytic_shortfall = 1 - sqrt(1 - 3 / n)
  )
}
cat("Normal mean / unknown scale calibration:\n")
print(
  do.call(rbind, lapply(c(12, 20, 4 / 0.13, 40, 80, 200), normal_scale)),
  digits = 7,
  row.names = FALSE
)

eps <- c(0.01, 0.025, 0.05, 0.10)
cat("\nThresholds for specified relative SD shortfall:\n")
print(
  data.frame(sd_tolerance = eps, bend_threshold = (4 / 3) * (2 * eps - eps^2)),
  digits = 7,
  row.names = FALSE
)

# Conditional mode y = 0 is perfectly straight, but conditional variance
# exp(a * x^2) contributes a quadratic volume term to the x marginal.
# The normalized x marginal is N(0, 1 / (1 - a)); the scanned marginal N(0, 1).
a <- 0.25
lp_volume <- function(th) {
  -0.5 * th[1]^2 - 0.5 * th[2]^2 * exp(-a * th[1]^2)
}
gr_volume <- function(th) {
  c(
    -th[1] + a * th[1] * th[2]^2 * exp(-a * th[1]^2),
    -th[2] * exp(-a * th[1]^2)
  )
}
volume_fit <- inla_vanilla(
  lp_volume,
  c(x = 0, y = 0),
  grad = gr_volume,
  vb = FALSE,
  silent = TRUE
)
volume_bend <- unname(volume_fit$diagnostics$marginals["x", "bend"])
volume_sd <- unname(volume_fit$summary["x", "SD"])
stopifnot(abs(volume_bend) < 1e-10, abs(volume_sd - 1) < 1e-4)
cat("\nZero-bend counterexample:\n")
print(
  data.frame(
    a = a,
    bend = volume_bend,
    sd_fit = volume_sd,
    sd_exact = 1 / sqrt(1 - a),
    sd_shortfall = 1 - volume_sd * sqrt(1 - a)
  ),
  digits = 7,
  row.names = FALSE
)
