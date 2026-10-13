# Unit tests of the skew-normal utilities in R/skew-normal.R and
# R/qsnorm-fast.R. Every reference is computed here, independently of the code
# under test: stats::integrate() of the density or of Owen's defining
# integral, closed forms (alpha = 0 and |alpha| = 1), or root finding on a
# score written out below. Quantile references invert psnorm(), which is
# itself checked against integrate() earlier in this file.

# Helpers ----------------------------------------------------------------------

# Standardised skew-normal density 2 phi(z) Phi(alpha z), written from scratch.
sn_density_ref <- function(z, alpha) {
  2 * stats::dnorm(z) * stats::pnorm(alpha * z)
}

# Log-density of SN(xi, omega, alpha), written from scratch.
sn_logdens_ref <- function(x, xi, omega, alpha) {
  z <- (x - xi) / omega
  log(2) +
    stats::dnorm(z, log = TRUE) +
    stats::pnorm(alpha * z, log.p = TRUE) -
    log(omega)
}

# Standardised distribution function by quadrature. The upper tail is
# integrated directly, so small upper-tail probabilities keep their relative
# accuracy.
sn_cdf_ref <- function(z, alpha, lower_tail = TRUE) {
  f <- function(x) sn_density_ref(x, alpha)
  if (lower_tail) {
    stats::integrate(f, -Inf, z, rel.tol = 1e-12, abs.tol = 0)$value
  } else {
    stats::integrate(f, z, Inf, rel.tol = 1e-12, abs.tol = 0)$value
  }
}

# Owen's T(h, a) from its defining integral.
owen_t_ref <- function(h, a) {
  if (a == 0) {
    return(0)
  }
  f <- function(x) exp(-h^2 * (1 + x^2) / 2) / (1 + x^2)
  stats::integrate(f, 0, a, rel.tol = 1e-13, abs.tol = 0)$value / (2 * pi)
}

# Mills ratio phi(u) / Phi(u), computed on the log scale.
mills_ref <- function(u) {
  exp(stats::dnorm(u, log = TRUE) - stats::pnorm(u, log.p = TRUE))
}

# Mode of the standardised skew-normal: the root of its score
# -z + alpha phi(alpha z) / Phi(alpha z).
sn_std_mode_ref <- function(alpha) {
  score <- function(z) -z + alpha * mills_ref(alpha * z)
  stats::uniroot(score, c(-3, 3), tol = 1e-14)$root
}

# SN parameters as the marginal scan sees them: the scan coordinate is centred
# at the mode and scaled by the Laplace SD, so the mode is at `mode` and the
# curvature of the log-density there is -1.
scan_sn_pars <- function(alpha, mode = 0) {
  z0 <- sn_std_mode_ref(alpha)
  u <- alpha * z0
  omega <- sqrt(1 + alpha^2 * mills_ref(u) * (u + mills_ref(u)))
  c(xi = mode - omega * z0, omega = omega, alpha = alpha)
}

# Reference quantile: the root of psnorm(q) = p, on the standardised scale.
sn_quantile_ref <- function(p, alpha) {
  vapply(
    p,
    function(pp) {
      stats::uniroot(
        function(z) psnorm(z, 0, 1, alpha) - pp,
        c(-40, 40),
        tol = 1e-13,
        maxiter = 500
      )$root
    },
    numeric(1)
  )
}

# Expect every element of a named vector of errors to be below `tol`, and on
# failure list the worst offenders by name.
expect_all_below <- function(err, tol) {
  bad <- sort(err[is.na(err) | err >= tol], decreasing = TRUE)
  shown <- utils::head(bad, 6)
  expect(
    length(bad) == 0L,
    sprintf(
      "%d of %d errors are not below %g. Worst: %s",
      length(bad),
      length(err),
      tol,
      paste(sprintf("[%s] %.2g", names(shown), shown), collapse = ", ")
    )
  )
  invisible(err)
}

grid_scan <- seq(-4, 4, length.out = 21)

# dsnorm ------------------------------------------------------------------------

sn_pars <- list(
  c(0, 1, 0),
  c(1, 2, 3),
  c(-1, 0.5, -5),
  c(0, 1, 20),
  c(2, 3, -20),
  c(0.5, 1.5, 0.7),
  c(0, 1, -1)
)

test_that("dsnorm integrates to one", {
  # Observed |integral - 1| <= 1.1e-16; integrate() itself reports ~1e-12.
  for (p in sn_pars) {
    total <- stats::integrate(
      function(x) dsnorm(x, p[1], p[2], p[3]),
      -Inf,
      Inf,
      rel.tol = 1e-12
    )$value
    expect_equal(total, 1, tolerance = 1e-10)
  }
})

test_that("dsnorm on the log scale agrees with the density", {
  x <- seq(-6, 6, by = 0.25)
  for (p in sn_pars) {
    d <- dsnorm(x, p[1], p[2], p[3])
    log_d <- dsnorm(x, p[1], p[2], p[3], log = TRUE)
    ok <- d > 1e-300
    expect_equal(log_d[ok], log(d[ok]), tolerance = 1e-12)
    expect_equal(log_d, sn_logdens_ref(x, p[1], p[2], p[3]), tolerance = 1e-12)
  }

  # The log-density stays finite where the density itself underflows
  expect_equal(dsnorm(-30, 0, 1, 5), 0)
  expect_true(is.finite(dsnorm(-30, 0, 1, 5, log = TRUE)))
  expect_equal(
    dsnorm(-30, 0, 1, 5, log = TRUE),
    sn_logdens_ref(-30, 0, 1, 5),
    tolerance = 1e-12
  )
})

test_that("dsnorm adds logC on the log scale and multiplies by exp(logC)", {
  x <- seq(-3, 3, by = 0.5)
  expect_equal(
    dsnorm(x, 1, 2, 3, logC = 0.7),
    exp(0.7) * dsnorm(x, 1, 2, 3),
    tolerance = 1e-14
  )
  expect_equal(
    dsnorm(x, 1, 2, 3, logC = 0.7, log = TRUE),
    0.7 + dsnorm(x, 1, 2, 3, log = TRUE),
    tolerance = 1e-14
  )
})

test_that("dsnorm with alpha = 0 is the normal density", {
  x <- seq(-5, 5, by = 0.1)
  expect_equal(dsnorm(x, 1, 2, 0), stats::dnorm(x, 1, 2), tolerance = 1e-14)
  expect_equal(
    dsnorm(x, 1, 2, 0, log = TRUE),
    stats::dnorm(x, 1, 2, log = TRUE),
    tolerance = 1e-14
  )
})

test_that("dsnorm is reflected about xi when alpha changes sign", {
  x <- seq(-5, 5, by = 0.1)
  expect_equal(
    dsnorm(x, 0.5, 1.5, -3),
    dsnorm(1 - x, 0.5, 1.5, 3),
    tolerance = 1e-12
  )
})

# owen_t ------------------------------------------------------------------------

test_that("owen_t matches its defining integral, including a > 1, a < 0, h < 0 and h = 0", {
  grid <- expand.grid(
    h = c(-6, -3, -1.5, -0.3, 0, 0.3, 1, 2.5, 5, 8),
    a = c(-25, -3, -1, -0.5, 0, 0.2, 0.7, 1, 1.5, 4, 30)
  )
  ref <- mapply(owen_t_ref, grid$h, grid$a)
  got <- owen_t(grid$h, grid$a)
  err <- abs(got - ref)
  names(err) <- sprintf("h=%g, a=%g", grid$h, grid$a)
  # Documented absolute accuracy ~1e-15; observed max 1.2e-16.
  expect_all_below(err, 1e-15)
})

test_that("owen_t reproduces its closed forms at h = 0, a = 0 and a = 1", {
  a <- c(-50, -3, -1, -0.2, 0, 0.5, 1, 2, 50)
  expect_equal(owen_t(0, a), atan(a) / (2 * pi), tolerance = 1e-14)

  h <- c(-4, -2, -0.5, 0, 0.5, 1, 3, 6)
  expect_equal(owen_t(h, 0), rep(0, length(h)))
  expect_lt(
    max(abs(owen_t(h, 1) - stats::pnorm(h) * stats::pnorm(-h) / 2)),
    1e-15
  )
})

# psnorm ------------------------------------------------------------------------

test_that("psnorm matches integrate() of the density, in both tails", {
  grid <- expand.grid(
    z = c(-4, -2, -1, -0.3, 0, 0.4, 1, 2.5, 4),
    alpha = c(-20, -4, -1, -0.5, 0, 0.5, 1, 2, 5, 20)
  )
  xi <- 0.5
  omega <- 2
  q <- xi + omega * grid$z
  lower_ref <- mapply(sn_cdf_ref, grid$z, grid$alpha)
  upper_ref <- mapply(sn_cdf_ref, grid$z, grid$alpha, lower_tail = FALSE)
  lower <- psnorm(q, xi, omega, grid$alpha)
  upper <- psnorm(q, xi, omega, grid$alpha, lower_tail = FALSE)
  labels <- sprintf("z=%g, alpha=%g", grid$z, grid$alpha)

  # Observed max absolute error 2.9e-14 (lower), 5.6e-16 (upper); the
  # reference is good to about 1e-12.
  expect_all_below(stats::setNames(abs(lower - lower_ref), labels), 1e-12)
  expect_all_below(stats::setNames(abs(upper - upper_ref), labels), 1e-12)
  expect_all_below(stats::setNames(abs(lower + upper - 1), labels), 1e-15)
})

test_that("psnorm stays in [0, 1] and absolutely accurate far into both tails", {
  grid <- expand.grid(
    z = c(5, 7, 8, 9, 12, 20, 40),
    alpha = c(-20, -5, -2, -1, -0.5, 0, 0.5, 1, 2, 5, 20)
  )
  upper <- psnorm(grid$z, 0, 1, grid$alpha, lower_tail = FALSE)
  lower <- psnorm(-grid$z, 0, 1, grid$alpha)
  expect_true(all(upper >= 0 & upper <= 1))
  expect_true(all(lower >= 0 & lower <= 1))

  upper_ref <- mapply(sn_cdf_ref, grid$z, grid$alpha, lower_tail = FALSE)
  lower_ref <- mapply(sn_cdf_ref, -grid$z, grid$alpha)
  labels <- sprintf("z=%g, alpha=%g", grid$z, grid$alpha)
  # Bounded by owen_t's documented absolute accuracy (~1e-15); observed max
  # 1.6e-16.
  expect_all_below(stats::setNames(abs(upper - upper_ref), labels), 1e-15)
  expect_all_below(stats::setNames(abs(lower - lower_ref), labels), 1e-15)
})

test_that("psnorm keeps its relative accuracy for small probabilities in the heavier tail", {
  # The heavier tail is the upper one for alpha >= 0 and the lower one for
  # alpha < 0; there 1 - F is between Phi(-|z|) and 2 Phi(-|z|), so a double
  # can hold it to full relative precision, and psnorm() computes the upper
  # tail by reflection precisely so that it does not cancel. Observed
  # relative error <= 5e-15 for |alpha| <= 1. For |alpha| > 1 the observed
  # error is 2.7e-10 at z = 5, 2.8e-8 at z = 6, 2.1e-5 at z = 7, 3.5e-2 at
  # z = 8 and 0.5 (a factor of 2) from z = 9 on: owen_t()'s a > 1 branch
  # forms Phi(h)/2 + Phi(ah)/2 - Phi(h) Phi(ah), a difference of numbers
  # close to one.
  grid <- expand.grid(
    z = c(5, 6, 7, 8, 9, 10),
    alpha = c(-20, -5, -2, -1, -0.5, 0, 0.5, 1, 2, 5, 20)
  )
  heavy_upper <- grid$alpha >= 0
  got <- ifelse(
    heavy_upper,
    psnorm(grid$z, 0, 1, grid$alpha, lower_tail = FALSE),
    psnorm(-grid$z, 0, 1, grid$alpha)
  )
  ref <- ifelse(
    heavy_upper,
    mapply(sn_cdf_ref, grid$z, grid$alpha, lower_tail = FALSE),
    mapply(sn_cdf_ref, -grid$z, grid$alpha)
  )
  err <- abs(got / ref - 1)
  names(err) <- sprintf("z=%g, alpha=%g", grid$z, grid$alpha)
  expect_all_below(err, 1e-9)
})

test_that("psnorm is 0 at -Inf and 1 at Inf", {
  expect_equal(psnorm(c(-Inf, Inf), 0, 1, 2), c(0, 1))
  expect_equal(psnorm(c(-Inf, Inf), 0, 1, -2, lower_tail = FALSE), c(1, 0))
  expect_equal(psnorm(c(-Inf, Inf), 0, 1, 0), c(0, 1))
})

test_that("psnorm returns NA for a non-positive or missing scale", {
  expect_equal(psnorm(0.3, 0, c(0, -1, NA, NaN), 2), rep(NA_real_, 4))
  expect_equal(psnorm(NA_real_, 0, 1, 2), NA_real_)

  # Valid entries next to invalid ones are still computed
  out <- psnorm(c(0.3, 0.3), 0, c(1, -1), 2)
  expect_equal(out[1], sn_cdf_ref(0.3, 2), tolerance = 1e-10)
  expect_true(is.na(out[2]))
})

# qsnorm_fast and qsnorm --------------------------------------------------------

p_wide <- c(
  1e-6,
  1e-5,
  1e-4,
  1e-3,
  0.01,
  0.025,
  0.1,
  0.25,
  0.5,
  0.75,
  0.9,
  0.975,
  0.99,
  0.999,
  1 - 1e-4,
  1 - 1e-5,
  1 - 1e-6
)
alpha_wide <- c(-20, -5, -2, -0.5, 0.3, 0.5, 2, 5, 20)

test_that("qsnorm_fast and qsnorm are exact where the quantile has a closed form", {
  # alpha = 0 is N(xi, omega^2); alpha = 1 has F = Phi^2 and alpha = -1 has
  # 1 - F = Phi(-z)^2.
  xi <- 0.4
  omega <- 1.7
  for (qfun in list(qsnorm_fast, qsnorm)) {
    expect_equal(
      qfun(p_wide, xi, omega, 0),
      stats::qnorm(p_wide, xi, omega),
      tolerance = 1e-12
    )
    expect_equal(
      qfun(p_wide, xi, omega, 1),
      xi + omega * stats::qnorm(sqrt(p_wide)),
      tolerance = 1e-12
    )
    expect_equal(
      qfun(p_wide, xi, omega, -1),
      xi - omega * stats::qnorm(sqrt(1 - p_wide)),
      tolerance = 1e-12
    )
    expect_equal(qfun(c(0, 1), xi, omega, 3), c(-Inf, Inf))
  }
})

test_that("qsnorm_fast scales, shifts and recycles its arguments", {
  p <- c(0.1, 0.5, 0.9)
  expect_equal(
    qsnorm_fast(p, 2, 3, -4),
    2 + 3 * qsnorm_fast(p, 0, 1, -4),
    tolerance = 1e-14
  )
  alpha <- c(-2, 0, 2)
  expect_equal(
    qsnorm_fast(p, 0, 1, alpha),
    vapply(1:3, function(i) qsnorm_fast(p[i], 0, 1, alpha[i]), numeric(1))
  )
  # SN(-alpha) is the mirror image of SN(alpha)
  expect_equal(
    qsnorm_fast(p, 0, 1, -2.5),
    -qsnorm_fast(1 - p, 0, 1, 2.5),
    tolerance = 1e-12
  )
})

test_that("qsnorm_fast inverts psnorm to its documented relative accuracy", {
  # R/qsnorm-fast.R documents a relative error below about 1e-7; the error
  # here is relative to max(1, |q|) on the standardised scale, so that it is
  # defined where q crosses 0, and the tolerance allows 10x headroom. Observed
  # errors are far larger: up to 2.7e-2 in the tails (p < 0.01 or > 0.99;
  # 6e-2 at alpha = 1.5, p = 1e-3), mostly from the Lambert W branch, and up
  # to 5e-3 in between, from the Taylor branch and the half-normal tail
  # switch at erf(qnorm(0.99) / sqrt(2) / alpha). Luu's original code sets
  # tol = 0.01: it is a starting value for iterative refinement, not a 1e-7
  # approximation.
  err <- NULL
  for (a in alpha_wide) {
    q_ref <- sn_quantile_ref(p_wide, a)
    e <- abs(qsnorm_fast(p_wide, 0, 1, a) - q_ref) / pmax(1, abs(q_ref))
    names(e) <- sprintf("alpha=%g, p=%g", a, p_wide)
    err <- c(err, e)
  }
  expect_all_below(err, 1e-6)
})

test_that("qsnorm inverts psnorm accurately at the probabilities summaries use", {
  # Observed max error 1.8e-10 over alpha in [-30, 30] at these p.
  probs <- c(0.025, 0.25, 0.5, 0.75, 0.975)
  err <- NULL
  for (a in seq(-30, 30, by = 0.5)) {
    q_ref <- sn_quantile_ref(probs, a)
    e <- abs(qsnorm(probs, 0, 1, a) - q_ref) / pmax(1, abs(q_ref))
    names(e) <- sprintf("alpha=%g, p=%g", a, probs)
    err <- c(err, e)
  }
  expect_all_below(err, 1e-9)

  # Location and scale carry through
  expect_equal(
    qsnorm(probs, -1, 2.5, 3),
    -1 + 2.5 * sn_quantile_ref(probs, 3),
    tolerance = 1e-10
  )
})

test_that("qsnorm inverts psnorm to near machine precision across (1e-6, 1 - 1e-6)", {
  # R/skew-normal.R documents qsnorm() as accurate "to near machine
  # precision". The tolerance allows for the conditioning of p near 1, where
  # the reference itself is good to about 1e-11. Observed: <= 7.7e-13 for
  # alpha = 0.3 (and <= 3.7e-12 for alpha in {0, -1, 1}), but up to 1.9e-4
  # in the light tail (p <= 1e-3 for alpha > 0) for every other alpha here
  # (2.9e-4 at alpha = 1.5): two Newton steps do not repair qsnorm_fast()'s
  # starting error of up to 6e-2 there.
  err <- NULL
  for (a in alpha_wide) {
    q_ref <- sn_quantile_ref(p_wide, a)
    e <- abs(qsnorm(p_wide, 0, 1, a) - q_ref) / pmax(1, abs(q_ref))
    names(e) <- sprintf("alpha=%g, p=%g", a, p_wide)
    err <- c(err, e)
  }
  expect_all_below(err, 1e-9)
})

# sn_moments and sn_mode --------------------------------------------------------

test_that("sn_moments gives the mean and SD found by integration", {
  # Observed errors <= 4.4e-16.
  for (p in sn_pars) {
    f <- function(x) dsnorm(x, p[1], p[2], p[3])
    mean_ref <- stats::integrate(
      function(x) x * f(x),
      -Inf,
      Inf,
      rel.tol = 1e-12
    )$value
    var_ref <- stats::integrate(
      function(x) (x - mean_ref)^2 * f(x),
      -Inf,
      Inf,
      rel.tol = 1e-12
    )$value
    mom <- sn_moments(p[1], p[2], p[3])
    expect_equal(mom$mean, mean_ref, tolerance = 1e-9)
    expect_equal(mom$sd, sqrt(var_ref), tolerance = 1e-9)
  }
})

test_that("sn_moments is vectorised over the parameters", {
  xi <- c(0, 1, -2)
  omega <- c(1, 2, 0.5)
  alpha <- c(0, 3, -7)
  mom <- sn_moments(xi, omega, alpha)
  for (i in 1:3) {
    one <- sn_moments(xi[i], omega[i], alpha[i])
    expect_equal(mom$mean[i], one$mean)
    expect_equal(mom$sd[i], one$sd)
  }
})

test_that("sn_mode finds the maximum of the density", {
  # Reference 1: the root of the score; reference 2: optimize() on the
  # density. Observed errors <= 9.7e-9 and <= 1.5e-8 SDs, the precision of
  # optimize().
  for (p in sn_pars) {
    mode <- sn_mode(p[1], p[2], p[3])
    sd <- sn_moments(p[1], p[2], p[3])$sd
    mode_score <- p[1] + p[2] * sn_std_mode_ref(p[3])
    mode_opt <- stats::optimize(
      function(x) sn_density_ref((x - p[1]) / p[2], p[3]),
      p[1] + c(-4, 4) * p[2],
      maximum = TRUE,
      tol = 1e-12
    )$maximum
    expect_lt(abs(mode - mode_score) / sd, 1e-7)
    expect_lt(abs(mode - mode_opt) / sd, 1e-7)
  }
})

test_that("sn_mode is xi when alpha = 0", {
  # Observed error 1.7e-8 omega.
  for (p in list(c(0, 1), c(3, 0.2), c(-5, 4))) {
    expect_lt(abs(sn_mode(p[1], p[2], 0) - p[1]) / p[2], 1e-7)
  }
})

# fit_skew_normal ---------------------------------------------------------------

fit_sn_check <- function(p, x = grid_scan) {
  y <- sn_logdens_ref(x, p[1], p[2], p[3])
  y <- y - max(y)
  fit <- fit_skew_normal(x, y)
  fit$y <- y
  fit
}

test_that("fit_skew_normal recovers xi, omega and alpha from exact skew-normal scans", {
  pars <- c(
    # As the scan sees them: mode 0, unit curvature at the mode
    lapply(c(-20, -5, -2, -0.5, 0.5, 2, 5, 20), scan_sn_pars),
    # Arbitrary location and scale
    list(
      c(-0.6, 1.3, 1.5),
      c(0.9, 1.4, -2),
      c(-1.2, 1.8, 4),
      c(1.1, 1.6, -4),
      c(0.3, 0.9, 0.5),
      c(-0.8, 1.5, 8)
    )
  )
  # Observed parameter errors <= 2e-14 and nmad <= 4.4e-16.
  for (p in pars) {
    fit <- fit_sn_check(p)
    expect_equal(
      c(fit$xi, fit$omega, fit$alpha),
      unname(p),
      tolerance = 1e-8
    )
    expect_lt(fit$nmad, 1e-10)
    # logC restores the normalisation removed by max(y) = 0
    expect_equal(
      dsnorm(grid_scan, fit$xi, fit$omega, fit$alpha, fit$logC, log = TRUE),
      fit$y,
      tolerance = 1e-8
    )
  }
})

test_that("fit_skew_normal fits mirror-image inputs with mirror-image parameters", {
  # On a grid symmetric about 0, rev(y) is the log-density of
  # SN(-xi, omega, -alpha). Observed: exact for alpha = +-2, but for
  # alpha = +-0.2 one of the pair recovers alpha and the other collapses to
  # alpha ~ 0 (nmad 3.8e-4). The trapezoidal weights in fit_skew_normal()
  # are off by one: the second grid point gets weight h / 2 instead of h, so
  # the moment-based starting values are not mirror-symmetric.
  for (a in c(-2, -0.2, 0.2, 2)) {
    p <- scan_sn_pars(a)
    y <- sn_logdens_ref(grid_scan, p[1], p[2], p[3])
    y <- y - max(y)
    fit <- fit_skew_normal(grid_scan, y)
    fit_rev <- fit_skew_normal(grid_scan, rev(y))
    expect_equal(
      c(fit_rev$xi, fit_rev$omega, fit_rev$alpha),
      c(-fit$xi, fit$omega, -fit$alpha),
      tolerance = 1e-6,
      label = sprintf("mirrored fit for alpha = %g", a)
    )
  }
})

test_that("fit_skew_normal recovers mildly skewed scans, centred or shifted", {
  # The volume correction adds gamma1 * t to the scan, which moves its mode
  # off t = 0, so shifted inputs are routine. The input is an exact
  # skew-normal, so the fit should be exact (as it is for |alpha| >= 0.5,
  # observed errors ~1e-14). Observed: 9 of these 24 cases collapse to a
  # Gaussian fit with |alpha| < 1e-4 and nmad up to 5.3e-3, because the
  # moment-based starting alpha has the wrong sign (the off-by-one
  # trapezoidal weight, and truncation of the scan at +-4) and nlminb() stops
  # at the stationary point that the direct parametrisation has at alpha = 0.
  cases <- expand.grid(
    alpha = c(-0.5, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.5),
    mode = c(-0.5, 0, 0.5)
  )
  err <- numeric(nrow(cases))
  for (i in seq_len(nrow(cases))) {
    fit <- fit_sn_check(scan_sn_pars(cases$alpha[i], cases$mode[i]))
    err[i] <- abs(fit$alpha - cases$alpha[i])
  }
  names(err) <- sprintf("alpha=%g, mode=%g", cases$alpha, cases$mode)
  expect_all_below(err, 1e-6)
})

test_that("fit_skew_normal gives alpha near 0 and the right moments for a Gaussian input", {
  # Near alpha = 0 the shape is weakly identified (the skewness is O(alpha^3)),
  # so alpha itself stops short of 0. Observed |alpha| <= 0.0063, moment
  # errors <= 1e-8 and nmad <= 1e-8.
  for (p in list(c(0, 1), c(0.4, 1.2), c(-0.7, 0.8), c(1, 1.5))) {
    y <- stats::dnorm(grid_scan, p[1], p[2], log = TRUE)
    y <- y - max(y)
    fit <- fit_skew_normal(grid_scan, y)
    mom <- sn_moments(fit$xi, fit$omega, fit$alpha)
    expect_lt(abs(fit$alpha), 0.02)
    expect_lt(abs(mom$mean - p[1]), 5e-8)
    expect_lt(abs(mom$sd - p[2]), 5e-8)
    expect_lt(fit$nmad, 5e-8)
  }
})

test_that("fit_skew_normal ignores points below threshold_log_drop", {
  # Six of the 21 points fall below -6 for this scan; pushing them further
  # down must not move the fit.
  p <- scan_sn_pars(5)
  y <- sn_logdens_ref(grid_scan, p[1], p[2], p[3])
  y <- y - max(y)
  low <- y < -6
  expect_gt(sum(low), 0)
  y[low] <- y[low] - 50
  fit <- fit_skew_normal(grid_scan, y)
  expect_equal(c(fit$xi, fit$omega, fit$alpha), unname(p), tolerance = 1e-8)
})
