# Unit tests of the skew-normal utilities in R/skew-normal.R and
# R/qsnorm-fast.R. Every reference is computed here, independently of the code
# under test: stats::integrate() of the density or of Owen's defining
# integral, closed forms (alpha = 0 and |alpha| = 1), root finding on a score
# or on a quadrature CDF written out below, and optim() on a least-squares
# objective written out below. Quantiles are checked by quadrature at the
# returned value, not by inverting psnorm(); only qsnorm_fast()'s coarse
# envelope uses a psnorm() root, and psnorm() is itself checked against
# integrate() first.
#
# Tolerances were calibrated by printing the actual errors first; each
# comment states the observed value.

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

# SD of the standardised skew-normal, 1 - 2 delta^2 / pi with
# delta = alpha / sqrt(1 + alpha^2) (checked against quadrature below).
sn_sd_std <- function(alpha) {
  sqrt(1 - 2 / pi * alpha^2 / (1 + alpha^2))
}

# Standardised distribution function by quadrature, lower or upper tail.
# Each tail is integrated directly, from its infinite end and split at 0
# (the density's bulk), so small tail probabilities keep their relative
# accuracy and probabilities near 1/2 their absolute accuracy.
sn_cdf_ref <- function(z, alpha, lower_tail = TRUE) {
  f <- function(x) sn_density_ref(x, alpha)
  integ <- function(lo, hi) {
    if (lo >= hi) {
      return(0)
    }
    stats::integrate(f, lo, hi, rel.tol = 1e-13, abs.tol = 0)$value
  }
  if (lower_tail) {
    integ(-Inf, min(z, 0)) + integ(min(z, 0), z)
  } else {
    integ(max(z, 0), Inf) + integ(z, max(z, 0))
  }
}

# Error of a standardised quantile q of probability p, from one Newton step
# on the log of the smaller tail: |log F(q) - log p| F(q) / f(q), with F by
# quadrature. To first order this is |F(q) - p| / f(q), and on the log scale
# it stays faithful when q is far out in a tail. Independent of psnorm(), and
# in standardised units. The target can be given as a tail probability
# (`tail`, on the side `upper`) when 1 - p would round.
sn_quantile_err <- function(q, p, alpha, upper = p > 0.5, tail = NULL) {
  if (is.null(tail)) {
    tail <- ifelse(upper, 1 - p, p)
  }
  mapply(
    function(qq, tt, up, a) {
      got <- sn_cdf_ref(qq, a, lower_tail = !up)
      abs(log(got) - log(tt)) * got / sn_density_ref(qq, a)
    },
    q,
    tail,
    upper,
    alpha
  )
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
# Used only for qsnorm_fast()'s coarse envelope.
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
  # Observed |integral - 1| <= 1.1e-16; the tolerance is what integrate() is
  # asked for.
  for (p in sn_pars) {
    total <- stats::integrate(
      function(x) dsnorm(x, p[1], p[2], p[3]),
      -Inf,
      Inf,
      rel.tol = 1e-12
    )$value
    expect_equal(total, 1, tolerance = 1e-12)
  }
})

test_that("dsnorm on the log scale agrees with the density and the reference", {
  # Observed: <= 3.5e-16 relative against log(density), exact against the
  # reference.
  x <- seq(-6, 6, by = 0.25)
  for (p in sn_pars) {
    d <- dsnorm(x, p[1], p[2], p[3])
    log_d <- dsnorm(x, p[1], p[2], p[3], log = TRUE)
    ok <- d > 1e-300
    expect_equal(log_d[ok], log(d[ok]), tolerance = 1e-14)
    expect_equal(log_d, sn_logdens_ref(x, p[1], p[2], p[3]), tolerance = 1e-14)
    expect_equal(
      d,
      sn_density_ref((x - p[1]) / p[2], p[3]) / p[2],
      tolerance = 1e-14
    )
  }

  # The log-density stays finite where the density itself underflows
  expect_equal(dsnorm(-30, 0, 1, 5), 0)
  expect_equal(
    dsnorm(-30, 0, 1, 5, log = TRUE),
    sn_logdens_ref(-30, 0, 1, 5),
    tolerance = 1e-14
  )
})

test_that("dsnorm adds logC on the log scale and multiplies by exp(logC)", {
  x <- seq(-3, 3, by = 0.5)
  expect_equal(
    dsnorm(x, 1, 2, 3, logC = 0.7),
    exp(0.7 + sn_logdens_ref(x, 1, 2, 3)),
    tolerance = 1e-14
  )
  expect_equal(
    dsnorm(x, 1, 2, 3, logC = 0.7, log = TRUE),
    0.7 + sn_logdens_ref(x, 1, 2, 3),
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
    tolerance = 1e-14
  )
})

# owen_t ------------------------------------------------------------------------

test_that("owen_t matches its defining integral, including a > 1, a < 0, h < 0 and h = 0", {
  grid <- expand.grid(
    h = c(-20, -12, -6, -3, -1.5, -0.3, 0, 0.3, 1, 2.5, 5, 8, 12, 20),
    a = c(-25, -3, -1, -0.5, 0, 0.2, 0.7, 1, 1.5, 4, 30)
  )
  ref <- mapply(owen_t_ref, grid$h, grid$a)
  got <- owen_t(grid$h, grid$a)
  labels <- sprintf("h=%g, a=%g", grid$h, grid$a)
  # Absolute: documented ~1e-15, observed max 4.2e-17.
  expect_all_below(stats::setNames(abs(got - ref), labels), 2e-16)
  # Relative, which psnorm()'s heavy tail needs: T(20, a) is about 1e-89.
  # Observed max 1.9e-14 (at h = 20), at the reference's own precision.
  nz <- ref != 0
  expect_all_below(
    stats::setNames(abs(got[nz] / ref[nz] - 1), labels[nz]),
    6e-14
  )
})

test_that("owen_t reproduces its closed forms at h = 0, a = 0 and a = 1", {
  # T(0, a) = atan(a) / (2 pi) and T(h, 1) = Phi(h) Phi(-h) / 2. Observed
  # 1.4e-17 absolute, 2.1e-15 relative (at h = 6).
  a <- c(-50, -3, -1, -0.2, 0, 0.5, 1, 2, 50)
  expect_equal(owen_t(0, a), atan(a) / (2 * pi), tolerance = 1e-15)

  h <- c(-4, -2, -0.5, 0, 0.5, 1, 3, 6)
  expect_identical(owen_t(h, 0), rep(0, length(h)))
  t1 <- stats::pnorm(h) * stats::pnorm(-h) / 2
  expect_lt(max(abs(owen_t(h, 1) - t1)), 1e-16)
  expect_lt(max(abs(owen_t(h, 1) / t1 - 1)), 6e-15)
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

  # Observed max absolute error 3.3e-16 (lower) and 2.2e-16 (upper), and
  # 1.1e-16 for lower + upper - 1.
  expect_all_below(stats::setNames(abs(lower - lower_ref), labels), 1e-15)
  expect_all_below(stats::setNames(abs(upper - upper_ref), labels), 1e-15)
  expect_all_below(stats::setNames(abs(lower + upper - 1), labels), 4.5e-16)
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
  # Bounded by rounding of Phi(-5) ~ 3e-7, the largest term here; observed
  # max 3.2e-22.
  expect_all_below(stats::setNames(abs(upper - upper_ref), labels), 1e-21)
  expect_all_below(stats::setNames(abs(lower - lower_ref), labels), 1e-21)
})

test_that("psnorm keeps its relative accuracy for small probabilities in the heavier tail", {
  # The heavier tail is the upper one for alpha >= 0 and the lower one for
  # alpha < 0; there the tail probability lies between Phi(-|z|) and
  # 2 Phi(-|z|), so a double holds it to full relative precision, and
  # psnorm() computes it without cancellation (owen_t()'s a > 1 identity in
  # upper-tail form). Observed relative error <= 3.2e-14, at z = 30, where
  # the tail is 1e-197.
  grid <- expand.grid(
    z = c(5, 6, 7, 8, 9, 10, 15, 20, 30),
    alpha = c(-20, -5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 5, 20)
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
  expect_all_below(err, 1e-13)
})

test_that("psnorm keeps enough relative accuracy in the lighter tail for qsnorm", {
  # In the lighter tail, F = Phi(z) - 2 T(z, alpha) is a difference of
  # numbers of size Phi(z), so relative accuracy is lost as F shrinks. It
  # still has to be good where qsnorm() works, down to p = 1e-6 and a little
  # beyond. For tail probabilities >= 1e-8, observed max 4.4e-11 (at
  # alpha = +-2, F = 1.8e-8).
  grid <- expand.grid(
    z = seq(0.25, 6, by = 0.25),
    alpha = c(-20, -5, -2, -1.5, -1, -0.5, -0.3, 0.3, 0.5, 1, 1.5, 2, 5, 20)
  )
  light_upper <- grid$alpha < 0
  got <- ifelse(
    light_upper,
    psnorm(grid$z, 0, 1, grid$alpha, lower_tail = FALSE),
    psnorm(-grid$z, 0, 1, grid$alpha)
  )
  ref <- ifelse(
    light_upper,
    mapply(sn_cdf_ref, grid$z, grid$alpha, lower_tail = FALSE),
    mapply(sn_cdf_ref, -grid$z, grid$alpha)
  )
  keep <- ref >= 1e-8
  expect_gt(sum(keep), 100)
  err <- abs(got / ref - 1)[keep]
  names(err) <- sprintf("z=%g, alpha=%g", grid$z, grid$alpha)[keep]
  expect_all_below(err, 1.5e-10)
})

test_that("psnorm is 0 at -Inf and 1 at Inf", {
  expect_identical(psnorm(c(-Inf, Inf), 0, 1, 2), c(0, 1))
  expect_identical(
    psnorm(c(-Inf, Inf), 0, 1, -2, lower_tail = FALSE),
    c(1, 0)
  )
  expect_identical(psnorm(c(-Inf, Inf), 0, 1, 0), c(0, 1))
})

test_that("psnorm returns NA for a non-positive or missing scale or shape", {
  expect_identical(psnorm(0.3, 0, c(0, -1, NA, NaN), 2), rep(NA_real_, 4))
  expect_identical(psnorm(NA_real_, 0, 1, 2), NA_real_)
  expect_identical(psnorm(0.3, 0, 1, c(NA, Inf)), rep(NA_real_, 2))

  # Valid entries next to invalid ones are still computed
  out <- psnorm(c(0.3, 0.3), 0, c(1, -1), 2)
  expect_equal(out[1], sn_cdf_ref(0.3, 2), tolerance = 1e-14)
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
alpha_wide <- c(-20, -5, -2, -1.5, -0.5, 0.3, 0.5, 1.5, 2, 5, 20)

test_that("qsnorm_fast and qsnorm are exact where the quantile has a closed form", {
  # alpha = 0 is N(xi, omega^2); alpha = 1 has F = Phi^2 and alpha = -1 has
  # 1 - F = Phi(-z)^2. Each reference is written in its well-conditioned
  # tail. Observed, in units of omega: qsnorm <= 2.2e-14; qsnorm_fast
  # <= 1.7e-11, since it takes qnorm(sqrt(p)), which is ill-conditioned near
  # p = 1.
  xi <- 0.4
  omega <- 1.7
  ref_pos <- xi +
    omega *
      ifelse(
        p_wide > 0.5,
        stats::qnorm((1 - p_wide) / (1 + sqrt(p_wide)), lower.tail = FALSE),
        stats::qnorm(sqrt(p_wide))
      )
  ref_neg <- xi -
    omega *
      ifelse(
        p_wide < 0.5,
        stats::qnorm(p_wide / (1 + sqrt(1 - p_wide)), lower.tail = FALSE),
        stats::qnorm(sqrt(1 - p_wide))
      )
  for (case in list(list(qsnorm_fast, 5e-11), list(qsnorm, 6e-14))) {
    qfun <- case[[1]]
    tol <- case[[2]]
    expect_equal(
      qfun(p_wide, xi, omega, 0),
      stats::qnorm(p_wide, xi, omega),
      tolerance = 1e-15
    )
    expect_lt(max(abs(qfun(p_wide, xi, omega, 1) - ref_pos) / omega), tol)
    expect_lt(max(abs(qfun(p_wide, xi, omega, -1) - ref_neg) / omega), tol)
    expect_identical(qfun(c(0, 1), xi, omega, 3), c(-Inf, Inf))
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

test_that("qsnorm_fast is within its documented envelope of the true quantile", {
  # R/qsnorm-fast.R documents qsnorm_fast() as a starting value whose error
  # reaches about 5e-2 in standardised units away from alpha = 0 and +-1
  # (Luu's original uses a tolerance of 0.01). Errors here are relative to
  # max(1, |q|). Observed: max 0.044 overall (alpha = +-1.5, p = 1e-3 or
  # 1 - 1e-3), from the Lambert W tail branch, and max 5.0e-3 for p in
  # [0.01, 0.99], from the Taylor branch. (On a finer grid of alpha the
  # overall max is 0.055, and 0.079 in absolute terms, at |alpha| ~ 1.25.)
  err <- NULL
  for (a in alpha_wide) {
    q_ref <- sn_quantile_ref(p_wide, a)
    e <- abs(qsnorm_fast(p_wide, 0, 1, a) - q_ref) / pmax(1, abs(q_ref))
    names(e) <- sprintf("alpha=%g, p=%g", a, p_wide)
    err <- c(err, e)
  }
  expect_all_below(err, 0.1)
  bulk <- rep(p_wide >= 0.01 & p_wide <= 0.99, length(alpha_wide))
  expect_all_below(err[bulk], 0.0125)
})

test_that("qsnorm_fast gives the same start in scalar and vector calls for tiny p", {
  # For alpha < 0 the function folds p to 1 - p, which is exactly 1 for
  # p < 1.1e-16. The special values for p = 0 and 1 are assigned before the
  # fold, so a scalar call returns -Inf and, in a vector call, the element
  # keeps its initial value 0. The true quantiles (heavy lower tail) are
  # about -8.5, -9.3 and -21.3. Observed: -Inf (scalar) and 0 (vector).
  p <- c(1e-17, 1e-20, 1e-100)
  for (a in c(-2, -5)) {
    scalar <- vapply(p, function(pp) qsnorm_fast(pp, 0, 1, a), numeric(1))
    vector <- qsnorm_fast(c(0.5, p), 0, 1, a)[-1]
    q_true <- sn_quantile_ref(p, a)
    expect_identical(vector, scalar)
    expect_all_below(
      stats::setNames(
        abs(vector - q_true) / abs(q_true),
        sprintf("alpha=%g, p=%g", a, p)
      ),
      0.1
    )
  }
})

test_that("qsnorm returns quantiles whose probability is p, across (1e-6, 1 - 1e-6)", {
  # R/skew-normal.R documents qsnorm() as accurate to near machine
  # precision. The error is |F(q) - p| / f(q) with F by quadrature on the
  # smaller tail, relative to max(1, |q|). Observed max 6.1e-13 (at alpha = 5,
  # p = 1e-6, the light tail); 4.4e-16 for alpha in {0, +-1}.
  alphas <- c(alpha_wide, -1, 0, 1)
  grid <- expand.grid(p = p_wide, alpha = alphas)
  q <- unlist(lapply(alphas, function(a) qsnorm(p_wide, 0, 1, a)))
  err <- sn_quantile_err(q, grid$p, grid$alpha) / pmax(1, abs(q))
  names(err) <- sprintf("alpha=%g, p=%g", grid$alpha, grid$p)
  expect_all_below(err, 2e-12)
})

test_that("qsnorm is accurate at the probabilities summaries use, over alpha in [-30, 30]", {
  probs <- c(0.025, 0.25, 0.5, 0.75, 0.975)
  alphas <- seq(-30, 30, by = 0.5)
  grid <- expand.grid(p = probs, alpha = alphas)
  q <- unlist(lapply(alphas, function(a) qsnorm(probs, 0, 1, a)))
  err <- sn_quantile_err(q, grid$p, grid$alpha) / pmax(1, abs(q))
  names(err) <- sprintf("alpha=%g, p=%g", grid$alpha, grid$p)
  # Observed max 7.4e-16.
  expect_all_below(err, 2.5e-15)

  # Location and scale carry through
  expect_equal(
    qsnorm(probs, -1, 2.5, 3),
    -1 + 2.5 * qsnorm(probs, 0, 1, 3),
    tolerance = 1e-15
  )
})

test_that("qsnorm recycles its arguments and returns NA or +-Inf at the edges", {
  p <- c(0.1, 0.5, 0.9)
  expect_equal(
    qsnorm(p, c(0, 1, 2), c(1, 2, 3), c(-2, 0, 5)),
    c(qsnorm(0.1, 0, 1, -2), qsnorm(0.5, 1, 2, 0), qsnorm(0.9, 2, 3, 5)),
    tolerance = 1e-15
  )
  expect_identical(
    qsnorm(c(-0.1, 1.1, NA, 0, 1), 0, 1, 2),
    c(NA, NA, NA, -Inf, Inf)
  )
  expect_identical(qsnorm(0.5, 0, c(0, -1, NA), 2), rep(NA_real_, 3))
  expect_identical(qsnorm(0.5, 0, 1, c(NA, Inf)), rep(NA_real_, 2))
})

test_that("qsnorm stays accurate far into the heavier tail", {
  # The heavier tail is the lower one for alpha < 0 (and the upper one for
  # alpha > 0, but p = 1 - 1e-20 is not representable). Probabilities in one
  # vector call, as summary() makes them. Observed error <= 1.7e-16 down to
  # p = 1e-100. At p = 1e-300 the true quantile is -37.05, but qsnorm()
  # returns -30.0: the vector start from qsnorm_fast() is 0 (see above), and
  # 30 Newton steps capped at one unit cannot reach it. (A scalar call starts
  # at -8 and returns -37.066, also short of convergence.)
  p <- c(0.5, 1e-10, 1e-20, 1e-50, 1e-100, 1e-300)
  for (a in c(-20, -5, -2, -0.5)) {
    q <- qsnorm(p, 0, 1, a)
    err <- sn_quantile_err(q, p, a) / pmax(1, abs(q))
    names(err) <- sprintf("alpha=%g, p=%g", a, p)
    expect_all_below(err, 5e-16)
  }
})

test_that("qsnorm does not run away deep in the lighter tail", {
  # In the lighter tail psnorm() cancels (see above), so qsnorm() cannot be
  # exact there. R/skew-normal.R says the iteration then stops "still
  # accurate to well within the noise of the distribution function". Where
  # psnorm() is still resolved, the error is small: observed <= 1e-5 SDs at
  # p = 1e-14 (tolerance 3e-5). Below that, psnorm() rounds to 0 within
  # about 0.1 SDs of the true quantile, which bounds what any refinement on
  # it can do (tolerance 0.25 SDs). Observed: 1.1e-5 and 5.5e-4 SDs for
  # alpha = 2 at p = 1e-16 and 1e-18, 1e-3 SDs for alpha = 5 at p = 1e-16,
  # but -47 SDs for alpha = 5 at p = 1e-18 and -41 and -50 SDs for
  # alpha = 20 at p = 1e-16 and 1e-18 (qsnorm(1e-18, 0, 1, 5) = -30.95, true
  # -1.627). Where psnorm() returns 0, qsn_std_logtail()'s safeguard
  # "step back towards the centre" sets step = +1 for a lower tail (-1 for
  # an upper), and q - step moves further out, every iteration up to the cap.
  sn_light_quantile_ref <- function(p, a) {
    # The lower tail for a > 0; F(z) <= 2 Phi(z) Phi(a z) brackets the root
    lo <- -sqrt(2 * log(1 / p) / (1 + a^2)) - 1
    stats::uniroot(
      function(z) log(sn_cdf_ref(z, a)) - log(p),
      c(lo, 0),
      tol = 1e-12
    )$root
  }
  cases <- expand.grid(p = c(1e-14, 1e-16, 1e-18), alpha = c(2, 5, 20))
  err <- vapply(
    seq_len(nrow(cases)),
    function(i) {
      p <- cases$p[i]
      a <- cases$alpha[i]
      abs(qsnorm(p, 0, 1, a) - sn_light_quantile_ref(p, a)) / sn_sd_std(a)
    },
    numeric(1)
  )
  names(err) <- sprintf("alpha=%g, p=%g", cases$alpha, cases$p)
  resolved <- cases$p >= 1e-14
  expect_all_below(err[resolved], 3e-5)
  expect_all_below(err[!resolved], 0.25)
})

# sn_score_map ------------------------------------------------------------------

test_that("sn_score_map maps normal scores to skew-normal quantiles within 1e-6 SDs", {
  # sampling() pushes each normal score z through this map, so it must give
  # Q(Phi(z)). R/skew-normal.R documents an error within 1e-6 posterior SDs
  # (5e-7 at |alpha| = 20) for |z| <= 7. Checked midway between the spline
  # knots, where the interpolation error peaks, by quadrature on the smaller
  # tail. Observed max 4.4e-7 SDs (alpha = -20); 1.8e-15 for alpha = 0.
  zk <- seq(-7, 7, length.out = 141)
  z_mid <- ((zk[-1] + zk[-141]) / 2)[seq(1, 140, by = 2)]
  xi <- 0.3
  omega <- 1.7
  err <- NULL
  for (a in c(-20, -5, -1.5, -0.3, 0, 0.3, 1, 2, 20)) {
    q <- (sn_score_map(xi, omega, a)(z_mid) - xi) / omega
    e <- sn_quantile_err(
      q,
      alpha = a,
      upper = z_mid > 0,
      tail = stats::pnorm(-abs(z_mid))
    ) /
      sn_sd_std(a)
    names(e) <- sprintf("alpha=%g, z=%g", a, z_mid)
    err <- c(err, e)
  }
  expect_all_below(err, 1e-6)
})

test_that("sn_score_map is monotone and continues linearly beyond its knots", {
  # Beyond |z| = 7 (tail probability 1.3e-12) the map is linear. It has no
  # documented accuracy there; observed error up to 2.5e-3 SDs at |z| = 8,
  # and these draws occur with probability 1e-15.
  for (a in c(-20, -1.5, 0, 2, 20)) {
    map <- sn_score_map(0, 1, a)
    z <- seq(-9, 9, by = 0.01)
    expect_true(all(diff(map(z)) > 0))
    # continuous at the ends of the knot range
    expect_lt(abs(map(7 + 1e-9) - map(7 - 1e-9)), 1e-8)
    expect_lt(abs(map(-7 + 1e-9) - map(-7 - 1e-9)), 1e-8)
    z_out <- c(-8, -7.5, 7.5, 8)
    err <- sn_quantile_err(
      map(z_out),
      alpha = a,
      upper = z_out > 0,
      tail = stats::pnorm(-abs(z_out))
    ) /
      sn_sd_std(a)
    names(err) <- sprintf("alpha=%g, z=%g", a, z_out)
    expect_all_below(err, 7.5e-3)
  }
})

# sn_moments and sn_mode --------------------------------------------------------

test_that("sn_moments gives the mean and SD found by integration", {
  # Observed relative errors <= 2.2e-16.
  for (p in c(sn_pars, list(c(0, 1, 100)))) {
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
    sd_ref <- sqrt(var_ref)
    expect_lt(abs(mom$mean - mean_ref) / sd_ref, 1e-13)
    expect_lt(abs(mom$sd / sd_ref - 1), 1e-13)
    expect_lt(abs(sn_sd_std(p[3]) / (sd_ref / p[2]) - 1), 1e-13)
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
  # Reference 1: the root of the score, to 1e-14. Reference 2: optimize() on
  # the density, which is itself only good to about 1e-8 relative. Observed
  # errors <= 9.7e-9 and <= 1.5e-8 SDs.
  for (p in c(sn_pars, list(c(0, 1, 100), c(0, 1, -100)))) {
    mode <- sn_mode(p[1], p[2], p[3])
    sd <- sn_moments(p[1], p[2], p[3])$sd
    mode_score <- p[1] + p[2] * sn_std_mode_ref(p[3])
    mode_opt <- stats::optimize(
      function(x) sn_density_ref((x - p[1]) / p[2], p[3]),
      p[1] + c(-4, 4) * p[2],
      maximum = TRUE,
      tol = 1e-12
    )$maximum
    expect_lt(abs(mode - mode_score) / sd, 3e-8)
    expect_lt(abs(mode - mode_opt) / sd, 5e-8)
  }
})

test_that("sn_mode is xi when alpha = 0", {
  # Observed error <= 6.7e-9 omega.
  for (p in list(c(0, 1), c(3, 0.2), c(-5, 4))) {
    expect_lt(abs(sn_mode(p[1], p[2], 0) - p[1]) / p[2], 2e-8)
  }
})

test_that("sn_mode is accurate when the location is large relative to the scale", {
  # The mode is equivariant: mode(xi, omega, alpha) = xi + omega z0(alpha).
  # sn_mode() maximises on the original scale with optimize(), whose
  # tolerance has a relative term sqrt(eps) |x| ~ 1.5e-8 |x|, so its error in
  # SDs grows like |xi| / omega. Observed: 4.3e-5 SDs at (100, 0.01),
  # 9.2e-4 at (1000, 0.001) and 0.24 SDs at (2000, 1e-4) even with
  # alpha = 0, where the mode is exactly xi. (The tolerance is the 3e-8 SDs
  # that sn_mode() reaches at unit scale, above.)
  pars <- list(
    c(100, 0.01, 3),
    c(1000, 0.001, 3),
    c(-500, 0.001, -2),
    c(2000, 1e-4, 0)
  )
  err <- vapply(
    pars,
    function(p) {
      ref <- p[1] + p[2] * sn_std_mode_ref(p[3])
      abs(sn_mode(p[1], p[2], p[3]) - ref) / (p[2] * sn_sd_std(p[3]))
    },
    numeric(1)
  )
  names(err) <- vapply(pars, paste, "", collapse = ", ")
  expect_all_below(err, 3e-8)
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
  # The input is an exact skew-normal, so the fit is exact up to the
  # optimiser's convergence. Observed parameter errors <= 1.8e-14, nmad
  # <= 8.9e-16 and log-density errors <= 2.8e-13 (on values down to -100).
  for (p in pars) {
    fit <- fit_sn_check(p)
    expect_lt(max(abs(c(fit$xi, fit$omega, fit$alpha) - p)), 1e-12)
    expect_lt(fit$nmad, 1e-14)
    # logC restores the normalisation removed by max(y) = 0
    expect_lt(
      max(abs(
        sn_logdens_ref(grid_scan, fit$xi, fit$omega, fit$alpha) +
          fit$logC -
          fit$y
      )),
      1e-11
    )
  }
})

test_that("fit_skew_normal fits mirror-image inputs with mirror-image parameters", {
  # On a grid symmetric about 0, rev(y) is the log-density of
  # SN(-xi, omega, -alpha). Observed max difference 1e-13 (alpha = +-0.05).
  for (a in c(-2, -0.2, -0.05, 0.05, 0.2, 2)) {
    p <- scan_sn_pars(a)
    y <- sn_logdens_ref(grid_scan, p[1], p[2], p[3])
    y <- y - max(y)
    fit <- fit_skew_normal(grid_scan, y)
    fit_rev <- fit_skew_normal(grid_scan, rev(y))
    expect_lt(
      max(abs(
        c(fit_rev$xi, fit_rev$omega, fit_rev$alpha) -
          c(-fit$xi, fit$omega, -fit$alpha)
      )),
      1e-12,
      label = sprintf("mirrored fit for alpha = %g", a)
    )
  }
})

test_that("fit_skew_normal recovers mildly skewed scans, centred or shifted", {
  # The volume correction adds gamma1 * t to the scan, which moves its mode
  # off t = 0, so shifted inputs are routine. Near alpha = 0 the shape is
  # weakly identified, and a start with the wrong sign of skew can stall at
  # the stationary point alpha = 0 (it did, for 9 of these 24 cases, before
  # the mirror-image start was added). Observed alpha errors <= 7e-14 and
  # nmad <= 8.9e-16.
  cases <- expand.grid(
    alpha = c(-0.5, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.5),
    mode = c(-0.5, 0, 0.5)
  )
  err <- numeric(nrow(cases))
  nmad <- numeric(nrow(cases))
  for (i in seq_len(nrow(cases))) {
    fit <- fit_sn_check(scan_sn_pars(cases$alpha[i], cases$mode[i]))
    err[i] <- abs(fit$alpha - cases$alpha[i])
    nmad[i] <- fit$nmad
  }
  labels <- sprintf("alpha=%g, mode=%g", cases$alpha, cases$mode)
  expect_all_below(stats::setNames(err, labels), 1e-12)
  expect_all_below(stats::setNames(nmad, labels), 1e-14)
})

test_that("fit_skew_normal gives alpha near 0 and the right moments for a Gaussian input", {
  # Near alpha = 0 the shape is weakly identified (the skewness is O(alpha^3)),
  # so alpha itself stops short of 0. Observed |alpha| <= 0.0064, moment
  # errors <= 1.4e-8 and nmad <= 1.3e-8.
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

test_that("fit_skew_normal handles uneven grids, as when scan points are dropped", {
  # fit_marginals() drops scan points where the log-posterior is not finite,
  # so the grid can have gaps or lose an end. Observed parameter errors
  # <= 6.7e-15.
  for (a in c(-5, -0.3, 0.3, 2)) {
    p <- scan_sn_pars(a, 0.3)
    for (drop in list(c(1, 2), c(20, 21), c(3, 8, 15), 1:4)) {
      fit <- fit_sn_check(p, grid_scan[-drop])
      expect_lt(max(abs(c(fit$xi, fit$omega, fit$alpha) - p)), 1e-12)
    }
  }
})

test_that("fit_skew_normal ignores points below threshold_log_drop", {
  # Six of the 21 points fall below -6 for this scan; pushing them further
  # down must not move the fit, and a stricter threshold still recovers it.
  p <- scan_sn_pars(5)
  y <- sn_logdens_ref(grid_scan, p[1], p[2], p[3])
  y <- y - max(y)
  low <- y < -6
  expect_gt(sum(low), 0)
  y_low <- y
  y_low[low] <- y_low[low] - 50
  fit <- fit_skew_normal(grid_scan, y_low)
  expect_lt(max(abs(c(fit$xi, fit$omega, fit$alpha) - p)), 1e-12)
  # nmad is over all points, including the dropped ones
  expect_gt(fit$nmad, 1e-4)

  fit <- fit_skew_normal(grid_scan, y, threshold_log_drop = -2)
  expect_lt(max(abs(c(fit$xi, fit$omega, fit$alpha) - p)), 1e-12)
})

test_that("fit_skew_normal keeps at least five points for a narrow scan", {
  # Only three points are above -6 here; with four parameters, the top five
  # points are kept so that the fit is still determined. Observed errors
  # <= 8.9e-16.
  for (omega in c(0.3, 0.2)) {
    p <- c(0.05, omega, 2)
    y <- sn_logdens_ref(grid_scan, p[1], p[2], p[3])
    y <- y - max(y)
    expect_lt(sum(y >= -6), 5)
    fit <- fit_skew_normal(grid_scan, y)
    expect_lt(max(abs(c(fit$xi, fit$omega, fit$alpha) - p)), 1e-12)
  }
})

test_that("fit_skew_normal returns the least-squares fit and its nmad for a non-skew-normal input", {
  # The scan of log(lambda) for lambda ~ Gamma(k, 1), centred at the mode
  # with unit curvature: y(t) = sqrt(k) t - k (exp(t / sqrt(k)) - 1). The
  # fit should minimise the weighted least-squares objective documented in
  # R/skew-normal.R (weights exp(y), points below the threshold dropped but
  # at least five kept), written out here. At the fit, a Newton step on that
  # objective, with central-difference gradient and Hessian, should be zero
  # up to the differencing noise. Observed: Newton steps <= 7.1e-9 (a 1e-4
  # shift in alpha gives a step of 1e-4), the Hessian is positive definite,
  # and nmad matches its definition to 5.1e-14 relative.
  ref_objective <- function(par, x, y, threshold) {
    keep <- y >= threshold | rank(-y, ties.method = "first") <= 5
    x <- x[keep]
    y <- y[keep]
    w <- exp(y) / sum(exp(y))
    ld <- sn_logdens_ref(x, par[1], exp(par[2]), par[3]) + par[4]
    sum(w * (y - ld)^2)
  }
  unit <- function(i, h) replace(numeric(4), i, h)
  fd_newton_step <- function(f, par) {
    h <- 1e-5
    g <- vapply(
      1:4,
      function(i) (f(par + unit(i, h)) - f(par - unit(i, h))) / (2 * h),
      numeric(1)
    )
    h <- 1e-4
    H <- matrix(0, 4, 4)
    for (i in 1:4) {
      for (j in 1:4) {
        ei <- unit(i, h)
        ej <- unit(j, h)
        H[i, j] <- (f(par + ei + ej) -
          f(par + ei - ej) -
          f(par - ei + ej) +
          f(par - ei - ej)) /
          (4 * h^2)
      }
    }
    list(step = solve(H, g), min_eigen = min(eigen(H, symmetric = TRUE)$values))
  }
  cases <- expand.grid(k = c(3, 50), sgn = c(1, -1), threshold = c(-6, -3))
  for (i in seq_len(nrow(cases))) {
    k <- cases$k[i]
    sgn <- cases$sgn[i]
    threshold <- cases$threshold[i]
    t <- sgn * grid_scan
    y <- sqrt(k) * t - k * (exp(t / sqrt(k)) - 1)
    y <- y - max(y)
    fit <- fit_skew_normal(grid_scan, y, threshold_log_drop = threshold)
    newton <- fd_newton_step(
      function(par) ref_objective(par, grid_scan, y, threshold),
      c(fit$xi, log(fit$omega), fit$alpha, fit$logC)
    )
    label <- sprintf("k = %g, sign %g, threshold %g", k, sgn, threshold)
    expect_lt(max(abs(newton$step)), 2e-8, label = label)
    expect_gt(newton$min_eigen, 0, label = label)
    # log(lambda) is skewed to the left
    expect_lt(fit$alpha * sgn, -0.9, label = label)

    nmad_ref <- max(abs(
      exp(y) -
        exp(
          sn_logdens_ref(grid_scan, fit$xi, fit$omega, fit$alpha) +
            fit$logC
        )
    )) /
      max(exp(y))
    expect_equal(fit$nmad, nmad_ref, tolerance = 1.5e-13, label = label)
    expect_gt(fit$nmad, 1e-3)
  }
})
