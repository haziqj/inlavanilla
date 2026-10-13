# Skew-normal utilities, ported from INLAvaan. Everything here is internal: the
# fitted marginals are summarised by the package itself, and users who want
# densities or quantiles on another scale go through sampling().

# Log or plain density of SN(xi, omega, alpha), with an optional log-constant.
dsnorm <- function(x, xi, omega, alpha, logC = 0, log = FALSE) {
  xx <- (x - xi) / omega
  if (isTRUE(log)) {
    log(2) +
      stats::dnorm(xx, log = TRUE) +
      stats::pnorm(alpha * xx, log.p = TRUE) -
      log(omega) +
      logC
  } else {
    2 / omega * stats::dnorm(xx) * stats::pnorm(alpha * xx) * exp(logC)
  }
}

# Gauss-Legendre nodes and weights on [0, 1] (Golub-Welsch). Evaluated once
# when the namespace is built, so Owen's T below costs one small matrix
# product per call.
gauss_legendre_01 <- function(n) {
  i <- seq_len(n - 1)
  b <- i / sqrt(4 * i^2 - 1)
  J <- matrix(0, n, n)
  J[cbind(i, i + 1)] <- b
  J[cbind(i + 1, i)] <- b
  e <- eigen(J, symmetric = TRUE)
  list(x = (e$values + 1) / 2, w = e$vectors[1, ]^2)
}
gl_owen <- gauss_legendre_01(48L)

# Owen's T(h, a) = (2 pi)^-1 int_0^a exp(-h^2 (1 + x^2) / 2) / (1 + x^2) dx,
# vectorised over h and a. It is reduced to h >= 0 and 0 <= a <= 1 by the
# symmetries T(-h, a) = T(h, a), T(h, -a) = -T(h, a) and, for a > 1,
# T(h, a) = Q(h)/2 + Q(ah)/2 - Q(h) Q(ah) - T(ah, 1/a), with Q(x) = 1 - Phi(x).
# (INLAvaan writes the last identity with Phi; it is the same algebraically,
# but as a difference of numbers close to one it loses the relative accuracy
# of small values of T, which the heavy tail of psnorm() needs.) On the
# reduced range 48 Gauss-Legendre nodes give about 1e-15 absolute accuracy.
owen_t <- function(h, a) {
  n <- max(length(h), length(a))
  h <- abs(rep_len(h, n))
  a <- rep_len(a, n)
  sgn <- sign(a)
  a <- abs(a)
  big <- a > 1
  ah <- a * h
  aa <- ifelse(big, 1 / a, a)
  hh <- ifelse(big, ah, h)
  X <- outer(aa, gl_owen$x)
  G <- exp(-0.5 * hh^2 * (1 + X^2)) / (1 + X^2)
  t_red <- as.numeric(G %*% gl_owen$w) * aa / (2 * pi)
  qh <- stats::pnorm(h, lower.tail = FALSE)
  qah <- stats::pnorm(ah, lower.tail = FALSE)
  sgn * ifelse(big, 0.5 * qh + 0.5 * qah - qh * qah - t_red, t_red)
}

# Skew-normal distribution function F(q) = Phi(z) - 2 T(z, alpha), with
# z = (q - xi) / omega. The upper tail is computed by reflection,
# 1 - F(z; alpha) = F(-z; -alpha), so a small tail probability never comes
# out as the difference of two numbers close to one. It is 0 and 1 at -Inf
# and Inf, and NA for a non-positive or missing scale.
psnorm <- function(q, xi = 0, omega = 1, alpha = 0, lower_tail = TRUE) {
  n <- max(length(q), length(xi), length(omega), length(alpha))
  q <- rep_len(q, n)
  xi <- rep_len(xi, n)
  omega <- rep_len(omega, n)
  alpha <- rep_len(alpha, n)
  z <- (q - xi) / omega
  if (!isTRUE(lower_tail)) {
    z <- -z
    alpha <- -alpha
  }
  valid <- !is.na(z) & is.finite(alpha) & !is.na(omega) & omega > 0
  ok <- valid & is.finite(z)
  out <- rep(NA_real_, n)
  out[ok] <- stats::pnorm(z[ok]) - 2 * owen_t(z[ok], alpha[ok])
  out[ok] <- pmin(pmax(out[ok], 0), 1)
  out[valid & z == -Inf] <- 0
  out[valid & z == Inf] <- 1
  out
}

# Skew-normal quantiles, solved to near machine precision. Each probability
# is turned into the log of its smaller tail, and qsn_std_logtail() solves for
# the standardised quantile.
qsnorm <- function(p, xi = 0, omega = 1, alpha = 0) {
  n <- max(length(p), length(xi), length(omega), length(alpha))
  p <- rep_len(p, n)
  xi <- rep_len(xi, n)
  omega <- rep_len(omega, n)
  alpha <- rep_len(alpha, n)
  out <- rep(NA_real_, n)
  ok <- !is.na(p) & p >= 0 & p <= 1 & is.finite(alpha) & omega > 0
  out[ok & p == 0] <- -Inf
  out[ok & p == 1] <- Inf
  mid <- which(ok & p > 0 & p < 1)
  if (length(mid) > 0L) {
    upper <- p[mid] > 0.5
    logp <- ifelse(upper, log1p(-p[mid]), log(p[mid]))
    z <- qsn_std_logtail(logp, alpha[mid], upper)
    out[mid] <- xi[mid] + omega[mid] * z
  }
  out
}

# Standardised skew-normal quantile from the log of a tail probability: solves
# log F(q) = logp, or log(1 - F(q)) = logp where `upper` is TRUE, for
# SN(0, 1, alpha). Newton's method on the log scale, from qsnorm_fast(). The
# skew-normal is log-concave, so log F and log(1 - F) are concave and
# monotone, and Newton's method converges monotonically from any start; steps
# are capped at one unit as a safeguard. Working on the log of the smaller
# tail keeps far-tail quantiles accurate. Deep in the light tail (tail
# probabilities around 1e-12 and below), psnorm() loses relative accuracy to
# cancellation and the iteration stops at the 30-step cap, still accurate to
# well within the noise of the distribution function there.
qsn_std_logtail <- function(logp, alpha, upper) {
  n <- length(logp)
  alpha <- rep_len(alpha, n)
  upper <- rep_len(upper, n)
  p <- exp(logp)
  q <- .qsn_std_fast(ifelse(upper, 1 - p, p), alpha)
  q[!is.finite(q)] <- ifelse(upper, 8, -8)[!is.finite(q)]
  active <- seq_len(n)
  for (it in 1:30) {
    qi <- q[active]
    up <- upper[active]
    tail <- numeric(length(active))
    if (any(!up)) {
      tail[!up] <- psnorm(qi[!up], 0, 1, alpha[active][!up])
    }
    if (any(up)) {
      tail[up] <- psnorm(qi[up], 0, 1, alpha[active][up], lower_tail = FALSE)
    }
    f <- dsnorm(qi, 0, 1, alpha[active])
    # d/dq log F = f / F and d/dq log(1 - F) = -f / (1 - F)
    step <- (log(tail) - logp[active]) * tail / ifelse(up, -f, f)
    # A tail that underflows means q is far out: step back towards the centre
    bad <- !is.finite(step)
    step[bad] <- ifelse(up[bad], -1, 1)
    step <- pmax(pmin(step, 1), -1)
    q[active] <- qi - step
    active <- active[abs(step) > 1e-10 * (1 + abs(qi))]
    if (length(active) == 0L) {
      break
    }
  }
  q
}

# The map z -> Q(Phi(z)) from normal scores to the quantiles of
# SN(xi, omega, alpha), for the copula draws. It is smooth and close to
# linear, so a cubic spline through knots solved by qsn_std_logtail() at
# z = -7, ..., 7 reproduces it to within 1e-6 posterior SDs (5e-7 at
# |alpha| = 20, far less for mild skewness), at the cost of one spline
# evaluation per draw. Beyond the knots, a probability of 3e-12, it continues
# linearly.
sn_score_map <- function(xi, omega, alpha, zmax = 7, n_knots = 141L) {
  zk <- seq(-zmax, zmax, length.out = n_knots)
  logp <- stats::pnorm(-abs(zk), log.p = TRUE)
  qk <- xi + omega * qsn_std_logtail(logp, alpha, upper = zk > 0)
  sp <- stats::splinefun(zk, qk, method = "fmm")
  slope <- sp(c(-zmax, zmax), deriv = 1)
  function(z) {
    out <- sp(z)
    lo <- z < -zmax
    hi <- z > zmax
    out[lo] <- qk[1] + slope[1] * (z[lo] + zmax)
    out[hi] <- qk[n_knots] + slope[2] * (z[hi] - zmax)
    out
  }
}

# Mean and standard deviation of SN(xi, omega, alpha), in closed form.
sn_moments <- function(xi, omega, alpha) {
  delta <- alpha / sqrt(1 + alpha^2)
  list(
    mean = xi + omega * delta * sqrt(2 / pi),
    sd = omega * sqrt(1 - 2 * delta^2 / pi)
  )
}

# Mode of SN(xi, omega, alpha). It has no closed form, but it lies between xi
# and the mean, so a bracket one SD wider on each side always contains it.
sn_mode <- function(xi, omega, alpha) {
  mom <- sn_moments(xi, omega, alpha)
  lo <- min(xi, mom$mean) - mom$sd
  hi <- max(xi, mom$mean) + mom$sd
  stats::optimize(
    function(x) dsnorm(x, xi, omega, alpha, log = TRUE),
    interval = c(lo, hi),
    maximum = TRUE,
    tol = 1e-8 * (hi - lo)
  )$maximum
}

mills_ratio <- function(u) {
  logphi <- stats::dnorm(u, log = TRUE)
  logPhi <- stats::pnorm(u, log.p = TRUE)
  # guard against -Inf in logPhi; cap to avoid exp overflow
  logPhi <- pmax(logPhi, -745) # ~exp(-745) ~ 5e-324
  exp(logphi - logPhi)
}

# Fit a skew-normal to log-density evaluations y = log f(x) at points x, by
# weighted least squares on the log scale. The weights are the density values
# raised to the temperature `temp`, so the fit concentrates where the mass is;
# points more than `threshold_log_drop` below the maximum are ignored. This is
# the R-INLA approach (gmrflib/fit-sn.c), as ported in INLAvaan, with the
# temperature fixed rather than estimated. `y` should be normalised so that
# max(y) = 0.
fit_skew_normal <- function(x, y, threshold_log_drop = -6, temp = 1) {
  # Trapezoidal integration weights, for moment-based starting values. (The
  # INLAvaan original gives the second point half its weight, which skews the
  # starting values of a symmetric scan.)
  n <- length(x)
  dx <- diff(x)
  w_int <- 0.5 * (c(dx, 0) + c(0, dx))

  p <- exp(y)
  m0 <- sum(w_int * p)
  m1 <- sum(w_int * p * x)
  m2 <- sum(w_int * p * x^2)
  if (m0 < .Machine$double.eps) {
    m0 <- 1 # nocov
  }
  mean_hat <- m1 / m0
  var_hat <- m2 / m0 - mean_hat^2
  sd_hat <- sqrt(max(.Machine$double.eps, var_hat))
  m3 <- sum(w_int * p * (x - mean_hat)^3)

  # Theoretical limit of the skew-normal skewness is 0.9952717
  skew_hat <- (m3 / m0) / (sd_hat^3)
  skew_hat <- max(min(skew_hat, 0.9), -0.9)

  # Skewness as a function of delta = alpha / sqrt(1 + alpha^2)
  f_skew_to_delta <- function(delta, target_skew) {
    num <- (4 - pi) / 2 * (delta * sqrt(2 / pi))^3
    den <- (1 - 2 * delta^2 / pi)^(3 / 2)
    num / den - target_skew
  }
  root <- tryCatch(
    stats::uniroot(
      f_skew_to_delta,
      interval = c(-0.995, 0.995),
      target_skew = skew_hat
    )$root,
    error = function(e) 0 # nocov
  )
  alpha_init <- root / sqrt(1 - root^2)
  omega_init <- sd_hat / sqrt(1 - 2 * root^2 / pi)
  xi_init <- mean_hat - omega_init * root * sqrt(2 / pi)

  # Fitting weights do not depend on the parameters, so compute them once.
  # Points below the threshold have zero weight and are dropped from the fit.
  keep <- y >= threshold_log_drop | rank(-y, ties.method = "first") <= 5
  x_all <- x
  y_all <- y
  x <- x[keep]
  y <- y[keep]
  w <- exp(temp * y)
  w <- w / sum(w)

  # Parameters: location, log inverse scale, shape, log-constant
  logdens <- function(param) {
    xx <- (x - param[1]) * exp(param[2])
    log(2) +
      stats::dnorm(xx, log = TRUE) +
      stats::pnorm(param[3] * xx, log.p = TRUE) +
      param[2] +
      param[4]
  }

  ob <- function(param) {
    sum(w * (y - logdens(param))^2)
  }

  gr <- function(param) {
    a <- param[3]
    s <- exp(param[2])
    xx <- s * (x - param[1])
    R <- mills_ratio(a * xx)
    r <- y - logdens(param)
    c(
      -2 * sum(w * r * (s * xx - a * s * R)),
      -2 * sum(w * r * (-xx^2 + a * xx * R + 1)),
      -2 * sum(w * r * (xx * R)),
      -2 * sum(w * r)
    )
  }

  hs <- function(param) {
    a <- param[3]
    s <- exp(param[2])
    xx <- s * (x - param[1])
    u <- a * xx
    R <- mills_ratio(u)
    Rprime <- -R * (u + R) # derivative of the Mills ratio
    r <- y - logdens(param)

    # First derivatives of the log-density
    L_mu <- s * xx - a * s * R
    L_ls <- -xx^2 + a * xx * R + 1
    L_a <- xx * R

    # Second derivatives (those involving logC vanish)
    L_mumu <- -s^2 + a^2 * s^2 * Rprime
    L_muls <- 2 * s * xx - a * s * R - a^2 * s * xx * Rprime
    L_mua <- -s * R - a * s * xx * Rprime
    L_lsls <- -2 * xx^2 + a * xx * R + a^2 * xx^2 * Rprime
    L_lsa <- xx * R + a * xx^2 * Rprime
    L_aa <- xx^2 * Rprime

    # Weighted Hessian: 2 * sum w * (grad L grad L^T - r * Hess L)
    H11 <- 2 * sum(w * (L_mu * L_mu - r * L_mumu))
    H12 <- 2 * sum(w * (L_mu * L_ls - r * L_muls))
    H13 <- 2 * sum(w * (L_mu * L_a - r * L_mua))
    H14 <- 2 * sum(w * L_mu)
    H22 <- 2 * sum(w * (L_ls * L_ls - r * L_lsls))
    H23 <- 2 * sum(w * (L_ls * L_a - r * L_lsa))
    H24 <- 2 * sum(w * L_ls)
    H33 <- 2 * sum(w * (L_a * L_a - r * L_aa))
    H34 <- 2 * sum(w * L_a)
    H44 <- 2 * sum(w)

    # fmt: skip
    matrix(c(H11, H12, H13, H14,
             H12, H22, H23, H24,
             H13, H23, H33, H34,
             H14, H24, H34, H44), nrow = 4, byrow = TRUE)
  }

  # Fit from the moment-based start and from its mirror image (the opposite
  # skew), and keep the better fit. The shape is weakly identified near zero,
  # where a skew-normal differs from a shifted Gaussian only at third order,
  # so a start with the wrong sign of skew can stall at alpha = 0 when the
  # skewness of a truncated scan is mild.
  d0 <- if (abs(root) < 0.2) 0.5 else root
  starts <- lapply(c(d0, -d0), function(d) {
    om <- sd_hat / sqrt(1 - 2 * d^2 / pi)
    c(mean_hat - om * d * sqrt(2 / pi), log(1 / om), d / sqrt(1 - d^2), 0)
  })
  if (abs(root) >= 0.2) {
    starts[[1]] <- c(xi_init, log(1 / omega_init), alpha_init, 0)
  }
  fits <- lapply(starts, function(st) stats::nlminb(st, ob, gr, hs))
  fit <- fits[[which.min(vapply(fits, `[[`, 0, "objective"))]]
  xi_hat <- fit$par[1]
  omega_hat <- exp(-fit$par[2])
  alpha_hat <- fit$par[3]
  logC_hat <- fit$par[4]

  # Normalised max abs deviation (NMAD) of the fitted density, over all points
  expy <- exp(y_all)
  fity <- dsnorm(x_all, xi_hat, omega_hat, alpha_hat, logC_hat)

  list(
    xi = xi_hat,
    omega = omega_hat,
    alpha = alpha_hat,
    logC = logC_hat,
    nmad = max(abs(expy - fity)) / max(expy)
  )
}
