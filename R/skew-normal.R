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
# T(h, a) = Phi(h)/2 + Phi(ah)/2 - Phi(h) Phi(ah) - T(ah, 1/a). On the
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
  ph <- stats::pnorm(h)
  pah <- stats::pnorm(ah)
  sgn * ifelse(big, 0.5 * ph + 0.5 * pah - ph * pah - t_red, t_red)
}

# Skew-normal distribution function F(q) = Phi(z) - 2 T(z, alpha), with
# z = (q - xi) / omega. The upper tail is computed by reflection,
# 1 - F(z; alpha) = F(-z; -alpha), so a small tail probability never comes
# out as the difference of two numbers close to one.
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
  ok <- is.finite(z) & is.finite(alpha) & omega > 0
  out <- rep(NA_real_, n)
  out[ok] <- stats::pnorm(z[ok]) - 2 * owen_t(z[ok], alpha[ok])
  out[ok] <- pmin(pmax(out[ok], 0), 1)
  out
}

# Skew-normal quantiles to near machine precision: the fast approximation of
# qsnorm_fast() (relative error below about 1e-7) polished by two Newton steps
# on psnorm(). Used for the reported quantiles; the copula draws use
# qsnorm_fast() directly, where its accuracy is ample.
qsnorm <- function(p, xi = 0, omega = 1, alpha = 0) {
  q <- qsnorm_fast(p, xi = xi, omega = omega, alpha = alpha)
  for (i in 1:2) {
    f <- dsnorm(q, xi = xi, omega = omega, alpha = alpha)
    step <- (psnorm(q, xi = xi, omega = omega, alpha = alpha) - p) / f
    ok <- is.finite(step) & f > 0
    q[ok] <- q[ok] - step[ok]
  }
  q
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
  # Trapezoidal integration weights, for moment-based starting values
  n <- length(x)
  w_int <- numeric(n)
  w_int[1] <- 0.5 * (x[2] - x[1])
  w_int[-1] <- 0.5 * (c(diff(x[-1]), 0) + c(0, diff(x[-n])))

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

  # Fitting weights do not depend on the parameters, so compute them once
  w <- exp(temp * y)
  w[y < threshold_log_drop] <- 0
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

  fit <- stats::nlminb(
    c(xi_init, log(1 / omega_init), alpha_init, 0),
    ob,
    gr,
    hs
  )
  xi_hat <- fit$par[1]
  omega_hat <- exp(-fit$par[2])
  alpha_hat <- fit$par[3]
  logC_hat <- fit$par[4]

  # Normalised max abs deviation (NMAD) of the fitted density
  expy <- exp(y)
  fity <- dsnorm(x, xi_hat, omega_hat, alpha_hat, logC_hat)

  list(
    xi = xi_hat,
    omega = omega_hat,
    alpha = alpha_hat,
    logC = logC_hat,
    nmad = max(abs(expy - fity)) / max(expy)
  )
}
