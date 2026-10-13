# Finite differences. Every step `h` is per coordinate and recycled to the
# length of `x`, so callers can scale the steps by the posterior SDs once those
# are known.

# Central-difference gradient of a scalar function: 2m calls of `fn`.
fd_grad <- function(fn, x, h) {
  m <- length(x)
  h <- rep_len(h, m)
  g <- numeric(m)
  for (j in seq_len(m)) {
    xp <- xm <- x
    xp[j] <- x[j] + h[j]
    xm[j] <- x[j] - h[j]
    g[j] <- (fn(xp) - fn(xm)) / (2 * h[j])
  }
  g
}

# One-sided (forward) difference gradient of a scalar function: m + 1 calls
# of `fn`, or m when f(x) is supplied. Its error is O(h), half a step times
# the curvature, so callers should correct for it (see inla_vanilla()).
fd_grad_forward <- function(fn, x, h, f0 = fn(x)) {
  m <- length(x)
  h <- rep_len(h, m)
  g <- numeric(m)
  for (j in seq_len(m)) {
    xp <- x
    xp[j] <- x[j] + h[j]
    g[j] <- (fn(xp) - f0) / h[j]
  }
  g
}

# Central-difference Jacobian of a gradient function: 2m calls of `gr`.
fd_jacobian <- function(gr, x, h) {
  m <- length(x)
  h <- rep_len(h, m)
  H <- matrix(0, m, m)
  for (j in seq_len(m)) {
    xp <- xm <- x
    xp[j] <- x[j] + h[j]
    xm[j] <- x[j] - h[j]
    H[, j] <- (gr(xp) - gr(xm)) / (2 * h[j])
  }
  H
}

# Hessian of a scalar function from function values alone. The diagonal uses
# the usual three-point second difference. Each off-diagonal entry reuses the
# axial points and needs only two new calls,
#   f(x + hi ei + hj ej) + f(x - hi ei - hj ej) - f(x +- hi ei) - f(x +- hj ej)
#     + 2 f(x) = 2 hi hj H_ij + O(h^4),
# so the error is O(h^2), as for the four-point formula, at m^2 + m + 1 calls
# of `fn` rather than 2m^2 + 1.
fd_hessian <- function(fn, x, h, f0 = fn(x)) {
  m <- length(x)
  h <- rep_len(h, m)
  fp <- fm <- numeric(m)
  H <- matrix(0, m, m)
  for (i in seq_len(m)) {
    xp <- xm <- x
    xp[i] <- x[i] + h[i]
    xm[i] <- x[i] - h[i]
    fp[i] <- fn(xp)
    fm[i] <- fn(xm)
    H[i, i] <- (fp[i] - 2 * f0 + fm[i]) / h[i]^2
  }
  if (m > 1) {
    for (i in seq_len(m - 1)) {
      for (j in (i + 1):m) {
        xpp <- xmm <- x
        xpp[c(i, j)] <- x[c(i, j)] + h[c(i, j)]
        xmm[c(i, j)] <- x[c(i, j)] - h[c(i, j)]
        H[i, j] <- H[j, i] <- (fn(xpp) +
          fn(xmm) -
          fp[i] -
          fm[i] -
          fp[j] -
          fm[j] +
          2 * f0) /
          (2 * h[i] * h[j])
      }
    }
  }
  H
}

# Second derivative of a scalar function along direction `v` at `x`, from
# three calls (two when f(x) is supplied).
fd_curvature <- function(fn, x, v, h, f0 = fn(x)) {
  (fn(x + h * v) - 2 * f0 + fn(x - h * v)) / h^2
}

# Format a number of seconds for messages.
format_time <- function(s) {
  if (s < 1) {
    paste0(round(1000 * s), "ms")
  } else if (s < 60) {
    paste0(format(round(s, 1), nsmall = 1), "s")
  } else {
    paste0(floor(s / 60), "m ", round(s %% 60), "s")
  }
}

elapsed <- function() {
  proc.time()[["elapsed"]]
}
