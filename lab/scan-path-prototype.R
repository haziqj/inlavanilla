# Run from the package root: Rscript lab/scan-path-prototype.R
# Prototype: does correcting the scan path towards the true conditional mode
# fix the small-n under-dispersion? Compares, for the Puromycin model:
#  exact (grid), current (linear CMP + linear volume), corrected height
#  (one or more Newton steps in the other parameters, Laplace curvature) with
#  the linear volume correction, and full Laplace.
devtools::load_all(".", quiet = TRUE)
d <- subset(Puromycin, state == "treated")
logpost <- function(theta) {
  Vm <- exp(theta[[1]])
  K <- exp(theta[[2]])
  sigma <- exp(theta[[3]])
  mu <- Vm * d$conc / (K + d$conc)
  sum(dnorm(d$rate, mu, sigma, log = TRUE)) +
    dnorm(theta[[1]], log(200), 0.5, log = TRUE) +
    dnorm(theta[[2]], log(0.1), 1, log = TRUE) +
    dnorm(sigma, 0, 20, log = TRUE) +
    log(2) +
    theta[[3]]
}
init <- c(log_Vm = log(200), log_K = log(0.1), log_sigma = log(10))
fit <- inla_vanilla(logpost, init, silent = TRUE)
gr <- function(x) fd_grad(logpost, x, 1e-5)
mode <- fit$mode
S <- fit$cov
H <- solve(S)
sds <- sqrt(diag(S))
m <- 3
L <- backsolve(chol(H), diag(m))
z <- seq(-4, 4, length.out = 21)

# exact marginal sd via grid
s <- fit$summary
axes <- lapply(1:3, function(j) {
  seq(
    s[j, "Mean"] - 7 * s[j, "SD"],
    s[j, "Mean"] + 7 * s[j, "SD"],
    length.out = 70
  )
})
grid <- as.matrix(expand.grid(axes))
lpg <- apply(grid, 1, logpost)
w <- exp(lpg - max(lpg))
exact_sd <- sapply(1:3, function(j) {
  p <- tapply(w, grid[, j], sum)
  p <- p / sum(p)
  x <- axes[[j]]
  mu <- sum(x * p)
  sqrt(sum((x - mu)^2 * p))
})

sd_from_scan <- function(y) {
  y <- y - max(y)
  f <- fit_skew_normal(z, y)
  sn_moments(f$xi, f$omega, f$alpha)$sd
}
res <- matrix(
  NA,
  3,
  5,
  dimnames = list(
    names(init),
    c("exact", "current", "newton1", "newton3", "full_laplace")
  )
)
res[, "exact"] <- exact_sd
for (j in 1:3) {
  vj <- S[, j] / sds[j]
  g1 <- gamma1_lp(
    j,
    logpost,
    mode,
    logpost(mode),
    sweep(S, 2, sds, "/"),
    L,
    whitened_trace_lp(logpost, mode, logpost(mode), L)
  )
  Hoo <- H[-j, -j, drop = FALSE]
  scan <- function(newton) {
    sapply(z, function(t) {
      x <- mode + vj * t
      for (k in seq_len(newton)) {
        x[-j] <- x[-j] + solve(Hoo, gr(x)[-j])
      }
      logpost(x) + g1 * t
    })
  }
  res[j, "current"] <- sds[j] * sd_from_scan(scan(0))
  res[j, "newton1"] <- sds[j] * sd_from_scan(scan(1))
  res[j, "newton3"] <- sds[j] * sd_from_scan(scan(3))
  # full Laplace: exact conditional mode and log det of conditional Hessian
  yfull <- sapply(z, function(t) {
    th <- mode[j] + sds[j] * t
    f <- function(xo) {
      x <- mode
      x[j] <- th
      x[-j] <- xo
      logpost(x)
    }
    o <- optim(
      mode[-j] + (vj * t)[-j],
      function(xo) -f(xo),
      method = "BFGS",
      control = list(reltol = 1e-12)
    )
    Hc <- -fd_hessian(f, o$par, 1e-4)
    -o$value - 0.5 * determinant(Hc)$modulus
  })
  res[j, "full_laplace"] <- sds[j] * sd_from_scan(yfull)
}
print(round(res, 4))
cat("relative error vs exact:\n")
print(round(sweep(res[, -1], 1, res[, 1], "/") - 1, 3))
