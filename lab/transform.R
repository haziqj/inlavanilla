library(inlavanilla)
n <- 100
x <- runif(n, 0, 5)
mu <- 0.5 + 2.5 * x
y <- rnorm(n, mu, sd = 1)

logpost <- function(theta, y, x) {
  b0 <- theta["b0"]
  b1 <- theta["b1"]
  log_tau <- theta["log_tau"]
  sd_tau <- sqrt(1 / exp(log_tau))
  mu <- b0 + b1 * x

  log_lik <- sum(dnorm(y, mu, sd = sd_tau, log = TRUE))
  log_pri <- dnorm(b0, 0, 5, log = TRUE) +
    dnorm(b1, 0, 5, log = TRUE) +
    dgamma(sd_tau, 1, 2, log = TRUE) +
    # Jacobian of sd = exp(-log_tau / 2): |d sd / d log_tau| = sd / 2
    log(sd_tau / 2)

  log_lik + log_pri
}

theta0 <- c("b0" = 1, "b1" = 1, "log_tau" = 0)
fit <- inla_vanilla(logpost, theta0, y = y, x = x)
print(fit)