# Draw from an approximate posterior

Draws from the approximate joint posterior by a Gaussian copula: the
skew-normal marginals of the fit are coupled through the correlation
matrix \\R\\ of the Laplace covariance. Each draw takes \\z \sim N(0,
R)\\, maps it to uniforms \\u_j = \Phi(z_j)\\ and then to each
parameter's skew-normal quantile, \\\vartheta_j = Q_j(u_j)\\.

## Usage

``` r
sampling(object, ...)

# S3 method for class 'inla_vanilla'
sampling(object, nsamp = 1000L, fun = NULL, ...)
```

## Arguments

- object:

  An object to draw from.

- ...:

  For the `inla_vanilla` method, further arguments passed to `fun`.

- nsamp:

  Number of draws.

- fun:

  Optional function applied to each draw. It receives one draw as a
  named numeric vector, like `logpost`, and returns a numeric vector of
  a fixed length (named, ideally). When `NULL`, the parameters
  themselves are returned.

## Value

A numeric matrix with one row per draw: the parameters on the
unconstrained scale, or the values of `fun`.

## Details

Draws never call the log-posterior, so they are cheap: large `nsamp`
costs little. Use `fun` for natural-scale or derived quantities, such as
`exp(log_sigma)` or a covariance matrix built from several parameters.

## Examples

``` r
set.seed(1)
y <- rnorm(30, mean = 5, sd = 2)
logpost <- function(theta, y) {
  sigma <- exp(theta[["log_sigma"]])
  sum(dnorm(y, theta[["mu"]], sigma, log = TRUE)) +
    dnorm(theta[["mu"]], 0, 10, log = TRUE) +
    dexp(sigma, log = TRUE) + theta[["log_sigma"]]
}
fit <- inla_vanilla(logpost, c(mu = 0, log_sigma = 0), y = y, verbose = FALSE)

draws <- sampling(fit, nsamp = 2000)
head(draws)
#>            mu log_sigma
#> [1,] 5.603356 0.5173617
#> [2,] 5.125348 0.5495693
#> [3,] 5.285843 0.5657669
#> [4,] 5.141380 0.6024758
#> [5,] 4.707959 0.8346202
#> [6,] 5.023140 0.7795313

# A derived quantity: the coefficient of variation sigma / mu
cv <- sampling(fit, nsamp = 2000, fun = function(theta) {
  c(cv = exp(theta[["log_sigma"]]) / theta[["mu"]])
})
quantile(cv, c(0.025, 0.5, 0.975))
#>      2.5%       50%     97.5% 
#> 0.2716513 0.3583291 0.4832763 
```
