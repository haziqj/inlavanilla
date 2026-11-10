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
sampling(object, nsamp = 1000L, fun = NULL, ..., transform = TRUE)
```

## Arguments

- object:

  An object to draw from.

- ...:

  For the `inla_vanilla` method, further arguments passed to `fun`. As
  in [`stats::optim()`](https://rdrr.io/r/stats/optim.html), a name that
  abbreviates `nsamp` or `fun` would be taken as that argument: a data
  argument named `n` would set the number of draws. Give such data a
  longer name, or wrap it in `fun`.

- nsamp:

  Number of draws.

- fun:

  Optional function applied to each draw. It receives one draw as a
  named numeric vector, like `logpost`, and returns a numeric vector of
  a fixed length (named, ideally). When `NULL`, the parameters
  themselves are returned.

- transform:

  Logical. Return the transformed parameters on their natural scale?
  Default `TRUE`. `FALSE` returns every parameter on the scale of
  `logpost`. It must be named in full, and has no effect on what `fun`
  receives.

## Value

A numeric matrix with one row per draw: the parameters, with the
transformed ones on their natural scale, or the values of `fun`.

## Details

For a fit with `marginals = "gaussian"` this is a draw from the joint
Gaussian \\N(\tilde\vartheta, \Omega)\\, centred at the VB mean (or at
the mode when `vb = "off"`).

Draws never call the log-posterior, so they are cheap: large `nsamp`
costs little. A parameter named for its scale (`log_sigma`, `logit_p`,
`atanh_rho`) is returned on its natural scale, in a column named
`sigma`, `p` or `rho` (see the section on natural-scale summaries in
[`inla_vanilla()`](http://haziqj.ml/inlavanilla/reference/inla_vanilla.md)).
Use `fun` for quantities of several parameters, such as a ratio or a
covariance matrix. `fun` is written as `logpost` is: it receives each
draw on the scale of `logpost`, under the names in `init`, whatever
`transform`.

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
fit <- inla_vanilla(logpost, c(mu = 0, log_sigma = 0), y = y, silent = TRUE)

# Draws of mu and sigma: the parameter named log_sigma comes back as sigma
draws <- sampling(fit, nsamp = 2000)
head(draws)
#>            mu    sigma
#> [1,] 5.603356 1.677596
#> [2,] 5.125348 1.732507
#> [3,] 5.285843 1.760798
#> [4,] 5.141380 1.826636
#> [5,] 4.707959 2.303939
#> [6,] 5.023140 2.180450

# The same draws as logpost sees them
head(sampling(fit, nsamp = 2000, transform = FALSE))
#>            mu log_sigma
#> [1,] 4.621822 0.5056868
#> [2,] 5.253538 0.3865430
#> [3,] 5.908196 0.7070464
#> [4,] 5.444282 0.5663662
#> [5,] 5.037834 0.4417673
#> [6,] 5.583696 0.5130578

# A derived quantity: the coefficient of variation sigma / mu
cv <- sampling(fit, nsamp = 2000, fun = function(theta) {
  c(cv = exp(theta[["log_sigma"]]) / theta[["mu"]])
})
quantile(cv, c(0.025, 0.5, 0.975))
#>      2.5%       50%     97.5% 
#> 0.2785661 0.3576390 0.4836230 
```
