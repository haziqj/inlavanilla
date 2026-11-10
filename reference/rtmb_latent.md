# Draw latent variables from an RTMB object

Builds the `fun` of
[`sampling()`](http://haziqj.ml/inlavanilla/reference/sampling.md) for a
model whose latent variables are the random effects of an RTMB or TMB
object. For each draw of the parameters from the fit, the function it
returns draws the random effects from their Gaussian approximation given
those parameters, so that together the draws follow the posterior of the
latent variables.

## Usage

``` r
rtmb_latent(obj)
```

## Arguments

- obj:

  An object from `RTMB::MakeADFun()` or
  [`TMB::MakeADFun()`](https://rdrr.io/pkg/TMB/man/MakeADFun.html),
  built with `random` naming the latent variables.

## Value

A function of the parameters, to pass as `fun` to
[`sampling()`](http://haziqj.ml/inlavanilla/reference/sampling.md). It
returns one draw of every random effect of `obj`, in the object's order,
named after the parameter it belongs to: `u[1]`, `u[2]` and so on for a
vector `u`.

## Details

An object built with `random` integrates its random effects out by a
Laplace approximation. Each call of `obj$fn()` finds their mode given
the parameters, and the Hessian there is the precision of a Gaussian
approximation of their distribution given the parameters. A draw is the
mode plus Gaussian noise with that precision. The approximation is exact
when the random effects are Gaussian given the parameters, as in a
linear mixed model, and close when each has plenty of data. It is
weakest for a random effect with few discrete observations: see the
article "Latent variables" on the package website.

The object must return a negative log-density, as RTMB and TMB expect,
and its parameters other than the random effects must be those of the
fit, in the same order and on the same scale. Both hold when the fit is
made from the object itself, with
`logpost = function(theta) -obj$fn(theta)` and `init = obj$par`. A fit
made from another function of the same model, such as a closed form with
the random effects integrated out by hand, can use an object built only
for this purpose.

The sparse Cholesky factor comes from the Matrix package, which RTMB and
TMB install. An object's tape does not survive
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html), so build the object
again in a new session.

## Examples

``` r
# Random intercepts u, integrated out by RTMB's Laplace approximation
set.seed(1)
group <- rep(1:8, each = 5)
y <- 2 + rnorm(8)[group] + rnorm(40, sd = 0.5)
nll <- function(p) {
  -sum(RTMB::dnorm(y, p$mu + p$u[group], exp(p$log_sigma), log = TRUE)) -
    sum(RTMB::dnorm(p$u, 0, exp(p$log_tau), log = TRUE)) -
    RTMB::dnorm(p$mu, 0, 10, log = TRUE) -
    RTMB::dnorm(p$log_tau, 0, 2, log = TRUE) -
    RTMB::dnorm(p$log_sigma, 0, 2, log = TRUE)
}
obj <- RTMB::MakeADFun(
  nll,
  list(mu = 0, log_tau = 0, log_sigma = 0, u = numeric(8)),
  random = "u",
  silent = TRUE
)
fit <- inla_vanilla(
  function(theta) -obj$fn(theta),
  obj$par,
  grad = function(theta) -as.vector(obj$gr(theta)),
  silent = TRUE
)

# One row per draw of the parameters, one column per random intercept
u <- sampling(fit, nsamp = 500, fun = rtmb_latent(obj))
round(colMeans(u), 2)
#>  u[1]  u[2]  u[3]  u[4]  u[5]  u[6]  u[7]  u[8] 
#> -0.61 -0.01 -0.65  1.06  0.29 -1.14  0.50  0.56 
```
