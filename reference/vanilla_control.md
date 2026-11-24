# Tuning of an INLAvanilla fit

Builds the `control` list of
[`inla_vanilla()`](https://haziqj.ml/inlavanilla/reference/inla_vanilla.md):
the limits of the mode search, of the VB correction and of the
skew-normal scans. The defaults suit most fits. A plain list with some
of these names works too, and is checked by the same rules:
`control = list(n_grid = 31)` is
`control = vanilla_control(n_grid = 31)`.

## Usage

``` r
vanilla_control(
  mode_maxit = 5000L,
  mode_maxeval = 5000L,
  mode_tol = 1e-10,
  vb_maxit = 25L,
  vb_tol = 0.001,
  n_grid = 21L,
  log_threshold = -6,
  nlminb = list()
)
```

## Arguments

- mode_maxit, mode_maxeval:

  The most iterations, and evaluations of `logpost`, of the mode search
  by [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html) (its
  `iter.max` and `eval.max`). When the search stops short of the mode,
  it is restarted once, with each parameter scaled by its curvature and
  the same limits.

- mode_tol:

  The relative tolerance of the mode search (the `rel.tol` of
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html)).

- vb_maxit, vb_tol:

  The iteration limit and the convergence tolerance, in posterior SDs,
  of the fixed-point iteration of the VB correction, with
  `vb_correction = "full"` and when a one-step correction has to
  iterate. A step over 1.5 posterior SDs is cut to 1.5. With
  `vb_correction = "full"`, if two steps in a row are cut, or it has not
  converged after `vb_maxit` steps,
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html) takes over.

- n_grid:

  The number of points of each skew-normal scan, spread over \\\pm 4\\
  posterior SDs. At least 7.

- log_threshold:

  Scan points whose log-posterior falls more than this below the maximum
  are left out of the skew-normal fit.

- nlminb:

  A named list of any other controls of
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html) for the mode
  search, such as `trace = 1`, which prints its iterations, or
  `abs.tol`. Its `iter.max`, `eval.max` and `rel.tol` stand in for
  `mode_maxit`, `mode_maxeval` and `mode_tol`, which must then be left
  at their defaults. It does not reach the fallback optimiser of the VB
  correction or the conditional optimisations of
  `marginal_scan = "laplace"`, which keep limits of their own.

## Value

A list with an entry for each argument, for the `control` argument of
[`inla_vanilla()`](https://haziqj.ml/inlavanilla/reference/inla_vanilla.md).

## See also

[`inla_vanilla()`](https://haziqj.ml/inlavanilla/reference/inla_vanilla.md).

## Examples

``` r
set.seed(1)
y <- rnorm(30, mean = 5, sd = 2)
logpost <- function(theta) {
  sigma <- exp(theta[["log_sigma"]])
  sum(dnorm(y, theta[["mu"]], sigma, log = TRUE)) +
    dnorm(theta[["mu"]], 0, 10, log = TRUE) +
    dexp(sigma, log = TRUE) + theta[["log_sigma"]]
}

# A finer scan of each marginal
fit <- inla_vanilla(
  logpost,
  c(mu = 0, log_sigma = 0),
  control = vanilla_control(n_grid = 41),
  silent = TRUE
)

# The same as a plain list
fit <- inla_vanilla(
  logpost,
  c(mu = 0, log_sigma = 0),
  control = list(n_grid = 41),
  silent = TRUE
)

# Print the iterations of the mode search
vanilla_control(nlminb = list(trace = 1))
#> $mode_maxit
#> [1] 5000
#> 
#> $mode_maxeval
#> [1] 5000
#> 
#> $mode_tol
#> [1] 1e-10
#> 
#> $vb_maxit
#> [1] 25
#> 
#> $vb_tol
#> [1] 0.001
#> 
#> $n_grid
#> [1] 21
#> 
#> $log_threshold
#> [1] -6
#> 
#> $nlminb
#> $nlminb$trace
#> [1] 1
#> 
#> 
```
