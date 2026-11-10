# Time taken by each stage of an INLAvanilla fit

The elapsed (wall-clock) time of each stage of
[`inla_vanilla()`](http://haziqj.ml/inlavanilla/reference/inla_vanilla.md),
in seconds, as INLAvaan's `timing()` gives them.

## Usage

``` r
timing(object, ...)

# Default S3 method
timing(object, ...)

# S3 method for class 'inla_vanilla'
timing(object, what = "total", ...)
```

## Arguments

- object:

  An object of class `inla_vanilla`.

- ...:

  Unused.

- what:

  The segments to return, by name, or `"all"` for every segment.
  Defaults to `"total"`.

## Value

A named numeric vector of class `inla_vanilla_timing`, in seconds. It
prints each time as
[`summary.inla_vanilla()`](http://haziqj.ml/inlavanilla/reference/summary.inla_vanilla.md)
does: in milliseconds below a second, then in seconds, then in minutes
and seconds.

## Details

The segments are:

- `init`: the checks of the arguments, with the calls of `logpost`,
  `grad` and `hess` at `init`.

- `laplace`: the posterior mode and the Hessian there.

- `vb`: the VB mean correction, close to zero with `vb = "off"`.

- `marginals`: the scans and the skew-normal fits, close to zero with
  `marginals = "gaussian"`.

- `diagnostics`: the diagnostics of the fit and its warning.

- `total`: the whole fit.

The segments are disjoint, and `total` is their sum. Each stage includes
its progress output, so with `verbose = TRUE` they can be a little
longer than the times printed during the fit, which time the work alone.
[`sampling()`](http://haziqj.ml/inlavanilla/reference/sampling.md) and
[`loo()`](http://haziqj.ml/inlavanilla/reference/loo.md) run after the
fit, so their time is not included.

INLAvaan also exports a `timing()` generic, an S4 one. When both
packages are attached, the one attached last masks the other. With
inlavanilla attached last, its `timing()` passes an INLAvaan fit on to
INLAvaan's. With INLAvaan attached last, call `inlavanilla::timing(fit)`
for an inlavanilla fit.

## See also

[`summary.inla_vanilla()`](http://haziqj.ml/inlavanilla/reference/summary.inla_vanilla.md),
which gives the time of each stage with the number of calls.

## Examples

``` r
set.seed(1)
y <- rnorm(30, mean = 5, sd = 2)
logpost <- function(theta, y) {
  sum(dnorm(y, theta[["mu"]], exp(theta[["log_sigma"]]), log = TRUE)) +
    dnorm(theta[["mu"]], 0, 10, log = TRUE) +
    dnorm(theta[["log_sigma"]], 0, 2, log = TRUE)
}
fit <- inla_vanilla(logpost, c(mu = 0, log_sigma = 0), y = y, silent = TRUE)

# The whole fit
timing(fit)
#> total 
#>  14ms 

# Every stage
timing(fit, what = "all")
#>        init     laplace          vb   marginals diagnostics       total 
#>         0ms         2ms         0ms        11ms         1ms        14ms 

# Some stages
timing(fit, what = c("vb", "marginals"))
#>        vb marginals 
#>       0ms      11ms 
```
