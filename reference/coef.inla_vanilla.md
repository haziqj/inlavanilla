# Posterior means of an INLAvanilla fit

The posterior mean of each parameter, from its skew-normal marginal, on
the unconstrained scale of `logpost`: the `Mean` column of
[`summary.inla_vanilla()`](http://haziqj.ml/inlavanilla/reference/summary.inla_vanilla.md).
These are the estimates INLAvaan's
[`coef()`](https://rdrr.io/r/stats/coef.html) gives. For a fit with
`marginals = "gaussian"` they are the VB-corrected means, or, with
`vb = "off"`, the mode. The joint posterior mode is `fit$mode`, and the
VB-corrected mean `fit$mean_vb`.

## Usage

``` r
# S3 method for class 'inla_vanilla'
coef(object, ...)
```

## Arguments

- object:

  An object of class `inla_vanilla`.

- ...:

  Unused.

## Value

A named numeric vector, with one posterior mean per parameter.
