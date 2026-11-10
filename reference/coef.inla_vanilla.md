# Posterior means of an INLAvanilla fit

The posterior mean of each parameter, from its skew-normal marginal: the
`Mean` column of the fit's summary table. A parameter named for its
scale (`log_sigma`) is on its natural scale, under its natural name
(`sigma`), unless the fit was made with `transform = FALSE`; the others
are on the scale of `logpost`. These are the estimates INLAvaan's
[`coef()`](https://rdrr.io/r/stats/coef.html) gives. For a fit with
`marginals = "gaussian"` and no transformed parameters they are the
VB-corrected means, or, with `vb = "off"`, the mode. The joint posterior
mode is `fit$mode`, and the VB-corrected mean `fit$mean_vb`, both on the
scale of `logpost`.

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
