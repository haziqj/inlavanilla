# Plot the posterior marginals of an INLAvanilla fit

Draws the marginal density of each parameter, on the unconstrained
scale, with the Gaussian marginal of the Laplace approximation for
comparison. The vertical line marks the posterior mean. For a fit with
`marginals = "gaussian"` the marginals are Gaussian: centred at the VB
mean, or, with `vb = FALSE`, the Laplace marginals themselves, drawn
without a second copy.

## Usage

``` r
# S3 method for class 'inla_vanilla'
plot(x, pars = NULL, laplace = TRUE, ...)
```

## Arguments

- x:

  An object of class `inla_vanilla`.

- pars:

  Parameters to plot, by name or position. Defaults to all.

- laplace:

  Logical. Overlay the Laplace (Gaussian) marginals? Default `TRUE`;
  ignored when the marginals are the Laplace ones.

- ...:

  Further graphical parameters passed to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html) for
  each panel, such as `xlab` or `main`; they replace the defaults (the
  parameter name as the title, no axis labels, no y axis).

## Value

`x`, invisibly.
