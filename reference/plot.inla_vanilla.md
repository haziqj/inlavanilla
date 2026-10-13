# Plot the posterior marginals of an INLAvanilla fit

Draws the skew-normal marginal density of each parameter, on the
unconstrained scale, with the Gaussian marginal of the Laplace
approximation for comparison. The vertical line marks the posterior
mean.

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

  Logical. Overlay the Laplace (Gaussian) marginals? Default `TRUE`.

- ...:

  Further graphical parameters passed to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html) for
  each panel, such as `xlab` or `main`; they replace the defaults (the
  parameter name as the title, no axis labels, no y axis).

## Value

`x`, invisibly.
