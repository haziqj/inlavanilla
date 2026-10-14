# Plot an INLAvanilla fit

By default (`type = "marg_pdf"`), draws the marginal density of each
parameter, on the unconstrained scale, with the Gaussian marginal of the
Laplace approximation for comparison. The vertical line marks the
posterior mean. For a fit with `marginals = "gaussian"` the marginals
are Gaussian: centred at the VB mean, or, with `vb = FALSE`, the Laplace
marginals themselves, drawn without a second copy.

## Usage

``` r
# S3 method for class 'inla_vanilla'
plot(
  x,
  pars = NULL,
  type = c("marg_pdf", "sn_fit", "sn_fit_log"),
  laplace = TRUE,
  ...
)
```

## Arguments

- x:

  An object of class `inla_vanilla`.

- pars:

  Parameters to plot, by name or position. Defaults to all.

- type:

  The plot: `"marg_pdf"` (the default) for the marginal densities,
  `"sn_fit"` for each skew-normal fit against its scan, or
  `"sn_fit_log"` for the same on a log scale, as INLAvaan names them.

- laplace:

  Logical. Overlay the Laplace (Gaussian) marginals? Default `TRUE`;
  ignored when the marginals are the Laplace ones, and by the `"sn_fit"`
  plots.

- ...:

  Further graphical parameters passed to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html) for
  each panel, such as `xlab` or `main`; they replace the defaults (the
  parameter name as the title, no axis labels, and for `"marg_pdf"` no y
  axis).

## Value

`x`, invisibly.

## Details

`type = "sn_fit"` shows each skew-normal fit against the scan it was
fitted to, as INLAvaan's `plot(type = "sn_fit")` does. The x axis is the
scan coordinate, in posterior SDs from the mode (the Laplace SDs, from
-4 to 4). "Raw" is the log-posterior along the scan and "Corrected" the
same with the volume correction `gamma1` added, each as a density scaled
to a peak of 1, and "Skew-normal fit" is the fitted skew-normal at the
height it was fitted at. The largest gap between the last two, as a
fraction of the peak, is the NMAD of
[`summary.inla_vanilla()`](http://haziqj.ml/inlavanilla/reference/summary.inla_vanilla.md).
`type = "sn_fit_log"` draws the same on a base-10 log scale, which shows
the tails. Both use the scans kept in the fit, with no further calls of
`logpost`, and need skew-normal marginals.
