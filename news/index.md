# Changelog

## inlavanilla 0.1.0

- Initial version.
  [`inla_vanilla()`](http://haziqj.ml/inlavanilla/reference/inla_vanilla.md)
  approximates the posterior of any log-posterior written on an
  unconstrained scale with the INLAvaan pipeline: a Laplace
  approximation, a VB mean correction (by a deterministic degree-5
  cubature rule), skew-normal marginals with the volume correction, and
  the log marginal likelihood.
  [`sampling()`](http://haziqj.ml/inlavanilla/reference/sampling.md)
  gives Gaussian-copula draws of the parameters or of any function of
  them. [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods
  report the fit and its diagnostics.
- [`loo()`](http://haziqj.ml/inlavanilla/reference/loo.md) and
  [`waic()`](http://haziqj.ml/inlavanilla/reference/loo.md) estimate
  predictive accuracy from the fit and a pointwise log-likelihood
  `loglik(theta, ...)`, by the deterministic second-order Taylor
  approximation of Alhyari et al. (2026): no refits and no draws. The
  results work with
  [`loo::loo_compare()`](https://mc-stan.org/loo/reference/loo_compare.html).
  A new vignette, “Model comparison with LOO and WAIC”, shows how. Both
  follow the fit’s `silent`, for the warnings of `loglik` too, and pass
  any other object on to the loo package’s functions of the same name.
- The VB correction (`vb`) and the skew-normal marginals
  (`marginals = "skewnormal"` or `"gaussian"`) can each be switched off,
  down to the Laplace approximation alone.
- Two website articles, “Speeding up with Rcpp” and “Using TMB”, fit a
  structural equation model with 39 parameters from a log-posterior
  written in R, in C++ with an analytic gradient, and as a TMB template.
- Three output levels: by default each stage prints one line;
  `verbose = TRUE` adds a full diagnostic of each stage (the optimiser’s
  report, how the Hessian was taken, the posterior correlations, the VB
  iterations, a line per marginal with the flagged values marked, and
  the calls each stage made); `silent = TRUE` prints nothing and muffles
  warnings, for simulations. The default output is what `verbose = TRUE`
  used to give, so code that passed `verbose = FALSE` to keep a fit
  quiet should now pass `silent = TRUE`. The warnings raised during a
  fit are kept in `fit$warnings` whatever the level.
- Progress bars follow the `cli.progress_show_after` option when it is
  set (half a second otherwise), so a knitr or Quarto setup chunk can
  keep them out of a rendered document.
- The marginal diagnostics include `scan_end_mass`, the mass each fitted
  skew-normal puts outside the scanned window, as in INLAvaan’s
  `diagnostics()`; [`summary()`](https://rdrr.io/r/base/summary.html)
  flags it above 0.05.
- A fit whose marginals are clearly off now ends with one warning, of
  class `inlavanilla_diagnostics_warning`: heavy tails, an NMAD above
  0.1, a scan-end mass above 0.05, a VB shift above 1 SD or the VB
  fallback, at INLAvaan’s levels. Healthy fits stay silent. A new
  vignette, “Diagnostics”, explains every check.
- The check of a user `grad` also compares it with finite differences
  one curvature scale from `init`, so a gradient that is wrong by a
  constant factor is caught even when `init` is at the mode.
