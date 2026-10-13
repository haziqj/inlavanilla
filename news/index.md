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
