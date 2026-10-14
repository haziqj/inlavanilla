# inlavanilla 0.1.0

* Initial version. `inla_vanilla()` approximates the posterior of any
  log-posterior written on an unconstrained scale with the INLAvaan pipeline:
  a Laplace approximation, a VB mean correction (by a deterministic degree-5
  cubature rule), skew-normal marginals with the volume correction, and the
  log marginal likelihood. `sampling()` gives Gaussian-copula draws of the
  parameters or of any function of them. `print()`, `summary()` and `plot()`
  methods report the fit and its diagnostics.
* `loo()` and `waic()` estimate predictive accuracy from the fit and a
  pointwise log-likelihood `loglik(theta, ...)`, by the deterministic
  second-order Taylor approximation of Alhyari et al. (2026): no refits and no
  draws. The results work with `loo::loo_compare()`. A new vignette,
  "Model comparison with LOO and WAIC", shows how. Both follow the fit's
  `silent`, for the warnings of `loglik` too, and pass any other object on to
  the loo package's functions of the same name.
* The VB correction (`vb`) and the skew-normal marginals (`marginals =
  "skewnormal"` or `"gaussian"`) can each be switched off, down to the
  Laplace approximation alone.
* The VB correction uses a degree-5 rule built on a regular simplex, with at
  most (m + 1)(m + 2) + 1 nodes, against 2m^2 + 1 before, all within 3
  posterior SDs. With a gradient it makes about half the calls at large m,
  with the same accuracy. Without one, it takes the expected gradient from
  values of `logpost` alone, by Stein's identity on the same nodes, instead of
  finite-difference gradients: at most (m + 1)(m + 2) calls per sweep instead
  of (m + 1)(2m^2 + 1), so it costs about as much as finding the mode. For the
  39-parameter SEM of the Rcpp article that is 8,200 calls instead of
  487,000. This VB mean is exact when the log-posterior is a polynomial of
  degree four in whitened coordinates (six with a gradient), and typically
  about 0.01 posterior SDs from exact VB.
* Two website articles, "Speeding up with Rcpp" and "Using TMB", fit a
  structural equation model with 39 parameters from a log-posterior written
  in R, in C++ with an analytic gradient, and as a TMB template.
* Three output levels: by default each stage prints one line; `verbose =
  TRUE` adds a full diagnostic of each stage (the optimiser's report, how the
  Hessian was taken, the posterior correlations, the VB iterations, a line per
  marginal with the flagged values marked, and the calls each stage made);
  `silent = TRUE` prints nothing and muffles warnings, for simulations. The
  default output is what `verbose = TRUE` used to give, so code that passed
  `verbose = FALSE` to keep a fit quiet should now pass `silent = TRUE`. The
  warnings raised during a fit are kept in `fit$warnings` whatever the level.
* Progress bars follow the `cli.progress_show_after` option when it is set
  (half a second otherwise), so a knitr or Quarto setup chunk can keep them
  out of a rendered document.
* The diagnostics follow INLAvaan's. Each parameter has the Newton step at
  the mode in posterior SDs (`mode_shift_sigma`), the VB shift on the scale
  of `logpost` (`vb_shift`) and in posterior SDs (`vb_shift_sigma`), with the
  Kullback-Leibler divergence of each shifted marginal (`kld`), and
  `unscanned_mass`, the mass each fitted skew-normal puts outside the scanned
  window; the VB columns are `NA` when `vb = FALSE`. The fit also records
  the condition number of the Hessian in posterior-SD units (`hess_cond`).
  `summary()` and `verbose = TRUE` show the unscanned mass.
* A fit whose checks look off ends with one warning, of class
  `inlavanilla_diagnostics_warning`, headed "Fit diagnostics flagged N
  potential issues", in INLAvaan's words and at its levels: a Newton step at
  the mode above 0.1 posterior SDs (whatever `nlminb()`'s return code, with
  its message and a hint), an NMAD above 0.1, a unscanned mass above 0.05, a
  VB shift above 1 SD, and a Hessian condition number above 1e8; heavy tails
  are inlavanilla's own check. Like INLAvaan, the fit only records the VB
  fallback. `summary()` flags the same checks at the same levels, and also
  a ridge bend above 0.25 nats, where the marginal SD may be about 10% too
  small, and the gap between the skew-normal and VB means; its
  lines list every flagged parameter, worst first, with its value. Healthy
  fits stay silent. A new vignette, "Diagnostics", explains every check,
  with citations of Jamil and Rue (2026a).
* `diagnostics(fit)` and `diagnostics(fit, type = "param")` return the
  checks as numbers, with the names, order and print formats of INLAvaan's
  `diagnostics()`, followed by inlavanilla's own. `plot(fit, type =
  "sn_fit")` and `type = "sn_fit_log"` draw each skew-normal fit against the
  scan it was fitted to, as INLAvaan's plot does, from the scans the fit now
  keeps in `fit$scan`.
* The check of a user `grad` also compares it with finite differences one
  curvature scale from `init`, so a gradient that is wrong by a constant
  factor is caught even when `init` is at the mode.
* `print()` names the approximation on its header rule ("Skew-normal + VB",
  "Skew-normal", "Laplace + VB" or "Laplace"), as `summary()` now does,
  instead of the run time, which `summary()` still reports. It ends with a
  line when the fit raised warnings and one when `summary()` flags any check,
  so a fit run with `silent = TRUE` does not look clean when it is not.
* `summary()` reports the 2.5%, 50% and 97.5% quantiles by default, as
  R-INLA does, instead of adding the quartiles, so its table fits in 80
  columns and matches `print()`. Pass `probs` for others.
* `print()` and `summary()` wrap long lines to the console width (the
  `width` option in R Markdown and Quarto), and split a table too wide for
  it into blocks of columns, as R does for a matrix.
