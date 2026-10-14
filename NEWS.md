# inlavanilla 0.1.1

* `inla_vanilla()` takes the Hessian of `logpost` as an optional `hess`,
  after `grad`, as `nlminb()` does. With it, `nlminb()` takes Newton steps to
  the mode, the Hessian at the mode is the supplied one, and the volume
  correction of each skew-normal marginal takes one `hess` call instead of
  2m + 2 `grad` calls. RTMB gives it with the gradient, as `obj$he`. The fit
  is the same, except that the mode is found more precisely, and on test
  models it took 13% to 35% less time. A supplied `hess` is checked at the
  mode against central differences of `grad`, and `fit$counts` counts its
  calls.
* The website article "Gradients with RTMB" is now "Gradients and Hessians
  with RTMB". It passes RTMB's Hessian too, and shows where the fit uses it,
  the time it saves and the more precise mode it finds.
* The website article "Diagnostic classification models" is now "Cognitive
  diagnostic models".

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
  degree four in whitened coordinates (six with a gradient). On test models
  it was a median 0.015 posterior SDs from exact VB, but up to about 0.1 on
  strongly skewed posteriors, where a gradient gives about 0.03.
* Two website articles, "Speeding up with Rcpp" and "Gradients with RTMB",
  fit a structural equation model with 39 parameters from a log-posterior
  written in R, in C++ with an analytic gradient, and in R for RTMB, which
  gives the gradient by automatic differentiation, and compare the routes
  for speed and accuracy.
* A website article, "Multilevel models", fits a linear mixed model with the
  random effects integrated out, checks its marginals against brms, and draws
  the random effects for a caterpillar plot.
* A website article, "Diagnostic classification models", fits the DINA model
  to the ECPE data with the latent classes summed out, checks its marginals
  against Stan, and classifies examinees from draws.
* A website article, "Correlation matrices", puts an LKJ prior on the
  correlation matrix of a multivariate normal, through its canonical partial
  correlations on the atanh scale, fits it to the Holzinger and Swineford
  tests, checks the marginals against Stan, and gives the reliability of each
  scale from draws.
* A website article, "Binary factor analysis", fits a one-factor model for
  binary items in the delta parameterisation, with the factor integrated out
  by Gauss-Hermite quadrature over the counts of the response patterns.
* A website article, "Diagnostic tests without a gold standard", fits the
  three latent class models of Koukounari et al. (2021) to five
  schistosomiasis tests: conditional independence, a person-level random
  effect integrated out by Gauss-Hermite quadrature, and a finite mixture with
  always-positive and always-negative groups. It compares their estimates with
  the published ones, and the models by their log marginal likelihoods.
* Three output levels: by default each stage prints one line; `verbose =
  TRUE` adds a full diagnostic of each stage (the optimiser's report, how the
  Hessian was taken, the posterior correlations, the VB iterations, a line per
  marginal with the flagged values marked, and the calls each stage made);
  `silent = TRUE` prints nothing and muffles warnings, for simulations. The
  default output is what `verbose = TRUE` used to give, so code that passed
  `verbose = FALSE` to keep a fit quiet should now pass `silent = TRUE`. The
  warnings raised during a fit are kept in `fit$warnings` whatever the level.
* While a stage runs, its line is the stage's sentence with a counter, as in
  INLAvaan ("Fitting 4/6 skew-normal marginals"), not a bar as wide as the
  console. After two seconds it ends with the time so far, or the time left
  where the work is known in advance (the marginals and `loo()`), in the
  place where the finished line shows the time taken. The lines, `loo()`'s
  included, follow the `cli.progress_show_after` option when it is set (half
  a second otherwise), so a knitr or Quarto setup chunk can keep them out of
  a rendered document.
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
  its message and a hint), an NMAD above 0.1, an unscanned mass above 0.05, a
  VB shift above 1 SD, and a Hessian condition number above 1e8; heavy tails
  are inlavanilla's own check. Like INLAvaan, the fit only records the VB
  fallback. `summary()` flags the same checks at the same levels, and also
  a ridge bend above 0.25 nats, where the marginal SD may be about 10% too
  small, and a mean gap (between the skew-normal and VB means) above 0.1
  SDs; its lines list every flagged parameter, worst first, with its value
  (names only for the heavy-tails check, which is yes or no).
  Healthy fits stay silent. A new vignette, "Diagnostics", explains every
  check, with citations of Jamil and Rue (2026a).
* `summary()` prints the fit first (the log marginal likelihood, and the
  calls and time each stage took), then the posterior marginals, then the
  diagnostics, each under its own heading. Its table of diagnostics (VB
  shift, NMAD, unscanned mass, ridge bend and mean gap) has numbers
  formatted as a tibble prints them and the flagged values in red. One line
  gives the largest value of each skew-normal check that passes, and each
  flagged check has a line of its own.
* The tables of `print()` and `summary()` number their rows, as a tibble
  does, and line up the decimal points down each column. The posterior
  marginals still give each parameter to three significant figures of its
  posterior SD, and a value that rounds to zero no longer prints as
  "-0.000".
* `print()` and `summary()` open with the log marginal likelihood and whether
  the fit converged, with the iterations and the time the fit took.
  Converged means that the computation reached its solution: no further
  Newton step at the mode, and no further fixed-point step at the VB mean,
  would move a parameter by more than 0.1 posterior SDs. The second is a new
  check, `vb_step_sigma` in `diagnostics(fit, type = "param")`: where the
  fallback optimiser of the VB correction lands was not checked before. A
  correction that has not settled is flagged and raises a warning.
* Printed output is never wider than 80 columns, or the console if it is
  narrower: the summaries, the verbose report and the fit-time warning wrap
  their long lines. The rendered documentation is 76 columns wide, so that
  no output scrolls sideways on the website.
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
* The VB correction no longer stops with "NA/NaN gradient evaluation" when
  `logpost` is `-Inf` at some of its nodes, which reach up to 3 posterior SDs
  from the mode, as near a boundary of the support. Its fallback optimiser
  now moves away from those nodes. If it finds a VB mean where `logpost` is
  finite at every node, the fit goes ahead. Otherwise, as for a posterior
  truncated within reach of the nodes, the fit stops with an error that says
  so and suggests a reparameterisation or `vb = FALSE`.
* With a gradient, a VB fallback that lands short of settled is rerun once
  from where it landed. That happens when the optimiser starts where a node
  fails, and it then reported convergence well short of the VB mean: 0.29
  posterior SDs on a one-parameter posterior with a boundary. Fits whose
  fallback settles are unchanged, and `verbose = TRUE` reports the rerun.
* A VB fixed-point step over 1.5 posterior SDs is now cut to 1.5, and the
  iteration goes on; only two such steps in a row hand over to the fallback
  optimiser. INLAvaan hands over at the first step over 1 SD, but on a skewed
  posterior the first step overshoots and the iteration then comes back. On
  a DINA model with 63 parameters and a gradient the handover cost 54,000
  calls to reach the point that six fixed-point steps reach.
