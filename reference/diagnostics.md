# Diagnostics of an INLAvanilla fit

The checks of each stage of the approximation, as numbers: a summary of
the whole fit, or a table with one row per parameter. They are the
checks of INLAvaan's `diagnostics()`, under the same names and in the
same order, followed by inlavanilla's own. The fit-time warning and the
flags of
[`summary.inla_vanilla()`](http://haziqj.ml/inlavanilla/reference/summary.inla_vanilla.md)
are computed from the same numbers, and
[`summary.inla_vanilla()`](http://haziqj.ml/inlavanilla/reference/summary.inla_vanilla.md)
explains each check and the level at which it is flagged.

Jamil and Rue (2026) use the VB shift and NMAD as fit diagnostics
(Sections 3.1 and 4.4), and the gap between the skew-normal and VB means
(`mean_gap`) as an internal check (Section 5.1). The mode shift, the
Hessian condition number and the unscanned mass come from INLAvaan's
`diagnostics()`, not from the paper. `vb_step_sigma`, `tails` and
`ridge_bend` are inlavanilla's own.

## Usage

``` r
diagnostics(object, ...)

# Default S3 method
diagnostics(object, ...)

# S3 method for class 'inla_vanilla'
diagnostics(object, type = c("global", "param"), ...)
```

## Arguments

- object:

  An object of class `inla_vanilla`.

- ...:

  Unused.

- type:

  `"global"` (the default) for the diagnostics of the whole fit, or
  `"param"` for those of each parameter.

## Value

For `type = "global"`, a named numeric vector of class
`inla_vanilla_diagnostics`. For `type = "param"`, a data frame of class
`inla_vanilla_diagnostics_param`. Their
[`print()`](https://rdrr.io/r/base/print.html) methods use INLAvaan's
formats.

## Details

A posterior SD here is the Laplace SD at the mode, \\s_j =
\sqrt{\Sigma\_{jj}}\\, as in INLAvaan; it is not the SD column of the
summary table, which is the skew-normal's. Values that a fit did not
compute are `NA`: the VB quantities with `vb = FALSE`, and the scan
quantities with `marginals = "gaussian"`.

**Global diagnostics** (`type = "global"`), a named numeric vector:

- `npar`: the number of parameters.

- `converged`: 1 if
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html) reported
  convergence (return code 0), 0 otherwise, as INLAvaan defines it. No
  check uses it: without a gradient
  [`nlminb()`](https://rdrr.io/r/stats/nlminb.html) can report false
  convergence at an accurate mode. `mode_shift_max` checks the mode
  instead.

- `iterations`: the iterations of
  [`nlminb()`](https://rdrr.io/r/stats/nlminb.html), in the run that was
  kept if it was restarted.

- `grad_inf`, `grad_inf_rel`, `grad_l2`: the largest absolute value, the
  largest relative value and the Euclidean norm of `grad` (see below).
  They depend on the scale of each parameter.

- `mode_shift_max`: the largest `mode_shift_sigma`. Unlike the gradient
  norms it is scale-free. It should be close to 0.

- `hess_cond`: the condition number of the Hessian at the mode, with
  each parameter in its own posterior-SD units: the ratio of the largest
  to the smallest eigenvalue of the posterior correlation matrix. Large
  values mean that some combination of the parameters is weakly
  identified. INLAvaan takes it on the raw scale of the parameters,
  where a mere change of units makes it large, so the two agree only
  when the posterior SDs are similar.

- `vb_applied`: 1 if the VB correction ran, 0 otherwise.

- `vb_shift_max`: the largest absolute `vb_shift_sigma`.

- `vb_kld_global`: the Kullback-Leibler divergence between the Laplace
  Gaussian and its VB-shifted copy, \\\frac12 \delta^\top H \delta\\,
  which is also the VB term of the log marginal likelihood.

- `kld_max`, `kld_mean`: the largest and the mean `kld`.

- `nmad_max`, `nmad_mean`: the largest and the mean `nmad`.

- `unscanned_mass_max`: the largest `unscanned_mass`.

and then inlavanilla's own:

- `vb_fallback`: 1 if the VB fixed-point iteration handed over to the
  fallback optimiser, 0 otherwise. Like INLAvaan, inlavanilla records it
  but flags nothing: the fallback is a route of the solver, not a
  measure of the approximation.

- `vb_step_max`: the largest `vb_step_sigma`. Like `mode_shift_max` it
  should be close to 0.

- `mean_gap_max`: the largest absolute `mean_gap`.

- `ridge_bend_max`: the largest `ridge_bend`.

- `n_heavy_tails`: the number of parameters whose `tails` is `FALSE`.

**Per-parameter diagnostics** (`type = "param"`), a data frame:

- `names`: the parameter.

- `grad`: the gradient of minus the log-posterior at the mode, as
  INLAvaan reports it: from `grad` when the fit has one, and otherwise
  from central differences of `logpost` a thousandth of a posterior SD
  apart. It should be close to 0.

- `grad_abs`: its absolute value.

- `grad_rel`: `grad_abs / (abs(par) + 1e-6)`, with `par` the VB mean
  (the mode when `vb = FALSE`), as in INLAvaan.

- `mode_shift_sigma`: the Newton step at the mode, in posterior SDs,
  \\\|(\Sigma g)\_j\| / s_j\\ with \\g\\ the gradient. It should be
  close to 0.

- `kld`: the Kullback-Leibler divergence between the Laplace Gaussian
  marginal and its VB-shifted copy, \\\delta_j^2 / (2\Sigma\_{jj})\\ =
  `vb_shift_sigma^2 / 2`. It is a function of the VB shift alone, and
  says nothing about the skew-normal fit.

- `vb_shift`: the VB shift of the mean, \\\delta_j\\, the VB mean minus
  the mode, on the scale of `logpost`.

- `vb_shift_sigma`: the same in posterior SDs.

- `nmad`: the NMAD of the skew-normal fit.

- `alpha`: the shape of the fitted skew-normal marginal.

- `unscanned_mass`: the mass the fitted skew-normal puts beyond the
  scanned window.

and then inlavanilla's own:

- `vb_step_sigma`: the step still to go at the VB mean, in posterior
  SDs: the fixed-point step there, the Newton step of the VB
  stationarity condition, as `mode_shift_sigma` is for the mode. After
  the fallback optimiser it is taken where the optimiser landed. On the
  fixed-point path it is the last step taken, which bounds the next,
  unless that step was 0.1 or more (a `vb_tol` above 0.1), when it is
  taken afresh. It should be close to 0 (`Inf` where the step is not
  finite).

- `mean_gap`, `gamma1`, `tails` and `ridge_bend`, which
  [`summary.inla_vanilla()`](http://haziqj.ml/inlavanilla/reference/summary.inla_vanilla.md)
  describes.

Some of INLAvaan's diagnostics are left out. `nsamp` counts posterior
draws, which
[`sampling()`](http://haziqj.ml/inlavanilla/reference/sampling.md) makes
on demand rather than during the fit. `grad_num` and `grad_diff` compare
an analytic gradient with finite differences at the mode, which would
take further calls; a supplied `grad` is checked at `init` and through
the symmetry of its Jacobian at the mode instead, with a warning if
either check fails. `vb_mcse_max`, `vb_mcse_mean` and `vb_mcse_sigma`
estimate the quadrature error of INLAvaan's Sobol points, and the
cubature rule of inlavanilla is deterministic (INLAvaan gives `NA` for
its own deterministic rule).

INLAvaan also exports a `diagnostics()` generic, an S4 one. When both
packages are attached, the one attached last masks the other. With
inlavanilla attached last, its `diagnostics()` passes an INLAvaan fit on
to INLAvaan's. With INLAvaan attached last, call
`inlavanilla::diagnostics(fit)` for an inlavanilla fit.

## References

Jamil, H., & Rue, H. (2026). *Approximate Bayesian inference for
structural equation models using integrated nested Laplace
approximations* (arXiv:2603.25690 \[stat.ME\]). arXiv.
<https://doi.org/10.48550/arXiv.2603.25690>

## See also

[`summary.inla_vanilla()`](http://haziqj.ml/inlavanilla/reference/summary.inla_vanilla.md)
for what each check means and when it is flagged, and
[`plot.inla_vanilla()`](http://haziqj.ml/inlavanilla/reference/plot.inla_vanilla.md)
with `type = "sn_fit"` for the skew-normal fits themselves.

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
diagnostics(fit)
#>               npar          converged         iterations           grad_inf 
#>                  2                  1                 10           1.96e-05 
#>       grad_inf_rel            grad_l2     mode_shift_max          hess_cond 
#>           3.13e-05           2.05e-05           2.53e-06           1.01e+00 
#>         vb_applied       vb_shift_max      vb_kld_global            kld_max 
#>                  1             0.2510             0.0315             0.0315 
#>           kld_mean           nmad_max          nmad_mean unscanned_mass_max 
#>             0.0158             0.0077             0.0039           6.58e-04 
#>        vb_fallback        vb_step_max       mean_gap_max     ridge_bend_max 
#>                  0           2.64e-04             0.0124             0.1324 
#>      n_heavy_tails 
#>                  0 

# One row per parameter
diagnostics(fit, type = "param")
#>       names grad grad_abs grad_rel mode_shift_sigma    kld vb_shift
#> 1        mu    0        0        0                0 0.0000  -0.0006
#> 2 log_sigma    0        0        0                0 0.0315   0.0323
#>   vb_shift_sigma   nmad   alpha unscanned_mass vb_step_sigma mean_gap
#> 1        -0.0017 0.0000 -0.2519       6.36e-05         0e+00  -0.0002
#> 2         0.2510 0.0077  1.3709       6.58e-04         3e-04   0.0124
#>    gamma1 tails ridge_bend
#> 1 -0.0002  TRUE     0.1324
#> 2  0.1282  TRUE     0.0000
```
