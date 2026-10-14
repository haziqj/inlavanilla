# INLA-style posterior approximation for any log-posterior

Approximates the posterior of a parameter vector from its log density
alone, using the approximation pipeline of INLAvaan (Jamil & Rue, 2026a,
2026b) with no latent field and nothing model-specific. The pipeline has
four stages:

1.  **Joint Gaussian (Laplace).** The posterior mode \\\vartheta^\*\\
    and the inverse Hessian \\\Omega\\ there give a Gaussian
    approximation \\N(\vartheta^\*, \Omega)\\ and the Laplace estimate
    of the log marginal likelihood.

2.  **VB mean correction.** The Gaussian is recentred at
    \\\tilde\vartheta = \vartheta^\* + \hat\delta\\, the mean that
    maximises its expected log-posterior, with the covariance unchanged.

3.  **Skew-normal marginals.** For each parameter, the log-posterior is
    scanned along the conditional mode path of the Gaussian, corrected
    for the changing volume of the other parameters, and fitted by a
    skew-normal \\SN(\hat\xi_j, \hat\omega_j, \hat\alpha_j)\\.

4.  **Gaussian copula draws.**
    [`sampling()`](http://haziqj.ml/inlavanilla/reference/sampling.md)
    couples the skew-normal marginals through the correlations of
    \\\Omega\\.

Stages 2 and 3 are optional and independent of each other. `vb = FALSE`
skips the VB correction, which leaves the marginals unchanged (they are
scanned around the mode) but gives the Laplace log marginal likelihood.
`marginals = "gaussian"` skips the scans and uses the Gaussian's own
marginals \\N(\tilde\vartheta_j, \Omega\_{jj})\\, centred at the VB mean
(or at the mode when `vb = FALSE`);
[`sampling()`](http://haziqj.ml/inlavanilla/reference/sampling.md) then
draws from the joint Gaussian. So `vb = FALSE, marginals = "gaussian"`
is the Laplace approximation alone, and `marginals = "gaussian"` alone
is the Laplace approximation with the VB correction.

## Usage

``` r
inla_vanilla(
  logpost,
  init,
  grad = NULL,
  hess = NULL,
  ...,
  vb = TRUE,
  marginals = c("skewnormal", "gaussian"),
  verbose = FALSE,
  silent = FALSE,
  control = list()
)
```

## Arguments

- logpost:

  A function `logpost(theta, ...)` returning the log-posterior density
  (up to a constant) at the named numeric vector `theta`, on the
  unconstrained scale. See the section on writing `logpost`.

- init:

  A named numeric vector of starting values. Its length sets the number
  of parameters and its names label them; unnamed values are called
  `theta[1]`, `theta[2]`, ..., and a repeated name is numbered in the
  same way.

- grad:

  Optional function `grad(theta, ...)` returning the gradient of
  `logpost` as a numeric vector. When `NULL` (the default), finite
  differences are used.

- hess:

  Optional function `hess(theta, ...)` returning the Hessian of
  `logpost`, the m x m matrix of its second derivatives (a sparse matrix
  from the Matrix package, and a number when m = 1, work too). It needs
  `grad`. When `NULL` (the default), finite differences of `grad` are
  used.

- ...:

  Further arguments passed to `logpost`, `grad` and `hess`, such as the
  data. As in [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html),
  they come after the gradient and the Hessian, and the arguments after
  them (`vb`, `marginals`, `verbose`, `silent`, `control`) must be given
  by their full names, so data arguments cannot partially match them.
  The arguments before them can be matched, as in
  [`nlminb()`](https://rdrr.io/r/stats/nlminb.html): a data argument
  named `g` or `gr` would be taken as `grad`, and one named `h` or `he`
  as `hess`, so give it a longer name.

- vb:

  Logical. Apply the VB mean correction? Default `TRUE`. It enters the
  log marginal likelihood, and the skew-normal means are checked against
  it (see
  [`summary.inla_vanilla()`](http://haziqj.ml/inlavanilla/reference/summary.inla_vanilla.md)).

- marginals:

  The marginal of each parameter: `"skewnormal"` (the default) fits a
  skew-normal to a scan of the log-posterior; `"gaussian"` uses the
  Gaussian marginal of the Laplace approximation, centred at the VB mean
  when `vb = TRUE`, with no scans. Abbreviations such as `"skew"` and
  `"gauss"` work too.

- verbose:

  Logical. Print a full diagnostic of each stage under its progress
  line? Default `FALSE`. See the section on output.

- silent:

  Logical. Print nothing at all, and muffle the warnings and messages
  raised during the fit? Default `FALSE`. The warnings are still
  recorded in `fit$warnings`. Cannot be combined with `verbose = TRUE`.

- control:

  A list of tuning parameters; see the section on control.

## Value

An object of class `inla_vanilla`, with methods for
[`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`plot()`](https://rdrr.io/r/graphics/plot.default.html),
[`sampling()`](http://haziqj.ml/inlavanilla/reference/sampling.md) and
[`diagnostics()`](http://haziqj.ml/inlavanilla/reference/diagnostics.md).
It is a list holding, among other things:

- `summary`: the posterior summary table (mean, SD, quantiles, mode).

- `mode`, `cov`: the Laplace mode and covariance.

- `mean_vb`: the VB-corrected mean (the mode when `vb = FALSE`).

- `sn`: the skew-normal parameters `xi`, `omega` and `alpha` per
  parameter (with `alpha = 0` for Gaussian marginals).

- `vb`, `marginals`, `verbose`, `silent`: the settings the fit used.

- `logml`: the log marginal likelihood, and `logml_laplace` its Laplace
  part.

- `diagnostics`: the checks of each stage, which
  [`diagnostics()`](http://haziqj.ml/inlavanilla/reference/diagnostics.md)
  returns under INLAvaan's names and
  [`summary.inla_vanilla()`](http://haziqj.ml/inlavanilla/reference/summary.inla_vanilla.md)
  explains. `optim` holds the return code, message and iterations of
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html), the gradient
  of minus the log-posterior at the mode (`grad`), the whitened gradient
  there (`grad_mode`), and the condition number of the Hessian in
  posterior-SD units (`hess_cond`). `vb` holds the iterations of the VB
  correction, whether it fell back to
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html) (`fallback`),
  and the Kullback-Leibler divergence of the shift (`kld_global`); it is
  `NULL` when `vb = FALSE`. `marginals` is a data frame with a row per
  parameter: `mode_shift_sigma`, `kld`, `vb_shift`, `vb_shift_sigma`,
  `nmad`, `alpha`, `unscanned_mass`, `vb_step_sigma`, `mean_gap`,
  `gamma1`, `tails` and `ridge_bend`, with `NA` for the checks the fit
  did not run.

- `scan`: the scans the skew-normals were fitted to, for
  `plot(fit, type = "sn_fit")`: the grid `z` in posterior SDs from the
  mode, the log-posterior at each point (`logpost`, a column per
  parameter), and the log scale factor `logC` of each fit on that grid.
  It is `NULL` for Gaussian marginals.

- `warnings`: the messages of the warnings raised during the fit, as a
  character vector of plain text, each message once (empty when there
  were none). They are recorded whether or not `silent` is `TRUE`.

## Writing `logpost`

The approximation works on an unconstrained scale: every parameter must
be free to take any real value. Reparameterise anything constrained, for
example a standard deviation \\\sigma \> 0\\ as \\\log\sigma\\, a
probability as its logit, or a correlation as its Fisher \\z\\. Means
and regression coefficients need no change.

Return `-Inf` wherever the density is zero or cannot be evaluated, such
as outside the support or where a matrix is not positive definite (`NA`
and `NaN` count the same). The scans and the VB nodes reach several
posterior SDs from the mode, and an error raised there stops the fit.

`logpost` is then the log density **of the parameters you pass**. When a
prior is stated on the natural scale, include the log-Jacobian of the
transformation: with \\\vartheta = \log\sigma\\ and a prior
\\p(\sigma)\\, the prior term is \\\log p(e^\vartheta) + \vartheta\\.
Leaving it out targets a different posterior, not a worse approximation
of the same one, and nothing downstream can detect it.

The choice of parameterisation also sets the quality of the
approximation, since each parameter gets its own skew-normal on the
scale you choose. Scales on which the posterior is close to Gaussian
(log standard deviations, Fisher-\\z\\ correlations, logits) work best.
Natural-scale or derived quantities are available through
`sampling(fit, fun = ...)`.

The log marginal likelihood is an estimate of \\\log p(y)\\ only when
`logpost` returns the full \\\log p(y, \vartheta)\\, including every
normalising constant of the likelihood and the prior. Otherwise it is
off by the missing constants.

## Accuracy

The skew-normal marginals capture skewness that the Laplace
approximation misses, and they are accurate to the extent that a
skew-normal fits each marginal (see the diagnostics in
[`summary.inla_vanilla()`](http://haziqj.ml/inlavanilla/reference/summary.inla_vanilla.md)).
Like INLAvaan, each marginal is scanned along a straight line through
the mode, and the changing volume of the other parameters enters only to
first order. When a parameter's spread depends strongly on another,
poorly determined parameter, its marginal comes out too narrow, by a
factor that vanishes as O(1/n). Examples are a regression coefficient or
a mean when the noise SD is uncertain. For a normal mean with unknown
variance the SD is about 11% too small at n = 12, 3% at n = 50, and
under 1% at n = 200. The `ridge_bend` diagnostic of
[`summary.inla_vanilla()`](http://haziqj.ml/inlavanilla/reference/summary.inla_vanilla.md)
flags this narrowing where it predicts about 10% or more.

A skew-normal also cannot follow a tail heavier than a Gaussian's. In a
structural equation model with 39 parameters and n = 75, the data could
not keep two residual variances away from zero, so on the log scale
their posteriors had a long, exponential left tail, set by the prior.
Checked against MCMC, those two skew-normal marginals were less than
half as wide as the posterior, with the 2.5% quantile about 2 posterior
SDs too high; the other 37 were within 7% in SD. On the variance scale,
where the tail near zero is compressed, the two were much closer: means
0.2 posterior SDs too high and SDs 10 to 15% too small. INLAvaan's
marginals of this model are the same as inlavanilla's: Jamil and Rue
(2026a, Section 4.4) find the same miss in the lower tail of a variance
near zero, the boundary pile-up that a unimodal skew-normal cannot
capture. No check is aimed at this. The `ridge_bend` diagnostic of
[`summary.inla_vanilla()`](http://haziqj.ml/inlavanilla/reference/summary.inla_vanilla.md)
flagged one of the two (0.44), but the other's ridge bend (0.11) was no
larger than that of marginals within 7% in SD, and the `tails` check,
which looks only 4 posterior SDs either side of the mode (the Laplace
SDs), missed both.

## Gradients and Hessians

`grad` is optional. Without it, derivatives come from finite differences
of `logpost`, and the number of calls grows with the square of the
number of parameters m. The mode and Hessian take about m^2 + 6m calls,
and the marginals about m(`n_grid` + 2m + 5). The VB correction takes at
most (m + 1)(m + 2) calls per sweep, over a few sweeps: it needs no
derivatives, because Stein's identity gives the expected gradient from
values of `logpost` alone. Its VB mean is then exact when the
log-posterior is a polynomial of degree four in whitened coordinates,
and up to degree six with `grad`. On test models it was about 0.01
posterior SDs from exact VB (0.003 with `grad`), and up to 0.1 on a
strongly skewed posterior (0.03 with `grad`), against a VB error of 0.01
to 0.2 SDs from the posterior mean itself. On the most skewed posteriors
the VB correction may then not settle, and the fit says it has not
converged (see
[`summary.inla_vanilla()`](http://haziqj.ml/inlavanilla/reference/summary.inla_vanilla.md));
`grad` is the remedy.

With `grad`, every stage is more accurate. The fit is also faster when a
gradient call costs about as much as a `logpost` call, as a hand-written
gradient usually does, since it then makes about half as many calls in
all. A gradient by automatic differentiation (RTMB, TMB) costs about
two, and the fit takes about as long as without one. The VB correction
takes at most (m + 1)(m + 2) + 1 gradient calls per sweep, at the same
nodes. Beyond a handful of parameters it and the marginals are the
slowest stages. `vb = FALSE` skips the VB correction (and its term in
the log marginal likelihood), and the skew-normal marginals do not
depend on it; `marginals = "gaussian"` saves the calls of the scans.

`hess` is optional too, and needs `grad`, as in
[`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html). With it,
[`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html) takes Newton
steps to the mode, the Hessian at the mode is the supplied one, and the
volume correction of each skew-normal marginal takes one `hess` call
instead of 2m + 2 `grad` calls. The VB correction and the scans do not
use it. The fit is the one `grad` alone gives, to the error of its
finite differences (about 1e-8, relative), except that the mode is found
more precisely, so `hess` saves time rather than adding accuracy. On
test models with an automatic-differentiation Hessian (RTMB) the fits
took 13% to 35% less time, the most where the Hessian was sparse.

A supplied `grad` is checked against finite differences at `init`, and
its Jacobian at the mode is checked for symmetry, with a warning if
either check fails. A supplied `hess` is checked at the mode against
central differences of `grad`, along two directions, with a warning if
it does not match. Any R function works for `logpost`, `grad` and
`hess`, including compiled code (Rcpp, TMB) and RTMB, which
differentiates R code. An RTMB function can return the log-posterior
itself, and its `obj$fn`, `obj$gr` and `obj$he` are then passed as they
are. For a TMB template that returns the negative log-posterior, pass
`function(x) -obj$fn(x)`, `function(x) -as.vector(obj$gr(x))` and
`function(x) -obj$he(x)`, with `init = obj$par`; its repeated names (one
per element of a parameter vector) are numbered, as `beta[1]`,
`beta[2]`, and so on. TMB gives no Hessian for a model with random
effects.

## Control

`control` is a list with any of:

- `n_grid`: number of grid points per skew-normal scan, spread over
  \\\pm 4\\ posterior SDs (default `21`).

- `log_threshold`: grid points whose log-posterior falls more than this
  below the maximum are ignored in the skew-normal fit (default `-6`).

- `vb_maxit`, `vb_tol`: iteration limit (default `25`) and convergence
  tolerance in posterior SDs (default `1e-3`) for the VB correction's
  fixed-point iteration. A step over 1.5 posterior SDs is cut to 1.5. If
  two steps in a row are cut, or it has not converged after `vb_maxit`
  steps, [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html) takes
  over.

- `nlminb`: a list of
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html) controls for
  the mode search. `iter.max` and `eval.max` default to `1000` and
  `2000`, and `trace = 1` prints its iterations.

## Output

By default each stage prints one line as it finishes, with the time it
took, and shows a progress bar while it runs for longer than half a
second, or than the option `cli.progress_show_after` when it is set. In
knitr or Quarto documents, where a bar cannot redraw in place, a larger
`cli.progress_show_after` in the setup chunk keeps bars out of the
output. Warnings, such as a gradient that does not match `logpost`, are
signalled as usual. A fit whose checks look off, such as a mode that was
not reached, ends with one warning that lists them, as INLAvaan's does
(see
[`summary.inla_vanilla()`](http://haziqj.ml/inlavanilla/reference/summary.inla_vanilla.md)
and
[`diagnostics()`](http://haziqj.ml/inlavanilla/reference/diagnostics.md)).

`verbose = TRUE` adds a full diagnostic of the fit, in grey lines under
each stage. It reports what the fit computes anyway, with no further
calls of `logpost` or `grad`:

- Before the mode search: the parameters, where the derivatives come
  from, `logpost` at `init`, and for a supplied `grad` its largest
  relative discrepancy from finite differences there.

- The mode: the message, iterations and evaluations of
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html), and whether
  it was restarted with the parameters rescaled; the log-posterior at
  the mode, the whitened gradient there (which decides the restart) and
  the largest Newton step in posterior SDs (which the fit-time warning
  checks); how the Hessian was taken (and, for a supplied `grad`, the
  asymmetry of its Jacobian, or, for a supplied `hess`, how far it is
  from central differences of `grad`); the largest posterior correlation
  and the smallest eigenvalue of the correlation matrix, which is near
  zero when some combination of the parameters is poorly identified, and
  the condition number of the Hessian in posterior-SD units; and the
  Laplace log marginal likelihood.

- The VB correction: the largest step, in posterior SDs, at each
  fixed-point iteration, whether Anderson acceleration or the fallback
  optimiser was used (and, after the fallback, the step still to go
  where it landed, and whether it was rerun from there), the largest
  shift, and the term it adds to the log marginal likelihood.

- The marginals: the shape, NMAD, ridge bend, tails, unscanned mass and
  volume correction `gamma1` of each parameter, with the values that
  [`summary.inla_vanilla()`](http://haziqj.ml/inlavanilla/reference/summary.inla_vanilla.md)
  flags marked.

- The calls of `logpost`, `grad` and `hess` in each stage, and in all,
  with the total time.

For the iteration trace of
[`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html) itself, pass
`control = list(nlminb = list(trace = 1))`.

`silent = TRUE` prints nothing, for simulations and other repeated fits:
no progress bars and no stage lines, and every warning and message
raised during the fit, including those of `logpost`, `grad` and `hess`,
is muffled. Errors are still raised. Output printed with
[`cat()`](https://rdrr.io/r/base/cat.html) or
[`print()`](https://rdrr.io/r/base/print.html), by `logpost`, `grad`,
`hess` or the trace of
[`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html), is not
captured.

Whatever the level, the warnings raised during the fit are recorded in
`fit$warnings`, so that a simulation can check them afterwards.

## References

Jamil, H., & Rue, H. (2026a). *Approximate Bayesian inference for
structural equation models using integrated nested Laplace
approximations* (arXiv:2603.25690 \[stat.ME\]). arXiv.
<https://doi.org/10.48550/arXiv.2603.25690>

Jamil, H., & Rue, H. (2026b). *Implementation and workflows for
INLA-based approximate Bayesian structural equation modelling*
(arXiv:2604.00671 \[stat.CO\]). arXiv.
<https://doi.org/10.48550/arXiv.2604.00671>

## See also

[`sampling()`](http://haziqj.ml/inlavanilla/reference/sampling.md) for
posterior draws and derived quantities,
[`summary.inla_vanilla()`](http://haziqj.ml/inlavanilla/reference/summary.inla_vanilla.md)
and
[`diagnostics()`](http://haziqj.ml/inlavanilla/reference/diagnostics.md)
for diagnostics.

## Examples

``` r
# Normal data with unknown mean and standard deviation, on the scale
# (mu, log sigma)
set.seed(1)
y <- rnorm(30, mean = 5, sd = 2)

logpost <- function(theta, y) {
  mu <- theta[["mu"]]
  sigma <- exp(theta[["log_sigma"]])
  sum(dnorm(y, mu, sigma, log = TRUE)) +
    dnorm(mu, 0, 10, log = TRUE) + # prior on mu
    dexp(sigma, rate = 1, log = TRUE) + # prior on sigma ...
    theta[["log_sigma"]] # ... and its Jacobian
}

fit <- inla_vanilla(logpost, init = c(mu = 0, log_sigma = 0), y = y)
#> ✔ Found the posterior mode in 10 iterations [3ms]
#> ✔ VB mean correction, max shift 0.24σ [1ms]
#> ✔ Fitted 2 skew-normal marginals [11ms]
fit
#> ── INLAvanilla ───────────────────────────────────────── Skew-normal + VB ──
#> 
#> Log marginal likelihood: -66.38 · Converged in 10 iterations (19ms)
#> 
#>               Mean     SD   2.5%    50%  97.5%   Mode
#> 1 mu         5.159  0.327  4.517  5.159  5.800  5.159
#> 2 log_sigma  0.616  0.129  0.377  0.610  0.885  0.599

# The posterior of sigma itself, through the draws
sigma <- sampling(fit, fun = function(theta) exp(theta[["log_sigma"]]))
quantile(sigma, c(0.025, 0.5, 0.975))
#>     2.5%      50%    97.5% 
#> 1.444354 1.832695 2.455869 

# A full diagnostic of each stage
fit <- inla_vanilla(
  logpost,
  init = c(mu = 0, log_sigma = 0),
  y = y,
  verbose = TRUE
)
#> ℹ 2 parameters: mu and log_sigma
#>   Derivatives by finite differences of logpost
#>   logpost at init -481.47
#> ✔ Found the posterior mode in 10 iterations [1ms]
#>   nlminb: relative convergence (4), after 10 iterations (14 function and 26
#>     gradient evaluations)
#>   logpost at the mode -65.06, whitened gradient 8.8e-08σ, Newton step
#>     8.8e-08σ (warning > 0.1σ)
#>   Hessian by second differences of logpost: 6 calls, and 12 more to set the
#>     steps
#>   Posterior correlation matrix: largest correlation -0.0042 (mu,
#>     log_sigma), smallest eigenvalue 1
#>   Hessian condition number 1, in SD units (warning > 1e+08)
#>   Laplace log marginal likelihood -66.41
#>   63 logpost calls
#> ✔ VB mean correction, max shift 0.24σ [2ms]
#>   Fixed-point steps, largest over the parameters (σ): 0.25, 0.0079, 2.3e-04
#>   Converged in 3 iterations, Anderson acceleration at 1 of them; no fallback
#>   Largest shift 0.24σ, for log_sigma (warning > 1σ)
#>   KL term 0.029, giving a log marginal likelihood of -66.38
#>   18 logpost calls
#> ✔ Fitted 2 skew-normal marginals [7ms]
#>              Shape     NMAD  Ridge bend  Tails  Unscanned    gamma1
#>   mu         -0.25  4.2e-05        0.13     ok    6.4e-05  -0.00014
#>   log_sigma   1.32   0.0076     5.5e-06     ok    0.00057      0.13
#>   ! marks what summary() flags: NMAD > 0.1, ridge bend > 0.25, heavy tails,
#>     unscanned mass > 0.05
#>   66 logpost calls
#> ℹ Done: 147 logpost calls in 27ms

# Nothing printed, as in a simulation; any warnings are kept in the fit
fit <- inla_vanilla(
  logpost,
  init = c(mu = 0, log_sigma = 0),
  y = y,
  silent = TRUE
)
fit$warnings
#> character(0)
```
