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

Stages 2 and 3 are optional and independent of each other. `vb = "off"`
skips the VB correction, which leaves the marginals unchanged (they are
scanned around the mode) but gives the Laplace log marginal likelihood.
`marginals = "gaussian"` skips the scans and uses the Gaussian's own
marginals \\N(\tilde\vartheta_j, \Omega\_{jj})\\, centred at the VB mean
(or at the mode when `vb = "off"`);
[`sampling()`](http://haziqj.ml/inlavanilla/reference/sampling.md) then
draws from the joint Gaussian. So `vb = "off", marginals = "gaussian"`
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
  vb = c("onestep", "full", "off"),
  marginals = c("skewnormal", "gaussian"),
  width = c("trace", "cmp", "tk"),
  tk_fallback = FALSE,
  cores = NULL,
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
  them (`vb`, `marginals`, `width`, `tk_fallback`, `cores`, `verbose`,
  `silent`, `control`) must be given by their full names, so data
  arguments cannot partially match them. A data argument cannot take one
  of those names. The arguments before them can be matched, as in
  [`nlminb()`](https://rdrr.io/r/stats/nlminb.html): a data argument
  named `g` or `gr` would be taken as `grad`, and one named `h` or `he`
  as `hess`, so give it a longer name.

- vb:

  The VB mean correction: `"onestep"` (the default) takes one Newton
  step from the mode, `"full"` iterates to the VB solution, and `"off"`
  (or `FALSE`) leaves it out. `TRUE` means the default, and
  abbreviations such as `"one"` work too. See the section on the VB
  correction. The correction enters the log marginal likelihood, and the
  skew-normal means are checked against it (see
  [`summary.inla_vanilla()`](http://haziqj.ml/inlavanilla/reference/summary.inla_vanilla.md)).

- marginals:

  The marginal of each parameter: `"skewnormal"` (the default) fits a
  skew-normal to a scan of the log-posterior; `"gaussian"` uses the
  Gaussian marginal of the Laplace approximation, centred at the VB mean
  unless `vb = "off"`, with no scans. Abbreviations such as `"skew"` and
  `"gauss"` work too.

- width:

  How each skew-normal scan accounts for the changing volume of the
  other parameters, the width term \\\gamma_j(t) = -\frac12
  \log\|H\_{-j}\|\\ (see the section on the width term). With m
  parameters and `n_grid` scan points, per parameter:

  - `"trace"` (the default): its linear term \\\gamma_j' t\\, as in
    INLAvaan, from 4 `grad` calls, and one more for all the parameters
    together (2m + 5 `logpost` calls without `grad`), on top of the
    `n_grid - 1` `logpost` calls of the scan.

  - `"cmp"`: the exact width term along the same straight scan, from the
    Hessian of the other parameters at each scan point: about
    `n_grid - 1` times 2(m - 1) `grad` calls (`n_grid - 1` `hess` calls
    with `hess`, or `n_grid - 1` times (m - 1)m `logpost` calls without
    `grad`).

  - `"tk"`: the Laplace marginal of Tierney and Kadane (1986). At each
    scan point the other parameters are moved to their conditional mode
    by [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html), started
    from the straight scan, and the height and the width term are taken
    there: the calls of `"cmp"` plus those of `n_grid - 1` conditional
    optimisations, typically 5 to 20 `logpost` and `grad` calls each
    with `grad`, and many more without.

  With one parameter the three are the same. Abbreviations work. Ignored
  for Gaussian marginals.

- tk_fallback:

  Logical. Rescan with `width = "tk"` every parameter whose `ridge_bend`
  is flagged (above 0.25, see
  [`summary.inla_vanilla()`](http://haziqj.ml/inlavanilla/reference/summary.inla_vanilla.md))?
  Default `FALSE`, so that a fit costs what `width` says and gives the
  same marginals whatever its diagnostics. The parameters rescanned are
  listed in `fit$fallback` and by
  [`summary.inla_vanilla()`](http://haziqj.ml/inlavanilla/reference/summary.inla_vanilla.md).
  Ignored for Gaussian marginals and with `width = "tk"`.

- cores:

  The number of CPU cores for the skew-normal scans, and for
  [`loo()`](http://haziqj.ml/inlavanilla/reference/loo.md) and
  [`waic()`](http://haziqj.ml/inlavanilla/reference/loo.md) on this fit.
  `NULL`, the default, follows
  [`future::plan()`](https://future.futureverse.org/reference/plan.html),
  and is serial if no plan is set. `1` is serial whatever the plan.
  Needs the future package (see the section on parallel computation).

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
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`plot()`](https://rdrr.io/r/graphics/plot.default.html),
[`sampling()`](http://haziqj.ml/inlavanilla/reference/sampling.md),
[`diagnostics()`](http://haziqj.ml/inlavanilla/reference/diagnostics.md)
and [`timing()`](http://haziqj.ml/inlavanilla/reference/timing.md). It
is a list holding, among other things:

- `summary`: the posterior summary table (mean, SD, quantiles, mode).

- `mode`, `cov`: the Laplace mode and covariance.

- `mean_vb`: the VB-corrected mean (the mode when `vb = "off"`).

- `sn`: the skew-normal parameters `xi`, `omega` and `alpha` per
  parameter (with `alpha = 0` for Gaussian marginals).

- `vb`, `marginals`, `width`, `tk_fallback`, `cores`, `verbose`,
  `silent`: the settings the fit used (`width` is `NA` for Gaussian
  marginals).

- `fallback`: the parameters that `tk_fallback` rescanned with
  `width = "tk"` (empty when there were none).

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
  correction, whether it stopped at the first step (`one_pass`), whether
  it fell back to
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html) (`fallback`),
  and the Kullback-Leibler divergence of the shift (`kld_global`); it is
  `NULL` when `vb = "off"`. `marginals` is a data frame with a row per
  parameter: `mode_shift_sigma`, `kld`, `vb_shift`, `vb_shift_sigma`,
  `nmad`, `alpha`, `unscanned_mass`, `vb_step_sigma`, `grad_gap`,
  `mean_gap`, `gamma1`, `tails` and `ridge_bend`, with `NA` for the
  checks the fit did not run, and `width`, the width method of each scan
  (see the section on the width term).

- `scan`: the scans the skew-normals were fitted to, for
  `plot(fit, type = "sn_fit")`: the grid `z` in posterior SDs from the
  mode, the log-posterior at each point (`logpost`, a column per
  parameter: along the CMP, or at the conditional modes for
  `width = "tk"`), the width term added to it (`width`, zero at the
  mode), and the log scale factor `logC` of each fit on that grid. It is
  `NULL` for Gaussian marginals.

- `workers`: the number of parallel jobs the scans ran on, 1 for a
  serial fit (see the section on parallel computation).

- `timing`: the elapsed time of each stage of the fit, in seconds, which
  [`timing()`](http://haziqj.ml/inlavanilla/reference/timing.md)
  returns.

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
flags this narrowing where it predicts about 10% or more. By default
(`width = "trace"`) this is the scan of INLAvaan; `width = "tk"`, or
`tk_fallback = TRUE` for the flagged parameters only, removes the
narrowing at a higher cost (see the section on the width term).

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

## The width term

Each skew-normal scan follows the conditional mode path (CMP) of the
Laplace Gaussian, \\\vartheta(t) = \vartheta^\* + t v_j\\ with \\v_j =
\Omega\_{\cdot j} / \sqrt{\Omega\_{jj}}\\, so that \\t\\ counts
posterior SDs. The log marginal is then about the height of the
log-posterior along the path plus the width term \\\gamma_j(t) =
-\frac12 \log\|H\_{-j}(t)\|\\, where \\H\_{-j}\\ is the negative Hessian
of the log-posterior with row and column j deleted: the log-volume of
the other parameters. `width` sets how it is taken, in three levels that
each remove one error:

- `"trace"`, the default and INLAvaan's, linearises it, \\\gamma_j(t)
  \approx \gamma_j' t\\, with the slope from whitened differences of the
  gradient (Jamil & Rue, 2026a, Corollary 3.3).

- `"cmp"` takes \\\gamma_j(t)\\ exactly along the same path. Its
  difference from `"trace"` is the error of the linearisation.

- `"tk"` moves the other parameters to their conditional mode at each
  point of the scan and takes the height and the width term there: the
  Laplace approximation of the marginal (Tierney & Kadane, 1986), INLA's
  "Laplace strategy". Its difference from `"cmp"` is the error of the
  straight path, which the `ridge_bend` diagnostic of
  [`summary.inla_vanilla()`](http://haziqj.ml/inlavanilla/reference/summary.inla_vanilla.md)
  predicts.

The skew-normal is fitted to the corrected heights in the same way for
all three, on the same grid. `tk_fallback = TRUE` keeps the cheap scan
for every parameter but those whose ridge bend is flagged, and rescans
those with `"tk"`. In `fit$diagnostics$marginals`, `width` gives the
method each parameter used and `gamma1` the slope of its width term at
the mode (for `"cmp"` and `"tk"`, a central difference across the mode).
The ridge bend of a parameter scanned with `"tk"` is still reported, but
not flagged.

## Gradients and Hessians

`grad` is optional. Without it, derivatives come from finite differences
of `logpost`, and the number of calls grows with the square of the
number of parameters m. The mode and Hessian take about m^2 + 6m calls,
and the marginals about m(`n_grid` + 2m + 6), with 2(m - 1) more for
each parameter whose ridge bend is above 0.1. The VB correction takes
(m + 1)(m + 2) calls, or that many per sweep with `vb = "full"`: it
needs no derivatives, because Stein's identity gives the expected
gradient from values of `logpost` alone. With `vb = "full"`, the VB mean
is then exact when the log-posterior is a polynomial of degree four in
whitened coordinates, and up to degree six with `grad`. On test models
it was about 0.01 posterior SDs from exact VB (0.003 with `grad`), and
up to 0.1 on a strongly skewed posterior (0.03 with `grad`), against a
VB error of 0.01 to 0.2 SDs from the posterior mean itself. On the most
skewed posteriors the VB correction may then not settle, and the fit
says it has not converged (see
[`summary.inla_vanilla()`](http://haziqj.ml/inlavanilla/reference/summary.inla_vanilla.md));
`grad` is the remedy.

With `grad`, every stage is more accurate. The fit is also faster when a
gradient call costs about as much as a `logpost` call, as a hand-written
gradient usually does, since it then makes far fewer calls in all. A
gradient by automatic differentiation (RTMB, TMB) costs about two, and
beyond a handful of parameters the fit is still several times faster
than without one (four times at m = 39). The default VB correction,
`vb = "onestep"`, takes 2m + 1 gradient calls, about as many as one
numerical Hessian, and `vb = "full"` takes at most (m + 1)(m + 2) + 1
per sweep, usually about four sweeps. Beyond a handful of parameters the
marginals are the slowest stage, along with the VB correction under
`vb = "full"`. `vb = "off"` skips the VB correction (and its term in the
log marginal likelihood), and the skew-normal marginals do not depend on
it; `marginals = "gaussian"` saves the calls of the scans.

`hess` is optional too, and needs `grad`, as in
[`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html). With it,
[`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html) takes Newton
steps to the mode, and the Hessian at the mode is the supplied one. The
VB correction does not use it, nor do the default marginals, whose
volume correction takes 4m + 1 `grad` calls for all the parameters
together, where a Hessian each would cost far more. `width = "cmp"` and
`width = "tk"` do use it. The fit is the one `grad` alone gives, to the
error of its finite differences (about 1e-8, relative), except that the
mode is found more precisely. A dense Hessian costs about m gradients,
so this takes a little longer: a structural equation model with 39
parameters and an automatic-differentiation Hessian (RTMB) took 0.18 s
with it and 0.17 s without.

The fit goes straight to the mode: at `init` it only checks that
`logpost`, `grad` and `hess` return finite values of the right shape. A
supplied `grad` is compared with `logpost` at the mode, from differences
the fit takes there anyway, so at no cost in calls, and its Jacobian
there is checked for symmetry, with a warning if either check fails (see
`grad_gap` in
[`summary.inla_vanilla()`](http://haziqj.ml/inlavanilla/reference/summary.inla_vanilla.md)).
A supplied `hess` is checked at the mode against central differences of
`grad`, along two directions, with a warning if it does not match;
`grad` is then not compared with `logpost`, since no differences of
`logpost` are taken. Any R function works for `logpost`, `grad` and
`hess`, including compiled code (Rcpp, TMB) and RTMB, which
differentiates R code. An RTMB function can return the log-posterior
itself, and its `obj$fn`, `obj$gr` and `obj$he` are then passed as they
are. For a TMB template that returns the negative log-posterior, pass
`function(x) -obj$fn(x)`, `function(x) -as.vector(obj$gr(x))` and
`function(x) -obj$he(x)`, with `init = obj$par`; its repeated names (one
per element of a parameter vector) are numbered, as `beta[1]`,
`beta[2]`, and so on. TMB gives no Hessian for a model with random
effects.

## The VB correction

The VB mean correction moves the centre of the Laplace Gaussian towards
the posterior, keeping its covariance. It maximises the expected
log-posterior under that Gaussian, so that at the solution the gradient,
averaged over the Gaussian, is zero.

`vb = "onestep"` takes the first Newton step towards that solution, with
the Laplace curvature, and stops. With a gradient, the average is taken
over the mode and \\\pm\sqrt 3\\ posterior SDs along each principal axis
of the covariance. Without one, it is one pass of Stein's identity on
the rule of `vb = "full"`. Most of the movement is in that step, and
exact VB is itself an approximation to the posterior mean with an error
of the same order. If the first step is over 2.5 posterior SDs, or a
node fails, the iteration carries on with capped steps, within the
number of calls that the mode and Hessian took.
[`diagnostics()`](http://haziqj.ml/inlavanilla/reference/diagnostics.md)
then has `vb_step_sigma`, the step still to go, which is `NA` after a
single step. After a single step the fit checks that it delivered what
it promised: `vb_gain_ratio` in
[`diagnostics()`](http://haziqj.ml/inlavanilla/reference/diagnostics.md)
is the rise in the expected log-posterior over the rise the step
assumed, and a ratio below 0.5 raises a warning. A large first step can
also overstate the VB term of the log marginal likelihood, so after a
shift above 1 posterior SD, which is flagged, compare models with care:
`vb = "full"` may be closer.

`vb = "full"` iterates to the solution, on a degree-5 rule with at most
(m + 1)(m + 2) + 1 nodes, and hands over to
[`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html) if it does not
converge.

## Parallel computation

Set `cores` to run the per-parameter scans of the skew-normal marginals
in parallel:

    cores <- future::availableCores(omit = 1)
    fit <- inla_vanilla(logpost, init, width = "tk", cores = cores)

`availableCores()` counts the cores the machine offers, and `omit = 1`
leaves one free. A fit never asks for more workers than the machine
offers. The results and the progress line are those of a serial fit. A
long scan uses all the `cores`, a shorter one fewer, and a scan expected
to take under a second stays serial. The mode, the Hessian and the VB
correction stay serial, so the gain is largest with `width = "tk"` or
`"cmp"`, where the scans are nearly the whole fit: about 7.5 times
faster on 12 cores on test models. With the default width the scans are
about half of a slow fit, and it was about 1.5 times at most.

`cores` needs the future package, which is optional. It uses forked
workers where those are safe, in a terminal on macOS or Linux, and
separate R sessions elsewhere: on Windows, in RStudio, and in Positron's
console, which refuses to fork. Forks share everything. Separate
sessions need the package installed, not just loaded with
`devtools::load_all()`, and are sent your functions and data. Compiled
code (Rcpp, TMB, RTMB) cannot be sent, so the fit first checks that a
worker gets the same value at the mode, and otherwise runs serially with
a message.

To use separate sessions with compiled code, rebuild it on each worker:
put whatever creates it in a function, and call that function again when
the process is not the one that started the fit. With RTMB and random
effects (without them, RTMB needs none of this):

    build <- function() {
      RTMB::MakeADFun(model, parameters, random = "u", silent = TRUE)
    }
    obj <- build()
    pid <- Sys.getpid()
    rebuild <- function() {
      if (Sys.getpid() != pid) {
        obj <<- build()
        pid <<- Sys.getpid()
      }
    }
    logpost <- function(theta) {
      rebuild()
      -obj$fn(theta)
    }
    grad <- function(theta) {
      rebuild()
      -obj$gr(theta)
    }

With Rcpp, `rebuild()` calls `Rcpp::sourceCpp(file, cacheDir = dir)`
instead, with a `dir` of your own so that the workers share one compiled
copy. The RTMB and Rcpp articles on the package website have worked
versions.

With `cores` left at `NULL` the fit follows a plan set with
[`future::plan()`](https://future.futureverse.org/reference/plan.html),
and runs serially if there is none, so a plan can apply to every fit in
a session. `cores = 1` runs serially whatever the plan.
[`loo()`](http://haziqj.ml/inlavanilla/reference/loo.md) and
[`waic()`](http://haziqj.ml/inlavanilla/reference/loo.md) use the
setting of the fit.

## Control

`control` is a list with any of:

- `n_grid`: number of grid points per skew-normal scan, spread over
  \\\pm 4\\ posterior SDs (default `21`).

- `log_threshold`: grid points whose log-posterior falls more than this
  below the maximum are ignored in the skew-normal fit (default `-6`).

- `vb_maxit`, `vb_tol`: iteration limit (default `25`) and convergence
  tolerance in posterior SDs (default `1e-3`) for the VB correction's
  fixed-point iteration, with `vb = "full"` and when a one-step
  correction has to iterate. A step over 1.5 posterior SDs is cut to
  1.5. With `vb = "full"`, if two steps in a row are cut, or it has not
  converged after `vb_maxit` steps,
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html) takes over.

- `nlminb`: a list of
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html) controls for
  the mode search. `iter.max` and `eval.max` both default to `5000`, and
  `trace = 1` prints its iterations.

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
  from, and `logpost` at `init`.

- The mode: the message, iterations and evaluations of
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html), and whether
  it was restarted with the parameters rescaled; the log-posterior at
  the mode, the whitened gradient there (which decides the restart) and
  the largest Newton step in posterior SDs (which the fit-time warning
  checks); how the Hessian was taken (and, for a supplied `grad`, the
  asymmetry of its Jacobian and its gap from `logpost`, or, for a
  supplied `hess`, how far it is from central differences of `grad`);
  the largest posterior correlation and the smallest eigenvalue of the
  correlation matrix, which is near zero when some combination of the
  parameters is poorly identified, and the condition number of the
  Hessian in posterior-SD units; and the Laplace log marginal
  likelihood.

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

Tierney, L., & Kadane, J. B. (1986). Accurate approximations for
posterior moments and marginal densities. *Journal of the American
Statistical Association*, *81*(393), 82–86.
<https://doi.org/10.1080/01621459.1986.10478240>

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
#> ✔ Found the posterior mode in 10 iterations [2ms]
#> ✔ VB mean correction, max shift 0.25σ [2ms]
#> ✔ Fitted 2 skew-normal marginals [5ms]
fit
#> ── INLAvanilla ──────────────────────────────── Skew-normal + one-step VB ──
#> 
#> Log marginal likelihood: -66.38 · Converged in 10 iterations (12ms)
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
#> ✔ Found the posterior mode in 10 iterations [2ms]
#>   nlminb: relative convergence (4), after 10 iterations (14 function and 26
#>     gradient evaluations)
#>   logpost at the mode -65.06, whitened gradient 9e-08σ, Newton step 9e-08σ
#>     (warning > 0.1σ)
#>   Hessian by second differences of logpost: 6 calls, and 12 more to set the
#>     steps
#>   Posterior correlation matrix: largest correlation -0.0042 (mu,
#>     log_sigma), smallest eigenvalue 1
#>   Hessian condition number 1, in SD units (warning > 1e+08)
#>   Laplace log marginal likelihood -66.41
#>   63 logpost calls
#> ✔ VB mean correction, max shift 0.25σ [2ms]
#>   Fixed-point steps, largest over the parameters (σ): 0.25
#>   One step, accepted (guard > 2.5σ)
#>   Gain ratio 0.95 (warning < 0.5)
#>   Largest shift 0.25σ, for log_sigma (warning > 1σ)
#>   KL term 0.031, giving a log marginal likelihood of -66.38
#>   16 logpost calls
#> ✔ Fitted 2 skew-normal marginals [6ms]
#>              Shape     NMAD  Ridge bend  Tails  Unscanned    gamma1
#>   mu         -0.25  4.2e-05        0.13     ok    6.4e-05  -0.00014
#>   log_sigma   1.32   0.0076     3.6e-06     ok    0.00057      0.13
#>   ! marks what summary() flags: NMAD > 0.1, ridge bend > 0.25, heavy tails,
#>     unscanned mass > 0.05
#>   68 logpost calls
#> ℹ Done: 147 logpost calls in 25ms

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
