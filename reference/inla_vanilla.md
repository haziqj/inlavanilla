# INLA-style posterior approximation for any log-posterior

Approximates the posterior of a parameter vector from its log density
alone, using the approximation pipeline of INLAvaan with no latent field
and nothing model-specific. The pipeline has four stages:

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

## Usage

``` r
inla_vanilla(
  logpost,
  init,
  grad = NULL,
  ...,
  vb = TRUE,
  verbose = TRUE,
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

- ...:

  Further arguments passed to `logpost` and `grad`, such as the data. As
  in [`stats::optim()`](https://rdrr.io/r/stats/optim.html) and
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html), they come
  after the gradient, and the arguments after them (`vb`, `verbose`,
  `control`) must be given by their full names, so data arguments cannot
  partially match them.

- vb:

  Logical. Apply the VB mean correction? Default `TRUE`. It enters the
  log marginal likelihood, and the skew-normal means are checked against
  it (see
  [`summary.inla_vanilla()`](http://haziqj.ml/inlavanilla/reference/summary.inla_vanilla.md)).

- verbose:

  Logical. Show progress? Default `TRUE`.

- control:

  A list of tuning parameters; see the section on control.

## Value

An object of class `inla_vanilla`, with methods for
[`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) and
[`sampling()`](http://haziqj.ml/inlavanilla/reference/sampling.md). It
is a list holding, among other things:

- `summary`: the posterior summary table (mean, SD, quantiles, mode).

- `mode`, `cov`: the Laplace mode and covariance.

- `mean_vb`: the VB-corrected mean (the mode when `vb = FALSE`).

- `sn`: the skew-normal parameters `xi`, `omega` and `alpha` per
  parameter.

- `logml`: the log marginal likelihood.

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
variance the SD is about 11\\ too small at n = 12, 3\\

## Gradients

`grad` is optional. Without it, derivatives come from finite differences
of `logpost`, and the number of calls grows quickly with the number of
parameters m. The mode and Hessian take about m^2 + 6m calls. The
marginals take about m(`n_grid` + 2m + 5). The VB correction takes about
(m + 1)(2m^2 + 1) calls per sweep, over a few sweeps, so it dominates
beyond a handful of parameters: at m = 20 a fit makes on the order of
50,000 calls. With `grad`, the VB correction takes 2m^2 + 1 gradient
calls per sweep instead, and the fit is faster and more accurate. For
more than about 10 parameters, supply `grad`, or set `vb = FALSE` to
skip the VB correction (and its term in the log marginal likelihood).

A supplied `grad` is checked against finite differences at `init`, and
its Jacobian at the mode is checked for symmetry, with a warning if
either check fails. Any R function works for `logpost` and `grad`,
including compiled code (Rcpp, TMB, RTMB). For a TMB object, which
minimises, pass `function(x) -obj$fn(x)` and
`function(x) -as.vector(obj$gr(x))`, with `init = obj$par`; its repeated
names (one per element of a parameter vector) are numbered, as
`beta[1]`, `beta[2]`, and so on.

## Control

`control` is a list with any of:

- `n_grid`: number of grid points per skew-normal scan, spread over
  \\\pm 4\\ posterior SDs (default `21`).

- `log_threshold`: grid points whose log-posterior falls more than this
  below the maximum are ignored in the skew-normal fit (default `-6`).

- `vb_maxit`, `vb_tol`: iteration limit (default `25`) and convergence
  tolerance in posterior SDs (default `1e-3`) for the VB correction.

- `nlminb`: a list of
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html) controls for
  the mode search. `iter.max` and `eval.max` default to `1000` and
  `2000`.

## References

Jamil, H. and Rue, H. (2026). Approximate Bayesian inference for
structural equation models using integrated nested Laplace
approximations.
[doi:10.48550/arXiv.2603.25690](https://doi.org/10.48550/arXiv.2603.25690)

Jamil, H. and Rue, H. (2026). Implementation and workflows for
INLA-based approximate Bayesian structural equation modelling.
[doi:10.48550/arXiv.2604.00671](https://doi.org/10.48550/arXiv.2604.00671)

## See also

[`sampling()`](http://haziqj.ml/inlavanilla/reference/sampling.md) for
posterior draws and derived quantities,
[`summary.inla_vanilla()`](http://haziqj.ml/inlavanilla/reference/summary.inla_vanilla.md)
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
#> ✔ Found the posterior mode in 10 iterations [5ms]
#> ✔ VB mean correction, max shift 0.24σ [2ms]
#> ✔ Fitted 2 skew-normal marginals [7ms]
fit
#> ── INLAvanilla ───────────────────────────────────────────────────────── 20ms ──
#> 2 parameters · skew-normal marginals · VB-corrected
#> 
#>             Mean     SD   2.5%    50%  97.5%   Mode
#> mu         5.159  0.327  4.517  5.159  5.800  5.159
#> log_sigma  0.616  0.129  0.377  0.610  0.885  0.599
#> 
#> Log marginal likelihood: -66.38

# The posterior of sigma itself, through the draws
sigma <- sampling(fit, fun = function(theta) exp(theta[["log_sigma"]]))
quantile(sigma, c(0.025, 0.5, 0.975))
#>     2.5%      50%    97.5% 
#> 1.444354 1.832695 2.455869 
```
