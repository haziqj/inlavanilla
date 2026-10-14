# Leave-one-out cross-validation and WAIC for an INLAvanilla fit

`loo()` estimates how well the model predicts new data by leave-one-out
cross-validation (LOO), and `waic()` by the widely applicable
information criterion. Both work from the fit alone, with no refitting
and no posterior draws, by a second-order Taylor expansion of each
unit's log-likelihood (Alhyari et al., 2026), as in INLAvaan.

Both need the log-likelihood split into one term per observation, or per
independent unit, which
[`inla_vanilla()`](http://haziqj.ml/inlavanilla/reference/inla_vanilla.md)
never sees. You supply it as `loglik`: write the log-likelihood terms
and the log-prior as separate functions, and `logpost` as their sum.

## Usage

``` r
loo(object, ...)

waic(object, ...)

# Default S3 method
loo(object, ...)

# Default S3 method
waic(object, ...)

# S3 method for class 'inla_vanilla'
loo(object, loglik, ...)

# S3 method for class 'inla_vanilla'
waic(object, loglik, ...)
```

## Arguments

- object:

  An object of class `inla_vanilla`.

- ...:

  Further arguments passed to `loglik`, such as the data. As in
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html), a name that
  abbreviates `object` or `loglik` (such as `l`) would be taken as that
  argument, so give such data a longer name. If the loo package is
  attached after inlavanilla, its generics `loo(x, ...)` and
  `waic(x, ...)` are the ones called, and they take an argument named
  `x` as the fit: then do not name data `x`, or call
  `inlavanilla::loo()`.

- loglik:

  A function `loglik(theta, ...)` returning the log-likelihood terms,
  one per observation or independent unit, at the named parameter vector
  `theta` (on the scale of `logpost`). Their sum plus the log-prior is
  `logpost`. Names, if any, label the units.

## Value

`loo()` returns an object of class `inla_vanilla_loo`, which also has
class `loo`: a list with

- `estimates`: a matrix with rows `elpd_loo`, `p_loo` and `looic` and
  columns `Estimate` and `SE`;

- `pointwise`: a matrix with one row per unit and columns `elpd_loo`
  (the log CPO), `p_loo`, `looic` and `k` (the unit's largest leverage);

- `order`: the order of the expansion used, 2 or 1;

- `first_order_units`: the units whose second-order term does not exist,
  if any;

- `theta`: the point of the expansion, and `calls`: the number of calls
  of `loglik`;

- `warnings`: the messages of any warnings raised.

`waic()` returns an object of class `inla_vanilla_waic` (also `waic` and
`loo`) with the same layout, with rows and columns `elpd_waic`, `p_waic`
and `waic`.

## Details

Write \\\ell_u(\vartheta) = \log p(y_u \mid \vartheta)\\ for the
log-likelihood of unit \\u\\. The leave-one-out predictive density of
unit \\u\\, its conditional predictive ordinate, is \$\$\mathrm{CPO}\_u
= p(y_u \mid y\_{-u}) = 1 / E\[1 / p(y_u \mid \vartheta) \mid y\],\$\$
an expectation over the full posterior. The posterior is taken as the
Gaussian \\N(\tilde\vartheta, \Omega)\\ of the fit, with the
VB-corrected mean (the mode when `vb = FALSE`) and the Laplace
covariance, and \\\ell_u\\ is expanded to second order about
\\\tilde\vartheta\\, with score \\s_u\\ and Hessian \\H_u\\. The
expectation is then a Gaussian integral: \$\$\log \mathrm{CPO}\_u =
\ell_u - \tfrac12 s_u^\top (\Omega^{-1} + H_u)^{-1} s_u + \tfrac12 \log
\|I + \Omega H_u\|.\$\$ The same expansion gives the log pointwise
predictive density \$\$\mathrm{lpd}\_u = \log E\[p(y_u \mid \vartheta)
\mid y\] = \ell_u + \tfrac12 s_u^\top (\Omega^{-1} - H_u)^{-1} s_u -
\tfrac12 \log \|I - \Omega H_u\|\$\$ and the WAIC penalty
\\p\_{\mathrm{waic},u} = \mathrm{Var}\[\ell_u(\vartheta)\] = s_u^\top
\Omega s_u + \tfrac12 \mathrm{tr}\[(\Omega H_u)^2\]\\. Then
\\\mathrm{elpd}\_{\mathrm{loo}} = \sum_u \log \mathrm{CPO}\_u\\,
\\p\_{\mathrm{loo}} = \sum_u (\mathrm{lpd}\_u - \log \mathrm{CPO}\_u)\\,
\\\mathrm{elpd}\_{\mathrm{waic}} = \sum_u (\mathrm{lpd}\_u -
p\_{\mathrm{waic},u})\\, and `looic` and `waic` are \\-2\\ times the
elpd. Each standard error is \\\sqrt{n \\ \mathrm{Var}}\\ of the
pointwise terms over the \\n\\ units, as in the loo package.

The scores and Hessians come from finite differences of `loglik` along
the principal axes of \\\Omega\\, a hundredth of a posterior SD apart:
\\m^2 + m + 1\\ calls of `loglik` for \\m\\ parameters, each returning
every unit.

**Existence and the leverage k.** For each unit, the eigenvalues \\k\\
of \\-\Omega H_u\\ are the shares of the posterior precision the unit
carries along its directions. The second-order \\\log \mathrm{CPO}\_u\\
exists only when every \\k \< 1\\: at \\k \ge 1\\, the rest of the data
and the prior cannot identify some combination of the parameters without
the unit, and \\E\[1 / p(y_u \mid \vartheta)\]\\ diverges. If any unit
fails, every LOO estimate falls back to the first-order \\\log
\mathrm{CPO}\_u = \ell_u - \tfrac12 s_u^\top \Omega s_u\\, with a
warning naming the units, rather than mixing orders or dropping units.
The largest \\k\\ of each unit is reported as `k`, a measure of its
influence. Likewise the second-order \\\mathrm{lpd}\_u\\ needs every \\k
\> -1\\. Where that fails, WAIC falls back to first order, at which it
equals the first-order LOO, and `p_loo` uses the first-order term for
that unit.

The expansion rests on the Gaussian approximation of the posterior, as
the VB correction does. A cross-check that uses the skew-normal
marginals is Pareto-smoothed importance sampling on copula draws:
`loo::loo(sampling(fit, 4000, fun = loglik, ...))`.

## Silent fits

On a fit made with `silent = TRUE`, `loo()` and `waic()` print nothing
and raise no warnings, their own or those of `loglik`; the warnings are
kept in the result's `warnings`, and a fallback to first order shows in
`order` and `first_order_units`.

## Other objects

For any other object, `loo()` and `waic()` call the loo package's
functions of the same name, if it is installed, so attaching inlavanilla
after loo (or a package that uses it, such as rstan or brms) does not
get in their way.

## Comparing models

The results have the layout of the loo package's objects, so
[`loo::loo_compare()`](https://mc-stan.org/loo/reference/loo_compare.html)
compares fits to the same units. By hand, with `d` the difference of two
fits' `pointwise[, "elpd_loo"]`, the difference in elpd is `sum(d)` with
standard error `sqrt(length(d)) * sd(d)`.

## References

Alhyari, M., Jamil, H., Montcho, H., & Rue, H. (2026). *Deterministic
leave-one-cluster-out cross-validation for multilevel Bayesian models*
(arXiv:2609.00670 \[stat.ME\]). arXiv.
<https://doi.org/10.48550/arXiv.2609.00670>

Vehtari, A., Gelman, A., & Gabry, J. (2017). Practical Bayesian model
evaluation using leave-one-out cross-validation and WAIC. *Statistics
and Computing*, *27*(5), 1413–1432.
<https://doi.org/10.1007/s11222-016-9696-4>

Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation
and widely applicable information criterion in singular learning theory.
*Journal of Machine Learning Research*, *11*, 3571–3594.
<https://jmlr.org/papers/v11/watanabe10a.html>

## Examples

``` r
# Poisson counts: is a log-linear trend better than a constant rate?
y <- c(2, 3, 6, 7, 8, 9, 10, 12, 15)
time <- seq_along(y) - 5

loglik <- function(theta, time, y) {
  dpois(y, exp(theta[["a"]] + theta[["b"]] * time), log = TRUE)
}
logprior <- function(theta) sum(dnorm(theta, 0, 5, log = TRUE))
logpost <- function(theta, time, y) {
  sum(loglik(theta, time, y)) + logprior(theta)
}

fit <- inla_vanilla(logpost, c(a = 0, b = 0), time = time, y = y)
#> ✔ Found the posterior mode in 9 iterations [1ms]
#> ✔ VB mean correction, max shift 0.12σ [2ms]
#> ✔ Fitted 2 skew-normal marginals [6ms]
loo(fit, loglik, time = time, y = y)
#> ── Leave-one-out cross-validation ─────────────────────────────────── 9 units ──
#> 
#>           Estimate   SE
#> elpd_loo     -19.2  0.8
#> p_loo          0.6  0.2
#> looic         38.3  1.6
#> 
#> ✔ Second-order expansion for every unit. Largest leverage k = 0.49, for unit 9.
waic(fit, loglik, time = time, y = y)
#> ── WAIC ───────────────────────────────────────────────────────────── 9 units ──
#> 
#>            Estimate   SE
#> elpd_waic     -19.1  0.8
#> p_waic          0.6  0.1
#> waic           38.3  1.6
#> 
#> ✔ Second-order expansion for every unit.
```
