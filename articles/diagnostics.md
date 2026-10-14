# Diagnostics

## What can go wrong

[`inla_vanilla()`](http://haziqj.ml/inlavanilla/reference/inla_vanilla.md)
approximates a posterior in stages, and each stage can fail in its own
way. The mode search can stop short. The Gaussian can be a poor centre,
so the VB correction moves far. A skew-normal can fail to match the
shape of a marginal. And each marginal is fitted to a scan of only \\\pm
4\\ posterior SDs around the mode, so whatever lies beyond is
extrapolated. A posterior SD in these checks is the Laplace SD at the
mode, \\\sqrt{\Sigma\_{jj}}\\, as in INLAvaan; it is not the SD column
of the summary table, which is the skew-normal’s.

Each stage is checked as it runs, at almost no cost: the checks reuse
what the fit has already computed.
[`summary()`](https://rdrr.io/r/base/summary.html) reports every check
and flags the values that look off, and
[`diagnostics()`](http://haziqj.ml/inlavanilla/reference/diagnostics.md)
returns them as numbers. Most problems also raise one warning at the end
of the fit. Printing the fit opens with the log marginal likelihood and
whether the fit converged, and ends with a line counting its warnings
and one counting the flagged checks, each shown only when there are any,
so even a fit run with `silent = TRUE` shows that something needs a
look.

Two of the checks say whether the computation reached its solution: the
mode check and the VB check. When both pass, the fit has converged. The
others judge how well the approximation fits the posterior, which a
converged fit can still do badly.

| Check | In `diagnostics(fit, type = "param")` | Flagged by [`summary()`](https://rdrr.io/r/base/summary.html) | Fit-time warning | From |
|----|----|----|----|----|
| Mode | `mode_shift_sigma` | above 0.1 SDs | the same | INLAvaan |
| VB settled | `vb_step_sigma` | above 0.1 SDs | the same | inlavanilla |
| Hessian | `hess_cond` (in `diagnostics(fit)`) | above \\10^8\\ | the same | INLAvaan |
| VB correction | `vb_shift_sigma` | above 1 SD | the same | INLAvaan |
| Skew-normal fit | `nmad` | above 0.1 | the same | INLAvaan |
| Unscanned mass | `unscanned_mass` | above 0.05 | the same | INLAvaan |
| Tails | `tails` | heavy | the same | inlavanilla |
| Ridge bend | `ridge_bend` | above 0.25 nats | no | inlavanilla |
| Mean gap | `mean_gap` | above 0.1 SDs | no | inlavanilla |

Most of these are the checks of INLAvaan’s
[`diagnostics()`](http://haziqj.ml/inlavanilla/reference/diagnostics.md),
under the same names. Jamil and Rue (2026) use the VB shift and NMAD as
fit diagnostics, and the gap between the skew-normal and VB means as an
internal check; the mode shift, the Hessian condition number and the
unscanned mass come from INLAvaan’s
[`diagnostics()`](http://haziqj.ml/inlavanilla/reference/diagnostics.md),
not from the paper. The levels of those five checks are INLAvaan’s, and
the fit-time warning uses its words. The VB settling and tails checks
are inlavanilla’s own, and so are the ridge bend and mean gap checks,
which only [`summary()`](https://rdrr.io/r/base/summary.html) flags.
Every level is a package default, chosen so that a healthy fit stays
silent, and not a calibrated cut-off: a tripped check is a prompt to
look again, and silence is not a certificate of accuracy. In the
simulations of Jamil and Rue (2026, Section 4.4), the VB shift and NMAD
track the loss of accuracy, and the paper reports a calibration of the
expected error against them in its Supplementary Material.

## A healthy fit

Poisson counts with a Gamma(2, 1) prior on the rate, on the log scale:

``` r

library(inlavanilla)

y <- c(2, 0, 3, 1, 4, 2, 1)
logpost_pois <- function(theta, y) {
  rate <- exp(theta[["log_rate"]])
  sum(dpois(y, rate, log = TRUE)) +
    dgamma(rate, 2, 1, log = TRUE) +
    theta[["log_rate"]] # the Jacobian
}
fit_pois <- inla_vanilla(logpost_pois, c(log_rate = 0), y = y, silent = TRUE)
summary(fit_pois)
#> ── INLAvanilla summary ───────────────────────────────── Skew-normal + VB ──
#> 
#> Log marginal likelihood: -12.354 · Converged in 7 iterations (18ms)
#> Laplace -12.362 + VB term 0.0082
#> 66 logpost calls: mode 6ms, VB 2ms, marginals 5ms
#> 
#> ── Posterior marginals ──
#> 
#>              Mean     SD   2.5%    50%  97.5%   Mode
#> 1 log_rate  0.595  0.263  0.047  0.607  1.078  0.631
#> 
#> ── Diagnostics ──
#> 
#>             VB shift    NMAD  Unscanned  Ridge bend  Mean gap
#> 1 log_rate     -0.13  0.0066    0.00038           0   -0.0016
#> 
#> ✔ Mode found in 7 iterations (Newton step 9.8e-08σ).
#> ✔ VB correction settled in 3 iterations: max shift 0.13σ.
#> ✔ Skew-normal fits: max NMAD 0.0066, unscanned mass 0.00038, ridge bend 0,
#>   mean gap 0.0016σ.
```

Every check passes. The marginal is mildly skewed (its skew-normal shape
`alpha`, in `diagnostics(fit_pois, type = "param")`, is about -1.4), the
skew-normal matches the scan to within 1% of its peak (NMAD), and its
mean agrees with the VB-corrected mean.

## The mode, the Hessian and the VB correction

The mode is checked by the Newton step there, in posterior SDs
(`mode_shift_sigma`): how far one step of Newton’s method would still
move each parameter. It is flagged above 0.1, INLAvaan’s level, whatever
the optimiser’s return code, because without a gradient the optimiser
can report failure at an accurate mode, and success short of one. When
the mode is off, the warning says so first, and suggests what to change.

`hess_cond` is the condition number of the Hessian at the mode, with
each parameter in posterior-SD units: the ratio of the largest to the
smallest eigenvalue of the posterior correlation matrix. Above \\10^8\\
some combination of the parameters is weakly identified. INLAvaan takes
the condition number on the raw scale of the parameters, where a mere
change of units can make it large; inlavanilla allows posterior SDs from
\\10^{-6}\\ to \\10^4\\, so it measures it in posterior SDs instead.

The VB correction recentres the Gaussian at the shift that minimises its
Kullback-Leibler divergence from the posterior (Jamil & Rue, 2026,
Section 3.1). `vb_shift_sigma` is how far it moved each mean from the
mode, in posterior SDs. `vb_shift` is the same shift on the scale of
`logpost`, and `kld`, half the square of `vb_shift_sigma`, is the
Kullback-Leibler divergence of each shifted marginal from the unshifted
one. A shift of a few tenths is normal for a skewed posterior; a shift
above 1 SD means that the posterior is far from the Laplace Gaussian.
When the fixed-point iteration hands over to the optimiser,
`fit$diagnostics$vb` records it, and
[`summary()`](https://rdrr.io/r/base/summary.html) and the verbose
output report it, but, as in INLAvaan, it is not flagged: it is a route
of the solver, not a sign of a poor approximation.

Where the solver lands is checked as the mode is. `vb_step_sigma` is the
step still to go at the VB mean, in posterior SDs: the fixed-point step
there, which is the Newton step of the condition the VB mean solves.
Above 0.1, the level of the mode check, the correction has not settled
and the fit has not converged. The VB mean is then off, and with it the
VB term of the log marginal likelihood and the expansions of
[`loo()`](http://haziqj.ml/inlavanilla/reference/loo.md) and
[`waic()`](http://haziqj.ml/inlavanilla/reference/loo.md). The
skew-normal marginals do not depend on it, since they are scanned around
the mode, and `vb = FALSE` leaves it out. Gaussian marginals
(`marginals = "gaussian"`) are centred on it, and `vb = FALSE` centres
them at the mode instead. Without a gradient, the usual cause is not the
solver but the VB mean itself: from `logpost` values alone it is exact
only for a log-posterior of degree four, and on a strongly skewed
posterior it can be about 0.1 SDs out. Supplying `grad` is the fix.

With a gradient the fallback takes one more step of its own. If it lands
short of settled, it is run once more from where it landed, much as the
mode search reruns its optimiser when the gradient at the mode is not
yet small. This can happen when the optimiser starts where a node fails:
the lowest node around the mode sits on a boundary of the support, say,
so its first value is a failure, and it then reports convergence too
early. The rerun costs nothing when the first run settles, and one more
run of the optimiser when it does not. It is not a check, so
[`summary()`](https://rdrr.io/r/base/summary.html) does not mention it,
but `verbose = TRUE` reports it: “Fell back to nlminb, as …, rerun once
from where it landed”. Without a gradient no rerun is made, since its
optimiser is steered by the same estimate of the expected gradient and
cannot settle any closer.

The skew-normal mean and the VB mean both estimate the posterior mean,
by different routes. Both recover the leading, \\O(n^{-1})\\ shift from
the mode to the mean, so they agree to \\O(n^{-2})\\ (Jamil & Rue, 2026,
Proposition 3.4). `mean_gap` is their gap in posterior SDs. A gap
signals that the higher-order terms of that expansion about the mode
matter, or numerical error (Section 5.1), and a gap above 0.1 is
flagged.

## The skew-normal fit

`nmad` is the largest gap between the fitted skew-normal and the scanned
log-posterior, on the density scale and relative to the peak (Jamil &
Rue, 2026, Section 4.4). A skew-normal has one mode and one direction of
skew, so it cannot follow a marginal with two bumps:

``` r

logpost_mix <- function(theta) {
  log(0.7 * dnorm(theta[[1]], 0, 1) + 0.3 * dnorm(theta[[1]], 3, 1))
}
fit_mix <- inla_vanilla(logpost_mix, c(theta = -0.2))
#> ✔ Found the posterior mode in 4 iterations [4ms]
#> ✔ VB mean correction, max shift 0.39σ [14ms]
#> ✔ Fitted 1 skew-normal marginal [3ms]
#> Warning: Fit diagnostics flagged 3 potential issues:
#> ✖ The fitted marginal deviates from the scanned posterior (NMAD > 0.1) for
#>   `theta` (0.20).
#> ✖ The fitted marginal puts more than 0.05 of its mass beyond the scanned
#>   window (4σ either side of the mode) for `theta` (0.09); its credible
#>   limits rely on extrapolation.
#> ✖ Heavy tails for `theta`: the density is still over 5% of its peak at the
#>   ends of the scanned window (4σ either side of the mode), and the
#>   skew-normal extrapolates beyond it.
#> ℹ Inspect with `diagnostics(fit)` and `diagnostics(fit, type = "param")`;
#>   `vignette("diagnostics", package = "inlavanilla")` explains each check.
```

The fit ends with one warning that lists every problem. The second bump
puts mass far beyond what the Laplace SD expects, so the skew-normal
fits the scan poorly, it stretches its own tail out of the scanned
window to cover the bump, and the tails are heavy.
[`diagnostics()`](http://haziqj.ml/inlavanilla/reference/diagnostics.md)
gives the numbers behind the warning, for the whole fit and for each
parameter:

``` r

diagnostics(fit_mix)
#>               npar          converged         iterations           grad_inf 
#>                  1                  1                  4           1.99e-08 
#>       grad_inf_rel            grad_l2     mode_shift_max          hess_cond 
#>           4.78e-08           1.99e-08           2.04e-08           1.00e+00 
#>         vb_applied       vb_shift_max      vb_kld_global            kld_max 
#>                  1             0.3921             0.0769             0.0769 
#>           kld_mean           nmad_max          nmad_mean unscanned_mass_max 
#>             0.0769             0.1963             0.1963           9.13e-02 
#>        vb_fallback        vb_step_max       mean_gap_max     ridge_bend_max 
#>                  0           6.87e-04             0.7627             0.0000 
#>      n_heavy_tails 
#>                  1
diagnostics(fit_mix, type = "param")
#>   names grad grad_abs grad_rel mode_shift_sigma    kld vb_shift
#> 1 theta    0        0        0                0 0.0769   0.4011
#>   vb_shift_sigma   nmad alpha unscanned_mass vb_step_sigma mean_gap gamma1
#> 1         0.3921 0.1963 4.546         0.0913         7e-04   0.7627      0
#>   tails ridge_bend
#> 1 FALSE          0
```

`plot(type = "sn_fit")` shows the fit against the scan it was fitted to,
as INLAvaan’s plot of the same name does: on the scan coordinate, in
posterior SDs from the mode, the log-posterior along the scan (“Raw”),
the same with the volume correction added (“Corrected”), and the fitted
skew-normal. NMAD is the largest gap between the last two.
`type = "sn_fit_log"` draws the same on a log scale.

``` r

plot(fit_mix, type = "sn_fit")
```

![The scanned log-posterior of the mixture, as a density scaled to a
peak of 1, with its second bump to the right, and the fitted
skew-normal, a single skewed bump that runs below the main peak and
above the dip between the
two.](diagnostics_files/figure-html/mixture-sn-fit-1.png)

Against the exact density, the skew-normal covers both modes but follows
neither:

``` r

# The skew-normal density, from the fit's xi, omega and alpha
dskewnorm <- function(x, xi, omega, alpha) {
  z <- (x - xi) / omega
  2 / omega * dnorm(z) * pnorm(alpha * z)
}
theta <- seq(-4, 7, length.out = 400)
exact <- 0.7 * dnorm(theta, 0, 1) + 0.3 * dnorm(theta, 3, 1)
sn <- fit_mix$sn[1, ]
fitted <- dskewnorm(theta, sn[["xi"]], sn[["omega"]], sn[["alpha"]])
par(mar = c(4, 4, 0.5, 0.5))
plot(
  theta, exact, type = "l", lwd = 2, las = 1, bty = "n",
  xlab = "theta", ylab = "density", ylim = c(0, 0.33)
)
lines(theta, fitted, col = "#eb6834", lwd = 2)
legend(
  "topright", c("exact", "skew-normal fit"),
  col = c("black", "#eb6834"), lwd = 2, bty = "n"
)
```

![The fitted skew-normal against the exact bimodal density. The
skew-normal is a single skewed bump that covers both modes of the exact
density but follows
neither.](diagnostics_files/figure-html/mixture-plot-1.png)

## Ridge bend

Each marginal is scanned along a straight line through the mode, the
conditional mode path (CMP) on which the Laplace Gaussian puts the other
parameters at their conditional mode given this one (Jamil & Rue, 2026,
Lemma 3.1). The straight line is a first-order approximation to a
possibly curved ridge of high posterior density (Section 3.2.1). When
the ridge bends away from the CMP, the scan misses log-density, and the
marginal comes out too narrow. The ridge bend, `ridge_bend`,
inlavanilla’s own check, is the log-density, in nats, that one Newton
step in the other parameters would recover at \\\pm 2\\ posterior SDs.

This happens when a parameter’s spread depends on another parameter that
the data pin down poorly. A normal mean with an unknown SD, from 12
observations, is the textbook case:

``` r

set.seed(4)
y12 <- rnorm(12, 10, 2)
logpost_norm <- function(theta, y) {
  sum(dnorm(y, theta[["mu"]], exp(theta[["log_sigma"]]), log = TRUE)) +
    dnorm(theta[["mu"]], 0, 100, log = TRUE) +
    dnorm(theta[["log_sigma"]], 0, 10, log = TRUE)
}
fit_norm <- inla_vanilla(
  logpost_norm, c(mu = 9, log_sigma = 0.5), y = y12, silent = TRUE
)
signif(fit_norm$diagnostics$marginals[, c("nmad", "ridge_bend")], 2)
#>              nmad ridge_bend
#> mu        1.9e-06    3.3e-01
#> log_sigma 1.1e-02    4.9e-07
```

[`summary()`](https://rdrr.io/r/base/summary.html) flags the ridge bend
of `mu`. The exact marginal SD of `mu`, from a grid, shows how much
narrower the skew-normal is:

``` r

grid <- expand.grid(
  mu = seq(7, 15, length.out = 401),
  log_sigma = seq(-1, 2.5, length.out = 401)
)
lp <- apply(grid, 1, function(g) logpost_norm(g, y12))
w <- exp(lp - max(lp))
w <- w / sum(w)
m <- sum(w * grid$mu)
sds <- c(
  exact = sqrt(sum(w * (grid$mu - m)^2)),
  skew_normal = fit_norm$summary["mu", "SD"]
)
sds
#>       exact skew_normal 
#>   0.6125832   0.5308595
```

The skew-normal SD of `mu` is 13% too small here. The gap shrinks as
\\1/n\\: for a normal mean with an unknown SD it is about 3% at 50
observations and under 1% at 200. For this model the ridge bend \\b\\
gives the gap: the SD is too small by \\1 - \sqrt{1 - 3b/4}\\, roughly
\\3b/8\\, and \\b\\ is close to \\4/n\\. From its ridge bend of 0.33,
the formula gives 13%.
[`summary()`](https://rdrr.io/r/base/summary.html) flags a ridge bend
above 0.25, where the formula gives 10%. In other models it is a rough
guide to how much: on test models it was mostly within a factor of 1.5
of the true gap, and it overstated the narrowing in logistic and Poisson
regressions.

## Tails and unscanned mass

Both checks are about the edge of the scanned window, but they look at
different things:

- `tails` reads the **log-posterior itself** at the two ends of the
  scan. It is flagged when either end is still within 3 nats of the
  peak, that is, when the density there is over 5% of its peak. A
  Gaussian is 8 nats down at \\\pm 4\\ SDs.
- `unscanned_mass` is the mass the **fitted skew-normal** puts beyond
  the window, where it is extrapolated (as INLAvaan’s
  [`diagnostics()`](http://haziqj.ml/inlavanilla/reference/diagnostics.md)
  reports it). A Gaussian marginal gives \\2\Phi(-4) = 6.3 \times
  10^{-5}\\; the larger it is, the more of the marginal is extrapolated.

Neither one can see the mass the true posterior has beyond the window.
The figure shows three one-parameter posteriors whose exact densities
are known. The x axis is in posterior SDs from the mode, so the scanned
window is always \\-4\\ to \\4\\ (the grey band).

![Three columns, one posterior each. Top row: the log density with the
scan points, the ends of the scan and a line 3 nats below the peak.
Bottom row: the density, with the mass outside the scanned window shaded
for the exact posterior and for the skew-normal fit. In the first column
both are tiny. In the second the scan ends sit above the 3-nat line, so
tails is flagged, and the exact posterior has much more mass outside the
window than the fit. In the third the scan ends are just below the line,
and the exact posterior again has more mass outside than the
fit.](diagnostics_files/figure-html/tails-figure-1.png)

In the top row, `tails` reads the open circles against the dashed line 3
nats below the peak. In the bottom row, the orange shading is
`unscanned_mass` and the grey shading is the exact posterior’s mass
beyond the window.

1.  **Well behaved.** The scan ends are 5.8 and 11.6 nats down, either
    side of a Gaussian’s 8, as the Poisson log rate is a little skewed,
    and neither the fit nor the posterior has much mass outside. Both
    checks are rightly quiet.
2.  **Heavy inside the window.** A Cauchy posterior is still high at the
    ends of the scan, so `tails` fires, and the fit warns.
    `unscanned_mass` stays under 0.05: the skew-normal has Gaussian
    tails, so whatever it fits, it puts little mass outside the window.
    It exceeds 0.05 only when the fit itself is very wide or very
    skewed, as for the mixture above.
3.  **Heavy beyond the window.** Here the data can only bound a variance
    from above, so its log has a long left tail that only the prior
    holds back. The tail flattens out beyond the window: the scan end is
    already more than 3 nats down, and the skew-normal puts little mass
    outside. Both checks are quiet, but the exact posterior has 0.025 of
    its mass outside, against 0.006 for the fit.

The third case is the known limit of a skew-normal at a boundary: no
check that reads only the scanned window can see a tail that starts
beyond it. Jamil and Rue (2026, Section 4.4) find the same miss in the
lower tail of a variance near zero, the boundary pile-up that a unimodal
skew-normal cannot capture. In the structural equation model of the
[Rcpp article](http://haziqj.ml/inlavanilla/articles/rcpp.md), two log
residual variances behave like this, more strongly. Their skew-normal
marginals are too narrow on the log scale, although on the variance
scale the error is much smaller. INLAvaan gives the same marginals
there.

## The warning

The fit-time warning is a condition of class
`inlavanilla_diagnostics_warning`, which is also an
`inlavanilla_warning`. Its headline counts the problems, as INLAvaan’s
does, and it names at most three parameters for each;
[`summary()`](https://rdrr.io/r/base/summary.html) lists them all. With
`silent = TRUE` nothing is printed, and every warning raised during the
fit is kept in `fit$warnings` as plain text, a line per bullet with no
wrapping, which suits simulations. Here it is wrapped to fit the page:

``` r

fit_mix_quiet <- inla_vanilla(logpost_mix, c(theta = -0.2), silent = TRUE)
length(fit_mix_quiet$warnings)
#> [1] 1
bullets <- strsplit(fit_mix_quiet$warnings, "\n")[[1]]
writeLines(strwrap(bullets, width = 76, exdent = 2))
#> Fit diagnostics flagged 3 potential issues:
#> ✖ The fitted marginal deviates from the scanned posterior (NMAD > 0.1) for
#>   `theta` (0.20).
#> ✖ The fitted marginal puts more than 0.05 of its mass beyond the scanned
#>   window (4σ either side of the mode) for `theta` (0.09); its credible
#>   limits rely on extrapolation.
#> ✖ Heavy tails for `theta`: the density is still over 5% of its peak at the
#>   ends of the scanned window (4σ either side of the mode), and the
#>   skew-normal extrapolates beyond it.
#> ℹ Inspect with `diagnostics(fit)` and `diagnostics(fit, type = "param")`;
#>   `vignette("diagnostics", package = "inlavanilla")` explains each check.
```

To handle the warning in code, catch its class:

``` r

fit <- withCallingHandlers(
  inla_vanilla(logpost_mix, c(theta = -0.2)),
  inlavanilla_diagnostics_warning = function(w) {
    message("Check this fit.")
    invokeRestart("muffleWarning")
  }
)
#> ✔ Found the posterior mode in 4 iterations [1ms]
#> ✔ VB mean correction, max shift 0.39σ [1ms]
#> ✔ Fitted 1 skew-normal marginal [4ms]
#> Check this fit.
```

## Checks on a supplied gradient

When you pass `grad`, it is compared with finite differences of
`logpost` at `init`, per posterior SD, and its Jacobian at the mode is
checked for symmetry, as the Hessian of a function must be. Either
failure gives a warning of its own, because a wrong gradient gives a
wrong mode and wrong marginals.

## References

Jamil, H., & Rue, H. (2026). *Approximate Bayesian inference for
structural equation models using integrated nested Laplace
approximations* (arXiv:2603.25690 \[stat.ME\]). arXiv.
<https://doi.org/10.48550/arXiv.2603.25690>
