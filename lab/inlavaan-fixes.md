# Fixes for INLAvaan found while building inlavanilla

For the agent working on `haziqj/INLAvaan`. Line numbers refer to INLAvaan at
commit `db52ad7`.

inlavanilla ports INLAvaan's skew-normal code, and its test suite (closed-form
references, quadrature, mutation testing) turned up the issues below. Every one
is fixed in inlavanilla. The fixed code is in `haziqj/inlavanilla`, branch
`claude/dazzling-ritchie-cj1doo`, in `R/skew-normal.R`, `R/qsnorm-fast.R` and
`R/sampling.R`, with tests in `tests/testthat/test-skew-normal.R` and
`tests/testthat/test-sampling.R`. Port the fixes and those tests. The
functions have the same names except where noted.

Items 1 to 5 are bugs that change results. Items 6 to 8 are latent or
low-impact. Items 9 to 12 are not bugs.

## 1. `qsnorm_fast()` is a starting value, not an accurate quantile (high)

**Where:**
- `R/qsnorm_fast.R:46` documents "< 10^-7 relative error". That is wrong.
- The function is used directly for results in four places:
  - `R/posterior-sampling.R:75`, the copula draws in `sample_params()`;
  - `R/posterior-marginals.R:146`, the reported skew-normal quantiles;
  - `R/posterior-sampling.R:752`, `summarise_samples_sn()`;
  - `R/norta.R:47` and `:51`, the NORTA adjustment.

**What:** away from alpha = 0 and +-1 the error reaches about 5e-2 in
standardised units. In the bulk it is 5e-3 to 1e-2; in the tails it is 1e-2 to
6e-2; at alpha = 1.5 and p = 1e-3 it is 6e-2. Luu's original (`sncdfinv`) uses
a tolerance of 0.01: it was meant as a starting value for refinement. In
posterior-SD units, copula draws and reported quantiles are off by about 0.01
SD in the bulk and up to about 0.04 SD in the light tail.

**Reproduce:**

```r
qref <- function(p, a) uniroot(function(z) psnorm(z, 0, 1, a) - p, c(-10, 10), tol = 1e-14)$root
qsnorm_fast(1e-3, 0, 1, 1.5) - qref(1e-3, 1.5)  # -0.061
qsnorm_fast(0.75, 0, 1, 2) - qref(0.75, 2)      #  0.0059
```

**Fix:**
- For reported quantiles, use `qsnorm()` from inlavanilla. It runs Newton's
  method on the log of the smaller tail probability to convergence
  (`qsn_std_logtail()`). The skew-normal is log-concave, so convergence is
  monotone from any start. Relative error is about 1e-12.
- For copula draws, use `sn_score_map()`. It builds a cubic spline through 141
  knots at z = -7..7, solved by `qsn_std_logtail()`, and maps normal scores
  z -> Q(Phi(z)) directly, with linear extension beyond the knots. Error is
  below 1e-6 posterior SDs, at about 5 ms per parameter, plus one spline
  evaluation per draw.
- Correct the documentation at `R/qsnorm_fast.R:46`.

## 2. Owen's T loses relative accuracy for a > 1, so `psnorm()`'s heavy tail is wrong (medium)

**Where:** `R/skew-normal.R:361-363`, `owen_t()`, the `big` branch.

**What:** `0.5 * Phi(h) + 0.5 * Phi(ah) - Phi(h) Phi(ah) - T(ah, 1/a)` is a
difference of numbers close to one when h is large. psnorm's heavier tail
(upper for alpha > 1, lower for alpha < -1) has the following relative error:

| z | Relative error |
|---|---|
| 5 | 2.7e-10 |
| 6 | 2.8e-8 |
| 7 | 2.1e-5 |
| 8 | 3.5e-2 |
| 9 and above | 0.5, a factor of 2 |

```r
f <- function(x) 2 * dnorm(x) * pnorm(2 * x)
psnorm(9, 0, 1, 2, lower_tail = FALSE)                   # 1.13e-19
integrate(f, 9, Inf, rel.tol = 1e-12, abs.tol = 0)$value # 2.26e-19
```

**Fix:** the same identity written with upper-tail probabilities
Q(x) = 1 - Phi(x):

```r
qh <- stats::pnorm(h, lower.tail = FALSE)
qah <- stats::pnorm(ah, lower.tail = FALSE)
sgn * ifelse(big, 0.5 * qh + 0.5 * qah - qh * qah - t_red, t_red)
```

## 3. Off-by-one trapezoid weights in `fit_skew_normal()` (medium)

**Where:** `R/skew-normal.R:58`.

**What:** `w_int[-1] <- 0.5 * (c(diff(x[-1]), 0) + c(0, diff(x[-length(x)])))`
gives the second grid point half its weight. On the scan grid
`seq(-4, 4, length.out = 21)` the weights are 0.2, 0.2, 0.4, ..., 0.4, 0.2 and
sum to 7.8, not 8. The moment-based starting values are then not
mirror-symmetric. Mode-centred exact skew-normal scans with alpha = -0.1 or
-0.2 collapse to alpha of about 0, while +0.1 and +0.2 are recovered to 1e-14.

**Fix:**

```r
dx <- diff(x)
w_int <- 0.5 * (c(dx, 0) + c(0, dx))
```

## 4. `fit_skew_normal()` stalls at alpha = 0 on mildly skewed scans (medium)

**Where:** `R/skew-normal.R:245-256`. There is a single `nlminb()` start.

**What:** the volume correction adds `gamma1 * t` to the scan, which moves its
mode off t = 0. For exact skew-normal inputs with |alpha| <= 0.5 and the mode
at +-0.5 on the grid, 9 of 24 cases returned a Gaussian (|alpha| < 1e-4,
nmad up to 5e-3). The moment-based starting skew, computed from a scan
truncated at +-4, can have the wrong sign. The shape is weakly identified near
zero (a skew-normal differs from a shifted Gaussian only at third order), so
nlminb stops there.

```r
x <- seq(-4, 4, length.out = 21)
y <- dsnorm(x, 0.2672, 1.0276, 0.3, log = TRUE); y <- y - max(y)
fit_skew_normal(x, y, temp = 1)[c("alpha", "nmad")]  # alpha ~ -4e-5, nmad 1.2e-3
```

**Fix:** fit from the moment start and from its mirror image (delta -> -delta,
or +-0.5 when |delta| < 0.2), and keep the fit with the smaller objective.
This roughly doubles the cost of the skew-normal fit, about 1 ms per
parameter.

Two related changes in inlavanilla's version are worth taking too. Points below
the threshold have zero weight, so they are dropped before the fit. This changes
nothing numerically, and it keeps a -Inf or NaN height out of the objective,
where 0 * Inf would give NaN.

## 5. `psnorm()` returns NA at -Inf and Inf (low)

**Where:** `R/skew-normal.R:383`. `ok <- is.finite(z) & ...` treats an infinite
argument as invalid.

**Fix:** return 0 at -Inf and 1 at Inf, as `stats::pnorm()` does. Keep NA for
a missing value or a non-positive scale.

## 6. `qsnorm_fast()` leaves some elements of a vector call at 0 (latent)

**Where:** `R/qsnorm_fast.R:113-114` and `:253`.

**What:**
- `res` is initialised to 0, and the special values for u = 0 and u = 1 are
  set before the alpha < 0 fold. Points outside `mask` (alpha of 0 or 1, or u
  folded to 0 or 1) therefore keep 0 in the main path.
- The early-return path (`!any(mask)`) does give them `z`. So
  `qsnorm_fast(c(0.2, 0.7), alpha = c(0, 2))` returns 0 for the first element,
  where `qnorm(0.2)` is right.
- For p < 1.1e-16 with alpha < 0, 1 - p rounds to 1: a scalar call gives -Inf
  and a vector call gives 0.

INLAvaan always calls it with a scalar alpha, so this is latent.

**Fix:** add `res[!mask] <- z[!mask]` after `res[mask] <- final_vals`. For tiny
p, start any refinement from the Gaussian tail; see `qsn_std_logtail()` in
inlavanilla, which starts from sqrt(-2 (log p - log 2)) on the heavy side and
sqrt(-2 log p / (1 + alpha^2)) on the light side.

## 7. Light-tail CDF by cancellation (low; only matters for refinement)

psnorm's light tail (lower for alpha > 0) is Phi(z) - 2 T(z, alpha), again a
difference of nearly equal numbers. Below about 1e-15 relative to Phi(z) it is
noise. If you add Newton refinement (item 1), use inlavanilla's
`log_sn_light_tail()` there. Since T(h, Inf) = Q(h)/2, the light tail is
(1/pi) times the integral over x > alpha of exp(-h^2 (1 + x^2)/2) / (1 + x^2),
which after the substitution x = alpha + u / (h^2 alpha) is a smooth
40-node Gauss-Laguerre integral. Quantiles are then exact to about 1e-13 SDs
down to p = 1e-40.

## 8. Mode and Hessian numerics at extreme scales (low for SEM)

- **Mode search.** `nlminb()` runs on raw parameters with no `scale`. With
  posterior SDs far from 1 (1e4, or 1e-6), it can stop up to a posterior SD
  short of the mode and report success (convergence 0), or report false
  convergence. SEM parameters are usually near unit scale, so the risk is
  small. inlavanilla checks the whitened gradient at the mode,
  max |L' g|. When it exceeds 1e-3 it reruns `nlminb(..., scale =
  sqrt(diag(H)))`. It warns from INLAvaan's Newton step in posterior SDs,
  |(Sigma g)_j| / sd_j above 0.1, rather than from nlminb's return code.
- **`fast_hessian()`** (`R/utils-optim.R:17`, used when
  `numerical_grad = TRUE`) takes second differences with a fixed absolute
  step of 1e-5. That is too small for rounding at posterior SDs above about
  1e3, and too large below about 1e-5. inlavanilla's `fd_scales()` sets
  per-coordinate steps from a cheap curvature search, and uses an
  off-diagonal formula with two evaluations per pair instead of four.
- **`density_mode()`** (`R/utils.R:257`). `optimize()`'s tolerance also has
  a sqrt(eps) |x| term, which swamps the scale when |location| >> scale (for
  example an intercept of 2000 with SD 1e-4: 0.24 SD error). Search on the
  standardised scale instead; see `sn_mode()` in inlavanilla.

## 9. Not a bug: VB integration rule (optional improvement)

`lab/vb-benchmark.R` in inlavanilla compares rules for the VB mean correction
against exact VB solutions on nine models. The fully symmetric degree-5
cubature rule (2m^2 + 1 nodes on the principal axes; weights
1 + (m^2 - 7m)/18 for the centre, (4 - m)/18 for the axes, 1/36 for the
pairs, scale sqrt(3)) was the most accurate on every model. Error is in
posterior SDs at the default tolerance:

| Model | m | Degree-5 | Sobol-64 (default) | Gauss-Hermite axis |
|---|---|---|---|---|
| Normal, n = 10 | 2 | 1.3e-4 | 0.019 | 0.019 |
| Logistic regression | 5 | 0.0026 | 0.0093 | 0.024 |
| Logistic regression | 10 | 0.0087 | 0.022 | 0.034 |
| Logistic regression | 30 | 0.016 | 0.058 | 0.051 |

The cost is 2m^2 + 1 gradient calls per sweep, against 65 for Sobol-64. It is
deterministic and needs no Sobol table.

## 10. Not a bug: small-n limit of the linear scan path (for the paper)

The CMP scan is a straight line, and the volume correction is first order.
When a parameter's conditional spread depends on another, poorly determined
parameter (a noise scale or a variance component), the marginal comes out too
narrow by a factor that vanishes as O(1/n).

- **Normal mean with unknown variance:** the SD is short by exactly
  sqrt(2(a_n - 1)/(2 a_n + 1)), about 11% at n = 12, 3% at n = 50 and 0.7% at
  n = 200.
- **Michaelis-Menten (Puromycin, n = 12):** 18% too narrow.
  - One Newton step towards the conditional mode at each scan point halves
    the gap (-8%).
  - The full Laplace marginal (conditional mode and conditional
    log-determinant at each point) brings it to -1%, at about 10 times the
    cost of the marginal stage.

- **Bollen's Political Democracy SEM (39 parameters, n = 75), the model of
  INLAvaan's SEM article with a mean structure, default priors:** checked
  against HMC (4 chains of 10,000 draws, R-hat 1.000) and against full
  Laplace marginals on a wide grid, which agree with each other.
  - The log residual variance of `x2` and the log disturbance variance of
    `dem65` have skew-normal SDs of 0.38 and 0.46 times the posterior SD.
    Their 2.5% quantiles are -3.46 against -6.98 and -3.60 against -7.34,
    two or more posterior SDs off; the upper quantiles agree.
  - Cause: the data cannot keep these variances away from zero, so on the
    log scale the left tail is exponential, with slope 0.50: the
    `gamma(1,.5)[sd]` prior with its Jacobian. A skew-normal cannot follow
    it, and the scan (4 Laplace SDs) stops before the tail flattens.
  - On the variance scale, the scale INLAvaan reports, the error is moderate.
    Against Stan NUTS on an independently written model with INLAvaan's
    default priors (`lab/poldem.stan` in inlavanilla), INLAvaan's own summary
    has x2~~x2 mean 0.132 against 0.114 (0.23 posterior SDs), SD 0.066
    against 0.077, 2.5% quantile 0.031 against 0.001; dem65~~dem65 mean 0.273
    against 0.236 (0.17 SDs). The other 37 parameters: means within 0.09 SDs,
    SD ratios 0.93-1.01. A validation against MCMC on the natural scale with
    a tolerance of about 0.25 SDs would pass this model.
  - The other 37 parameters are within 7% in SD (median 3%), 0.08 SDs in the
    mean and 0.21 SDs in the quantiles.
  - inlavanilla's fit matches INLAvaan's here (modes within 2e-4 SDs, log
    marginal likelihood -1686.85 against -1686.86), so INLAvaan's marginals
    of these two parameters will be off the same way. Worth checking any
    validation against MCMC on the log scale, where it shows most.

See `lab/DESIGN.md`, `lab/scan-path-prototype.R`, `lab/sem-marginal-check.R`
and `lab/sem-hmc.R` in inlavanilla.

## 11. Not a bug: rename `scan_end_mass` to `unscanned_mass`

In October 2026 the owner renamed this diagnostic in inlavanilla so that its
name says what it measures: the mass of the fitted marginal that no scan
point supports. INLAvaan is to follow, so that `diagnostics()` uses the same
names in both packages.

- Per parameter: `scan_end_mass` becomes `unscanned_mass`.
- Global: `scan_end_mass_max` becomes `unscanned_mass_max`.
- The printed name is "Unscanned", and prose says "unscanned mass".
- The level, 0.05, and the definition are unchanged.

inlavanilla also renamed two checks of its own that INLAvaan does not have:
`bend` became `ridge_bend` and `sn_vs_vb` became `mean_gap`.

## 12. Suggestion: check where the VB fallback lands

inlavanilla's VB correction, ported from INLAvaan, accepted the fallback
optimiser's result without checking it. It now measures the step still to go
at the VB mean, in posterior SDs: the fixed-point step there, which is the
Newton step of the VB stationarity condition, as `mode_shift_sigma` is for
the mode. On the fixed-point path it is the last step, below the tolerance;
after the fallback it costs one more sweep. Above 0.1, the mode check's
level, the correction has not settled: the fit warns, `summary()` flags it,
and `print()` says the fit has not converged. `diagnostics()` reports it as
`vb_step_sigma` per parameter and `vb_step_max` globally. INLAvaan's VB mean
feeds its marginal likelihood and its LOO in the same way, so the same check
would help there.
