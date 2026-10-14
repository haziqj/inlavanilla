# Research: calibrating the scan-path bend threshold

Date: 2026-10-13, with the cross-model check of 2026-10-14. Reproduce the
calculations with `Rscript lab/bend-threshold.R` and
`Rscript lab/bend-calibration.R` from the package root.

The diagnostic was renamed `ridge_bend` on 2026-10-14, after this note was
written. The note keeps its old name, `bend`. The scripts use the new one.

## Decision (2026-10-14)

`summary()` flags `bend` above **0.25 nats**, up from 0.1. The level comes
from an accuracy target chosen first: a marginal SD about 10% too small. In
the normal-mean reference family below, a 10% shortfall is a bend of 0.253.
The cross-model check at the end of this note shows that bend predicts the
scan-path narrowing well enough for that target to carry over, and that 0.1,
a 3.8% target, is too tight to separate the marginals that are off from
those that are not. The flag still certifies nothing: the zero-bend
counterexample below stands, and the SEM's `logv_x2`, which is badly off for
another reason, is no longer flagged.

The original recommendation of 2026-10-13, kept for the record:

> There is no universal, literature-established cutoff for this diagnostic.
> The existing 0.1 is a reasonable screening heuristic, but not an error
> bound. A principled calibration starts by choosing an accuracy target and
> a reference family. For a 5% relative marginal SD shortfall in the
> normal-mean/unknown-scale reference family derived below, the threshold is
> 0.13 nats. Do not simply increase the package default on this evidence
> alone: the SEM investigation in `lab/DESIGN.md` includes a badly inaccurate
> marginal with bend about 0.11, and increasing the cutoff would hide that
> flag. Keep 0.1 until calibration across the intended model families
> justifies a change.

The cross-model check is that calibration. The SEM marginal with bend 0.11
is `logv_x2`, and bend flagged it by chance (see below).

For a diagnostic that directly answers whether marginal uncertainty is
materially affected, compare SDs and credible limits from the current scan
with a more accurate conditional-mode and conditional-volume calculation.
Flag a chosen relative SD change (for example, 5%) and separately assess
credible-limit changes. Even that comparison is sensitivity to a more
accurate approximation, not a guarantee against unvisited tails or modes.

The research of 2026-10-13 changed no package behaviour; the decision above
did.

## What the package actually calculates

`R/marginals.R::scan_bend()` computes

\[
B_j=\max_{t\approx-2,+2}\frac12 g_{-j}(t)^T C_j g_{-j}(t),\qquad
C_j=\Sigma_{-j,-j}-\frac{\Sigma_{-j,j}\Sigma_{j,-j}}{\Sigma_{jj}}.
\]

Here the gradient is evaluated on the straight scan path and the covariance
is held at its joint-mode Laplace value. The code selects the scan grid points
nearest to -2 and +2, so a custom grid need not evaluate exactly at those
coordinates. The SD unit is the Laplace SD at the mode.

This is the log-density improvement predicted by a quadratic model in the
nuisance parameters. It resembles half the squared Newton decrement, but
uses the Hessian from the joint mode, not the local conditional Hessian at
the scan point. It need not equal the actual log-density gain after a step,
or the gain at the true conditional optimum. Boyd and Vandenberghe explain
the Newton-decrement interpretation in their
[optimization slides, slide 9.23](https://web.stanford.edu/~boyd/cvxbook/bv_cvxslides.pdf).

In this checkout, `R/diagnostics.R` sets `summary_levels$bend = 0.1`.
`R/methods.R::check_flags()` uses it for summaries; the fit-time warning omits
bend. The originating commit, `a3795be`, gave a normal-scale example and a
rough relation to SD error, but no formal calibration of the cutoff.

An estimated gain of 0.1 means an estimated unnormalized height ratio
`exp(0.1) = 1.1052`. It does **not** mean 10% marginal SD error, a KL divergence
of 0.1, or a probability of approximation failure. There is no chi-squared
reference distribution for this deterministic quantity. A Wilks cutoff or a
Pinsker bound cannot be applied to it just because its unit is nats.

## What the literature supports

- [Jamil and Rue, Section 3.2](https://arxiv.org/html/2603.25690v1) describes
  the straight conditional mean path as an approximation to a curved ridge.
  It also separates conditional peak height from the conditional volume
  contribution, minus half a log determinant. This explains what bend measures
  and why it cannot capture every marginal error; it does not supply this
  package-specific cutoff. The linked text is explicitly version 1.
- [McDonald and Campbell, Section 5](https://arxiv.org/html/2411.01697v1)
  calibrates a different Laplace diagnostic using Student-t reference
  functions at the boundary of an acceptable error. One proposed boundary
  is 5% error in the integral. This supports the *calibration strategy*, not
  a 5% SD rule or a 0.13 bend cutoff. Their diagnostic uses Bayesian
  quadrature and requires choices of grid and hyperparameters.
- [Kasprzak, Giordano and Broderick](https://arxiv.org/abs/2209.14992)
  develops finite-sample Laplace error bounds that can control means and
  variances under stated assumptions. Such guarantees require additional
  analysis; they are not consequences of two local gradient evaluations.
- [Giordano, Broderick and Jordan](https://jmlr.org/papers/v19/17-670.html)
  relates derivatives of posterior expectations under log-density
  perturbations to covariances. This supplies a useful local sensitivity
  interpretation for the derivation below.

The numerical calibrations and counterexample below are our derivations for
this code, not thresholds claimed by those papers.

## An exact reference-family calibration

Take normal observations with unknown mean mu and SD sigma, write
eta = log(sigma), and use a flat prior in (mu, eta), equivalently a density
proportional to 1/sigma in (mu, sigma). This improper reference prior gives
a proper posterior when n > 1 and the residual sum of squares S is positive;
the marginal variance exists when n > 3.

Writing d = mu - sample_mean, the log posterior is

\[
\ell(\mu,\eta)=-n\eta-\tfrac12 e^{-2\eta}(S+nd^2)+\text{constant}.
\]

At the joint mode, sigma_hat^2 = S/n. The Laplace variance of mu is S/n^2,
and the conditional variance of eta is 1/(2n). On the straight mu scan at
t Laplace SDs, the eta gradient is t^2, giving

\[
B(t)=\frac{t^4}{4n},\qquad B(2)=\frac4n.
\]

The straight scan is exactly Gaussian in mu, and its linear volume
correction vanishes by symmetry. Therefore the package's fitted marginal
SD is sqrt(S)/n. Integrating eta exactly gives

\[
p(\mu\mid y)\propto(S+nd^2)^{-n/2},\qquad
\operatorname{Var}(\mu\mid y)=\frac{S}{n(n-3)}.
\]

Consequently, the SD shortfall **relative to the exact SD** is

\[
\epsilon=1-\frac{s_{\rm scan}}{s_{\rm exact}}
=1-\sqrt{1-\frac3n}
=1-\sqrt{1-\frac34 B}.
\]

Solving for a specified tolerance gives

\[
\boxed{B_{\rm crit}(\epsilon)=\frac43(2\epsilon-\epsilon^2).}
\]

| Allowed SD shortfall | Reference-family bend cutoff |
| --- | ---: |
| 1% | 0.02653 |
| 2.5% | 0.06583 |
| 5% | 0.13000 |
| 10% | 0.25333 |

Thus 0.1 corresponds to 3.823% SD shortfall in this family. This explains
why the existing number is of a sensible order, while exposing its implicit
accuracy tolerance. The 5% tolerance is still a user or package policy
choice; mathematics translates it rather than choosing it.

The benchmark's true conditional-mode height gain at t = 2 is
`2 - n/2 * log(1 + 4/n)`, not `4/n`. This also illustrates why the actual
diagnostic must be described as a *predicted* gain.

The accompanying script fits the package with analytic gradients and
`vb = FALSE`, and checks its SD against the formula. Independent quadrature
of the exact Student-t marginal confirms the variance:

| n | bend | SD shortfall |
| ---: | ---: | ---: |
| 12 | 0.33333 | 13.397% |
| 20 | 0.20000 | 7.805% |
| 40 | 0.10000 | 3.823% |
| 80 | 0.05000 | 1.893% |
| 200 | 0.02000 | 0.753% |

The script also evaluates n = 4/0.13 as a continuous interpolation of the
reference kernel, obtaining exactly 5% shortfall. This is not a claim that
a sample can contain a fractional number of observations. Different priors
or parameterizations can change the calibration.

## Why the local coefficient is approximately 3/8

At an accurately located joint mode, the nuisance gradient vanishes.
Its first derivative along the Gaussian conditional mode path also
vanishes. For a smooth, regular posterior its leading term is therefore
quadratic in t, so the predicted height loss is quartic:

\[
b(t)\approx B t^4/16.
\]

For a nearly standard normal scanned marginal, the first-order change in
variance under this log-density perturbation is

\[
\delta\operatorname{Var}(t)
\approx\operatorname{Cov}_{N(0,1)}(t^2,b(t))
=\frac B{16}(15-3)=\frac{3B}{4}.
\]

The relative SD change is consequently about 3B/8. A 5% tolerance gives
0.1333 at first order, consistent with the exact 0.13 calibration above.

This is a perturbative calculation: exp(-t^2/2 + B*t^4/16) is not a proper
density for positive B. The quartic describes a local correction, not an
expression to exponentiate over the whole real line. The conversion also
assumes a nearly Gaussian baseline, small correction, negligible omitted
volume effects, and broadly similar behavior on the two sides. The maximum
of two sampled losses cannot verify those assumptions.

## A zero-bend counterexample

Consider the proper joint density with log kernel

\[
\ell(x,y)=-\tfrac12x^2-\tfrac12 e^{-a x^2}y^2,\qquad 0<a<1.
\]

The joint mode is (0,0), the Hessian there is the identity, and the
conditional mode of y is always zero. The x scan therefore has bend
exactly zero at both probe points (indeed everywhere). Its linear volume
correction is also zero. The package returns an N(0,1) marginal for x.

However, integrating y contributes exp(a*x^2/2), so the exact x marginal is
N(0, 1/(1-a)). With a = 0.25, the fitted SD is 1 but the exact SD is
1.154701: **13.397% too small despite zero bend**. Both joint second moments
are finite at this value of a. The script confirms this package behavior.

The missed term is quadratic conditional volume, not scan-path curvature.
This is why no positive threshold on bend alone can certify marginal SD
accuracy. Sampling only two locations also cannot rule out failures farther
into the tails or additional modes.

## Implications for implementation

1. Retain bend as a cheap continuous indicator of scan-path curvature.
   Describe it as a predicted log-density gain, and label any cutoff as a
   screening convention with an explicit calibration family.
2. A reference-calibrated option can map a relative SD tolerance epsilon
   to `(4/3) * (2*epsilon - epsilon^2)`. At 5%, this is 0.13. It must not
   promise a general error bound or replace the other diagnostics.
3. A more informative, optional accuracy check would recompute conditional
   modes and conditional Hessian determinants over a sufficiently wide
   marginal grid, normalize the resulting marginal, and compare its SD
   and quantiles to those currently reported. Retain the grid-based
   distribution for the comparison rather than forcing both versions into
   the same skew-normal family, which could hide a shared fitting error.
4. Calibrate and evaluate that check against exact integration or accurate
   MCMC across normal-scale, nonlinear regression, skewed, heavy-tailed,
   near-boundary and volume-curvature examples. Measure missed inaccurate
   fits as well as unnecessary flags, with held-out model settings.

A simple replacement of 0.1 by another constant cannot supply that missing
validation.

## Cross-model check (2026-10-14)

`lab/bend-calibration.R` compares bend with the error of each skew-normal
marginal on 51 parameters: 12 small models against a grid over +-8 Laplace
SDs, and the 39 parameters of the Political Democracy SEM against HMC
(`lab/sem-hmc.R`). "Predicted" is the normal-family shortfall
`1 - sqrt(1 - 3 * bend / 4)`; "cover" is the exact mass inside the
skew-normal 95% interval.

| Model | Parameter | Bend | Predicted | SD short | Cover |
| --- | --- | ---: | ---: | ---: | ---: |
| Poisson, n = 8 (interface tests) | slope | 0.139 | 5.3% | 3.3% | 94.1% |
| Logistic, mtcars (main vignette) | `b0` | 0.144 | 5.5% | 4.5% | 93.7% |
| Normal mean, n = 30 | `mu` | 0.133 | 5.1% | 5.1% | 93.5% |
| Normal mean, n = 20 | `mu` | 0.200 | 7.8% | 7.8% | 92.7% |
| Weibull, n = 15 | `log_scale` | 0.294 | 11.7% | 11.0% | 91.8% |
| Logistic, n = 30 (simulated) | `b0` | 0.308 | 12.3% | 4.1% | 93.7% |
| Normal mean, n = 12 | `mu` | 0.333 | 13.4% | 13.4% | 91.1% |
| Logistic, n = 30 (simulated) | `b2` | 0.346 | 13.9% | 9.1% | 92.4% |
| Puromycin, n = 12 | `log_Vm`, `log_K` | 0.35 | 14% | 17.8% | 89.7% |

The other small-model parameters have bend below 0.09 and SDs within about
5% (one of them 5% too wide). On the SEM:

| Parameter | Bend | Predicted | SD short | Cover |
| --- | ---: | ---: | ---: | ---: |
| 33 others | 0.015 to 0.094 | under 3.6% | -0.1% to 6.8% | 93.2% to 95.0% |
| `b` | 0.108 | 4.1% | 6.6% | 93.3% |
| `logv_x2` | 0.114 | 4.4% | 59.9% | 81.0% |
| `z_y2_y6` | 0.115 | 4.4% | 3.1% | 94.3% |
| `z_y4_y8` | 0.127 | 4.9% | 6.7% | 93.6% |
| `logv_x1` | 0.141 | 5.4% | 6.1% | 93.7% |
| `logv_dem65` | 0.436 | 18.0% | 53.0% | 81.1% |

Findings:

- **Bend predicts the scan-path narrowing.** The normal-family formula is
  exact for the normal mean and within about a factor of 1.5 in most other
  models. It overstates the narrowing in GLMs (most of all for `b0` at
  n = 30, 12% against 4%) and understates it for Puromycin, where the linear
  volume correction adds to the error (14% against 18%).
- **0.1 is too tight to separate anything.** It is a 3.8% target in the
  reference family. Of the 16 parameters it flags, 4 are under 5% short,
  including the intercept of the main vignette, whose 95% interval holds
  93.7% of the exact mass. On the SEM, 4 of its 6 flags are 3% to 7% short,
  no different from the 33 unflagged parameters (up to 6.8%). For scale, an
  MCMC run with an effective sample size of 400 estimates an SD only to
  about 3.5%.
- **0.25 flags every scan-path narrowing of 10% or more** (Weibull, normal
  n = 12, Puromycin, `logv_dem65`), with one false alarm under 5% (`b0` at
  n = 30) and one near miss (`b2` at n = 30, 9.1%).
- **`logv_x2` was flagged by chance.** It is 60% too narrow because of the
  boundary tail of a variance near zero, which bend does not measure; its
  bend of 0.114 sits between `b` (0.108, 6.6% short) and `z_y2_y6` (0.115,
  3.1% short). At 0.25 no check flags it. The fix is a check aimed at that
  failure, such as the slope of the scan at its ends against the Gaussian
  slope (`lab/DESIGN.md`, "Later"), not a bend level kept low to catch it.
