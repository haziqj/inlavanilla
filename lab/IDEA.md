# Task: scaffold and build the R package **INLAvanilla** in haziqj/inlavanilla

## What it is
INLAvanilla runs the INLAvaan approximation pipeline for any user-supplied log-density, with no latent field and nothing SEM-specific. People can try "non-latent INLA" on any model, and it doubles as a well-tested reference engine. Tagline: *INLA-style skew-normal marginals for any log-posterior.*

Keep it MINIMAL. Make the core pipeline correct and well tested before adding any features.

## Primary sources (read before writing maths code)
- INLAvaan source: github.com/haziqj/INLAvaan. Port and adapt the existing implementations of the joint Gaussian, VB mean correction, skew-normal marginal fitting (including the volume correction), and copula sampling. Do not re-derive them from generic INLA knowledge.
- Papers: arXiv 2603.25690 and 2604.00671 (Lemma 3.1, Corollary 3.3, Proposition 3.4 are referenced below).
- If the INLAvaan code and anything in this prompt disagree, flag it and ask. Do not guess.

## User interface (minimum)
The user supplies **one** of:
- `loglik(theta)` + `logprior(theta)`, or
- `logpost(theta)` (unnormalised).

Optionally they also supply `grad` (gradient of the log-posterior; fall back to numerical gradients, e.g. numDeriv), `init` (starting values; required if dimension can't be inferred), and parameter names. They can also declare per-parameter transforms (`"identity"`, `"log"` for positive params, `"tanh"`/artanh for correlations), so the user works on the natural scale while the engine works on the unconstrained scale ϑ.

Main entry point, something like:
```r
fit <- inla_vanilla(logpost = ..., init = ..., grad = NULL, transform = NULL, ...)
```

## Pipeline (mirror the INLAvaan diagram)
1. **1a. Joint Gaussian (Laplace):** optimise to the mode ϑ*, Ω = H(ϑ*)⁻¹. Marginal likelihood: log π(y) ≈ log π(y, ϑ*) − log π̃_G(ϑ* | y). This only applies when loglik + logprior are given separately, or when the user asserts that logpost is properly scaled; document this clearly.
2. **1b. VB mean correction:** ϑ̃ = ϑ* + δ̂, covariance unchanged. Make this optional, on by default if INLAvaan does so.
3. **2. Skew-normal marginals** for j = 1..m:
   - (a) height: scan h_j(t) = log π(ϑ* + t v_j | y) along the CMP direction (Lem. 3.1)
   - (b) volume: efficient correction + t γ′_j (Cor. 3.3)
   - (c) fit SN(ξ̂_j, ω̂_j, α̂_j) to p_j(ϑ_j(t)) ∝ exp(h_j(t)), with ϑ_j(t) = ϑ*_j + t √Ω_jj
   - Back-transform densities via change of variables: log → f(log ψ)/ψ; tanh → f(artanh ρ)/(1 − ρ²).
4. **3. Gaussian copula draws:** R = cor(Ω); z ~ N(0, R); u = Φ(z); ϑ_j = Q_SN^(j)(u_j). Return B draws on both scales.

## Outputs (S3 class, minimal)
- `summary()`: per-parameter mean, sd, quantiles, mode from the SN marginals (natural scale), plus log marginal likelihood
- `marginals()`: density/quantile functions per parameter
- `sample()` / `draws()`: copula draws, optionally `g(theta)` for derived quantities
- `plot()`: marginal densities; overlaying MCMC samples, if supplied, is a nice-to-have
- `print()`

Out of scope for v0.1: factor scores, PPP/DIC, multiple-chain anything, compiled code.

## Tests (testthat) — the main point of the package
Use closed-form or easily computed posteriors as ground truth:
- Normal–normal conjugate (Gaussian approximation should be exact; SN α ≈ 0)
- Normal with unknown variance (log transform; compare to inverse-gamma marginal)
- Beta–binomial on logit scale (skewed; compare SN marginal to exact Beta)
- Gamma–Poisson (log scale)
- A bivariate case with a correlation parameter (tanh transform)
- Check the marginal likelihood against analytic values for conjugate cases
- Check analytic grad vs numerical grad give the same fit
- Check that copula draws reproduce the marginal quantiles and R

Use tolerances appropriate to approximation error, not exact equality. Where an approximation is inexact, test it against a long Stan/MCMC-free reference: numerical quadrature in 1–2D.

## Conventions
- devtools / usethis / testthat (3e), roxygen2 with markdown, MIT or GPL-3 to match INLAvaan
- Base R + minimal Imports (numDeriv, sn if appropriate; avoid heavy deps)
- A short README with one worked example (beta–binomial is good) and one Quarto vignette
- GitHub Actions R-CMD-check
- Small, reviewable commits. Open a draft PR rather than pushing large changes to main.

## Order of work
1. Read the INLAvaan source and the papers; write a short plan mapping INLAvaan functions to INLAvanilla functions, and show it to me before coding.
2. Scaffold the package skeleton.
3. Implement 1a and its tests, then 1b, 2, 3, each with tests before moving on.
4. Add the S3 methods, README, and vignette.