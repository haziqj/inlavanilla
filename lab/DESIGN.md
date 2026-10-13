# INLAvanilla design decisions

Decisions agreed after `IDEA.md` was written. Where the two disagree, this file wins.

## Settled

- **One `logpost`.** The user supplies `logpost(theta, ...)` returning log p(y, θ). There is no `loglik` / `logprior` split. The log marginal likelihood is valid only when `logpost` includes every normalising constant, and the docs must say so.
- **Pure R.** No Rcpp and no `src/`. Users bring compiled code by passing any R-callable function, e.g. `Rcpp::cppFunction()`, TMB's `obj$fn` / `obj$gr` (negate them: TMB minimises), RTMB, or BridgeStan.
- **No numDeriv, no sn.** Port INLAvaan's own finite differences (`fast_grad`, `fast_jacobian`, `fast_hessian`, now `fd_*` with per-coordinate steps) and skew-normal code (`fit_skew_normal`, `dsnorm`, `psnorm`, `qsnorm_fast`, `owen_t`). Imports are limited to base packages (stats, graphics, utils) plus cli.
- **v1 has no transform layer.** `logpost` is written on the unconstrained scale ϑ ∈ ℝᵐ. Including any Jacobian is the user's job: a prior on σ with ϑ = log σ needs `+ theta["log_sigma"]`. The docs must stress this, because forgetting it gives the wrong posterior, not merely a worse approximation. `summary()` reports on the ϑ scale. Natural-scale or derived quantities, such as σ, Σ, or continuous-time SEM drift matrices, come from `sampling(fit, fun = ...)`, and the user summarises the draws. `summary()` takes no `fun` argument for now. Copula draws never call `logpost`, so `nsamp` can be large at almost no cost. The reason for this choice: real constraints are often joint and model-specific (covariance matrices, stable drift matrices), and no generic transform layer can express them.
- **S3, not S7.** S4 and R6 are ruled out. The fit is a classed list built by one internal constructor, `new_inla_vanilla()`. `summary()` returns a `summary.inla_vanilla` object with its own `print()` method. S7 would add an Imports entry and an `@` workaround for R < 4.3, for features we don't use (validated properties, inheritance, multiple dispatch).
- **cli is a required Imports entry.** It is used for formatted messages, warnings and errors (`cli_abort()` etc.) and for progress bars and steps, as in INLAvaan. It may be dropped at a later optimisation stage, but not before. Add it to `DESCRIPTION` with the first code that uses it, to avoid an unused-Imports NOTE.
- **VB mean correction on by default.** It feeds the log marginal likelihood (+½ δᵀHδ, as in INLAvaan) and later LOO/WAIC. It is also a diagnostic of the skew-normal fit, because skew-normal means should agree with ϑ̃ (Prop. 3.4). In INLAvaan the skew-normal scan is centred at the mode ϑ*, so the VB shift does not move the marginals themselves.
- **Stages 2 and 3 are independent switches** (agreed after v0.1.0): `vb = TRUE` and `marginals = c("skewnormal", "gaussian")`. They are not a ladder: the skew-normal scans are centred at the mode, so VB never moves them, and skew-normal without VB is the cheap path when there is no gradient (VB costs (m + 1)(2m² + 1) calls per sweep without one, about 230,000 calls at m = 30 against 6,900 for the rest of the fit). `vb = FALSE, marginals = "gaussian"` is the Laplace approximation alone; `marginals = "gaussian"` alone is the VB-centred Gaussian. The value may be abbreviated (`"skew"`, `"gauss"`), as with `match.arg()`. Gaussian marginals are stored as skew-normals with α = 0, so every method works unchanged, and `sampling()` then draws from the joint Gaussian exactly.
- **Three output levels, `verbose = FALSE` and `silent = FALSE` by default** (agreed after v0.1.0: `verbose` is a full diagnostic of the fit, and `silent` silences everything, for simulations). The default is the old `verbose = TRUE` output, one line per stage, and must stay as fast. Progress bars respect the user's `cli.progress_show_after` (0.5 s when unset); rendered documents keep them out through options in their setup chunk, as INLAvaan's README does, rather than through a terminal check in the package. `verbose = TRUE` adds grey detail lines under each stage, from what the fit has already computed (no extra calls): nlminb's message, iterations and evaluations and whether it was restarted, the whitened gradient, how the Hessian was taken, the largest correlation and the smallest eigenvalue of the correlation matrix, the VB steps per iteration, Anderson and fallback use, a table of the marginal diagnostics with the `summary()` flags marked, and the calls per stage. nlminb's own trace stays with `control = list(nlminb = list(trace = 1))`. `silent = TRUE` prints nothing and muffles every warning and message raised during the fit, including those of `logpost` and `grad`; errors still stop it. Both `TRUE` is an error. Whatever the level, each distinct warning is recorded as plain text in `fit$warnings`, and the fit stores `verbose` and `silent`, so later methods (a planned `loo()`) can follow `fit$silent`.
- **Posterior draws use `sampling()` only**, matching INLAvaan's exported generic. There is no `draws()` alias, and no `sample()` method, since `base::sample` is not generic.
- **One method per stage**, keeping INLAvaan's defaults:
  - skew-normal marginals
  - `nlminb`
  - the `"shortcut"` volume correction
  - for the VB expectation, the fully symmetric degree-5 cubature rule (2m² + 1 nodes on the principal axes of Ω), replacing INLAvaan's Sobol points. It is exact when the log-posterior is a degree-6 polynomial in whitened coordinates, and deterministic. In `lab/vb-benchmark.R` it was the most accurate rule on every model tested, from 1 to 30 parameters: 10 to 150 times better than 64 Sobol points in 1–2D, and 0.016 against 0.058 SDs at m = 30. The cost is 2m² + 1 gradient calls per sweep, so the docs advise supplying `grad` for larger m.
- **Not ported:**
  - the `asymgaus`, `marggaus` and `sampling` marginals
  - NORTA
  - the `ucminf` and `optim` optimisers
  - Gauss–Hermite VB
  - parallel cores
  - PPP, DIC, LOO and WAIC (v0.1)

## Open

Nothing at present.

## Findings to discuss

- **Linear scan path at small n.** The skew-normal scan follows the conditional mode path (CMP) of the Laplace Gaussian, which is a straight line (Lemma 3.1). When another parameter's conditional mode bends along that path, the scanned heights are too low in the tails. A common case is a noise scale, where σ̂²(K) ∝ SSR(K) bends quadratically. The marginal then comes out too narrow, by a factor that shrinks as O(1/n).
  - Example: Michaelis–Menten on the treated `Puromycin` data (n = 12; parameters log Vm, log K, log σ).
  - Exact marginal SDs from a 3D grid: 0.0382 and 0.1537 for log Vm and log K.
  - Skew-normal: 0.0315 and 0.127, both about 18% too narrow. Laplace is the same.
  - A full Laplace marginal with exact conditional modes and conditional Hessians recovers 0.0383 and 0.1541. So the gap comes from the linear path, not from the skew-normal fit or the volume correction.
  - log σ itself is captured almost exactly: skew-normal mean 2.4733 against an exact 2.4734, while Laplace gives 2.3324.
  - The test suite found the same effect independently:
    - For the normal mean under a normal–inverse-gamma prior, the skew-normal SD is short by exactly √(2(aₙ − 1)/(2aₙ + 1)), a Gaussian against the true Student-t.
    - For the log scale of a Weibull model (n ≈ 15), both the skew-normal and Laplace SDs are about 11% too small, and a full Laplace marginal is within 0.1%.
  - Prototype remedies on Puromycin (`lab/scan-path-prototype.R`), SD error for log Vm and log K:

    | Scan | SD error |
    |---|---|
    | Current | −17.7% |
    | One Newton step towards the conditional mode (Laplace curvature) | −8.1% |
    | Three Newton steps | −5.7% |
    | Full Laplace: exact conditional mode and conditional log-determinant at each scan point | −1.0% |

    So the bent path accounts for about two thirds of the gap, and the linear volume correction for the rest.
  - **Decision.** Keep the CMP and the shortcut volume correction, which make the method fast. Document the limit in the Accuracy section of `?inla_vanilla`, and flag it with the cheap `bend` diagnostic: two gradient calls per parameter, flagged above 0.1 nats.
  - The full-Laplace remedy is INLA's "Laplace strategy" for marginals. It needs a conditional optimisation and a conditional Hessian at each of the 21 scan points: about 21(m + 3) gradient calls per parameter, against 2m + 2 now, so roughly 10× the cost of the marginal stage. It is a change of method and needs a decision.
  - In this example the SN − VB diagnostic flagged log σ (0.13 SDs). That gap came from the VB mean being off, not the skew-normal, so the diagnostic is documented as "either may be off".

## Later

- **Heavy-tailed marginals on the SEM (known limitation, kept for now).** Recorded after v0.1.0 for future work. On the Political Democracy SEM of `vignettes/articles/rcpp.qmd` (39 parameters, n = 75), checked against HMC (4 chains of 10,000 draws, R-hat 1.000, tail ESS over 8,000; `lab/sem-hmc.R`) and against full Laplace marginals on a wide grid (`lab/sem-marginal-check.R`), which agree with each other:
  - `logv_x2` and `logv_dem65` are badly off. The skew-normal SD is 0.38 and 0.46 of the posterior SD, and the 2.5% quantile is −3.46 against −6.98 and −3.60 against −7.34 (2.4 and 2.0 posterior SDs). The upper quantiles agree. The data cannot keep these residual variances away from zero, so the log-density's left tail is linear with slope 0.50: the Gamma(1, 0.5) prior on the SD, with its Jacobian. No skew-normal can follow an exponential tail. On the variance scale the skew-normal overstates the mean by 16–18% and puts the 2.5% quantile at 0.03 instead of 0.001.
  - The other 37 parameters: SD ratios 0.93–1.00 (median 0.97), means within 0.08 posterior SDs, quantiles within 0.21. This is the straight-scan narrowing below.
  - On the natural scale (INLAvaan reports variances), against Stan NUTS on an independently written model (`lab/poldem.stan`, 4 × 10,000 draws, R-hat 1.00, no divergences): x2~~x2 mean 0.23 posterior SDs too high, SD ratio 0.85, 2.5% quantile 0.031 against 0.001; dem65~~dem65 mean 0.17 SDs too high, SD ratio 0.90. The other 37 natural-scale parameters: means within 0.09 SDs, SD ratios 0.93–1.01, quantiles within 0.28 SDs. So the error is large on the unconstrained scale that inlavanilla reports, and moderate on the natural scale that INLAvaan reports.
  - INLAvaan (0.3.2.9006, `asem(..., meanstructure = TRUE)`, default priors) gives the same marginals as inlavanilla to about 3 significant figures (x2~~x2 mean 0.1322 against 0.1323, 2.5% quantile 0.0314 against 0.0315), so matching INLAvaan does not rule this out; both share the limitation.
  - Diagnostics: `bend` (> 0.1) flags both bad parameters (0.11 and 0.44), but its size does not rank how far off a marginal is. `tails` passes everywhere, because at ±4 Laplace SDs the drop is still over 3 nats: the tail flattens only further out.
  - An earlier comparison against full Laplace marginals on a ±4.5 SD grid gave "about 25% too narrow" for `logv_dem65`: the grid cut off the tail. Use a wide grid.
  - Candidate fixes: a tail check from the slope of the scan at its ends (no extra calls), compared with the Gaussian slope; a wider or adaptive scan for flagged parameters; the full Laplace marginal on an adaptive grid for flagged parameters only (about 10 times the cost of the marginal stage); one Newton step towards the conditional mode at each scan point (halved the narrowing on Puromycin, `lab/scan-path-prototype.R`).
- **LOO/WAIC (planned for v1.1).** Port INLAvaan's deterministic second-order Taylor LOO and WAIC (Alhyari, Jamil, Montcho and Rue, 2026, "Deterministic leave-one-cluster-out cross-validation for multilevel Bayesian models", arXiv:2609.00670; `R/loo.R` and `R/method-loo.R` in INLAvaan): per-unit score s_u and Hessian H_u of the pointwise log-likelihood at the VB mean, with the Laplace covariance Ω; log CPO_u = ℓ_u − ½ s_uᵀ(Ω⁻¹ + H_u)⁻¹ s_u + ½ log|I + ΩH_u|, lpd_u likewise with −H_u, p_waic,u = s_uᵀΩs_u + ½ tr[(H_uΩ)²]. No draws and no Monte Carlo error. It needs a user function `loglik(theta, ...)` returning the vector of per-unit log-likelihoods; where it is supplied (`loo(fit, loglik, ...)` or `inla_vanilla(..., loglik = )`) is still to be decided. Shaping the results like `loo` package objects would let `loo::loo_compare()` compare models without a dependency. The sampling route already works today: `loo::loo(sampling(fit, 4000, fun = loglik, ...))`.
- INLAvaan and inlavanilla both export `sampling`, so loading both masks one of them. If INLAvaan later uses inlavanilla as its engine, the generic should live here.
- **Transform layer, as a wrapper around the engine.** Given `transform`, build an unconstrained `logpost(ginv(ϑ)) + log|det J|` and its gradient, run the unchanged v1 engine, then post-process the summaries. Size-1 blocks map the skew-normal back analytically; joint blocks go through the draws. It would support:
  - named elementwise transforms (`"log"`, `"logit"`, `"tanh"`);
  - custom elementwise transforms `list(g, ginv)`;
  - joint blocks `list(pars, type = "cov")` or `list(pars, g, ginv)`, with log|det J| taken numerically when not supplied. A built-in `"cov"`/`"cor"` would follow INLAvaan: log sds plus atanh of (partial) correlations.

  A structured `init` such as `list(Sigma = diag(p))` could sit on top. None of this changes the engine.
- A TMB adaptor: with random effects, TMB's `fn` is already the Laplace-integrated π(y | ϑ).

## API sketch (v0.1)

```r
fit <- inla_vanilla(logpost, init, grad = NULL, ..., vb = TRUE,
                    marginals = c("skewnormal", "gaussian"), verbose = FALSE,
                    silent = FALSE, control = list())
fit; summary(fit); plot(fit)
sampling(fit, nsamp = 1000, fun = NULL)
```
