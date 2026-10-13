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

## Later

- LOO/WAIC need pointwise log-likelihoods. This would be an optional extra argument (e.g. `loglik_i`) and does not reopen the `logpost` decision.
- INLAvaan and inlavanilla both export `sampling`, so loading both masks one of them. If INLAvaan later uses inlavanilla as its engine, the generic should live here.
- **Transform layer, as a wrapper around the engine.** Given `transform`, build an unconstrained `logpost(ginv(ϑ)) + log|det J|` and its gradient, run the unchanged v1 engine, then post-process the summaries. Size-1 blocks map the skew-normal back analytically; joint blocks go through the draws. It would support:
  - named elementwise transforms (`"log"`, `"logit"`, `"tanh"`);
  - custom elementwise transforms `list(g, ginv)`;
  - joint blocks `list(pars, type = "cov")` or `list(pars, g, ginv)`, with log|det J| taken numerically when not supplied. A built-in `"cov"`/`"cor"` would follow INLAvaan: log sds plus atanh of (partial) correlations.

  A structured `init` such as `list(Sigma = diag(p))` could sit on top. None of this changes the engine.
- A TMB adaptor: with random effects, TMB's `fn` is already the Laplace-integrated π(y | ϑ).

## API sketch (v0.1)

```r
fit <- inla_vanilla(logpost, init, grad = NULL, vb = TRUE, control = list(), ...)
fit; summary(fit); plot(fit)
sampling(fit, nsamp = 1000, fun = NULL)
```
