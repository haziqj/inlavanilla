# INLAvanilla design decisions

Decisions agreed after `IDEA.md` was written. Where the two disagree, this file wins.

## Settled

- **One `logpost`.** The user supplies `logpost(theta, ...)` on the natural scale, returning log p(y, θ). There is no `loglik` / `logprior` split. The log marginal likelihood is valid only when `logpost` includes every normalising constant, and the docs must say so.
- **Pure R.** No Rcpp and no `src/`. Users bring compiled code by passing any R-callable function, e.g. `Rcpp::cppFunction()`, TMB's `obj$fn` / `obj$gr` (negate them: TMB minimises), RTMB, or BridgeStan.
- **No numDeriv, no sn.** Port INLAvaan's own finite differences (`fast_grad`, `fast_jacobian`, `fast_hessian`) and skew-normal code (`fit_skew_normal`, `dsnorm`, `psnorm`, `qsnorm_fast`, `owen_t`). Imports are limited to base packages (stats, graphics, utils) plus cli.
- **Transforms** cover identity, log, logit and tanh. The interface is still open (see below).
- **S3, not S7.** S4 and R6 are ruled out. The fit is a classed list built by one internal constructor, `new_inla_vanilla()`. `summary()` returns a `summary.inla_vanilla` object with its own `print()` method. S7 would add an Imports entry and an `@` workaround for R < 4.3, for features we don't use (validated properties, inheritance, multiple dispatch).
- **cli is a required Imports entry.** It is used for formatted messages, warnings and errors (`cli_abort()` etc.) and for progress bars and steps, as in INLAvaan. It may be dropped at a later optimisation stage, but not before. Add it to `DESCRIPTION` with the first code that uses it, to avoid an unused-Imports NOTE.
- **VB mean correction on by default.** It feeds the log marginal likelihood (+½ δᵀHδ, as in INLAvaan) and later LOO/WAIC. It is also a diagnostic of the skew-normal fit, because skew-normal means should agree with ϑ̃ (Prop. 3.4). In INLAvaan the skew-normal scan is centred at the mode ϑ*, so the VB shift does not move the marginals themselves.
- **One method per stage**, keeping INLAvaan's defaults:
  - skew-normal marginals
  - `nlminb`
  - the `"shortcut"` volume correction
  - Sobol QMC for VB
- **Not ported:**
  - the `asymgaus`, `marggaus` and `sampling` marginals
  - NORTA
  - the `ucminf` and `optim` optimisers
  - Gauss–Hermite VB
  - parallel cores
  - PPP, DIC, LOO and WAIC (v0.1)

## Open

- Transform interface: named transforms (`transform = c(sigma = "log")`) vs bounds (`lower = c(sigma = 0)`).
- Posterior draws: `sampling()` to match INLAvaan's exported generic (`base::sample` is not generic, so `sample()` would mask base R). Is a `draws()` alias wanted?

## Later

- LOO/WAIC need pointwise log-likelihoods. This would be an optional extra argument (e.g. `loglik_i`) and does not reopen the `logpost` decision.
- INLAvaan and inlavanilla both export `sampling`, so loading both masks one of them. If INLAvaan later uses inlavanilla as its engine, the generic should live here.
- A TMB adaptor: with random effects, TMB's `fn` is already the Laplace-integrated π(y | ϑ).

## API sketch (v0.1)

```r
fit <- inla_vanilla(logpost, init, grad = NULL, <transforms>, vb = TRUE,
                    control = list(), ...)
fit; summary(fit); plot(fit)
sampling(fit, nsamp = 1000, fun = NULL)
```
