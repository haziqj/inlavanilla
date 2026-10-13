# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project status

`inlavanilla` (styled "INLAvanilla" in prose) is a small, pure-R package that exports `inla_vanilla()`, `sampling()`, `loo()` and `waic()`, with S3 `print`, `summary` and `plot` methods. The original specification is `lab/IDEA.md`, and `lab/DESIGN.md` records the decisions made since. Read both before starting any feature work. Where they disagree, `DESIGN.md` wins. Together they define the API, the maths pipeline and the test plan. Test posteriors with a natural-scale closed form (e.g. the beta–binomial) are written as unconstrained `logpost`s with the Jacobian included, and compared through `sampling(fit, fun = ...)` or the ϑ-scale density. `lab/` is excluded from the package build (`.Rbuildignore`), so it is a safe place for design notes.

Tagline: *INLA-style skew-normal marginals for any log-posterior.* The package runs the INLAvaan approximation pipeline on any user-supplied log-density, with no latent field and nothing SEM-specific. The README and vignette open with the package owner's wording: "If you can write down the log-likelihood plus the log-prior, `{inlavanilla}` gives you approximate Bayesian inference for the parameters of your model."

## Commands

Standard devtools workflow, run from the repo root in R:

```r
devtools::load_all()                 # load the package for interactive work
devtools::document()                 # regenerate man/*.Rd and NAMESPACE from roxygen
devtools::test()                     # run all testthat tests
devtools::test(filter = "laplace")   # run tests/testthat/test-laplace.R only
testthat::test_file("tests/testthat/test-laplace.R")
devtools::check()                    # full R CMD check (what CI runs)
devtools::build_readme()             # re-knit README.Rmd -> README.md
pkgdown::build_site()                # local pkgdown site (docs/ is gitignored)
```

Format with Air (`air format .`). The `format-suggest` workflow runs Air on every PR and posts suggestions, so format before pushing, but don't restyle code unrelated to your change.

Generated files that must not be edited by hand: `NAMESPACE` and `man/*.Rd` (roxygen2 with markdown, edit the `#'` comments instead) and `README.md` (edit `README.Rmd`, then rebuild and commit any new `man/figures/` images).

## CI

- `R-CMD-check.yaml` runs on macOS, Windows, and Ubuntu with R devel, release and oldrel-1 through oldrel-4. Avoid base-R features that are too new for oldrel-4.
- `test-coverage.yaml` runs covr and uploads to Codecov. The Codecov checks are informational only.
- `pkgdown.yaml` builds the site and deploys it to `gh-pages` on pushes to main.

## Architecture

Entry point: `inla_vanilla(logpost, init, grad = NULL, ..., vb = TRUE, marginals = c("skewnormal", "gaussian"), verbose = FALSE, silent = FALSE, control = list())` in `R/inla-vanilla.R`. As in `optim()` and `nlminb()`, the dots come after the gradient, so `vb`, `marginals`, `verbose`, `silent` and `control` must be named in full and data arguments cannot partially match them. `marginals` is checked by `check_choice()`, which accepts unambiguous abbreviations in any case (`"skew"`, `"Gauss"`), as `match.arg()` does, but raises the package's own error. `inla_vanilla()` checks the arguments and sets up the condition handlers; the stages run in the internal `vanilla_fit()`.

- **Input.** The user supplies a single `logpost(theta, ...)` on the unconstrained scale ϑ ∈ ℝᵐ, with any Jacobian already included; `...` is passed through as in `optim()`.
- **No transform layer in v1.** Natural-scale and derived quantities come from `sampling(fit, fun = ...)`. A later transform layer is planned as a wrapper that builds an unconstrained `logpost` for the unchanged engine (see `lab/DESIGN.md`), so keep the engine free of transform logic.
- **Wrappers.** `inla_vanilla()` binds `...` once in `logpost_x()` and `grad_x()`, and `vanilla_fit()` wraps them in closures (`lp`, `gr`) that name the parameters, count calls, record the last parameter value and tick progress. Every stage receives these closures, never the user functions. Never forward `...` to internal helpers, because a user's data argument could partially match one of their formals (the same applies to `fun` in `sampling()`). `with_context()` reports an error from `logpost` or `grad` with the stage and the parameter value.

The stages run in order, one file each:

1. **Laplace** (`R/laplace.R`, `laplace_fit()`).
   - The mode comes from `nlminb`. If the whitened gradient at the mode is still above 1e-3, `nlminb` reruns with `scale = sqrt(diag(H))`.
   - Finite-difference steps are set per coordinate by `fd_scales()`, a cheap curvature search. With a gradient the Hessian takes central differences of it at 1e-4 of that scale; without one, `fd_hessian()` at 0.01 of it. This keeps posterior SDs from 1e-6 to 1e4 working.
   - It returns R (upper Cholesky of the precision), Σ, L (L Lᵀ = Σ), the whitened gradient at the mode (which, not `nlminb`'s return code, decides the convergence warning), the Jacobian asymmetry (a check on a user `grad`), and the Laplace log marginal likelihood.
2. **VB** (`R/vb.R`).
   - `vb_correct()` is INLAvaan's fixed-point iteration with Anderson(1) acceleration and an `nlminb` fallback.
   - The Gaussian expectation uses `vb_nodes()`, the fully symmetric degree-5 cubature rule (2m²+1 nodes on Σ's principal axes). This replaced INLAvaan's Sobol points after the benchmark in `lab/vb-benchmark.R`.
   - `fit$logml` = Laplace + ½ δᵀHδ.
3. **Skew-normal marginals** (`R/marginals.R`, `fit_marginals()`).
   - Each parameter is scanned on `seq(-4, 4)` posterior SDs along the conditional mode path, corrected by `gamma1_*`, and fitted by `fit_skew_normal()`.
   - The volume correction has two routes that must agree: `gamma1_grad()` is INLAvaan's shortcut, using the gradient; `gamma1_lp()` uses second differences of `logpost` when there is no gradient, at O(m) calls per parameter instead of O(m²).
   - The marginals are scanned around the mode ϑ*, so VB does not move them.
   - `scan_bend()` computes the diagnostic `bend`: the log-density, in nats, that the straight scan misses at ±2 SDs, from one Newton step in the other parameters. `summary()` flags it above 0.1. It catches the documented small-n narrowing.
   - `marginals = "gaussian"` skips this stage: `sn` holds the Gaussian marginals N(mean_vb, Σ_jj) as skew-normals with `alpha = 0`, the scan diagnostics are `NA`, and `print()`, `summary()`, `plot()` and `sampling()` branch on `is_gaussian()` or `alpha == 0`. VB and the marginals are independent switches, not a ladder: VB never moves the skew-normal marginals.
   - The scan is a fixed ±4 SDs, as in INLAvaan. `tails` records a scan whose ends stay above `log_threshold`. `scan_end_mass()` gives INLAvaan's `scan_end_mass`, the fitted mass outside that window; `summary()` flags it above 0.05.
   - `R/diagnostics.R`: `warn_diagnostics()` raises one `inlavanilla_diagnostics_warning` at the end of `vanilla_fit()` (before the fit object is built, so `fit$warnings` records it) when a marginal is clearly off: heavy `tails`, NMAD > 0.1, `scan_end_mass` > 0.05, VB shift > 1 SD or the VB fallback (`warn_levels`, INLAvaan's levels). `bend` and the SN-VB gap are flagged in `summary()` only. Healthy fits in every document must stay silent. Extending the scan into heavy tails was tried and rejected, because it hurt the fit in other cases.
4. **Copula draws** (`R/sampling.R`). `sampling()` is an S3 generic, matching INLAvaan's exported `sampling()`. Draws go through `sn_score_map()`, a spline from normal scores to skew-normal quantiles; do not use `qsnorm_fast()` directly, because it is only a starting value, accurate to about 5e-2.

Supporting files:

- **Skew-normal utilities** (`R/skew-normal.R`, `R/qsnorm-fast.R`): ported from INLAvaan. The Apache-licensed `qsnorm_fast()` keeps its attribution header. All of them are internal.
- **Methods** (`R/methods.R`): `print()` and `summary()` print via `cli::cat_*` to stdout.
- **LOO and WAIC** (`R/loo.R`): `loo(fit, loglik, ...)` and `waic()`, the deterministic second-order Taylor method of Alhyari, Jamil, Montcho and Rue (2026, arXiv:2609.00670), as in INLAvaan: no refits and no draws.
  - `loglik(theta, ...)` returns one log-likelihood term per unit; `logpost` is their sum plus the log-prior, and the fit never sees `loglik`.
  - `taylor_pointwise()` takes each unit's score and Hessian by central differences in the whitened coordinates of N(mean_vb, Σ), 0.01 SD apart: m² + m + 1 calls of `loglik`. If any unit's leverage k (an eigenvalue of −ΣH_u) reaches 1, every unit falls back to first order, with a warning; there is no mixing of orders.
  - The results have the loo package's layout (classes `loo` and `waic`, `estimates`, `pointwise` with one `elpd_*` column), so `loo::loo_compare()` works with no dependency. The methods are also registered on loo's generics (`loo_x_inla_vanilla()`, `waic_x_inla_vanilla()`), and `loo.default()` and `waic.default()` forward other objects to loo when it is installed.
  - The generics take `object` first, not `x`, so data named `x` reaches `loglik`. They follow `fit$silent`, for `loglik`'s own warnings too, and record each distinct warning in the result's `warnings`.
- **Progress** (`R/progress.R`): three output levels, held in the `prog` environment from `progress_new()`.
  - Default (`verbose = FALSE, silent = FALSE`): each stage has a cli progress bar that appears after the user's `cli.progress_show_after` seconds, or 0.5 s when that option is unset (never override a value the user set: setup chunks of rendered documents use it to keep bars out), and ends with a plain styled `message()` line. cli's glue-based helpers (`cli_progress_step()`, `cli_alert_*()`) cost several ms per call, as much as a whole small fit, so avoid them in the fitting path, and keep the default path free of anything built only for `verbose`.
  - `verbose = TRUE`: the `report_*()` functions add indented grey detail lines under each stage, behind `if (prog$detail)`. They use only what the fit has already computed (`lap$restart`, `lap$hessian`, `vbc$history`, the call counters), never further calls of `logpost` or `grad`.
  - `silent = TRUE`: nothing is printed, and the handlers in `inla_vanilla()` muffle every warning and message raised during the fit. `verbose = TRUE, silent = TRUE` is an error.
  - Whatever the level, the handler records each distinct warning in `prog$warnings`, which becomes `fit$warnings`, as plain text: `warn()` attaches its message unwrapped and uncoloured as `text`. The fit also stores `verbose` and `silent`.
- **Errors and warnings**: raise them with `abort()` and `warn()` (`R/utils.R`), which format with cli and signal with base R, giving the classes `inlavanilla_error` and `inlavanilla_warning`. Do not use `cli::cli_abort()` or `cli::cli_warn()`: they call rlang, which cli only suggests.
- **Non-ASCII output**: user-visible non-ASCII characters (σ, ·) go through `sd_unit()` and `sep_dot()`, which fall back to ASCII. R source files must stay ASCII.

Out of scope for v0.1: factor scores, PPP and DIC, multiple chains, compiled code, and INLAvaan's alternative methods (other marginal types, NORTA, other optimisers, Sobol or Gauss–Hermite VB, parallel cores). Do not add alternative methods; there is one method per stage.

## Source-of-truth rules for the maths

- Port and adapt the implementations in the INLAvaan source (github.com/haziqj/INLAvaan): joint Gaussian, VB correction, SN fitting with the volume correction, and copula sampling. Do not re-derive them from generic INLA knowledge. The references are the papers arXiv 2603.25690 and 2604.00671.
- If INLAvaan's code and `lab/IDEA.md` / `lab/DESIGN.md` disagree, flag the disagreement and ask the user. Do not guess.
- The pipeline to lift is `R/inlavaan.R` from mode-finding through copula sampling. It reaches the model only through the `joint_lp()` / `joint_lp_grad()` closures and the per-parameter transform functions `g`, `ginv`, `g_prime` and `ginv_prime` (`R/partable.R`). Its helpers are all pure R: `R/utils-optim.R`, `R/utils-inlavaan.R` (`compute_gamma1j`), `R/skew-normal.R`, `R/qsnorm_fast.R`, `R/sobol.R` (with the `SobolOwen` table in `R/sysdata.rda`), `R/posterior-marginals.R` and `R/posterior-sampling.R`.

## Testing philosophy

The tests are described as "the main point of the package." They use closed-form or cheaply computed posteriors as ground truth: normal–normal conjugate, normal with unknown variance (log), beta–binomial on the logit scale, gamma–Poisson (log), a bivariate case with a correlation (tanh), analytic marginal likelihoods, analytic vs numerical gradients, and copula draws reproducing the marginal quantiles and R. Use tolerances that fit the approximation error, not exact equality. Where there is no closed form, compare against 1–2D numerical quadrature rather than MCMC.

## Dependencies and housekeeping

- Imports should stay limited to base packages (stats, graphics, utils) plus cli, which the user requires for messages and progress bars. Do not add numDeriv, sn or Rcpp.
- The license is GPL (>= 3).
- Add a `NEWS.md` bullet for user-facing changes, as `CONTRIBUTING.md` asks.
- Write references and in-text citations in APA 7 style, by hand, in roxygen, the README and the vignettes: "Jamil, H., & Rue, H. (2026a). *Title in sentence case* (arXiv:2603.25690 [stat.ME]). arXiv. <https://doi.org/...>" (the arXiv identifier and primary category in brackets after the title, as Zotero's APA style gives; escape the square brackets as `\[stat.ME\]` in roxygen), cited as "(Jamil & Rue, 2026a)" or "Jamil and Rue (2026a)". The a/b suffixes appear only where both 2026 papers are cited. In roxygen use `--` for page ranges (it renders as an en dash) and `<https://doi.org/...>` rather than `\doi{}`. `DESCRIPTION` keeps CRAN's `<doi:...>` form.
- The vignette (`vignettes/inlavanilla.qmd`) is Quarto, as in INLAvaan; CI installs Quarto, and R-CMD-check installs the package itself so that Quarto's R session finds it on Windows.

## Coding principles

- **tidyverse style.** Follow the [tidyverse style guide](https://style.tidyverse.org). Format R code with Air and respect its `# fmt: skip` directives. Use descriptive `snake_case` names.
- **tinyverse: no unnecessary dependencies.** Prefer base R. A new `Imports` entry must be justified by substantial functionality that would be unreasonable to replicate, and must not bring in long dependency chains.
- **Works out of the box everywhere.** The package must install cleanly on macOS, Linux and Windows with no required system libraries. Optional system libraries (OpenMP, OpenBLAS) may be used for performance but must degrade gracefully. Never assume Homebrew or apt is available.

## Workflow conventions

Make small, reviewable commits. Open a draft PR rather than pushing large changes to main.
