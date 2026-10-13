# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project status

`inlavanilla` (styled "INLAvanilla" in prose) is a freshly scaffolded R package (usethis template). There is no `R/` directory, no tests and no exported functions yet. The original specification is `lab/IDEA.md`, and `lab/DESIGN.md` records the decisions made since. Read both before starting any feature work. Where they disagree, `DESIGN.md` wins. Together they define the API, the maths pipeline, the test plan and the required order of work. Test posteriors with a natural-scale closed form (e.g. the beta–binomial) are written as unconstrained `logpost`s with the Jacobian included, and compared through `sampling(fit, fun = ...)` or the ϑ-scale density. `lab/` is excluded from the package build (`.Rbuildignore`), so it is a safe place for design notes.

Tagline: *INLA-style skew-normal marginals for any log-posterior.* The package runs the INLAvaan approximation pipeline on any user-supplied log-density, with no latent field and nothing SEM-specific.

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

## Planned architecture (from `lab/IDEA.md` and `lab/DESIGN.md`)

Main entry point: `inla_vanilla(logpost, init, grad = NULL, vb = TRUE, control = list(), ...)`. The user supplies a single `logpost(theta, ...)` on the unconstrained scale ϑ ∈ ℝᵐ, with any Jacobian already included, and `...` is passed through as in `optim()`. v1 has no transform layer: natural-scale and derived quantities come from `sampling(fit, fun = ...)`. A later transform layer is planned as a wrapper that builds an unconstrained `logpost` for the unchanged engine (see `lab/DESIGN.md`). Keep the engine free of transform logic. When `grad` is missing, gradients fall back to INLAvaan's own finite differences, not numDeriv.

The pipeline has four stages. Implement them in this order, each with tests passing before moving on:

1. **1a. Joint Gaussian (Laplace).** Find the mode ϑ* and set Ω = H(ϑ*)⁻¹. The log marginal likelihood is log π(y, ϑ*) − log π̃_G(ϑ* | y). It is only meaningful when `logpost` includes every normalising constant, and the docs must say so. With VB on, INLAvaan adds ½ δᵀHδ.
2. **1b. VB mean correction.** ϑ̃ = ϑ* + δ̂, with the covariance unchanged. It is on by default. In INLAvaan's skew-normal path it does not move the marginals, which are scanned around ϑ*. It feeds the log marginal likelihood and serves as a check on the skew-normal fit (Prop. 3.4).
3. **2. Skew-normal marginals per parameter.** Scan the height along the CMP direction (Lemma 3.1), add the volume correction (Corollary 3.3), then fit SN(ξ, ω, α). In v1 the marginals stay on the ϑ scale; back-transforming densities belongs to the later transform layer.
4. **3. Gaussian copula draws.** Set R = cor(Ω), draw z ~ N(0, R), take u = Φ(z) and apply the SN quantile per margin. Return draws on the ϑ scale, or of `fun(theta)` when `fun` is given.

Outputs: `print()`, `summary()` (ϑ-scale moments, quantiles and mode, plus the log marginal likelihood), `plot()`, and `sampling()` for copula draws, optionally of a derived `g(theta)`. The name `sampling()` matches INLAvaan's exported generic; `base::sample` is not generic. The object system is S3 (not S7, S4 or R6), and the fit is a classed list built by `new_inla_vanilla()`.

Out of scope for v0.1: factor scores, PPP/DIC/LOO/WAIC, multiple chains, compiled code, and INLAvaan's alternative methods (other marginal types, NORTA, other optimisers, Gauss–Hermite VB, parallel cores).

## Source-of-truth rules for the maths

- Port and adapt the implementations in the INLAvaan source (github.com/haziqj/INLAvaan): joint Gaussian, VB correction, SN fitting with the volume correction, and copula sampling. Do not re-derive them from generic INLA knowledge. The references are the papers arXiv 2603.25690 and 2604.00671.
- If INLAvaan's code and `lab/IDEA.md` / `lab/DESIGN.md` disagree, flag the disagreement and ask the user. Do not guess.
- The pipeline to lift is `R/inlavaan.R` from mode-finding through copula sampling. It reaches the model only through the `joint_lp()` / `joint_lp_grad()` closures and the per-parameter transform functions `g`, `ginv`, `g_prime` and `ginv_prime` (`R/partable.R`). Its helpers are all pure R: `R/utils-optim.R`, `R/utils-inlavaan.R` (`compute_gamma1j`), `R/skew-normal.R`, `R/qsnorm_fast.R`, `R/sobol.R` (with the `SobolOwen` table in `R/sysdata.rda`), `R/posterior-marginals.R` and `R/posterior-sampling.R`.
- The first step in `lab/IDEA.md` is a written plan mapping INLAvaan functions to INLAvanilla functions, shown to the user before any coding.

## Testing philosophy

The tests are described as "the main point of the package." They use closed-form or cheaply computed posteriors as ground truth: normal–normal conjugate, normal with unknown variance (log), beta–binomial on the logit scale, gamma–Poisson (log), a bivariate case with a correlation (tanh), analytic marginal likelihoods, analytic vs numerical gradients, and copula draws reproducing the marginal quantiles and R. Use tolerances that fit the approximation error, not exact equality. Where there is no closed form, compare against 1–2D numerical quadrature rather than MCMC.

## Known scaffold inconsistencies

- `DESCRIPTION` still has template placeholders (Title, Description, Authors@R).
- Imports should stay limited to base packages (stats, graphics, utils) plus cli, which the user requires for messages and progress bars. Do not add numDeriv, sn or Rcpp.
- The license is GPL (>= 3).
- `CONTRIBUTING.md` asks for a `NEWS.md` bullet for user-facing changes, but there is no `NEWS.md` yet.

## Coding principles

- **tidyverse style.** Follow the [tidyverse style guide](https://style.tidyverse.org). Format R code with Air and respect its `# fmt: skip` directives. Use descriptive `snake_case` names.
- **tinyverse: no unnecessary dependencies.** Prefer base R. A new `Imports` entry must be justified by substantial functionality that would be unreasonable to replicate, and must not bring in long dependency chains.
- **Works out of the box everywhere.** The package must install cleanly on macOS, Linux and Windows with no required system libraries. Optional system libraries (OpenMP, OpenBLAS) may be used for performance but must degrade gracefully. Never assume Homebrew or apt is available.

## Workflow conventions

Make small, reviewable commits. Open a draft PR rather than pushing large changes to main.
