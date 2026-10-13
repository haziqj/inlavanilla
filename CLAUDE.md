# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project status

`inlavanilla` (styled "INLAvanilla" in prose) is a freshly scaffolded R package (usethis template). There is no `R/` directory, no tests and no exported functions yet. The full specification lives in `lab/IDEA.md`. Read it before starting any feature work: it defines the API, the maths pipeline, the test plan and the required order of work. `lab/` is excluded from the package build (`.Rbuildignore`), so it is a safe place for design notes.

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

## Planned architecture (from `lab/IDEA.md`)

Main entry point: `inla_vanilla(logpost = , init = , grad = NULL, transform = NULL, ...)`. The user supplies either `loglik` + `logprior` or a single unnormalised `logpost`. Per-parameter transforms (`"identity"`, `"log"`, `"tanh"`) let the user work on the natural scale while the engine works on the unconstrained scale ϑ. Gradients fall back to numerical ones (numDeriv) when `grad` is missing.

The pipeline has four stages. Implement them in this order, each with tests passing before moving on:

1. **1a. Joint Gaussian (Laplace).** Find the mode ϑ* and set Ω = H(ϑ*)⁻¹. The log marginal likelihood is log π(y, ϑ*) − log π̃_G(ϑ* | y). It is only meaningful when `loglik`/`logprior` are given separately or `logpost` is asserted to be properly normalised, and the docs must say so.
2. **1b. VB mean correction.** ϑ̃ = ϑ* + δ̂, with the covariance unchanged. It is optional, and the default should match INLAvaan.
3. **2. Skew-normal marginals per parameter.** Scan the height along the CMP direction (Lemma 3.1), add the volume correction (Corollary 3.3), then fit SN(ξ, ω, α). Back-transform the densities with the Jacobian (log → f(log ψ)/ψ; tanh → f(artanh ρ)/(1 − ρ²)).
4. **3. Gaussian copula draws.** Set R = cor(Ω), draw z ~ N(0, R), take u = Φ(z) and apply the SN quantile per margin. Return draws on both the unconstrained and natural scales.

Outputs form a minimal S3 class with `print()`, `summary()` (natural-scale moments, quantiles and mode, plus the log marginal likelihood), `marginals()`, `sample()`/`draws()` (optionally of a derived `g(theta)`), and `plot()`.

Out of scope for v0.1: factor scores, PPP/DIC, multiple chains, and compiled code.

## Source-of-truth rules for the maths

- Port and adapt the implementations in the INLAvaan source (github.com/haziqj/INLAvaan): joint Gaussian, VB correction, SN fitting with the volume correction, and copula sampling. Do not re-derive them from generic INLA knowledge. The references are the papers arXiv 2603.25690 and 2604.00671.
- If INLAvaan's code and `lab/IDEA.md` disagree, flag the disagreement and ask the user. Do not guess.
- The first step in `lab/IDEA.md` is a written plan mapping INLAvaan functions to INLAvanilla functions, shown to the user before any coding.

## Testing philosophy

The tests are described as "the main point of the package." They use closed-form or cheaply computed posteriors as ground truth: normal–normal conjugate, normal with unknown variance (log), beta–binomial on the logit scale, gamma–Poisson (log), a bivariate case with a correlation (tanh), analytic marginal likelihoods, analytic vs numerical gradients, and copula draws reproducing the marginal quantiles and R. Use tolerances that fit the approximation error, not exact equality. Where there is no closed form, compare against 1–2D numerical quadrature rather than MCMC.

## Known scaffold inconsistencies

- `DESCRIPTION` still has template placeholders (Title, Description, Authors@R).
- `DESCRIPTION` declares `LinkingTo: Rcpp` / `Imports: Rcpp` and a `src/` directory exists, but the spec puts compiled code out of scope. Clarify with the user before adding C++ or removing Rcpp. Planned Imports are numDeriv and possibly sn; keep dependencies minimal.
- The license is GPL (>= 3).
- `CONTRIBUTING.md` asks for a `NEWS.md` bullet for user-facing changes, but there is no `NEWS.md` yet.

## Workflow conventions

Make small, reviewable commits. Open a draft PR rather than pushing large changes to main. Follow the tidyverse style guide.
