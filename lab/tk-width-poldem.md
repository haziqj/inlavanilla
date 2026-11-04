# The width term on the Political Democracy SEM

Date: 2026-11-01. Reproduce with `Rscript lab/tk-width-poldem.R` from the
package root (about 1 minute, most of it HMC; the draws are cached in
`lab/tk-width-poldem_cache/`, which is gitignored). Figures are in
`lab/tk-width-poldem-figures/`.

## Question

Each skew-normal marginal is fitted to height(t) + width(t) along the
conditional mode path (CMP) θ(t) = θ* + t vⱼ, where the width term
γⱼ(t) = −½ log|H₋ⱼ| is the log-volume of the other parameters. Three levels
of approximation, each removing one error, are now `inla_vanilla(width = )`:

1. `"tk"`: Tierney–Kadane. Re-optimise θ₋ⱼ at each grid point (warm start at
   the CMP point) and take height and width at the conditional mode.
2. `"cmp"`: no re-optimisation. Take the exact log-determinant of H₋ⱼ at each
   CMP grid point.
3. `"trace"` (default, INLAvaan's): the linear term γⱼ′t.

Trace against cmp is the linearisation error. Cmp against tk is the error of
the straight path. `tk_fallback = TRUE` uses trace everywhere except where
the ridge bend is flagged (above 0.25), and rescans those parameters with tk.

## Setup

- Model, parameterisation and priors as in `vignettes/articles/rcpp.qmd`
  (INLAvaan's defaults), with its C++ log-posterior and analytic gradient.
  The script extracts the C++ from the article. 39 parameters, n = 75.
- Every fit uses `grad`, the VB correction and the default 21-point grid on
  ±4 Laplace SDs. The skew-normal fit is the same for all methods.
- **Reference:** HMC in the whitened Laplace coordinates (the sampler of
  `lab/sem-hmc.R`). It ran 4 chains of 2,000 warm-up and 25,000 kept draws
  (100,000 draws, 24 s). Max R-hat 1.0002, min bulk ESS 29,792, min tail ESS
  21,257, and 4 rejected trajectories with an energy error over 1,000.
  `lab/DESIGN.md` records that this sampler agrees with Stan NUTS on
  `lab/poldem.stan`. Stan was not rerun here.
- **Metrics,** per parameter on the unconstrained scale that inlavanilla
  reports:
  - **JS%:** 100 × JS / log 2, with the MCMC KDE taken as in INLAvaan's
    `compare_mcmc()`. The grid covers both densities, so an MCMC tail beyond
    a narrow skew-normal counts. `compare_mcmc()` uses the approximation's
    own grid, so it would understate those tails. JS is invariant under the
    monotone maps to the natural scale, apart from KDE error.
  - **Δmean:** |mean − MCMC mean| / MCMC SD. This is the metric behind
    INLAvaan's ablation figure of 0.035 (skew-normal row, natural scale,
    blavaan reference, 42 parameters). It is also reported on the natural
    scale (variances, residual correlations).
  - **KS distance and SD ratio.**
- **MCMC noise floor of JS:** chains 1–2 against chains 3–4. It is 0.025%
  on average, and at most 0.075%.
- **Costs:** the marginal stage only, as the median of 5 fits. Calls are the
  fit's counts minus those of a `marginals = "gaussian"` fit with the same
  mode and VB correction. One `grad_cpp()` call takes about 10 µs, so these
  times are mostly R overhead. With an expensive gradient, time scales with
  the calls.

## Results

Aggregates are means over parameters. The bulk is the 37 parameters other
than `logv_x2` and `logv_dem65`, the two boundary variances whose log-scale
posteriors have long exponential left tails (`lab/DESIGN.md`).

| Method | JS% (all 39) | JS% median | KS | Δmean | \|SD ratio − 1\| | Δmean, natural scale | JS% (bulk 37) | Δmean (bulk) | \|SD ratio − 1\| (bulk) | Marginal time (s) | grad calls | logpost calls |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| trace | 0.380 | 0.057 | 0.019 | 0.032 | 0.061 | 0.025 | 0.081 | 0.016 | 0.034 | 0.08 | 3,198 | 780 |
| cmp | 0.253 | 0.062 | 0.016 | 0.026 | 0.057 | 0.019 | 0.087 | 0.015 | 0.036 | 0.53 | 59,358 | 780 |
| tk | 0.177 | 0.016 | 0.011 | 0.022 | 0.027 | 0.015 | 0.025 | 0.012 | 0.006 | 0.63 | 64,038 | 5,488 |
| trace + fallback | 0.291 | 0.057 | 0.017 | 0.028 | 0.057 | 0.022 | 0.081 | 0.016 | 0.034 | 0.08 | 4,773 | 937 |

The JS% between the methods themselves (mean over 39, median) is 0.093 and
0.003 for trace against cmp, and 0.055 and 0.035 for cmp against tk. Trace
against cmp is 2.14 for `logv_dem65` and 1.05 for `logv_x2`. Cmp against tk
is at most 0.36 (`z_y4_y8`).

**Flags.** The ridge bend flags only `logv_dem65` (0.44). The fallback
rescanned only it, for 1,575 extra gradient calls (+49%). Every TK
conditional optimisation succeeded (no grid point lost).

### Per parameter

JS% against HMC, the mean error in MCMC SDs, and the skew-normal SD over the
MCMC SD. The fallback column equals trace except for `logv_dem65`.

| Parameter | Ridge bend | JS trace | JS cmp | JS tk | JS fallback | Δmean trace | Δmean tk | SD ratio trace | SD ratio cmp | SD ratio tk |
|---|---|---|---|---|---|---|---|---|---|---|
| nu_x1 | 0.05 | 0.014 | 0.018 | 0.007 | 0.014 | 0.002 | 0.003 | 0.99 | 0.99 | 1.00 |
| nu_x2 | 0.05 | 0.031 | 0.038 | 0.011 | 0.031 | 0.005 | 0.005 | 0.98 | 0.98 | 1.00 |
| nu_x3 | 0.05 | 0.033 | 0.043 | 0.007 | 0.033 | 0.005 | 0.006 | 0.98 | 0.97 | 1.00 |
| nu_y1 | 0.05 | 0.045 | 0.033 | 0.010 | 0.045 | 0.002 | 0.002 | 0.97 | 0.98 | 1.00 |
| nu_y2 | 0.04 | 0.060 | 0.055 | 0.008 | 0.060 | 0.004 | 0.004 | 0.97 | 0.97 | 1.00 |
| nu_y3 | 0.05 | 0.051 | 0.048 | 0.011 | 0.051 | 0.002 | 0.002 | 0.97 | 0.97 | 1.00 |
| nu_y4 | 0.05 | 0.071 | 0.065 | 0.020 | 0.071 | 0.005 | 0.006 | 0.96 | 0.97 | 1.00 |
| nu_y5 | 0.05 | 0.053 | 0.035 | 0.016 | 0.053 | 0.000 | 0.000 | 0.97 | 0.98 | 1.00 |
| nu_y6 | 0.05 | 0.069 | 0.058 | 0.011 | 0.069 | 0.002 | 0.003 | 0.96 | 0.97 | 1.00 |
| nu_y7 | 0.05 | 0.069 | 0.052 | 0.012 | 0.069 | 0.004 | 0.004 | 0.96 | 0.97 | 0.99 |
| nu_y8 | 0.05 | 0.057 | 0.043 | 0.010 | 0.057 | 0.004 | 0.004 | 0.97 | 0.97 | 1.00 |
| lam_x2 | 0.04 | 0.182 | 0.076 | 0.066 | 0.182 | 0.076 | 0.054 | 0.95 | 0.99 | 0.99 |
| lam_x3 | 0.08 | 0.185 | 0.142 | 0.016 | 0.185 | 0.026 | 0.003 | 0.93 | 0.94 | 0.99 |
| a | 0.08 | 0.116 | 0.117 | 0.012 | 0.116 | 0.022 | 0.002 | 0.95 | 0.95 | 0.99 |
| b | 0.11 | 0.162 | 0.179 | 0.022 | 0.162 | 0.019 | 0.003 | 0.94 | 0.93 | 0.99 |
| c | 0.08 | 0.159 | 0.098 | 0.035 | 0.159 | 0.029 | 0.009 | 0.94 | 0.95 | 0.98 |
| beta_dem60_ind60 | 0.08 | 0.048 | 0.057 | 0.018 | 0.048 | 0.009 | 0.019 | 0.97 | 0.97 | 1.00 |
| beta_dem65_ind60 | 0.07 | 0.053 | 0.023 | 0.021 | 0.053 | 0.010 | 0.015 | 0.97 | 0.99 | 1.00 |
| beta_dem65_dem60 | 0.04 | 0.023 | 0.025 | 0.031 | 0.023 | 0.013 | 0.013 | 1.00 | 1.02 | 1.02 |
| logv_ind60 | 0.04 | 0.114 | 0.100 | 0.016 | 0.114 | 0.039 | 0.024 | 0.96 | 0.96 | 1.00 |
| logv_dem60 | 0.06 | 0.109 | 0.087 | 0.008 | 0.109 | 0.020 | 0.009 | 0.95 | 0.96 | 1.00 |
| **logv_dem65** | **0.44** | 5.937 | 2.795 | 2.475 | 2.475 | 0.339 | 0.167 | 0.46 | 0.61 | 0.64 |
| logv_x1 | 0.14 | 0.319 | 0.202 | 0.186 | 0.319 | 0.073 | 0.087 | 0.94 | 0.98 | 1.00 |
| logv_x2 | 0.11 | 5.875 | 3.853 | 3.482 | 5.875 | 0.328 | 0.230 | 0.40 | 0.48 | 0.51 |
| logv_x3 | 0.01 | 0.063 | 0.049 | 0.070 | 0.063 | 0.052 | 0.058 | 1.00 | 1.01 | 1.01 |
| logv_y1 | 0.03 | 0.101 | 0.077 | 0.042 | 0.101 | 0.002 | 0.008 | 0.93 | 0.94 | 0.96 |
| logv_y2 | 0.02 | 0.018 | 0.034 | 0.010 | 0.018 | 0.006 | 0.010 | 0.99 | 0.98 | 1.00 |
| logv_y3 | 0.02 | 0.009 | 0.018 | 0.008 | 0.009 | 0.000 | 0.004 | 0.99 | 0.99 | 1.00 |
| logv_y4 | 0.04 | 0.180 | 0.148 | 0.051 | 0.180 | 0.027 | 0.014 | 0.93 | 0.94 | 0.99 |
| logv_y5 | 0.03 | 0.039 | 0.037 | 0.025 | 0.039 | 0.011 | 0.015 | 0.98 | 0.98 | 0.99 |
| logv_y6 | 0.02 | 0.027 | 0.033 | 0.008 | 0.027 | 0.003 | 0.002 | 0.98 | 0.98 | 1.00 |
| logv_y7 | 0.02 | 0.014 | 0.015 | 0.015 | 0.014 | 0.009 | 0.014 | 1.00 | 1.00 | 1.01 |
| logv_y8 | 0.03 | 0.133 | 0.128 | 0.032 | 0.133 | 0.026 | 0.016 | 0.95 | 0.95 | 0.99 |
| z_y1_y5 | 0.05 | 0.034 | 0.062 | 0.017 | 0.034 | 0.016 | 0.017 | 0.98 | 0.97 | 1.00 |
| z_y2_y4 | 0.09 | 0.037 | 0.146 | 0.009 | 0.037 | 0.014 | 0.002 | 0.98 | 0.94 | 1.00 |
| z_y2_y6 | 0.12 | 0.044 | 0.190 | 0.007 | 0.044 | 0.005 | 0.001 | 0.97 | 0.93 | 0.99 |
| z_y3_y7 | 0.06 | 0.056 | 0.103 | 0.007 | 0.056 | 0.011 | 0.006 | 0.97 | 0.95 | 1.00 |
| z_y4_y8 | 0.13 | 0.200 | 0.449 | 0.072 | 0.200 | 0.026 | 0.001 | 0.93 | 0.89 | 0.98 |
| z_y6_y8 | 0.08 | 0.032 | 0.134 | 0.005 | 0.032 | 0.015 | 0.003 | 0.98 | 0.95 | 1.00 |

### Figures

- `width-curves.png`: γⱼ(t) − γⱼ(0) for trace, cmp and tk, per parameter.
- `height-curves.png`: the log-posterior along the CMP and at the conditional
  modes (TK).
- `marginals-vs-hmc.png`: the three skew-normal marginals over the HMC KDE.
- `js-by-parameter.png`: JS% per parameter and method, with each
  parameter's MCMC noise floor.

## Conclusions

- **The linearisation error (trace against cmp) is negligible for 37 of the
  39 parameters.** The median JS% between the two marginals is 0.003, against
  an error against MCMC of about 0.06. It matters only for the two boundary
  variances, where the width curve bends strongly (JS% 2.1 and 1.1 between
  the methods), and there cmp moves the marginal towards MCMC (JS% 5.9 → 2.8
  and 5.9 → 3.9). Elsewhere the curves part only beyond about ±2 SDs. There,
  the exact CMP width is often worse than the line: for the residual
  correlations, cmp has JS% up to 0.45 against 0.20 for trace. The linear
  term partly cancels the straight-path error, and the exact CMP width
  removes that cancellation.
- **The straight-path error (cmp against tk) is small per parameter but
  systematic, and it accounts for nearly all the error left in the bulk.**
  Trace and cmp marginals are 0–7% too narrow on all 37 bulk parameters (SD
  ratios 0.93–1.00 and 0.89–1.02). TK puts them at 0.96–1.02, with a mean
  |SD ratio − 1| of 0.006 against 0.034. Its bulk JS% (0.025) equals the
  MCMC noise floor (0.025), so at this MCMC precision TK cannot be told apart
  from HMC there. Bulk Δmean barely moves (0.016 → 0.012): the error is in
  the spread, not the location.
- **The overall error against MCMC is dominated by the skew-normal, not the
  width term.** `logv_dem65` and `logv_x2` give 80% of trace's summed JS (11.8
  of 14.8). TK roughly halves their JS (5.9 → 2.5 and 3.5), but
  their SD ratios stay at 0.64 and 0.51, because a skew-normal fitted on
  ±4 Laplace SDs cannot follow the exponential left tail. The Laplace
  marginal itself follows it (`lab/sem-marginal-check.R`), so a wider or
  adaptive grid, or a non-skew-normal fit, would be the next step for these
  two, not a better width.
- **The ridge bend diagnostic is a weak trigger here.** It flags only
  `logv_dem65` (0.44). The fallback fixes part of that error (JS 5.9 → 2.5,
  overall JS 0.38 → 0.29) for +49% gradient calls. It misses `logv_x2`
  (bend 0.11) and the uniform 3–7% narrowing of the bulk, where bends are
  0.01–0.14. So `tk_fallback = TRUE` is a cheap partial fix, and
  `width = "tk"` is the fix.
- **Cost.** On the marginal stage, per fit: trace takes 3,198 gradient and
  780 logpost calls (0.08 s). Cmp takes 59,358 gradient calls, 19× trace,
  for (n_grid − 1) · 2(m − 1) per parameter (0.53 s). TK takes 64,038
  gradient and 5,488 logpost calls (0.63 s). Its conditional optimisations
  are cheap (about 6 gradient and 7 logpost calls per grid point, from the
  CMP warm start in whitened coordinates), so nearly all its cost is the
  2(m − 1)-call Hessian per point, the same as cmp's. With this C++
  gradient (about 10 µs), TK is 8× the time of the trace marginals and under
  1 s in all. With an expensive gradient, expect about 17× the marginal
  stage (69,526 calls against 3,978). Since TK costs barely more than cmp, cmp is not worth using except
  as a diagnostic of the linearisation.
- **Comparison with INLAvaan's ablation (Δmean 0.035).** Trace here gives a
  natural-scale Δmean of 0.025, and TK 0.015. These are not directly
  comparable: the reference is HMC rather than blavaan, there are 39
  parameters rather than 42, and the grid differs. The ordering is what
  carries over.

## Caveats

- One data set and one model. The bulk improvement from TK is small in
  absolute terms (JS 0.08 → 0.025%), and only the HMC precision of 100,000
  draws resolves it.
- The JS grid differs from `compare_mcmc()`'s (see Setup). JS on the natural
  scale was not computed separately.
- The four rejected HMC trajectories are counted with a crude rule
  (|ΔH| > 1,000 or a non-finite gradient). They did not affect R-hat or ESS.
