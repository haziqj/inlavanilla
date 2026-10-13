# Summarise an INLAvanilla fit

Posterior summaries from the skew-normal marginals, on the unconstrained
scale of `logpost`, together with diagnostics of each stage of the
approximation.

## Usage

``` r
# S3 method for class 'inla_vanilla'
summary(object, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), ...)
```

## Arguments

- object:

  An object of class `inla_vanilla`.

- probs:

  Probabilities for the reported quantiles.

- ...:

  Unused.

## Value

An object of class `summary.inla_vanilla`: a list with the summary table
`table`, the log marginal likelihood `logml` (and its Laplace part
`logml_laplace`), the per-parameter diagnostics `diagnostics`, and
details of the optimiser, the VB correction, the number of calls and the
timing.

## Details

The diagnostics, and the levels at which they are flagged, are:

- **Mode.** The largest gradient at the mode in whitened units: how far,
  in posterior SDs, a Newton step would still move it. The mode is
  accepted when this is below 0.01, whatever the return code of
  [`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html), whose
  message is shown otherwise.

- **VB correction.** The shift of each mean from the mode, in posterior
  SDs (`vb_shift`). It is flagged when it exceeds 1 SD, or when the
  fixed-point iteration failed and the fallback optimiser ran: both mean
  the posterior is far from the Laplace Gaussian.

- **Tails.** Each parameter is scanned over \\\pm 4\\ Laplace SDs.
  `tails` is `FALSE` (and flagged) when the density at an end of the
  scan is still over 5% of its peak: the Laplace SD badly understates a
  heavy tail there, and the skew-normal extrapolates it.

- **Scan path.** `bend` is the log-density, in nats, that the straight
  scan misses at \\\pm 2\\ posterior SDs: the gain from one Newton step
  towards the conditional mode of the other parameters. It is flagged
  above 0.1. A large value means the other parameters' conditional mode
  bends away from the line, as a poorly determined noise scale does for
  a regression coefficient at small n, and the marginal comes out too
  narrow, by roughly 0.35 to 0.5 times `bend` in relative SD (see the
  Accuracy section of
  [`inla_vanilla()`](http://haziqj.ml/inlavanilla/reference/inla_vanilla.md)).

- **Skew-normal fits.** `nmad`, the largest gap between the fitted and
  the scanned density relative to its peak (flagged above 0.05), and
  `sn_vs_vb`, the gap between the skew-normal mean and the VB-corrected
  mean in posterior SDs. Both means target the posterior mean, so they
  should agree closely (to \\O(n^{-2})\\; Jamil & Rue, 2026, Proposition
  3.4). A gap above 0.1 SDs is flagged: the posterior is far from
  Gaussian along that parameter, and either approximation may be off.

- **Scan-end mass.** `scan_end_mass`, the mass each fitted skew-normal
  puts outside the scanned window (the mode \\\pm 4\\ Laplace SDs),
  where the marginal is extrapolated, as INLAvaan's `diagnostics()`
  reports it. A Gaussian marginal gives \\2\Phi(-4) = 6.3 \times
  10^{-5}\\; the larger it is, the more of the marginal is extrapolated.
  It is flagged above 0.05.

Some problems also raise one warning at the end of the fit, of class
`inlavanilla_diagnostics_warning`: heavy tails, an NMAD above 0.1, a
scan-end mass above 0.05, a VB shift above 1 SD or the VB fallback.
These levels match INLAvaan's fit-time warning; the flags here are more
sensitive. With `silent = TRUE` the warning is only recorded in
`fit$warnings`.
[`vignette("diagnostics", package = "inlavanilla")`](http://haziqj.ml/inlavanilla/articles/diagnostics.md)
explains every check with examples.

With `vb = FALSE` there is no VB mean, and the VB columns are left out.
With `marginals = "gaussian"` there are no scans, so only the VB shift
is reported.

## References

Jamil, H., & Rue, H. (2026). *Approximate Bayesian inference for
structural equation models using integrated nested Laplace
approximations* (arXiv:2603.25690 \[stat.ME\]). arXiv.
<https://doi.org/10.48550/arXiv.2603.25690>
